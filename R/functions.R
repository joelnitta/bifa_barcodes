#' Check monophyly of species in a tree
#'#'
#' @param tree Input phylogeny; list of class phylo.
#' @param outgroup Name of sinlge tip to use as outgroup.
#' @param sp_delim Character used to separate species name from voucher data
#' in tree tip labels
#'
#' @return Dataframe (tibble) with monophyly status for each species
check_monophy <- function(
  tree, outgroup = "Zygnema_circumcarinatum", sp_delim = "__") {

  assertthat::assert_that(assertthat::is.string(outgroup))
  assertthat::assert_that(assertthat::is.string(sp_delim))

  # Get list of outgroup taxa
  og_taxa <- get_og_taxa()

  # Root tree
  tree_rooted <-
    phytools::reroot(
        tree = tree,
        node.number = which(tree$tip.label == outgroup)
      )

  # Make dataframe mapping tips to species names
  taxonomy_data <- tibble(tip = tree_rooted$tip.label) %>%
    separate_wider_delim(
      tip,
      sp_delim,
      too_many = "drop", too_few = "align_start",
      names = c("species", "voucher"), cols_remove = FALSE) %>%
    select(tip, species) %>%
    as.data.frame()

  monophy_res <- MonoPhy::AssessMonophyly(tree_rooted, taxonomy = taxonomy_data)

  monophy_res[[1]]$result %>%
    rownames_to_column("species") %>%
    as_tibble() %>%
    filter(!species %in% og_taxa) %>%
    janitor::clean_names()

}

# Convert DNA sequences in DNAbin format to DNAStringSet
dnabin_to_dnass <- function(dnabin) {
  temp_file <- tempfile(fileext = ".fasta")
  ape::write.FASTA(dnabin, temp_file)
  seq <- Biostrings::readDNAStringSet(temp_file)
  fs::file_delete(temp_file)
  seq
}

analyze_dist <- function(alignment) {

  dist_matrix <-
    alignment %>%
    dnabin_to_dnass() %>%
    DECIPHER::DistanceMatrix()

  diag(dist_matrix) <- NA
  dist_matrix[upper.tri(dist_matrix)] <- NA

  dist_matrix %>%
    as.data.frame() %>%
    rownames_to_column("acc_1") %>%
    as_tibble() %>%
    pivot_longer(names_to = "acc_2", values_to = "dist", -acc_1) %>%
    filter(acc_1 != acc_2) %>%
    filter(!is.na(dist)) %>%
    separate_wider_regex(
      acc_1, c(sp_1 = "^.*?", "__", voucher_1 = ".*")
    ) %>%
    separate_wider_regex(
      acc_2, c(sp_2 = "^.*?", "__", voucher_2 = ".*")
    ) %>%
    # assert that all vouchers are different
    mutate(voucher_check = voucher_1 != voucher_2) %>%
    assert(in_set(TRUE), voucher_check) %>%
    select(-voucher_check) %>%
    mutate(
      comp_type = if_else(sp_1 == sp_2, "intra", "inter")
    )
}

write_phy <- function(x, file, ...) {
  ips::write.phy(x = x, file = file, ...)
  return(file)
}

#' Run IQ-TREE
#'
#' For details, see http://www.iqtree.org/doc/
#' 
#' Requires docker to be installed
#'
#' @param alignment DNA alignment to use for phylogenetic analysis. Must be
#'   matrix (i.e., aligned sequences) of class DNAbin
#' @param aln_path Path to DNA alignment. Either alignment, or aln_path must be
#'   provided, but not both
#' @param tree_path Optional; path to tree(s) written out by IQ-TREE,
#'   useful if this differs from default alignment name.
#' @param wd Path to working directory. The alignment and IQ-TREE intermediate
#'   files and results will be written here.
#' @param bb Optional; number of ultrafast bootstrap replicates to run.
#' @param nt Optional; number of cores to use. Set to "AUTO" to determine
#'   automatically.
#' @param alrt Optional; number of SH-aLRT tests to run.
#' @param m Optional; specify model. If no model is given, ModelTest will be run
#'   to identify the best model for the data.
#' @param redo Logical; should the analysis be redone from scratch if output
#'   from previous runs is present?
#' @param spp Path to partition file.
#' @param seed Optional; Specify a random number seed to reproduce a previous
#'   run.
#' @param echo Logical; should STDERR be written to the screen?
#' @param other_args Other arguments to pass to IQ tree; must be entered as a
#'   character vector with the name and value of each argument separately. For
#'   example, c("-pers", "0.2", "-nstop", "500").
#' @param ... Additional arguments; not used by this function.
#'
#' @return List; a single phylogenetic tree (list of class "phylo")
iqtree_docker <- function(alignment_path = NULL, wd = getwd(),
                   tree_path = NULL,
                   bb = NULL, nt = NULL, alrt = NULL, m = NULL, redo = FALSE,
                   spp = NULL,
                   seed = NULL,
                   echo = FALSE,
                   other_args = NULL, ...) {

  alignment_dir <- fs::path_dir(alignment_path) %>%
    fs::path_abs()
  alignment_file <- fs::path_file(alignment_path)
  out_tree <- fs::path_ext_set(alignment_path, ".phy.contree")

  # Set up arguments
  iqtree_arguments <- c(
    "-s", alignment_file,
    if (!is.null(bb)) "-bb",
    bb,
    if (!is.null(alrt)) "-alrt",
    alrt,
    if (!is.null(nt)) "-nt",
    nt,
    if (!is.null(m)) "-m",
    m,
    if (!is.null(seed)) "-seed",
    seed,
    if (!is.null(spp)) "-spp",
    fs::path_abs(spp),
    if (isTRUE(redo)) "-redo",
    other_args
  )

  babelwhale::run(
    container_id = "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0",
    command = "iqtree2",
    args = iqtree_arguments,
    volumes = glue::glue("{alignment_dir}:/wd"),
    workspace = "/wd"
  )

  ape::read.tree(out_tree)

}

#' Run IQ-TREE
#'
#' For details, see http://www.iqtree.org/doc/
#' 
#' Requires iqtree2 to be installed
#'
#' @param alignment DNA alignment to use for phylogenetic analysis. Must be
#'   matrix (i.e., aligned sequences) of class DNAbin
#' @param aln_path Path to DNA alignment. Either alignment, or aln_path must be
#'   provided, but not both
#' @param tree_path Optional; path to tree(s) written out by IQ-TREE,
#'   useful if this differs from default alignment name.
#' @param wd Path to working directory. The alignment and IQ-TREE intermediate
#'   files and results will be written here.
#' @param bb Optional; number of ultrafast bootstrap replicates to run.
#' @param nt Optional; number of cores to use. Set to "AUTO" to determine
#'   automatically.
#' @param alrt Optional; number of SH-aLRT tests to run.
#' @param m Optional; specify model. If no model is given, ModelTest will be run
#'   to identify the best model for the data.
#' @param redo Logical; should the analysis be redone from scratch if output
#'   from previous runs is present?
#' @param spp Path to partition file.
#' @param seed Optional; Specify a random number seed to reproduce a previous
#'   run.
#' @param echo Logical; should STDERR be written to the screen?
#' @param other_args Other arguments to pass to IQ tree; must be entered as a
#'   character vector with the name and value of each argument separately. For
#'   example, c("-pers", "0.2", "-nstop", "500").
#' @param ... Additional arguments; not used by this function.
#'
#' @return List; either a single phylogenetic tree (list of class "phylo"),
#' or list containing phylogenetic trees
#'
#' @examples
#' \dontrun{
#' library(ape)
#' temp_dir <- fs::dir_create(tempdir(), "blah")
#' data(woodmouse)
#' # Rapid boot-strap tree with 1000 replicates on best-fitting model
#' tree <- iqtree(woodmouse, temp_dir, bb = 1000, echo = TRUE)
#' plot(tree)
#' # Check the optimum number of cores to use for GTR+I+G model
#' iqtree(
#'   woodmouse,
#'   temp_dir,
#'   m = "GTR+I+G", nt = "AUTO", echo = TRUE, redo = TRUE)
#' fs::dir_delete(temp_dir)
#' }
iqtree <- function(alignment = NULL, wd = getwd(),
                   aln_path = NULL,
                   tree_path = NULL,
                   bb = NULL, nt = NULL, alrt = NULL, m = NULL, redo = FALSE,
                   spp = NULL,
                   seed = NULL,
                   echo = FALSE,
                   other_args = NULL, ...) {

  assertthat::assert_that(
    !is.null(alignment) | !is.null(aln_path),
    msg = "Either alignment or aln_path must be provided, but not both")

  assertthat::assert_that(
    is.null(alignment) | is.null(aln_path),
    msg = "Either alignment or aln_path must be provided, but not both")

  assertthat::assert_that(assertthat::is.dir(wd))

  assertthat::assert_that(is.logical(echo))

  assertthat::assert_that(is.logical(redo))

  if (!is.null(bb))
    assertthat::assert_that(assertthat::is.number(bb))

  if (!is.null(alrt))
    assertthat::assert_that(assertthat::is.number(alrt))

  if (!is.null(nt))
    assertthat::assert_that(
      assertthat::is.number(nt) | assertthat::is.string(nt))

  if (!is.null(m))
    assertthat::assert_that(assertthat::is.string(m))

  if (!is.null(spp))
    assertthat::assert_that(assertthat::is.readable(spp))

  if (!is.null(seed))
    assertthat::assert_that(assertthat::is.number(seed))

  wd <- fs::path_norm(wd)

  # check that iqtree is installed and on the PATH
  tryCatch({
    processx::run("iqtree2", "-h", echo = FALSE)
  }, warning = function(w) {
    stop("iqtree not installed and on path")
  }, error = function(e) {
    stop("iqtree not installed and on path")
  }, finally = {
    TRUE
  })

  # Write alignment to working directory in phylip format if alignment
  # is provided via R as DNAbin
  if (is.null(aln_path)) {
    assertthat::assert_that(inherits(alignment, "DNAbin"),
                            msg = "alignment must be of class 'DNAbin'")
    assertthat::assert_that(
      is.matrix(alignment),
      msg = "alignment must be a matrix (not a list of unaligned sequences)")
   
    aln_path <- fs::path(wd, deparse(substitute(alignment))) %>%
      fs::path_ext_set("phy")
   
    phangorn::write.phyDat(alignment, aln_path, format = "phylip")
  }

  assertthat::assert_that(assertthat::is.readable(aln_path))

  # Set up arguments
  iqtree_arguments <- c(
    "-s", fs::path_abs(aln_path),
    if (!is.null(bb)) "-bb",
    bb,
    if (!is.null(alrt)) "-alrt",
    alrt,
    if (!is.null(nt)) "-nt",
    nt,
    if (!is.null(m)) "-m",
    m,
    if (!is.null(seed)) "-seed",
    seed,
    if (!is.null(spp)) "-spp",
    fs::path_abs(spp),
    if (isTRUE(redo)) "-redo",
    other_args
  )

  # Run iqtree command
  processx::run(
    "iqtree2",
    iqtree_arguments, wd = wd, echo = echo,
    # Include env variable as workaround for initial parsimony analysis
    # using all cores
    # https://github.com/iqtree/iqtree2/issues/18
    env = c("current", OMP_NUM_THREADS = "1"))

  # Read in resulting tree(s)
  # Default: use default treefile if tree_path not provided
  if (is.null(tree_path)) {
    tree_path <- fs::path(wd, deparse(substitute(alignment))) %>%
      fs::path_ext_set(".phy.treefile")
  }

  # Return single tree if only one file in tree_path
  if(length(tree_path) == 1) {
    assertthat::assert_that(assertthat::is.readable(tree_path))
    res <- ape::read.tree(tree_path)
  }

  # Return list of trees if multiple files in tree_path
  if(length(tree_path) > 1) {
    # Set up results list to have same
    # names as tree_path
    res <- vector(length = length(tree_path))
    names(res) <- names(tree_path)
    res <- as.list(res)
    for(i in seq_along(tree_path)) {
      assertthat::assert_that(assertthat::is.readable(tree_path[[i]]))
      # although we check if file is readable
      # beware that ape will crash R if it is not a newick file!
      res[[i]] <- ape::read.tree(tree_path[[i]])
    }
  }

  return(res)

}

align_seqs <- function(seqs) {
  # Align
  alignment <- ips::mafft(
    x = seqs,
    options = "--adjustdirection",
    exec = system("which mafft", intern = TRUE)
  )

  # remove _R_ appended by mafft
  seq_names <- rownames(alignment)
  seq_names <- stringr::str_remove_all(seq_names, "_R_$")
  rownames(alignment) <- seq_names

  alignment

}

get_og_taxa <- function() {
  # get outgroup sequences used by ftol
  ftolr::accessions_wide %>%
    filter(outgroup) %>%
    pull(species)
}

align_with_og <- function(fern_seqs, marker) {
  
  og_taxa <- get_og_taxa()

  ftol_seqs <- ftolr::ft_seqs(loci = marker, aligned = FALSE, del_gaps = TRUE)

  og_seqs <- ftol_seqs[names(ftol_seqs) %in% og_taxa]

  # Combine outgroup seqs and fern seqs, align
  align_seqs(c(og_seqs, fern_seqs))
}

drop_indets <- function(seq_list) {
  to_drop <- str_detect(names(seq_list), "sp.__")
  seq_list[!to_drop]
}
