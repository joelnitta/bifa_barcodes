#' Check monophyly of species in a tree
#'#'
#' @param tree Input phylogeny; list of class phylo.
#' @param outgroup Name of sinlge tip to use as outgroup.
#' @param sp_delim Character used to separate species name from voucher data
#' in tree tip labels
#' @param dataset Name of dataset analyzed
#'
#' @return Dataframe (tibble) with monophyly status for each species
check_monophy <- function(
  tree, outgroup = "Zygnema_circumcarinatum", sp_delim = "__", dataset) {

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
    janitor::clean_names() %>%
    mutate(dataset = dataset)

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



#' Write DNA seqeunces to PHYLIP
#'
#' @param x Object of class DNAbin.
#' @param file Path to write sequences.
#' @param ... Other args not used by this function.
#'
#' @return Path to written file
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

# Drop sequences that are not identified to species
drop_indets <- function(seq_list) {
  to_drop <- str_detect(names(seq_list), "sp.__")
  seq_list[!to_drop]
}

#' Drop species complexes and hybrids from sequences
#'
#' @param seq_list DNA sequences; List of class DNAbin.
#' @param sp_complex Vector of species names that are species complexes.
#' @param hybrid_taxa Vector of specimens that are hybrids.
#'
#' @return DNA sequences with species complexes and hybrid taxa removed
drop_complex_hybrids <- function(seq_list, sp_complex, hybrid_taxa) {
  is_sp_complex <- str_detect(
    names(seq_list), paste(sp_complex, collapse = "|"))
  is_hybrid <- names(seq_list) %in% hybrid_taxa
  assertthat::assert_that(
    isTRUE(sum(is_hybrid) <= length(hybrid_taxa))
  )
  to_drop <- is_sp_complex | is_hybrid
  seq_list[!to_drop]
}

#' Build a BLAST database.
#'
#' This is a wrapper for makeblastdb.
#'
#' @param in_seqs Character vector of length one; the path to the fasta
#' file containing the sequences to be included in the database.
#' @param db_type Character vector of length one; "nucl" for DNA or "prot"
#' for amino acids (proteins).
#' @param out_name Character vector of length one; name of BLAST database
#' to be created (optional). This will be used to name all database files;
#' if omitted, the name of the `in_seqs` file will be used instead.
#' @param title Character vector of length one; title of BLAST database
#' to be created (optional).
#' @param parse_seqids Logical; should the makeblastdb flag
#' "parse_seqids" be used?
#' @param wd Character vector of length one; working directory. The blast
#' database will be made here.
#' @param echo Logical; should standard error and output be printed?
#' @param ... Additional other arguments. Not used by this function, but
#' meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return A series of files starting with \code{out_name} and ending in
#' .phr, .pin, .pog, .psd, .psi, psq (for proteins) or .nhr, .nin, .nog,
#' .nsd, .nsi, and .nsq (for DNA) that constitute the BLAST database in
#' the working directory.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @examples
#' library(ape)
#' data(woodmouse)
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#' list.files(temp_dir)
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   title = "test db",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#' list.files(temp_dir)
#' fs::file_delete(temp_dir)
#' @export
build_blast_db <- function(in_seqs,
                            db_type = "nucl",
                            out_name = NULL,
                            title = NULL,
                            parse_seqids = FALSE,
                            echo = TRUE,
                            wd, ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(in_seqs))
  assertthat::assert_that(assertthat::is.string(db_type))
  assertthat::assert_that(assertthat::is.string(wd))
  assertthat::assert_that(is.logical(parse_seqids))
  assertthat::assert_that(assertthat::is.string(title) | is.null(title))
  assertthat::assert_that(assertthat::is.string(out_name) | is.null(out_name))

  assertthat::assert_that(db_type %in% c("nucl", "prot"))

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  in_seqs <- fs::path_abs(in_seqs)
  assertthat::assert_that(assertthat::is.readable(in_seqs))

  # Prepare arguments
  parse_seqids <- if(isTRUE(parse_seqids)) "-parse_seqids" else NULL
  title <- if(!is.null(title)) c("-title", title) else NULL
  out_name <- if(!is.null(out_name)) c("-out", out_name) else NULL

  arguments <- c("-in", in_seqs,
                 "-dbtype", db_type,
                 parse_seqids,
                 out_name,
                 title)

  # run command
  processx::run("makeblastdb", arguments, wd = wd, echo = echo)

}

write_fasta <- function(x, file, ...) {
  ape::write.FASTA(x = x, file = file, ...)
  file
}

#' Make a blast database and return paths to the created files
#'
#' @param seqs_path Path to sequences
#' @param blast_db_dir Folder to write BLAST database
#' @param out_name Name of BLAST database
#'
#' @return Paths to components of BLAST database.
#'   Externally, database will be created
#'
make_blast_db <- function(seqs_path, blast_db_dir, out_name) {

  # Define output files
  out_files <- fs::path(blast_db_dir, out_name) %>%
    paste0(c(".nhr", ".nin", ".nsq"))

  # Delete any existing output
  for (i in seq_along(out_files)) {
    if(fs::file_exists(out_files[[i]])) {
      fs::file_delete(out_files[[i]])
    }
  }

  # Create blast DB (side-effect)
  build_blast_db(
    seqs_path,
    title = out_name,
    out_name = out_name,
    parse_seqids = FALSE,
    wd = blast_db_dir)

  # Return path to output file
  out_files

}

#' Run a blastn query.
#'
#' This is a wrapper for blastn.
#'
#' @param query Character vector of length one; the path to the fasta
#' file to use as the query sequence(s).
#' @param database Character vector of length one; the name of the blast
#' database.
#' @param out_file Character vector of length one; the name to use for
#' the results file.
#' @param outfmt Character vector of length one; value to pass to
#' \code{blastn} \code{outfmt} argument. Default = "6".
#' @param other_args Character vector; other arguments to pass on to
#' \code{blastn}.
#' Must be formatted so that each argument name and its value are
#' separate, consecutive elements of the vector, e.g.,
#' \code{c("-evalue", 10, "-num_threads", 1)}.
#' The argument name must be preceded by a hyphen.
#' For a list of options, run \code{blastn -help}.
#' @param wd Character vector of length one; working directory. The blast
#' search will be conducted here.
#' @param echo Logical; should standard error and output be printed?
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return A tab-separated text file with the results of the blastn
#' search, named with the value of \code{out_file}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @examples
#' library(ape)
#'
#' # Make temp dir for storing files
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Write out ape::woodmouse dataset as DNA
#' data(woodmouse)
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#'
#' # Make blast database
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   db_type = "nucl",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#'
#' # Blast the original sequences against the database
#' blast_n(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   database = "wood",
#'   out_file = "blastn_results",
#'   wd = temp_dir,
#'   echo = TRUE
#' )
#'
#' # Take a look at the results.
#' readr::read_tsv(
#'   fs::path(temp_dir, "blastn_results"),
#'   col_names = FALSE
#'   )
#'
#' # Cleanup.
#' fs::file_delete(temp_dir)
#' @export
blast_n <- function(query,
                     database,
                     out_file = NULL,
                     outfmt = "6",
                     other_args = NULL,
                     echo = TRUE,
                     wd,
                     ...) {

  # Check input

  assertthat::assert_that(assertthat::is.string(query))
  assertthat::assert_that(assertthat::is.string(database))
  assertthat::assert_that(assertthat::is.string(out_file) | is.null(out_file))
  assertthat::assert_that(assertthat::is.string(outfmt))
  assertthat::assert_that(is.character(other_args) | is.null(other_args))
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.string(wd))
  assertthat::assert_that(
    length(other_args) > 1 | is.null(other_args),
    msg = "other_args not formatted correctly.")

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  query <- fs::path_abs(query)
  assertthat::assert_that(assertthat::is.readable(query))

  # modify arguments

  arguments <- c("-query", query,
                 "-db", database,
                 "-outfmt", outfmt,
                 "-out", out_file,
                 other_args)

  # run command
  processx::run("blastn", arguments, wd = wd, echo = echo)

  # return output file
  fs::path(wd, out_file)

}

load_blast_tsv <- function(blast_out_path, dataset_name) {
  # Read in sequences.
  # BLAST doesn't output column headers, so we need to specify
  # (make sure they match correctly first!)
  fmt6_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  # Read in BLAST output
  readr::read_tsv(
    blast_out_path,
    col_names = fmt6_cols,
    col_types = "ccdddddddddd" # two ID cols are char, rest is numeric
  ) %>%
  mutate(dataset = dataset_name)
}

prep_blast_res <- function(blast_res_raw, q_seqs) {

  # Get lengths (bp) of query seqs
  q_seqs <- ape::del.gaps(q_seqs)
  seq_lengths <-
    tibble(
      qseqid = names(q_seqs),
      q_seq_len = map_dbl(q_seqs, length)
    )

  blast_res_raw %>%
    left_join(seq_lengths, by = "qseqid", relationship = "many-to-one") %>%
    assert(not_na, q_seq_len) %>%
    mutate(
      q_species = str_split(qseqid, "__") %>% map_chr(1)) %>%
    mutate(
      s_species = str_split(sseqid, "__") %>% map_chr(1)) %>%
    filter(length >= 0.95 * q_seq_len) %>%
    group_by(q_species) %>%
    mutate(n_indiv = n_distinct(qseqid)) %>%
    ungroup() %>%
    mutate(
      comp_type = case_when(
        qseqid == sseqid ~ "self",
        q_species == s_species ~ "intra",
        q_species != s_species ~ "inter"
      )
    )
}

get_species_cutoff <- function(blast_res, dataset_select) {
  blast_res %>%
    filter(comp_type == "intra") %>%
    filter(dataset == dataset_select) %>%
    summarize(
      mean_pident = mean(pident),
      min_pident = min(pident),
      .groups = "drop") %>%
    mutate(dataset = dataset_select) %>%
    pivot_longer(
      names_to = "cutoff_type",
      -dataset
    )
}

test_blast <- function(
  blast_res,
  dataset_select, cutoff_table, cutoff_type_select = "mean") {

  # Obtain value to use for infraspecific cutoff
  intra_cutoff <-
    cutoff_table %>%
    filter(str_detect(dataset, "no_hybrids")) %>%
    filter(str_detect(dataset, dataset_select)) %>%
    filter(str_detect(cutoff_type, cutoff_type_select)) %>%
    pull(value)

  blast_res %>%
    filter(dataset == dataset_select) %>%
    mutate(
      above_cutoff = pident > intra_cutoff,
      match_other = q_species != s_species,
      fail = above_cutoff & match_other) %>%
    group_by(q_species) %>%
    summarize(fail = any(fail), .groups = "drop") %>%
    mutate(dataset = dataset_select)
}

summarize_blast_fail_rate <- function(blast_test_res_raw) {
  blast_test_res_raw %>%
    group_by(dataset) %>%
    count(fail) %>%
    ungroup()
}
