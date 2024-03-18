check_monophy <- function(tree) {

}

analyze_dist <- function(dna_ss) {
  dna_matrix <- DECIPHER::DistanceMatrix(dna_ss)

  diag(dna_matrix) <- NA
  dna_matrix[upper.tri(dna_matrix)] <- NA

  dna_matrix %>%
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
iqtree <- function(alignment_path = NULL, wd = getwd(),
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

align_with_og <- function(fern_seqs, marker) {
  # get outgroup sequences used by ftol
  og_taxa <-
    ftolr::accessions_wide %>%
    filter(outgroup) %>%
    pull(species)

  ftol_seqs <- ftolr::ft_seqs(loci = marker, aligned = FALSE)

  og_seqs <- ftol_seqs[og_taxa]

  # Combine outgroup seqs and fern seqs
  combined_seqs <- c(og_seqs, fern_seqs)

  # Align
  alignment <- ips::mafft(
    x = combined_seqs,
    options = "--adjustdirection",
    exec = system("which mafft", intern = TRUE)
  )

  # remove _R_ appended by mafft
  seq_names <- rownames(alignment)
  seq_names <- stringr::str_remove_all(seq_names, "_R_$")
  rownames(alignment) <- seq_names

  alignment
}
