# Filter data ----

# Drop sequences that are not identified to species
drop_indets <- function(seq_tbl) {
  seq_tbl %>%
    # drop `sp.`
    filter(!str_detect(voucher, "_sp_")) %>%
    # drop `aff.`
    filter(!str_detect(voucher, "_aff_")) %>%
    # drop `cf.`
    filter(!str_detect(voucher, "_cf_"))
}

#' Drop species complexes and hybrids from sequences
#'
#' @param seq_tbl Tibble of DNA sequences; `seq` is list of class DNAbin.
#' @param hybrid_taxa Tibble of specimens that are hybrids.
#' @param sp_complex Tibble of species names that are species complexes.
#'
#' @return DNA sequences with species complexes and hybrid taxa removed
drop_complex_hybrids <- function(seq_tbl, hybrid_taxa, sp_complex) {
  seq_tbl %>%
    anti_join(hybrid_taxa, by = "voucher") %>%
    anti_join(sp_complex, by = "species")
}

# Metadata ----

add_family <- function(seqs, ppgi) {
  seqs %>%
    mutate(
      genus = str_split_i(species, " ", 1)
    ) %>%
    left_join(
      select(ppgi, genus, family),
      by = "genus"
    ) %>%
    assert(not_na, family) %>%
    select(-genus)
}

load_rbcl_og_seqs <- function() {

  og_taxa <- get_og_taxa()

  ftolr::accessions_wide

  ftol_seqs <- ftolr::ft_seqs(loci = "rbcL", aligned = FALSE, del_gaps = TRUE)

  rbcl <- ftol_seqs[names(ftol_seqs) %in% og_taxa] %>%
    dnabin_to_seqtbl(name_col = "species") %>%
    left_join(
      select(ftolr::accessions_wide, species, voucher = plastome),
      by = "species"
    ) %>%
    mutate(
      voucher = glue::glue("{species}_{voucher}"),
      species = str_replace_all(species, "_", " "),
      family = "outgroup",
      marker = "rbcl"
    )
  rbcl_short <- mutate(rbcl, marker = "rbcl_short")
  bind_rows(rbcl, rbcl_short)

}

# Sequence alignment ----

remove_mafft_r <- function(alignment) {
  seq_names <- rownames(alignment)
  seq_names <- stringr::str_remove_all(seq_names, "_R_$")
  seq_names <- stringr::str_remove_all(seq_names, "^_R_")
  rownames(alignment) <- seq_names
  return(alignment)
}

align_seqs <- function(seqs) {
  # Align
  alignment <- ips::mafft(
    x = seqs,
    options = "--adjustdirection",
    exec = system("which mafft", intern = TRUE)
  )

  # remove _R_ appended by mafft
  remove_mafft_r(alignment)
}

align_with_og <- function(fern_seqs, marker) {

  og_taxa <- get_og_taxa()

  ftol_seqs <- ftolr::ft_seqs(loci = marker, aligned = FALSE, del_gaps = TRUE)

  og_seqs <- ftol_seqs[names(ftol_seqs) %in% og_taxa]

  # Combine outgroup seqs and fern seqs, align
  align_seqs(c(og_seqs, fern_seqs))
}

#' Align sequences in a tibble
#'
#' @param seqtbl Tibble containing one DNA sequence per row
#' in a list-column
#' @param name_col Name of column with sequence name
#' @param seq_col Name of column with sequences (list-column)
#'
#' @return Dataframe with aligned sequences. New column logical column "reversed" will
#' be appended indicating if sequence was reversed when aligning or not.
#'
align_seqs_tbl <- function(seqs_tbl, name_col = "accession", seq_col = "seq") {
  
  # Need to use quasiquotation for name_col
  name_col_sym <- ensym(name_col)
  
  # Extract sequences, convert to ape DNAbin list
  seqs <-
    seqs_tbl %>%
    # Name column must be unique values
    assert(is_uniq, !!name_col_sym) %>%
    seqtbl_to_dnabin(name_col = name_col, seq_col = seq_col)

  # Set aside metadata
  seqs_data <-
    select(seqs_tbl, -all_of(seq_col)) %>%
    # Metadata can't already contain column name "reversed"
    verify(!"reversed" %in% colnames(.))

  # Align sequences
  alignment <- ips::mafft(
    x = seqs,
    options = "--adjustdirection",
    exec = "/usr/bin/mafft")

  # Join metadata back to aligned sequences
  # - make tibble of which seqs were reversed
  reversed_seqs_tbl <-
    alignment %>%
    dnabin_to_seqtbl(name_col = name_col, seq_col = seq_col) %>%
    # remove '_R_' from species names inserted by mafft
    select(!!name_col_sym) %>%
    mutate(reversed = str_detect(!!name_col_sym, "_R_")) %>%
    mutate(across(-reversed, ~str_remove_all(., "_R_")))

  # - fix names if mafft changed them
  rownames(alignment) <- str_remove_all(rownames(alignment), "_R_")

  # Join metadata back to aligned sequences
  alignment %>%
    dnabin_to_seqtbl(name_col = name_col, seq_col = seq_col) %>%
    left_join(reversed_seqs_tbl, by = name_col) %>%
    left_join(seqs_data, by = name_col) %>%
    assert(not_na, reversed)
}

concatenate_trnlf_rbcl <- function(
  rbcl_aligned_to_cat,
  trnlf_aligned_to_cat) {

  trnlf_aligned_to_cat_dnabin <-
    trnlf_aligned_to_cat %>%
    seqtbl_to_dnabin(name_col = "voucher") %>%
    as.matrix()

  rbcl_aligned_to_cat_dnabin <-
    rbcl_aligned_to_cat %>%
    seqtbl_to_dnabin(name_col = "voucher") %>%
    as.matrix()

  metadata <- bind_rows(
    select(rbcl_aligned_to_cat, voucher, species, taxon_sampling, family),
    select(trnlf_aligned_to_cat, voucher, species, taxon_sampling, family)
  ) %>%
    unique() %>%
    assert(is_uniq, voucher) %>%
    mutate(marker = "rbcl_trnlf") %>%
    assert(not_na, everything())

  ape::cbind.DNAbin(
    rbcl_aligned_to_cat_dnabin,
    trnlf_aligned_to_cat_dnabin,
    fill.with.gaps = TRUE) %>%
    dnabin_to_seqtbl(name_col = "voucher") %>%
    left_join(metadata, by = "voucher", relationship = "one-to-one") %>%
    assert(not_na, everything())
  
}

#' Merge subalignments with MAFFT
#'
#' Modified from ips::mafft.merge(), but allows adding
# "singleton" seqs that are not part of any subalignment
#'
#' @param subMSA List of sub-alignments to merge; each one a matrix of class "DNAbin"
#' @param other_seqs List of class "DNAbin"; "singleton" sequences that don't belong
#' to any sub-alignment, but should be merged with the sub-alignments.
#' @param method Name of method for MAFFT to use
#' @param gt List of class "phylo"; guide tree
#' @param thread Number of threads to use
#' @param exec Path to MAFFT executable
#' @param quiet Logical; should MAFFT output be printed to screen?
#' @param adjustdirection Logical; should MAFFT attempt to automatically adjust sequence direction?
#'
#' @return Matrix of class "DNAbin"; the merged alignment
#'
mafft_merge <- function(subMSA, other_seqs, method = "auto", gt, thread = -1, exec, quiet = TRUE, adjustdirection = FALSE)
{
    quiet <- ifelse(quiet, "--quiet", "")
    adjustdirection <- ifelse(adjustdirection, "--adjustdirection", "")
    method <- match.arg(method, c("auto", "localpair", "globalpair",
        "genafpair", "parttree", "retree 1", "retree 2"))
    method <- paste("--", method, sep = "")
    if (missing(gt)) {
        gt <- ""
    }
    else {
        phylo2mafft(gt)
        gt <- "--treein tree.mafft"
    }
    n <- sapply(subMSA, nrow)
    subMSAtable <- vector(length = length(n))
    init <- 0
    for (i in seq_along(n)) {
        nn <- 1:n[i] + init
        init <- max(nn)
        subMSAtable[i] <- paste(nn, collapse = " ")
    }
    subMSA <- lapply(subMSA, as.list)
    subMSA <- c(subMSA, list(other_seqs))
    subMSA <- do.call(c, subMSA)
    names(subMSA) <- gsub("^.+[.]", "", names(subMSA))
    fns <- vector(length = 3)
    for (i in seq_along(fns)) fns[i] <- tempfile(pattern = "mafft",
        tmpdir = tempdir())
    write(subMSAtable, fns[1])
    ips::write.fas(subMSA, fns[2])
    call.mafft <- paste(exec, method, "--merge", fns[1], quiet,
      adjustdirection,
      gt, "--thread", thread, fns[2], ">", fns[3])
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- length(scan(fns[3], what = "c", quiet = TRUE))
    if (res != 0) {
        res <- ape::read.FASTA(fns[3])
        if (length(unique(sapply(res, length))) == 1) {
            res <- as.matrix(res)
        }
    }
    unlink(fns[file.exists(fns)])
    return(res)
}

merge_trnlf <- function(trnlf_seqs, trnlf_aligned_family_grouped) {
  seqtbl <- trnlf_aligned_family_grouped %>%
    assert(not_na, voucher) %>%
    assert(is_uniq, voucher)

  trnlf_singletons <- trnlf_seqs %>%
    assert(not_na, voucher) %>%
    assert(is_uniq, voucher) %>%
    anti_join(seqtbl, by = "voucher") %>%
    seqtbl_to_dnabin(name_col = "voucher")

  sub_msa_list <-
    seqtbl %>%
    group_by(family) %>%
    nest() %>%
    mutate(seqs = map(data, ~seqtbl_to_dnabin(.x, name_col = "voucher"))) %>%
    pull(seqs) %>%
    map(as.matrix)

  mafft_merge(
    sub_msa_list, trnlf_singletons, thread = 2,
    exec = "/usr/bin/mafft", adjustdirection = TRUE) %>%
    # remove _R_ from sequence names
    remove_mafft_r() %>%
    # convert to tbl
    dnabin_to_seqtbl(name_col = "voucher") %>%
    # add metadata
    left_join(
      select(trnlf_seqs, voucher, species, marker, family, taxon_sampling),
      by = "voucher",
      relationship = "one-to-one"
    ) %>%
    assert(not_na, everything())
}

# Barcoding gap test ----

calc_barcode_dist <- function(seqs_aligned) {

  # skip calculation if less than two rows (singleton family)
  if (nrow(seqs_aligned) < 2) {
    return(tibble())
  }

  taxon_sampling <- unique(seqs_aligned$taxon_sampling)
  assertthat::assert_that(assertthat::is.string(taxon_sampling))

  marker <- unique(seqs_aligned$marker)
  assertthat::assert_that(assertthat::is.string(marker))

  seqs_aligned_as_list <- seqtbl_to_dnabin(
    seqs_aligned, name_col = "voucher")

  seq_lens <- seqs_aligned_as_list %>%
    map_dbl(length) %>%
    unique()

  assertthat::assert_that(
    isTRUE(length(seq_lens) == 1),
    msg = "Aligned seqs not all the same length."
  )

  # Calculate distances as matrix
  seqs_dist <- 
    seqs_aligned_as_list %>%
    as.matrix() %>%
    ape::dist.dna(model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)
  # set diagonal (self-matches) to NA
  diag(seqs_dist) <- NaN
  # set upper triangle (reverse direction match) to NA
  seqs_dist[upper.tri(seqs_dist)] <- NaN

  # convert to dataframe
  as.data.frame(seqs_dist) %>%
    rownames_to_column("voucher_1") %>%
    pivot_longer(names_to = "voucher_2", values_to = "dist", -voucher_1) %>%
    filter(!is.na(dist)) %>%
    as_tibble() %>%
    # add taxon_sampling and marker
    mutate(taxon_sampling = taxon_sampling, marker = marker) %>%
    # extract species name and match type
    mutate(
      species_1 = voucher_to_species(voucher_1),
      species_2 = voucher_to_species(voucher_2),
      comp_type = case_when(
        voucher_1 == voucher_2 ~ "self",
        species_1 == species_2 ~ "intra",
        species_1 != species_2 ~ "inter"
      )
    ) %>%
    select(
      voucher_1, voucher_2, species_1, species_2,
      dist, comp_type, marker, taxon_sampling) %>%
    assert(not_na, everything())
}

test_barcode_gap <- function(barcode_dist) {

  # skip calculation if less than two rows (singleton family)
  if (nrow(barcode_dist) < 2) {
    return(tibble())
  }

  # skip calculation if lack both "intra" and "inter" comparisons
  # (need to compare inter vs. intra)
  if (n_distinct(barcode_dist$comp_type) < 2) {
    return(tibble())
  }

  min_inter_dist <-
    barcode_dist %>%
    filter(comp_type == "inter") %>%
    select(species_1, species_2, dist, taxon_sampling, marker) %>%
    pivot_longer(
      names_to = "side", values_to = "species", contains("species")) %>%
    group_by(species, taxon_sampling, marker) %>%
    slice_min(order_by = dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(min_inter_dist = dist, taxon_sampling, marker, species)

  max_intra_dist <-
    barcode_dist %>%
    filter(comp_type == "intra") %>%
    select(species_1, species_2, dist, taxon_sampling, marker) %>%
    pivot_longer(
      names_to = "side", values_to = "species", contains("species")) %>%
    group_by(species, taxon_sampling, marker) %>%
    slice_max(order_by = dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(max_intra_dist = dist, taxon_sampling, marker, species)

  inner_join(
    min_inter_dist,
    max_intra_dist,
    by = c("species", "taxon_sampling", "marker"),
    relationship = "one-to-one") %>%
  assert(not_na, everything()) %>%
  mutate(fail = min_inter_dist < max_intra_dist) %>%
  select(marker, taxon_sampling, species, min_inter_dist, max_intra_dist, fail)

}

# Monophyly test ----

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
    if (is.null(aln_path)) {
      tree_path <- fs::path(wd, deparse(substitute(alignment))) %>%
        fs::path_ext_set(".phy.treefile")
    } else {
      tree_path <- fs::path(fs::path_abs(aln_path)) %>%
        fs::path_ext_set(".phy.treefile")
    }
  }

  # Return single tree if only one file in tree_path
  if (length(tree_path) == 1) {
    assertthat::assert_that(assertthat::is.readable(tree_path))
    res <- ape::read.tree(tree_path)
  }

  # Return list of trees if multiple files in tree_path
  if (length(tree_path) > 1) {
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

#' Check monophyly of species in a tree
#'#'
#' @param tree Input phylogeny; list of class phylo.
#' @param outgroup Name of sinlge tip to use as outgroup.
#' @param sp_delim Character used to separate species name from voucher data
#' in tree tip labels
#' @param aln_file Name of alignment file analyzed
#'
#' @return Dataframe (tibble) with monophyly status for each species
check_monophy <- function(
  tree, aln_file) {

  sp_delim <- "__"
  marker <- str_extract(aln_file, "rbcl_trnlf|rbcl_short|trnlf|rbcl")
  taxon_sampling <- str_extract(aln_file, "all|no_hybrid")

  assertthat::assert_that(assertthat::is.string(sp_delim))
  assertthat::assert_that(assertthat::is.string(marker))
  assertthat::assert_that(assertthat::is.string(taxon_sampling))

  # Get list of outgroup taxa
  og_taxa <- get_og_taxa()

  # Root on Zygnema for rbcL, Equisetum otherwise
  outgroup <- if_else(
    str_detect(marker, "rbcl"),
    tree$tip.label[str_detect(tree$tip.label, "Zygnema")][1],
    tree$tip.label[str_detect(tree$tip.label, "Equisetum")][1]
    )

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
    mutate(taxon_sampling = taxon_sampling) %>%
    mutate(marker = marker)

}

# Blast test ----

#' Combine rbcL and trnL-F sequences
#'
#' Any gaps are removed before combining (pasting) sequences
#'
#' @param seqs Tibble with DNA sequences as a list-column called "seq".
#' @param taxon_sampling_select Taxon sampling to use; must be "all" or 'no_hybrid'.
#'
#' @return Tibble with DNA sequences as a list-column called "seq"; the seq
#' column includes the rbcL and trnLF sequences pasted together into a single
#' sequence.
combine_rbcl_trnlf_seqs <- function(seqs, taxon_sampling_select) {

  rbcl <-
    seqs %>%
    filter(marker == "rbcl", taxon_sampling == taxon_sampling_select) %>%
    mutate(seq = map(seq, ape::del.gaps)) %>%
    select(rbcl_seq = seq, voucher)

  trnlf <-
    seqs %>%
    filter(marker == "trnlf", taxon_sampling == taxon_sampling_select) %>%
    mutate(seq = map(seq, ape::del.gaps)) %>%
    select(trnlf_seq = seq, voucher)

  inner_join(
      rbcl, trnlf, by = "voucher", relationship = "one-to-one"
    ) %>%
    rowwise() %>%
    mutate(seq = paste_seqs(rbcl_seq, trnlf_seq)) %>%
    select(seq, voucher) %>%
    ungroup() %>%
    mutate(marker = "rbcl_trnlf", taxon_sampling = taxon_sampling_select)
}

write_seqtbl_to_blast_fasta <- function(seqs_for_blast_test, dir) {
  marker <- seqs_for_blast_test %>% pull(marker) %>% unique()
  taxon_sampling <- seqs_for_blast_test %>% pull(taxon_sampling) %>% unique()
  assertthat::assert_that(assertthat::is.string(marker))
  assertthat::assert_that(assertthat::is.string(taxon_sampling))

  file <- glue::glue("{marker}_{taxon_sampling}") %>%
    fs::path_ext_set(".fasta") %>%
    fs::path(dir, .)

  seqs_for_blast_test %>%
    seqtbl_to_dnabin(name_col = "voucher") %>%
    ape::del.gaps() %>%
    ape::write.FASTA(file = file)

  file

}

format_blast_db_name <- function(file_name) {
  marker <- str_extract(file_name, "rbcl_trnlf|rbcl_short|trnlf|rbcl")
  taxon_sampling <- str_extract(file_name, "all|no_hybrid")
  assertthat::assert_that(assertthat::is.string(marker))
  assertthat::assert_that(assertthat::is.string(taxon_sampling))
  paste(marker, taxon_sampling, sep = "_")
}

format_blast_output_name <- function(file_name) {
  marker <- str_extract(file_name, "rbcl_trnlf|rbcl_short|trnlf|rbcl")
  taxon_sampling <- str_extract(file_name, "all|no_hybrid")
  assertthat::assert_that(assertthat::is.string(marker))
  assertthat::assert_that(assertthat::is.string(taxon_sampling))
  paste(marker, taxon_sampling, "blastn_results.tsv", sep = "_")
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
blast_n <- function(
  query,
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

load_blast_tsv <- function(blast_out_path) {
  # Read in sequences.
  # BLAST doesn't output column headers, so we need to specify
  # (make sure they match correctly first!)
  fmt6_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  marker <- str_extract(blast_out_path, "rbcl_trnlf|rbcl_short|trnlf|rbcl")
  taxon_sampling <- str_extract(blast_out_path, "all|no_hybrid")
  assertthat::assert_that(assertthat::is.string(marker))
  assertthat::assert_that(assertthat::is.string(taxon_sampling))

  # Read in BLAST output
  readr::read_tsv(
    blast_out_path,
    col_names = fmt6_cols,
    col_types = "ccdddddddddd" # two ID cols are char, rest is numeric
  ) %>%
  mutate(taxon_sampling = taxon_sampling) %>%
  mutate(marker = marker)
}

prep_blast_res <- function(blast_res_raw, seqs_for_blast_test) {

  q_seqs <- seqtbl_to_dnabin(seqs_for_blast_test, name_col = "voucher")

  # Get lengths (bp) of query seqs
  q_seqs <- ape::del.gaps(q_seqs)
  seq_lengths <-
    tibble(
      qseqid = names(q_seqs),
      q_seq_len = map_dbl(q_seqs, length)
    )

  blast_res_raw %>%
    # Add query seq lengths
    left_join(seq_lengths, by = "qseqid", relationship = "many-to-one") %>%
    assert(not_na, q_seq_len) %>%
    # Filter to only those matching 95% of original sequence length
    filter(length >= 0.95 * q_seq_len) %>%
    # Parse out query species and subject species
    # while ignoring infraspecific taxa
    mutate(
      # '__' separates taxon and voucher
      q_taxon = str_split_i(qseqid, "__", 1),
      s_taxon = str_split_i(sseqid, "__", 1),
      # '_' separates genus, specific epithet, and infraspecific epithet
      # eg: "Asplenium_wilfordii_var._densum"
      q_genus = str_split_i(q_taxon, "_", 1),
      s_genus = str_split_i(s_taxon, "_", 1),
      q_specific_epithet = str_split_i(q_taxon, "_", 2),
      s_specific_epithet = str_split_i(s_taxon, "_", 2),
      q_species = paste(q_genus, q_specific_epithet),
      s_species = paste(s_genus, s_specific_epithet)
    ) %>%
    group_by(q_species) %>%
    mutate(n_indiv = n_distinct(qseqid)) %>%
    ungroup() %>%
    # Classify comparison types
    mutate(
      comp_type = case_when(
        qseqid == sseqid ~ "self",
        q_species == s_species ~ "intra",
        q_species != s_species ~ "inter"
      )
    )
}

get_species_cutoff <- function(blast_res) {
  blast_res %>%
    filter(comp_type == "intra") %>%
    summarize(
      mean_pident = mean(pident),
      min_pident = min(pident),
      marker = unique(marker),
      taxon_sampling = unique(taxon_sampling),
      .groups = "drop") %>%
    pivot_longer(
      names_to = "cutoff_type",
      -c(taxon_sampling, marker)
    )
}


#' Choose a species-level cutoff value for BLAST barcode test
#'
#' Also groups the dataframe for looping in the targets pipeline
#'
#' @param cutoff_table Dataframe; tibble with potential cutoff values to use
#' (including mean and minimum interspecific identities).
#'
#' @return Dataframe with single cutoff value selected for each combination of
#' marker and taxon sampling
select_and_group_cutoff <- function(cutoff_table) {
  cutoff_select <-
    cutoff_table %>%
    filter(taxon_sampling == "no_hybrid", cutoff_type == "mean_pident") %>%
    select(marker, taxon_sampling, cutoff_use = value)

  # Use the mean pident from the no-hybrids taxon sampling for each marker
  cutoff_select %>%
    bind_rows(
      mutate(cutoff_select, taxon_sampling = "all")
    ) %>%
    group_by(marker, taxon_sampling) %>%
    targets::tar_group()
}

test_blast <- function(
  blast_res,
  cutoff_table) {

  # Obtain value to use for infraspecific cutoff
  intra_cutoff <-
    cutoff_table %>%
    pull(cutoff_use)

  assertthat::assert_that(assertthat::is.number(intra_cutoff))

  blast_res %>%
    mutate(
      above_cutoff = pident > intra_cutoff,
      match_other = q_species != s_species,
      fail = above_cutoff & match_other) %>%
    group_by(q_species) %>%
    summarize(
      fail = any(fail),
      marker = unique(marker),
      taxon_sampling = unique(taxon_sampling),
      .groups = "drop")
}

summarize_blast_fail_rate <- function(blast_test_res_raw) {
  blast_test_res_raw %>%
    group_by(marker, taxon_sampling) %>%
    count(fail) %>%
    ungroup()
}

# I/O ----

read_fasta_to_tbl <- function(
  path, ...) {
  ape::read.FASTA(path) %>%
  dnabin_to_seqtbl(., name_col = "voucher", seq_col = "seq") %>%
  # drop periods in sequence names
  mutate(voucher = str_remove_all(voucher, "\\.")) %>%
  # add species column
  mutate(species = voucher_to_species(voucher)) %>%
  mutate(...)
}

read_txt_to_tbl <- function(path, col_name = "voucher") {
  x <- readr::read_lines(path)
  tibble::tibble("{col_name}" := x)
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

write_fasta <- function(x, file, ...) {
  ape::write.FASTA(x = x, file = file, ...)
  file
}

write_seqtbl_to_phy <- function(
  seqtbl, file, name_col = "accession", seq_col = "seq") {
  seqtbl %>%
    seqtbl_to_dnabin(name_col = name_col, seq_col = seq_col) %>%
    write_phy(file = file)
}

write_seqtbl_aln_to_phy <- function(seqtbl, dir) {

  marker <- unique(seqtbl$marker)
  taxon_sampling <- unique(seqtbl$taxon_sampling)
  assertthat::assert_that(assertthat::is.string(marker))
  assertthat::assert_that(assertthat::is.string(taxon_sampling))

  seq_lens <-
    seqtbl %>%
    mutate(seq_len = length(seq)) %>%
    pull(seq_len) %>%
    unique()

  assertthat::assert_that(
    isTRUE(length(seq_lens) == 1),
    msg = "Aligned seqs not all the same length.")

  file <- glue::glue("{marker}_{taxon_sampling}") %>%
    fs::path_ext_set(".phy") %>%
    fs::path(dir, .)

  seqtbl %>%
    seqtbl_to_dnabin(name_col = "voucher") %>%
    as.matrix() %>%
    write_phy(file = file)
}

# Utils ----

get_og_taxa <- function() {
  # get outgroup sequences used by ftol
  ftolr::accessions_wide %>%
    filter(outgroup) %>%
    pull(species)
}

subset_seqs_to_family <- function(seqs) {
  seqs %>%
    add_count(family, taxon_sampling) %>%
    filter(n > 1) %>%
    group_by(family, taxon_sampling) %>%
    mutate(n_species = n_distinct(species)) %>%
    filter(n_species > 1) %>%
    ungroup()
}

voucher_to_species <- function(voucher, sep = " ") {
  # '__' separates taxon and voucher
  taxon <- str_split_i(voucher, "__", 1)
  # '_' separates genus, specific epithet, and infraspecific epithet
  # eg: "Asplenium_wilfordii_var._densum"
  genus <- str_split_i(voucher, "_", 1)
  specific_epithet <- str_split_i(voucher, "_", 2)
  paste0(genus, sep, specific_epithet)
}

#' Convert DNA sequences from list of class DNAbin to tibble
#'
#' @param dnabin List or matrix of class DNAbin
#' @param name_col Name of column to use for sequence name
#' @param seq_col Name of column to use for sequences (list-column)
#' @return Tibble with one row per sequence
#'
dnabin_to_seqtbl <- function(dnabin, name_col = "accession", seq_col = "seq") {
  # Convert to a list in case DNA sequences are aligned (in matrix)
  dnabin <- as.list(dnabin)

  # Check input
  assertthat::assert_that(inherits(dnabin, "DNAbin"))
  assertthat::assert_that(assertthat::is.string(name_col))
  assertthat::assert_that(assertthat::is.string(seq_col))

  tibble::tibble(
    "{seq_col}" := split(dnabin, 1:length(dnabin)),
    "{name_col}" := names(dnabin))
}

#' Convert a list-column of DNA sequences in a tibble to a list of class DNAbin
#'
#' @param seqtbl Tibble containing one DNA sequence per row
#' in a list-column
#' @param name_col Name of column with sequence name
#' @param seq_col Name of column with sequences (list-column)
#' @return List of class DNAbin
#'
seqtbl_to_dnabin <- function(seqtbl, name_col = "accession", seq_col = "seq") {
  require(ape)

  # Check input
  assertthat::assert_that(inherits(seqtbl, "tbl"))
  assertthat::assert_that(assertthat::is.string(name_col))
  assertthat::assert_that(name_col %in% colnames(seqtbl))

  # Extract sequences from metadata and rename
  seqs_dnabin <- do.call(c, seqtbl[[seq_col]])
  names(seqs_dnabin) <- seqtbl[[name_col]]

  # Make sure that went OK
  assertthat::assert_that(is.list(seqs_dnabin))
  assertthat::assert_that(inherits(seqs_dnabin, "DNAbin"))
  assertthat::assert_that(all(names(seqs_dnabin) == seqtbl[[name_col]]))

  seqs_dnabin
}

#' Paste two sequences end to end
#'
#' x and y must have the same sequence name
#'
#' @param x First sequence; list of class DNAbin of length 1.
#' @param y Second sequence; list of class DNAbin of length 1.
#'
#' @return List of class DNAbin of length 1: x and y pasted together end to end
#' @examples 
#' paste_seqs(
#'   set_names(as.list(woodmouse[1, ]), "a"),
#'   set_names(as.list(woodmouse[2, ]), "a")
#' )
paste_seqs <- function(x, y) {
  x_char <- as.character(x)[[1]] %>% paste(collapse = "")
  y_char <- as.character(y)[[1]] %>% paste(collapse = "")
  assertthat::assert_that(
    isTRUE(names(x) == names(y))
  )
  seq_name <- names(x)
  
  combined_char <- paste0(x_char, y_char)
  res <- combined_char %>%
    strsplit("") %>%
    ape::as.DNAbin() %>%
    set_names(seq_name)
  list(seq = res)
}
