source("R/packages.R")
source("R/functions.R")


tar_plan(
  # Load sequences ---
  tar_file_read(
    rbcl_seqs_all,
    "data_raw/BIFA_rbcL_analyses_240313.fasta",
    ape::read.FASTA(!!.x),
  ),
  tar_file_read(
    trnlf_seqs_all,
    "data_raw/BIFA_trnLF_analyses_240315.fasta",
    ape::read.FASTA(!!.x),
  ),
  # Drop specimens not identified to species
  rbcl_seqs = drop_indets(rbcl_seqs_all),
  trnlf_seqs = drop_indets(trnlf_seqs_all),
  # Align ingroup only ----
  rbcl_fern_align = align_seqs(rbcl_seqs),
  trnlf_fern_align = align_seqs(trnlf_seqs),
  # Analyze genetic distances ----
  rbcl_dist = analyze_dist(rbcl_fern_align),
  trnlf_dist = analyze_dist(trnlf_fern_align),
  # Add outgroups and align ----
  # only do rbcL, since trnLF sequences are too diverged to align
  rbcl_with_og = align_with_og(rbcl_seqs, marker = "rbcL"),
  # Build tree ----
  rbcl_tree = iqtree(
    alignment = rbcl_with_og,
    wd = "_targets/user/iqtree",
    m = "MFP", # test model followed by ML analysis
    bb = 1000,
    seed = 20240301,
    redo = TRUE,
    other_args = c(
      "-mset", "GTR", # only test GTR models
      "-mrate", "E,I,G,I+G", # don't test free-rate models
      "-t", "PARS",
      "-nt", "AUTO"
    )
  ),
  rbcl_monophyly = check_monophy(rbcl_tree),
  # Write report ----
  tarchetypes::tar_quarto(
    report
  )
)