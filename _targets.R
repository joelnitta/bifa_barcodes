source("R/packages.R")
source("R/functions.R")


tar_plan(
  # Load data ---
  tar_file(
    rbcl_file,
    "data_raw/BIFA_rbcL_analyses_240313.fasta"
  ),
  tar_file(
    trnlf_file,
    "data_raw/BIFA_trnLF_analyses_240315.fasta"
  ),
  rbcl = ape::read.FASTA(rbcl_file),
  rbcl_ss = Biostrings::readDNAStringSet(rbcl_file),
  trnlf = ape::read.FASTA(trnlf_file),
  trnlf_ss = Biostrings::readDNAStringSet(trnlf_file),
  # Analyze genetic distances ----
  rbcl_dist = analyze_dist(rbcl_ss),
  trnlf_dist = analyze_dist(trnlf_ss),
  # Add outgroups and align ----
  rbcl_with_og = align_with_og(rbcl, marker = "rbcL"),
  tar_file(
    rbcl_align_file,
    write_phy(rbcl_with_og, "_targets/user/iqtree/rbcl.phy")
  ),
  rbcl_tree = iqtree(
    alignment = rbcl_align_file,
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
  # rbcl_monophyly = check_monophy(rbcl_tree),
  # trnlf_monophyly = check_monophy(trnlf_tree),
  tarchetypes::tar_quarto(
    report
  )
)