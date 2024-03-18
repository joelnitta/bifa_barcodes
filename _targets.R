source("R/packages.R")
source("R/functions.R")


tar_plan(
  # Load data ---
  tar_file_read(
    rbcl,
    "data_raw/BIFA_rbcL_analyses_240313.fasta",
    ape::read.FASTA(!!.x),
  ),
  tar_file_read(
    trnlf,
    "data_raw/BIFA_trnLF_analyses_240315.fasta",
    ape::read.FASTA(!!.x),
  ),
  # Align ingroup only ----
  trnlf_fern_align = align_seqs(trnlf),
  rbcl_fern_align = align_seqs(rbcl),
  # Analyze genetic distances ----
  rbcl_dist = analyze_dist(rbcl_fern_align),
  trnlf_dist = analyze_dist(trnlf_fern_align),
  # Add outgroups and align ----
  rbcl_with_og = align_with_og(rbcl, marker = "rbcL"),
  tar_file(
    rbcl_align_file,
    write_phy(rbcl_with_og, "_targets/user/iqtree/rbcl.phy")
  ),
  # TODO: can't get outgroup seqs for trnLF
  trnlf_with_og = align_with_og(trnlf, marker = "trnL-trnF"),
  tar_file(
    trnlf_align_file,
    write_phy(rbcl_with_og, "_targets/user/iqtree/trnlf.phy")
  ),
  # Build tree ----
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
  trnlf_tree = iqtree(
    alignment = trnlf_align_file,
    m = "MFP",
    bb = 1000,
    seed = 20240301,
    redo = TRUE,
    other_args = c(
      "-mset", "GTR",
      "-mrate", "E,I,G,I+G",
      "-t", "PARS",
      "-nt", "AUTO"
    )
  ),
  # rbcl_monophyly = check_monophy(rbcl_tree),
  # trnlf_monophyly = check_monophy(trnlf_tree),
  # Write report ----
  tarchetypes::tar_quarto(
    report
  )
)