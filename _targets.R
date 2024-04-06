source("R/packages.R")
source("R/functions.R")

tar_plan(
  # Load data ---
  # - rbcL sequences
  tar_file_read(
    rbcl_seqs_all,
    "data_raw/BIFA_rbcL_analyses_240404.fasta",
    ape::read.FASTA(!!.x),
  ),
  # - trnLF sequences
  tar_file_read(
    trnlf_seqs_all,
    "data_raw/BIFA_trnLF_analyses_240404.fasta",
    ape::read.FASTA(!!.x),
  ),
  # - hybrid taxa
  tar_file_read(
    hybrid_taxa,
    "data_raw/hybrid_taxa.txt",
    read_lines(!!.x)
  ),
  # - species complexes
  tar_file_read(
    sp_complex,
    "data_raw/species_complex.txt",
    read_lines(!!.x)
  ),

  # Filter data ---
  # - Drop specimens not identified to species
  rbcl_seqs = drop_indets(rbcl_seqs_all),
  trnlf_seqs = drop_indets(trnlf_seqs_all),
  # - Also drop species complexes and hybrids
  rbcl_seqs_no_hybrids = drop_complex_hybrids(
    rbcl_seqs, sp_complex, hybrid_taxa),
  trnlf_seqs_no_hybrids = drop_complex_hybrids(
    trnlf_seqs, sp_complex, hybrid_taxa),

  # Monophyly test ----
  # Only do rbcL, since trnLF sequences are too diverged to align
  # - Prep lists for looping
  seqs_for_phy_analysis = list(rbcl_seqs, rbcl_seqs_no_hybrids),
  seqs_phy_names = c("rbcl", "rbcl_no_hybrids"),
  # - Align sequences with outgroup and write to file
  tar_target(
    seqs_for_phy_aligned,
    align_with_og(seqs_for_phy_analysis[[1]], marker = "rbcL"),
    pattern = map(seqs_for_phy_analysis)
  ),
  tar_target(
    alignment_files,
    write_phy(
      seqs_for_phy_aligned,
      paste0("_targets/user/iqtree/", seqs_phy_names, ".phy")
    ),
    pattern = map(seqs_for_phy_aligned, seqs_phy_names)
  ),
  # - Build tree
  tar_target(
    barcode_tree,
    iqtree(
      aln_path = alignment_files[[2]],
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
    pattern = map(alignment_files)
  ),
  # - Analyze monophyly
  tar_target(
    monophyly_res,
    check_monophy(barcode_tree, dataset = seqs_phy_names),
    pattern = map(barcode_tree, seqs_phy_names)
  ),

  # BLAST test ----
  # - create lists for looping
  seqs_list = list(
    rbcl_seqs, trnlf_seqs, rbcl_seqs_no_hybrids, trnlf_seqs_no_hybrids),
  seqs_names = c(
    "rbcl", "trnlf", "rbcl_no_hybrids", "trnlf_no_hybrids"
  ),
  # - Write out sequences for BLAST db
  tar_target(
    seqs_for_blast_db,
    write_fasta(
      ape::del.gaps(seqs_list[[1]]),
      paste0("_targets/user/blast/", seqs_names, ".fasta")
    ),
    format = "file",
    pattern = map(seqs_list, seqs_names)
  ),
  # - Create BLAST databases
  tar_target(
    blast_db,
    make_blast_db(
      seqs_path = seqs_for_blast_db,
      blast_db_dir = "_targets/user/blast/",
      out_name = seqs_names),
    format = "file",
    pattern = map(seqs_for_blast_db, seqs_names)
  ),
  # - Run BLAST query
  tar_target(
    blast_res_tsv,
    blast_n(
      query = seqs_for_blast_db,
      database = seqs_names,
      out_file = paste0(seqs_names, "_blastn_results.tsv"),
      wd = "_targets/user/blast/",
      depends = blast_db
    ),
    format = "file",
    pattern = map(seqs_for_blast_db, seqs_names, blast_db)
  ),
  tar_target(
    blast_res_raw,
    load_blast_tsv(blast_res_tsv, seqs_names),
    pattern = map(blast_res_tsv, seqs_names)
  ),
  # Filter out matches at < 95% of sequence length,
  # flag intra- vs. inter-specific matches
  tar_target(
    blast_res,
    prep_blast_res(blast_res_raw, seqs_list[[1]]),
    pattern = map(blast_res_raw, seqs_list)
  ),
  tar_target(
    cutoff_table,
    get_species_cutoff(blast_res, seqs_names),
    pattern = map(blast_res, seqs_names)
  ),
  tar_target(
    blast_test_res_raw,
    test_blast(
      blast_res,
      dataset_select = seqs_names,
      cutoff_table
    ),
    pattern = map(seqs_names)
  ),
  blast_fail_rate_summary = summarize_blast_fail_rate(blast_test_res_raw),
  # Write report ----
  tarchetypes::tar_quarto(
    report
  )
)