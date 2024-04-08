source("R/packages.R")
source("R/functions.R")

tar_option_set(
  workspace_on_error = TRUE
)

tar_plan(
  # Load data ---
  # - rbcL sequences
  tar_file_read(
    rbcl_seqs_all,
    "data_raw/BIFA_rbcL_analyses_240407.fasta",
    ape::read.FASTA(!!.x),
  ),
  # - rbcL sequences (short)
  tar_file_read(
    rbcl_seqs_all_short,
    "data_raw/BIFA_rbcL_analyses_240407_short.fasta",
    ape::read.FASTA(!!.x),
  ),
  # - trnLF sequences
  tar_file_read(
    trnlf_seqs_all,
    "data_raw/BIFA_trnLF_analyses_240408.fasta",
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
  # - PPG I classification system (modified)
  tar_file_read(
    ppgi,
    "data_raw/ppgi_taxonomy_mod.csv",
    read_csv(!!.x)
  ),

  # Filter data ---
  # - Drop specimens not identified to species
  rbcl_seqs = drop_indets(rbcl_seqs_all),
  rbcl_seqs_short = drop_indets(rbcl_seqs_all_short),
  trnlf_seqs = drop_indets(trnlf_seqs_all),
  # - Also drop species complexes and hybrids
  rbcl_seqs_no_hybrids = drop_complex_hybrids(
    rbcl_seqs, sp_complex, hybrid_taxa),
  rbcl_seqs_no_hybrids_short = drop_complex_hybrids(
    rbcl_seqs_short, sp_complex, hybrid_taxa),
  trnlf_seqs_no_hybrids = drop_complex_hybrids(
    trnlf_seqs, sp_complex, hybrid_taxa),
  # - create lists for looping
  seqs_list = list(
    rbcl_seqs, rbcl_seqs_no_hybrids,
    trnlf_seqs, trnlf_seqs_no_hybrids,
    rbcl_seqs_short, rbcl_seqs_no_hybrids_short),
  seqs_names = c(
    "rbcl", "rbcl_no_hybrids",
    "trnlf", "trnlf_no_hybrids",
    "rbcl_short", "rbcl_no_hybrids_short"
  ),

  # Barcoding gap ---
  # - Subset each set of sequences to families with >1 species per family
  tar_target(
    seqs_by_family,
    subset_seqs_to_family(
      seqs_list[[1]],
      seqs_names,
      ppgi
    ),
    pattern = map(seqs_list, seqs_names)
  ),
  # - Align seqs within each family
  seqs_by_family_for_align = seqs_by_family, # need this to map by family
  tar_target(
    seqs_by_family_aligned,
    align_seqs(seqs_by_family_for_align[[1]]),
    pattern = map(seqs_by_family_for_align),
    iteration = "list"
  ),
  # - Calculate barcode marker distances within each family
  tar_target(
    barcode_dist,
    calc_barcode_dist(seqs_by_family_aligned),
    pattern = map(seqs_by_family_aligned)
  ),

  # Monophyly test ----
  # Only do rbcL, since trnLF sequences are too diverged to align
  # - Prep lists for looping
  seqs_for_phy_analysis = list(
    rbcl_seqs, rbcl_seqs_no_hybrids,
    rbcl_seqs_short, rbcl_seqs_no_hybrids_short
    ),
  seqs_phy_names = c(
    "rbcl", "rbcl_no_hybrids",
    "rbcl_short", "rbcl_no_hybrids_short"),
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
      aln_path = alignment_files,
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
    pattern = map(alignment_files),
    iteration = "list"
  ),
  # - Analyze monophyly
  tar_target(
    monophyly_res,
    check_monophy(barcode_tree, dataset = seqs_phy_names),
    pattern = map(barcode_tree, seqs_phy_names)
  ),

  # BLAST test ----
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
