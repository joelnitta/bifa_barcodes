source("R/packages.R")
source("R/functions.R")

tar_option_set(
  workspace_on_error = TRUE,
  controller = crew_controller_local(workers = 8)
)

tar_plan(
  # Load data ---
  # DNA sequences
  tar_file(
    fern_seq_paths,
    c(
      rbcl = "data_raw/BIFA_rbcL_analyses_240407.fasta",
      rbcl_short = "data_raw/BIFA_rbcL_analyses_240407_barcoding_short.fasta",
      trnlf = "data_raw/BIFA_trnLF_analyses_240408.fasta"
    )
  ),
  fern_seq_paths_map = fern_seq_paths, # work-around to map over fern_seq_paths
  markers = c("rbcl", "rbcl_short", "trnlf"),
  tar_target(
    fern_seqs_raw,
    read_fasta_to_tbl(fern_seq_paths_map, marker = markers),
    pattern = map(fern_seq_paths_map, markers)
  ),
  outgroup_seqs = load_rbcl_og_seqs(),
  seqs_raw = bind_rows(
    add_family(fern_seqs_raw, ppgi),
    outgroup_seqs),
  # - hybrid taxa
  tar_file_read(
    hybrid_taxa,
    "data_raw/hybrid_taxa.txt",
    read_txt_to_tbl(!!.x, col_name = "voucher")
  ),
  # - species complexes
  tar_file_read(
    sp_complex,
    "data_raw/species_complex.txt",
    read_txt_to_tbl(!!.x, col_name = "species")
  ),
  # - PPG I classification system (modified)
  tar_file_read(
    ppgi,
    "data_raw/ppgi_taxonomy_mod.csv",
    read_csv(!!.x)
  ),
  # - Drop specimens not identified to species
  seqs_all_with_indets = seqs_raw %>%
    mutate(taxon_sampling = "all_with_indet"),
  seqs_all = drop_indets(seqs_raw) %>%
    mutate(taxon_sampling = "all"),
  # - Also drop species complexes and hybrids
  seqs_no_hybrid = drop_complex_hybrids(
    seqs_all,
    hybrid_taxa,
    sp_complex
  ) %>% mutate(taxon_sampling = "no_hybrid"),
  # - Final unaligned sequences
  seqs = bind_rows(seqs_all_with_indets, seqs_all, seqs_no_hybrid) %>%
    assert(not_na, everything()),

  # DNA alignment ---
  # - trnlf: subset each set of sequences to families with >1 species per family
  tar_target(
    trnlf_to_align,
    subset_seqs_to_family(filter(seqs, marker == "trnlf")) %>%
      group_by(family, taxon_sampling) %>%
      tar_group,
    iteration = "group"
  ),
  # trnlf: align within family and taxon_sampling
  tar_target(
    trnlf_aligned_family,
    align_seqs_tbl(trnlf_to_align, name_col = "voucher"),
    pattern = map(trnlf_to_align),
    iteration = "vector"
  ),
  # rbcl: align without outgroups for gap test
  tar_target(
    rbcl_to_align_for_gap_test,
    filter(seqs, family != "outgroup", str_detect(marker, "rbcl")) %>%
      subset_seqs_to_family() %>%
      group_by(marker, taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    rbcl_aligned_for_gap_test,
    align_seqs_tbl(rbcl_to_align_for_gap_test, name_col = "voucher"),
    pattern = map(rbcl_to_align_for_gap_test),
    iteration = "vector"
  ),
  # rbcl: align with outgroup for monophyly test
  tar_target(
    rbcl_to_align_for_monophy_test,
    filter(seqs, str_detect(marker, "rbcl")) %>%
      group_by(marker, taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    rbcl_aligned_for_monophy_test,
    align_seqs_tbl(rbcl_to_align_for_monophy_test, name_col = "voucher"),
    pattern = map(rbcl_to_align_for_monophy_test),
    iteration = "vector"
  ),
  # trnlf: merge subalignments
  tar_target(
    trnlf_seqs,
    filter(seqs, marker == "trnlf") %>%
      group_by(taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    trnlf_aligned_family_grouped,
    trnlf_aligned_family %>%
      group_by(taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    trnlf_aligned,
    merge_trnlf(trnlf_seqs, trnlf_aligned_family_grouped),
    pattern = map(trnlf_seqs, trnlf_aligned_family_grouped)
  ),
  # concatenate rbcl and trnlf: with outgroups (for monophy test)
  tar_target(
    trnlf_aligned_to_cat,
    trnlf_aligned %>%
      group_by(taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    rbcl_aligned_to_cat_for_monophy_test,
    rbcl_aligned_for_monophy_test %>%
      filter(marker == "rbcl") %>%
      group_by(taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    rbcl_trnlf_aligned_for_monophy_test,
    concatenate_trnlf_rbcl(
      rbcl_aligned_to_cat_for_monophy_test,
      trnlf_aligned_to_cat
    ),
    pattern = map(rbcl_aligned_to_cat_for_monophy_test, trnlf_aligned_to_cat)
  ),
  # concatenate rbcl and trnlf: without outgroups (for gap test)
  tar_target(
    rbcl_aligned_to_cat_for_gap_test,
    rbcl_aligned_for_gap_test %>%
      filter(marker == "rbcl") %>%
      group_by(taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    rbcl_trnlf_aligned_for_gap_test,
    concatenate_trnlf_rbcl(
      rbcl_aligned_to_cat_for_gap_test,
      trnlf_aligned_to_cat
    ),
    pattern = map(rbcl_aligned_to_cat_for_gap_test, trnlf_aligned_to_cat)
  ),

  # Barcode gap test ----
  # - Assemble aligned sequences to test
  tar_target(
    seqs_aligned_for_gap_test,
    bind_rows(
      trnlf_aligned,
      rbcl_aligned_for_gap_test,
      rbcl_trnlf_aligned_for_gap_test) %>%
      group_by(marker, taxon_sampling, family) %>%
      tar_group(),
    iteration = "group"
  ),
  # - Calculate barcode marker distances within each family
  tar_target(
    barcode_dist,
    calc_barcode_dist(seqs_aligned_for_gap_test),
    pattern = map(seqs_aligned_for_gap_test)
  ),
  # - Compare maximum intraspecific vs. minimum interspecific distances
  tar_target(
    barcode_gap_res_raw,
    test_barcode_gap(barcode_dist),
    pattern = map(barcode_dist)
  ),

  # Monophyly test ----
  # Select alignments to test
  tar_target(
    alignments_for_monophy_test,
    bind_rows(
      rbcl_aligned_for_monophy_test,
      rbcl_trnlf_aligned_for_monophy_test,
      trnlf_aligned
    ) %>%
    group_by(marker, taxon_sampling) %>%
    tar_group(),
    iteration = "group"
  ),
  # Write alignments to file
  tar_target(
    alignment_files,
    write_seqtbl_aln_to_phy(
      alignments_for_monophy_test,
      dir = "_targets/user/iqtree/"
    ),
    pattern = map(alignments_for_monophy_test)
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
    check_monophy(barcode_tree, alignment_files),
    pattern = map(barcode_tree, alignment_files)
  ),

  # BLAST test ----
  taxon_sampling_names = c("all", "no_hybrid"),
  tar_target(
    rbcl_trnlf_pasted_seqs,
    combine_rbcl_trnlf_seqs(seqs, taxon_sampling_select = taxon_sampling_names),
    pattern = map(taxon_sampling_names)
  ),
  # Select sequences to test
  tar_target(
    seqs_for_blast_test,
    bind_rows(
      rbcl_trnlf_pasted_seqs,
      seqs,
    ) %>%
    group_by(marker, taxon_sampling) %>%
    tar_group(),
    iteration = "group"
  ),
  # Write seqs to file
  tar_target(
    seq_files_for_blast_db,
    write_seqtbl_to_blast_fasta(
      seqs_for_blast_test,
      dir = "_targets/user/blast/"
    ),
    pattern = map(seqs_for_blast_test)
  ),

  # - Create BLAST databases
  tar_target(
    blast_db,
    make_blast_db(
      seqs_path = seq_files_for_blast_db,
      blast_db_dir = "_targets/user/blast/",
      out_name = format_blast_db_name(seq_files_for_blast_db)),
    format = "file",
    pattern = map(seq_files_for_blast_db)
  ),
  # - Run BLAST query
  tar_target(
    blast_res_tsv,
    blast_n(
      query = seq_files_for_blast_db,
      database = format_blast_db_name(seq_files_for_blast_db),
      out_file = format_blast_output_name(seq_files_for_blast_db),
      wd = "_targets/user/blast/",
      depends = blast_db
    ),
    format = "file",
    pattern = map(seq_files_for_blast_db, blast_db)
  ),
  tar_target(
    blast_res_raw,
    load_blast_tsv(blast_res_tsv),
    pattern = map(blast_res_tsv)
  ),
  # Filter out matches at < 95% of sequence length,
  # flag intra- vs. inter-specific matches
  tar_target(
    blast_res_raw_grouped,
    blast_res_raw %>%
      group_by(marker, taxon_sampling) %>%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    blast_res,
    prep_blast_res(blast_res_raw_grouped, seqs_for_blast_test),
    pattern = map(blast_res_raw_grouped, seqs_for_blast_test)
  ),
  # Determine marker-specific species level cutoffs
  tar_target(
    cutoff_table,
    get_species_cutoff(blast_res),
    pattern = map(blast_res)
  ),
  tar_target(
    cutoff_table_grouped,
    select_and_group_cutoff(cutoff_table),
    iteration = "group"
  ),
  tar_target(
    blast_test_res_raw,
    test_blast(
      blast_res,
      cutoff_table_grouped
    ),
    pattern = map(blast_res, cutoff_table_grouped)
  ),
  blast_fail_rate_summary = summarize_blast_fail_rate(blast_test_res_raw),
  # Write report ----
  tar_file(
    references,
    "reports/references.yaml"
  ),
  tarchetypes::tar_quarto(
    report,
    quiet = FALSE
  )
)
