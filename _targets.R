source("R/packages.R")

tar_plan(
  tar_file_read(
    rbcl,
    "data_raw/BIFA_rbcL_analyses_240313.fasta",
    ape::read.FASTA(!!.x)
  ),
  tarchetypes::tar_quarto(
    report
  )
)