#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# sc-mix: scMixology, five cell lines on 10x (Tian et al. 2019).
# Truth: cell line.

out <- dataset_args()

sce <- with_phase("download", {
  url <- "https://github.com/LuyiTian/sc_mixology/raw/refs/heads/master/data/sincell_with_class_5cl.RData"
  raw_path <- file.path(tempdir(), basename(url))
  options(timeout = max(600, getOption("timeout")))
  download.file(url, destfile = raw_path)
  env <- new.env()
  load(raw_path, envir = env)
  unlink(raw_path)
  env$sce_sc_10x_5cl_qc
})

write_dataset(
  sce,
  truth = colData(sce)$cell_line,
  qc = make_qc_df(
    nFeature_min = 200, nFeature_max = 6200,
    nCount_max = 60000,
    percent_mt_max = 10
  ),
  out
)
