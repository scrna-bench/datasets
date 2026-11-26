#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# 1.3m: 10x 1.3 million mouse brain cells (TENxBrainData).
# Truth: library index, only so the full pipeline runs; ARIs on this dataset
# have no biological meaning.

suppressPackageStartupMessages({
  library(TENxBrainData)
})

out <- dataset_args()

sce <- with_phase("download", TENxBrainData())

sce <- with_phase("build", {
  rownames(sce) <- rowData(sce)$Symbol
  sce
})

write_dataset(
  sce,
  truth = colData(sce)$Library,
  qc = make_qc_df(
    nFeature_min = 200, nFeature_max = 2500,
    nCount_max = 4000,
    percent_mt_max = 5
  ),
  out
)
