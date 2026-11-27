#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# cb: cord blood CITE-seq (SingleCellMultiModal), human genes only.
# Truth: provided cell type annotation.

suppressPackageStartupMessages({
  library(SingleCellMultiModal)
})

out <- dataset_args()

sce <- with_phase("download", {
  CITEseq(
    DataType = "cord_blood", modes = "*", dry.run = FALSE, version = "1.0.0",
    DataClass = "SingleCellExperiment"
  )
})

sce <- with_phase("build", {
  # keep human genes and drop the "HUMAN_" prefix
  sce <- sce[grep("^HUMAN_", rownames(sce)), ]
  rownames(sce) <- sub("^HUMAN_", "", rownames(sce))
  sce
})

write_dataset(
  sce,
  truth = colData(sce)$celltype,
  qc = make_qc_df(
    nFeature_min = 200, nFeature_max = 2500,
    nCount_max = 4000,
    percent_mt_max = 5
  ),
  out
)
