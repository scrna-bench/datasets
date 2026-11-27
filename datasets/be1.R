#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# be1: GEO GSE243665, one 10x run per sample, concatenated.
# Truth: sample.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(GEOquery)
})

out <- dataset_args()

raw_dir <- with_phase("download", {
  gse_id <- "GSE243665"
  getGEOSuppFiles(GEO = gse_id, baseDir = tempdir())
  file.path(tempdir(), gse_id)
})

sce <- with_phase("build", {
  # one 10x-style folder per sample, with the GEO/sample prefix stripped from
  # the file names so read10xCounts finds matrix/features/barcodes
  files <- list.files(raw_dir, pattern = "\\.(mtx|tsv)\\.gz$", full.names = TRUE)
  samples <- unique(regmatches(basename(files), regexpr("_[^_]+_", basename(files))))
  samples <- gsub("_", "", samples)
  for (s in samples) {
    sample_dir <- file.path(raw_dir, s)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    for (f in files[grepl(paste0("_", s, "_"), basename(files))]) {
      file.rename(f, file.path(sample_dir, sub("^[^_]+_[^_]+_", "", basename(f))))
    }
  }

  sce <- do.call(cbind, lapply(samples, function(s) {
    read10xCounts(samples = file.path(raw_dir, s), sample.names = s, col.names = TRUE)
  }))
  unlink(raw_dir, recursive = TRUE)
  rownames(sce) <- make.unique(rowData(sce)$Symbol)
  sce
})

write_dataset(
  sce,
  truth = colData(sce)$Sample,
  qc = make_qc_df(
    nFeature_min = 200, nFeature_max = 5000,
    nCount_max = 25000,
    percent_mt_max = 5
  ),
  out
)
