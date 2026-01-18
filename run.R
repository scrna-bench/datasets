#!/usr/bin/env Rscript

library(argparse)
library(zellkonverter)
library(TENxBrainData)
library(SingleCellMultiModal)
library(SingleCellExperiment)
library(DropletUtils)
library(GEOquery)
library(stringr)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

parser$add_argument(
  "--output_dir", "-o",
  dest = "output_dir", type = "character",
  help = "output directory where files will be saved",
  default = getwd(), required = TRUE
)
parser$add_argument(
  "--name", "-n",
  dest = "name", type = "character",
  help = "name of the module",
  required = TRUE
)
parser$add_argument(
  "--dataset_name",
  dest = "dataset_name", type = "character",
  help = "name of the dataset",
  choices = c("sc-mix", "be1", "cb", "1.3m"), required = TRUE
)

args <- parser$parse_args()

h5ad_path <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
clusters_truth_path <- file.path(
  args$output_dir, paste0(args$name, ".clusters_truth.tsv")
)

if (args$dataset_name == "sc-mix") {
  url <- "https://github.com/LuyiTian/sc_mixology/raw/refs/heads/master/data/sincell_with_class_5cl.RData"
  bn <- basename(url)
  raw_path <- file.path(args$output_dir, bn)
  if (!file.exists(raw_path)) {
    download.file(url, destfile = raw_path)
  }

  load(raw_path)
  sce <- sce_sc_10x_5cl_qc
  colData(sce)$clusters.truth <- colData(sce)$cell_line
  file.remove(raw_path)
} else if (args$dataset_name == "be1") {
  gse_id <- "GSE243665"
  getGEOSuppFiles(
    GEO = gse_id,
    baseDir = args$output_dir
  )

  raw_dir <- file.path(args$output_dir, gse_id)
  files <- list.files(
    raw_dir,
    pattern = "\\.(mtx|tsv)\\.gz$", full.names = TRUE
  )
  samples <- unique(str_match(basename(files), "_(.*?)_")[, 2])

  for (s in samples) {
    sample_dir <- file.path(raw_dir, s)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

    s_files <- files[grepl(paste0("_", s, "_"), basename(files))]
    file.copy(s_files, file.path(raw_dir, s), overwrite = TRUE)

    for (f in s_files) {
      bn <- basename(f)
      new_name <- sub("^[^_]+_[^_]+_", "", bn)
      dest <- file.path(sample_dir, new_name)
      file.rename(f, dest)
    }
  }

  sce_list <- lapply(samples, function(sample) {
    read10xCounts(
      samples = file.path(raw_dir, sample),
      sample.names = sample,
      col.names = TRUE
    )
  })
  sce <- do.call(cbind, sce_list)
  file.remove(raw_dir)
  colData(sce)$clusters.truth <- colData(sce)$Sample
} else if (args$dataset_name == "cb") {
  sce <- CITEseq(
    DataType = "cord_blood", modes = "*", dry.run = FALSE, version = "1.0.0",
    DataClass = "SingleCellExperiment"
  )

  gene_m <- grep("^HUMAN_", rownames(sce), value = TRUE)
  sce <- sce[gene_m, ]
  rownames(sce) <- sub("^HUMAN_", "", rownames(sce))

  colData(sce)$clusters.truth <- colData(sce)$celltype
} else if (args$dataset_name == "1.3m") {
  sce <- TENxBrainData()
  rownames(sce) <- rowData(sce)$Symbol
}

writeH5AD(sce, file = h5ad_path, compression = "gzip")
write.table(
  data.frame(
    cell_id = colnames(sce),
    truths = sce$clusters.truth
  ),
  clusters_truth_path,
  sep = "\t", quote = F, row.names = F
)
