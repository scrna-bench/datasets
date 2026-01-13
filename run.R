#!/usr/bin/env Rscript

library(argparse)
library(zellkonverter)
library(TENxBrainData)
library(SingleCellMultiModal)

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

if (args$dataset_name == "sc-mix") {
  url <- "https://github.com/LuyiTian/sc_mixology/raw/refs/heads/master/data/sincell_with_class_5cl.RData"
  bn <- basename(url)

  raw_path <- file.path(args$output_dir, bn)

  if (!file.exists(raw_path)) {
    download.file(url, destfile = raw_path)
  }

  load(raw_path)
  writeH5AD(sce_sc_10x_5cl_qc, file = h5ad_path, compression = "gzip")
  file.remove(raw_path)
} else if (args$dataset_name == "cb") {
  sce <- CITEseq(
    DataType = "cord_blood", modes = "*", dry.run = FALSE, version = "1.0.0",
    DataClass = "SingleCellExperiment"
  )

  gene_m <- grep("^HUMAN_", rownames(sce), value = TRUE)
  sce <- sce[gene_m, ]
  rownames(sce) <- sub("^HUMAN_", "", rownames(sce))

  writeH5AD(sce, file = h5ad_path, compression = "gzip")
} else if (args$dataset_name == "1.3m") {
  sce <- TENxBrainData()
  writeH5AD(sce, file = h5ad_path, compression = "gzip")
}
