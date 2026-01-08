#!/usr/bin/env Rscript

library(argparse)
library(zellkonverter)

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
  required = TRUE,
)
parser$add_argument(
  "--dataset_name",
  choices = c("sc-mix"), type = "character",
  help = "name of the dataset",
  required = TRUE
)

args <- parser$parse_args()

if (args$dataset_generator == "sc-mix") {
  url <- "https://github.com/LuyiTian/sc_mixology/raw/refs/heads/master/data/sincell_with_class_5cl.RData"
  bn <- basename(url)

  raw_path <- file.path(args$output_dir, bn)
  h5ad_path <- file.path(args$output_dir, args$name)

  if (!file.exists(raw_path)) {
    download.file(url, destfile = raw_path)
  }

  load(raw_path)

  writeH5AD(sce_sc_10x_5cl_qc, file = h5ad_path)
}
