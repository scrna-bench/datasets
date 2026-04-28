#!/usr/bin/env Rscript
# test loading of dependencies for run.R
suppressPackageStartupMessages({
  library(argparse)
  library(anndataR)
  library(TENxBrainData)
  library(SingleCellMultiModal)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(GEOquery)
  library(stringr)
})

message("OK")
