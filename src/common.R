# Shared by every dataset entrypoint: argument handling, QC threshold table,
# and writing the stage outputs. Source with chdir = TRUE.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

source("cli.R")
source("timing.R")

# Parse the CLI, start the event log, and return the output paths.
dataset_args <- function() {
  args <- parse_cli()
  out <- list(
    h5ad = output_path(args, "data_h5ad", ".h5ad"),
    truth = output_path(args, "data_clusters_truth", ".clusters_truth.tsv"),
    truth_num = output_path(args, "data_clusters_truth_num", ".clusters_truth_num.txt"),
    events = args$outputs$data_events %||% file.path(args$output_dir, "obkit-events.jsonl")
  )
  timing_init(out$events)
  out$work_dir <- dirname(out$h5ad)
  out
}

# Suggested QC thresholds, stored in metadata(sce)$qc_thresholds and used by
# the pipelines' "manual" filter.
make_qc_df <- function(
  nFeature_min = NA_real_, nFeature_max = NA_real_,
  nCount_min = NA_real_, nCount_max = NA_real_,
  percent_mt_min = NA_real_, percent_mt_max = NA_real_
) {
  data.frame(
    metric = c("nFeature", "nCount", "percent.mt"),
    min = c(nFeature_min, nCount_min, percent_mt_min),
    max = c(nFeature_max, nCount_max, percent_mt_max),
    stringsAsFactors = FALSE
  )
}

# Keep sparse counts only, attach truth labels and QC thresholds, drop cells
# without a label, then write the h5ad, the truth table and the number of
# truth clusters.
write_dataset <- function(sce, truth, qc, out) {
  sce <- with_phase("finalize", {
    for (a in setdiff(assayNames(sce), "counts")) assay(sce, a) <- NULL
    counts(sce) <- as(counts(sce), "dgCMatrix")
    colData(sce)$clusters.truth <- truth
    metadata(sce) <- list(qc_thresholds = qc)
    sce[, !is.na(colData(sce)$clusters.truth)]
  })

  with_phase("write", {
    adata <- anndataR::as_AnnData(sce, x_mapping = "counts")
    anndataR::write_h5ad(adata, path = out$h5ad)
    write.table(
      data.frame(cell_id = colnames(sce), truths = sce$clusters.truth),
      out$truth,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    writeLines(as.character(length(unique(sce$clusters.truth))), con = out$truth_num)
  })
  invisible(sce)
}
