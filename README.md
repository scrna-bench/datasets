# Datasets

Datasets for the scRNA pipelines benchmark, one omnibenchmark (api 0.8)
entrypoint per dataset:

| entrypoint | dataset | truth labels |
|---|---|---|
| `sc_mix` | scMixology, five cell lines (10x) | cell line |
| `be1` | GEO GSE243665 | sample |
| `cb` | cord blood CITE-seq (SingleCellMultiModal) | cell type |
| `tenx_1_3m` | 10x 1.3M mouse brain (TENxBrainData) | library (no biological meaning) |

Each writes `data_h5ad` (sparse counts, QC thresholds in `uns`),
`data_clusters_truth`, `data_clusters_truth_num` and `data_events`
(obkit phase events). Shared code lives in `src/`.

## Environment

`envs/r/` is an [rv](https://github.com/A2-ai/rv) project. CI
(`.github/workflows/containers.yml`) builds it into
`ghcr.io/scrna-bench/datasets:r-<sha>`, which the benchmark plan references.
