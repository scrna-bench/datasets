# Datasets

Datasets for the scRNA pipelines benchmark.

## Setup

```sh
pixi install
pixi run check
```

`pixi run check` triggers the activation script, which installs any Bioconductor data packages that require post-link steps (pixi does not run conda post-link scripts automatically). Both steps are guarded and only run when needed.

## Conda environment export

The activation script also exports the resolved environment to `envs/datasets.yml` whenever `pixi.lock` is newer than the file. To export manually:

```sh
pixi run export-env
```

The environment is named after the repo root folder (e.g. `datasets`).

## Usage

```sh
pixi run Rscript run.R --output_dir <dir> --name <name> --dataset_name <sc-mix|be1|cb|1.3m>
```
