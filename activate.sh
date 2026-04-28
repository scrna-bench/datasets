#!/bin/bash
# pixi does not run conda post-link scripts; run them here if needed
_bioc_post_link() {
    local pkg=$1 rname=$2
    if [ ! -d "$CONDA_PREFIX/lib/R/library/$rname" ]; then
        PREFIX=$CONDA_PREFIX bash "$CONDA_PREFIX/bin/.bioconductor-${pkg}-post-link.sh"
    fi
}

_bioc_post_link tenxbraindata TENxBrainData
_bioc_post_link singlecellmultimodal SingleCellMultiModal

# re-export only when pixi.lock is newer than the env file
_env_file="$PIXI_PROJECT_ROOT/envs/datasets.yml"
if [ ! -f "$_env_file" ] || [ "$PIXI_PROJECT_ROOT/pixi.lock" -nt "$_env_file" ]; then
    mkdir -p "$PIXI_PROJECT_ROOT/envs"
    _repo=$(basename "$(git -C "$PIXI_PROJECT_ROOT" rev-parse --show-toplevel)")
    pixi workspace export conda-environment --name "$_repo" "$_env_file"
fi
