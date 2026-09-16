#!/bin/bash
# Submit per-variant DeepLIFT/SHAP attributions for every cell type.
#
# One GPU job per cell type, each averaging over all folds. The array bound is
# derived from the cell-type list rather than hardcoded, so adding a cell type
# here cannot silently leave it unscored.
#
# Usage:  bash run_shap.sh <variants.tsv> <out_dir> [celltype ...]
#
# <variants.tsv> is the five-column variant-scorer table written by
# ../score/1_prepare_caqtl_variants.sh.
#
#   conda activate chrombpnet && bash run_shap.sh variants.tsv <dir>

set -euo pipefail

source "$(dirname "$0")/../../config.sh"

if [[ $# -lt 2 ]]; then
    echo "usage: run_shap.sh <variants.tsv> <out_dir> [celltype ...]" >&2
    exit 1
fi

VARIANTS=$1
OUT_DIR=$2
shift 2

if [[ $# -gt 0 ]]; then
    CELLTYPES=("$@")
else
    CELLTYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
fi

if [[ ! -f "${VARIANTS}" ]]; then
    echo "no such variant table: ${VARIANTS}" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

for CT in "${CELLTYPES[@]}"; do
    if [[ -f "${OUT_DIR}/${CT}_ism.h5" ]]; then
        echo "[skip] ${CT}: attributions exist"
        continue
    fi

    sbatch --job-name="shap_${CT}" \
           --partition="${SLURM_GPU_PARTITION}" \
           --gres="${SLURM_GPU_GRES}" \
           --cpus-per-task=4 \
           --mem=32G \
           --time=12:00:00 \
           --output="${OUT_DIR}/${CT}_shap.log" \
           --wrap "python ${SCRIPT_DIR}/run_shap.py \
               --celltype ${CT} \
               --variants ${VARIANTS} \
               --out_dir ${OUT_DIR}"

    echo "submitted: ${CT}"
done

echo "${#CELLTYPES[@]} cell types -> ${OUT_DIR}"
