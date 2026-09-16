#!/bin/bash
# Score variants with every ChromBPNet fold model.
#
# Submits a SLURM array over (cell type x fold). Each task scores the whole
# variant list with one fold's bias-free model; 4_ensemble.py then averages.
#
# Usage:  bash 3_score_folds.sh {caqtl|eqtl}
#
#   caqtl  one variant list shared by all cell types
#   eqtl   a per-cell-type list, written by 2_prepare_eqtl_variants.py
#
#   conda activate chrombpnet && bash 3_score_folds.sh caqtl

set -euo pipefail

source "$(dirname "$0")/../../config.sh"

if [[ $# -lt 1 ]]; then
    echo "usage: 3_score_folds.sh {caqtl|eqtl}" >&2
    exit 1
fi
MODE=$1
if [[ "${MODE}" != "caqtl" && "${MODE}" != "eqtl" ]]; then
    echo "mode must be 'caqtl' or 'eqtl', got '${MODE}'" >&2
    exit 1
fi

CELLTYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
N_TASKS=$(( ${#CELLTYPES[@]} * CHROMBPNET_N_FOLDS ))

OUT_BASE="${WORK_DIR}/chrombpnet/${MODE}/scores_folds"
mkdir -p "${OUT_BASE}"

sbatch --job-name="cbpnet_score_${MODE}" \
       --partition="${SLURM_GPU_PARTITION}" \
       --gres="${SLURM_GPU_GRES}" \
       --mem=32G \
       --array="0-$(( N_TASKS - 1 ))%12" \
       --output="${OUT_BASE}/job_%A_%a.out" \
       --error="${OUT_BASE}/job_%A_%a.err" \
       --export=ALL,MODE="${MODE}",OUT_BASE="${OUT_BASE}" \
       "$(cd "$(dirname "$0")" && pwd)/_score_one.sh"
