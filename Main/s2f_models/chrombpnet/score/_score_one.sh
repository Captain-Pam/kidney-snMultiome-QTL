#!/bin/bash
# One array task of 3_score_folds.sh: score one variant list with one
# (cell type, fold) ChromBPNet model. Not meant to be run directly.

set -euo pipefail

source "$(dirname "$0")/../../config.sh"

CELLTYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)

CT=${CELLTYPES[$(( SLURM_ARRAY_TASK_ID / CHROMBPNET_N_FOLDS ))]}
FOLD=$(( SLURM_ARRAY_TASK_ID % CHROMBPNET_N_FOLDS ))

if [[ "${MODE}" == "caqtl" ]]; then
    VARIANTS="${WORK_DIR}/chrombpnet/caqtl/caqtl_variants.tsv"
    EXTRA=()
else
    VARIANTS="${WORK_DIR}/chrombpnet/eqtl/variants/${CT}_variants.tsv"
    # The eQTL lists are small; without this variant-scorer's default
    # chunking produces an empty background distribution.
    EXTRA=(-nc 1)
fi

MODEL="${CHROMBPNET_MODEL_DIR}/${CT}_fold_${FOLD}/models/chrombpnet_nobias.h5"
PEAKS="${CHROMBPNET_PEAK_DIR}/${CT}_peaks.bed"
OUT_DIR="${OUT_BASE}/${CT}/fold_${FOLD}"
OUT_PREFIX="${OUT_DIR}/${CT}_fold_${FOLD}"

if [[ ! -f "${VARIANTS}" ]]; then
    echo "[skip] ${CT} fold_${FOLD}: no variant list at ${VARIANTS}"
    exit 0
fi
if [[ -f "${OUT_PREFIX}.variant_scores.tsv" ]]; then
    echo "[skip] ${CT} fold_${FOLD}: already scored"
    exit 0
fi

mkdir -p "${OUT_DIR}"

echo "cell type: ${CT}  fold: ${FOLD}"
echo "model:     ${MODEL}"

# Fail fast rather than silently falling back to CPU, which would take days.
python - <<'PY'
import tensorflow as tf
gpus = tf.config.list_physical_devices("GPU")
if not gpus:
    raise SystemExit("No GPU detected - aborting rather than running on CPU")
print(f"GPU: {[g.name for g in gpus]}")
PY

# -sc original    the five-column variant schema written by steps 1 and 2
# -p <peaks>      normalise scores against this cell type's peak distribution,
#                 which is what makes the quantile columns comparable across
#                 cell types
# -n 10           10 shuffled-background variants per real variant
python "${VARIANT_SCORER_SRC}/variant_scoring.py" \
    -l "${VARIANTS}" \
    -g "${CHROMBPNET_FASTA}" \
    -m "${MODEL}" \
    -o "${OUT_PREFIX}" \
    -s "${GENOME_CHROM_SIZES}" \
    -sc original \
    -p "${PEAKS}" \
    -n 10 \
    --no_hdf5 \
    "${EXTRA[@]}"

echo "done: ${CT} fold_${FOLD}"
