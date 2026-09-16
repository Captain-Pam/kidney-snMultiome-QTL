#!/bin/bash
# One array task of 3_eval_folds.sh. Not meant to be run directly.

set -euo pipefail

source "$(dirname "$0")/../config.sh"

N_DATA_FOLDS=16

# task -> (model fold, genome, data fold)
MODEL_FOLD=$(( SLURM_ARRAY_TASK_ID / (2 * N_DATA_FOLDS) ))
GENOME_IDX=$(( (SLURM_ARRAY_TASK_ID % (2 * N_DATA_FOLDS)) / N_DATA_FOLDS ))
DATA_FOLD=$(( SLURM_ARRAY_TASK_ID % N_DATA_FOLDS ))

MODEL_DIR="${CERBERUS_MODEL_DIR}/f${MODEL_FOLD}c0"
PARAMS="${MODEL_DIR}/train/params.json"
MODEL="${MODEL_DIR}/train/model_best.pth"

# Dataset 0 is the human head, dataset 1 the mouse head, matching the order the
# two data directories were passed to hound_train_folds.
if [[ ${GENOME_IDX} -eq 0 ]]; then
    DATA_DIR="${CERBERUS_DATA_HG38}"
else
    DATA_DIR="${CERBERUS_DATA_MM10}"
fi

OUT_DIR="${MODEL_DIR}/eval${GENOME_IDX}/fold${DATA_FOLD}"

if [[ -f "${OUT_DIR}/acc.txt" ]]; then
    echo "[skip] already evaluated: ${OUT_DIR}/acc.txt"
    exit 0
fi

mkdir -p "${OUT_DIR}"

echo "model fold ${MODEL_FOLD}, dataset ${GENOME_IDX}, data fold ${DATA_FOLD}"

# --rc averages predictions over the forward and reverse-complement strand,
# matching how the model is used for variant scoring.
hound_eval \
    --dataset "${GENOME_IDX}" \
    --split "fold${DATA_FOLD}" \
    --rc \
    -o "${OUT_DIR}" \
    "${PARAMS}" \
    "${MODEL}" \
    "${DATA_DIR}"

# Mark the model's own held-out split, which is what the accuracy figures read.
if [[ ${DATA_FOLD} -eq ${MODEL_FOLD} ]]; then
    ln -sfn "fold${DATA_FOLD}" "${MODEL_DIR}/eval${GENOME_IDX}/test"
fi
