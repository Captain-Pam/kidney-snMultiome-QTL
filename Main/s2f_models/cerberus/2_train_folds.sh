#!/bin/bash
# Fine-tune Cerberus on the kidney tracks: 8 folds, human + mouse heads jointly.
#
# Each fold f<I>c0 holds out data fold I as test and fold I+1 as validation, and
# is initialised from the matching fold of the pretrained Cerberus foundation
# model. Keeping the fold index aligned matters: initialising fold I from a
# checkpoint that saw fold I during pretraining would leak test data.
#
# hound_train_folds writes one submission script per fold and submits them.
#
#   conda activate baskerville && bash 2_train_folds.sh

set -euo pipefail

source "$(dirname "$0")/../config.sh"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "${CERBERUS_MODEL_DIR}"

for (( FOLD=0; FOLD<CERBERUS_N_FOLDS; FOLD++ )); do
    FOLD_DIR="${CERBERUS_MODEL_DIR}/f${FOLD}c0"
    PRETRAINED="${CERBERUS_PRETRAINED_DIR}/f${FOLD}c0/train/model_best.pth"

    if [[ ! -f "${PRETRAINED}" ]]; then
        echo "[skip] fold ${FOLD}: no pretrained checkpoint at ${PRETRAINED}"
        continue
    fi

    mkdir -p "${FOLD_DIR}"

    # params.json is shared; only the initialisation checkpoint differs per fold.
    sed "s|PRETRAINED_MODEL|${PRETRAINED}|" \
        "${SCRIPT_DIR}/params.json" > "${FOLD_DIR}/params.json"
done

# -f 8   number of folds; fold I is the test split for model f<I>c0
# the two data directories become dataset 0 (human) and dataset 1 (mouse),
# which is the --head / --dataset numbering used by every downstream script
hound_train_folds \
    -f "${CERBERUS_N_FOLDS}" \
    -o "${CERBERUS_MODEL_DIR}" \
    "${SCRIPT_DIR}/params.json" \
    "${CERBERUS_DATA_HG38}" \
    "${CERBERUS_DATA_MM10}"

echo "Training submitted for ${CERBERUS_N_FOLDS} folds -> ${CERBERUS_MODEL_DIR}"
