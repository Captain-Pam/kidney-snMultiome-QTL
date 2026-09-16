#!/bin/bash
# One array task of 2_score_folds.sh: one fold, one chunk of variants.
# Not meant to be run directly.

set -euo pipefail

source "$(dirname "$0")/../../../config.sh"

FOLD=$(( SLURM_ARRAY_TASK_ID / N_CHUNKS ))
CHUNK=$(( SLURM_ARRAY_TASK_ID % N_CHUNKS ))

PARAMS="${CERBERUS_MODEL_DIR}/f${FOLD}c0/train/params.json"
MODEL="${CERBERUS_MODEL_DIR}/f${FOLD}c0/train/model_best.pth"
OUT_DIR="${OUT_BASE}/local_f${FOLD}c0/chunk_${CHUNK}"

mkdir -p "${OUT_DIR}"

# hound_snp records how far it got, so an interrupted array can be resubmitted.
if [[ -f "${OUT_DIR}/scores.h5" ]] && python - "${OUT_DIR}/scores.h5" <<'PY'
import sys, h5py
with h5py.File(sys.argv[1]) as f:
    sys.exit(0 if f["progress_status"][()].decode() == "completed" else 1)
PY
then
    echo "[skip] fold ${FOLD} chunk ${CHUNK}: complete"
    exit 0
fi

echo "fold ${FOLD}, chunk ${CHUNK}"

# --head 0            human head
# --stats logSUM      log of summed predicted coverage
# --local_window 1024 restrict the sum to +/-512 bp around the variant. This
#                     only takes effect because ${TARGETS_HUMAN_LOCAL} carries a
#                     `window` column set to "local"; with the plain targets
#                     file the flag is silently ignored.
# --rc --shifts 0     average forward and reverse-complement, no offsets
hound_snp \
    --head 0 \
    --rc \
    --shifts "0" \
    --stats "logSUM" \
    --mix_dtype bfloat16 \
    --local_window 1024 \
    --index_start $(( CHUNK * CHUNK_SIZE )) \
    --index_end $(( CHUNK * CHUNK_SIZE + CHUNK_SIZE )) \
    -f "${CERBERUS_FASTA}" \
    -t "${TARGETS_HUMAN_LOCAL}" \
    -o "${OUT_DIR}" \
    "${PARAMS}" \
    "${MODEL}" \
    "${VCF}"

echo "done: fold ${FOLD} chunk ${CHUNK}"
