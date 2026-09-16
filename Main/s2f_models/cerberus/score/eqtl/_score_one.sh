#!/bin/bash
# One array task of 3_score_folds.sh: one fold, one chunk of variants.
# Not meant to be run directly.

set -euo pipefail

source "$(dirname "$0")/../../../config.sh"

FOLD=$(( SLURM_ARRAY_TASK_ID / N_CHUNKS ))
CHUNK=$(( SLURM_ARRAY_TASK_ID % N_CHUNKS ))

PARAMS="${CERBERUS_MODEL_DIR}/f${FOLD}c0/train/params.json"
MODEL="${CERBERUS_MODEL_DIR}/f${FOLD}c0/train/model_best.pth"
OUT_DIR="${OUT_BASE}/f${FOLD}c0/chunk_${CHUNK}"

mkdir -p "${OUT_DIR}"

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

# --stats covgene/logSED  sum predicted RNA coverage over the gene's exons and
#                         take the log ratio between alleles
# -g ${EQTL_GTF}          restricts which gene each variant is paired with; see
#                         1_prepare_gtf.py for why the full annotation is wrong
# note: the plain targets file, not the _local one — the exon set defines the
# region here, so no local window is wanted
hound_snp \
    --head 0 \
    --rc \
    --shifts "0" \
    --stats "covgene/logSED" \
    --mix_dtype bfloat16 \
    --index_start $(( CHUNK * CHUNK_SIZE )) \
    --index_end $(( CHUNK * CHUNK_SIZE + CHUNK_SIZE )) \
    -f "${CERBERUS_FASTA}" \
    -g "${EQTL_GTF}" \
    -t "${TARGETS_HUMAN}" \
    -o "${OUT_DIR}" \
    "${PARAMS}" \
    "${MODEL}" \
    "${VCF}"

echo "done: fold ${FOLD} chunk ${CHUNK}"
