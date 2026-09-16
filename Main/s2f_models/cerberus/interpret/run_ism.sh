#!/bin/bash
# In-silico mutagenesis around variants, with Cerberus.
#
# For each variant, every base in a window is substituted with all three
# alternatives and the sequence re-scored. The contribution shown for a position
# is the score of the base actually present minus the mean of the three
# alternatives, which is computed at plotting time from these raw scores.
#
# Usage:  bash run_ism.sh {atac|rna} <variants.vcf> <out_dir> [fold]
#
#   atac  logSUM over a 1,024 bp window, 50 bp mutagenised
#         -> accessibility attributions
#   rna   covgene/logSED over the target gene's exons, 100 bp mutagenised
#         -> expression attributions
#
# With no fold argument this submits one job per fold, so the results can be
# averaged with ../score/ensemble_folds.py:
#
#   python ../score/ensemble_folds.py --fold_dir '<out_dir>/f{fold}c0' \
#       --out_dir '<out_dir>/ensemble' --keys ref/cov/logSUM alt/cov/logSUM
#
#   conda activate baskerville && bash run_ism.sh atac variants.vcf <dir>

set -euo pipefail

source "$(dirname "$0")/../../config.sh"

if [[ $# -lt 3 ]]; then
    echo "usage: run_ism.sh {atac|rna} <variants.vcf> <out_dir> [fold]" >&2
    exit 1
fi

MODE=$1
VCF=$2
OUT_ROOT=$3
FOLD=${4:-}

if [[ "${MODE}" != "atac" && "${MODE}" != "rna" ]]; then
    echo "mode must be 'atac' or 'rna', got '${MODE}'" >&2
    exit 1
fi

SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"

if [[ -z "${FOLD}" ]]; then
    mkdir -p "${OUT_ROOT}"
    for (( F=0; F<CERBERUS_N_FOLDS; F++ )); do
        sbatch --job-name="ism_${MODE}_f${F}" \
               --partition="${SLURM_GPU_PARTITION}" \
               --gres="${SLURM_GPU_GRES}" \
               --mem=32G \
               --time=12:00:00 \
               --output="${OUT_ROOT}/ism_f${F}.out" \
               --error="${OUT_ROOT}/ism_f${F}.err" \
               "${SELF}" "${MODE}" "${VCF}" "${OUT_ROOT}" "${F}"
    done
    echo "submitted ${CERBERUS_N_FOLDS} folds -> ${OUT_ROOT}"
    exit 0
fi

PARAMS="${CERBERUS_MODEL_DIR}/f${FOLD}c0/train/params.json"
MODEL="${CERBERUS_MODEL_DIR}/f${FOLD}c0/train/model_best.pth"
OUT_DIR="${OUT_ROOT}/f${FOLD}c0"
mkdir -p "${OUT_DIR}"

if [[ "${MODE}" == "atac" ]]; then
    # The 11-track ATAC subset keeps the output small; it must carry the
    # window=local column for --local_window to take effect.
    hound_ism_snp \
        --head 0 \
        --rc \
        --shifts "0" \
        --stats "logSUM" \
        --local_window 1024 \
        --mix_dtype bfloat16 \
        -l 50 \
        -f "${CERBERUS_FASTA}" \
        -t "${TARGETS_ATAC_SUBSET}" \
        -o "${OUT_DIR}" \
        "${PARAMS}" "${MODEL}" "${VCF}"
else
    hound_ism_snp \
        --head 0 \
        --rc \
        --shifts "0" \
        --stats "covgene/logSED" \
        --mix_dtype bfloat16 \
        -l 100 \
        -f "${CERBERUS_FASTA}" \
        -g "${EQTL_GTF}" \
        -t "${TARGETS_HUMAN}" \
        -o "${OUT_DIR}" \
        "${PARAMS}" "${MODEL}" "${VCF}"
fi

echo "done: ${MODE} fold ${FOLD} -> ${OUT_DIR}"
