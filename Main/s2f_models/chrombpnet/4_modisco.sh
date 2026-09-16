#!/bin/bash
# De novo motif discovery from ChromBPNet profile contribution scores.
#
# TF-MoDISco clusters recurring high-attribution subsequences ("seqlets") into
# motifs, then the report step matches them against a known-motif database.
#
#   conda activate chrombpnet && bash 4_modisco.sh

set -euo pipefail

source "$(dirname "$0")/../config.sh"

CELLTYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)

CONTRIBS_DIR="${WORK_DIR}/chrombpnet/contribs_bw"
OUT_DIR="${WORK_DIR}/chrombpnet/modisco"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "${OUT_DIR}"

for CT in "${CELLTYPES[@]}"; do
    H5="${CONTRIBS_DIR}/${CT}.profile_scores.h5"
    CT_OUT="${OUT_DIR}/${CT}"
    RESULTS="${CT_OUT}/modisco_results.h5"

    if [[ ! -f "${H5}" ]]; then
        echo "[skip] ${CT}: no contribution scores (run 3_contribs_bw.sh)"
        continue
    fi
    if [[ -f "${RESULTS}" ]]; then
        echo "[skip] ${CT}: modisco output exists"
        continue
    fi

    mkdir -p "${CT_OUT}"

    sbatch --job-name="modisco_${CT}" \
           --partition="${SLURM_CPU_PARTITION}" \
           --cpus-per-task=8 \
           --mem=64G \
           --time=12:00:00 \
           --output="${CT_OUT}/modisco.log" \
           --wrap "set -euo pipefail
python ${SCRIPT_DIR}/h5_to_npz.py \
    --h5 ${H5} \
    --ohe ${CT_OUT}/ohe.npz \
    --shap ${CT_OUT}/shap.npz

modisco motifs \
    -s ${CT_OUT}/ohe.npz \
    -a ${CT_OUT}/shap.npz \
    -n 50000 \
    -w 500 \
    -o ${RESULTS}

rm -f ${CT_OUT}/ohe.npz ${CT_OUT}/shap.npz

modisco report \
    -i ${RESULTS} \
    -o ${CT_OUT}/report/ \
    -m ${MODISCO_REPORT_MEME} \
    -l"

    echo "submitted: ${CT}"
done
