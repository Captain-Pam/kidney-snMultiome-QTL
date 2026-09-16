#!/bin/bash
# Genome-wide DeepLIFT/SHAP contribution scores over each cell type's peaks.
#
# Produces the per-peak contribution tracks that TF-MoDISco consumes
# (4_modisco.sh). Variant-level attributions for the example loci are computed
# separately, per variant, by interpret/run_shap.py.
#
# Uses the fold-0 model; the peak set is the filtered one the pipeline wrote
# during training, so regions match what the model was trained on.
#
#   conda activate chrombpnet && bash 3_contribs_bw.sh

set -euo pipefail

source "$(dirname "$0")/../config.sh"

CELLTYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)

OUT_DIR="${WORK_DIR}/chrombpnet/contribs_bw"
mkdir -p "${OUT_DIR}"

for CT in "${CELLTYPES[@]}"; do
    MODEL="${CHROMBPNET_MODEL_DIR}/${CT}_fold_0/models/chrombpnet_nobias.h5"
    REGIONS="${CHROMBPNET_MODEL_DIR}/${CT}_fold_0/auxiliary/filtered.peaks.bed"
    PREFIX="${OUT_DIR}/${CT}"

    if [[ ! -f "${MODEL}" || ! -f "${REGIONS}" ]]; then
        echo "[skip] ${CT}: missing model or filtered peaks"
        continue
    fi
    if [[ -f "${PREFIX}.profile_scores.h5" ]]; then
        echo "[skip] ${CT}: contributions exist"
        continue
    fi

    sbatch --job-name="contribs_${CT}" \
           --partition="${SLURM_GPU_PARTITION}" \
           --gres="${SLURM_GPU_GRES}" \
           --cpus-per-task=4 \
           --mem=32G \
           --output="${OUT_DIR}/${CT}.log" \
           --wrap "chrombpnet contribs_bw \
               -m ${MODEL} \
               -r ${REGIONS} \
               -g ${CHROMBPNET_FASTA} \
               -c ${GENOME_CHROM_SIZES} \
               -op ${PREFIX} \
               -pc counts profile"

    echo "submitted: ${CT}"
done
