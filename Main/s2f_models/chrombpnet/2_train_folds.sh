#!/bin/bash
# Train bias-factorized ChromBPNet models: 12 cell types x 5 chromosome folds.
#
# Submits one GPU job per (cell type, fold). Each job runs the full chrombpnet
# pipeline, which trains the accessibility component on top of a frozen,
# pre-trained Tn5 bias model and writes chrombpnet_nobias.h5 — the bias-free
# model every downstream stage uses.
#
#   conda activate chrombpnet && bash 2_train_folds.sh

set -euo pipefail

source "$(dirname "$0")/../config.sh"

CELLTYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)

mkdir -p "${CHROMBPNET_MODEL_DIR}"

for (( FOLD=0; FOLD<CHROMBPNET_N_FOLDS; FOLD++ )); do
    FOLD_JSON="${SPLITS_DIR}/fold_${FOLD}.json"
    # The Tn5 bias model is pre-trained on HepG2 (ENCODE ENCSR291GJU) and held
    # fixed during training. It is fold-matched so that no bias-model training
    # chromosome leaks into this fold's test set.
    BIAS="${BIAS_MODEL_DIR}/fold_${FOLD}/model.bias.fold_${FOLD}.ENCSR291GJU.h5"

    for CT in "${CELLTYPES[@]}"; do
        PEAKS="${CHROMBPNET_PEAK_DIR}/${CT}_peaks.bed"
        NONPEAKS="${CHROMBPNET_PEAK_DIR}/${CT}_nonpeaks_negatives.bed"
        FRAG="${FRAG_DIR}/${CT}_merged_fragments.tsv.gz"
        OUTDIR="${CHROMBPNET_MODEL_DIR}/${CT}_fold_${FOLD}"

        if [[ ! -f "${PEAKS}" || ! -f "${NONPEAKS}" || ! -f "${FRAG}" ]]; then
            echo "[skip] ${CT} fold_${FOLD}: missing inputs (run 1_prepare_peaks.sh)"
            continue
        fi
        if [[ -f "${OUTDIR}/models/chrombpnet_nobias.h5" ]]; then
            echo "[skip] ${CT} fold_${FOLD}: model exists"
            continue
        fi

        mkdir -p "${OUTDIR}"

        sbatch --job-name="cbpnet_${CT}_f${FOLD}" \
               --partition="${SLURM_GPU_PARTITION}" \
               --gres="${SLURM_GPU_GRES}" \
               --cpus-per-task=4 \
               --mem=32G \
               --time=48:00:00 \
               --output="${OUTDIR}/train.out" \
               --error="${OUTDIR}/train.err" \
               --wrap "chrombpnet pipeline \
                   -ifrag ${FRAG} \
                   -d ATAC \
                   -g ${CHROMBPNET_FASTA} \
                   -c ${GENOME_CHROM_SIZES} \
                   -p ${PEAKS} \
                   -n ${NONPEAKS} \
                   -fl ${FOLD_JSON} \
                   -b ${BIAS} \
                   -o ${OUTDIR}"

        echo "submitted: ${CT} fold_${FOLD}"
    done
done
