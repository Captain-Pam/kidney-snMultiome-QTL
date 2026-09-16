#!/bin/bash
# Call peaks and prepare GC-matched background regions for ChromBPNet.
#
# Runs on CPU. One pass over all cell types:
#   1. MACS2 peak calling on the pseudobulk fragment file
#   2. blacklist filtering
#   3. chrombpnet prep nonpeaks
#
#   conda activate chrombpnet && bash 1_prepare_peaks.sh

set -euo pipefail

source "$(dirname "$0")/../config.sh"

CELLTYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)

mkdir -p "${CHROMBPNET_PEAK_DIR}" "${CHROMBPNET_MACS2_DIR}"

# ── Peak calling ──────────────────────────────────────────────────────────────
for CT in "${CELLTYPES[@]}"; do
    PEAKS="${CHROMBPNET_PEAK_DIR}/${CT}_peaks.bed"
    FRAG="${FRAG_DIR}/${CT}_merged_fragments.tsv.gz"
    TMPBED="${CHROMBPNET_MACS2_DIR}/${CT}_frags_tmp.bed"

    if [[ -f "${PEAKS}" ]]; then
        echo "  [skip] ${CT}: peaks exist"
        continue
    fi
    if [[ ! -f "${FRAG}" ]]; then
        echo "  [skip] ${CT}: no fragment file at ${FRAG}"
        continue
    fi

    echo "  Calling peaks: ${CT}"
    zcat "${FRAG}" | grep -v "^#" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > "${TMPBED}"

    # ATAC settings: no model building, shift back by half the extension so the
    # 200 bp window is centred on the Tn5 insertion site. Permissive p-value —
    # ChromBPNet wants a generous peak set, not a calibrated one.
    macs2 callpeak \
        -t "${TMPBED}" \
        -f BED \
        -g hs \
        -n "${CT}" \
        --outdir "${CHROMBPNET_MACS2_DIR}" \
        --nomodel \
        --shift -100 \
        --extsize 200 \
        --keep-dup all \
        -p 0.01

    rm -f "${TMPBED}"

    # Drop blacklisted regions and non-canonical contigs. The blacklist is
    # extended by 1057 bp on each side so that no training window can reach
    # into a blacklisted region: the model sees 2114 bp of input context.
    bedtools intersect -v \
        -a "${CHROMBPNET_MACS2_DIR}/${CT}_peaks.narrowPeak" \
        -b "${BLACKLIST_EXT}" \
        | grep -P "^chr([0-9]+|X|Y|M)\t" \
        > "${PEAKS}"

    echo "  ${CT}: $(wc -l < "${PEAKS}") peaks"
done

# ── Background regions ────────────────────────────────────────────────────────
# GC-matched non-peak regions, sampled once against the fold-0 split and reused
# for every fold.
for CT in "${CELLTYPES[@]}"; do
    PEAKS="${CHROMBPNET_PEAK_DIR}/${CT}_peaks.bed"
    PREFIX="${CHROMBPNET_PEAK_DIR}/${CT}_nonpeaks"

    if [[ ! -f "${PEAKS}" ]]; then
        echo "  [skip] ${CT}: no peaks"
        continue
    fi
    if [[ -f "${PREFIX}_negatives.bed" ]]; then
        echo "  [skip] ${CT}: nonpeaks exist"
        continue
    fi

    rm -rf "${PREFIX}_auxiliary"   # stale output from an interrupted run

    echo "  Preparing nonpeaks: ${CT}"
    chrombpnet prep nonpeaks \
        -g "${CHROMBPNET_FASTA}" \
        -p "${PEAKS}" \
        -c "${GENOME_CHROM_SIZES}" \
        -st 500 \
        -fl "${SPLITS_DIR}/fold_0.json" \
        -br "${GENOME_BLACKLIST}" \
        -o "${PREFIX}"
done

echo "Peak preparation complete. Next: 2_train_folds.sh"
