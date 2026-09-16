#!/bin/bash
# Split one multiome sample into per-cell-type RNA BAMs and ATAC fragments.
#
# Worked example for a single sample. The full dataset was produced by running
# this over every sample x cell type, then merging across samples (step 2).
# The merged outputs are released, so this is here to document how they were
# made rather than as a step you need to run. See README.md.
#
# Usage:  bash 1_split_sample_bam.sh <sample_id> <barcode_dir> <outdir>
#
# <barcode_dir> holds one headerless CSV per cell type, <cell_type>.csv, listing
# the cell barcodes of that sample assigned to that cell type.

set -euo pipefail

source "$(dirname "$0")/../config.sh"

SAMPLE=$1
BARCODE_DIR=$2
OUTDIR=$3

# Cell Ranger ARC output for this sample.
GEX_BAM="${CELLRANGER_DIR}/${SAMPLE}/outs/gex_possorted_bam.bam"
ATAC_FRAGMENTS="${CELLRANGER_DIR}/${SAMPLE}/outs/atac_fragments.tsv.gz"

mkdir -p "${OUTDIR}/split_bams" "${OUTDIR}/split_fragments"

for BARCODE_FILE in "${BARCODE_DIR}"/*.csv; do
    CT=$(basename "${BARCODE_FILE}" .csv)

    # ── RNA: subset the gene-expression BAM to this cell type ─────────────────
    FILTER_BAM="${OUTDIR}/split_bams/${CT}_filter.bam"
    FINAL_BAM="${OUTDIR}/split_bams/${CT}_filter_final.bam"

    subset-bam_linux -b "${GEX_BAM}" -c "${BARCODE_FILE}" -o "${FILTER_BAM}" --cores 32
    samtools index "${FILTER_BAM}"

    # Keep only confidently mapped, deduplicated transcriptomic reads. Cell
    # Ranger flags these with xf:i:25; without this filter the coverage is
    # dominated by PCR duplicates.
    samtools view -h "${FILTER_BAM}" \
        | awk 'BEGIN{OFS="\t"} /^@/ {print; next} ($0 ~ /xf:i:25/) {print}' \
        | samtools view -b -o "${FINAL_BAM}"
    samtools index "${FINAL_BAM}"
    rm -f "${FILTER_BAM}" "${FILTER_BAM}.bai"

    # ── ATAC: subset the fragment file to this cell type ──────────────────────
    FRAG_TSV="${OUTDIR}/split_fragments/${CT}_atac_fragments.tsv"

    python "$(dirname "$0")/subset_fragments.py" \
        "${ATAC_FRAGMENTS}" "${BARCODE_FILE}" "${FRAG_TSV}"
    bgzip -f "${FRAG_TSV}"

    echo "Done: ${SAMPLE} ${CT}"
done
