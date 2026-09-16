#!/bin/bash
# Merge per-sample splits into per-cell-type pseudobulk coverage.
#
# Takes the per-sample, per-cell-type outputs of 1_split_sample_bam.sh for one
# cell type, merges them across samples, and writes bigWig + .w5 coverage:
#
#   <ct>_atac.{bw,w5}    Tn5 insertions from the merged fragment file
#   <ct>_rna+.{bw,w5}    forward-strand RNA coverage
#   <ct>_rna-.{bw,w5}    reverse-strand RNA coverage
#
# These are the released files; see README.md. Run once per cell type, or once
# per cell type x donor to get the per-donor tracks used for genotype-stratified
# coverage plots.
#
# Usage:  bash 2_make_coverage.sh <cell_type> <split_root> <outdir>
#
# <split_root> contains one directory per sample, each holding the split_bams/
# and split_fragments/ subdirectories written by step 1.

set -euo pipefail

source "$(dirname "$0")/../config.sh"

CT=$1
SPLIT_ROOT=$2
OUTDIR=$3

FRAG_OUT="${OUTDIR}/bam_frag"
BW_OUT="${OUTDIR}/bw_w5"
mkdir -p "${FRAG_OUT}" "${BW_OUT}"

# ── ATAC ──────────────────────────────────────────────────────────────────────
MERGED_FRAG="${FRAG_OUT}/${CT}_merged_fragments.tsv.gz"
TMP_FRAG="${FRAG_OUT}/${CT}_merged_fragments_unsorted.tsv"

: > "${TMP_FRAG}"
for FRAG in "${SPLIT_ROOT}"/*/split_fragments/"${CT}"_atac_fragments.tsv.gz; do
    [[ -f "${FRAG}" ]] || continue
    zcat "${FRAG}" | grep -v "^#" >> "${TMP_FRAG}"
done

sort -k1,1 -k2,2n "${TMP_FRAG}" | bgzip -c > "${MERGED_FRAG}"
tabix -p bed "${MERGED_FRAG}"
rm -f "${TMP_FRAG}"

scatac_fragment_tools bigwig \
    -i "${MERGED_FRAG}" \
    -c "${GENOME_CHROM_SIZES}" \
    -o "${BW_OUT}/${CT}_atac.bw" \
    -x -v

# ── RNA ───────────────────────────────────────────────────────────────────────
MERGED_BAM="${FRAG_OUT}/${CT}_merged_filtered.bam"

mapfile -t BAMS < <(ls "${SPLIT_ROOT}"/*/split_bams/"${CT}"_filter_final.bam 2>/dev/null)
samtools merge -@ 32 -o - "${BAMS[@]}" | samtools sort -@ 32 -o "${MERGED_BAM}"
samtools index "${MERGED_BAM}"

# deeptools names strands by the read orientation, which is reversed relative to
# the transcript for this library chemistry: --filterRNAstrand forward gives the
# minus-strand transcript track and vice versa. The outputs are named for the
# transcript strand, which is what the model targets expect.
bamCoverage --filterRNAstrand forward --binSize 1 --skipNAs -p 32 \
    -b "${MERGED_BAM}" -o "${BW_OUT}/${CT}_rna-.bw"
bamCoverage --filterRNAstrand reverse --binSize 1 --skipNAs -p 32 \
    -b "${MERGED_BAM}" -o "${BW_OUT}/${CT}_rna+.bw"

# ── bigWig -> .w5 ─────────────────────────────────────────────────────────────
# .w5 is the HDF5 coverage format hound_data reads directly.
for TRACK in atac rna+ rna-; do
    bw_w5.py "${BW_OUT}/${CT}_${TRACK}.bw" "${BW_OUT}/${CT}_${TRACK}.w5"
done

echo "Done: ${CT}"
