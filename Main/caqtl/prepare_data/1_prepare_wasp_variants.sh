#!/usr/bin/env bash
# Purpose: Build WASP SNP-index and haplotype HDF5 files from phased hg38 VCFs.
# Inputs:  data/genotype/phased/chr<1-22>.vcf.gz and hg38 chromosome metadata.
# Outputs: prepare_data/output/wasp_variants/{samples.txt,snp_tab.h5,snp_index.h5,haplotypes.h5}.
# Run:     bash prepare_data/1_prepare_wasp_variants.sh

set -euo pipefail

WASP_DIR="${WASP_DIR:-software/WASP}"
PHASED_DIR="data/genotype/phased"
CHROM_INFO="data/reference/chromInfo.hg38.txt.gz"
OUT_DIR="prepare_data/output/wasp_variants"

mkdir -p "$OUT_DIR"

first_vcf="${PHASED_DIR}/chr1.vcf.gz"
[[ -f "$first_vcf" ]] || { echo "Error: missing ${first_vcf}" >&2; exit 1; }

bcftools query -l "$first_vcf" > "${OUT_DIR}/samples.txt"

"${WASP_DIR}/snp2h5/snp2h5" \
    --chrom "$CHROM_INFO" \
    --format vcf \
    --haplotype "${OUT_DIR}/haplotypes.h5" \
    --snp_index "${OUT_DIR}/snp_index.h5" \
    --snp_tab "${OUT_DIR}/snp_tab.h5" \
    --samples "${OUT_DIR}/samples.txt" \
    "${PHASED_DIR}"/chr{1..22}.vcf.gz

echo "WASP variant files written to ${OUT_DIR}."
