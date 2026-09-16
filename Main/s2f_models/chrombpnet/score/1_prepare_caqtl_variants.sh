#!/bin/bash
# Convert the caQTL VCF to the flat table variant-scorer reads.
#
# The VCF is built once, upstream, by
# ../../cerberus/score/caqtl/1_prepare_variants.py, so both models score exactly
# the same variant set.
#
# variant-scorer's "original" schema is five headerless columns:
#   chr  pos  variant_id  allele1  allele2
#
# Usage:  bash 1_prepare_caqtl_variants.sh

set -euo pipefail

source "$(dirname "$0")/../../config.sh"

VCF="${WORK_DIR}/caqtl/caqtl_full_fdr0.1_inPeaks.vcf"
OUT="${WORK_DIR}/chrombpnet/caqtl/caqtl_variants.tsv"

mkdir -p "$(dirname "${OUT}")"

grep -v "^#" "${VCF}" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5}' > "${OUT}"

echo "Variants written: $(wc -l < "${OUT}") -> ${OUT}"
