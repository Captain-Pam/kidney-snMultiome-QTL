#!/usr/bin/env bash

# Purpose:
#   Run one pairwise GWAS-eQTL SMR/HEIDI analysis using hg38 inputs.
# Inputs:
#   1. GWAS summary-statistics file prepared for SMR.
#   2. eQTL BESD prefix with matching .besd, .esi, and .epi files.
#   3. hg38 LD-reference prefix with matching .bed, .bim, and .fam files.
#   4. Output prefix.
#   5. Maximum allele-frequency difference accepted by SMR.
# Output:
#   SMR/HEIDI result files under the requested output prefix.
# Example:
#   SMR_BIN=/path/to/smr bash 1_run_smr.sh \
#     trait.smr.tsv cell_type_eqtl reference_hg38 results/trait_cell_type 0.2

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <gwas_summary> <eqtl_besd_prefix> <hg38_ld_prefix> <output_prefix> <diff_freq>" >&2
  exit 1
fi

GWAS_SUMMARY="$1"
BEQTL_PREFIX="$2"
LD_REF="$3"
OUT_PREFIX="$4"
DIFF_FREQ="$5"
SMR_BIN="${SMR_BIN:-smr}"

# Check the essential input files before starting SMR.
[[ -f "${GWAS_SUMMARY}" ]] || { echo "Missing GWAS file: ${GWAS_SUMMARY}" >&2; exit 1; }
for suffix in besd esi epi; do
  [[ -f "${BEQTL_PREFIX}.${suffix}" ]] || { echo "Missing eQTL file: ${BEQTL_PREFIX}.${suffix}" >&2; exit 1; }
done
for suffix in bed bim fam; do
  [[ -f "${LD_REF}.${suffix}" ]] || { echo "Missing LD-reference file: ${LD_REF}.${suffix}" >&2; exit 1; }
done
command -v "${SMR_BIN}" >/dev/null 2>&1 || { echo "SMR executable not found: ${SMR_BIN}" >&2; exit 1; }

mkdir -p "$(dirname "${OUT_PREFIX}")"

# Use the study-wide eQTL instrument threshold and the caller-supplied frequency filter.
"${SMR_BIN}" \
  --bfile "${LD_REF}" \
  --gwas-summary "${GWAS_SUMMARY}" \
  --beqtl-summary "${BEQTL_PREFIX}" \
  --peqtl-smr 5e-8 \
  --diff-freq "${DIFF_FREQ}" \
  --out "${OUT_PREFIX}"
