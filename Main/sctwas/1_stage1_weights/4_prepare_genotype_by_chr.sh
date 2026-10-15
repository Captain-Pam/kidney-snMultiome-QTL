#!/usr/bin/env bash
# Purpose: Prepare per-chromosome QCed genotype files for 250 kb scTWAS Stage 1.
# Inputs:
#   1. Path to the public Main/sctwas directory.
#   2. QCed donor genotype PLINK bfile prefix.
# Outputs:
#   - 1_stage1_weights/2_data/genotype/chr<1-22>.bed/.bim/.fam with MAF >= 0.05.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <sctwas_workflow_dir> <genotype_bfile_prefix>" >&2
  exit 1
fi

SCTWAS_DIR="$1"
GENO="$2"
PLINK_BIN="${PLINK_BIN:-plink}"
OUT_DIR="${SCTWAS_DIR}/1_stage1_weights/2_data/genotype"

for suffix in bed bim fam; do
  [[ -f "${GENO}.${suffix}" ]] || { echo "Missing genotype file: ${GENO}.${suffix}" >&2; exit 1; }
done
command -v "${PLINK_BIN}" >/dev/null 2>&1 || { echo "PLINK executable not found: ${PLINK_BIN}" >&2; exit 1; }
mkdir -p "$OUT_DIR"

echo "=== Prepare per-chromosome genotype files (MAF >= 0.05) ==="
for chr in $(seq 1 22); do
  echo "[chr${chr}] PLINK subset"
  "$PLINK_BIN" --bfile "$GENO" --chr "$chr" --maf 0.05 --make-bed --out "$OUT_DIR/chr${chr}" \
    > "$OUT_DIR/chr${chr}.prepare.log" 2>&1
  printf "[chr%-2s] SNPs: %s\n" "$chr" "$(wc -l < "$OUT_DIR/chr${chr}.bim")"
done

echo "Genotype preparation complete: $OUT_DIR"
