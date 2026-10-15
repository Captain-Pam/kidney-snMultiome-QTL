#!/usr/bin/env bash
# Purpose: Build per-chromosome 1KGP EUR LD reference files restricted to the union of 250 kb scTWAS cis windows.
# Inputs:
#   1. Path to the public Main/sctwas directory.
#   2. 1KGP PLINK bfile prefix.
#   3. EUR keep file.
# Outputs:
#   - 2_stage2_focus/2_data/ldref/cis_ranges_250kb_union.txt.
#   - 2_stage2_focus/2_data/ldref/chr<chr>.bed/.bim/.fam.

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <sctwas_workflow_dir> <1kgp_bfile_prefix> <eur_keep_file>" >&2
  exit 1
fi

SCTWAS_DIR="$1"
KGP="$2"
EUR_KEEP="$3"
PLINK_BIN="${PLINK_BIN:-plink}"
STAGE1_DIR="$SCTWAS_DIR/1_stage1_weights"
STAGE2_DIR="$SCTWAS_DIR/2_stage2_focus"
LD_DIR="$STAGE2_DIR/2_data/ldref"
RANGE_OUT="$LD_DIR/cis_ranges_250kb_union.txt"

mkdir -p "$LD_DIR"
shopt -s nullglob
RANGE_FILES=("$STAGE1_DIR"/2_data/gene_manifests/*.cis_ranges.txt)
if [[ "${#RANGE_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: no cis range files found under $STAGE1_DIR/2_data/gene_manifests" >&2
  exit 1
fi
if [[ ! -f "$EUR_KEEP" ]]; then
  echo "ERROR: missing EUR keep file: $EUR_KEEP" >&2
  exit 1
fi
for suffix in bed bim fam; do
  [[ -f "${KGP}.${suffix}" ]] || { echo "ERROR: missing 1KGP file: ${KGP}.${suffix}" >&2; exit 1; }
done
command -v "$PLINK_BIN" >/dev/null 2>&1 || { echo "ERROR: PLINK executable not found: $PLINK_BIN" >&2; exit 1; }

cat "${RANGE_FILES[@]}" | sort -k1,1n -k2,2n -k3,3n -u > "$RANGE_OUT"
cp "$EUR_KEEP" "$LD_DIR/eur.keep"

echo "=== Build 1KGP EUR LD references for 250 kb cis-window union ==="
echo "Range files: ${#RANGE_FILES[@]}"
echo "Union ranges: $(wc -l < "$RANGE_OUT")"

for cc in $(seq 1 22); do
  echo "[chr${cc}] start"
  "$PLINK_BIN" \
    --bfile "$KGP" \
    --keep "$LD_DIR/eur.keep" \
    --chr "$cc" \
    --extract range "$RANGE_OUT" \
    --maf 0.05 \
    --make-bed \
    --out "$LD_DIR/chr${cc}" \
    > "$LD_DIR/chr${cc}.plink.stdout" \
    2> "$LD_DIR/chr${cc}.plink.stderr"
  echo "[chr${cc}] SNPs: $(wc -l < "$LD_DIR/chr${cc}.bim")"
done

echo "LD reference build complete: $LD_DIR"
