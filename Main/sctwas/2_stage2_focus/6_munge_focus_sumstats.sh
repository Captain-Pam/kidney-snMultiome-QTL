#!/usr/bin/env bash
# Purpose: Munge standardized GWAS summary statistics for FOCUS finemapping.
# Inputs:
#   1. Path to the public Main/sctwas directory.
#   2. GWAS trait name from gwas_manifest.tsv.
# Outputs:
#   - 2_stage2_focus/2_data/sumstats_by_trait/<trait>/<trait>.focus.sumstats.gz.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <sctwas_workflow_dir> <trait>" >&2
  exit 1
fi

SCTWAS_DIR="$1"
TRAIT="$2"
FOCUS_BIN="${FOCUS_BIN:-focus}"
WD="$SCTWAS_DIR/2_stage2_focus"
INPUT="$WD/2_data/sumstats_by_trait/$TRAIT/${TRAIT}.focus_input.tsv"
OUTPUT_PREFIX="$WD/2_data/sumstats_by_trait/$TRAIT/${TRAIT}.focus"

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: missing FOCUS input: $INPUT" >&2
  exit 1
fi

command -v "$FOCUS_BIN" >/dev/null 2>&1 || { echo "ERROR: FOCUS executable not found: $FOCUS_BIN" >&2; exit 1; }

echo "=== FOCUS munge ==="
echo "Trait: $TRAIT"
echo "Input: $INPUT"
echo "Output prefix: $OUTPUT_PREFIX"

"$FOCUS_BIN" munge "$INPUT" \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --signed-sumstats Z,0 \
  --N-col N \
  -o "$OUTPUT_PREFIX"

if [[ ! -f "${OUTPUT_PREFIX}.sumstats.gz" ]]; then
  echo "ERROR: expected munged file was not created: ${OUTPUT_PREFIX}.sumstats.gz" >&2
  exit 1
fi
echo "Munged FOCUS sumstats: ${OUTPUT_PREFIX}.sumstats.gz"
