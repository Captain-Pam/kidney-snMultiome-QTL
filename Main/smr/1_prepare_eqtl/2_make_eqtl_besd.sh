#!/usr/bin/env bash
set -euo pipefail

# Purpose:
#   Convert one cleaned SMR eQTL qfile into a dense BESD bundle.
# Inputs:
#   1) Cleaned SMR qfile.
#   2) Output prefix for the .besd, .esi, and .epi files.
# Outputs:
#   A dense SMR BESD bundle at the requested output prefix.
# Example:
#   SMR_BIN=/path/to/smr bash 2_make_eqtl_besd.sh \
#     celltype.qfile.tsv output/celltype

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <qfile.tsv> <output_prefix>" >&2
  exit 1
fi

QFILE="$1"
OUT_PREFIX="$2"
SMR_BIN="${SMR_BIN:-smr}"

if [[ ! -f "${QFILE}" ]]; then
  echo "Input qfile not found: ${QFILE}" >&2
  exit 1
fi
if ! command -v "${SMR_BIN}" >/dev/null 2>&1; then
  echo "SMR executable not found: ${SMR_BIN}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUT_PREFIX}")"

echo "Building dense BESD: ${OUT_PREFIX}"
"${SMR_BIN}" \
  --qfile "${QFILE}" \
  --make-besd-dense \
  --out "${OUT_PREFIX}"

echo "Completed dense BESD: ${OUT_PREFIX}"
