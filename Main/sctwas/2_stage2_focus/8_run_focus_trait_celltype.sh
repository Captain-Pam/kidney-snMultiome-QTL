#!/usr/bin/env bash
# Purpose: Run FOCUS finemapping for one GWAS trait and one cell type across autosomes.
# Inputs:
#   1. Path to the public Main/sctwas directory.
#   2. GWAS trait name.
#   3. Cell type name.
#   4. Optional number of parallel chromosome jobs.
#   - LD_BLOCK_BED must point to the hg38 EUR LD-block BED file.
# Outputs:
#   - 2_stage2_focus/3_outputs/<trait>/focus/<celltype>/finemap_chr<chr>.focus.tsv.
#   - 2_stage2_focus/3_outputs/<trait>/focus/<celltype>/<celltype>.focus_all.tsv.

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 <sctwas_workflow_dir> <trait> <celltype> [parallel_jobs]" >&2
  exit 1
fi

SCTWAS_DIR="$1"
TRAIT="$2"
CELLTYPE="$3"
PARALLEL_JOBS="${4:-4}"
FOCUS_BIN="${FOCUS_BIN:-focus}"
LD_BLOCK_BED="${LD_BLOCK_BED:?Set LD_BLOCK_BED to the hg38 EUR LD-block BED file.}"
WD="$SCTWAS_DIR/2_stage2_focus"
GWAS="$WD/2_data/sumstats_by_trait/$TRAIT/${TRAIT}.focus.sumstats.gz"
DB="$WD/3_outputs/focus_db/${CELLTYPE}.db"
LD_DIR="$WD/2_data/ldref"
OUT_DIR="$WD/3_outputs/$TRAIT/focus/$CELLTYPE"

if [[ ! -f "$GWAS" ]]; then
  echo "ERROR: missing munged FOCUS GWAS: $GWAS" >&2
  exit 1
fi
if [[ ! -f "$DB" ]]; then
  echo "ERROR: missing FOCUS DB: $DB" >&2
  echo "Build it first with 7_build_focus_db_celltype.R." >&2
  exit 1
fi
if [[ ! -f "$LD_BLOCK_BED" ]]; then
  echo "ERROR: missing FOCUS LD block BED: $LD_BLOCK_BED" >&2
  exit 1
fi
command -v "$FOCUS_BIN" >/dev/null 2>&1 || { echo "ERROR: FOCUS executable not found: $FOCUS_BIN" >&2; exit 1; }
[[ "$PARALLEL_JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: parallel_jobs must be a positive integer." >&2; exit 1; }

mkdir -p "$OUT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "=== FOCUS finemap ==="
echo "Trait: $TRAIT"
echo "Cell type: $CELLTYPE"
echo "Parallel chromosomes: $PARALLEL_JOBS"

FAIL=0
PIDS=()
for cc in $(seq 1 22); do
  LD_PREFIX="$LD_DIR/chr${cc}"
  if [[ ! -f "${LD_PREFIX}.bed" ]]; then
    echo "ERROR: missing LD reference for chr${cc}: ${LD_PREFIX}.bed" >&2
    exit 1
  fi
  "$FOCUS_BIN" finemap "$GWAS" "$LD_PREFIX" "$DB" \
    --chr "$cc" \
    --locations "$LD_BLOCK_BED" \
    --tissue "$CELLTYPE" \
    --p-threshold 5e-8 \
    -o "$OUT_DIR/finemap_chr${cc}" \
    > "$OUT_DIR/finemap_chr${cc}.runlog" \
    2>&1 &
  PIDS+=("$!")
  if [[ "${#PIDS[@]}" -ge "$PARALLEL_JOBS" ]]; then
    for pid in "${PIDS[@]}"; do
      wait "$pid" || FAIL=1
    done
    PIDS=()
  fi
done

for pid in "${PIDS[@]}"; do
  wait "$pid" || FAIL=1
done
if [[ "$FAIL" -ne 0 ]]; then
  echo "ERROR: one or more FOCUS chromosome jobs failed. Check $OUT_DIR/*.runlog." >&2
  exit 1
fi

for cc in $(seq 1 22); do
  if [[ ! -f "$OUT_DIR/finemap_chr${cc}.focus.tsv" ]]; then
    echo "ERROR: missing expected FOCUS output: $OUT_DIR/finemap_chr${cc}.focus.tsv" >&2
    exit 1
  fi
done

OUT_ALL="$OUT_DIR/${CELLTYPE}.focus_all.tsv"
FOCUS_HEADER=$'block\tens_gene_id\tens_tx_id\tmol_name\ttissue\tref_name\ttype\tchrom\ttx_start\ttx_stop\tblock_genes\ttrait\tinference_pop1\tinter_z_pop1\tcv.R2_pop1\tcv.R2.pval_pop1\tldregion_pop1\ttwas_z_pop1\tpips_pop1\tin_cred_set_pop1'
TMP_OUT="${OUT_ALL}.tmp"
printf "%s\n" "$FOCUS_HEADER" > "$TMP_OUT"
NONEMPTY_CHR=0
EMPTY_CHR=0
for cc in $(seq 1 22); do
  CHR_OUT="$OUT_DIR/finemap_chr${cc}.focus.tsv"
  if [[ -s "$CHR_OUT" ]]; then
    FIRST_LINE=$(head -n 1 "$CHR_OUT")
    if [[ "$FIRST_LINE" != block$'\t'* ]]; then
      echo "ERROR: unexpected FOCUS header in $CHR_OUT" >&2
      exit 1
    fi
    tail -n +2 "$CHR_OUT" >> "$TMP_OUT"
    NONEMPTY_CHR=$((NONEMPTY_CHR + 1))
  else
    EMPTY_CHR=$((EMPTY_CHR + 1))
  fi
done
mv "$TMP_OUT" "$OUT_ALL"
echo "FOCUS combined rows: $(( $(wc -l < "$OUT_ALL") - 1 ))"
echo "FOCUS non-empty chromosome files: $NONEMPTY_CHR"
echo "FOCUS empty chromosome files: $EMPTY_CHR"
echo "FOCUS output: $OUT_ALL"
