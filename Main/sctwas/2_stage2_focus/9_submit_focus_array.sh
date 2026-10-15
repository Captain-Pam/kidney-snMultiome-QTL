#!/usr/bin/env bash
# Purpose: Submit FOCUS finemapping trait x cell type jobs as a Slurm array.
# Inputs:
#   - SCTWAS_DIR and LD_BLOCK_BED environment variables.
#   - 2_stage2_focus/2_data/stage2_tasks.tsv from 3_make_stage2_task_table.R.
#   - Munged FOCUS GWAS files and per-cell-type FOCUS DBs.
# Outputs:
#   - Slurm logs in the current code directory.
#   - 2_stage2_focus/3_outputs/<trait>/focus/<celltype>/<celltype>.focus_all.tsv.

#SBATCH --job-name=sctwas_focus
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=focus_%A_%a.out
#SBATCH --error=focus_%A_%a.err

set -euo pipefail

SCTWAS_DIR="${SCTWAS_DIR:?Set SCTWAS_DIR to the public Main/sctwas directory.}"
LD_BLOCK_BED="${LD_BLOCK_BED:?Set LD_BLOCK_BED to the hg38 EUR LD-block BED file.}"
FOCUS_BIN="${FOCUS_BIN:-focus}"
export LD_BLOCK_BED FOCUS_BIN
STAGE_DIR="$SCTWAS_DIR/2_stage2_focus"
TASK_FILE="$STAGE_DIR/2_data/stage2_tasks.tsv"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: submit this script with sbatch --array=1-N." >&2
  exit 1
fi
if [[ ! -f "$TASK_FILE" ]]; then
  echo "ERROR: missing task table: $TASK_FILE" >&2
  exit 1
fi

N_TASKS=$(awk 'END {print NR - 1}' "$TASK_FILE")
if [[ "$SLURM_ARRAY_TASK_ID" -gt "$N_TASKS" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID exceeds task count $N_TASKS." >&2
  exit 1
fi

IFS=$'\t' read -r TASK_ID TRAIT CELLTYPE READY_SUMSTATS READY_WEIGHTS < <(
  awk -v id="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS=OFS="\t"} NR == id + 1 {print $1, $2, $3, $4, $5}' "$TASK_FILE"
)

echo "=== FOCUS array task ==="
echo "Task ID: ${TASK_ID}"
echo "Trait: ${TRAIT}"
echo "Cell type: ${CELLTYPE}"
echo "Parallel chromosomes: ${SLURM_CPUS_PER_TASK}"

bash "$STAGE_DIR/8_run_focus_trait_celltype.sh" \
  "$SCTWAS_DIR" \
  "$TRAIT" \
  "$CELLTYPE" \
  "$SLURM_CPUS_PER_TASK"
