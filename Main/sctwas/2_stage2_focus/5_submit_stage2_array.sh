#!/usr/bin/env bash
# Purpose: Submit scTWAS Stage 2 trait x cell type jobs as a Slurm array.
# Inputs:
#   - SCTWAS_DIR environment variable pointing to Main/sctwas.
#   - 2_stage2_focus/2_data/stage2_tasks.tsv from 3_make_stage2_task_table.R.
# Outputs:
#   - Slurm logs in the current code directory.
#   - 2_stage2_focus/3_outputs/<trait>/stage2/<celltype>/<celltype>.twas.tsv.

#SBATCH --job-name=sctwas_s2
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=stage2_%A_%a.out
#SBATCH --error=stage2_%A_%a.err

set -euo pipefail

SCTWAS_DIR="${SCTWAS_DIR:?Set SCTWAS_DIR to the public Main/sctwas directory.}"
STAGE_DIR="${SCTWAS_DIR}/2_stage2_focus"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
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

echo "=== scTWAS Stage 2 ==="
echo "Task ID: ${TASK_ID}"
echo "Trait: ${TRAIT}"
echo "Cell type: ${CELLTYPE}"
echo "Stage2 sumstats ready flag: ${READY_SUMSTATS}"
echo "Weights ready flag: ${READY_WEIGHTS}"
echo "Cores: ${SLURM_CPUS_PER_TASK}"

if [[ "$READY_SUMSTATS" != "TRUE" || "$READY_WEIGHTS" != "TRUE" ]]; then
  echo "ERROR: task inputs are not marked ready in $TASK_FILE." >&2
  exit 1
fi

"$RSCRIPT_BIN" "$STAGE_DIR/4_run_stage2_trait_celltype.R" \
  --workflow-dir "$SCTWAS_DIR" \
  --trait "$TRAIT" \
  --celltype "$CELLTYPE" \
  --cores "$SLURM_CPUS_PER_TASK"
