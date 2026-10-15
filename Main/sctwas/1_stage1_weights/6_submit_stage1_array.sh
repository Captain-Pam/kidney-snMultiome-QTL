#!/usr/bin/env bash
# Purpose: Submit scTWAS Stage 1 training for 12 cell types x 22 chromosomes as a Slurm array.
# Inputs:
#   - SCTWAS_DIR environment variable pointing to Main/sctwas.
#   - 1_stage1_weights/2_data/celltypes.tsv.
#   - Prepared SCT objects, qcovar files, gene manifests, and per-chromosome genotype files.
# Outputs:
#   - Slurm logs in the current code directory.
#   - 1_stage1_weights/3_outputs/weights/<celltype>/<gene>.wgt.RDat.
#   - 1_stage1_weights/3_outputs/summaries/stage1/<celltype>.chr<chr>.stage1.tsv.

#SBATCH --job-name=sctwas_s1_250kb
#SBATCH --array=1-264
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=stage1_%A_%a.out
#SBATCH --error=stage1_%A_%a.err

set -euo pipefail

SCTWAS_DIR="${SCTWAS_DIR:?Set SCTWAS_DIR to the public Main/sctwas directory.}"
STAGE_DIR="${SCTWAS_DIR}/1_stage1_weights"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
PLINK_BIN="${PLINK_BIN:-plink}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: submit this script with sbatch so SLURM_ARRAY_TASK_ID is defined." >&2
  exit 1
fi

CELLTYPE=$(awk 'NR>1 {print $1}' "$STAGE_DIR/2_data/celltypes.tsv" | sed -n "$(( (SLURM_ARRAY_TASK_ID - 1) / 22 + 1 ))p")
CHR=$(( (SLURM_ARRAY_TASK_ID - 1) % 22 + 1 ))

echo "=== scTWAS Stage 1 ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Cell type: ${CELLTYPE}"
echo "Chromosome: ${CHR}"
echo "Cores: ${SLURM_CPUS_PER_TASK}"

"$RSCRIPT_BIN" "$STAGE_DIR/5_run_stage1_celltype_chr.R" \
  --workflow-dir "$SCTWAS_DIR" \
  --celltype "$CELLTYPE" \
  --chr "$CHR" \
  --cores "$SLURM_CPUS_PER_TASK" \
  --plink "$PLINK_BIN"
