#!/bin/bash
# Evaluate every trained fold on every data fold, for both genomes.
#
# Submits an array of (model fold x genome x data fold) jobs. Each writes an
# acc.txt of per-track Pearson r and r^2. The accuracy figures read only the
# diagonal — model f<I>c0 evaluated on data fold I, its held-out test split —
# but the full grid is cheap and lets you check for fold-specific artefacts.
#
#   conda activate baskerville && bash 3_eval_folds.sh

set -euo pipefail

source "$(dirname "$0")/../config.sh"

N_DATA_FOLDS=16
N_GENOMES=2
N_TASKS=$(( CERBERUS_N_FOLDS * N_GENOMES * N_DATA_FOLDS ))

mkdir -p "${CERBERUS_MODEL_DIR}/logs"

sbatch --job-name=cerberus_eval \
       --partition="${SLURM_GPU_PARTITION}" \
       --gres="${SLURM_GPU_GRES}" \
       --cpus-per-task=4 \
       --mem=30G \
       --time=24:00:00 \
       --array="0-$(( N_TASKS - 1 ))" \
       --output="${CERBERUS_MODEL_DIR}/logs/eval_%A_%a.out" \
       --error="${CERBERUS_MODEL_DIR}/logs/eval_%A_%a.err" \
       "$(cd "$(dirname "$0")" && pwd)/_eval_one.sh"
