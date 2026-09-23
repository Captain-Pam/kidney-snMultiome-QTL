#!/usr/bin/env bash
# Purpose: Run one lead-mode RASQUAL covariate-tuning task.
# Inputs:  tune/tuning_tasks.tsv and SLURM_ARRAY_TASK_ID.
# Outputs: tune/output/<cell_type>/<model>/ chromosome-level RASQUAL results.
# Run:     sbatch --array=0-<N-1> tune/2_run_tuning_array.sh

#SBATCH --job-name=caqtl_tune
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=3-00:00:00
#SBATCH --partition=genoa-std-mem
#SBATCH --output=tune/logs/tune_%A_%a.out
#SBATCH --error=tune/logs/tune_%A_%a.err

set -euo pipefail

export TASKS_FILE="tune/tuning_tasks.tsv"
export OUTPUT_ROOT="tune/output"
bash run/2_run_rasqual_array.sh
