#!/usr/bin/env bash
# Purpose: Submit the 60-chunk Open4Gene analysis as a Slurm array.
# Inputs: INPUT_DIR with prepared chunk data; OUTPUT_DIR for chunk results.
# Outputs: res_<chunk>.tsv and status_<chunk>.tsv files in OUTPUT_DIR, plus Slurm logs here.
# Example: INPUT_DIR=/path/to/input OUTPUT_DIR=/path/to/output sbatch Main/peak_gene_linkage/2_submit_open4gene_array.sh

#SBATCH --job-name=open4gene
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --time=5-00:00:00
#SBATCH --account=ksusztak-nephrobase
#SBATCH --array=1-60%10
#SBATCH --output=open4gene_%A_%a.out
#SBATCH --error=open4gene_%A_%a.err

set -euo pipefail

# Require the prepared input and output locations from the submission environment.
: "${INPUT_DIR:?Set INPUT_DIR to the prepared input directory}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to the result directory}"

# Slurm runs a spool copy, so use the directory from which sbatch was invoked.
# Fall back to the script directory when this file is run directly.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR=$SLURM_SUBMIT_DIR
else
  SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fi
RUNNER="${SCRIPT_DIR}/1_run_open4gene_chunk.R"
RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}
PYTHON_BIN=${PYTHON_BIN:-}
chunk_name=$(printf 'chunk_%03d' "$SLURM_ARRAY_TASK_ID")

# Report the assigned chunk and resources before starting the analysis.
printf 'Current chunk: %s\nTask ID: %s\nNode: %s\nCPUs: %s\nMemory: 200G\n' \
  "$chunk_name" "$SLURM_ARRAY_TASK_ID" "${SLURMD_NODENAME:-unknown}" "${SLURM_CPUS_PER_TASK:-20}"

# Build arguments safely and include the optional Python executable only when set.
args=(
  "$RUNNER"
  --chunk-name "$chunk_name"
  --input-dir "$INPUT_DIR"
  --output-dir "$OUTPUT_DIR"
  --batches 20
  --cores "${SLURM_CPUS_PER_TASK:-20}"
)
if [[ -n "$PYTHON_BIN" ]]; then
  args+=(--python "$PYTHON_BIN")
fi

"$RSCRIPT_BIN" "${args[@]}"
