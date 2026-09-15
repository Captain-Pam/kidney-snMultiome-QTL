#!/usr/bin/env bash
#
# Purpose:
#   Remove ambient RNA from one Cell Ranger gene-expression matrix using
#   CellBender remove-background.
#
# Input:
#   1. A Cell Ranger matrix directory containing matrix.mtx.gz.
#   2. An output prefix, including the destination directory and sample name.
#
# Output:
#   <output_prefix>.h5: CellBender-filtered gene-expression matrix.
#   <output_prefix>.tar.gz: CellBender checkpoint archive.
#
#SBATCH --job-name=cellbender
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --gres=gpu:1
#SBATCH --output=cellbender_%j.out
#SBATCH --error=cellbender_%j.err

set -euo pipefail

# Initial CellBender parameters.
# Values were adjusted for individual samples based on CellBender diagnostic results.
readonly LEARNING_RATE="0.0001"
readonly EPOCHS="150"
readonly FPR="0.01"
readonly LOW_COUNT_THRESHOLD="100"

# Read sample-specific paths from the submission command.
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_matrix_directory> <output_prefix>" >&2
    exit 2
fi

readonly INPUT_MATRIX_DIR="$1"
readonly OUTPUT_PREFIX="$2"
readonly CELLBENDER_BIN="${CELLBENDER_BIN:-cellbender}"

if [[ ! -f "${INPUT_MATRIX_DIR}/matrix.mtx.gz" ]]; then
    echo "ERROR: matrix.mtx.gz not found in ${INPUT_MATRIX_DIR}" >&2
    exit 1
fi

if ! command -v "${CELLBENDER_BIN}" >/dev/null 2>&1; then
    echo "ERROR: CellBender executable not found: ${CELLBENDER_BIN}" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUTPUT_PREFIX}")"

echo "[$(date --iso-8601=seconds)] Starting CellBender."

"${CELLBENDER_BIN}" remove-background \
    --input "${INPUT_MATRIX_DIR}" \
    --output "${OUTPUT_PREFIX}.h5" \
    --cuda \
    --learning-rate "${LEARNING_RATE}" \
    --epochs "${EPOCHS}" \
    --fpr "${FPR}" \
    --low-count-threshold "${LOW_COUNT_THRESHOLD}" \
    --checkpoint "${OUTPUT_PREFIX}.tar.gz"

echo "[$(date --iso-8601=seconds)] CellBender completed."
