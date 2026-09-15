#!/usr/bin/env bash
#
# Purpose:
#   Align and quantify paired Gene Expression and Chromatin Accessibility reads
#   from one single-nucleus multiome library using Cell Ranger ARC.
#
# Input:
#   1. A run identifier.
#   2. A Cell Ranger ARC libraries CSV file.
#   3. A Cell Ranger ARC-compatible reference directory.
#   4. An output directory.
#
# Output:
#   <output_directory>/<run_id>/outs: Cell Ranger ARC count outputs.
#
#SBATCH --job-name=cellranger_arc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=80G
#SBATCH --output=cellranger_arc_%j.out
#SBATCH --error=cellranger_arc_%j.err

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <run_id> <libraries_csv> <reference_directory> <output_directory>" >&2
    exit 2
fi

readonly RUN_ID="$1"
readonly LIBRARIES_CSV_INPUT="$2"
readonly REFERENCE_DIR_INPUT="$3"
readonly OUTPUT_DIR_INPUT="$4"
readonly CELLRANGER_ARC_BIN="${CELLRANGER_ARC_BIN:-cellranger-arc}"

# Cell Ranger ARC resource parameters used in the study.
readonly LOCAL_CORES="30"
readonly LOCAL_MEMORY_GB="80"

# Check the required executable and input paths.
if ! command -v "${CELLRANGER_ARC_BIN}" >/dev/null 2>&1; then
    echo "ERROR: Cell Ranger ARC executable not found: ${CELLRANGER_ARC_BIN}" >&2
    exit 1
fi

if [[ ! -f "${LIBRARIES_CSV_INPUT}" ]]; then
    echo "ERROR: Libraries CSV not found: ${LIBRARIES_CSV_INPUT}" >&2
    exit 1
fi

if [[ ! -d "${REFERENCE_DIR_INPUT}" ]]; then
    echo "ERROR: Reference directory not found: ${REFERENCE_DIR_INPUT}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR_INPUT}"

# Resolve input paths before changing to the output directory.
readonly LIBRARIES_CSV="$(realpath "${LIBRARIES_CSV_INPUT}")"
readonly REFERENCE_DIR="$(realpath "${REFERENCE_DIR_INPUT}")"
readonly OUTPUT_DIR="$(realpath "${OUTPUT_DIR_INPUT}")"

echo "[$(date --iso-8601=seconds)] Starting Cell Ranger ARC for ${RUN_ID}."

cd "${OUTPUT_DIR}"

"${CELLRANGER_ARC_BIN}" count \
    --id="${RUN_ID}" \
    --reference="${REFERENCE_DIR}" \
    --libraries="${LIBRARIES_CSV}" \
    --localcores="${LOCAL_CORES}" \
    --localmem="${LOCAL_MEMORY_GB}"

echo "[$(date --iso-8601=seconds)] Cell Ranger ARC completed for ${RUN_ID}."
