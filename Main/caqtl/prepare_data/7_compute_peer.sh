#!/usr/bin/env bash
# Purpose: Compute candidate PEER factors from normalized chr1 peak matrices.
# Inputs:  prepare_data/output/sample_qc/<cell_type>/<cell_type>.chr1_peak.bed.gz.
# Outputs: prepare_data/output/peer/<cell_type>_peer<k>.PEER_covariates.txt.
# Run:     bash prepare_data/7_compute_peer.sh

set -euo pipefail

GTEX_QTL_DIR="${GTEX_QTL_DIR:-software/gtex-pipeline/qtl}"
INPUT_ROOT="prepare_data/output/sample_qc"
OUTPUT_ROOT="prepare_data/output/peer"
CELL_TYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
PEER_COUNTS=(1 2 3 4 5 10)

mkdir -p "$OUTPUT_ROOT" prepare_data/logs

for cell_type in "${CELL_TYPES[@]}"; do
    input_bed="${INPUT_ROOT}/${cell_type}/${cell_type}.chr1_peak.bed.gz"
    [[ -f "$input_bed" ]] || { echo "Error: missing PEER input ${input_bed}" >&2; exit 1; }

    sbatch \
        --job-name="peer_${cell_type}" \
        --output="prepare_data/logs/peer_${cell_type}.out" \
        --error="prepare_data/logs/peer_${cell_type}.err" \
        --time=24:00:00 \
        --mem=32G \
        --wrap="set -euo pipefail; for k in ${PEER_COUNTS[*]}; do Rscript ${GTEX_QTL_DIR}/src/run_PEER.R ${input_bed} ${OUTPUT_ROOT}/${cell_type}_peer\${k} \${k}; done"
done
