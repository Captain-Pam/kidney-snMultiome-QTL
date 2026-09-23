#!/usr/bin/env bash
# Purpose: Merge RASQUAL chromosome files in fixed chr1-chr22 order.
# Inputs:  run/output/<cell_type>/<run_label>/rasqual_*.txt.
# Outputs: run/merged/<cell_type>/<cell_type>_<run_label>.txt.
# Run:     bash run/3_merge_results.sh

set -euo pipefail

INPUT_ROOT="run/output"
OUTPUT_ROOT="run/merged"
CELL_TYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
RUN_LABELS=(lead_observed lead_perm1 lead_perm2 full_observed)

for cell_type in "${CELL_TYPES[@]}"; do
    out_dir="${OUTPUT_ROOT}/${cell_type}"
    mkdir -p "$out_dir"
    for run_label in "${RUN_LABELS[@]}"; do
        output_file="${out_dir}/${cell_type}_${run_label}.txt"
        : > "$output_file"
        for chromosome_number in {1..22}; do
            chromosome="chr${chromosome_number}"
            input_file="${INPUT_ROOT}/${cell_type}/${run_label}/rasqual_${cell_type}_${chromosome}_final_peer2_pc0.txt"
            [[ -f "$input_file" ]] || { echo "Error: missing ${input_file}" >&2; exit 1; }
            cat "$input_file" >> "$output_file"
        done
        echo "Merged ${output_file}."
    done
done
