#!/usr/bin/env bash
# Purpose: Run one cell-type/chromosome/mode RASQUAL task from a task table.
# Inputs:  TASKS_FILE, SLURM_ARRAY_TASK_ID and ordered RASQUAL inputs.
# Outputs: OUTPUT_ROOT/<cell_type>/<run_label>/rasqual_*.txt.
# Run:     sbatch --array=0-<N-1> run/2_run_rasqual_array.sh

#SBATCH --job-name=caqtl_rasqual
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=3-00:00:00
#SBATCH --partition=genoa-std-mem
#SBATCH --output=run/logs/rasqual_%A_%a.out
#SBATCH --error=run/logs/rasqual_%A_%a.err

set -euo pipefail

RASQUAL_BIN="${RASQUAL_BIN:-software/rasqual/bin/rasqual}"
RASQUAL_ROOT="prepare_data/output/rasqual"
TASKS_FILE="${TASKS_FILE:-run/final_tasks.tsv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-run/output}"
task_id="${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"

record="$(awk -F '\t' -v line="$((task_id + 2))" 'NR == line {print; exit}' "$TASKS_FILE")"
[[ -n "$record" ]] || { echo "Error: no task row for index ${task_id}" >&2; exit 1; }
IFS=$'\t' read -r listed_id cell_type chromosome mode permutation model_label covariate_bin run_label <<< "$record"
[[ "$listed_id" == "$task_id" ]] || { echo "Error: task IDs are not contiguous" >&2; exit 1; }

data_dir="${RASQUAL_ROOT}/${cell_type}/data"
peak_file="${RASQUAL_ROOT}/${cell_type}/peaks/${chromosome}.peaks.tsv"
vcf_file="${data_dir}/${chromosome}.rasqual.plink_qc.vcf.gz"
sample_file="${RASQUAL_ROOT}/${cell_type}/samples.qc.txt"
counts_bin="${data_dir}/counts.mtx.bin"
size_bin="${data_dir}/size_factors.mtx.bin"
out_dir="${OUTPUT_ROOT}/${cell_type}/${run_label}"
output_file="${out_dir}/rasqual_${cell_type}_${chromosome}_${model_label}.txt"

for required_file in "$peak_file" "$vcf_file" "$sample_file" "$counts_bin" "$size_bin" "$covariate_bin"; do
    [[ -f "$required_file" ]] || { echo "Error: missing ${required_file}" >&2; exit 1; }
done
mkdir -p "$out_dir"
: > "$output_file"

mode_options=()
permutation_options=()
[[ "$mode" == "lead" ]] && mode_options=(-t)
[[ "$mode" == "lead" || "$mode" == "full" ]] || { echo "Error: invalid mode ${mode}" >&2; exit 1; }
[[ "$permutation" == "perm" ]] && permutation_options=(-r)
[[ "$permutation" == "none" || "$permutation" == "perm" ]] || { echo "Error: invalid permutation ${permutation}" >&2; exit 1; }

sample_count="$(wc -l < "$sample_file")"
threads="${SLURM_CPUS_PER_TASK:-8}"

echo "Running ${cell_type} ${chromosome} ${run_label}."
while IFS=$'\t' read -r peak_chr peak_start peak_end peak_id rasqual_index; do
    peak_midpoint=$(((peak_start + peak_end) / 2))
    cis_start=$((peak_midpoint - 10000))
    ((cis_start < 1)) && cis_start=1
    cis_end=$((peak_midpoint + 10000))
    n_cis_variants="$(tabix "$vcf_file" "${peak_chr}:${cis_start}-${cis_end}" | wc -l)"
    n_peak_variants="$(tabix "$vcf_file" "${peak_chr}:${peak_start}-${peak_end}" | wc -l)"
    ((n_cis_variants > 0)) || continue

    tabix "$vcf_file" "${peak_chr}:${cis_start}-${cis_end}" \
        | "$RASQUAL_BIN" \
            -y "$counts_bin" \
            -k "$size_bin" \
            -x "$covariate_bin" \
            -l "$n_cis_variants" \
            -m "$n_peak_variants" \
            -n "$sample_count" \
            -j "$rasqual_index" \
            -s "$peak_start" \
            -e "$peak_end" \
            -f "$peak_id" \
            -z --forced --n-threads "$threads" \
            "${mode_options[@]}" "${permutation_options[@]}" \
            >> "$output_file"
done < <(tail -n +2 "$peak_file")

echo "RASQUAL task completed: ${output_file}"
