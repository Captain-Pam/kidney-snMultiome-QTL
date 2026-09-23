#!/usr/bin/env bash
# Purpose: Build final observed, permutation and full-mode RASQUAL tasks.
# Inputs:  ordered RASQUAL inputs under prepare_data/output/rasqual/.
# Outputs: run/final_tasks.tsv.
# Run:     bash run/1_make_tasks.sh

set -euo pipefail

RASQUAL_ROOT="prepare_data/output/rasqual"
TASK_FILE="run/final_tasks.tsv"
CELL_TYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
RUN_SPECS=(
    "lead:none:lead_observed"
    "lead:perm:lead_perm1"
    "lead:perm:lead_perm2"
    "full:none:full_observed"
)

printf 'task_id\tcell_type\tchromosome\tmode\tpermutation\tmodel_label\tcovariate_bin\trun_label\n' > "$TASK_FILE"
task_count=0
for cell_type in "${CELL_TYPES[@]}"; do
    covariate="${RASQUAL_ROOT}/${cell_type}/data/covar_final_peer2_pc0.bin"
    for chromosome_number in {1..22}; do
        chromosome="chr${chromosome_number}"
        for run_spec in "${RUN_SPECS[@]}"; do
            IFS=: read -r mode permutation run_label <<< "$run_spec"
            printf '%d\t%s\t%s\t%s\t%s\tfinal_peer2_pc0\t%s\t%s\n' \
                "$task_count" "$cell_type" "$chromosome" "$mode" "$permutation" \
                "$covariate" "$run_label" >> "$TASK_FILE"
            task_count=$((task_count + 1))
        done
    done
done

echo "Wrote ${task_count} final tasks to ${TASK_FILE}."
echo "Submit with: sbatch --array=0-$((task_count - 1)) run/2_run_rasqual_array.sh"
