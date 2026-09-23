#!/usr/bin/env bash
# Purpose: Build lead-mode RASQUAL tasks for PEER-factor and genotype-PC tuning.
# Inputs:  ordered RASQUAL inputs under prepare_data/output/rasqual/.
# Outputs: tune/tuning_tasks.tsv.
# Run:     bash tune/1_make_tuning_tasks.sh

set -euo pipefail

RASQUAL_ROOT="prepare_data/output/rasqual"
TASK_FILE="tune/tuning_tasks.tsv"
CELL_TYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
PEER_COUNTS=(0 1 2 3 4 5 10)
PC_COUNTS=(1 2 3 4 5 10 15)

printf 'task_id\tcell_type\tchromosome\tmode\tpermutation\tmodel_label\tcovariate_bin\trun_label\n' > "$TASK_FILE"
task_count=0

for cell_type in "${CELL_TYPES[@]}"; do
    data_dir="${RASQUAL_ROOT}/${cell_type}/data"
    for chromosome_number in {1..22}; do
        chromosome="chr${chromosome_number}"
        for n_peer in "${PEER_COUNTS[@]}"; do
            model="peer${n_peer}_pc0"
            covariate="${data_dir}/covar_base_peer${n_peer}_pc0.bin"
            printf '%d\t%s\t%s\tlead\tnone\t%s\t%s\t%s\n' \
                "$task_count" "$cell_type" "$chromosome" "$model" "$covariate" "$model" \
                >> "$TASK_FILE"
            task_count=$((task_count + 1))
        done
        for n_pc in "${PC_COUNTS[@]}"; do
            model="peer2_pc${n_pc}"
            covariate="${data_dir}/covar_base_peer2_pc${n_pc}.bin"
            printf '%d\t%s\t%s\tlead\tnone\t%s\t%s\t%s\n' \
                "$task_count" "$cell_type" "$chromosome" "$model" "$covariate" "$model" \
                >> "$TASK_FILE"
            task_count=$((task_count + 1))
        done
    done
done

echo "Wrote ${task_count} tuning tasks to ${TASK_FILE}."
echo "Submit with: sbatch --array=0-$((task_count - 1)) tune/2_run_tuning_array.sh"
