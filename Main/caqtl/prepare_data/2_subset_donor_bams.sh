#!/usr/bin/env bash
# Purpose: Extract one donor's ATAC reads from a library-level Cell Ranger ARC BAM.
# Inputs:  data/manifests/donor_bams.tsv and data/barcodes/donor/<donor>.txt.
# Outputs: prepare_data/output/donor_bams/<donor>/<donor>.bam and index.
# Run:     sbatch --array=0-<N-1> prepare_data/2_subset_donor_bams.sh

#SBATCH --job-name=caqtl_subset
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --partition=genoa-std-mem
#SBATCH --output=prepare_data/subset_%A_%a.out
#SBATCH --error=prepare_data/subset_%A_%a.err

set -euo pipefail

MANIFEST="data/manifests/donor_bams.tsv"
BARCODE_DIR="data/barcodes/donor"
OUT_DIR="prepare_data/output/donor_bams"
task_id="${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"

record="$(awk -F '\t' -v line="$((task_id + 2))" 'NR == line {print; exit}' "$MANIFEST")"
[[ -n "$record" ]] || { echo "Error: no manifest row for task ${task_id}" >&2; exit 1; }
IFS=$'\t' read -r donor_id bam_path _ <<< "$record"

barcode_file="${BARCODE_DIR}/${donor_id}.txt"
donor_dir="${OUT_DIR}/${donor_id}"
output_bam="${donor_dir}/${donor_id}.bam"

[[ -f "$bam_path" ]] || { echo "Error: missing library BAM ${bam_path}" >&2; exit 1; }
[[ -f "$barcode_file" ]] || { echo "Error: missing donor barcode list ${barcode_file}" >&2; exit 1; }
mkdir -p "$donor_dir"

echo "Extracting donor ${donor_id} from ${bam_path}."
subset-bam \
    --bam "$bam_path" \
    --cell-barcodes "$barcode_file" \
    --out-bam "$output_bam" \
    --cores "${SLURM_CPUS_PER_TASK:-4}"
samtools index -@ "${SLURM_CPUS_PER_TASK:-4}" "$output_bam"

echo "Donor BAM completed: ${output_bam}"
