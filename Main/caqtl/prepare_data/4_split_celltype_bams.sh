#!/usr/bin/env bash
# Purpose: Split one WASP-corrected donor BAM into the 12 analysis cell types.
# Inputs:  corrected donor BAM and data/barcodes/cell_type/<donor>/<cell_type>.txt.
# Outputs: prepare_data/output/celltype_bams/<cell_type>/<donor>.bam and index.
# Run:     sbatch --array=0-<N-1> prepare_data/4_split_celltype_bams.sh

#SBATCH --job-name=caqtl_split
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --partition=genoa-std-mem
#SBATCH --output=prepare_data/split_%A_%a.out
#SBATCH --error=prepare_data/split_%A_%a.err

set -euo pipefail

MANIFEST="data/manifests/donor_bams.tsv"
WASP_ROOT="prepare_data/output/wasp"
BARCODE_ROOT="data/barcodes/cell_type"
OUT_ROOT="prepare_data/output/celltype_bams"
CELL_TYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
task_id="${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"

donor_id="$(awk -F '\t' -v line="$((task_id + 2))" 'NR == line {print $1; exit}' "$MANIFEST")"
[[ -n "$donor_id" ]] || { echo "Error: no donor for task ${task_id}" >&2; exit 1; }
input_bam="${WASP_ROOT}/${donor_id}/${donor_id}.wasp.final.bam"
[[ -f "$input_bam" ]] || { echo "Error: missing corrected BAM ${input_bam}" >&2; exit 1; }

for cell_type in "${CELL_TYPES[@]}"; do
    barcode_file="${BARCODE_ROOT}/${donor_id}/${cell_type}.txt"
    [[ -f "$barcode_file" ]] || { echo "Error: missing barcode file ${barcode_file}" >&2; exit 1; }
    out_dir="${OUT_ROOT}/${cell_type}"
    output_bam="${out_dir}/${donor_id}.bam"
    mkdir -p "$out_dir"

    echo "Splitting ${donor_id}: ${cell_type}."
    subset-bam \
        --bam "$input_bam" \
        --cell-barcodes "$barcode_file" \
        --out-bam "$output_bam" \
        --cores "${SLURM_CPUS_PER_TASK:-4}"
    samtools index -@ "${SLURM_CPUS_PER_TASK:-4}" "$output_bam"
done

echo "Cell-type BAM splitting completed for ${donor_id}."
