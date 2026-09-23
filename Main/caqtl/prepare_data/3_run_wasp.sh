#!/usr/bin/env bash
# Purpose: Correct allele-dependent mapping in one donor ATAC BAM with WASP.
# Inputs:  donor BAMs, WASP HDF5 files and the hg38 BWA reference.
# Outputs: prepare_data/output/wasp/<donor>/<donor>.wasp.final.bam, index and flagstat.
# Run:     sbatch --array=0-<N-1> prepare_data/3_run_wasp.sh

#SBATCH --job-name=caqtl_wasp
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --partition=genoa-lrg-mem
#SBATCH --output=prepare_data/logs/wasp_%A_%a.out
#SBATCH --error=prepare_data/logs/wasp_%A_%a.err

set -euo pipefail

WASP_DIR="${WASP_DIR:-software/WASP}"
REFERENCE_FASTA="data/reference/genome.fa"
MANIFEST="data/manifests/donor_bams.tsv"
H5_DIR="prepare_data/output/wasp_variants"
DONOR_BAM_DIR="prepare_data/output/donor_bams"
OUT_ROOT="prepare_data/output/wasp"
task_id="${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"

donor_id="$(awk -F '\t' -v line="$((task_id + 2))" 'NR == line {print $1; exit}' "$MANIFEST")"
[[ -n "$donor_id" ]] || { echo "Error: no donor for task ${task_id}" >&2; exit 1; }

input_bam="${DONOR_BAM_DIR}/${donor_id}/${donor_id}.bam"
out_dir="${OUT_ROOT}/${donor_id}"
final_bam="${out_dir}/${donor_id}.wasp.final.bam"

[[ -f "$input_bam" ]] || { echo "Error: missing donor BAM ${input_bam}" >&2; exit 1; }
[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: missing reference FASTA ${REFERENCE_FASTA}" >&2; exit 1; }
mkdir -p "$out_dir"

echo "Step 1/5: identify SNP-overlapping reads for ${donor_id}."
python "${WASP_DIR}/mapping/find_intersecting_snps.py" \
    --is_paired_end \
    --is_sorted \
    --output_dir "$out_dir" \
    --snp_tab "${H5_DIR}/snp_tab.h5" \
    --snp_index "${H5_DIR}/snp_index.h5" \
    --haplotype "${H5_DIR}/haplotypes.h5" \
    --samples "$donor_id" \
    "$input_bam"

echo "Step 2/5: remap allele-swapped read pairs with BWA-MEM."
bwa mem -t "${SLURM_CPUS_PER_TASK:-8}" "$REFERENCE_FASTA" \
    "${out_dir}/${donor_id}.remap.fq1.gz" \
    "${out_dir}/${donor_id}.remap.fq2.gz" \
    | samtools sort -@ "${SLURM_CPUS_PER_TASK:-8}" -o "${out_dir}/realigned.bam" -

echo "Step 3/5: retain reads that remap consistently."
python "${WASP_DIR}/mapping/filter_remapped_reads.py" \
    "${out_dir}/${donor_id}.to.remap.bam" \
    "${out_dir}/realigned.bam" \
    "${out_dir}/filtered_keep.bam"

echo "Step 4/5: merge corrected reads and remove duplicate-flagged alignments."
samtools merge -u - \
    "${out_dir}/${donor_id}.keep.bam" \
    "${out_dir}/filtered_keep.bam" \
    | samtools sort -@ "${SLURM_CPUS_PER_TASK:-8}" - \
    | samtools view -b -F 0x400 -o "$final_bam"

echo "Step 5/5: index and summarize the corrected BAM."
samtools index -@ "${SLURM_CPUS_PER_TASK:-8}" "$final_bam"
samtools flagstat -@ "${SLURM_CPUS_PER_TASK:-8}" "$final_bam" \
    > "${out_dir}/${donor_id}.wasp.final.flagstat"

echo "WASP correction completed for ${donor_id}."
