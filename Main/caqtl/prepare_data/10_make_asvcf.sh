#!/usr/bin/env bash
# Purpose: Create allele-specific VCFs in the canonical donor order.
# Inputs:  phased VCFs, samples.qc.txt and WASP-corrected cell-type BAMs.
# Outputs: prepare_data/output/rasqual/<cell_type>/data/chr<1-22>.rasqual.vcf.gz.
# Run:     sbatch --array=0-263 prepare_data/10_make_asvcf.sh

#SBATCH --job-name=caqtl_asvcf
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=1-00:00:00
#SBATCH --partition=genoa-std-mem
#SBATCH --output=prepare_data/asvcf_%A_%a.out
#SBATCH --error=prepare_data/asvcf_%A_%a.err

set -euo pipefail

RASQUAL_DIR="${RASQUAL_DIR:-software/rasqual}"
PHASED_DIR="data/genotype/phased"
BAM_ROOT="prepare_data/output/celltype_bams"
RASQUAL_ROOT="prepare_data/output/rasqual"
CELL_TYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)
task_id="${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"

cell_type_index=$((task_id / 22))
chromosome_number=$((task_id % 22 + 1))
cell_type="${CELL_TYPES[$cell_type_index]:?Array index must be between 0 and 263}"
chromosome="chr${chromosome_number}"

out_dir="${RASQUAL_ROOT}/${cell_type}"
data_dir="${out_dir}/data"
sample_file="${out_dir}/samples.qc.txt"
bam_list="${out_dir}/bam.list.${chromosome}.txt"
reference_vcf="${data_dir}/${chromosome}.ref.vcf.gz"
asvcf="${data_dir}/${chromosome}.rasqual.vcf.gz"

mkdir -p "$data_dir"
: > "$bam_list"
while IFS= read -r sample_id; do
    bam="${BAM_ROOT}/${cell_type}/${sample_id}.bam"
    [[ -f "$bam" ]] || { echo "Error: missing BAM ${bam}" >&2; exit 1; }
    printf '%s\n' "$bam" >> "$bam_list"
done < "$sample_file"

bcftools view \
    --samples-file "$sample_file" \
    --output-type z \
    --output "$reference_vcf" \
    "${PHASED_DIR}/${chromosome}.vcf.gz"
tabix -f -p vcf "$reference_vcf"

bcftools query -l "$reference_vcf" > "${data_dir}/${chromosome}.vcf.samples.txt"
if ! cmp -s "$sample_file" "${data_dir}/${chromosome}.vcf.samples.txt"; then
    echo "Error: VCF and samples.qc.txt orders differ for ${cell_type} ${chromosome}" >&2
    exit 1
fi

echo "Creating ASVCF for ${cell_type} ${chromosome}."
bash "${RASQUAL_DIR}/src/ASVCF/createASVCF.sh" paired_end \
    "$bam_list" "$reference_vcf" "$asvcf" atac
tabix -f -p vcf "$asvcf"
