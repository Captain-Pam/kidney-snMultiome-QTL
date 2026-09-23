#!/usr/bin/env bash
# Purpose: Restrict allele-specific VCFs to allele-aware WGS-QC variants.
# Inputs:  chr*.rasqual.vcf.gz and data/genotype/qc_variants.txt.
# Outputs: chr*.rasqual.plink_qc.vcf.gz and tabix indexes.
# Run:     bash prepare_data/11_filter_variants.sh

set -euo pipefail

RASQUAL_ROOT="prepare_data/output/rasqual"
QC_VARIANTS="data/genotype/qc_variants.txt"
CELL_TYPES=(CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT)

[[ -s "$QC_VARIANTS" ]] || { echo "Error: missing QC variant list ${QC_VARIANTS}" >&2; exit 1; }

for cell_type in "${CELL_TYPES[@]}"; do
    data_dir="${RASQUAL_ROOT}/${cell_type}/data"
    for chromosome_number in {1..22}; do
        chromosome="chr${chromosome_number}"
        input_vcf="${data_dir}/${chromosome}.rasqual.vcf.gz"
        output_vcf="${data_dir}/${chromosome}.rasqual.plink_qc.vcf.gz"
        [[ -f "$input_vcf" ]] || { echo "Error: missing ${input_vcf}" >&2; exit 1; }

        echo "Filtering ${cell_type} ${chromosome}."
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou "$input_vcf" \
            | bcftools view -i "ID=@${QC_VARIANTS}" -Oz -o "$output_vcf"
        tabix -f -p vcf "$output_vcf"
    done
done
