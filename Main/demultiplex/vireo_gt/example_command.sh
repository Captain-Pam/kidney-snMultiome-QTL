#!/bin/bash
# Example cellSNP-lite + vireo command for one pooled library (HK2385_HK2814, donors HK2385 + HK2814), RNA BAM.
# One such job is run per library; the Slurm loop generator is not shown.
# The donor VCF here is chr-renamed (chr1 -> 1) beforehand: see preprocess_vcf/rename.sh
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=HK2385_HK2814

conda activate cellSNP

# 1. subset the (chr-renamed) WGS VCF to the two donors
bcftools view preprocess_vcf/all_multiome.chr_rename.vcf.gz \
  | bcftools view -s HK2385,HK2814 \
  | bcftools filter -e 'GT[0]="ref" && GT[1]="ref"' \
  -Oz -o donor_filter.vcf.gz
tabix -p vcf donor_filter.vcf.gz

# 2. pile up informative sites over the called barcodes
cellsnp-lite \
  -s path/to/cellranger/HK2385_HK2814/HK2385_HK2814_output/outs/gex_possorted_bam.bam \
  -b path/to/cellbender/HK2385_HK2814_cell_barcodes.csv \
  -O cellsnp_out -R donor_filter.vcf.gz \
  -p 16 --minMAF 0.1 --minCOUNT 20 --gzip

# 3. restrict the donor VCF to the pileup sites and assign donors
bcftools view donor_filter.vcf.gz -R cellsnp_out/cellSNP.base.vcf.gz -Oz -o donor_filter_cellsnp.vcf.gz
vireo -c cellsnp_out -d donor_filter_cellsnp.vcf.gz -o vireo_out
