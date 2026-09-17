#!/bin/bash
# Example cellSNP-lite + vireo command for one pooled library (HK2385_HK2814, donors HK2385 + HK2814), ATAC BAM.
# One such job is run per library; the Slurm loop generator is not shown.
# Reuses the per-donor VCF produced by the RNA vireo step (vireo_gt).
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=HK2385_HK2814

conda activate cellSNP

# 1. pile up informative sites on the ATAC BAM (no UMIs, permissive thresholds)
cellsnp-lite \
  -s path/to/cellranger/HK2385_HK2814/HK2385_HK2814_output/outs/atac_possorted_bam.bam \
  -b path/to/atac_barcodes/HK2385_HK2814.csv \
  -O cellsnp_out \
  -R ../../vireo_gt/outputs/HK2385_HK2814/donor_filter.vcf.gz \
  -p 16 --minMAF 0 --minCOUNT 2 --UMItag None --gzip

# 2. restrict the donor VCF to the pileup sites and assign donors
bcftools view ../../vireo_gt/outputs/HK2385_HK2814/donor_filter.vcf.gz \
  -R cellsnp_out/cellSNP.base.vcf.gz -Oz -o donor_filter_cellsnp.vcf.gz
vireo -c cellsnp_out -d donor_filter_cellsnp.vcf.gz -o vireo_out
