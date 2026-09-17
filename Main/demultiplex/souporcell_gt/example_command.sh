#!/bin/bash
# Example souporcell command for one pooled library (HK2385_HK2814, donors HK2385 + HK2814), RNA BAM.
# One such job is run per library; the Slurm loop generator is not shown.
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=HK2385_HK2814

conda activate souporcell

# reuse the per-donor VCF prepared by the demuxlet step
zcat ../../demuxlet/outputs/HK2385_HK2814/donor_filter.reorder.vcf.gz > donor.vcf

# genotype-aware clustering seeded with the two known donors (k = 2)
souporcell_pipeline.py \
  -i path/to/cellranger/HK2385_HK2814/HK2385_HK2814_output/outs/gex_possorted_bam.bam \
  -b path/to/cellbender/HK2385_HK2814_cell_barcodes.csv \
  -f path/to/hg38.fa \
  --known_genotypes donor.vcf \
  --known_genotypes_sample_names HK2385 HK2814 \
  -t 16 -o output -k 2

rm donor.vcf
