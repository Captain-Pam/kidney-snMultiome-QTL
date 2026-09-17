#!/bin/bash
# Example demuxlet command for one pooled library (HK2385_HK2814, donors HK2385 + HK2814).
# One such job is run per library; the Slurm loop generator is not shown.
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --job-name=dem_HK2385_HK2814

# 1. subset the WGS VCF to the two pooled donors and drop uninformative / low-quality sites
bcftools view path/to/wgs_donors.vcf.gz \
  | bcftools view -s HK2385,HK2814 \
  | bcftools filter -e 'GT[0]="ref" && GT[1]="ref"' \
  | bcftools view -e 'GT[0]="AA" && GT[1]="AA"' \
  | bcftools view -e 'GT[0]="RA" && GT[1]="RA"' \
  | bcftools view -e 'FORMAT/DP[0]<5 | FORMAT/DP[1]<5 | AF<0.1' \
  -Oz -o donor_filter.vcf.gz

# 2. reorder VCF contigs to match the BAM header (required by popscle), then index
./sort_vcf_same_as_bam.sh \
  path/to/cellranger/HK2385_HK2814/HK2385_HK2814_output/outs/gex_possorted_bam.bam \
  donor_filter.vcf.gz z > donor_filter.reorder.vcf.gz
tabix -p vcf donor_filter.reorder.vcf.gz

# 3. assign RNA barcodes to donors
demuxlet --alpha 0 --alpha 0.5 \
  --group-list path/to/cellbender/HK2385_HK2814_cell_barcodes.csv \
  --field GT \
  --sam path/to/cellranger/HK2385_HK2814/HK2385_HK2814_output/outs/gex_possorted_bam.bam \
  --vcf donor_filter.reorder.vcf.gz \
  --out results
