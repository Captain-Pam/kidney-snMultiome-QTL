bcftools annotate --rename-chrs chr_rename.txt path/to/wgs_donors.vcf.gz -Oz -o all_multiome.chr_rename.vcf.gz

tabix -p vcf all_multiome.chr_rename.vcf.gz
