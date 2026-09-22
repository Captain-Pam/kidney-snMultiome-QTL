#!/bin/bash
# Step 4 - Normalized expression BEDs (GTEx TMM pipeline).
#
# Applies GTEx expression thresholds and TMM normalization to produce the
# phenotype BED consumed by tensorQTL.
#
# Inputs (per cell type, from step 3):
#   prepare_data/sample_qc/<ct>/{tpm.gct, raw.gct, sample_lookup.txt}
#   prepare_data/sample_qc/gencode_v46_rev.gtf
# Output:
#   prepare_data/bed/<ct>.expression.bed.gz
#
# Requires the GTEx qtl pipeline (https://github.com/broadinstitute/gtex-pipeline).
# Run from the eqtl/ directory.

qc_dir=prepare_data/sample_qc
gtf=${qc_dir}/gencode_v46_rev.gtf
gtex=gtex-pipeline/qtl/src   # path to the GTEx pipeline scripts

mkdir -p prepare_data/bed

for ct in CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT; do
    python ${gtex}/eqtl_prepare_expression.py \
        ${qc_dir}/${ct}/tpm.gct \
        ${qc_dir}/${ct}/raw.gct \
        ${gtf} \
        ${qc_dir}/${ct}/sample_lookup.txt \
        prepare_data/bed/${ct} \
        --tpm_threshold 0.1 \
        --count_threshold 6 \
        --sample_frac_threshold 0.2 \
        --normalization_method tmm > prepare_data/bed/${ct}.log 2>&1
    echo ${ct}
done
