#!/bin/bash
# Step 5 - PEER factors.
#
# Computes k = 5, 10, 15, 20 PEER factors per cell type (one Slurm job each);
# the tuning stage later picks the optimal k.
#
# Input:   prepare_data/bed/<ct>.expression.bed.gz
# Output:  prepare_data/peer/<ct>_peer<k>.PEER_covariates.txt
#
# Requires the GTEx qtl pipeline run_PEER.R and a PEER R environment.
# Run from the eqtl/ directory.

gtex=gtex-pipeline/qtl/src   # path to the GTEx pipeline scripts
mkdir -p prepare_data/peer

for ct in CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT; do
    sbatch --job-name=peer_${ct} --output=prepare_data/peer/${ct}.out \
           --time=24:00:00 --mem=32G \
           --wrap="source ~/.bashrc && conda activate peer && \
                   for k in 5 10 15 20; do \
                       Rscript ${gtex}/run_PEER.R prepare_data/bed/${ct}.expression.bed.gz \
                               prepare_data/peer/${ct}_peer\${k} \${k}; \
                   done"
done
