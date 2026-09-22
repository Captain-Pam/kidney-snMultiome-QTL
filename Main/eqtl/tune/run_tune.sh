#!/bin/bash
# Tuning - choose the number of genotype PCs and PEER factors.
#
# Metric: number of eGenes at permutation p < 0.05 (from src/map_cis.py).
#   - genotype PCs: swept on a single representative cell type (PTS), PEER fixed at 5
#   - PEER factors: swept per cell type, no genotype PCs
#
# Outputs:
#   tune/tune_pc/cis_df_pc_<n>.csv
#   tune/tune_peer/<ct>_peer<k>.csv
# Then run:  python tune/select_hyperparams.py
#
# Run from the eqtl/ directory.
conda activate tensorqtl

geno=data/genotype/final_QCed
prep=prepare_data
gpu="--partition=gpu --gres=gpu:1 -c 4 --mem=23000 --time=1-00:00:00"

#############################
# tune genotype PCs (on PTS) #
#############################
ct=PTS
mkdir -p tune/tune_pc

# n_pc = 0 (no genotype PCs)
python src/map_cis.py --geno ${geno} \
    --rna ${prep}/bed/${ct}.expression.bed.gz \
    --cov ${prep}/peer/${ct}_peer5.PEER_covariates.txt,${prep}/covs/covs_${ct}.txt \
    --out tune/tune_pc/cis_df_pc_0.csv

for n_pc in 1 2 3 4 5 10 15 20; do
    python src/map_cis.py --geno ${geno} \
        --rna ${prep}/bed/${ct}.expression.bed.gz \
        --cov ${prep}/peer/${ct}_peer5.PEER_covariates.txt,${prep}/pca/PC_${n_pc}.txt,${prep}/covs/covs_${ct}.txt \
        --out tune/tune_pc/cis_df_pc_${n_pc}.csv
done

##############################
# tune PEER (per cell type) #
##############################
mkdir -p tune/tune_peer
for ct in CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT; do
    for n_peer in 5 10 15 20; do
        sbatch ${gpu} --job-name=peer_${ct}_${n_peer} --output=tune/tune_peer/${ct}_${n_peer}.out \
            --wrap="source ~/.bashrc && conda activate tensorqtl && \
                    python src/map_cis.py --geno ${geno} \
                        --rna ${prep}/bed/${ct}.expression.bed.gz \
                        --cov ${prep}/peer/${ct}_peer${n_peer}.PEER_covariates.txt,${prep}/covs/covs_${ct}.txt \
                        --out tune/tune_peer/${ct}_peer${n_peer}.csv"
    done
done
