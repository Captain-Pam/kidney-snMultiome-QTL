#!/bin/bash
# Final cis-eQTL run.
#
# Uses a fixed 5 PEER factors for every cell type (the value chosen after
# tuning) and no genotype PCs. Runs two passes with the same covariates:
#   - map_cis.py      permutation pass  (window +/-250 kb, MAF >= 0.1) -> q-values
#   - map_nominal.py  nominal pass      (window +/-500 kb, all pairs)
# After the jobs finish, run:  python run/assemble_parquet.py
#
# Run from the eqtl/ directory.
conda activate tensorqtl

geno=data/genotype/final_QCed
prep=prepare_data
n_peer=5
gpu="--partition=gpu --gres=gpu:1 -c 4 --mem=23000 --time=3-00:00:00"

mkdir -p run/outputs

for ct in CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT; do
    rna=${prep}/bed/${ct}.expression.bed.gz
    cov=${prep}/peer/${ct}_peer${n_peer}.PEER_covariates.txt,${prep}/covs/covs_${ct}.txt

    sbatch ${gpu} --job-name=cis_${ct} --output=run/outputs/${ct}.cis.out \
        --wrap="source ~/.bashrc && conda activate tensorqtl && \
                python src/map_cis.py --geno ${geno} --rna ${rna} --cov ${cov} \
                    --out run/outputs/${ct}_cis_df.csv"

    sbatch ${gpu} --job-name=nom_${ct} --output=run/outputs/${ct}.nominal.out \
        --wrap="source ~/.bashrc && conda activate tensorqtl && \
                python src/map_nominal.py --geno ${geno} --rna ${rna} --cov ${cov} \
                    --out run/outputs/${ct}"
done
