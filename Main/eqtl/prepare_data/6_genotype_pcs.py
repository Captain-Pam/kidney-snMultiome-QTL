#!/usr/bin/env python
"""Step 6 - Genotype principal components (for tuning).

Splits a precomputed genotype PCA eigenvector file into covariate tables with
1..20 PCs. Genotype PCs are only evaluated during tuning; the final run uses
none (PC = 0 was optimal).

The eigenvector file is produced beforehand with plink, e.g.:
    plink2 --bfile data/genotype/final_QCed --pca 20 \
           --out data/genotype/genotype_pca

Input:   data/genotype/genotype_pca.eigenvec   (FID IID PC1 ... PC20)
Outputs: prepare_data/pca/PC_<n>.txt           (PCs x samples), n in 1..20

Run from the eqtl/ directory.
"""
import os

import pandas as pd

os.makedirs('prepare_data/pca', exist_ok=True)

ev = pd.read_csv('data/genotype/genotype_pca.eigenvec', sep=r'\s+', header=None)
ev.index = ev.iloc[:, 1]          # IID
ev = ev.iloc[:, 2:].T             # PCs x samples
ev.index = ['PC%d' % i for i in range(1, ev.shape[0] + 1)]

for n in [1, 2, 3, 4, 5, 10, 15, 20]:
    ev.iloc[:n, :].to_csv('prepare_data/pca/PC_%d.txt' % n, sep='\t')
