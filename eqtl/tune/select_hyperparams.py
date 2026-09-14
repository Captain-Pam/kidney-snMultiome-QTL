#!/usr/bin/env python
"""Select optimal hyperparameters from the tuning sweeps.

Optimum = the value maximizing the number of eGenes (permutation p < 0.05).

Inputs:
  tune/tune_pc/cis_df_pc_<n>.csv          genotype-PC sweep (one cell type)
  tune/tune_peer/<ct>_peer<k>.csv         PEER sweep (per cell type)
Output:
  tune/peer_optimal.csv                   cell type -> chosen number of PEER factors

Run from the eqtl/ directory.
"""
import glob
import os

import pandas as pd


def n_egenes(path):
    return (pd.read_csv(path, index_col=0)['pval_perm'] < 0.05).sum()


# genotype PCs (single cell type) -> reported only
pcs = [0, 1, 2, 3, 4, 5, 10, 15, 20]
pc_counts = pd.Series(
    {n: n_egenes('tune/tune_pc/cis_df_pc_%d.csv' % n) for n in pcs}
)
print('genotype PC sweep (eGenes):')
print(pc_counts)
print('optimal n_pc =', pc_counts.idxmax())

# PEER factors (per cell type)
peers = [5, 10, 15, 20]
cts = sorted(os.path.basename(p).replace('_peer5.csv', '')
             for p in glob.glob('tune/tune_peer/*_peer5.csv'))

choice = {}
for ct in cts:
    counts = pd.Series({k: n_egenes('tune/tune_peer/%s_peer%d.csv' % (ct, k)) for k in peers})
    choice[ct] = counts.idxmax()

pd.Series(choice, name='n_peer').to_frame().to_csv('tune/peer_optimal.csv')
print('\nchosen PEER factors:')
print(pd.Series(choice))
