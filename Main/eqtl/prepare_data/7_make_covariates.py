#!/usr/bin/env python
"""Step 7 - Build the covariate tables (one per cell type).

Covariates: Gender, Age (from metadata) plus per-cell-type technical terms
averaged across cells within each (cell type x sample):
  log_nCount_RNA = log2(mean nCount_RNA + 1)
  log_cell_count = log2(number of cells + 1)
  percent_mito, percent_ribo

Inputs:
  data/metadata/multiome_metadata.txt   (index = sample; cols Gender, Age)
  data/rna_clean.h5ad                    (obs: multivi_cell_type, demultiplex_sample,
                                          nCount_RNA, percent_mito, percent_ribo)
Outputs:
  prepare_data/covs/covs_<ct>.txt        covariates x samples

Run from the eqtl/ directory.
"""
import os

import anndata
import numpy as np
import pandas as pd

os.makedirs('prepare_data/covs', exist_ok=True)

meta = pd.read_csv('data/metadata/multiome_metadata.txt', sep='\t', index_col=0)
meta = meta.loc[:, ['Gender', 'Age']]

ad = anndata.read_h5ad('data/rna_clean.h5ad')
obs = ad.obs.loc[:, ['multivi_cell_type', 'demultiplex_sample',
                     'nCount_RNA', 'percent_mito', 'percent_ribo']]
means = obs.groupby(['multivi_cell_type', 'demultiplex_sample']).mean()
sizes = obs.groupby(['multivi_cell_type', 'demultiplex_sample']).size()
df = pd.concat([means, sizes.rename('cell_count')], axis=1).reset_index()

for ct in df['multivi_cell_type'].unique():
    d = df.loc[df['multivi_cell_type'] == ct, :].copy()
    d.index = d['demultiplex_sample']
    d['log_nCount_RNA'] = np.log2(d['nCount_RNA'] + 1)
    d['log_cell_count'] = np.log2(d['cell_count'] + 1)
    d = d.loc[:, ['log_nCount_RNA', 'log_cell_count', 'percent_mito', 'percent_ribo']]

    ovlp = meta.index.intersection(d.index)
    covs = pd.concat([meta.loc[ovlp, :], d.loc[ovlp, :]], axis=1)
    covs.T.to_csv('prepare_data/covs/covs_%s.txt' % ct, sep='\t')
    print(ct)
