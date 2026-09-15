#!/usr/bin/env python
"""Step 1 - Pseudobulk expression from the cleaned snMultiome object.

Reads the QC'd RNA AnnData and, for every (cell type x sample), sums raw counts
into a pseudobulk matrix and counts the number of cells.

Input:
  data/rna_clean.h5ad   (obs: multivi_cell_type, demultiplex_sample)
Outputs:
  prepare_data/anndata/counts/<ct>.csv   genes x samples, summed raw counts
  prepare_data/anndata/cellcounts.csv    samples x cell types, number of cells

Run from the eqtl/ directory.
"""
import os

import anndata
import numpy as np
import pandas as pd

ad = anndata.read_h5ad('data/rna_clean.h5ad')

cts = ad.obs['multivi_cell_type'].cat.categories
donors = sorted(ad.obs['demultiplex_sample'].unique())

os.makedirs('prepare_data/anndata/counts', exist_ok=True)

cellcount = pd.DataFrame(0, index=donors, columns=cts)
for ct in cts:
    ad_ct = ad[ad.obs['multivi_cell_type'] == ct, :]
    n = ad_ct.obs['demultiplex_sample'].value_counts()
    cellcount.loc[n.index, ct] = n

    samples = sorted(ad_ct.obs['demultiplex_sample'].unique())
    out = pd.DataFrame(np.nan, index=ad_ct.var_names, columns=samples)
    for s in samples:
        tmp = ad_ct[ad_ct.obs['demultiplex_sample'] == s, :]
        out.loc[:, s] = np.array(tmp.X.sum(axis=0)).flatten()
    out.to_csv('prepare_data/anndata/counts/%s.csv' % ct)
    print(ct)

cellcount.to_csv('prepare_data/anndata/cellcounts.csv')
