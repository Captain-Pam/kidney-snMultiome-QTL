#!/usr/bin/env python
"""Assemble the final cis-eQTL tables.

Concatenates the per-chromosome nominal pairs, filters to common variants in
the cis window, adds allele info and the gene-level significance from the
permutation pass.

Inputs (per cell type, from run/run_eqtl.sh):
  run/outputs/<ct>.cis_qtl_pairs.chr<1..22>.parquet   nominal pairs
  run/outputs/<ct>_cis_df.csv                          gene-level q-values
  data/genotype/final_QCed.bim                         REF/ALT alleles
Outputs:
  run/all_parquets/<ct>.parquet   all filtered pairs (+ egene_qval, thresholds)
  run/fdr/<ct>.parquet            significant subset (egene_qval < 0.1)

Run from the eqtl/ directory.
"""
import glob
import os

import pandas as pd

bim = pd.read_csv('data/genotype/final_QCed.bim', header=None, sep='\t')
bim.columns = ['V%d' % (i + 1) for i in range(6)]
bim.index = bim['V2']

cts = sorted(os.path.basename(p).replace('_cis_df.csv', '')
             for p in glob.glob('run/outputs/*_cis_df.csv'))

os.makedirs('run/all_parquets', exist_ok=True)
os.makedirs('run/fdr', exist_ok=True)

for ct in cts:
    out = pd.concat(
        [pd.read_parquet('run/outputs/%s.cis_qtl_pairs.chr%d.parquet' % (ct, i)) for i in range(1, 23)],
        axis=0,
    )
    out = out.loc[(out['af'] >= 0.1) & (out['start_distance'].abs() < 250000), :]

    out['A1'] = bim.loc[out['variant_id'], 'V5'].values
    out['A2'] = bim.loc[out['variant_id'], 'V6'].values

    map_cis = pd.read_csv('run/outputs/%s_cis_df.csv' % ct, index_col=0)
    out['egene_qval'] = map_cis.loc[out['phenotype_id'], 'qval'].values
    out['pval_nominal_threshold'] = map_cis.loc[out['phenotype_id'], 'pval_nominal_threshold'].values

    out.to_parquet('run/all_parquets/%s.parquet' % ct)
    out.loc[out['egene_qval'] < 0.1, :].to_parquet('run/fdr/%s.parquet' % ct)
    print(ct)
