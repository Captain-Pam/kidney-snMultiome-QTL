"""
Combine per-caller demultiplexing results into a single consensus assignment.

For every pooled library (2 donors per lane) this collects the singlet /
doublet calls from all five callers:

    DMXT_*   demuxlet          (RNA)
    SPrna_*  souporcell        (RNA)
    VRrna_*  vireo             (RNA)
    SPatac_* souporcell        (ATAC)
    VRatac_* vireo             (ATAC)

Each caller contributes four columns: Sing (assigned donor), Doub
(SNG/DBL/UNK), Sing_LLK (donor1 vs donor2 log-likelihood) and Doub_LLK.
The consensus singlet / doublet label is the per-cell majority vote across the
five callers. ATAC-only barcodes (no RNA call) fall back to the vireo-ATAC call.

Two libraries (HK2524_HK3043, HK3039_HK3007) turned out to contain a single
donor and are assigned directly.

Paths below point at the original run layout; edit them for your environment.
"""

import os
import re
import glob

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# per-caller output directories (each holds outputs/<sample>/...)
# ---------------------------------------------------------------------------
demuxlet_dir = 'path/to/demuxlet'
souporcell_dir = 'path/to/souporcell_gt'
vireo_dir = 'path/to/vireo_gt'
souporcell_atac_dir = 'path/to/souporcell_atac'
vireo_atac_dir = 'path/to/vireo_atac'

cellbender_bc = 'path/to/cellbender/%s_cell_barcodes.csv'
atac_bc = 'path/to/atac_barcodes/%s.csv'

all_samples = sorted(re.sub('.*/', '', i)
                     for i in glob.glob('%s/outputs/*' % demuxlet_dir))


def donor_pair(sample):
    """Return the two donor IDs for a pooled library (with hand-fixed swaps)."""
    if sample == 'HK2518_HK2763':
        return 'HK2518', 'HK2762'
    if sample == 'HK3106_HK2598':
        return 'HK3106', 'HK2648'
    return tuple(sample.split('_'))


# ---------------------------------------------------------------------------
# 1. harmonize each caller into Sing / Doub / Sing_LLK / Doub_LLK columns
# ---------------------------------------------------------------------------
for samples in all_samples:
    g1, g2 = donor_pair(samples)
    os.makedirs('outputs/%s' % samples, exist_ok=True)

    # --- demuxlet (RNA) ---
    demuxlet = pd.read_csv('%s/outputs/%s/results.best' % (demuxlet_dir, samples),
                           sep='\t', index_col=0)
    demuxlet['DBL'] = np.where(demuxlet['BEST'].str.contains('DBL'), 'DBL', 'SNG')
    a = pd.concat([demuxlet['SNG.LLK2'][demuxlet['SNG.2ND'] == g1],
                   demuxlet['SNG.LLK1'][demuxlet['SNG.1ST'] == g1]])
    b = pd.concat([demuxlet['SNG.LLK2'][demuxlet['SNG.2ND'] == g2],
                   demuxlet['SNG.LLK1'][demuxlet['SNG.1ST'] == g2]])
    demuxlet['SNG_LLK'] = a[demuxlet.index] - b[demuxlet.index]
    demuxlet = demuxlet.loc[:, ['SNG.1ST', 'DBL', 'SNG_LLK', 'LLK12']]
    demuxlet.columns = 'DMXT_' + pd.Index(['Sing', 'Doub', 'Sing_LLK', 'Doub_LLK'])

    # --- souporcell (RNA) ---
    sp_rna = pd.read_csv('%s/outputs/%s/output/clusters.tsv' % (souporcell_dir, samples),
                         sep='\t', index_col=0)
    sp_rna['Sing'] = np.where(sp_rna['cluster0'] - sp_rna['cluster1'] > 0, g1, g2)
    sp_rna['Doub'] = (sp_rna['status']
                      .str.replace('singlet', 'SNG')
                      .str.replace('doublet', 'DBL')
                      .str.replace('unassigned', 'UNK'))
    sp_rna['Sing_LLK'] = sp_rna['cluster0'] - sp_rna['cluster1']
    sp_rna['Doub_LLK'] = sp_rna['log_prob_doublet']
    sp_rna = sp_rna.loc[:, ['Sing', 'Doub', 'Sing_LLK', 'Doub_LLK']]
    sp_rna.columns = 'SPrna_' + sp_rna.columns

    # --- vireo (RNA) ---
    vir_rna1 = pd.read_csv('%s/outputs/%s/vireo_out/prob_singlet.tsv.gz' % (vireo_dir, samples),
                           sep='\t', index_col=0)
    vir_rna = pd.read_csv('%s/outputs/%s/vireo_out/donor_ids.tsv' % (vireo_dir, samples),
                          sep='\t', index_col=0)
    vir_rna['SNG_LLK'] = vir_rna1[g1] - vir_rna1[g2]
    vir_rna['DBL'] = 'SNG'
    vir_rna.loc[vir_rna['donor_id'] == 'unassigned', 'DBL'] = 'UNK'
    vir_rna.loc[vir_rna['donor_id'] == 'doublet', 'DBL'] = 'DBL'
    vir_rna = vir_rna.loc[:, ['best_singlet', 'DBL', 'SNG_LLK', 'prob_doublet']]
    vir_rna.columns = ['Sing', 'Doub', 'Sing_LLK', 'Doub_LLK']
    vir_rna['Doub_LLK'] = np.log10(vir_rna['Doub_LLK'])
    vir_rna.columns = 'VRrna_' + vir_rna.columns

    rna_genotype = pd.concat([demuxlet,
                              sp_rna.loc[demuxlet.index, :],
                              vir_rna.loc[demuxlet.index, :]], axis=1)

    # --- souporcell (ATAC) ---
    sp_atac = pd.read_csv('%s/outputs/%s/output/clusters.tsv' % (souporcell_atac_dir, samples),
                          sep='\t', index_col=0)
    sp_atac['Sing'] = np.where(sp_atac['cluster0'] - sp_atac['cluster1'] > 0, g1, g2)
    sp_atac['Doub'] = (sp_atac['status']
                       .str.replace('singlet', 'SNG')
                       .str.replace('doublet', 'DBL')
                       .str.replace('unassigned', 'UNK'))
    sp_atac['Sing_LLK'] = sp_atac['cluster0'] - sp_atac['cluster1']
    sp_atac['Doub_LLK'] = sp_atac['log_prob_doublet']
    sp_atac = sp_atac.loc[:, ['Sing', 'Doub', 'Sing_LLK', 'Doub_LLK']]
    sp_atac.columns = 'SPatac_' + sp_atac.columns

    # --- vireo (ATAC) ---
    vir_atac1 = pd.read_csv('%s/outputs/%s/vireo_out/prob_singlet.tsv.gz' % (vireo_atac_dir, samples),
                            sep='\t', index_col=0)
    vir_atac = pd.read_csv('%s/outputs/%s/vireo_out/donor_ids.tsv' % (vireo_atac_dir, samples),
                           sep='\t', index_col=0)
    vir_atac['SNG_LLK'] = vir_atac1[g1] - vir_atac1[g2]
    vir_atac['DBL'] = 'SNG'
    vir_atac.loc[vir_atac['donor_id'] == 'unassigned', 'DBL'] = 'UNK'
    vir_atac.loc[vir_atac['donor_id'] == 'doublet', 'DBL'] = 'DBL'
    vir_atac = vir_atac.loc[:, ['best_singlet', 'DBL', 'SNG_LLK', 'prob_doublet']]
    vir_atac.columns = ['Sing', 'Doub', 'Sing_LLK', 'Doub_LLK']
    vir_atac['Doub_LLK'] = np.log10(vir_atac['Doub_LLK'])
    vir_atac.columns = 'VRatac_' + vir_atac.columns

    atac_genotype = pd.concat([sp_atac, vir_atac.loc[sp_atac.index, :]], axis=1)

    all_genotype = pd.concat([rna_genotype, atac_genotype], axis=1)
    rna_genotype.to_csv('outputs/%s/rna.csv' % samples)
    atac_genotype.to_csv('outputs/%s/atac.csv' % samples)
    all_genotype.to_csv('outputs/%s/all.csv' % samples)
    print(samples)


# ---------------------------------------------------------------------------
# 2. majority-vote consensus across the five callers
# ---------------------------------------------------------------------------
for i in all_samples:
    a = pd.read_csv('outputs/%s/all.csv' % i, index_col=0)
    out = pd.DataFrame(np.nan, index=a.index, columns=['singlet', 'doublet'])

    sing = a.loc[:, ['DMXT_Sing', 'SPrna_Sing', 'VRrna_Sing', 'SPatac_Sing', 'VRatac_Sing']]
    out['singlet'] = sing.mode(axis=1).iloc[:, 0]

    doub = a.loc[:, ['DMXT_Doub', 'SPrna_Doub', 'VRrna_Doub', 'SPatac_Doub', 'VRatac_Doub']]
    out['doublet'] = doub.mode(axis=1).iloc[:, 0]

    # ATAC-only barcodes have no RNA call -> use vireo-ATAC
    atac_only = a['SPrna_Sing'].isna()
    out.loc[atac_only, 'singlet'] = a.loc[atac_only, 'VRatac_Sing']
    out.loc[atac_only, 'doublet'] = a.loc[atac_only, 'VRatac_Doub']

    out.to_csv('outputs/%s/final.csv' % i)
    print(i)


# ---------------------------------------------------------------------------
# 3. single-donor libraries assigned directly
# ---------------------------------------------------------------------------
for i, donor in [('HK2524_HK3043', 'HK2524'), ('HK3039_HK3007', 'HK3039')]:
    bc_rna = pd.read_csv(cellbender_bc % i, index_col=0, header=None)
    bc_atac = pd.read_csv(atac_bc % i, index_col=0, header=None)
    all_bc = bc_rna.index.union(bc_atac.index)
    all_bc.name = None
    out = pd.DataFrame(np.nan, index=all_bc, columns=['singlet', 'doublet'])
    out['singlet'] = donor
    out['doublet'] = 'SNG'
    os.makedirs('outputs/%s' % i, exist_ok=True)
    out.to_csv('outputs/%s/final.csv' % i)


# ---------------------------------------------------------------------------
# 4. concatenate all libraries into one barcode -> donor table
# ---------------------------------------------------------------------------
samples = [re.sub('.*/', '', i) for i in glob.glob('outputs/*')]
out = []
for i in samples:
    a = pd.read_csv('outputs/%s/final.csv' % i, index_col=0)
    a.index = a.index.astype(str) + '-' + i
    out += [a]
pd.concat(out).to_csv('all_genotype.csv')
