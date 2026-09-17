"""
QC: pairwise agreement between demultiplexing callers.

For every pooled library, computes the fraction of cells on which pairs of
callers agree on the singlet donor and on the singlet/doublet status, over
cells that (a) have a call from every caller and (b) pass RNA metainfo QC.
Produces box plots of the agreement distributions across libraries.

Run after `combine.py` (reads its `outputs/<sample>/all.csv`).
"""

import glob
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rna_metainfo = 'path/to/rna_metainfo/%s.csv'

all_samples = sorted(re.sub('.*/', '', i) for i in glob.glob('outputs/*'))

pairs = ['DMXT_SPrna', 'DMXT_VRrna', 'SPrna_VRrna', 'SPrna_SPatac', 'VRrna_VRatac']
singlet = pd.DataFrame(np.nan, index=all_samples, columns=pairs)
doublet = pd.DataFrame(np.nan, index=all_samples, columns=pairs)


def agree(a, col1, col2):
    return (a[col1] == a[col2]).sum() / a.shape[0]


for samples in all_samples:
    rna = pd.read_csv(rna_metainfo % samples, index_col=0)
    a = pd.read_csv('outputs/%s/all.csv' % samples, index_col=0)
    a = a.loc[a.isna().sum(axis=1) == 0, :]
    a = a.loc[a.index.isin(rna.index), :]

    for c1, c2, name in [('DMXT', 'SPrna', 'DMXT_SPrna'),
                         ('DMXT', 'VRrna', 'DMXT_VRrna'),
                         ('SPrna', 'VRrna', 'SPrna_VRrna'),
                         ('SPrna', 'SPatac', 'SPrna_SPatac'),
                         ('VRrna', 'VRatac', 'VRrna_VRatac')]:
        singlet.loc[samples, name] = agree(a, '%s_Sing' % c1, '%s_Sing' % c2)
        doublet.loc[samples, name] = agree(a, '%s_Doub' % c1, '%s_Doub' % c2)


f, axs = plt.subplots(ncols=2, figsize=(10, 3))
singlet.plot.box(ax=axs[0])
axs[0].set_ylabel('fraction agreement')
axs[0].set_title('singlet prediction')
doublet.plot.box(ax=axs[1])
axs[1].set_ylabel('fraction agreement')
axs[1].set_title('doublet prediction')
f.tight_layout()
f.savefig('caller_agreement.png', dpi=150)
