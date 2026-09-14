#!/usr/bin/env python
"""cis-eQTL permutation pass (tensorQTL map_cis) + genome-wide FDR.

Outputs one CSV: per-gene top association, empirical permutation p-value
(`pval_perm`), nominal significance threshold (`pval_nominal_threshold`) and
Storey q-value (`qval`). Used both for hyperparameter tuning and the final run.
"""
import argparse

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, post


def load_covariates(cov_files, samples):
    dfs = []
    for f in cov_files.split(','):
        c = pd.read_csv(f, sep='\t', index_col=0).T
        dfs.append(c.loc[samples, :])
    return pd.concat(dfs, axis=1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geno", required=True, help="plink genotype prefix")
    p.add_argument("--rna", required=True, help="normalized expression BED")
    p.add_argument("--cov", required=True, help="comma-separated covariate files")
    p.add_argument("--out", required=True, help="output CSV")
    p.add_argument("--window", type=int, default=250000, help="cis window (+/- bp)")
    p.add_argument("--maf", type=float, default=0.1, help="minor allele freq threshold")
    p.add_argument("--fdr", type=float, default=0.1, help="genome-wide FDR")
    args = p.parse_args()

    print("using device:", "cuda" if torch.cuda.is_available() else "cpu")

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.rna)

    pr = genotypeio.PlinkReader(args.geno)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    variant_df['chrom'] = 'chr' + variant_df['chrom'].astype(str)
    variant_df['chrom'] = variant_df['chrom'].str.replace('chr23', 'chrX')

    cov_df = load_covariates(args.cov, phenotype_df.columns)

    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                         cov_df, window=args.window, maf_threshold=args.maf)
    post.calculate_qvalues(cis_df, fdr=args.fdr, qvalue_lambda=0.85)
    cis_df.to_csv(args.out)


if __name__ == "__main__":
    main()
