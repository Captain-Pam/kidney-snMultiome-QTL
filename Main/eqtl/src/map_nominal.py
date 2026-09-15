#!/usr/bin/env python
"""cis-eQTL nominal pass (tensorQTL map_nominal): every variant-gene pair.

Writes one parquet per chromosome: <out>.cis_qtl_pairs.chr<N>.parquet.
Run with a wider 500 kb window and no MAF filter; the final tables are trimmed
to +/-250 kb and MAF >= 0.1 during assembly.
"""
import argparse

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geno", required=True, help="plink genotype prefix")
    p.add_argument("--rna", required=True, help="normalized expression BED")
    p.add_argument("--cov", required=True, help="comma-separated covariate files")
    p.add_argument("--out", required=True, help="output prefix (incl. directory)")
    p.add_argument("--window", type=int, default=500000, help="cis window (+/- bp)")
    p.add_argument("--maf", type=float, default=0.0, help="minor allele freq threshold")
    args = p.parse_args()

    print("using device:", "cuda" if torch.cuda.is_available() else "cpu")

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.rna)

    pr = genotypeio.PlinkReader(args.geno)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    variant_df['chrom'] = 'chr' + variant_df['chrom'].astype(str)
    variant_df['chrom'] = variant_df['chrom'].str.replace('chr23', 'chrX')

    dfs = []
    for f in args.cov.split(','):
        c = pd.read_csv(f, sep='\t', index_col=0).T
        dfs.append(c.loc[phenotype_df.columns, :])
    cov_df = pd.concat(dfs, axis=1)

    for i in range(1, 23):
        chrom = 'chr%d' % i
        cis.map_nominal(genotype_df, variant_df,
                        phenotype_df.loc[phenotype_pos_df['chr'] == chrom],
                        phenotype_pos_df.loc[phenotype_pos_df['chr'] == chrom],
                        args.out, covariates_df=cov_df,
                        window=args.window, maf_threshold=args.maf)
        print(chrom)


if __name__ == "__main__":
    main()
