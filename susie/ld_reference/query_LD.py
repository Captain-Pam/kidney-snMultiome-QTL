#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import os
import scipy.sparse as sparse

def load_ld_npz(ld_prefix):
    
    #load the SNPs metadata
    gz_file = '%s.gz'%(ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps['A1'] + '.' + df_ld_snps['A2']
        
    #load the LD matrix
    npz_file = '%s.npz'%(ld_prefix)
    if not os.path.exists(npz_file):
        npz_file = '%s.npz2'%(ld_prefix)
    try: 
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s'%(npz_file))

    #create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps


def main():
    parser = argparse.ArgumentParser(description="Query LD matrix for SNPs. #path: ")
    parser.add_argument("--ld_prefix", type=str, help="Prefix for LD matrix files.")
    parser.add_argument("--input_file", type=str, help="Path to input SNPs CSV file.", default='tmp_snps.csv')
    parser.add_argument("--output_file", type=str, help="Path to save the output LD matrix CSV file.", default='ld_m_tmp.csv')
    
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    ld_prefix = args.ld_prefix
    
    tmp_snps = pd.read_csv(input_file)
    df_R, df_ld_snps = load_ld_npz(ld_prefix)
    df_ld_snps.index = 'chr' + df_ld_snps.index.str.replace('.', ':')
    df_R.index = df_R.columns = df_ld_snps.index
    
    df_ld_snps = df_ld_snps.loc[~df_ld_snps['BP'].duplicated(),:] # remove multi-allelic site
    df_ld_snps = df_ld_snps.loc[df_ld_snps.index.isin(tmp_snps['hg19_id1']) | df_ld_snps.index.isin(tmp_snps['hg19_id2']), :]
    ld_m = df_R.loc[df_ld_snps.index, df_ld_snps.index]
    ld_m.to_csv(output_file)

if __name__ == "__main__":
    main()
