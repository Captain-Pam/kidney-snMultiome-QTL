#!/usr/bin/env python
"""Build per-cell-type variant lists of fine-mapped eQTL SNPs for ChromBPNet.

ChromBPNet is accessibility-only, so it cannot score a variant's effect on a
gene directly. Scoring the eQTL variants with it anyway gives the accessibility
counterpart of the Cerberus expression score, which is what the supplementary
comparison against eQTL effect size uses.

Allele convention, matching ../../cerberus/score/eqtl/2_prepare_vcf.py:

    match_type 'direct' -> REF = A2 (hg38 reference), ALT = A1 (effect allele)
    match_type 'swap'   -> REF = A1,                  ALT = A2

Model scores are log(ALT/REF), so for 'swap' rows the eQTL slope is negated
(`slope_adj`) to put both measures in the same ALT-vs-REF orientation. Getting
this backwards silently flips the sign of the correlation for those variants.

Outputs, under ${WORK_DIR}/chrombpnet/eqtl/variants/:
    <ct>_variants.tsv   headerless chr, pos, variant_id, ref, alt
    eqtl_cs_meta.tsv     per-row metadata for the comparison plots
"""
import os

import pandas as pd

CHR_OK = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

# The correlation panels are restricted to confidently fine-mapped variants.
PIP_THRESH = 0.5


def main():
    eqtl_tsv = os.environ["EQTL_SUSIE_TSV"]
    out_dir = os.path.join(os.environ["WORK_DIR"], "chrombpnet", "eqtl", "variants")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Reading {eqtl_tsv}")
    df = pd.read_csv(
        eqtl_tsv,
        sep="\t",
        compression="gzip",
        usecols=[
            "celltype", "variant_id", "phenotype_id", "A1", "A2",
            "match_type", "slope", "variable_prob", "cs",
        ],
    )
    print(f"  {len(df):,} rows")

    # cs == -1 means the variant fell outside every credible set.
    df = df[df["cs"] != -1].copy()
    print(f"  {len(df):,} credible-set rows")

    split = df["variant_id"].str.split(":", expand=True)
    df["chrom"] = split[0]
    df["pos"] = split[1].astype(int)
    df = df[df["chrom"].isin(CHR_OK)].copy()

    is_swap = df["match_type"] == "swap"
    df["ref"] = df["A2"].where(~is_swap, df["A1"])
    df["alt"] = df["A1"].where(~is_swap, df["A2"])

    # ChromBPNet scores substitutions only.
    df = df[(df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)].copy()
    print(f"  {len(df):,} credible-set SNP rows")

    df["slope_adj"] = df["slope"].where(~is_swap.loc[df.index], -df["slope"])

    # Encode the alleles into the id so that multi-allelic positions stay
    # distinguishable when scores are joined back on.
    df["scorer_vid"] = (
        df["chrom"] + ":" + df["pos"].astype(str) + ":" + df["ref"] + ":" + df["alt"]
    )

    meta_cols = [
        "celltype", "variant_id", "scorer_vid", "phenotype_id", "chrom", "pos",
        "ref", "alt", "match_type", "slope", "slope_adj", "variable_prob", "cs",
    ]
    meta_path = os.path.join(out_dir, "eqtl_cs_meta.tsv")
    df[meta_cols].to_csv(meta_path, sep="\t", index=False)
    print(f"  wrote {meta_path} ({len(df):,} rows)")

    df = df[df["variable_prob"] >= PIP_THRESH].copy()
    print(f"\nPer-cell-type unique SNPs (PIP >= {PIP_THRESH}):")
    for celltype, group in df.groupby("celltype"):
        unique = (
            group.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
            .sort_values(["chrom", "pos"])
        )
        out_path = os.path.join(out_dir, f"{celltype}_variants.tsv")
        unique[["chrom", "pos", "scorer_vid", "ref", "alt"]].to_csv(
            out_path, sep="\t", header=False, index=False
        )
        print(f"  {celltype:12s} {len(unique):6,d}")

    total = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    print(f"\n  unique across all cell types: {len(total):,}")


if __name__ == "__main__":
    main()
