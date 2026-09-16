#!/usr/bin/env python
"""Write a VCF of fine-mapped eQTL credible-set variants for scoring.

Keeps variants assigned to a SuSiE credible set (`cs != -1`), deduplicates and
sorts. The PIP >= 0.5 cut used in the concordance analysis is applied later, at
plotting time, so the scored set stays available for the PIP-threshold sweeps.

Allele convention. The eQTL table records A1 as the effect allele, which is not
always the hg38 reference; `match_type` says which way round it is:

    direct -> REF = A2, ALT = A1
    swap   -> REF = A1, ALT = A2

The VCF is written in hg38 REF/ALT orientation, so model scores are
log(ALT/REF). Comparisons against the eQTL slope must negate the slope for
`swap` rows to match — see ../../../figures/eqtl_effect_concordance.py.

Usage:
    python 2_prepare_vcf.py --out <work>/eqtl/eqtl_susie_cs.vcf
"""
import argparse
import os

import pandas as pd

CHR_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHR_ORDER["chrX"] = 23
CHR_ORDER["chrY"] = 24


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="output VCF")
    parser.add_argument("--eqtl_tsv", default=os.environ.get("EQTL_SUSIE_TSV"))
    parser.add_argument(
        "--min_pip",
        type=float,
        default=0.0,
        help="optional PIP floor; the published run scored every credible-set variant",
    )
    args = parser.parse_args()

    print(f"Reading {args.eqtl_tsv}")
    eqtl = pd.read_csv(
        args.eqtl_tsv,
        sep="\t",
        compression="gzip",
        usecols=[
            "celltype", "variant_id", "phenotype_id",
            "A1", "A2", "match_type", "variable_prob", "cs",
        ],
    )
    print(f"  {len(eqtl):,} rows")

    credible = eqtl[eqtl["cs"] != -1].copy()
    print(f"  {len(credible):,} in a credible set")

    if args.min_pip > 0:
        credible = credible[credible["variable_prob"] >= args.min_pip]
        print(f"  {len(credible):,} with PIP >= {args.min_pip}")

    split = credible["variant_id"].str.split(":", expand=True)
    credible["chrom"] = split[0]
    credible["pos"] = split[1].astype(int)
    credible = credible[credible["chrom"].isin(CHR_ORDER)]
    print(f"  {len(credible):,} on canonical chromosomes")

    is_swap = credible["match_type"] == "swap"
    credible["ref"] = credible["A2"].where(~is_swap, credible["A1"])
    credible["alt"] = credible["A1"].where(~is_swap, credible["A2"])

    # A variant can appear for several cell types and genes; it only needs
    # scoring once, and the gene pairing comes from the GTF at scoring time.
    unique = credible.drop_duplicates(subset=["chrom", "pos", "ref", "alt"]).copy()
    print(f"  {len(unique):,} unique variants")

    unique["_order"] = unique["chrom"].map(CHR_ORDER)
    unique = unique.sort_values(["_order", "pos"])

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=2_prepare_vcf.py\n")
        handle.write(
            "##INFO=<ID=GENE,Number=1,Type=String,"
            'Description="eQTL target gene (phenotype_id)">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in unique.iterrows():
            handle.write(
                f"{row['chrom']}\t{row['pos']}\t{row['variant_id']}\t"
                f"{row['ref']}\t{row['alt']}\t.\tPASS\t"
                f"GENE={row['phenotype_id']}\n"
            )

    print(f"\n{len(unique):,} variants -> {args.out}")
    print(f"For 1000-variant chunks, score with --array=0-{len(unique) // 1000}")


if __name__ == "__main__":
    main()
