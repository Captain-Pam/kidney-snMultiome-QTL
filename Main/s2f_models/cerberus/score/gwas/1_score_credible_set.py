#!/usr/bin/env python
"""Score a GWAS credible set with Cerberus: accessibility and expression effects.

For each variant in a locus's credible set, reports two alternate-minus-
reference effects in the chosen cell type:

    atac_logSUM   log2 ratio of predicted ATAC coverage summed over a
                  +/-512 bp window around the variant
    rna_logSED    log2 ratio of predicted RNA+ coverage summed over the
                  target gene's exons

These are the same statistics `hound_snp` computes in bulk for the caQTL and
eQTL sets; doing it directly here keeps the score tied to one named gene and one
cell type, which is what a locus plot needs.

Usage:
    python 1_score_credible_set.py \
        --variants cs_variants.tsv --gene FGF5 \
        --atac_track 21 --rna_track 22 --out cs_borzoi_scores.tsv

`--variants` is a TSV with columns: SNP, chrom, pos, ref, alt, pip.
Track indices refer to the strand-collapsed output ordering; the script asserts
the track descriptions match `--celltype` before scoring, so a wrong index
fails loudly rather than producing plausible numbers for the wrong cell type.
"""
import argparse
import os
import sys

import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cerberus_model import (  # noqa: E402
    CerberusFold,
    exon_intervals,
    log_ratio,
    mask_intervals,
)

# Matches hound_snp --local_window 1024.
LOCAL_BP = 512


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", required=True, help="credible-set TSV")
    parser.add_argument("--gene", required=True, help="target gene symbol")
    parser.add_argument("--out", required=True)
    parser.add_argument("--celltype", default="PTS")
    parser.add_argument("--atac_track", type=int, required=True)
    parser.add_argument("--rna_track", type=int, required=True)
    parser.add_argument("--fold", type=int, default=0)
    args = parser.parse_args()

    model_dir = os.environ["CERBERUS_MODEL_DIR"]
    fold_dir = os.path.join(model_dir, f"f{args.fold}c0", "train")

    fold = CerberusFold(
        os.path.join(fold_dir, "params.json"),
        os.path.join(fold_dir, "model_best.pth"),
        os.environ["TARGETS_HUMAN"],
    )
    fold.check_track(args.atac_track, f"ATAC:{args.celltype}")
    fold.check_track(args.rna_track, f"RNA+:{args.celltype}")
    print(f"ATAC track: {fold.describe(args.atac_track)}")
    print(f"RNA  track: {fold.describe(args.rna_track)}")

    exons = exon_intervals(os.environ["EQTL_GTF"], args.gene)
    print(f"{args.gene}: {len(exons)} exons")

    fasta = pysam.FastaFile(os.environ["CERBERUS_FASTA"])
    variants = pd.read_csv(args.variants, sep="\t")

    rows = []
    for _, variant in variants.iterrows():
        ref_cov, alt_cov, centers = fold.predict_alleles(
            fasta, variant.chrom, int(variant.pos), variant.ref, variant.alt
        )

        local = (centers >= variant.pos - LOCAL_BP) & (centers <= variant.pos + LOCAL_BP)
        atac_ref = ref_cov[local, args.atac_track].sum()
        atac_alt = alt_cov[local, args.atac_track].sum()

        exonic = mask_intervals(centers, exons)
        rna_ref = ref_cov[exonic, args.rna_track].sum()
        rna_alt = alt_cov[exonic, args.rna_track].sum()

        atac_effect = log_ratio(atac_ref, atac_alt)
        rna_effect = log_ratio(rna_ref, rna_alt)

        print(
            f"{str(variant.SNP):14s} {variant.chrom}:{variant.pos} "
            f"{variant.ref}>{variant.alt} PIP={variant.pip:.3f}  "
            f"ATAC={atac_effect:+.4f}  RNA={rna_effect:+.4f}"
        )

        rows.append(
            {
                "SNP": variant.SNP,
                "chrom": variant.chrom,
                "pos": int(variant.pos),
                "ref": variant.ref,
                "alt": variant.alt,
                "pip": variant.pip,
                "celltype": args.celltype,
                "gene": args.gene,
                "atac_logSUM": atac_effect,
                "rna_logSED": rna_effect,
                "atac_ref": atac_ref,
                "atac_alt": atac_alt,
                "rna_ref": rna_ref,
                "rna_alt": rna_alt,
            }
        )

    fasta.close()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"\n{len(rows)} variants -> {args.out}")


if __name__ == "__main__":
    main()
