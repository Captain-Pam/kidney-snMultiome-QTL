#!/usr/bin/env python
"""Cache predicted reference/alternate coverage around variants.

The example-locus figures draw the predicted coverage profile itself rather than
a single summary score, so the profile is computed once here and written to an
npz the plotting scripts read. That keeps the figures runnable on a machine with
no GPU and no model checkpoint.

Output, one file per variant: <out_dir>/<chrom>_<pos>_<ref>_<alt>.npz
    pred_ref     (bins, tracks)  predicted coverage, reference allele
    pred_alt     (bins, tracks)  predicted coverage, alternate allele
    bin_centers  (bins,)         genomic coordinate of each bin
    track_descriptions           targets `description` per track, for labelling

Usage:
    python 2_predict_tracks.py --variants variants.tsv --out_dir <dir> \
        [--targets <targets.txt>]

`--variants` is a TSV with columns chrom, pos, ref, alt.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cerberus_model import CerberusFold  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument(
        "--targets",
        default=None,
        help="targets file; defaults to the full human set. Pass the ATAC subset "
        "to keep the cached arrays small when only accessibility is plotted.",
    )
    parser.add_argument("--fold", type=int, default=0)
    args = parser.parse_args()

    targets_path = args.targets or os.environ["TARGETS_HUMAN"]
    fold_dir = os.path.join(
        os.environ["CERBERUS_MODEL_DIR"], f"f{args.fold}c0", "train"
    )

    fold = CerberusFold(
        os.path.join(fold_dir, "params.json"),
        os.path.join(fold_dir, "model_best.pth"),
        targets_path,
    )
    descriptions = np.array(fold.targets["description"].tolist(), dtype=object)

    fasta = pysam.FastaFile(os.environ["CERBERUS_FASTA"])
    os.makedirs(args.out_dir, exist_ok=True)

    variants = pd.read_csv(args.variants, sep="\t")
    for _, variant in variants.iterrows():
        key = f"{variant.chrom}_{int(variant.pos)}_{variant.ref}_{variant.alt}"
        out_path = os.path.join(args.out_dir, f"{key}.npz")
        if os.path.exists(out_path):
            print(f"  [skip] {out_path}")
            continue

        print(f"  predicting {key}")
        pred_ref, pred_alt, centers = fold.predict_alleles(
            fasta, variant.chrom, int(variant.pos), variant.ref, variant.alt
        )

        np.savez_compressed(
            out_path,
            pred_ref=pred_ref.astype(np.float32),
            pred_alt=pred_alt.astype(np.float32),
            bin_centers=centers,
            track_descriptions=descriptions,
        )

    fasta.close()
    print(f"done -> {args.out_dir}")


if __name__ == "__main__":
    main()
