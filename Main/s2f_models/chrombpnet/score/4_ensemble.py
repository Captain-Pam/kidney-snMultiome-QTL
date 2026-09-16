#!/usr/bin/env python
"""Average ChromBPNet variant scores across folds.

Every numeric column is averaged over the folds; the identifying columns are
taken from fold 0. The folds differ only in their train/validation/test
chromosome split, so averaging them is the model-ensembling step, not an
aggregation over different variants — hence the check that all folds scored the
same variants in the same order.

Usage:
    python 4_ensemble.py --scores_folds_dir <...>/scores_folds \
                         --out_dir <...>/scores_ensemble
"""
import argparse
import os

import numpy as np
import pandas as pd

CELL_TYPES = [
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT",
]
META_COLS = ["chr", "pos", "variant_id", "allele1", "allele2"]

# variant-scorer writes the real variants and the shuffled background to
# separate files with the same schema.
SUFFIXES = ["variant_scores.tsv", "variant_scores.shuffled.tsv"]


def ensemble_celltype(celltype, scores_folds_dir, out_dir, n_folds, suffix):
    frames = []
    for fold in range(n_folds):
        path = os.path.join(
            scores_folds_dir, celltype, f"fold_{fold}",
            f"{celltype}_fold_{fold}.{suffix}",
        )
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        frames.append(pd.read_csv(path, sep="\t"))

    reference = frames[0]
    for fold, frame in enumerate(frames[1:], start=1):
        if not (frame["variant_id"].values == reference["variant_id"].values).all():
            raise ValueError(
                f"{celltype} fold_{fold} scored a different variant order than fold_0"
            )

    numeric_cols = [c for c in reference.columns if c not in META_COLS]
    stacked = np.stack([frame[numeric_cols].values for frame in frames], axis=0)
    averaged = np.nanmean(stacked, axis=0)

    ensemble = reference[META_COLS].copy()
    for i, col in enumerate(numeric_cols):
        ensemble[col] = averaged[:, i]

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{celltype}.{suffix}")
    ensemble.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores_folds_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument(
        "--n_folds", type=int, default=int(os.environ.get("CHROMBPNET_N_FOLDS", 5))
    )
    parser.add_argument("--celltypes", nargs="+", default=CELL_TYPES)
    args = parser.parse_args()

    for celltype in args.celltypes:
        print(f"{celltype}:")
        ct_out = os.path.join(args.out_dir, celltype)
        for suffix in SUFFIXES:
            try:
                ensemble_celltype(
                    celltype, args.scores_folds_dir, ct_out, args.n_folds, suffix
                )
            except FileNotFoundError as exc:
                print(f"  skipping {suffix}: missing {exc}")


if __name__ == "__main__":
    main()
