#!/usr/bin/env python
"""Collect ChromBPNet held-out accuracy into two small tables for plotting.

Reads each trained run's predictions and its own coverage bigWig, so the
observed counts are exactly what the model was trained against. No GPU and no
model loading needed.

    pearson_r.csv     cell_type, fold, pearsonr, n_peaks
    scatter_data.csv  cell_type, fold, observed, predicted   (examples only)

Peaks with log1p(observed counts) at or below `--log_count_threshold` are
excluded. Near-empty peaks carry almost no signal to predict, and including
them inflates the correlation by stretching the range.

Usage:
    python 5_export_metrics.py --out_dir <dir> [--examples TAL:0]
"""
import argparse
import csv
import os

import h5py
import numpy as np
import pyBigWig
from scipy.stats import pearsonr

CELL_TYPES = [
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT",
]

# ChromBPNet predicts total counts over the central 1 kb of its input.
WINDOW_HALF = 500


def load_counts(model_dir, celltype, fold, threshold):
    """Observed and predicted log counts at this fold's held-out peaks."""
    run_dir = os.path.join(model_dir, f"{celltype}_fold_{fold}")
    pred_h5 = os.path.join(run_dir, "evaluation", "chrombpnet_predictions.h5")
    bigwig = os.path.join(run_dir, "auxiliary", "data_unstranded.bw")

    with h5py.File(pred_h5, "r") as handle:
        chroms = handle["coords/coords_chrom"][:].astype(str)
        centers = handle["coords/coords_center"][:]
        predicted = handle["predictions/logcounts"][:]

    coverage = pyBigWig.open(bigwig)
    chrom_sizes = dict(coverage.chroms())
    observed = np.zeros(len(centers), dtype=np.float64)
    for i, (chrom, center) in enumerate(zip(chroms, centers)):
        start = max(0, int(center) - WINDOW_HALF)
        end = min(chrom_sizes.get(chrom, int(center) + WINDOW_HALF),
                  int(center) + WINDOW_HALF)
        observed[i] = np.nansum(coverage.values(chrom, start, end, numpy=True))
    coverage.close()

    observed = np.log1p(observed)
    keep = observed > threshold
    return observed[keep], predicted[keep]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument(
        "--examples",
        nargs="+",
        default=["TAL:0"],
        help="cell_type:fold pairs to export per-peak values for",
    )
    parser.add_argument(
        "--n_folds", type=int, default=int(os.environ.get("CHROMBPNET_N_FOLDS", 5))
    )
    parser.add_argument("--log_count_threshold", type=float, default=2.0)
    args = parser.parse_args()

    model_dir = os.environ["CHROMBPNET_MODEL_DIR"]
    os.makedirs(args.out_dir, exist_ok=True)

    examples = {}
    for item in args.examples:
        celltype, fold = item.split(":")
        examples.setdefault(celltype, set()).add(int(fold))

    pearson_rows = []
    scatter_rows = []

    for celltype in CELL_TYPES:
        for fold in range(args.n_folds):
            run_dir = os.path.join(model_dir, f"{celltype}_fold_{fold}")
            if not os.path.isdir(run_dir):
                print(f"  missing, skipping: {run_dir}")
                continue

            observed, predicted = load_counts(
                model_dir, celltype, fold, args.log_count_threshold
            )
            r = pearsonr(observed, predicted)[0]
            pearson_rows.append(
                {
                    "cell_type": celltype,
                    "fold": fold,
                    "pearsonr": r,
                    "n_peaks": len(observed),
                }
            )

            if fold in examples.get(celltype, ()):
                for obs, pred in zip(observed, predicted):
                    scatter_rows.append(
                        {
                            "cell_type": celltype,
                            "fold": fold,
                            "observed": obs,
                            "predicted": pred,
                        }
                    )

        folds = [row["pearsonr"] for row in pearson_rows if row["cell_type"] == celltype]
        if folds:
            print(f"  {celltype:12s} mean r = {np.mean(folds):.4f}")

    pearson_path = os.path.join(args.out_dir, "pearson_r.csv")
    with open(pearson_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["cell_type", "fold", "pearsonr", "n_peaks"]
        )
        writer.writeheader()
        writer.writerows(pearson_rows)
    print(f"{len(pearson_rows)} rows -> {pearson_path}")

    scatter_path = os.path.join(args.out_dir, "scatter_data.csv")
    with open(scatter_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["cell_type", "fold", "observed", "predicted"]
        )
        writer.writeheader()
        writer.writerows(scatter_rows)
    print(f"{len(scatter_rows)} rows -> {scatter_path}")


if __name__ == "__main__":
    main()
