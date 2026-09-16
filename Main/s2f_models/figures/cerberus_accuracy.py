#!/usr/bin/env python
"""Cerberus held-out prediction accuracy per cell type.

Pearson r between predicted and observed coverage on each fold's own held-out
test split, for the tracks from this study's multiome dataset. RNA is the mean
of the plus- and minus-strand tracks; ATAC is a single track. Each box is the
eight folds, ordered by median ATAC r.

Only the diagonal of the evaluation grid is read — model f<I>c0 evaluated on
data fold I. Evaluating a model on a fold it trained on would report inflated
accuracy.

Usage:
    python cerberus_accuracy.py [--out_dir <dir>] [--dataset Susztak_CKD_multiome]
"""
import argparse
import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figlib import CT_LABELS, cfg, model_ct, setup_style

RNA_BOX, RNA_EDGE, RNA_DOT = "#AED6F1", "#2A6DB5", "#1A5276"
ATAC_BOX, ATAC_EDGE, ATAC_DOT = "#A9DFBF", "#1E8449", "#145A32"

CT_SPACING = 1.4   # distance between cell-type centres
OFFSET = 0.24      # RNA drawn left of centre, ATAC right
BOX_WIDTH = 0.45


def load_test_folds(model_dir, n_folds, dataset_index=0):
    """Concatenate each fold's accuracy on its own held-out split."""
    frames = []
    for fold in range(n_folds):
        path = os.path.join(
            model_dir, f"f{fold}c0", f"eval{dataset_index}", f"fold{fold}", "acc.txt"
        )
        if not os.path.exists(path):
            raise SystemExit(f"Missing evaluation output: {path}")
        frame = pd.read_csv(path, sep="\t")
        frame["model_fold"] = fold
        frames.append(frame)

    accuracy = pd.concat(frames, ignore_index=True)
    # `description` is "<assay>:<cell_type>:<dataset>"; see data/3_make_targets.py.
    accuracy[["modality", "celltype", "dataset"]] = accuracy["description"].str.split(
        ":", n=2, expand=True
    )
    return accuracy


def draw_boxes(ax, data, positions, facecolor, edgecolor, dotcolor, rng):
    ax.boxplot(
        data,
        positions=positions,
        widths=BOX_WIDTH,
        patch_artist=True,
        medianprops=dict(color="black", linewidth=1.0),
        boxprops=dict(facecolor=facecolor, edgecolor=edgecolor, linewidth=0.5),
        whiskerprops=dict(color=edgecolor, linewidth=0.5),
        capprops=dict(color=edgecolor, linewidth=0.5),
        flierprops=dict(marker="", linestyle="none"),
        showfliers=False,
        zorder=2,
    )
    for position, values in zip(positions, data):
        jitter = rng.uniform(-0.10, 0.10, size=len(values))
        ax.scatter(
            position + jitter, values, s=5, color=dotcolor,
            alpha=0.85, zorder=3, linewidths=0,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", default=None)
    parser.add_argument("--dataset", default="Susztak_CKD_multiome")
    args = parser.parse_args()

    setup_style()
    out_dir = args.out_dir or cfg("FIGURE_DIR")
    os.makedirs(out_dir, exist_ok=True)

    n_folds = int(cfg("CERBERUS_N_FOLDS", "8"))
    accuracy = load_test_folds(cfg("CERBERUS_MODEL_DIR"), n_folds)

    subset = accuracy[accuracy["dataset"] == args.dataset].copy()
    if subset.empty:
        raise SystemExit(
            f"No tracks for dataset {args.dataset!r}. Available: "
            f"{sorted(accuracy['dataset'].unique())}"
        )
    subset["label"] = subset["celltype"].map(lambda c: CT_LABELS.get(model_ct(c), c))

    rna = (
        subset[subset["modality"].isin(["RNA+", "RNA-"])]
        .groupby(["label", "model_fold"])["pearsonr"]
        .mean()
        .reset_index()
    )
    rna_by_label = {label: group["pearsonr"].values for label, group in rna.groupby("label")}

    atac = subset[subset["modality"] == "ATAC"]
    atac_by_label = {
        label: group["pearsonr"].values for label, group in atac.groupby("label")
    }

    labels = sorted(
        atac_by_label,
        key=lambda label: np.median(atac_by_label[label]),
        reverse=True,
    )
    for label in labels:
        print(
            f"  {label:12s} ATAC median r={np.median(atac_by_label[label]):.3f}  "
            f"RNA median r={np.median(rna_by_label.get(label, [np.nan])):.3f}"
        )

    fig, ax = plt.subplots(figsize=(3.3, 0.8))

    rng = np.random.default_rng(42)
    centres = np.arange(len(labels)) * CT_SPACING
    rna_values = [rna_by_label[label] for label in labels]
    atac_values = [atac_by_label[label] for label in labels]

    draw_boxes(ax, rna_values, centres - OFFSET, RNA_BOX, RNA_EDGE, RNA_DOT, rng)
    draw_boxes(ax, atac_values, centres + OFFSET, ATAC_BOX, ATAC_EDGE, ATAC_DOT, rng)

    ax.set_xticks(centres)
    ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=5)
    ax.set_xlim(-0.7, centres[-1] + 0.7)

    everything = [v for values in rna_values + atac_values for v in values]
    ax.set_ylim(max(0, min(everything) - 0.03), min(1.0, max(everything) + 0.03))
    ax.set_ylabel("Pearson r", fontsize=6)
    ax.tick_params(axis="y", labelsize=5, width=0.5, length=2)
    ax.tick_params(axis="x", width=0.5, length=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)

    ax.set_title(
        f"Cerberus sequence prediction accuracy\n({n_folds} folds, {args.dataset})",
        fontsize=6, pad=3,
    )

    ax.legend(
        handles=[
            Line2D([0], [0], color="black", lw=1.0, label="Median"),
            plt.scatter([], [], s=5, color=RNA_DOT, alpha=0.85, label="RNA (avg ±)"),
            plt.scatter([], [], s=5, color=ATAC_DOT, alpha=0.85, label="ATAC"),
        ],
        fontsize=5, loc="lower left", frameon=True, framealpha=0.85,
        handletextpad=0.4, borderpad=0.4, labelspacing=0.3,
    )

    stem = os.path.join(out_dir, "cerberus_accuracy")
    fig.savefig(f"{stem}.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {stem}.pdf / .png")


if __name__ == "__main__":
    main()
