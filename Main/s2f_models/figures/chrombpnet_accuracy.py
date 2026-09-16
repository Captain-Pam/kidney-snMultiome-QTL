#!/usr/bin/env python
"""ChromBPNet held-out prediction accuracy.

Left: Pearson r between predicted and observed log counts at held-out peaks,
per cell type, across folds. Right: the per-peak values for one example fold as
a density scatter.

Reads the two CSVs written by `../chrombpnet/5_export_metrics.py`, so this runs
without a GPU or any model files.

Usage:
    python chrombpnet_accuracy.py [--metrics_dir <dir>] [--out_dir <dir>]
"""
import argparse
import os

import numpy as np
import pandas as pd

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

from figlib import CT_LABELS, cfg, qtl_ct, setup_style

BOX_FACE, BOX_EDGE, DOT_COLOR = "#AED6F1", "#2A6DB5", "#1A5276"
EXAMPLE_DOT = "#C0392B"


def density_scatter(ax, x, y, bins=200, cmap="YlOrRd", size=1):
    """Scatter coloured by local 2-D point density, drawn densest-last.

    Rasterized so the PDF holds a bitmap rather than tens of thousands of
    vector points.
    """
    counts, xedges, yedges = np.histogram2d(x, y, bins=bins)
    xi = np.searchsorted(xedges[1:-1], x)
    yi = np.searchsorted(yedges[1:-1], y)
    density = counts[xi, yi]
    order = density.argsort()
    return ax.scatter(
        x[order], y[order], c=density[order], s=size,
        cmap=cmap, linewidths=0, rasterized=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics_dir", default=None)
    parser.add_argument("--out_dir", default=None)
    args = parser.parse_args()

    setup_style()
    metrics_dir = args.metrics_dir or os.path.join(
        cfg("WORK_DIR"), "chrombpnet", "metrics"
    )
    out_dir = args.out_dir or cfg("FIGURE_DIR")
    os.makedirs(out_dir, exist_ok=True)

    pearson = pd.read_csv(os.path.join(metrics_dir, "pearson_r.csv"))
    scatter = pd.read_csv(os.path.join(metrics_dir, "scatter_data.csv"))

    by_ct = {ct: group for ct, group in pearson.groupby("cell_type")}
    order = sorted(by_ct, key=lambda ct: by_ct[ct]["pearsonr"].mean(), reverse=True)
    for celltype in order:
        print(f"  {celltype:12s} mean r = {by_ct[celltype]['pearsonr'].mean():.4f}")

    # The fold drawn on the right is highlighted in the box plot on the left.
    example_ct = scatter["cell_type"].iloc[0]
    example_fold = int(scatter["fold"].iloc[0])

    fig = plt.figure(figsize=(3.6, 1.0))
    grid = gridspec.GridSpec(1, 2, width_ratios=[1.4, 1], wspace=0.45)

    # ── Per-cell-type accuracy ────────────────────────────────────────────────
    ax_box = fig.add_subplot(grid[0])
    rng = np.random.default_rng(42)
    positions = np.arange(len(order))
    values = [by_ct[ct]["pearsonr"].values for ct in order]

    ax_box.boxplot(
        values,
        positions=positions,
        widths=0.45,
        patch_artist=True,
        medianprops=dict(color="black", linewidth=1.0),
        boxprops=dict(facecolor=BOX_FACE, edgecolor=BOX_EDGE, linewidth=0.5),
        whiskerprops=dict(color=BOX_EDGE, linewidth=0.5),
        capprops=dict(color=BOX_EDGE, linewidth=0.5),
        flierprops=dict(marker="", linestyle="none"),
        showfliers=False,
        zorder=2,
    )

    for position, celltype in zip(positions, order):
        folds = by_ct[celltype]
        jitter = rng.uniform(-0.14, 0.14, size=len(folds))
        colors = [
            EXAMPLE_DOT
            if (celltype == example_ct and fold == example_fold)
            else DOT_COLOR
            for fold in folds["fold"].values
        ]
        ax_box.scatter(
            position + jitter, folds["pearsonr"].values,
            s=5, c=colors, alpha=0.85, zorder=3, linewidths=0,
        )

    ax_box.set_xticks(positions)
    ax_box.set_xticklabels(
        [CT_LABELS.get(qtl_ct(ct), ct) for ct in order],
        rotation=90, ha="center", fontsize=5,
    )
    ax_box.set_ylabel("Pearson r (log counts)", fontsize=6)
    ax_box.set_title(
        f"ChromBPNet counts prediction accuracy\n"
        f"({pearson['fold'].nunique()} folds, peaks with log counts > 2)",
        fontsize=6, pad=3,
    )

    flat = pearson["pearsonr"].values
    ax_box.set_ylim(max(0, flat.min() - 0.03), min(1.0, flat.max() + 0.03))
    ax_box.tick_params(axis="y", labelsize=5, width=0.5, length=2)
    ax_box.tick_params(axis="x", width=0.5, length=2)
    ax_box.spines[["top", "right"]].set_visible(False)
    for spine in ["left", "bottom"]:
        ax_box.spines[spine].set_linewidth(0.5)

    ax_box.legend(
        handles=[
            Line2D([0], [0], color="black", lw=1.0, label="Median"),
            plt.scatter([], [], s=5, color=DOT_COLOR, alpha=0.85, label="Individual fold"),
            plt.scatter([], [], s=5, color=EXAMPLE_DOT, alpha=0.85, label="Shown at right"),
        ],
        fontsize=5, loc="upper right", frameon=True, framealpha=0.8,
        handletextpad=0.3, borderpad=0.4, labelspacing=0.3,
    )

    # Shrink the box panel vertically so it does not tower over the scatter.
    fig.canvas.draw()
    box_position = ax_box.get_position()
    height = box_position.height * 0.70
    ax_box.set_position(
        [
            box_position.x0,
            box_position.y0 + (box_position.height - height) / 2,
            box_position.width,
            height,
        ]
    )

    # ── Example fold ──────────────────────────────────────────────────────────
    ax = fig.add_subplot(grid[1])
    observed = scatter["observed"].values
    predicted = scatter["predicted"].values
    r = pearsonr(observed, predicted)[0]

    points = density_scatter(ax, observed, predicted)
    colorbar = plt.colorbar(points, ax=ax, pad=0.02, shrink=0.85)
    colorbar.ax.tick_params(labelsize=5, width=0.4, length=1.5)
    colorbar.set_label("Density", fontsize=5)
    colorbar.outline.set_linewidth(0.4)

    limits = [
        min(observed.min(), predicted.min()) - 0.2,
        max(observed.max(), predicted.max()) + 0.2,
    ]
    ax.plot(limits, limits, ls=":", lw=0.6, color="grey", alpha=0.8)
    ax.set_xlim(limits)
    ax.set_ylim(limits)

    ax.set_xlabel("Observed log counts", fontsize=6)
    ax.set_ylabel("Predicted log counts", fontsize=6)
    ax.set_title(
        f"{CT_LABELS.get(qtl_ct(example_ct), example_ct)}, fold {example_fold}\n"
        f"Pearson r = {r:.4f}  (n = {len(observed):,})",
        fontsize=6, pad=3,
    )
    ax.tick_params(labelsize=5, width=0.5, length=2)
    ax.spines[["top", "right"]].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(0.5)

    stem = os.path.join(out_dir, "chrombpnet_accuracy")
    fig.savefig(f"{stem}.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {stem}.pdf / .png")


if __name__ == "__main__":
    main()
