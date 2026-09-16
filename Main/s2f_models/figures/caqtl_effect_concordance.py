#!/usr/bin/env python
"""Predicted versus measured caQTL effect, for both models.

Two rows of scatter panels — ChromBPNet `logfc` on top, Cerberus `logSUM`
below — against the RASQUAL allele-specific effect size, plus a summary bar of
Pearson r across all 12 cell types.

One point per ATAC peak: the peak's variants are typically in tight LD with
indistinguishable measured effects, so each peak contributes only its
top-scoring variant (`figlib.select_lead`). RASQUAL effect sizes run on (0, 1)
centred at 0.5, so the reference lines are drawn at x = 0.5 and y = 0.

Usage:
    python caqtl_effect_concordance.py [--out_dir <dir>]
"""
import argparse
import os

import numpy as np
from scipy import stats

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from caqtl_data import CELL_TYPES, load_merged
from figlib import (
    CT_COLORS,
    CT_LABELS,
    add_ols,
    add_ref_lines,
    cfg,
    select_lead,
    setup_style,
    style_ax,
)

# The four examples shown in the caQTL-examples figure, outlined here so the
# reader can place them in the overall distribution.
EXAMPLE_SNPS = {
    74669588: "PT",     # chr18:74669588  PT   HNF1A
    59107669: "TAL",    # chr1:59107669   TAL  MEF2A/B
    181152560: "IC",    # chr5:181152560  IC   FOXS1/D1
    107693236: "EC",    # chr4:107693236  EC   FOXO1::ELK1/3
}

# Cell types given their own scatter panel, ordered by r at plot time.
SCATTER_CTS = ["PT", "TAL", "IC", "EC"]

DOT_SIZE = 0.8
DOT_ALPHA = 0.25
POINT_COLOR = "#3498DB"


def plot_summary_r(ax, r_by_ct, color, xlabel):
    """Horizontal bar of Pearson r per cell type, ascending."""
    order = sorted(r_by_ct, key=lambda c: r_by_ct[c])
    values = np.array([r_by_ct[c] for c in order])
    positions = np.arange(len(order))

    ax.barh(positions, values, height=0.6, color=color, linewidth=0, zorder=2)
    ax.set_yticks(positions)
    ax.set_yticklabels([CT_LABELS.get(c, c) for c in order], fontsize=5)
    ax.tick_params(axis="x", labelsize=5, width=0.4, length=2, pad=1.5)
    ax.tick_params(axis="y", width=0, length=0, pad=2)
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=3, prune="both"))

    # Zoom the axis to the observed range so differences between cell types
    # stay readable, rather than anchoring at zero.
    pad = max(0.02, (values.max() - values.min()) * 0.20)
    xmin = max(0.0, values.min() - pad)
    xmax = min(1.0, values.max() + pad)
    ax.set_xlim(xmin, xmax)
    ax.axvline(xmin, color="lightgrey", lw=0.4, zorder=1)

    ax.set_xlabel(xlabel, fontsize=5, labelpad=2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)


def draw_scatter_row(fig, grid, row, per_ct, r_by_ct, highlights, order, ylabel):
    for column, celltype in enumerate(order):
        ax = fig.add_subplot(grid[row, column + 1])
        x, y = per_ct[celltype]
        r = r_by_ct.get(celltype, float("nan"))

        ax.scatter(
            x, y, s=DOT_SIZE, color=POINT_COLOR, alpha=DOT_ALPHA,
            linewidths=0, rasterized=True, zorder=2,
        )
        add_ols(ax, x, y)
        add_ref_lines(ax, 0.5, 0)

        for hx, hy, color in highlights.get(celltype, []):
            ax.scatter(
                hx, hy, s=8, color=color,
                edgecolors="black", linewidths=0.4, zorder=5,
            )

        ax.set_title(
            f"{CT_LABELS[celltype]}\nr={r:.2f}, n={len(x):,}",
            fontsize=6, pad=2.5, linespacing=1.4,
        )
        ax.set_xlabel("RASQUAL effect size", fontsize=6, labelpad=2)
        ax.set_ylabel(ylabel if column == 0 else "", fontsize=6, labelpad=2)
        style_ax(ax)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", default=None)
    args = parser.parse_args()

    setup_style()

    out_dir = args.out_dir or cfg("FIGURE_DIR")
    os.makedirs(out_dir, exist_ok=True)

    print("Loading caQTL scores...")
    merged = load_merged()

    chrombpnet_lead, cerberus_lead = {}, {}
    chrombpnet_highlights, cerberus_highlights = {}, {}
    r_chrombpnet, r_cerberus = {}, {}

    for celltype in CELL_TYPES:
        subset = merged[merged["cell_type"] == celltype].copy()
        cbp_col = f"logfc_{celltype}"
        cer_col = f"logSUM_{celltype}"

        for column, store, correlations in (
            (cbp_col, chrombpnet_lead, r_chrombpnet),
            (cer_col, cerberus_lead, r_cerberus),
        ):
            lead = select_lead(subset, column)
            x = lead["Effect_size"].values.astype(float)
            y = lead[column].values.astype(float)
            finite = np.isfinite(x) & np.isfinite(y)
            store[celltype] = (x[finite], y[finite])
            if finite.sum() >= 5:
                correlations[celltype] = stats.pearsonr(x[finite], y[finite])[0]

        for snp_pos, example_ct in EXAMPLE_SNPS.items():
            if example_ct != celltype:
                continue
            rows = merged[
                (merged["SNP_position"] == snp_pos)
                & (merged["cell_type"] == celltype)
            ]
            if not len(rows):
                continue
            row = rows.iloc[0]
            effect = float(row["Effect_size"])
            color = CT_COLORS[celltype]
            for column, store in ((cbp_col, chrombpnet_highlights),
                                  (cer_col, cerberus_highlights)):
                value = float(row[column])
                if np.isfinite(effect) and np.isfinite(value):
                    store.setdefault(celltype, []).append((effect, value, color))

    order = sorted(SCATTER_CTS, key=lambda c: r_chrombpnet.get(c, 0), reverse=True)
    for celltype in order:
        print(
            f"  {CT_LABELS[celltype]:12s} "
            f"r_chrombpnet={r_chrombpnet.get(celltype, float('nan')):.3f} "
            f"n={len(chrombpnet_lead[celltype][0]):,}   "
            f"r_cerberus={r_cerberus.get(celltype, float('nan')):.3f} "
            f"n={len(cerberus_lead[celltype][0]):,}"
        )

    fig = plt.figure(figsize=(5.1, 2.31))
    grid = gridspec.GridSpec(
        2, len(order) + 1, hspace=0.55, wspace=0.4,
        height_ratios=[1, 1], width_ratios=[1.5] + [1] * len(order),
    )

    ax_top = fig.add_subplot(grid[0, 0])
    plot_summary_r(
        ax_top, r_chrombpnet, POINT_COLOR,
        "Pearson r\n(ChromBPNet logfc\nvs RASQUAL)",
    )
    ax_top.set_title("All cell types", fontsize=6, pad=2.5, fontweight="bold")

    ax_bottom = fig.add_subplot(grid[1, 0])
    plot_summary_r(
        ax_bottom, r_cerberus, POINT_COLOR,
        "Pearson r\n(Cerberus logSUM\nvs RASQUAL)",
    )

    draw_scatter_row(
        fig, grid, 0, chrombpnet_lead, r_chrombpnet,
        chrombpnet_highlights, order, "ChromBPNet logfc",
    )
    draw_scatter_row(
        fig, grid, 1, cerberus_lead, r_cerberus,
        cerberus_highlights, order, "Cerberus logSUM",
    )

    stem = os.path.join(out_dir, "caqtl_effect_concordance")
    fig.savefig(f"{stem}.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {stem}.pdf / .png")


if __name__ == "__main__":
    main()
