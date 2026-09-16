#!/usr/bin/env python
"""Cerberus logSUM versus ChromBPNet logfc, across cell types.

Two models with very different architectures and receptive fields — 786 kb of
context versus 2,114 bp — scored on the same variants. Points are coloured by
the measured RASQUAL effect size, so it is visible that the models agree most
where the measured effect is strongest.

Unlike the effect-concordance figure this uses **all** scored variants, not one
lead per peak: the question here is whether the two models agree variant by
variant, and there is no measured quantity to de-duplicate against.

Usage:
    python caqtl_model_concordance.py [--out_dir <dir>]
"""
import argparse
import os

import numpy as np
from scipy import stats

import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from caqtl_data import CELL_TYPES, load_merged
from figlib import (
    CT_LABELS,
    add_ols,
    add_ref_lines,
    cfg,
    setup_style,
    style_ax,
)

DOT_SIZE = 0.8
DOT_ALPHA = 0.25
BAR_COLOR = "#888888"

# Cell type given its own scatter panel.
FOCUS_CT = "PT"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", default=None)
    args = parser.parse_args()

    setup_style()
    out_dir = args.out_dir or cfg("FIGURE_DIR")
    os.makedirs(out_dir, exist_ok=True)

    print("Loading caQTL scores...")
    merged = load_merged()

    per_ct = {}
    correlations = {}
    for celltype in CELL_TYPES:
        subset = merged[merged["cell_type"] == celltype]
        chrombpnet = subset[f"logfc_{celltype}"].values.astype(float)
        cerberus = subset[f"logSUM_{celltype}"].values.astype(float)
        effect = subset["Effect_size"].values.astype(float)

        finite = np.isfinite(chrombpnet) & np.isfinite(cerberus) & np.isfinite(effect)
        if finite.sum() < 5:
            continue
        per_ct[celltype] = (chrombpnet[finite], cerberus[finite], effect[finite])
        correlations[celltype] = stats.pearsonr(
            chrombpnet[finite], cerberus[finite]
        )[0]

    for celltype, r in sorted(correlations.items(), key=lambda kv: -kv[1]):
        print(f"  {CT_LABELS[celltype]:12s} r={r:.3f}  n={len(per_ct[celltype][0]):,}")

    # Diverging colour scale centred on 0.5, the RASQUAL no-imbalance point.
    # Percentile limits keep a few extreme variants from washing out the rest.
    all_effects = np.concatenate([per_ct[c][2] for c in per_ct])
    colormap = mcolors.LinearSegmentedColormap.from_list(
        "effect_div",
        [(0.0, "#2166AC"), (0.5, "#D8D8D8"), (1.0, "#B2182B")],
        N=256,
    )
    norm = mcolors.TwoSlopeNorm(
        vmin=np.percentile(all_effects, 2),
        vcenter=0.5,
        vmax=np.percentile(all_effects, 98),
    )

    fig = plt.figure(figsize=(0.75, 2.31))
    grid = gridspec.GridSpec(
        2, 2, height_ratios=[1, 1], width_ratios=[1, 0.07],
        hspace=0.55, wspace=0.12,
    )

    # ── Per-cell-type correlation ─────────────────────────────────────────────
    ax_bar = fig.add_subplot(grid[0, :])
    order = sorted(correlations, key=lambda c: correlations[c])
    values = np.array([correlations[c] for c in order])
    positions = np.arange(len(order))

    ax_bar.barh(positions, values, height=0.6, color=BAR_COLOR, linewidth=0, zorder=2)
    ax_bar.set_yticks(positions)
    ax_bar.set_yticklabels([CT_LABELS.get(c, c) for c in order], fontsize=5)
    ax_bar.tick_params(axis="x", labelsize=5, width=0.4, length=2, pad=1.5)
    ax_bar.tick_params(axis="y", width=0, length=0, pad=2)
    ax_bar.xaxis.set_major_locator(plt.MaxNLocator(nbins=3, prune="both"))

    pad = max(0.02, (values.max() - values.min()) * 0.20)
    xmin = max(0.0, values.min() - pad)
    ax_bar.set_xlim(xmin, min(1.0, values.max() + pad))
    ax_bar.axvline(xmin, color="lightgrey", lw=0.4, zorder=1)
    ax_bar.set_xlabel(
        "Pearson r\n(Cerberus logSUM vs ChromBPNet logfc)", fontsize=5, labelpad=2
    )
    ax_bar.set_title("All cell types", fontsize=6, pad=3, fontweight="bold")
    ax_bar.spines[["top", "right", "left"]].set_visible(False)
    ax_bar.spines["bottom"].set_linewidth(0.4)

    # ── Focus scatter ─────────────────────────────────────────────────────────
    ax = fig.add_subplot(grid[1, 0])
    ax_cbar = fig.add_subplot(grid[1, 1])

    chrombpnet, cerberus, effect = per_ct[FOCUS_CT]
    ax.scatter(
        chrombpnet, cerberus, s=DOT_SIZE, c=effect, cmap=colormap, norm=norm,
        alpha=DOT_ALPHA, linewidths=0, rasterized=True, zorder=2,
    )
    add_ols(ax, chrombpnet, cerberus)
    add_ref_lines(ax, 0, 0)

    ax.set_title(
        f"{CT_LABELS[FOCUS_CT]}\nr={correlations[FOCUS_CT]:.2f}, n={len(chrombpnet):,}",
        fontsize=6, pad=2.5, linespacing=1.4,
    )
    ax.set_xlabel("ChromBPNet logfc", fontsize=6, labelpad=2)
    ax.set_ylabel("Cerberus logSUM", fontsize=6, labelpad=2)
    style_ax(ax)

    mappable = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    mappable.set_array([])
    colorbar = fig.colorbar(mappable, cax=ax_cbar)
    colorbar.set_label("RASQUAL\neffect size", fontsize=5, labelpad=2)
    colorbar.ax.tick_params(labelsize=4, width=0.4, length=2)
    colorbar.outline.set_linewidth(0.4)

    stem = os.path.join(out_dir, "caqtl_model_concordance")
    fig.savefig(f"{stem}.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {stem}.pdf / .png")


if __name__ == "__main__":
    main()
