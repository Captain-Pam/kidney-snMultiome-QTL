#!/usr/bin/env python
"""Cerberus logSED versus fine-mapped eQTL effect size.

One point per variant-gene pair among confidently fine-mapped credible-set
variants (PIP >= 0.5), coloured by cell type. logSED is the log ratio of
predicted RNA coverage summed over the eGene's exons between alleles.

The eQTL slope has already been aligned to hg38 ALT-vs-REF orientation when the
data is loaded — see `eqtl_data.load_eqtl`.

Usage:
    python eqtl_effect_concordance.py [--out_dir <dir>] [--pip 0.5]
"""
import argparse
import os

import numpy as np
from scipy import stats

import matplotlib.pyplot as plt

from eqtl_data import CELL_TYPES, load_matched
from figlib import CT_COLORS, CT_LABELS, cfg, qtl_ct, setup_style

DOT_SIZE = 6
DOT_ALPHA = 0.55

# Cell type called out separately: the largest cell population, and the one the
# example loci are drawn from.
FOCUS_CT = "PTS"

# The two examples shown in the eQTL-examples figure, outlined here.
EXAMPLE_VARIANTS = [
    ("chr7:94392163", "EC"),   # COL1A2-AS1
    ("chr19:9324196", "IC"),   # ZNF559
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", default=None)
    parser.add_argument("--pip", type=float, default=0.5, help="minimum PIP")
    args = parser.parse_args()

    setup_style()
    out_dir = args.out_dir or cfg("FIGURE_DIR")
    os.makedirs(out_dir, exist_ok=True)

    data = load_matched(pip_threshold=args.pip)

    slope = data["slope"].values.astype(float)
    logsed = data["logSED"].values.astype(float)
    r_all = stats.pearsonr(slope, logsed)[0]

    focus = data[data["celltype"] == FOCUS_CT]
    focus_slope = focus["slope"].values.astype(float)
    focus_logsed = focus["logSED"].values.astype(float)
    r_focus = stats.pearsonr(focus_slope, focus_logsed)[0]

    print(
        f"  overall r={r_all:.3f} n={len(data):,}   "
        f"{FOCUS_CT} r={r_focus:.3f} n={len(focus):,}"
    )

    fig, ax = plt.subplots(figsize=(2.016, 2.016))

    # Draw the focus cell type last so its points are not buried.
    order = [ct for ct in CELL_TYPES if ct != FOCUS_CT] + [FOCUS_CT]
    for celltype in order:
        subset = data[data["celltype"] == celltype]
        if not len(subset):
            continue
        key = qtl_ct(celltype)
        ax.scatter(
            subset["slope"], subset["logSED"],
            s=DOT_SIZE, alpha=DOT_ALPHA, linewidths=0,
            color=CT_COLORS[key], label=CT_LABELS[key],
            rasterized=True, zorder=3,
        )

    for values_x, values_y, color, zorder in (
        (slope, logsed, "grey", 4),
        (focus_slope, focus_logsed, CT_COLORS[qtl_ct(FOCUS_CT)], 5),
    ):
        fit_slope, intercept = np.polyfit(values_x, values_y, 1)
        line_x = np.array([values_x.min(), values_x.max()])
        ax.plot(
            line_x, fit_slope * line_x + intercept,
            ls=":", lw=0.8, color=color, alpha=0.9, zorder=zorder,
        )

    for variant_id, celltype in EXAMPLE_VARIANTS:
        row = data[
            (data["variant_id"] == variant_id) & (data["celltype"] == celltype)
        ]
        if len(row):
            ax.scatter(
                row["slope"], row["logSED"], s=18,
                color=CT_COLORS[qtl_ct(celltype)],
                edgecolors="black", linewidths=0.6, zorder=10,
            )

    ax.axhline(0, color="lightgrey", lw=0.4, zorder=1)
    ax.axvline(0, color="lightgrey", lw=0.4, zorder=1)

    ax.text(
        0.04, 0.97,
        f"Overall: r={r_all:.2f}, n={len(data):,}\n"
        f"{FOCUS_CT}: r={r_focus:.2f}, n={len(focus):,}",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=5.5, linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8),
    )

    legend = ax.legend(
        fontsize=4.5, ncol=2, loc="lower right", framealpha=0.8, frameon=True,
        handletextpad=0.3, columnspacing=0.5, borderpad=0.4,
        labelspacing=0.3, edgecolor="#cccccc",
    )
    legend.get_frame().set_linewidth(0.4)

    ax.set_xlabel("eQTL β", fontsize=6, labelpad=2)
    ax.set_ylabel("Cerberus logSED", fontsize=6, labelpad=2)
    ax.set_title(
        f"eQTL β vs Cerberus logSED\n(PIP ≥ {args.pip})",
        fontsize=6, pad=3, linespacing=1.4,
    )

    ax.spines[["top", "right"]].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(0.4)
    ax.tick_params(labelsize=6, width=0.4, length=2, pad=1.5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="both"))
    ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="both"))

    stem = os.path.join(out_dir, "eqtl_effect_concordance")
    fig.savefig(f"{stem}.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {stem}.pdf / .png")


if __name__ == "__main__":
    main()
