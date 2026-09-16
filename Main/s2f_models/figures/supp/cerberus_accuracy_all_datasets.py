#!/usr/bin/env python
"""
plot_supp_panel_e.py

Supplementary panel E: Cerberus (fine-tuned Cerberus) sequence prediction accuracy
(Pearson R) across 8 test folds, for ALL datasets in the fine-tune — human head
(eval0: Susztak, Ledru, Muto PKD, Muto AKI) and mouse head (eval1: White, Chen,
Muto PKD, Muto IRI, Kirita). RNA (RNA+/RNA- pooled) and ATAC boxes per dataset;
each point is one track, Pearson r averaged across the 8 test folds.
Style matches plot_supp_panel_d.py.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figlib import (  # noqa: E402
    ALT_COLOR, ATAC_IDX, CT_COLORS, CT_LABELS, NT_COLORS, REF_COLOR,
    add_ols, add_ref_lines, cfg, model_ct, qtl_ct, select_lead,
    setup_style, style_ax,
)

setup_style()

# ── Paths ──────────────────────────────────────────────────────────────────────
EVAL_BASE = cfg("CERBERUS_MODEL_DIR")
OUT_PDF   = os.path.join(cfg("FIGURE_DIR"), "cerberus_accuracy_all_datasets.pdf")
OUT_PNG   = os.path.join(cfg("FIGURE_DIR"), "cerberus_accuracy_all_datasets.png")

# ── Load 8 test folds, human (eval0) + mouse (eval1) ───────────────────────────
records = []
for species, evald in [("Human", "eval0"), ("Mouse", "eval1")]:
    for fi in range(8):
        df = pd.read_csv(f"{EVAL_BASE}/f{fi}c0/{evald}/fold{fi}/acc.txt", sep="\t")
        df["model_fold"] = fi
        df["species"]    = species
        records.append(df)
all_df = pd.concat(records, ignore_index=True)
all_df[["modality", "celltype", "dataset"]] = (
    all_df["description"].str.split(":", n=2, expand=True)
)
# modality group: pool RNA strands
all_df["mod_group"] = all_df["modality"].map(
    {"ATAC": "ATAC", "RNA+": "RNA", "RNA-": "RNA"}
)

# per-track accuracy = mean Pearson r across the 8 folds
per_track = (all_df
             .groupby(["species", "dataset", "mod_group", "identifier"])["pearsonr"]
             .mean()
             .reset_index())

# ── Dataset order + display labels (human first, then mouse) ────────────────────
DATASETS = [
    ("Human", "Susztak_CKD_multiome",   "Susztak\nCKD"),
    ("Human", "Ledru_Control_multiome", "Ledru\nControl"),
    ("Human", "Muto_PKD_paired",        "Muto\nPKD"),
    ("Human", "Muto_AKI_paired",        "Muto\nAKI"),
    ("Mouse", "White_CKD_multiome",     "White\nCKD"),
    ("Mouse", "Chen_timeseries_multiome","Chen\ntime-series"),
    ("Mouse", "Muto_PKD_multiome",      "Muto\nPKD"),
    ("Mouse", "Muto_IRI_ATAC",          "Muto\nIRI"),
    ("Mouse", "Kirita_IRI_RNA",         "Kirita\nIRI"),
]

def vals(species, dataset, mod):
    m = ((per_track["species"] == species) & (per_track["dataset"] == dataset)
         & (per_track["mod_group"] == mod))
    return per_track.loc[m, "pearsonr"].values

# ── Colors (match supp_panel_d) ────────────────────────────────────────────────
RNA_BOX  = "#AED6F1";  RNA_EDGE  = "#2A6DB5";  RNA_DOT  = "#1A5276"
ATAC_BOX = "#A9DFBF";  ATAC_EDGE = "#1E8449";  ATAC_DOT = "#145A32"

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5.2, 1.15))

DS_SPACING = 1.4
OFFSET     = 0.24
W          = 0.45
rng = np.random.default_rng(42)

n       = len(DATASETS)
centres = np.arange(n) * DS_SPACING

def draw_box(ax, data, pos, box_fc, box_ec, dot_color, width):
    if len(data) == 0:
        return
    ax.boxplot(
        [data], positions=[pos], widths=width, patch_artist=True,
        medianprops=dict(color="black", linewidth=1.0),
        boxprops=dict(facecolor=box_fc, edgecolor=box_ec, linewidth=0.5),
        whiskerprops=dict(color=box_ec, linewidth=0.5),
        capprops=dict(color=box_ec, linewidth=0.5),
        flierprops=dict(marker="", linestyle="none"),
        showfliers=False, zorder=2,
    )
    jit = rng.uniform(-0.09, 0.09, size=len(data))
    ax.scatter(pos + jit, data, s=3, color=dot_color,
               alpha=0.6, zorder=3, linewidths=0)

all_vals = []
for i, (sp, ds, _lab) in enumerate(DATASETS):
    rv = vals(sp, ds, "RNA")
    av = vals(sp, ds, "ATAC")
    draw_box(ax, rv, centres[i] - OFFSET, RNA_BOX,  RNA_EDGE,  RNA_DOT,  W)
    draw_box(ax, av, centres[i] + OFFSET, ATAC_BOX, ATAC_EDGE, ATAC_DOT, W)
    all_vals.extend(list(rv) + list(av))

# ── Species divider + labels ───────────────────────────────────────────────────
n_human = sum(1 for sp, _, _ in DATASETS if sp == "Human")
div_x   = (centres[n_human - 1] + centres[n_human]) / 2
ax.axvline(div_x, color="0.6", lw=0.5, ls=(0, (3, 3)), zorder=1)

# ── Axes ───────────────────────────────────────────────────────────────────────
ax.set_xticks(centres)
ax.set_xticklabels([lab for _, _, lab in DATASETS], rotation=0, ha="center", fontsize=5)
ax.set_xlim(-0.7, centres[-1] + 0.7)

y_lo = max(0.0, min(all_vals) - 0.03)
y_hi = min(1.0, max(all_vals) + 0.05)
ax.set_ylim(y_lo, y_hi)
ax.set_ylabel("Pearson r", fontsize=6)
ax.tick_params(axis="y", labelsize=5, width=0.5, length=2)
ax.tick_params(axis="x", width=0.5, length=2)
ax.spines[["top", "right"]].set_visible(False)
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)

# species headers
ax.text((centres[0] + centres[n_human - 1]) / 2, y_hi, "Human",
        ha="center", va="bottom", fontsize=6, fontweight="bold")
ax.text((centres[n_human] + centres[-1]) / 2, y_hi, "Mouse",
        ha="center", va="bottom", fontsize=6, fontweight="bold")

ax.set_title(
    "Cerberus sequence prediction accuracy across fine-tune datasets\n"
    "(8 test folds; each point = one track, mean r over folds)",
    fontsize=6, pad=10,
)

legend_handles = [
    Line2D([0], [0], color="black", lw=1.0, label="Median"),
    plt.scatter([], [], s=4, color=RNA_DOT,  alpha=0.85, label="RNA (RNA+/RNA-)"),
    plt.scatter([], [], s=4, color=ATAC_DOT, alpha=0.85, label="ATAC"),
]
ax.legend(handles=legend_handles, fontsize=5, loc="lower left",
          frameon=True, framealpha=0.85, handletextpad=0.4,
          borderpad=0.4, labelspacing=0.3)

fig.savefig(OUT_PDF, dpi=600, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")

# ── Report medians ─────────────────────────────────────────────────────────────
print("\nMedian Pearson r by dataset:")
for sp, ds, lab in DATASETS:
    rv = vals(sp, ds, "RNA"); av = vals(sp, ds, "ATAC")
    rmed = f"{np.median(rv):.3f}" if len(rv) else "  -  "
    amed = f"{np.median(av):.3f}" if len(av) else "  -  "
    print(f"  {sp:5s} {ds:26s}  RNA={rmed} (n={len(rv):3d})  ATAC={amed} (n={len(av):3d})")
