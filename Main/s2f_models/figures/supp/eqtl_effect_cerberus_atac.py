#!/usr/bin/env python
"""
merged_panel_h_borzoi_atac.py

Same style/layout as merged_panel_h_v2.py, but y = Cerberus ATAC-track local logSUM
(8-fold ensemble, matched cell type) instead of Cerberus RNA logSED.
eQTL beta is swap-aligned (slope_adj), identical convention to panel_h.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ── Font / PDF ────────────────────────────────────────────────────────────────
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figlib import (  # noqa: E402
    ALT_COLOR, ATAC_IDX, CT_COLORS, CT_LABELS, NT_COLORS, REF_COLOR,
    add_ols, add_ref_lines, cfg, model_ct, qtl_ct, select_lead,
    setup_style, style_ax,
)

setup_style()

# ── Paths ──────────────────────────────────────────────────────────────────────
MERGED = os.path.join(cfg("WORK_DIR"), "eqtl", "atac_logSUM", "eqtl_atac_logSUM_ensemble_merged.tsv")
OUT_PDF = os.path.join(cfg("FIGURE_DIR"), "eqtl_effect_cerberus_atac.pdf")
OUT_PNG = os.path.join(cfg("FIGURE_DIR"), "eqtl_effect_cerberus_atac.png")
YCOL = "borzoi_atac"
PIP_THRESH = 0.5

CELL_TYPES = ["CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
              "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT"]
CT_DISPLAY = {"CNT_CD_PC": "CNT/CD/PC", "DCT": "DCT", "DTL_ATL": "DTL/ATL",
              "EC": "EC", "IC": "IC", "Immune": "Immune", "PEC": "PEC",
              "PTS": "PTS", "Podocyte": "Podocyte", "Stromal": "Stromal",
              "TAL": "TAL", "injPT": "injPT"}
DOT_SIZE, DOT_ALPHA = 6, 0.55

# ── Load ───────────────────────────────────────────────────────────────────────
df = pd.read_csv(MERGED, sep="\t")
sub = df[df["variable_prob"] >= PIP_THRESH].dropna(subset=[YCOL, "slope_adj"]).copy()

x_all = sub["slope_adj"].values.astype(float)
y_all = sub[YCOL].values.astype(float)
r_all, _ = stats.pearsonr(x_all, y_all)

sub_pts = sub[sub["celltype"] == "PTS"]
r_pts, _ = stats.pearsonr(sub_pts["slope_adj"].values.astype(float),
                          sub_pts[YCOL].values.astype(float))
print(f"n={len(sub)}  overall r={r_all:.3f}  PTS r={r_pts:.3f}  n_PTS={len(sub_pts)}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(2.016, 2.016))

ct_order = [ct for ct in CELL_TYPES if ct != "PTS"] + ["PTS"]
for ct in ct_order:
    s = sub[sub["celltype"] == ct]
    if len(s) == 0:
        continue
    ax.scatter(s["slope_adj"], s[YCOL], s=DOT_SIZE, alpha=DOT_ALPHA, linewidths=0,
               color=CT_COLORS[ct], label=CT_DISPLAY[ct], rasterized=True, zorder=3)

m, b = np.polyfit(x_all, y_all, 1)
x_line = np.array([x_all.min(), x_all.max()])
ax.plot(x_line, m * x_line + b, ls=":", lw=0.8, color="grey", alpha=0.9, zorder=4)

x_pts = sub_pts["slope_adj"].values.astype(float)
y_pts = sub_pts[YCOL].values.astype(float)
m_pts, b_pts = np.polyfit(x_pts, y_pts, 1)
x_line_pts = np.array([x_pts.min(), x_pts.max()])
ax.plot(x_line_pts, m_pts * x_line_pts + b_pts, ls=":", lw=0.8,
        color=CT_COLORS["PTS"], alpha=0.9, zorder=5)

ax.axhline(0, color="lightgrey", lw=0.4, zorder=1)
ax.axvline(0, color="lightgrey", lw=0.4, zorder=1)

ax.text(0.04, 0.97,
        f"Overall: r={r_all:.2f}, n={len(sub):,}\nPTS: r={r_pts:.2f}, n={len(sub_pts):,}",
        transform=ax.transAxes, ha="left", va="top", fontsize=5.5, linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))

leg = ax.legend(fontsize=4.5, ncol=2, loc="lower right", framealpha=0.8, frameon=True,
                handletextpad=0.3, columnspacing=0.5, borderpad=0.4, labelspacing=0.3,
                edgecolor="#cccccc")
leg.get_frame().set_linewidth(0.4)

ax.set_xlabel("eQTL β", fontsize=6, labelpad=2)
ax.set_ylabel("Cerberus ATAC logSUM", fontsize=6, labelpad=2)
ax.set_title(f"eQTL β vs Cerberus ATAC logSUM\n(PIP ≥ {PIP_THRESH})",
             fontsize=6, pad=3, linespacing=1.4)

ax.spines[["top", "right"]].set_visible(False)
for sp in ["left", "bottom"]:
    ax.spines[sp].set_linewidth(0.4)
ax.tick_params(labelsize=6, width=0.4, length=2, pad=1.5)
ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="both"))
ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="both"))

fig.savefig(OUT_PDF, dpi=600, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_PDF}")
