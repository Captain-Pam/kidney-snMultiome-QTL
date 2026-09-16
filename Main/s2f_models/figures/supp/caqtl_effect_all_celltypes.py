#!/usr/bin/env python
"""
merged_panel_e_supp.py

Supplementary version: 3 rows × 12 cell types (all), sorted by n.
No summary column — all cell types visible in scatter.
"""
import os, glob
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
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

# ── Paths ─────────────────────────────────────────────────────────────────────
BORZOI_DIR  = os.path.join(cfg("WORK_DIR"), "caqtl", "scores", "logSUM", "local_ensemble")
MERGED_TSV  = os.path.join(cfg("WORK_DIR"), "chrombpnet", "caqtl", "caqtl_scores_merged.tsv.gz")
OUT_PDF     = os.path.join(cfg("FIGURE_DIR"), "caqtl_effect_all_celltypes.pdf")
OUT_PNG     = os.path.join(cfg("FIGURE_DIR"), "caqtl_effect_all_celltypes.png")

# ── Cell types ────────────────────────────────────────────────────────────────
RASQUAL_CELLTYPES = [
    "CNT_CD_PC","DCT","DTL_ATL","EC","IC","Immune","PEC",
    "PT","Podocyte","Stromal","TAL","injPT",
]

EXAMPLE_COLORS = {
    38888600:  "#D62728",
    101281470: "#FF7F0E",
    11862718:  "#2CA02C",
    2029767:   "#9467BD",
}
EXAMPLE_SNPS = {
    38888600:  ["PT"],
    101281470: ["CNT_CD_PC", "TAL"],
    11862718:  ["EC"],
    2029767:   ["IC"],
}

DOT_SIZE     = 0.8
DOT_ALPHA    = 0.25
RASQUAL_COLOR = "#3498DB"

# ── Helpers ───────────────────────────────────────────────────────────────────
# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading Cerberus scores...")
borzoi_rows, borzoi_raws = [], []
for cd in sorted(glob.glob(f"{BORZOI_DIR}/chunk_*"),
                 key=lambda p: int(p.split("_")[-1])):
    h5p = f"{cd}/scores.h5"
    if not os.path.exists(h5p): continue
    with h5py.File(h5p) as h5:
        borzoi_rows.append(pd.DataFrame({
            "chrom": h5["chr"][:].astype(str), "pos": h5["pos"][:],
            "ref":   h5["ref_allele"][:].astype(str),
            "alt":   h5["alt_allele"][:].astype(str),
        }))
        borzoi_raws.append(h5["cov/logSUM"][:].astype(np.float32))
borzoi_df  = pd.concat(borzoi_rows, ignore_index=True)
borzoi_raw = np.concatenate(borzoi_raws, axis=0)
borzoi_df["_key"] = (borzoi_df["chrom"]+"_"+borzoi_df["pos"].astype(str)
                     +"_"+borzoi_df["ref"]+"_"+borzoi_df["alt"])
for ct, idxs in ATAC_IDX.items():
    borzoi_df[f"logSUM_{ct}"] = borzoi_raw[:, idxs].mean(axis=1)
borzoi_df = borzoi_df.set_index("_key")

print("Loading merged table...")
merged = pd.read_csv(MERGED_TSV, sep="\t")
merged["_key"] = (merged["Chromosome"]+"_"+merged["SNP_position"].astype(str)
                  +"_"+merged["Ref_allele"]+"_"+merged["Alt_allele"])
for ct in RASQUAL_CELLTYPES:
    merged[f"logSUM_{ct}"] = merged["_key"].map(borzoi_df[f"logSUM_{ct}"])

# ── Build per-CT data ─────────────────────────────────────────────────────────
cbpnet_lead, borzoi_lead, borzoi_vs_cbpnet = {}, {}, {}
highlight_pts_cbpnet = {}
highlight_pts_borzoi = {}
r_cbpnet, r_borzoi, r_bvsc = {}, {}, {}

for ct in RASQUAL_CELLTYPES:
    sub = merged[merged["cell_type"] == ct].copy()

    cbp_col = f"logfc_{ct}"
    bor_col = f"logSUM_{ct}"

    lead_c = select_lead(sub, cbp_col)
    x = lead_c["Effect_size"].values.astype(float)
    y = lead_c[cbp_col].values.astype(float)
    mask = np.isfinite(x) & np.isfinite(y)
    cbpnet_lead[ct] = (x[mask], y[mask])
    if mask.sum() >= 5:
        r_cbpnet[ct], _ = stats.pearsonr(x[mask], y[mask])

    lead_b = select_lead(sub, bor_col)
    x = lead_b["Effect_size"].values.astype(float)
    y = lead_b[bor_col].values.astype(float)
    mask = np.isfinite(x) & np.isfinite(y)
    borzoi_lead[ct] = (x[mask], y[mask])
    if mask.sum() >= 5:
        r_borzoi[ct], _ = stats.pearsonr(x[mask], y[mask])

    ok = sub[cbp_col].notna() & sub[bor_col].notna() & sub["Effect_size"].notna()
    sub_ok = sub[ok]
    x = sub_ok[cbp_col].values.astype(float)
    y = sub_ok[bor_col].values.astype(float)
    c = sub_ok["Effect_size"].values.astype(float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
    if mask.sum() >= 5:
        borzoi_vs_cbpnet[ct] = (x[mask], y[mask], c[mask])
        r_bvsc[ct], _ = stats.pearsonr(x[mask], y[mask])

    for snp_pos, cts in EXAMPLE_SNPS.items():
        if ct not in cts: continue
        color = EXAMPLE_COLORS[snp_pos]
        rows = merged[(merged["SNP_position"] == snp_pos) & (merged["cell_type"] == ct)]
        if len(rows):
            row = rows.iloc[0]
            hx  = float(row["Effect_size"])
            hy1 = float(row[cbp_col])
            hy2 = float(row[bor_col])
            if np.isfinite(hx) and np.isfinite(hy1):
                highlight_pts_cbpnet.setdefault(ct, []).append((hx, hy1, color))
            if np.isfinite(hx) and np.isfinite(hy2):
                highlight_pts_borzoi.setdefault(ct, []).append((hx, hy2, color))

# All 12 cell types, sorted by n (most data first)
ct_ordered = sorted(borzoi_vs_cbpnet, key=lambda c: -len(borzoi_vs_cbpnet[c][0]))
for ct in ct_ordered:
    print(f"  {CT_LABELS[ct]:12s}  cbpnet n={len(cbpnet_lead[ct][0]):,}"
          f"  borzoi n={len(borzoi_lead[ct][0]):,}"
          f"  vs n={len(borzoi_vs_cbpnet[ct][0]):,}")

# ── Colormap for row 3 ────────────────────────────────────────────────────────
all_c = np.concatenate([borzoi_vs_cbpnet[ct][2] for ct in ct_ordered])
vmin  = np.percentile(all_c, 2)
vmax  = np.percentile(all_c, 98)
cmap  = mcolors.LinearSegmentedColormap.from_list(
    "effect_div",
    [(0.0,"#2166AC"),(0.5,"#D8D8D8"),(1.0,"#B2182B")], N=256)
norm  = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.5, vmax=vmax)

# ── Layout: 3 rows × 12 scatter + colorbar ───────────────────────────────────
NCOLS = len(ct_ordered)
fig = plt.figure(figsize=(11.0, 4.29))
gs  = gridspec.GridSpec(3, NCOLS + 1,
                        hspace=0.55, wspace=0.35,
                        height_ratios=[1, 1, 1],
                        width_ratios=[1]*NCOLS + [0.04])

# ── Row 1: ChromBPNet logfc vs RASQUAL ───────────────────────────────────────
for ci, ct in enumerate(ct_ordered):
    ax = fig.add_subplot(gs[0, ci])
    x, y = cbpnet_lead[ct]
    r    = r_cbpnet.get(ct, float("nan"))

    ax.scatter(x, y, s=DOT_SIZE, color=RASQUAL_COLOR, alpha=DOT_ALPHA,
               linewidths=0, rasterized=True, zorder=2)
    add_ols(ax, x, y)
    add_ref_lines(ax, 0.5, 0)

    if ct in highlight_pts_cbpnet:
        from itertools import groupby as _gb
        pts = highlight_pts_cbpnet[ct]
        for hcolor, grp in _gb(sorted(pts, key=lambda t: t[2]), key=lambda t: t[2]):
            gpts = list(grp)
            ax.scatter([p[0] for p in gpts], [p[1] for p in gpts],
                       s=3, color=hcolor, edgecolors="none", linewidths=0, zorder=5)

    ax.set_title(f"{CT_LABELS[ct]}\nr={r:.2f}, n={len(x):,}",
                 fontsize=5, pad=2.0, linespacing=1.4)
    ax.set_xlabel("RASQUAL ES", fontsize=5, labelpad=2)
    ax.set_ylabel("ChromBPNet\nlogfc" if ci == 0 else "", fontsize=5, labelpad=2)
    style_ax(ax)

# ── Row 2: Cerberus logSUM vs RASQUAL ──────────────────────────────────────────
for ci, ct in enumerate(ct_ordered):
    ax = fig.add_subplot(gs[1, ci])
    x, y = borzoi_lead[ct]
    r    = r_borzoi.get(ct, float("nan"))

    ax.scatter(x, y, s=DOT_SIZE, color=RASQUAL_COLOR, alpha=DOT_ALPHA,
               linewidths=0, rasterized=True, zorder=2)
    add_ols(ax, x, y)
    add_ref_lines(ax, 0.5, 0)

    if ct in highlight_pts_borzoi:
        from itertools import groupby as _gb
        pts = highlight_pts_borzoi[ct]
        for hcolor, grp in _gb(sorted(pts, key=lambda t: t[2]), key=lambda t: t[2]):
            gpts = list(grp)
            ax.scatter([p[0] for p in gpts], [p[1] for p in gpts],
                       s=3, color=hcolor, edgecolors="none", linewidths=0, zorder=5)

    ax.set_title(f"{CT_LABELS[ct]}\nr={r:.2f}, n={len(x):,}",
                 fontsize=5, pad=2.0, linespacing=1.4)
    ax.set_xlabel("RASQUAL ES", fontsize=5, labelpad=2)
    ax.set_ylabel("Cerberus\nlogSUM" if ci == 0 else "", fontsize=5, labelpad=2)
    style_ax(ax)

# ── Row 3: Cerberus logSUM vs ChromBPNet logfc ─────────────────────────────────
for ci, ct in enumerate(ct_ordered):
    ax = fig.add_subplot(gs[2, ci])
    x, y, c = borzoi_vs_cbpnet[ct]
    r        = r_bvsc.get(ct, float("nan"))

    ax.scatter(x, y, s=DOT_SIZE, c=c, cmap=cmap, norm=norm,
               alpha=0.4, linewidths=0, rasterized=True, zorder=2)
    add_ols(ax, x, y)
    add_ref_lines(ax, 0, 0)

    ax.set_title(f"{CT_LABELS[ct]}\nr={r:.2f}, n={len(x):,}",
                 fontsize=5, pad=2.0, linespacing=1.4)
    ax.set_xlabel("ChromBPNet\nlogfc", fontsize=5, labelpad=2)
    ax.set_ylabel("Cerberus\nlogSUM" if ci == 0 else "", fontsize=5, labelpad=2)
    style_ax(ax)

# ── Colorbar for row 3 ────────────────────────────────────────────────────────
cbar_ax = fig.add_subplot(gs[2, NCOLS])
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, cax=cbar_ax)
cb.set_label("Effect size", fontsize=5, labelpad=3)
cb.ax.tick_params(labelsize=4, width=0.4, length=2)
cb.outline.set_linewidth(0.4)

# ── Save ──────────────────────────────────────────────────────────────────────
fig.savefig(OUT_PDF, dpi=600, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")
