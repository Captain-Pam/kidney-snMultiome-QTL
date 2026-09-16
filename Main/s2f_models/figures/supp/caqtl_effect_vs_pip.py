"""
ChromBPNet and Cerberus variant effect scores vs caQTL SuSiE finemapping PIP

Data sources:
  ChromBPNet: scores_ensemble/{CT}/{CT}.variant_scores.tsv  → abs_logfc
  Cerberus:     caqtl_scoring/scores/logSUM/local_ensemble/   → abs(logSUM) on ATAC track
  PIP:        5_susie/2_caQTL/caqtl_kgp_hg38_res_susie_merged_all_celltypes.tsv.gz
              column: variable_prob (SuSiE PIP per variant × cell type)

Output: <FIGURE_DIR>/caqtl_effect_vs_pip.pdf
"""

import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figlib import (  # noqa: E402
    ALT_COLOR, ATAC_IDX, CT_COLORS, CT_LABELS, NT_COLORS, REF_COLOR,
    add_ols, add_ref_lines, cfg, model_ct, qtl_ct, select_lead,
    setup_style, style_ax,
)

setup_style()

# ── Paths ──────────────────────────────────────────────────────────────────────
CBPNET_DIR  = os.path.join(cfg("WORK_DIR"), "chrombpnet", "caqtl", "scores_ensemble")
BORZOI_DIR  = os.path.join(cfg("WORK_DIR"), "caqtl", "scores", "logSUM", "local_ensemble")
SUSIE_PATH  = cfg("CAQTL_SUSIE_TSV")
OUT_PDF     = os.path.join(cfg("FIGURE_DIR"), "caqtl_effect_vs_pip.pdf")

# ── Cell types & Cerberus ATAC track indices ─────────────────────────────────────
CT_ORDER = ["CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC",
            "Immune", "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT"]

BORZOI_ATAC_IDX = {
    "CNT_CD_PC": 0,  "DCT":      2,  "DTL_ATL": 4,  "EC":      6,
    "IC":         8,  "Immune":  10,  "PEC":    12,  "PTS":    14,
    "Podocyte":  16,  "Stromal": 18,  "TAL":    20,  "injPT":  22,
}

# ── PIP bins ───────────────────────────────────────────────────────────────────
PIP_BINS = [0,   0.01, 0.1, 0.5, 0.9, 1.001]
PIP_LABS = ["0–0.01", "0.01–0.1", "0.1–0.5", "0.5–0.9", ">0.9"]
BIN_COLORS = ["#E8E8E8", "#C5D9EE", "#7EB3D8", "#3A78B5", "#1A3D6B"]

# ── Load SuSiE PIPs ────────────────────────────────────────────────────────────
print("Loading SuSiE PIPs …")
susie = pd.read_csv(SUSIE_PATH, sep="\t",
                    usecols=["celltype", "variant_id", "variable_prob"])
susie = susie.rename(columns={"celltype": "cell_type", "variable_prob": "pip"})
susie = susie.dropna(subset=["pip"])
print(f"  {len(susie):,} variant×CT rows")
print(f"  PIP > 0.5: {(susie.pip > 0.5).sum():,}")

# ── Load ChromBPNet scores ─────────────────────────────────────────────────────
print("Loading ChromBPNet scores …")
cbp_dfs = []
for ct in CT_ORDER:
    f = os.path.join(CBPNET_DIR, ct, f"{ct}.variant_scores.tsv")
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, sep="\t", usecols=["variant_id", "abs_logfc"])
    df["cell_type"] = ct
    cbp_dfs.append(df[["variant_id", "cell_type", "abs_logfc"]])

cbp = pd.concat(cbp_dfs, ignore_index=True)
print(f"  {len(cbp):,} scored variants")

cbp_merged = cbp.merge(susie, on=["variant_id", "cell_type"], how="inner")
cbp_merged = cbp_merged.dropna(subset=["abs_logfc", "pip"])
print(f"  Matched: {len(cbp_merged):,}")

# ── Load Cerberus scores ─────────────────────────────────────────────────────────
print("Loading Cerberus scores …")
chunks = sorted(glob.glob(os.path.join(BORZOI_DIR, "chunk_*/scores.h5")),
                key=lambda p: int(p.split("chunk_")[1].split("/")[0]))

snp_ids_all, logsum_all = [], []
for h5_path in chunks:
    with h5py.File(h5_path, "r") as h5:
        snp_ids_all.append(h5["snp"][:].astype(str))
        logsum_all.append(h5["cov/logSUM"][:].astype(np.float32))

snp_ids = np.concatenate(snp_ids_all)   # chr:pos:ref:alt
logsum  = np.concatenate(logsum_all)    # (N, 124)
print(f"  {len(snp_ids):,} variants, {logsum.shape[1]} tracks")

# Build variant_id → row index
borzoi_idx = pd.Series(np.arange(len(snp_ids)), index=snp_ids)
borzoi_idx = borzoi_idx[~borzoi_idx.index.duplicated(keep="first")]

print("Extracting per-CT Cerberus logSUM …")
bor_rows = []
for ct in CT_ORDER:
    if ct not in BORZOI_ATAC_IDX:
        continue
    track_idx = BORZOI_ATAC_IDX[ct]
    sub = susie[susie["cell_type"] == ct].copy()
    in_borzoi = sub["variant_id"].isin(borzoi_idx.index)
    matched = sub[in_borzoi].copy()
    if matched.empty:
        continue
    row_idxs = borzoi_idx.loc[matched["variant_id"].values].values
    matched["borzoi_abs_logsum"] = np.abs(logsum[row_idxs, track_idx])
    bor_rows.append(matched[["variant_id", "cell_type", "pip", "borzoi_abs_logsum"]])

bor_merged = pd.concat(bor_rows, ignore_index=True).dropna()
print(f"  Matched: {len(bor_merged):,}")

# ── Bin by PIP ─────────────────────────────────────────────────────────────────
cbp_merged["pip_bin"] = pd.cut(cbp_merged["pip"], bins=PIP_BINS,
                                labels=PIP_LABS, include_lowest=True)
bor_merged["pip_bin"] = pd.cut(bor_merged["pip"], bins=PIP_BINS,
                                labels=PIP_LABS, include_lowest=True)

cbp_groups = [cbp_merged.loc[cbp_merged["pip_bin"] == lb, "abs_logfc"].values
              for lb in PIP_LABS]
bor_groups = [bor_merged.loc[bor_merged["pip_bin"] == lb, "borzoi_abs_logsum"].values
              for lb in PIP_LABS]
cbp_ns = [len(g) for g in cbp_groups]
bor_ns = [len(g) for g in bor_groups]

# Statistics
cbp_rho, cbp_p = stats.spearmanr(cbp_merged["pip"], cbp_merged["abs_logfc"])
bor_rho, bor_p = stats.spearmanr(bor_merged["pip"], bor_merged["borzoi_abs_logsum"])
cbp_kw_p = stats.kruskal(*[g for g in cbp_groups if len(g) > 1])[1]
bor_kw_p = stats.kruskal(*[g for g in bor_groups if len(g) > 1])[1]

print(f"ChromBPNet  ρ={cbp_rho:.3f} p={cbp_p:.2e}  KW p={cbp_kw_p:.2e}")
print(f"Cerberus      ρ={bor_rho:.3f} p={bor_p:.2e}  KW p={bor_kw_p:.2e}")
for lb, g in zip(PIP_LABS, cbp_groups):
    print(f"  CBP {lb:10s}  n={len(g):6d}  median={np.median(g):.4f}" if len(g) else f"  CBP {lb} n=0")

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(6.0, 3.2),
                          gridspec_kw=dict(wspace=0.38, left=0.11, right=0.97,
                                           top=0.86, bottom=0.22))

def draw_boxplot(ax, groups, ns, title, ylabel, rho, kw_p):
    bp = ax.boxplot(groups, positions=range(len(PIP_LABS)),
                    widths=0.55, patch_artist=True,
                    medianprops=dict(color="black", lw=1.5),
                    whiskerprops=dict(lw=0.8, color="#555"),
                    capprops=dict(lw=0.8, color="#555"),
                    flierprops=dict(marker=".", markersize=1.5, alpha=0.3,
                                   markeredgewidth=0, color="#888"),
                    showfliers=True)
    for patch, c in zip(bp["boxes"], BIN_COLORS):
        patch.set_facecolor(c)
        patch.set_alpha(0.9)
    ax.set_xticks(range(len(PIP_LABS)))
    ax.set_xticklabels(PIP_LABS, rotation=30, ha="right", fontsize=7.5)
    ax.set_xlabel("SuSiE PIP", fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_title(f"{title}\nρ={rho:.2f}  KW p={kw_p:.1e}", fontsize=8.5, pad=3)
    ymax = ax.get_ylim()[1]
    for i, n in enumerate(ns):
        ax.text(i, ymax * 0.97, f"n={n:,}", ha="center", va="top",
                fontsize=5.5, color="#444")

draw_boxplot(axes[0], cbp_groups, cbp_ns,
             "ChromBPNet", "|log FC|", cbp_rho, cbp_kw_p)
draw_boxplot(axes[1], bor_groups, bor_ns,
             "Cerberus", "|log SUM|", bor_rho, bor_kw_p)

fig.suptitle("Model variant effect vs caQTL finemapping PIP (SuSiE)", fontsize=10, y=0.97)
fig.savefig(OUT_PDF, dpi=200, bbox_inches="tight")
print(f"Saved: {OUT_PDF}")
