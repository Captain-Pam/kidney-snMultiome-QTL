"""
Cell-type specificity of ChromBPNet and Cerberus caQTL predictions,
stratified by whether the variant falls within the caQTL peak.

Heatmap: rows = caQTL cell type (SuSiE PIP > 0.5 variants),
         cols = model cell type,
         value = fold-enrichment over null background (max PIP < 0.01 in ALL CTs).

Layout: 2 × 2
  row 0: variant IN the caQTL peak    (ChromBPNet | Cerberus)
  row 1: variant OUTSIDE the caQTL peak (ChromBPNet | Cerberus)

Output: <FIGURE_DIR>/caqtl_celltype_specificity.pdf
"""

import glob, os
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
CBPNET_DIR = os.path.join(cfg("WORK_DIR"), "chrombpnet", "caqtl", "scores_ensemble")
BORZOI_DIR = os.path.join(cfg("WORK_DIR"), "caqtl", "scores", "logSUM", "local_ensemble")
SUSIE_PATH = cfg("CAQTL_SUSIE_TSV")
OUT_PDF    = os.path.join(cfg("FIGURE_DIR"), "caqtl_celltype_specificity.pdf")

# ── Cell types & Cerberus track indices ─────────────────────────────────────────
CT_ORDER = ["PTS", "injPT", "TAL", "DTL_ATL", "DCT", "CNT_CD_PC",
            "IC", "EC", "Podocyte", "PEC", "Stromal", "Immune"]
CT_LABELS = {
    "PTS": "PTS", "injPT": "injPT", "TAL": "TAL", "DTL_ATL": "DTL/ATL",
    "DCT": "DCT", "CNT_CD_PC": "CNT/CD/PC", "IC": "IC", "EC": "EC",
    "Podocyte": "Podocyte", "PEC": "PEC", "Stromal": "Stromal", "Immune": "Immune",
}
CT_LAB = [CT_LABELS[ct] for ct in CT_ORDER]

BORZOI_ATAC_IDX = {
    "CNT_CD_PC": 0,  "DCT":  2,  "DTL_ATL": 4,  "EC":      6,
    "IC":         8,  "Immune": 10, "PEC":  12,  "PTS":    14,
    "Podocyte":  16,  "Stromal": 18, "TAL": 20,  "injPT":  22,
}

PIP_THRESH = 0.5
BG_PIP_MAX = 0.01   # background: max PIP < this across ALL CTs

# ── Load SuSiE PIPs ────────────────────────────────────────────────────────────
print("Loading SuSiE PIPs …")
susie = pd.read_csv(SUSIE_PATH, sep="\t",
                    usecols=["celltype", "variant_id", "variable_prob",
                             "peak_id", "position"])
susie = susie.rename(columns={"celltype": "cell_type", "variable_prob": "pip"})

# Determine whether variant falls within peak boundaries
pk = susie["peak_id"].str.rsplit("-", n=2, expand=True)
susie["peak_start"] = pk[1].astype(int)
susie["peak_end"]   = pk[2].astype(int)
susie["in_peak"] = ((susie["position"] >= susie["peak_start"]) &
                    (susie["position"] <= susie["peak_end"]))

# Background: variants with max PIP < BG_PIP_MAX across all CT × peak entries
max_pip_per_var = susie.groupby("variant_id")["pip"].max()
bg_variants = set(max_pip_per_var.index[max_pip_per_var < BG_PIP_MAX])
print(f"  Background variants (max PIP < {BG_PIP_MAX} in all CTs): {len(bg_variants):,}")

fg = susie[susie["pip"] > PIP_THRESH]
print(f"  PIP > {PIP_THRESH}: {len(fg):,} variant×CT pairs")
print(f"    in peak:     {fg['in_peak'].sum():,}")
print(f"    outside peak:{(~fg['in_peak']).sum():,}")

# Filter CT_ORDER: keep only CTs with >= 10 variants at PIP > 0.9
# (matched to ChromBPNet scored variants — proxy: count from SuSiE joined with cbp_wide later;
#  use SuSiE counts as pre-filter since cbp coverage is ~100%)
n_high09 = (susie[susie["pip"] > 0.9]
            .groupby("cell_type")["variant_id"]
            .nunique())
CT_ORDER = [ct for ct in CT_ORDER if n_high09.get(ct, 0) >= 10]
CT_LAB   = [CT_LABELS[ct] for ct in CT_ORDER]
print(f"\n  CTs retained (≥10 variants with PIP>0.9): {CT_ORDER}")
print(f"  CTs removed:  {[ct for ct in ['PTS','injPT','TAL','DTL_ATL','DCT','CNT_CD_PC','IC','EC','Podocyte','PEC','Stromal','Immune'] if ct not in CT_ORDER]}")

# ── Load ChromBPNet scores (all CTs → wide) ───────────────────────────────────
print("\nLoading ChromBPNet scores …")
cbp_wide = None
for ct in CT_ORDER:
    f = os.path.join(CBPNET_DIR, ct, f"{ct}.variant_scores.tsv")
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, sep="\t", usecols=["variant_id", "abs_logfc"])
    df = df.rename(columns={"abs_logfc": f"cbp_{ct}"}).set_index("variant_id")
    cbp_wide = df if cbp_wide is None else cbp_wide.join(df, how="outer")
cbp_wide = cbp_wide.reset_index()
print(f"  {len(cbp_wide):,} variants × {len(CT_ORDER)} CT scores")

# Background scores (subset of cbp_wide)
cbp_bg = cbp_wide[cbp_wide["variant_id"].isin(bg_variants)]
bg_means_cbp = {ct: cbp_bg[f"cbp_{ct}"].mean() for ct in CT_ORDER}
print(f"  Background variants in ChromBPNet: {len(cbp_bg):,}")

# ── Load Cerberus scores ─────────────────────────────────────────────────────────
print("\nLoading Cerberus scores …")
chunks = sorted(glob.glob(os.path.join(BORZOI_DIR, "chunk_*/scores.h5")),
                key=lambda p: int(p.split("chunk_")[1].split("/")[0]))

sorted_tracks = sorted(BORZOI_ATAC_IDX[ct] for ct in CT_ORDER)
track_to_col  = {t: i for i, t in enumerate(sorted_tracks)}
snp_ids_all, logsum_all = [], []
for h5_path in chunks:
    with h5py.File(h5_path, "r") as h5:
        snp_ids_all.append(h5["snp"][:].astype(str))
        logsum_all.append(h5["cov/logSUM"][:, sorted_tracks].astype(np.float32))

snp_ids = np.concatenate(snp_ids_all)
logsum  = np.concatenate(logsum_all)   # (N, 12)
borzoi_idx = pd.Series(np.arange(len(snp_ids)), index=snp_ids)
borzoi_idx = borzoi_idx[~borzoi_idx.index.duplicated(keep="first")]
print(f"  {len(snp_ids):,} variants, {logsum.shape[1]} ATAC tracks")

# Background Cerberus
bg_in_bor = [v for v in bg_variants if v in borzoi_idx.index]
bg_row_idxs = borzoi_idx.loc[bg_in_bor].values
bg_means_bor = {
    ct: float(np.abs(logsum[bg_row_idxs, track_to_col[BORZOI_ATAC_IDX[ct]]]).mean())
    for ct in CT_ORDER
}
print(f"  Background variants in Cerberus: {len(bg_in_bor):,}")

# ── Build enrichment matrix for a given foreground subset ─────────────────────
def build_matrices(fg_susie):
    """fg_susie: rows from susie filtered to the foreground of interest."""
    # ChromBPNet
    cbp_fg = fg_susie.merge(cbp_wide, on="variant_id", how="inner")
    cbp_mat = np.full((len(CT_ORDER), len(CT_ORDER)), np.nan)
    cbp_ns  = np.zeros(len(CT_ORDER), dtype=int)
    for i, caqtl_ct in enumerate(CT_ORDER):
        sub = cbp_fg[cbp_fg["cell_type"] == caqtl_ct]
        cbp_ns[i] = len(sub)
        for j, model_ct in enumerate(CT_ORDER):
            vals = sub[f"cbp_{model_ct}"].dropna()
            if len(vals) >= 3:
                cbp_mat[i, j] = vals.mean() / bg_means_cbp[model_ct]

    # Cerberus
    fg_bor = fg_susie[fg_susie["variant_id"].isin(borzoi_idx.index)].copy()
    if len(fg_bor):
        row_idxs = borzoi_idx.loc[fg_bor["variant_id"].values].values
        for j, ct in enumerate(CT_ORDER):
            col_j = track_to_col[BORZOI_ATAC_IDX[ct]]
            fg_bor[f"bor_{ct}"] = np.abs(logsum[row_idxs, col_j])

    bor_mat = np.full((len(CT_ORDER), len(CT_ORDER)), np.nan)
    bor_ns  = np.zeros(len(CT_ORDER), dtype=int)
    for i, caqtl_ct in enumerate(CT_ORDER):
        sub = fg_bor[fg_bor["cell_type"] == caqtl_ct]
        bor_ns[i] = len(sub)
        for j, model_ct in enumerate(CT_ORDER):
            vals = sub[f"bor_{model_ct}"].dropna()
            if len(vals) >= 3:
                bor_mat[i, j] = vals.mean() / bg_means_bor[model_ct]

    return cbp_mat, bor_mat, cbp_ns, bor_ns

# Two foreground strata
high_pip = susie[susie["pip"] > PIP_THRESH]
strata = [
    (high_pip[high_pip["in_peak"]],       "in peak"),
    (high_pip[~high_pip["in_peak"]],      "outside peak"),
]

# ── Figure ─────────────────────────────────────────────────────────────────────
n_ct = len(CT_ORDER)
fig, axes = plt.subplots(2, 2, figsize=(11, 9.5),
                          gridspec_kw=dict(wspace=0.35, hspace=0.55,
                                           left=0.09, right=0.97,
                                           top=0.90, bottom=0.12))

def draw_heatmap(ax, mat, title, ns):
    # Row-wise z-score: within each caQTL CT row, centre at row mean
    # so colour shows which model CT scores highest for that caQTL CT
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        row_mean = np.nanmean(mat, axis=1, keepdims=True)
        row_std  = np.nanstd(mat,  axis=1, keepdims=True)
    row_std[~np.isfinite(row_std) | (row_std == 0)] = 1
    mat_z = (mat - row_mean) / row_std

    # Symmetric colour limits
    zlim = max(np.nanpercentile(np.abs(mat_z[np.isfinite(mat_z)]), 98), 1.0)
    im = ax.imshow(mat_z, cmap="RdYlBu_r", vmin=-zlim, vmax=zlim, aspect="auto")

    ax.set_xticks(range(n_ct))
    ax.set_xticklabels(CT_LAB, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_ct))
    ax.set_yticklabels([f"{lab} (n={ns[k]})" for k, lab in enumerate(CT_LAB)], fontsize=7)
    ax.set_xlabel("Model cell type", fontsize=8)
    ax.set_ylabel("caQTL cell type (SuSiE PIP > 0.5)", fontsize=8)
    ax.set_title(title, fontsize=9, pad=4)

    # Annotate diagonal and highest off-diagonal cell per row
    for k in range(n_ct):
        v = mat[k, k]
        if np.isfinite(v):
            ax.text(k, k, f"{v:.1f}×", ha="center", va="center",
                    fontsize=5, color="black")

        row_vals = mat[k].copy()
        row_vals[k] = np.nan
        j_max = int(np.nanargmax(row_vals)) if np.any(np.isfinite(row_vals)) else -1
        if j_max >= 0 and np.isfinite(row_vals[j_max]):
            ax.text(j_max, k, f"{row_vals[j_max]:.1f}×", ha="center", va="center",
                    fontsize=5, color="black")

    cb = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label("Row z-score\n(diagonal text = actual fold-enrichment)", fontsize=6.5)
    cb.ax.tick_params(labelsize=7)
    cb.ax.axhline(0, color="black", lw=0.8, ls="--")

for row, (fg, label) in enumerate(strata):
    print(f"\nBuilding matrices: PIP > {PIP_THRESH}, {label} …")
    cbp_mat, bor_mat, cbp_ns, bor_ns = build_matrices(fg)
    draw_heatmap(axes[row, 0], cbp_mat, f"ChromBPNet — {label}", cbp_ns)
    draw_heatmap(axes[row, 1], bor_mat,  f"Cerberus — {label}",     bor_ns)
    print("  Diagonal:")
    for i, ct in enumerate(CT_ORDER):
        print(f"    {ct:12s}  n={cbp_ns[i]:4d}  "
              f"CBP={cbp_mat[i,i]:.2f}×  Cerberus={bor_mat[i,i]:.2f}×")

n_bg_cbp = len(cbp_bg)
fig.suptitle(
    f"Cell-type specificity of caQTL predictions (SuSiE PIP > {PIP_THRESH})\n"
    f"Fold-enrichment over null variants (max PIP < 0.01 in all CTs; n={n_bg_cbp:,} null variants)\n"
    "stratified by whether the variant falls within the caQTL peak",
    fontsize=9.5, y=0.97)

fig.savefig(OUT_PDF, dpi=200, bbox_inches="tight")
print(f"\nSaved: {OUT_PDF}")
