#!/usr/bin/env python3
"""A proximal-tubule caQTL disrupting a cell-type-specific HNF1A motif.

chr16:89,622,209 A>G, shown across six kidney cell types. Left: observed
pseudobulk ATAC coverage stratified by donor genotype. Middle: per-donor in-peak
fragment RPM by genotype. Right: ChromBPNet contribution scores expressed
relative to the across-cell-type mean, which isolates the cell-type-differential
component of the attribution, with the HNF1A motif beneath.

Writes several figures; the combined one is the main panel.

Needs donor genotypes and per-donor fragment files, both controlled access.
See the README.

Usage:
    python chr16_caqtl_locus.py [--out_dir <dir>]
"""
import argparse
import os, re, subprocess
from multiprocessing import Pool
import numpy as np
import pandas as pd
import pysam
import h5py
import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

from figlib import (
    ALT_COLOR, NT_COLORS, REF_COLOR,
    align_pwm, cfg, parse_pwm, setup_style,
)

setup_style()

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--out_dir", default=None)
_args = _parser.parse_args()

# ── Config ────────────────────────────────────────────────────────────────────
WORK_DIR = cfg("WORK_DIR")
OUT_DIR  = _args.out_dir or cfg("FIGURE_DIR")

# Donor genotypes and per-donor fragments: controlled access.
VCF_FILE = cfg("GENOTYPE_VCF")
FRAG_DIR = cfg("FRAG_DIR_INDIVIDUAL")

# ChromBPNet attributions at this variant, per cell type.
CBPNET_DIR = os.path.join(WORK_DIR, "chrombpnet", "caqtl", "shap")
CBPNET_ALL_CTS_DIR = os.path.join(CBPNET_DIR, "chr16_all_cts")

# Cerberus ISM, averaged over folds where available.
_ISM_DIR      = os.path.join(WORK_DIR, "caqtl", "chr16_ism")
_BOR_ENSEMBLE = os.path.join(_ISM_DIR, "ensemble", "scores.h5")
_BOR_SINGLE   = os.path.join(_ISM_DIR, "scores.h5")
BORZOI_ISM_H5 = _BOR_ENSEMBLE if os.path.exists(_BOR_ENSEMBLE) else _BOR_SINGLE

VKEY      = "chr16_89622209_A_G"
CHR       = "chr16"
VAR_POS   = 89622209
PEAK_S    = 89621958   # chr16-89621958-89622460 (lead peak, q=2.33e-09)
PEAK_E    = 89622460

COV_HALF  = 5000          # ±5 kb coverage window
BIN_SIZE  = 10            # 10 bp bins (matches original downsample-by-10)
ISM_HALF  = 25            # ±25 bp ISM logos

COV_START = VAR_POS - COV_HALF
COV_END   = VAR_POS + COV_HALF
N_BINS    = (COV_END - COV_START) // BIN_SIZE   # 1000

# Our 6 CTs; fragment files use PTS not PT
CELL_TYPES = ["PT", "injPT", "CNT_CD_PC", "TAL", "EC", "IC"]
CT_LABELS  = {"PT": "PTS", "injPT": "injPT", "CNT_CD_PC": "CNT/CD/PC",
              "TAL": "TAL",  "EC": "EC",      "IC": "IC"}
CT_FRAG    = {"PT": "PTS",  "injPT": "injPT", "CNT_CD_PC": "CNT_CD_PC",
              "TAL": "TAL",  "EC": "EC",       "IC": "IC"}
BORZOI_IDX = {"CNT_CD_PC": 0, "EC": 3, "IC": 4, "TAL": 8, "injPT": 9, "PT": 10}

GENO_COLORS = {"ref": "#EF3B2C", "het": "#4292C6", "alt": "#084594"}
GENO_ORDER  = ["ref", "het", "alt"]   # back→front order for area overlap

MOTIF_ID = "MA0046.2"
MOTIF_NM = "HNF1A"

# ── 1. Genotypes ──────────────────────────────────────────────────────────────
print("Loading genotypes …")
vcf = pysam.VariantFile(VCF_FILE)
geno_groups = {"ref": [], "het": [], "alt": []}
for rec in vcf.fetch(CHR, VAR_POS - 1, VAR_POS):
    for s in vcf.header.samples:
        gt = rec.samples[s]["GT"]
        if gt == (0, 0):        geno_groups["ref"].append(s)
        elif set(gt) == {0, 1}: geno_groups["het"].append(s)
        elif gt == (1, 1):      geno_groups["alt"].append(s)
vcf.close()
for g, ss in geno_groups.items():
    print(f"  {g}: {len(ss)} samples")

# ── 2. RPM: precompute total fragment counts per file ─────────────────────────
TOTAL_COUNTS_CACHE = os.path.join(OUT_DIR, "frag_total_counts.csv")

def _count_file(fp):
    """Sum count column (col 5) across all rows of one fragment file."""
    try:
        result = subprocess.run(
            f"zcat {fp} | awk '{{s+=$5}} END {{print s}}'",
            shell=True, capture_output=True, text=True, timeout=120
        )
        return fp, int(result.stdout.strip())
    except Exception:
        return fp, 0

def load_total_counts():
    """Load or compute total fragment counts for all relevant files."""
    files = [
        os.path.join(FRAG_DIR, f"{CT_FRAG[ct]}_{s}_fragments.tsv.gz")
        for ct in CELL_TYPES
        for samples in geno_groups.values()
        for s in samples
    ]
    files = [f for f in files if os.path.exists(f)]
    files = list(set(files))

    if os.path.exists(TOTAL_COUNTS_CACHE):
        cached = pd.read_csv(TOTAL_COUNTS_CACHE, index_col=0)["total"].to_dict()
        missing = [f for f in files if f not in cached]
    else:
        cached, missing = {}, files

    if missing:
        print(f"  Computing total counts for {len(missing)} files (parallel) …")
        with Pool(16) as pool:
            results = pool.map(_count_file, missing)
        for fp, cnt in results:
            cached[fp] = cnt
        df = pd.DataFrame(list(cached.items()), columns=["file", "total"]).set_index("file")
        df.to_csv(TOTAL_COUNTS_CACHE)
        print(f"  Cached to {TOTAL_COUNTS_CACHE}")

    return cached

print("Loading/computing RPM totals …")
total_counts = load_total_counts()

# ── 3. Coverage from fragment files (RPM normalized) ──────────────────────────
def sample_coverage(ct_frag, sample, chrom, start, end, bin_size):
    """RPM-normalised per-bin coverage; returns array or None."""
    fp = os.path.join(FRAG_DIR, f"{ct_frag}_{sample}_fragments.tsv.gz")
    if not os.path.exists(fp):
        return None
    total = total_counts.get(fp, 0)
    if total == 0:
        return None
    n_bins = (end - start) // bin_size
    cov = np.zeros(n_bins)
    try:
        tbx = pysam.TabixFile(fp)
        for row in tbx.fetch(chrom, start, end):
            cols = row.split("\t")
            fs, fe = int(cols[1]), int(cols[2])
            cnt = int(cols[4]) if len(cols) > 4 else 1
            b0 = max(0, (fs - start) // bin_size)
            b1 = min(n_bins, (fe - start + bin_size - 1) // bin_size)
            cov[b0:b1] += cnt
        tbx.close()
    except (ValueError, KeyError):
        pass
    return cov / total * 1_000_000   # RPM

def group_coverage(ct, samples):
    """Average RPM coverage across samples; returns (array_or_None, n_valid)."""
    ct_frag = CT_FRAG[ct]
    cov_sum = np.zeros(N_BINS)
    n = 0
    for s in samples:
        c = sample_coverage(ct_frag, s, CHR, COV_START, COV_END, BIN_SIZE)
        if c is not None:
            cov_sum += c
            n += 1
    return ((cov_sum / n), n) if n > 0 else (None, 0)

# ── 2b. Fragment-based coverage for arbitrary window (for combined figure) ────
def group_coverage_bw(ct, samples, chrom, start, end, n_bins):
    """Same RPM normalization as sample_coverage(), but for an arbitrary window."""
    bin_size = (end - start) // n_bins
    ct_frag = CT_FRAG[ct]
    cov_sum = np.zeros(n_bins)
    n = 0
    for s in samples:
        c = sample_coverage(ct_frag, s, chrom, start, end, bin_size)
        if c is not None:
            cov_sum += c
            n += 1
    return (cov_sum / n, n) if n > 0 else (None, 0)

print("Computing coverage …")
bin_centers = COV_START + (np.arange(N_BINS) + 0.5) * BIN_SIZE
cov = {}      # cov[ct][geno]  = array or None
cov_n = {}    # cov_n[ct][geno] = n samples contributing
for ct in CELL_TYPES:
    cov[ct] = {}; cov_n[ct] = {}
    for geno, samples in geno_groups.items():
        arr, n = group_coverage(ct, samples)
        cov[ct][geno]   = arr
        cov_n[ct][geno] = n
        mx = f"{arr.max():.2f}" if arr is not None else "n/a"
        print(f"  {CT_LABELS[ct]:10s} {geno:3s}: n={n}  max={mx}")

# ── 3. Per-sample peak RPM (for boxplot) ─────────────────────────────────────
print("Computing peak RPM per sample …")
peak_counts = {}
for ct in CELL_TYPES:
    ct_frag = CT_FRAG[ct]
    peak_counts[ct] = {g: [] for g in GENO_ORDER}
    for geno, samples in geno_groups.items():
        for s in samples:
            c = sample_coverage(ct_frag, s, CHR, PEAK_S, PEAK_E, BIN_SIZE)
            if c is not None:
                peak_counts[ct][geno].append(c.sum())

# ── 4. ChromBPNet SHAP ────────────────────────────────────────────────────────
print("Loading ChromBPNet SHAP …")
cbpnet = {}
ref_oh_any = None
for ct in CELL_TYPES:
    p = os.path.join(CBPNET_DIR, f"{ct}_ism.h5")
    if not os.path.exists(p):
        cbpnet[ct] = None; continue
    with h5py.File(p) as f:
        ids = list(f["variant_ids"][:].astype(str))
        if VKEY not in ids: cbpnet[ct] = None; continue
        i = ids.index(VKEY)
        delta  = f["ref_contrib"][i].mean(axis=-1) - f["alt_contrib"][i].mean(axis=-1)
        ref_oh = f["ref_onehot"][i]
    cbpnet[ct] = pd.DataFrame(ref_oh * delta[:, np.newaxis], columns=list("ACGT"))
    if ref_oh_any is None: ref_oh_any = ref_oh
    print(f"  {CT_LABELS[ct]}: OK")

# ── 5. Cerberus ISM ─────────────────────────────────────────────────────────────
print(f"Loading Cerberus ISM … ({os.path.basename(BORZOI_ISM_H5)})")
borzoi_ism = {}
with h5py.File(BORZOI_ISM_H5) as h5:
    labels   = list(h5["label"][:].astype(str))
    ref_sc   = h5["ref/cov/logSUM"][:].astype(np.float32)
    alt_sc   = h5["alt/cov/logSUM"][:].astype(np.float32)
    ref_seqs = h5["ref/seqs"][:].astype(np.int8)
vi = labels.index(VKEY)
r_sc, a_sc, r_seq = ref_sc[vi], alt_sc[vi], ref_seqs[vi]
for ct in CELL_TYPES:
    tidx = BORZOI_IDX.get(ct)
    if tidx is None: borzoi_ism[ct] = None; continue
    df = pd.DataFrame(0.0, index=range(50), columns=list("ACGT"))
    for pos in range(50):
        ridx   = int(np.argmax(r_seq[:, pos]))
        others = [j for j in range(4) if j != ridx]
        df.iloc[pos, ridx] = float(-np.mean(a_sc[pos, others, tidx] - r_sc[pos, others, tidx]))
    borzoi_ism[ct] = df
    print(f"  {CT_LABELS[ct]}: sum_abs={df.abs().values.sum():.1f}")

# ── 5b. Difference from cross-CT mean attribution ─────────────────────────────
def mean_logo(logo_dict):
    vals = [v.values for v in logo_dict.values() if v is not None]
    if not vals: return None
    return pd.DataFrame(np.mean(vals, axis=0), columns=list("ACGT"))

# ChromBPNet SHAP: mean across all CTs in chr16_all_cts/ (one h5 per CT,
# all folds averaged). Falls back to the 6 plotted CTs if not yet generated.
cbpnet_all = {}
if os.path.isdir(CBPNET_ALL_CTS_DIR):
    for fname in sorted(os.listdir(CBPNET_ALL_CTS_DIR)):
        if not fname.endswith("_shap.h5"): continue
        ct_key = fname[:-len("_shap.h5")]
        p = os.path.join(CBPNET_ALL_CTS_DIR, fname)
        with h5py.File(p) as f:
            ids = list(f["variant_ids"][:].astype(str))
            if VKEY not in ids: continue
            i = ids.index(VKEY)
            delta  = f["ref_contrib"][i].mean(axis=-1) - f["alt_contrib"][i].mean(axis=-1)
            ref_oh = f["ref_onehot"][i]
        cbpnet_all[ct_key] = pd.DataFrame(ref_oh * delta[:, np.newaxis], columns=list("ACGT"))
    print(f"ChromBPNet mean: loaded {len(cbpnet_all)} CTs from chr16_all_cts/")
if cbpnet_all:
    cbpnet_mean = mean_logo(cbpnet_all)
else:
    print("chr16_all_cts/ not ready — using 6 plotted CTs for mean")
    cbpnet_mean = mean_logo(cbpnet)

# Cerberus ISM: mean across ALL tracks in the h5, not just the 6 plotted CTs
n_tracks = r_sc.shape[-1]
borzoi_mean_arr = np.zeros((50, 4))
for tidx in range(n_tracks):
    for pos in range(50):
        ridx   = int(np.argmax(r_seq[:, pos]))
        others = [j for j in range(4) if j != ridx]
        borzoi_mean_arr[pos, ridx] += float(
            -np.mean(a_sc[pos, others, tidx] - r_sc[pos, others, tidx]))
borzoi_mean = pd.DataFrame(borzoi_mean_arr / n_tracks, columns=list("ACGT"))
print(f"Cerberus mean computed over {n_tracks} tracks")

cbpnet_diff = {
    ct: (cbpnet[ct] - cbpnet_mean
         if cbpnet[ct] is not None and cbpnet_mean is not None else None)
    for ct in CELL_TYPES
}
borzoi_diff = {
    ct: (borzoi_ism[ct] - borzoi_mean
         if borzoi_ism[ct] is not None and borzoi_mean is not None else None)
    for ct in CELL_TYPES
}

# ── 6. PWM ────────────────────────────────────────────────────────────────────
raw_pwm = parse_pwm(MOTIF_ID)
if ref_oh_any is not None:
    pwm_df, is_rc, bpos = align_pwm(
        ref_oh_any, raw_pwm, var_pos=None, n_bases=50, ic_mode="simple")
    print(f"PWM: {MOTIF_NM} → pos {bpos}, RC={is_rc}")
else:
    pwm_df = pd.DataFrame(np.zeros((50, 4)), columns=list("ACGT"))
    is_rc = False

# ── 7. Y-limits ───────────────────────────────────────────────────────────────
def ylim_for(logos):
    vals = [v for v in logos.values() if v is not None]
    if not vals:
        return (0, 0.01)
    all_vals = np.concatenate([v.values for v in vals])
    ylo, yhi = all_vals.min(), all_vals.max()
    margin = (yhi - ylo) * 0.12 + 0.001
    return (ylo - margin if ylo < 0 else 0, yhi + margin)

cbp_ylim      = ylim_for(cbpnet)
bor_ylim      = ylim_for(borzoi_ism)
cbp_diff_ylim = ylim_for(cbpnet_diff)
bor_diff_ylim = ylim_for(borzoi_diff)

cov_ymax = 0
for ct in CELL_TYPES:
    for g in GENO_ORDER:
        y = cov[ct][g]
        if y is not None: cov_ymax = max(cov_ymax, y.max())
cov_ymax *= 1.1

# per-CT y-limits (for second figure)
cov_ymax_ct = {}
for ct in CELL_TYPES:
    mx = max((cov[ct][g].max() for g in GENO_ORDER if cov[ct][g] is not None), default=0)
    cov_ymax_ct[ct] = mx * 1.1 if mx > 0 else 1.0

# ── 8. Drawing helpers ────────────────────────────────────────────────────────
def _ax_style(ax, ytick, labelsize=5.5):
    ax.set_xticks([])
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.tick_params(axis="y", labelsize=labelsize, width=0.4, length=1.5, pad=1)
    ax.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
    if not ytick:
        ax.set_yticklabels([]); ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)

def draw_ism(ax, ldf, ylim, ytick, labelsize=5.5, na_fs=5, var_line=False):
    ylo, yhi = ylim
    if ldf is not None:
        logomaker.Logo(ldf, color_scheme=NT_COLORS, ax=ax,
                       flip_below=True, baseline_width=0, font_name="DejaVu Sans")
        for p in ax.patches: p.set_zorder(3)
        ax.axhline(0, color="#bbb", lw=0.35, zorder=1)
    else:
        ax.text(0.5, 0.5, "n/a", transform=ax.transAxes,
                ha="center", va="center", fontsize=na_fs, color="#ccc")
    if var_line:
        ax.axvline(25, color="#555", lw=0.6, ls="--", zorder=4, alpha=0.8)
    else:
        ax.axvspan(24.5, 25.5, color="#d0e8ff", alpha=0.55, zorder=0, lw=0)
    ax.set_xlim(-0.5, 49.5); ax.set_ylim(ylo, yhi)
    _ax_style(ax, ytick, labelsize=labelsize)

def draw_cov(ax, x, cov_dict, n_dict, ymax, ytick):
    has = False
    for geno in GENO_ORDER:
        y = cov_dict.get(geno)
        if y is None: continue
        has = True
        c = GENO_COLORS[geno]
        ax.fill_between(x, 0, y, color=c, alpha=0.45, linewidth=0)
        ax.plot(x, y, color=c, linewidth=0.5)
    if not has:
        ax.text(0.5, 0.5, "n/a", transform=ax.transAxes,
                ha="center", va="center", fontsize=5, color="#ccc")
    ax.axvspan(PEAK_S, PEAK_E, color="#e8e8e8", alpha=0.4, zorder=0)
    ax.axvline(VAR_POS, color="#555", lw=0.6, ls="--", zorder=4, alpha=0.8)
    ax.set_xlim(COV_START, COV_END)
    ax.set_ylim(0, ymax)
    _ax_style(ax, ytick)
    # per-row n legend (always shown)
    GT_LABEL = {"ref": "0/0", "het": "0/1", "alt": "1/1"}
    patches = [mpatches.Patch(facecolor=GENO_COLORS[g], alpha=0.7,
                              label=f"{GT_LABEL[g]} (n={n_dict.get(g, 0)})")
               for g in GENO_ORDER if n_dict.get(g, 0) > 0]
    ax.legend(handles=patches, fontsize=4, frameon=False,
              loc="upper right", handlelength=1,
              labelspacing=0.1, handletextpad=0.3, borderpad=0.2)

# ── 9. Figure: tracks + ISM ───────────────────────────────────────────────────
N_CT      = len(CELL_TYPES)
H_RATIOS  = [1.0] * N_CT + [0.5]   # CT rows + PWM row

fig = plt.figure(figsize=(10, 4.0))
fig.suptitle("chr16:89622209:A:G  (HNF1A caQTL)",
             fontsize=7.5, fontweight="bold", y=1.01)

outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.22,
                          width_ratios=[2, 3],
                          top=0.93, bottom=0.04, left=0.08, right=0.99)

# — ISM section (ChromBPNet | Cerberus) —
ism_outer = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], wspace=0.06)
for mi, (logo_dict, yabs_v, lbl) in enumerate([
    (cbpnet_diff,   cbp_diff_ylim, "ChromBPNet\nSHAP"),
    (borzoi_diff,   bor_diff_ylim, "Cerberus\nISM"),
]):
    gs = gridspec.GridSpecFromSubplotSpec(
        N_CT + 1, 1, subplot_spec=ism_outer[mi],
        height_ratios=H_RATIOS, hspace=0.06)
    for ri, ct in enumerate(CELL_TYPES):
        ax = fig.add_subplot(gs[ri])
        draw_ism(ax, logo_dict.get(ct), yabs_v, ytick=True)
        if mi == 0:
            ax.set_ylabel(CT_LABELS[ct], fontsize=6, rotation=0,
                          ha="right", va="center", labelpad=34)
        if ri == 0:
            ax.annotate(lbl, xy=(0.5, 1), xycoords="axes fraction",
                        xytext=(0, 16), textcoords="offset points",
                        fontsize=5, color="#555", ha="center", va="bottom",
                        annotation_clip=False)
    # PWM row
    ax_pwm = fig.add_subplot(gs[N_CT])
    if pwm_df is not None:
        logomaker.Logo(pwm_df, color_scheme=NT_COLORS, ax=ax_pwm,
                       flip_below=False, baseline_width=0, font_name="DejaVu Sans")
        for p in ax_pwm.patches: p.set_zorder(3)
    ax_pwm.axvspan(24.5, 25.5, color="#d0e8ff", alpha=0.55, zorder=0, lw=0)
    ax_pwm.set_xlim(-0.5, 49.5); ax_pwm.set_xticks([])
    ax_pwm.spines[["top", "right", "bottom"]].set_visible(False)
    ax_pwm.spines["left"].set_linewidth(0.4)
    ax_pwm.tick_params(axis="y", labelsize=5.5, width=0.4, length=1.5, pad=1)
    ax_pwm.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
    if mi == 0:
        ax_pwm.set_ylabel("PWM\n(IC)", fontsize=5, rotation=0,
                          ha="right", va="center", labelpad=34)
    if mi == 1:
        tag = " (−)" if is_rc else ""
        ax_pwm.text(0.5, -0.22, MOTIF_NM + tag,
                    transform=ax_pwm.transAxes, ha="center", va="top",
                    fontsize=5.5, style="italic", color="#333")

# — Coverage section —
cov_gs = gridspec.GridSpecFromSubplotSpec(
    N_CT + 1, 1, subplot_spec=outer[1],
    height_ratios=H_RATIOS, hspace=0.06)
for ri, ct in enumerate(CELL_TYPES):
    ax = fig.add_subplot(cov_gs[ri])
    draw_cov(ax, bin_centers, cov[ct], cov_n[ct], cov_ymax, ytick=True)
    if ri == 0:
        ax.annotate("ATAC (observed, RPM)", xy=(0.5, 1),
                    xycoords="axes fraction", xytext=(0, 16),
                    textcoords="offset points",
                    fontsize=5, color="#555", ha="center", va="bottom",
                    annotation_clip=False)
    # x-axis on bottom CT row
    if ri == N_CT - 1:
        ax.set_xticks([COV_START, VAR_POS, COV_END])
        ax.set_xticklabels(
            [f"{COV_START//1000}k", f"{VAR_POS}", f"{COV_END//1000}k"],
            fontsize=4.5)
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_linewidth(0.4)
        ax.tick_params(axis="x", which="both", length=2, width=0.4, pad=1,
                       labelsize=4.5)

ax_blank = fig.add_subplot(cov_gs[N_CT])
ax_blank.set_visible(False)

fig.text(0.005, 0.52, "ISM contribution (CT − mean, ref − alt)",
         va="center", rotation="vertical", fontsize=5)
fig.text(0.49, 0.52, "ATAC coverage (RPM)",
         va="center", rotation="vertical", fontsize=5)

out1 = os.path.join(OUT_DIR, "chr16_tracks_ism.pdf")
fig.savefig(out1, dpi=300, bbox_inches="tight")
print(f"Saved: {out1}")
fig.savefig(out1.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# ── 10. Figure: tracks + ISM, per-CT y-limits ─────────────────────────────────
fig_ct = plt.figure(figsize=(10, 4.0))
fig_ct.suptitle("chr16:89622209:A:G  (HNF1A caQTL)  —  per-CT y-scale",
                fontsize=7.5, fontweight="bold", y=1.01)

outer_ct = gridspec.GridSpec(1, 2, figure=fig_ct, wspace=0.22,
                             width_ratios=[2, 3],
                             top=0.93, bottom=0.04, left=0.08, right=0.99)

ism_outer_ct = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_ct[0], wspace=0.06)
for mi, (logo_dict, yabs_v, lbl) in enumerate([
    (cbpnet_diff,   cbp_diff_ylim, "ChromBPNet\nSHAP"),
    (borzoi_diff,   bor_diff_ylim, "Cerberus\nISM"),
]):
    gs = gridspec.GridSpecFromSubplotSpec(
        N_CT + 1, 1, subplot_spec=ism_outer_ct[mi],
        height_ratios=H_RATIOS, hspace=0.06)
    for ri, ct in enumerate(CELL_TYPES):
        ax = fig_ct.add_subplot(gs[ri])
        draw_ism(ax, logo_dict.get(ct), yabs_v, ytick=True)
        if mi == 0:
            ax.set_ylabel(CT_LABELS[ct], fontsize=6, rotation=0,
                          ha="right", va="center", labelpad=34)
        if ri == 0:
            ax.annotate(lbl, xy=(0.5, 1), xycoords="axes fraction",
                        xytext=(0, 16), textcoords="offset points",
                        fontsize=5, color="#555", ha="center", va="bottom",
                        annotation_clip=False)
    ax_pwm = fig_ct.add_subplot(gs[N_CT])
    if pwm_df is not None:
        logomaker.Logo(pwm_df, color_scheme=NT_COLORS, ax=ax_pwm,
                       flip_below=False, baseline_width=0, font_name="DejaVu Sans")
        for p in ax_pwm.patches: p.set_zorder(3)
    ax_pwm.axvspan(24.5, 25.5, color="#d0e8ff", alpha=0.55, zorder=0, lw=0)
    ax_pwm.set_xlim(-0.5, 49.5); ax_pwm.set_xticks([])
    ax_pwm.spines[["top", "right", "bottom"]].set_visible(False)
    ax_pwm.spines["left"].set_linewidth(0.4)
    ax_pwm.tick_params(axis="y", labelsize=5.5, width=0.4, length=1.5, pad=1)
    ax_pwm.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
    if mi == 0:
        ax_pwm.set_ylabel("PWM\n(IC)", fontsize=5, rotation=0,
                          ha="right", va="center", labelpad=34)
    if mi == 1:
        tag = " (−)" if is_rc else ""
        ax_pwm.text(0.5, -0.22, MOTIF_NM + tag,
                    transform=ax_pwm.transAxes, ha="center", va="top",
                    fontsize=5.5, style="italic", color="#333")

cov_gs_ct = gridspec.GridSpecFromSubplotSpec(
    N_CT + 1, 1, subplot_spec=outer_ct[1],
    height_ratios=H_RATIOS, hspace=0.06)
for ri, ct in enumerate(CELL_TYPES):
    ax = fig_ct.add_subplot(cov_gs_ct[ri])
    draw_cov(ax, bin_centers, cov[ct], cov_n[ct], cov_ymax_ct[ct], ytick=True)
    if ri == 0:
        ax.annotate("ATAC (observed, RPM, per-CT scale)", xy=(0.5, 1),
                    xycoords="axes fraction", xytext=(0, 16),
                    textcoords="offset points",
                    fontsize=5, color="#555", ha="center", va="bottom",
                    annotation_clip=False)
    if ri == N_CT - 1:
        ax.set_xticks([COV_START, VAR_POS, COV_END])
        ax.set_xticklabels(
            [f"{COV_START//1000}k", f"{VAR_POS}", f"{COV_END//1000}k"],
            fontsize=4.5)
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_linewidth(0.4)
        ax.tick_params(axis="x", which="both", length=2, width=0.4, pad=1,
                       labelsize=4.5)

ax_blank_ct = fig_ct.add_subplot(cov_gs_ct[N_CT])
ax_blank_ct.set_visible(False)

fig_ct.text(0.005, 0.52, "ISM contribution (CT − mean)",
            va="center", rotation="vertical", fontsize=5)
fig_ct.text(0.49, 0.52, "ATAC coverage (RPM)",
            va="center", rotation="vertical", fontsize=5)

out1ct = os.path.join(OUT_DIR, "chr16_tracks_ism_perct.pdf")
fig_ct.savefig(out1ct, dpi=300, bbox_inches="tight")
print(f"Saved: {out1ct}")
fig_ct.savefig(out1ct.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
plt.close(fig_ct)

# ── 12. Figure: boxplots (peak accessibility per sample per genotype) ──────────
N_CT = len(CELL_TYPES)
fig2, axes = plt.subplots(1, N_CT, figsize=(3.0 * N_CT, 2.8), sharey=True)
fig2.suptitle("chr16:89622209:A:G  —  peak ATAC RPM per sample",
              fontsize=7, fontweight="bold", y=1.02)

for ci, ct in enumerate(CELL_TYPES):
    ax = axes[ci]
    data = [np.array(peak_counts[ct][g]) for g in GENO_ORDER]
    labels_box = ["0/0", "0/1", "1/1"]
    colors_box  = [GENO_COLORS[g] for g in GENO_ORDER]

    bp = ax.boxplot(data, patch_artist=True, widths=0.5,
                    medianprops=dict(color="white", linewidth=1.5),
                    whiskerprops=dict(linewidth=0.7),
                    capprops=dict(linewidth=0.7),
                    flierprops=dict(marker="o", markersize=2, alpha=0.5))
    for patch, col in zip(bp["boxes"], colors_box):
        patch.set(facecolor=col, alpha=0.7)

    for gi, (gdata, col) in enumerate(zip(data, colors_box)):
        jitter = np.random.default_rng(gi).uniform(-0.15, 0.15, len(gdata))
        ax.scatter(np.full(len(gdata), gi + 1) + jitter, gdata,
                   color=col, s=8, alpha=0.6, zorder=3, linewidths=0)

    ax.set_title(CT_LABELS[ct], fontsize=7)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(labels_box, fontsize=6)
    ax.set_xlabel("genotype", fontsize=6)
    ax.tick_params(axis="y", labelsize=5.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.5)
    if ci == 0:
        ax.set_ylabel("peak RPM", fontsize=6)

fig2.tight_layout()
out2 = os.path.join(OUT_DIR, "chr16_boxplot.pdf")
fig2.savefig(out2, dpi=300, bbox_inches="tight")
print(f"Saved: {out2}")
fig2.savefig(out2.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
plt.close(fig2)

# ── 13. Combined figure: bigwig coverage | boxplot | ChromBPNet SHAP (CT-mean) ─
print("Building combined figure …")

# DPEP1 promoter — MANE Select transcript ENST00000690203.1 (gencode v48, + strand).
# The caQTL peak is an intronic enhancer ~8.6 kb downstream of the DPEP1 TSS. Show a
# narrow window spanning the DPEP1 promoter (left) and the enhancer peak (right).
DPEP1_STRAND = "+"
DPEP1_TSS    = 89613642            # MANE TSS = nearest promoter to the peak
DPEP1_EXON1  = (89613642, 89613719)

PEAK_MID  = (PEAK_S + PEAK_E) // 2      # 89,622,209 (= variant / peak summit)
GENE_DIST = PEAK_MID - DPEP1_TSS        # promoter → peak distance (bp)

COV2_START = 89612500
COV2_END   = 89623500
COV2_BINS  = (COV2_END - COV2_START) // BIN_SIZE
cov2_x     = COV2_START + (np.arange(COV2_BINS) + 0.5) * BIN_SIZE

# Compute genotype-stratified bigwig coverage over 2kb window
print("  Computing bigwig coverage …")
cov_bw   = {}
cov_bw_n = {}
for ct in CELL_TYPES:
    cov_bw[ct] = {}
    cov_bw_n[ct] = {}
    for geno, samples in geno_groups.items():
        arr, n = group_coverage_bw(ct, samples, CHR, COV2_START, COV2_END, COV2_BINS)
        cov_bw[ct][geno]   = arr
        cov_bw_n[ct][geno] = n

cov_bw_ymax = max(
    (cov_bw[ct][g].max() for ct in CELL_TYPES for g in GENO_ORDER
     if cov_bw[ct][g] is not None),
    default=1.0
) * 1.1

# shared boxplot y-limit across all CTs
box_ymax = max(
    v for ct in CELL_TYPES for g in GENO_ORDER for v in peak_counts[ct][g]
) * 1.15

def draw_boxplot_row(ax, data_dict, ymax, bottom_row=False):
    data = [np.array(data_dict[g]) for g in GENO_ORDER]
    bp = ax.boxplot(
        data, patch_artist=True, widths=0.55,
        medianprops=dict(color="white", linewidth=1.2),
        whiskerprops=dict(linewidth=0.6), capprops=dict(linewidth=0.6),
        flierprops=dict(marker="o", markersize=1.5, alpha=0, linestyle="none"),
    )
    for patch, g in zip(bp["boxes"], GENO_ORDER):
        patch.set(facecolor=GENO_COLORS[g], alpha=0.7)
    for med_line, gkey in zip(bp["medians"], GENO_ORDER):
        h = GENO_COLORS[gkey].lstrip("#")
        rv, gv, bv = int(h[:2],16)/255, int(h[2:4],16)/255, int(h[4:],16)/255
        med_line.set(color=(rv*0.55, gv*0.55, bv*0.55), linewidth=1.5)
    for gi, g in enumerate(GENO_ORDER):
        gd = np.array(data_dict[g])
        jitter = np.random.default_rng(gi).uniform(-0.1, 0.1, len(gd))
        ax.scatter(np.full(len(gd), gi + 1) + jitter, gd,
                   color="#222", s=0.8, alpha=0.6, zorder=3, linewidths=0)
    ax.set_xlim(0.3, 3.7)
    ax.set_ylim(-ymax * 0.06, ymax)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.4)
    ax.tick_params(axis="y", labelsize=6, width=0.4, length=1.5, pad=1)
    ax.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
    if bottom_row:
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(["0/0", "0/1", "1/1"], fontsize=6,
                           rotation=40, ha="right", rotation_mode="anchor")
        ax.tick_params(axis="x", width=0.4, length=1.5, pad=1)
    else:
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)

H_RATIOS = [1.0] * N_CT + [0.5]
fig3 = plt.figure(figsize=(4.46, 3.60))   # 30% larger, fonts stay 6pt

outer3 = gridspec.GridSpec(
    1, 3, figure=fig3,
    width_ratios=[1.1, 0.36, 1.05],
    wspace=0.26, top=0.83, bottom=0.13, left=0.16, right=0.99,
)

# ── Coverage (left) — fragment RPM, genotype-stratified, 2kb window ──────────
cov_gs3 = gridspec.GridSpecFromSubplotSpec(
    N_CT + 1, 1, subplot_spec=outer3[0],
    height_ratios=H_RATIOS, hspace=0.06)

ISM_WIN    = 512   # ±512 bp = 1024 bp ChromBPNet/Cerberus local scoring window
ISM_WIN_S  = VAR_POS - ISM_WIN
ISM_WIN_E  = VAR_POS + ISM_WIN

for ri, ct in enumerate(CELL_TYPES):
    ax = fig3.add_subplot(cov_gs3[ri])
    # model scoring window (lightest, drawn first)
    ax.axvspan(ISM_WIN_S, ISM_WIN_E, color="#cce5ff", alpha=0.35, zorder=0, lw=0)
    # peak region
    ax.axvspan(PEAK_S, PEAK_E, color="#d4edda", alpha=0.55, zorder=1, lw=0)
    for geno in GENO_ORDER:
        y = cov_bw[ct].get(geno)
        if y is None: continue
        c = GENO_COLORS[geno]
        ax.fill_between(cov2_x, 0, y, color=c, alpha=0.45, linewidth=0, zorder=2)
        ax.plot(cov2_x, y, color=c, linewidth=0.5, zorder=3)
    ax.axvline(VAR_POS, color="#555", lw=0.6, ls="--", zorder=4, alpha=0.8)
    ax.set_xlim(COV2_START, COV2_END)
    ax.set_ylim(0, cov_bw_ymax)
    _ax_style(ax, ytick=True, labelsize=6)
    GT_LABEL = {"ref": "0/0", "het": "0/1", "alt": "1/1"}
    ax.set_ylabel(CT_LABELS[ct], fontsize=6, rotation=0,
                  ha="right", va="center", labelpad=3)
    if ri == 0:
        ax.set_title("ATAC coverage (RPM)", fontsize=6, color="#555", pad=3)
        # variant marker (title removed)
        fig3.text(0.5, 1.04, "chr16:89,622,209  A>G", ha="center", va="top",
                  fontsize=6, fontweight="bold", color="#333")
        # single genotype legend + peak/SHAP-window legend
        geno_h = [mpatches.Patch(facecolor=GENO_COLORS[g], alpha=0.7,
                  label=f"{GT_LABEL[g]} (n={cov_bw_n['PT'].get(g,0)})") for g in GENO_ORDER]
        win_h = [mpatches.Patch(facecolor="#d4edda", label="peak"),
                 mpatches.Patch(facecolor="#cce5ff", label="SHAP window")]
        fig3.legend(handles=geno_h, ncol=3, fontsize=6, frameon=False,
                    loc="upper left", bbox_to_anchor=(0.15, 0.99),
                    columnspacing=1.0, handlelength=1.0, handletextpad=0.3)
        fig3.legend(handles=win_h, ncol=1, fontsize=6, frameon=False,
                    loc="upper right", bbox_to_anchor=(0.995, 1.0),
                    labelspacing=0.3, handlelength=1.0, handletextpad=0.3)
    # genomic x-axis is drawn on the gene-annotation row below, not on the CT rows

# ── DPEP1 promoter + promoter→peak distance (bottom row of coverage column) ───
ax_gene = fig3.add_subplot(cov_gs3[N_CT])
ax_gene.set_xlim(COV2_START, COV2_END)
ax_gene.set_ylim(0, 1)
# peak / scoring-window / variant markers, for continuity with the coverage rows
ax_gene.axvspan(ISM_WIN_S, ISM_WIN_E, color="#cce5ff", alpha=0.35, zorder=0, lw=0)
ax_gene.axvspan(PEAK_S, PEAK_E, color="#d4edda", alpha=0.55, zorder=1, lw=0)
ax_gene.axvline(VAR_POS, color="#555", lw=0.6, ls="--", zorder=4, alpha=0.8)

WSPAN = COV2_END - COV2_START
# DPEP1 promoter glyph (bent arrow, + strand → pointing right)
PY_BASE = 0.55
ax_gene.plot([DPEP1_TSS, DPEP1_TSS], [PY_BASE, PY_BASE + 0.34],
             color="#333", lw=0.9, zorder=3)
ax_gene.annotate("", xy=(DPEP1_TSS + WSPAN * 0.06, PY_BASE + 0.34),
                 xytext=(DPEP1_TSS, PY_BASE + 0.34),
                 arrowprops=dict(arrowstyle="-|>", color="#333", lw=0.9,
                                 mutation_scale=6), zorder=3)
# exon 1 block
ax_gene.add_patch(mpatches.Rectangle(
    (DPEP1_EXON1[0], PY_BASE - 0.12), max(DPEP1_EXON1[1] - DPEP1_EXON1[0], WSPAN * 0.006),
    0.24, facecolor="#333", edgecolor="none", zorder=3))
ax_gene.text(DPEP1_TSS + WSPAN * 0.075, PY_BASE + 0.34, "DPEP1",
             fontsize=6, style="italic", ha="left", va="center", color="#333")

# promoter → peak distance annotation
DY = 0.24
ax_gene.annotate("", xy=(PEAK_MID, DY), xytext=(DPEP1_TSS, DY),
                 arrowprops=dict(arrowstyle="<|-|>", color="#777", lw=0.6,
                                 mutation_scale=5), zorder=3)
ax_gene.text((DPEP1_TSS + PEAK_MID) / 2, DY + 0.06,
             f"{GENE_DIST/1000:.1f} kb", fontsize=5.5, ha="center", va="bottom",
             color="#555")

# genomic coordinate axis
ax_gene.set_yticks([])
ax_gene.set_xticks([COV2_START, VAR_POS, COV2_END])
ax_gene.set_xticklabels(
    [f"{COV2_START:,}", f"{VAR_POS:,}", f"{COV2_END:,}"],
    fontsize=6, rotation=40, ha="right", rotation_mode="anchor")
ax_gene.spines[["top", "right", "left"]].set_visible(False)
ax_gene.spines["bottom"].set_linewidth(0.4)
ax_gene.tick_params(axis="x", length=2, width=0.4, pad=1, labelsize=6)

# ── Boxplots (middle) — shared y-limit ───────────────────────────────────────
box_gs3 = gridspec.GridSpecFromSubplotSpec(
    N_CT + 1, 1, subplot_spec=outer3[1],
    height_ratios=H_RATIOS, hspace=0.06)

for ri, ct in enumerate(CELL_TYPES):
    ax = fig3.add_subplot(box_gs3[ri])
    draw_boxplot_row(ax, peak_counts[ct], box_ymax, bottom_row=(ri == N_CT - 1))
    if ri == 0:
        ax.set_title("peak RPM", fontsize=6, color="#555", pad=3)

ax_blank_box = fig3.add_subplot(box_gs3[N_CT])
ax_blank_box.set_visible(False)

# ── ChromBPNet SHAP (CT − mean), right ───────────────────────────────────────
shap_gs3 = gridspec.GridSpecFromSubplotSpec(
    N_CT + 1, 1, subplot_spec=outer3[2],
    height_ratios=H_RATIOS, hspace=0.06)

for ri, ct in enumerate(CELL_TYPES):
    ax = fig3.add_subplot(shap_gs3[ri])
    draw_ism(ax, cbpnet_diff.get(ct), cbp_diff_ylim, ytick=True, labelsize=6, na_fs=6,
             var_line=True)
    if ri == 0:
        ax.set_title("ChromBPNet SHAP\n(CT − mean, ref − alt)", fontsize=6, color="#555", pad=3)

# PWM row — genomic range marked on x-axis (same ±25 bp window as the SHAP logos)
ax_pwm3 = fig3.add_subplot(shap_gs3[N_CT])
if pwm_df is not None:
    logomaker.Logo(pwm_df, color_scheme=NT_COLORS, ax=ax_pwm3,
                   flip_below=False, baseline_width=0, font_name="DejaVu Sans")
    for p in ax_pwm3.patches: p.set_zorder(3)
ax_pwm3.axvline(25, color="#555", lw=0.6, ls="--", zorder=4, alpha=0.8)
ax_pwm3.set_xlim(-0.5, 49.5)
ax_pwm3.set_xticks([0, 25, 49])
ax_pwm3.set_xticklabels([f"{VAR_POS-25:,}", f"{VAR_POS:,}", f"{VAR_POS+24:,}"],
                        fontsize=6, rotation=40, ha="right", rotation_mode="anchor")
ax_pwm3.spines[["top", "right"]].set_visible(False)
ax_pwm3.spines["left"].set_linewidth(0.4); ax_pwm3.spines["bottom"].set_linewidth(0.4)
ax_pwm3.tick_params(axis="y", labelsize=6, width=0.4, length=1.5, pad=1)
ax_pwm3.tick_params(axis="x", length=2, width=0.4, pad=1, labelsize=6)
ax_pwm3.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
ax_pwm3.set_ylabel("PWM (IC)", fontsize=6, rotation=90,
                   ha="center", va="bottom", labelpad=2)
tag = " (−)" if is_rc else ""
ax_pwm3.text(1.0, 1.0, MOTIF_NM + tag, transform=ax_pwm3.transAxes,
             ha="right", va="bottom", fontsize=6, style="italic", color="#333")

out3 = os.path.join(OUT_DIR, "chr16_combined_figure.pdf")
fig3.savefig(out3, dpi=300, bbox_inches="tight")
print(f"Saved: {out3}")
fig3.savefig(out3.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
plt.close(fig3)

# ── Merged (all-sample) coverage for all 12 CTs, raw + RPM ──────────
print("Building merged coverage figure …")

MERGED_FRAG_DIR = cfg("FRAG_DIR")
ALL_CTS = ["CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
           "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT"]
ALL_CT_LABELS = {
    "CNT_CD_PC": "CNT/CD/PC", "DCT": "DCT", "DTL_ATL": "DTL/ATL",
    "EC": "EC", "IC": "IC", "Immune": "Immune", "PEC": "PEC",
    "PTS": "PTS", "Podocyte": "Podocyte", "Stromal": "Stromal",
    "TAL": "TAL", "injPT": "injPT",
}
MERGED_TOTALS = {
    "CNT_CD_PC": 111751421, "DCT": 49387919, "DTL_ATL": 76407656,
    "EC": 86728422, "IC": 49805670, "Immune": 31253488,
    "PEC": 10263317, "PTS": 269942307, "Podocyte": 10946558,
    "Stromal": 55137511, "TAL": 327268693, "injPT": 50164952,
}
CT_COLORS_12 = {
    "CNT_CD_PC": "#b5bd61", "DCT": "#006fa6", "DTL_ATL": "#ffb3c6",
    "EC": "#d62728", "IC": "#aa40fc", "Immune": "#279e68",
    "PEC": "#6a3a4c", "PTS": "#e377c2", "Podocyte": "#ff7f0e",
    "injPT": "#950046", "Stromal": "#8c6d31", "TAL": "#aec7e8",
}

def merged_coverage(ct, chrom, start, end, bin_size):
    fp = os.path.join(MERGED_FRAG_DIR, f"{ct}_merged_fragments.tsv.gz")
    n_bins = (end - start) // bin_size
    cov = np.zeros(n_bins)
    try:
        tbx = pysam.TabixFile(fp)
        for row in tbx.fetch(chrom, start, end):
            cols = row.split("\t")
            fs, fe = int(cols[1]), int(cols[2])
            cnt = int(cols[4]) if len(cols) > 4 else 1
            b0 = max(0, (fs - start) // bin_size)
            b1 = min(n_bins, (fe - start + bin_size - 1) // bin_size)
            cov[b0:b1] += cnt
        tbx.close()
    except (ValueError, KeyError):
        pass
    return cov

M4_START = COV2_START
M4_END   = COV2_END
M4_BINS  = (M4_END - M4_START) // BIN_SIZE
m4_x     = M4_START + (np.arange(M4_BINS) + 0.5) * BIN_SIZE

raw_cov = {}
rpm_cov = {}
for ct in ALL_CTS:
    raw = merged_coverage(ct, CHR, M4_START, M4_END, BIN_SIZE)
    raw_cov[ct] = raw
    rpm_cov[ct] = raw / MERGED_TOTALS[ct] * 1_000_000

# Two panels side by side: raw | RPM
fig4, axes4 = plt.subplots(
    len(ALL_CTS), 2,
    figsize=(5.5, len(ALL_CTS) * 0.55 + 0.4),
    sharex="col",
)
fig4.subplots_adjust(hspace=0.06, wspace=0.35)

raw_ymax = max(raw_cov[ct].max() for ct in ALL_CTS if raw_cov[ct].max() > 0)
rpm_ymax = max(rpm_cov[ct].max() for ct in ALL_CTS if rpm_cov[ct].max() > 0)

for ri, ct in enumerate(ALL_CTS):
    color = CT_COLORS_12[ct]
    for ci, (cov_arr, ymax, col_lbl) in enumerate([
        (raw_cov[ct], raw_ymax, "Raw counts"),
        (rpm_cov[ct], rpm_ymax, "RPM"),
    ]):
        ax = axes4[ri, ci]
        ax.fill_between(m4_x, cov_arr, color=color, lw=0, alpha=0.8)
        ax.plot(m4_x, cov_arr, color=color, lw=0.4)
        ax.axvspan(max(PEAK_S, M4_START), min(PEAK_E, M4_END),
                   color="#ffe8a0", alpha=0.35, zorder=0, lw=0)
        ax.axvline(VAR_POS, color="#CC3333", lw=0.5, ls="--", zorder=3)
        ax.set_xlim(M4_START, M4_END)
        ax.set_ylim(0, ymax * 1.08)
        ax.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
        ax.tick_params(axis="y", labelsize=4.5, width=0.35, length=1.5, pad=1)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_linewidth(0.35)
        if ri < len(ALL_CTS) - 1:
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
        else:
            ax.tick_params(axis="x", labelsize=4.5, width=0.35, length=1.5, pad=1)
            ax.set_xticks([M4_START, VAR_POS, M4_END])
            ax.set_xticklabels(
                [f"{M4_START:,}", str(VAR_POS), f"{M4_END:,}"],
                fontsize=4, rotation=15, ha="right")
        if ci == 0:
            ax.set_ylabel(ALL_CT_LABELS[ct], fontsize=5, rotation=0,
                          ha="right", va="center", labelpad=30)

axes4[0, 0].set_title("Raw counts", fontsize=5.5, pad=3)
axes4[0, 1].set_title("RPM", fontsize=5.5, pad=3)
fig4.text(0.5, 1.01, f"chr16:{M4_START:,}–{M4_END:,}  (variant {VAR_POS:,} A→G)",
          ha="center", va="bottom", fontsize=5.5, transform=fig4.transFigure)

out4 = os.path.join(OUT_DIR, "chr16_merged_coverage.pdf")
fig4.savefig(out4, dpi=300, bbox_inches="tight")
print(f"Saved: {out4}")
plt.close(fig4)

print("Done.")
