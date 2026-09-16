#!/usr/bin/env python
"""eGFR GWAS locus: credible set, model scores, and the disrupted motif.

The FGF5 locus, where the eGFR GWAS signal colocalizes with both a proximal
tubule caQTL and an eQTL. Rows, top to bottom:

    eGFR association over +/-200 kb, credible-set variants highlighted
    SuSiE fine-mapping PIP
    predicted effect of each credible-set variant, all three statistics
    ChromBPNet and Cerberus predicted accessibility on both alleles
    their attributions, and the aligned motif
    Cerberus predicted RNA across the gene body, gene model, RNA attributions

The motif is not specified here: it is chosen by matching the JASPAR database
against the alternate-allele attribution profile, so the panel shows what the
models actually implicate.

Usage:
    python fgf5_gwas_locus.py [--out_dir <dir>]
"""
import argparse
import os, re
import numpy as np
import pandas as pd
import h5py
import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

from figlib import (
    ALT_COLOR, CT_COLORS, NT_COLORS, REF_COLOR,
    align_pwm_to_cwm, cfg, parse_gtf_gene, parse_pwm, setup_style,
)

setup_style()

# ── Paths ──────────────────────────────────────────────────────────────────────
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--out_dir", default=None)
_parser.add_argument("--locus", default="rs16998073_chr4_80684341_81684341",
                     help="SuSiE locus directory name under GWAS_SUSIE_DIR")
_args = _parser.parse_args()

WORK_DIR = cfg("WORK_DIR")
OUT_DIR  = _args.out_dir or cfg("FIGURE_DIR")
GWAS_DIR = os.path.join(WORK_DIR, "gwas")

# Fine-mapping output for this locus (controlled access; see the README).
LOCUS_DIR  = os.path.join(cfg("GWAS_SUSIE_DIR"), _args.locus)
GTF_FILE   = cfg("EQTL_GTF")
ATAC_NPZ   = os.path.join(GWAS_DIR, "cerberus_pred", "chr4_80261400_T_C.npz")
RNA_NPZ    = os.path.join(GWAS_DIR, "cerberus_pred", "chr4_80261400_T_C_PTS.npz")
ATAC_ISM   = os.path.join(GWAS_DIR, "cerberus_ism", "atac", "scores.h5")
RNA_ISM    = os.path.join(GWAS_DIR, "cerberus_ism", "rna", "scores.h5")
CBP_PRED   = os.path.join(GWAS_DIR, "chrombpnet_pred", "chr4_80261400_T_C.npz")
CBP_SHAP   = os.path.join(GWAS_DIR, "chrombpnet_shap", "PT_ism.h5")
CBP_VKEY   = "chr4_80261400_T_C"
CBP_COV_FLANK = 500   # ±bp shown for ChromBPNet coverage (panel G COV_WIN_BP)
CS_BORZOI  = os.path.join(GWAS_DIR, "cs_cerberus_scores.tsv")
CS_CBP     = os.path.join(GWAS_DIR, "cs_chrombpnet_scores.tsv")
OUT_PDF    = os.path.join(OUT_DIR, "fgf5_gwas_locus.pdf")
OUT_PNG    = os.path.join(OUT_DIR, "fgf5_gwas_locus.png")

# ── Constants ─────────────────────────────────────────────────────────────────
VAR_CHROM     = "chr4"
VAR_POS       = 80261400     # rs12509595, T->C
VAR_REF       = "T"
VAR_ALT       = "C"
WIDE_FLANK    = 200_000
ATAC_PRED_FLANK = 1_000
RNA_PRED_START  = 80255000
RNA_PRED_END    = 80345000
PEAK_START    = 80261167     # PTS_peak_125977 (MACS2)
PEAK_END      = 80261480

ATAC_ISM_TIDX = 10   # PTS ATAC in targets_panel_h.txt (11 tracks)
RNA_ISM_TIDX  = 15   # PTS RNA+ in rna_ism/targets_cov.txt (124 tracks)
RNA_ISM_GENE  = b"ENSG00000138675.17"   # FGF5
ATAC_VAR_IDX  = 24   # variant index in 50-bp ATAC ISM window
RNA_VAR_IDX   = 49   # variant index in 100-bp RNA ISM window

MOTIF_DB   = cfg("MOTIF_DB")
MOTIF_ID   = "MA0067.2"   # PAX2 (represents PAX2/8 paralogs) — kidney master TF
MOTIF_NAME = "PAX2/8"

CS_COLOR   = "#E07B39"
GRAY_DASH  = dict(color="#888888", lw=0.6, ls="--")
HIGHLIGHT  = "#ffffb3"
PEAK_SHADE = "#d9d9d9"

# ── ISM logo helpers (match panels G and I) ───────────────────────────────────
def atac_ism_logo(seqs_oh, logstat, tidx):
    """Panel G style: contribution = actual - mean(others)."""
    L = seqs_oh.shape[1]; oh_T = seqs_oh.T
    ls = logstat[:, :, tidx].astype(float); arr = np.zeros((L, 4))
    for p in range(L):
        idx = int(np.argmax(oh_T[p])); others = [j for j in range(4) if j != idx]
        arr[p, idx] = ls[p, idx] - np.mean(ls[p, others])
    return pd.DataFrame(arr, columns=list("ACGT"))

def rna_ism_logo(seqs_oh, logstat, tidx):
    """Panel I style: contribution = -mean(others) at actual base."""
    L = seqs_oh.shape[1]; oh_T = seqs_oh.T
    ls = logstat[:, :, tidx].astype(float); arr = np.zeros((L, 4))
    for p in range(L):
        idx = int(np.argmax(oh_T[p])); others = [j for j in range(4) if j != idx]
        arr[p, idx] = -float(np.mean(ls[p, others]))
    return pd.DataFrame(arr, columns=list("ACGT"))

def draw_logo(ax, logo_df, ylim, var_pos_win, n_bases):
    logomaker.Logo(logo_df, color_scheme=NT_COLORS, ax=ax,
                   flip_below=True, baseline_width=0, font_name="DejaVu Sans")
    for p in ax.patches: p.set_zorder(3)
    ax.axvline(var_pos_win, color="#888", lw=0.5, ls="--", zorder=1)
    ax.axhline(0, color="#bbbbbb", lw=0.4, zorder=1)
    ax.set_xlim(-0.5, n_bases - 0.5)
    ax.set_ylim(*ylim)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=5, length=2)
    ax.spines[["top", "right", "bottom"]].set_visible(False)

# ── PWM helpers (match panel I) ───────────────────────────────────────────────
def draw_pwm(ax, pwm_df, var_pos_win, n_bases, rc, label):
    logomaker.Logo(pwm_df, color_scheme=NT_COLORS, ax=ax,
                   baseline_width=0, font_name="DejaVu Sans")
    ax.axvline(var_pos_win, color="#888", lw=0.5, ls="--", zorder=1)
    ax.set_xlim(-0.5, n_bases - 0.5)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    tag = " (−)" if rc else ""
    ax.text(0.5, -0.12, f"{label}{tag}", transform=ax.transAxes, ha="center",
            va="top", fontsize=5, style="italic", color="#333")

# ── Load GWAS data ────────────────────────────────────────────────────────────
print("Loading GWAS data...")
ss = pd.read_csv(os.path.join(LOCUS_DIR, "sumstats.tsv"), sep="\t")
pip_df = pd.read_csv(os.path.join(LOCUS_DIR, "pip.tsv"), sep="\t")
wide_lo, wide_hi = VAR_POS - WIDE_FLANK, VAR_POS + WIDE_FLANK

ss_win = ss[(ss["hg38_chr"] == 4) & ss["hg38_bp"].between(wide_lo, wide_hi)].copy()
ss_win["neglog10p"] = -np.log10(ss_win["P"].clip(lower=1e-300))
pip_win = pip_df[(pip_df["hg38_chr"] == 4) & pip_df["hg38_bp"].between(wide_lo, wide_hi)].copy()
cs_pos = set(pip_win.loc[pip_win["cs"] != -1, "hg38_bp"])
print(f"  Variants: {len(ss_win)}, credible set: {len(cs_pos)}")

# ── Credible-set model effect scores (ChromBPNet / Cerberus ATAC / Cerberus RNA) ──
csb = pd.read_csv(CS_BORZOI, sep="\t")
csc = pd.read_csv(CS_CBP, sep="\t")
cs_scores = csb.merge(csc[["SNP", "chrombpnet_log2fc"]], on="SNP").sort_values("pos")
CS_LO, CS_HI = cs_scores["pos"].min() - 1500, cs_scores["pos"].max() + 1500

# ── Load predictions ──────────────────────────────────────────────────────────
print("Loading predictions...")
atac_npz = np.load(ATAC_NPZ); rna_npz = np.load(RNA_NPZ)
bins_start = int(atac_npz["bins_start"]); bin_size = int(atac_npz["bin_size"])
n_bins = len(atac_npz["PT_ref"])
bin_centers = bins_start + (np.arange(n_bins) + 0.5) * bin_size
atac_ref, atac_alt = atac_npz["PT_ref"], atac_npz["PT_alt"]
rna_ref,  rna_alt  = rna_npz["pred_ref"], rna_npz["pred_alt"]
atac_lo, atac_hi = VAR_POS - ATAC_PRED_FLANK, VAR_POS + ATAC_PRED_FLANK
atac_mask = (bin_centers >= atac_lo) & (bin_centers <= atac_hi)
rna_mask  = (bin_centers >= RNA_PRED_START) & (bin_centers <= RNA_PRED_END)

# ── Load ISM ─────────────────────────────────────────────────────────────────
print("Loading ISM...")
with h5py.File(ATAC_ISM) as h5:
    atac_ref_oh   = h5["ref/seqs"][0]   # (4, 50)
    atac_logo_ref = atac_ism_logo(atac_ref_oh, h5["ref/cov/logSUM"][0], ATAC_ISM_TIDX)
    atac_logo_alt = atac_ism_logo(h5["alt/seqs"][0], h5["alt/cov/logSUM"][0], ATAC_ISM_TIDX)
a_max = max(atac_logo_ref.values.max(), atac_logo_alt.values.max())
a_min = min(atac_logo_ref.values.min(), atac_logo_alt.values.min())
atac_ylim = (a_min - abs(a_max - a_min) * 0.05, a_max + abs(a_max - a_min) * 0.05)

with h5py.File(RNA_ISM) as h5:
    gene_ids = h5["gene_ids"][:]
    pidx = int(np.where(gene_ids == RNA_ISM_GENE)[0][0])
    snp_i = int(h5["snp_idx"][pidx])
    rna_ref_oh   = h5["ref/seqs"][snp_i]   # (4, 100)
    rna_logo_ref = rna_ism_logo(rna_ref_oh, h5["ref/covgene/logSED"][pidx], RNA_ISM_TIDX)
    rna_logo_alt = rna_ism_logo(h5["alt/seqs"][snp_i], h5["alt/covgene/logSED"][pidx], RNA_ISM_TIDX)
r_max = max(rna_logo_ref.values.max(), rna_logo_alt.values.max())
r_min = min(rna_logo_ref.values.min(), rna_logo_alt.values.min())
rna_ylim = (r_min - abs(r_max - r_min) * 0.05, r_max + abs(r_max - r_min) * 0.05)

# ── Load ChromBPNet predictions + SHAP contributions (PT) ─────────────────────
print("Loading ChromBPNet...")
cbp_npz   = np.load(CBP_PRED)
cbp_ostart = int(cbp_npz["out_start"][0]); cbp_olen = int(cbp_npz["out_len"][0])
cbp_x      = cbp_ostart + np.arange(cbp_olen) + 0.5
cbp_ref, cbp_alt = cbp_npz["PT_ref"], cbp_npz["PT_alt"]
cbp_lo, cbp_hi = VAR_POS - CBP_COV_FLANK, VAR_POS + CBP_COV_FLANK

with h5py.File(CBP_SHAP) as h5:
    ids = list(h5["variant_ids"][:].astype(str))
    si  = ids.index(CBP_VKEY)
    shap_ref_oh = h5["ref_onehot"][si]   # (50, 4)
    shap_alt_oh = h5["alt_onehot"][si]
    shap_ref_ct = h5["ref_contrib"][si]  # (50, 4)
    shap_alt_ct = h5["alt_contrib"][si]
shap_logo_ref = pd.DataFrame(shap_ref_oh * shap_ref_ct, columns=list("ACGT"))
shap_logo_alt = pd.DataFrame(shap_alt_oh * shap_alt_ct, columns=list("ACGT"))
SHAP_VAR_IDX  = int(np.where((shap_ref_oh != shap_alt_oh).any(axis=1))[0][0])
s_max = max(shap_logo_ref.values.max(), shap_logo_alt.values.max())
s_min = min(shap_logo_ref.values.min(), shap_logo_alt.values.min())
shap_ylim = (s_min - abs(s_max - s_min) * 0.05, s_max + abs(s_max - s_min) * 0.05)

# Align NR2C2 PWM by matching its IC profile to each panel's ALT contribution
# scores (CWM / alt ISM), the TomTom convention — not genomic log-odds.
pwm_prob = parse_pwm(MOTIF_ID)
pwm_prob = pwm_prob / pwm_prob.sum(axis=1, keepdims=True)
atac_pwm_df, atac_rc, _ = align_pwm_to_cwm(atac_logo_alt.values, pwm_prob, 50)
rna_pwm_df,  rna_rc, _  = align_pwm_to_cwm(rna_logo_alt.values,  pwm_prob, 100)
shap_pwm_df, shap_rc, _ = align_pwm_to_cwm(shap_logo_alt.values, pwm_prob, 50)
print(f"  {MOTIF_NAME} CWM-aligned: ATAC rc={atac_rc}, RNA rc={rna_rc}, SHAP rc={shap_rc}")

print("Parsing GTF...")
fgf5_start, fgf5_end, fgf5_exons, _ = parse_gtf_gene(GTF_FILE, "FGF5")
cfap_start, cfap_end, cfap_exons, _ = parse_gtf_gene(GTF_FILE, "CFAP299")

# exonic coverage sums over FGF5 (panel I v2 Σ legend)
def exon_sums(ref, alt, xf, exons):
    r = s = 0.0
    for es, ee in exons:
        m = (xf >= es) & (xf <= ee)
        r += ref[m].sum(); s += alt[m].sum()
    return r, s
rna_rsum, rna_asum = exon_sums(rna_ref, rna_alt, bin_centers, fgf5_exons)

# ── Manual figure layout ──────────────────────────────────────────────────────
# Manhattan + ATAC coverage rows are 40% shorter than before.
# ATAC section is two columns: ChromBPNet (left) and Cerberus (right).
FIG_W      = 4.2
MARGIN_TOP = 0.20   # inches (room for column titles)
MARGIN_BOT = 0.34   # inches (for x-label)
FULL_L, FULL_W   = 0.15, 0.82
COL_L_L, COL_W   = 0.155, 0.355   # ChromBPNet (left) column
COL_L_R          = 0.585          # Cerberus (right) column
# three side-by-side columns for the credible-set prioritization row
CS3_L  = [0.155, 0.445, 0.735]
CS3_W  = 0.235

# (key, height_in, mode, gap_below_in)   mode: "full" | "atac2" | "cs3"
ROWS = [
    ("manh",     0.41, "full",  0.05),
    ("pip",      0.45, "full",  0.34),
    ("csprior",  0.31, "cs3",   0.30),
    ("atac_cov", 0.57, "atac2", 0.05),
    ("atac_ref", 0.42, "atac2", 0.03),
    ("atac_alt", 0.42, "atac2", 0.03),
    ("atac_pwm", 0.18, "atac2", 0.34),
    ("rna_pred", 0.80, "full",  0.02),
    ("gene",     0.28, "full",  0.30),
    ("rna_ism_ref", 0.42, "full", 0.03),
    ("rna_ism_alt", 0.42, "full", 0.03),
    ("rna_pwm",     0.18, "full", 0.00),
]
FIG_H = MARGIN_TOP + MARGIN_BOT + sum(h for _, h, _, _ in ROWS) + sum(g for *_, g in ROWS)

fig = plt.figure(figsize=(FIG_W, FIG_H))
A = {}
y_cursor = FIG_H - MARGIN_TOP
for key, h, mode, gap in ROWS:
    bottom = y_cursor - h
    if mode == "full":
        A[key] = fig.add_axes([FULL_L, bottom / FIG_H, FULL_W, h / FIG_H])
    elif mode == "cs3":  # three side-by-side columns
        for j in range(3):
            A[f"{key}_{j}"] = fig.add_axes([CS3_L[j], bottom / FIG_H, CS3_W, h / FIG_H])
    else:  # atac2: left (ChromBPNet) + right (Cerberus)
        A[key + "_l"] = fig.add_axes([COL_L_L, bottom / FIG_H, COL_W, h / FIG_H])
        A[key + "_r"] = fig.add_axes([COL_L_R, bottom / FIG_H, COL_W, h / FIG_H])
    y_cursor = bottom - gap

# Column titles over the ATAC section
_cov_top = (A["atac_cov_l"].get_position().y1)
fig.text(COL_L_L + COL_W / 2, _cov_top + 0.004, "ChromBPNet",
         ha="center", va="bottom", fontsize=6.5, color="#333333")
fig.text(COL_L_R + COL_W / 2, _cov_top + 0.004, "Cerberus",
         ha="center", va="bottom", fontsize=6.5, color="#333333")

# ── Manhattan ─────────────────────────────────────────────────────────────────
ax = A["manh"]
non_cs   = ss_win[~ss_win["hg38_bp"].isin(cs_pos)]
focal    = ss_win[ss_win["hg38_bp"] == VAR_POS]
other_cs = ss_win[ss_win["hg38_bp"].isin(cs_pos) & (ss_win["hg38_bp"] != VAR_POS)]
ax.scatter(non_cs["hg38_bp"], non_cs["neglog10p"], s=3, color="#9a9a9a",
           alpha=0.6, rasterized=True, linewidths=0)
ax.scatter(other_cs["hg38_bp"], other_cs["neglog10p"], s=9, color=CS_COLOR,
           alpha=0.95, zorder=4, edgecolors="white", linewidths=0.3)
# focal: same size as CS dots, black outline, drawn on top
ax.scatter(focal["hg38_bp"], focal["neglog10p"], s=9, color=CS_COLOR,
           edgecolors="black", linewidths=0.7, zorder=6)
ax.axhline(-np.log10(5e-8), color="#d62728", lw=0.7, ls="--")
ax.axvline(VAR_POS, **GRAY_DASH)
ax.set_ylabel("−log$_{10}$(P)", fontsize=6.5)
ax.set_ylim(0, ss_win["neglog10p"].max() * 1.08)
ax.set_xlim(wide_lo, wide_hi)
ax.tick_params(axis="both", labelsize=6, length=2); ax.set_xticks([])
ax.spines[["top", "right", "bottom"]].set_visible(False)
ax.text(0.98, 0.92, "eGFR GWAS", ha="right", va="top", fontsize=6,
        transform=ax.transAxes, color="#333333")

# ── PIP ───────────────────────────────────────────────────────────────────────
ax = A["pip"]
ncp = pip_win[pip_win["cs"] == -1]
fp  = pip_win[pip_win["hg38_bp"] == VAR_POS]
ocp = pip_win[(pip_win["cs"] != -1) & (pip_win["hg38_bp"] != VAR_POS)]
ax.scatter(ncp["hg38_bp"], ncp["pip"], s=3, color="#9a9a9a",
           alpha=0.5, rasterized=True, linewidths=0)
ax.scatter(ocp["hg38_bp"], ocp["pip"], s=9, color=CS_COLOR,
           alpha=0.95, zorder=4, edgecolors="white", linewidths=0.3)
# focal: same size as CS dots, black outline, drawn on top
ax.scatter(fp["hg38_bp"], fp["pip"], s=9, color=CS_COLOR,
           edgecolors="black", linewidths=0.7, zorder=6)
ax.axvline(VAR_POS, **GRAY_DASH)
ax.set_ylabel("PIP", fontsize=6.5)
ax.set_ylim(0, 1.0); ax.set_yticks([0, 0.5, 1.0]); ax.set_xlim(wide_lo, wide_hi)
ax.tick_params(axis="both", labelsize=6, length=2)
ax.set_xticks([80.10e6, 80.20e6, 80.30e6, 80.40e6])
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:.2f}"))
ax.tick_params(axis="x", labelsize=6, length=2)
ax.spines[["top", "right"]].set_visible(False)

# ── Credible-set model prioritization: 3 side-by-side lollipop panels ──────────
# zoom to the ~20 kb credible-set region; one panel per model (own y-scale).
cs_models = [
    ("csprior_0", "chrombpnet_log2fc",      "ChromBPNet",  "log$_2$FC"),
    ("csprior_1", "atac_logSUM", "Cerberus ATAC", "ΔlogSUM"),
    ("csprior_2", "rna_logSED",  "Cerberus RNA",  "ΔlogSED"),
]
for key, col, mname, metric in cs_models:
    ax = A[key]
    vals = cs_scores[col].values
    ax.axhline(0, color="#cccccc", lw=0.4, zorder=1)
    ax.axvline(VAR_POS, **GRAY_DASH)   # dashed gray marks the highlighted lead
    for _, r in cs_scores.iterrows():
        lead = (r["pos"] == VAR_POS)
        if not lead:   # lead's stem replaced by the dashed gray line
            ax.plot([r["pos"], r["pos"]], [0, r[col]], color=CS_COLOR, lw=0.9, zorder=2)
        ax.scatter([r["pos"]], [r[col]], s=16 if lead else 9, color=CS_COLOR,
                   edgecolors="black" if lead else "white",
                   linewidths=0.7 if lead else 0.3, zorder=3)
    ax.set_xlim(CS_LO, CS_HI)
    lo, hi = min(0.0, vals.min()), max(0.0, vals.max())
    span = hi - lo
    ax.set_ylim(lo - span * 0.08, hi + span * 0.16)
    ax.set_title(mname, fontsize=6, pad=2)
    ax.set_ylabel(metric, fontsize=5, labelpad=1)
    ax.tick_params(axis="y", labelsize=4.5, width=0.4, length=2, pad=1)
    ax.yaxis.set_major_locator(MaxNLocator(3, prune="both"))
    ax.set_xticks([80.246e6, 80.258e6])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:.3f}"))
    ax.tick_params(axis="x", labelsize=4.5, length=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_linewidth(0.4); ax.spines["bottom"].set_linewidth(0.4)
fig.text(CS3_L[0], A["csprior_0"].get_position().y1 + 0.013,
         "Predicted effect across the 7 credible-set variants  (lead rs12509595 outlined in black)",
         ha="left", va="bottom", fontsize=5.2, color="#555555")

# ── ATAC section: ChromBPNet (left) + Cerberus (right), panel G coverage style ──
def draw_cov_panelG(ax, x, ref, alt, var_pos, xlim, ylabel=None, peak=False):
    if peak:
        ax.axvspan(PEAK_START, PEAK_END, color=PEAK_SHADE, alpha=0.8, zorder=0, lw=0)
    ax.axvline(var_pos, color="#888", lw=0.5, ls="--", zorder=2)
    mask = (x >= xlim[0]) & (x <= xlim[1])
    xp, rp, ap = x[mask], ref[mask], alt[mask]
    if rp.sum() <= ap.sum():
        ax.plot(xp, ap, color=ALT_COLOR, lw=0.4, zorder=3); ax.plot(xp, rp, color=REF_COLOR, lw=0.4, zorder=4)
    else:
        ax.plot(xp, rp, color=REF_COLOR, lw=0.4, zorder=3); ax.plot(xp, ap, color=ALT_COLOR, lw=0.4, zorder=4)
    ax.set_xlim(*xlim); ax.set_ylim(bottom=0); ax.set_xticks([])
    ax.spines[["top", "right", "bottom"]].set_visible(False); ax.spines["left"].set_linewidth(0.4)
    ax.tick_params(axis="y", labelsize=4.5, width=0.4, length=2, pad=1)
    ax.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=5.5, rotation=0, ha="right", va="center", labelpad=2)
    ax.legend(handles=[Line2D([0],[0], color=REF_COLOR, lw=0.8, label=f"ref ({VAR_REF})"),
                       Line2D([0],[0], color=ALT_COLOR, lw=0.8, label=f"alt ({VAR_ALT})")],
              fontsize=4.5, loc="upper right", frameon=False,
              handlelength=0.8, handleheight=0.55, borderpad=0, labelspacing=0.2)

# Row 1: coverage
draw_cov_panelG(A["atac_cov_l"], cbp_x, cbp_ref, cbp_alt, VAR_POS, (cbp_lo, cbp_hi),
                ylabel="pred\n(PTS ATAC)", peak=True)
draw_cov_panelG(A["atac_cov_r"], bin_centers, atac_ref, atac_alt, VAR_POS, (atac_lo, atac_hi),
                peak=True)
A["atac_cov_l"].text((PEAK_START + PEAK_END) / 2, 1.0, "peak",
                     transform=A["atac_cov_l"].get_xaxis_transform(),
                     ha="center", va="bottom", fontsize=4.5, color="#777777", clip_on=False)

# Rows 2-3: contribution (SHAP) / ISM logos
draw_logo(A["atac_ref_l"], shap_logo_ref, shap_ylim, SHAP_VAR_IDX, 50)
A["atac_ref_l"].set_ylabel(f"contrib\nref ({VAR_REF})", fontsize=5.5, rotation=0, ha="right", va="center", labelpad=2)
draw_logo(A["atac_alt_l"], shap_logo_alt, shap_ylim, SHAP_VAR_IDX, 50)
A["atac_alt_l"].set_ylabel(f"contrib\nalt ({VAR_ALT})", fontsize=5.5, rotation=0, ha="right", va="center", labelpad=2)
draw_logo(A["atac_ref_r"], atac_logo_ref, atac_ylim, ATAC_VAR_IDX, 50)
draw_logo(A["atac_alt_r"], atac_logo_alt, atac_ylim, ATAC_VAR_IDX, 50)

# Row 4: HNF4A PWM under each column
draw_pwm(A["atac_pwm_l"], shap_pwm_df, SHAP_VAR_IDX, 50, shap_rc, MOTIF_NAME)
A["atac_pwm_l"].set_ylabel("motif", fontsize=5.5, rotation=0, ha="right", va="center", labelpad=2)
draw_pwm(A["atac_pwm_r"], atac_pwm_df, ATAC_VAR_IDX, 50, atac_rc, MOTIF_NAME)

# ── RNA predicted (full) — panel I (merged_panel_i_v2.py) coverage style ──────
ax = A["rna_pred"]
xpr, rpr, apr = bin_centers[rna_mask], rna_ref[rna_mask], rna_alt[rna_mask]
ytop = max(rpr.max(), apr.max()) * 1.15
# larger-sum allele drawn on bottom (panel I convention)
if rpr.sum() >= apr.sum():
    ax.plot(xpr, rpr, color=REF_COLOR, lw=0.5, zorder=3); ax.plot(xpr, apr, color=ALT_COLOR, lw=0.5, zorder=4)
else:
    ax.plot(xpr, apr, color=ALT_COLOR, lw=0.5, zorder=3); ax.plot(xpr, rpr, color=REF_COLOR, lw=0.5, zorder=4)
ax.legend(
    handles=[Line2D([0],[0], color=REF_COLOR, lw=0.8, label=f"ref ({VAR_REF})  Σ={rna_rsum:.1f}"),
             Line2D([0],[0], color=ALT_COLOR, lw=0.8, label=f"alt ({VAR_ALT})  Σ={rna_asum:.1f}")],
    fontsize=3.8, loc="upper right", frameon=False,
    handlelength=0.8, handleheight=0.55, borderpad=0, labelspacing=0.2)
ax.axvline(VAR_POS, color="#888", lw=0.7, ls="--", zorder=7)
ax.text(VAR_POS, ytop * 0.96, f"{VAR_CHROM}:{VAR_POS:,} {VAR_REF}→{VAR_ALT}",
        ha="left", va="top", fontsize=3.8, color="#888")
ax.set_xlim(RNA_PRED_START, RNA_PRED_END); ax.set_ylim(0, ytop); ax.set_xticks([])
_ymax = max(rpr.max(), apr.max())
ax.set_yticks([0, _ymax]); ax.set_yticklabels(["0", f"{_ymax:.1f}"])
ax.tick_params(axis="y", labelsize=4.5, width=0.4, length=2, pad=1)
ax.text(0.01, 0.03, f"{VAR_CHROM}:{RNA_PRED_START:,}–{RNA_PRED_END:,}",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=3.8, color="#888")
ax.spines[["top", "right", "bottom"]].set_visible(False)
ax.spines["left"].set_linewidth(0.4)
ax.set_ylabel("Cerberus pred\n(PTS RNA+)", fontsize=5.5, rotation=0, ha="right",
              va="center", labelpad=2)

# ── Gene annotation (full) ────────────────────────────────────────────────────
ax = A["gene"]
ax.set_xlim(RNA_PRED_START, RNA_PRED_END); ax.set_ylim(0, 1)
def draw_gene(g_start, g_end, exons, label, color="#333333"):
    if g_end < RNA_PRED_START or g_start > RNA_PRED_END: return
    s = max(g_start, RNA_PRED_START); e = min(g_end, RNA_PRED_END)
    ax.plot([s, e], [0.5, 0.5], color=color, lw=1.0, solid_capstyle="butt")
    for es, ee in exons:
        if ee < RNA_PRED_START or es > RNA_PRED_END: continue
        ax.add_patch(Rectangle((max(es, RNA_PRED_START), 0.5 - 0.19),
                               min(ee, RNA_PRED_END) - max(es, RNA_PRED_START), 0.38,
                               color=color, zorder=3))
    ax.text(s + (RNA_PRED_END - RNA_PRED_START) * 0.01, 0.80, label,
            fontsize=5.5, va="bottom", color=color)
draw_gene(fgf5_start, fgf5_end, fgf5_exons, "FGF5")
draw_gene(cfap_start, cfap_end, cfap_exons[:12], "CFAP299", color="#555555")
ax.axvline(VAR_POS, **GRAY_DASH)
ax.spines[["top", "right", "left"]].set_visible(False); ax.set_yticks([])
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:.3f}"))
ax.tick_params(axis="x", labelsize=5.5, length=2)

# ── RNA ISM (full) ────────────────────────────────────────────────────────────
for key, logo_df, label in [("rna_ism_ref", rna_logo_ref, f"Ref ({VAR_REF})"),
                            ("rna_ism_alt", rna_logo_alt, f"Alt ({VAR_ALT})")]:
    ax = A[key]
    draw_logo(ax, logo_df, rna_ylim, RNA_VAR_IDX, 100)
    ax.set_ylabel(label, fontsize=5.5, rotation=0, ha="right", va="center", labelpad=2)
A["rna_ism_ref"].set_title("RNA ISM / FGF5 exons (PTS)", fontsize=6, pad=2, loc="left")
draw_pwm(A["rna_pwm"], rna_pwm_df, RNA_VAR_IDX, 100, rna_rc, MOTIF_NAME)
A["rna_pwm"].set_ylabel("motif", fontsize=5.5, rotation=0, ha="right", va="center", labelpad=2)

# ── Save ──────────────────────────────────────────────────────────────────────
print("Saving...")
os.makedirs(OUT_DIR, exist_ok=True)
fig.savefig(OUT_PDF, dpi=200)
fig.savefig(OUT_PNG, dpi=200)
print(f"Saved: {OUT_PDF}\nSaved: {OUT_PNG}")
