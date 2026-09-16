#!/usr/bin/env python
"""eQTL examples: association, fine-mapping, predicted expression and motif.

One column per example variant-gene pair. Rows, top to bottom:

    eQTL association, -log10 p over +/-100 kb
    SuSiE fine-mapping PIP over the same window
    Cerberus predicted RNA coverage on both alleles, with exonic sums
    gene model
    Cerberus ISM, reference allele
    Cerberus ISM, alternate allele
    the aligned JASPAR motif

The two examples run in opposite directions: at one the reference allele
predicts higher expression, at the other the alternate does.

Usage:
    python eqtl_examples.py [--out_dir <dir>]
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
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, FuncFormatter

from figlib import (
    ALT_COLOR, CT_COLORS, HIGHLIGHT, NT_COLORS, REF_COLOR,
    align_pwm, cfg, ism_to_logo, parse_gtf_gene, parse_pwm, setup_style,
)

setup_style()

# ── Paths ──────────────────────────────────────────────────────────────────────
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--out_dir", default=None)
_args = _parser.parse_args()

WORK_DIR = cfg("WORK_DIR")
OUT_DIR  = _args.out_dir or cfg("FIGURE_DIR")

# ISM over the concordant fine-mapped variants, averaged across folds.
ISM_H5   = os.path.join(WORK_DIR, "eqtl", "examples", "cerberus_ism",
                        "ensemble", "scores.h5")
# Cached predicted coverage, one npz per example.
PRED_DIR = os.path.join(WORK_DIR, "eqtl", "examples", "cerberus_pred")
GTF_FILE = cfg("EQTL_GTF")
EQTL_TSV = cfg("EQTL_SUSIE_TSV")
OUT_PDF  = os.path.join(OUT_DIR, "eqtl_examples.pdf")
OUT_PNG  = os.path.join(OUT_DIR, "eqtl_examples.png")

MANHATTAN_FLANK = 100_000

# ── Per-example constants ──────────────────────────────────────────────────────
COL1A2AS1 = dict(
    chrom="chr7", pos=94392163, ref="T", alt="C",
    gene_name="ENSG00000285090", gene_id="ENSG00000285090",
    phenotype_id="ENSG00000285090",
    ct="EC", logsed=0.604, slope=0.974, pip=0.981,
    ism_track=7,
    motif_id="MA0517.1", motif_name="STAT1::STAT2",
    npz=os.path.join(PRED_DIR, "chr7_94392163_EC.npz"),
    flank=2000,
)
ZNF559 = dict(
    chrom="chr19", pos=9324196, ref="C", alt="T",
    gene_name="ZNF559", gene_id="ENSG00000188321",
    phenotype_id="ZNF559",
    ct="IC", logsed=-0.5146, slope=-1.133, pip=0.822,
    ism_track=9,
    motif_id="MA0095.3", motif_name="YY1",
    npz=os.path.join(PRED_DIR, "chr19_9324196_IC.npz"),
    flank=2000,
    gene_start=9323772, gene_end=9351162, strand="+",
    exons=[(9324179,9324228),(9324695,9324780),(9337796,9337858),
           (9338494,9338582),(9339193,9339319),(9341102,9341184),
           (9341695,9345871)],
)

# ── Helpers ────────────────────────────────────────────────────────────────────
# ── Load eQTL data ─────────────────────────────────────────────────────────────
print("Loading eQTL data...")
eqtl = pd.read_csv(EQTL_TSV, sep="\t", compression="gzip",
                   usecols=["celltype","variant_id","phenotype_id","position",
                             "pval_nominal","variable_prob","cs"])

def load_locus(ex):
    sub = eqtl[(eqtl["phenotype_id"] == ex["phenotype_id"]) &
               (eqtl["celltype"] == ex["ct"]) &
               (eqtl["position"] >= ex["pos"] - MANHATTAN_FLANK) &
               (eqtl["position"] <= ex["pos"] + MANHATTAN_FLANK)].copy()
    sub["neglog10p"] = -np.log10(sub["pval_nominal"].clip(lower=1e-300))
    return sub

locus_d = load_locus(COL1A2AS1)
locus_z = load_locus(ZNF559)
print(f"  COL1A2-AS1 locus: {len(locus_d)} variants")
print(f"  ZNF559 locus: {len(locus_z)} variants")

# ── Load ISM ───────────────────────────────────────────────────────────────────
print("Loading ISM data...")
with h5py.File(ISM_H5) as h5:
    labels   = h5["label"][:].astype(str)
    starts   = h5["start"][:]
    gene_ids = h5["gene_ids"][:].astype(str)
    snp_idx  = h5["snp_idx"][:]
    gene_idx = h5["gene_idx"][:]
    pair_lut = {(int(snp_idx[p]), int(gene_idx[p])): p for p in range(len(snp_idx))}
    gid_base = np.array([g.split(".")[0] for g in gene_ids])

    def _load(ex, label):
        si  = int(np.where(labels == label)[0][0])
        gi  = int(np.where(gid_base == ex["gene_id"])[0][0])
        pi  = pair_lut[(si, gi)]
        oh  = h5["ref"]["seqs"][si].T.astype(float)
        ir  = h5["ref"]["covgene"]["logSED"][pi, :, :, ex["ism_track"]]
        ia  = h5["alt"]["covgene"]["logSED"][pi, :, :, ex["ism_track"]]
        voff = ex["pos"] - 1 - int(starts[si])
        return oh, ir, ia, int(starts[si]), voff

    oh_d, ir_d, ia_d, ws_d, vo_d = _load(COL1A2AS1, f"chr7:{COL1A2AS1['pos']}")
    oh_z, ir_z, ia_z, ws_z, vo_z = _load(ZNF559,    f"chr19:{ZNF559['pos']}")

lr_d = ism_to_logo(ir_d, oh_d); la_d = ism_to_logo(ia_d, oh_d)
lr_z = ism_to_logo(ir_z, oh_z); la_z = ism_to_logo(ia_z, oh_z)

# ── PWM alignment ─────────────────────────────────────────────────────────────
print("Aligning PWMs...")
pwm_d, rc_d, _ = align_pwm(oh_d, parse_pwm(COL1A2AS1["motif_id"]), vo_d, 100)
pwm_z, rc_z, _ = align_pwm(oh_z, parse_pwm(ZNF559["motif_id"]), vo_z, 100)
print(f"  STAT1::STAT2 RC={rc_d}  YY1 RC={rc_z}")

# ── Predicted coverage ─────────────────────────────────────────────────────────
print("Loading predictions...")
def _load_npz(ex):
    npz = np.load(ex["npz"], allow_pickle=True)
    bs  = int(npz["bins_start"]); bsz = int(npz["bin_size"])
    xf  = bs + np.arange(len(npz["pred_ref"])) * bsz + bsz // 2
    return npz, xf

npz_d, xf_d = _load_npz(COL1A2AS1)
npz_z, xf_z = _load_npz(ZNF559)

# ── Gene annotation + coverage window ─────────────────────────────────────────
print("Parsing GTF (COL1A2-AS1)...")
tx_s_d, tx_e_d, exons_d, strand_d = parse_gtf_gene(GTF_FILE, COL1A2AS1["gene_name"])
win_s_d = max(1, min(tx_s_d, COL1A2AS1["pos"]) - COL1A2AS1["flank"])
win_e_d = max(tx_e_d, COL1A2AS1["pos"]) + COL1A2AS1["flank"]

tx_s_z, tx_e_z = ZNF559["gene_start"], ZNF559["gene_end"]
exons_z = ZNF559["exons"]; strand_z = ZNF559["strand"]
win_s_z = min(tx_s_z, ZNF559["pos"]) - ZNF559["flank"]
win_e_z = max(tx_e_z, ZNF559["pos"]) + ZNF559["flank"]

def exon_sums(npz, xf, exons):
    r = s = 0.0
    for es, ee in exons:
        m = (xf >= es) & (xf <= ee)
        r += npz["pred_ref"].astype(float)[m].sum()
        s += npz["pred_alt"].astype(float)[m].sum()
    return r, s

rs_d, as_d = exon_sums(npz_d, xf_d, exons_d)
rs_z, as_z = exon_sums(npz_z, xf_z, exons_z)

# ── Figure ─────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(7.8, 2.8))
gs  = gridspec.GridSpec(7, 2,
                        height_ratios=[0.65, 0.35, 1.2, 0.22, 0.6, 0.6, 0.24],
                        hspace=0.25, wspace=0.38)

ax_man  = [fig.add_subplot(gs[0, c]) for c in range(2)]
ax_pip  = [fig.add_subplot(gs[1, c], sharex=ax_man[c]) for c in range(2)]
ax_cov  = [fig.add_subplot(gs[2, c]) for c in range(2)]
ax_gene = [fig.add_subplot(gs[3, c], sharex=ax_cov[c]) for c in range(2)]
ax_ir   = [fig.add_subplot(gs[4, c]) for c in range(2)]
ax_ia   = [fig.add_subplot(gs[5, c], sharex=ax_ir[c]) for c in range(2)]
ax_pw   = [fig.add_subplot(gs[6, c], sharex=ax_ir[c]) for c in range(2)]

examples = [
    dict(ex=COL1A2AS1, locus=locus_d,
         npz=npz_d, xf=xf_d, win_s=win_s_d, win_e=win_e_d,
         tx_s=tx_s_d, tx_e=tx_e_d, exons=exons_d, strand=strand_d,
         rs=rs_d, as_=as_d,
         lr=lr_d, la=la_d, voff=vo_d, win_start=ws_d,
         pwm=pwm_d, rc=rc_d,
         ax_man=ax_man[0], ax_pip=ax_pip[0],
         ax_cov=ax_cov[0], ax_gene=ax_gene[0],
         ax_ir=ax_ir[0],   ax_ia=ax_ia[0], ax_pw=ax_pw[0]),
    dict(ex=ZNF559, locus=locus_z,
         npz=npz_z, xf=xf_z, win_s=win_s_z, win_e=win_e_z,
         tx_s=tx_s_z, tx_e=tx_e_z, exons=exons_z, strand=strand_z,
         rs=rs_z, as_=as_z,
         lr=lr_z, la=la_z, voff=vo_z, win_start=ws_z,
         pwm=pwm_z, rc=rc_z,
         ax_man=ax_man[1], ax_pip=ax_pip[1],
         ax_cov=ax_cov[1], ax_gene=ax_gene[1],
         ax_ir=ax_ir[1],   ax_ia=ax_ia[1], ax_pw=ax_pw[1]),
]

for D in examples:
    ex      = D["ex"]
    locus   = D["locus"]
    ct_col  = CT_COLORS[ex["ct"]]
    man_xlim = (ex["pos"] - MANHATTAN_FLANK, ex["pos"] + MANHATTAN_FLANK)

    # ── Row 0: Manhattan ──────────────────────────────────────────────────────
    a_man = D["ax_man"]
    bg    = locus[locus["cs"] != 1]
    cs1   = locus[locus["cs"] == 1]
    a_man.scatter(bg["position"],  bg["neglog10p"],  s=1,  color="#bbbbbb",
                  alpha=0.6, lw=0, rasterized=True, zorder=2)
    a_man.scatter(cs1["position"], cs1["neglog10p"], s=4,  color=ct_col,
                  alpha=0.85, lw=0, zorder=3)
    focal_man = locus[locus["position"] == ex["pos"]]
    if len(focal_man):
        a_man.scatter(focal_man["position"], focal_man["neglog10p"],
                      s=18, color=ct_col, edgecolors="black",
                      linewidths=0.5, zorder=10)
    a_man.axvline(ex["pos"], color="#888", lw=0.5, ls="--", zorder=1)
    a_man.set_xlim(*man_xlim)
    a_man.set_ylim(bottom=0)
    a_man.spines[["top","right","bottom"]].set_visible(False)
    a_man.spines["left"].set_linewidth(0.4)
    a_man.tick_params(axis="x", labelbottom=False, bottom=False)
    a_man.tick_params(axis="y", labelsize=4, width=0.4, length=2, pad=1)
    a_man.yaxis.set_major_locator(MaxNLocator(3, prune="both"))
    a_man.set_ylabel("−log₁₀(p)", fontsize=4, labelpad=2)
    a_man.set_title(
        f"{ex['gene_name']}  {ex['ct']}  logSED={ex['logsed']:+.3f}  "
        f"β={ex['slope']:+.3f}  PIP={ex['pip']:.3f}",
        fontsize=5, pad=2)

    # ── Row 1: PIP ────────────────────────────────────────────────────────────
    a_pip = D["ax_pip"]
    a_pip.scatter(bg["position"],  bg["variable_prob"],  s=1, color="#bbbbbb",
                  alpha=0.6, lw=0, rasterized=True, zorder=2)
    a_pip.scatter(cs1["position"], cs1["variable_prob"], s=4, color=ct_col,
                  alpha=0.85, lw=0, zorder=3)
    focal_pip = locus[locus["position"] == ex["pos"]]
    if len(focal_pip):
        a_pip.scatter(focal_pip["position"], focal_pip["variable_prob"],
                      s=18, color=ct_col, edgecolors="black",
                      linewidths=0.5, zorder=10)
    a_pip.axvline(ex["pos"], color="#888", lw=0.5, ls="--", zorder=1)
    a_pip.set_xlim(*man_xlim)
    a_pip.set_ylim(0, 1.08)
    a_pip.spines[["top","right"]].set_visible(False)
    a_pip.spines[["left","bottom"]].set_linewidth(0.4)
    a_pip.tick_params(axis="x", labelsize=4, width=0.4, length=2, pad=1)
    a_pip.tick_params(axis="y", labelsize=4, width=0.4, length=2, pad=1)
    a_pip.yaxis.set_major_locator(MaxNLocator(3, prune="upper"))
    a_pip.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:.2f}"))
    a_pip.xaxis.set_major_locator(MaxNLocator(5, prune="both"))
    a_pip.set_ylabel("PIP", fontsize=4, labelpad=2)
    a_pip.set_xlabel("Position (Mb)", fontsize=4, labelpad=1)

    # ── Row 2: coverage ───────────────────────────────────────────────────────
    a_cov = D["ax_cov"]
    mask  = (D["xf"] >= D["win_s"]) & (D["xf"] <= D["win_e"])
    xp    = D["xf"][mask]
    rp    = D["npz"]["pred_ref"].astype(float)[mask]
    ap    = D["npz"]["pred_alt"].astype(float)[mask]
    ytop  = max(rp.max(), ap.max()) * 1.15
    if rp.sum() >= ap.sum():
        a_cov.plot(xp, rp, color=REF_COLOR, lw=0.5, zorder=3)
        a_cov.plot(xp, ap, color=ALT_COLOR, lw=0.5, zorder=4)
    else:
        a_cov.plot(xp, ap, color=ALT_COLOR, lw=0.5, zorder=3)
        a_cov.plot(xp, rp, color=REF_COLOR, lw=0.5, zorder=4)
    a_cov.legend(
        handles=[Line2D([0],[0], color=REF_COLOR, lw=0.8, label=f"ref ({ex['ref']})  Σ={D['rs']:.1f}"),
                 Line2D([0],[0], color=ALT_COLOR, lw=0.8, label=f"alt ({ex['alt']})  Σ={D['as_']:.1f}")],
        fontsize=3.8, loc="upper right", frameon=False,
        handlelength=0.8, handleheight=0.55, borderpad=0, labelspacing=0.2)
    a_cov.set_ylim(0, ytop)
    a_cov.axvline(ex["pos"], color="#888", lw=0.7, ls="--", zorder=7)
    a_cov.text(ex["pos"], ytop*0.96,
               f"{ex['chrom']}:{ex['pos']:,} {ex['ref']}→{ex['alt']}",
               ha="left", va="top", fontsize=3.8, color="#888")
    a_cov.set_xlim(D["win_s"], D["win_e"])
    _ymax = max(rp.max(), ap.max())
    a_cov.set_yticks([0, _ymax])
    a_cov.set_yticklabels(["0", f"{_ymax:.1f}"])
    a_cov.tick_params(axis="y", labelsize=4.5, width=0.4, length=2, pad=1)
    a_cov.tick_params(axis="x", labelbottom=False, bottom=False)
    a_cov.text(0.01, 0.03, f"{ex['chrom']}:{D['win_s']:,}–{D['win_e']:,}",
               transform=a_cov.transAxes, ha="left", va="bottom",
               fontsize=3.8, color="#888")
    a_cov.spines[["top","right","bottom"]].set_visible(False)
    a_cov.spines["left"].set_linewidth(0.4)
    a_cov.set_ylabel(f"Cerberus pred\n({ex['ct']} RNA+)", fontsize=4.5, labelpad=2)

    # ── Row 3: gene annotation ────────────────────────────────────────────────
    a_gene = D["ax_gene"]
    a_gene.set_xlim(D["win_s"], D["win_e"])
    a_gene.set_ylim(0, 1)
    a_gene.axvline(ex["pos"], color="#888", lw=0.7, ls="--", zorder=5)
    iy = 0.5
    a_gene.plot([D["tx_s"], D["tx_e"]], [iy, iy], color="#333", lw=0.4, zorder=2)
    for es, ee in D["exons"]:
        a_gene.add_patch(Rectangle((es, iy-0.25), ee-es, 0.5,
                                   fc="#333", ec="none", zorder=3))
    arw = max((D["tx_e"] - D["tx_s"]) // 8, 3000)
    for apos in range(D["tx_s"] + arw, D["tx_e"], arw):
        dx = 1 if D["strand"] == "+" else -1
        a_gene.annotate("", xy=(apos+dx, iy), xytext=(apos-dx, iy),
            arrowprops=dict(arrowstyle="-|>", color="#333", lw=0.4,
                            mutation_scale=3, shrinkA=0, shrinkB=0), zorder=4)
    a_gene.text((D["tx_s"]+D["tx_e"])/2, iy+0.32, ex["gene_name"],
                ha="center", va="bottom", fontsize=5, style="italic", color="#333")
    a_gene.set_yticks([]); a_gene.spines[:].set_visible(False)
    a_gene.tick_params(axis="x", labelbottom=False, bottom=False)

    # ── ISM shared scale ──────────────────────────────────────────────────────
    _abs = max(abs(D["lr"].values.max()), abs(D["lr"].values.min()),
               abs(D["la"].values.max()), abs(D["la"].values.min()))
    ra_ylim  = (-_abs * 1.15, _abs * 1.15)
    ism_xlim = (-0.5, 99.5)
    voff     = D["voff"]

    def _style_ism(ax, ylim, ylabel):
        ax.axvspan(voff-0.5, voff+0.5, color=HIGHLIGHT, alpha=0.85, zorder=0, lw=0)
        ax.axvline(voff, color="#888", lw=0.5, ls="--", zorder=1)
        ax.axhline(0, color="#bbb", lw=0.35, zorder=1)
        ax.set_xlim(*ism_xlim); ax.set_ylim(*ylim); ax.set_xticks([])
        ax.yaxis.set_major_locator(MaxNLocator(3, prune="both"))
        ax.tick_params(axis="y", labelsize=4.5, width=0.4, length=2, pad=1)
        ax.spines[["top","right","bottom"]].set_visible(False)
        ax.spines["left"].set_linewidth(0.4)
        ax.set_ylabel(ylabel, fontsize=4.5, labelpad=2)

    # Rows 4-5: ISM ref / alt
    for ax_ism, logo, lbl in [
        (D["ax_ir"], D["lr"], f"ISM ref\n({ex['ct']} RNA+)"),
        (D["ax_ia"], D["la"], f"ISM alt\n({ex['ct']} RNA+)"),
    ]:
        logomaker.Logo(logo, color_scheme=NT_COLORS, ax=ax_ism,
                       flip_below=True, baseline_width=0, font_name="DejaVu Sans", zorder=2)
        _style_ism(ax_ism, ra_ylim, lbl)
    D["ax_ir"].set_title(f"{ex['chrom']}:{D['win_start']+1:,}–{D['win_start']+100:,}",
                         fontsize=4, pad=2)

    # Row 6: motif PWM
    a_pw = D["ax_pw"]
    rc_tag = " (−)" if D["rc"] else ""
    a_pw.axvspan(voff-0.5, voff+0.5, color=HIGHLIGHT, alpha=0.85, zorder=0, lw=0)
    a_pw.axvline(voff, color="#888", lw=0.5, ls="--", zorder=1)
    logomaker.Logo(D["pwm"], color_scheme=NT_COLORS, ax=a_pw,
                   flip_below=False, baseline_width=0, font_name="DejaVu Sans", zorder=2)
    a_pw.text(0.5, -0.28, f"{ex['motif_name']}{rc_tag}",
              transform=a_pw.transAxes, ha="center", va="top",
              fontsize=4.5, style="italic", color="#333")
    a_pw.set_ylabel(f"{ex['motif_name']}\n(bits)", fontsize=4.5, labelpad=2)
    a_pw.set_xlim(*ism_xlim); a_pw.set_ylim(bottom=0)
    a_pw.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
    a_pw.tick_params(axis="y", labelsize=4.5, width=0.4, length=2, pad=1)
    a_pw.spines[["top","right","bottom"]].set_visible(False)
    a_pw.spines["left"].set_linewidth(0.4)
    a_pw.set_xticks([])

# ── Reposition ISM rows to 60 % of each column width ─────────────────────────
fig.canvas.draw()
for c in range(2):
    pos_gene = ax_gene[c].get_position()
    pos_r    = ax_ir[c].get_position()
    pos_a    = ax_ia[c].get_position()
    pos_p    = ax_pw[c].get_position()
    ism_w = pos_r.width * 0.60
    ism_x = pos_r.x0 + (pos_r.width - ism_w) / 2
    gap   = 0.018
    bot_r = pos_gene.y0 - gap - pos_r.height
    bot_a = bot_r - gap * 0.5 - pos_a.height
    bot_p = bot_a - gap        - pos_p.height
    ax_ir[c].set_position([ism_x, bot_r, ism_w, pos_r.height])
    ax_ia[c].set_position([ism_x, bot_a, ism_w, pos_a.height])
    ax_pw[c].set_position([ism_x, bot_p, ism_w, pos_p.height])

os.makedirs(OUT_DIR, exist_ok=True)
fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")
