#!/usr/bin/env python
"""caQTL examples: predicted coverage, attributions and the disrupted motif.

One column per example variant, each in its own target cell type. Rows, top to
bottom:

    ChromBPNet predicted accessibility   reference red, alternate blue
    ChromBPNet SHAP, reference allele
    ChromBPNet SHAP, alternate allele
    Cerberus predicted accessibility
    Cerberus ISM, reference allele
    Cerberus ISM, alternate allele
    the aligned JASPAR motif

Both models' attributions are shown on both alleles so the loss of signal at the
variant is visible rather than inferred.

Examples are listed in a TSV (`--examples`) so they can be changed without
editing this file. Columns:

    variant_key  chrom_pos_ref_alt, matching the prediction/ISM filenames
    label        what to print in the column header
    motif_id     JASPAR accession
    motif_name   display name
    chrom
    pos          1-based
    celltype     target cell type, QTL naming

Usage:
    python caqtl_examples.py --examples examples_caqtl.tsv [--out_dir <dir>]
"""
import argparse
import os

import h5py
import logomaker
import numpy as np
import pandas as pd

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from figlib import (
    ALT_COLOR,
    ATAC_SUBSET_IDX,
    CT_COLORS,
    CT_LABELS,
    HIGHLIGHT,
    NT_COLORS,
    REF_COLOR,
    align_pwm,
    cfg,
    parse_pwm,
    setup_style,
)

ISM_SHADE = "#e0e0e0"   # band marking the mutagenised window on coverage tracks

COV_WIN_BP = 500        # bases either side of the variant shown in coverage
ISM_WIN_BP = 25         # half-width of the shaded ISM window
ISM_LEN = 50            # bases in the attribution windows

# The variant's index within each model's saved window. They differ by one
# because the two tools centre their windows differently; getting this wrong
# shifts every logo by a base.
CHROMBPNET_VAR_POS = 25
CERBERUS_VAR_POS = 24

ROW_LABELS = [
    "ChromBPNet\npredicted",
    "ChromBPNet\nSHAP (ref)",
    "ChromBPNet\nSHAP (alt)",
    "Cerberus\npredicted",
    "Cerberus ISM\n(ref)",
    "Cerberus ISM\n(alt)",
    "PWM",
]

COV_H, ISM_H, PWM_H = 0.35, 0.343, 0.315


def load_chrombpnet_shap(shap_dir, celltypes):
    """Per-cell-type SHAP attributions, keyed by variant."""
    loaded = {}
    for celltype in celltypes:
        path = os.path.join(shap_dir, f"{celltype}_ism.h5")
        if not os.path.exists(path):
            print(f"  missing ChromBPNet SHAP for {celltype}: {path}")
            continue
        with h5py.File(path, "r") as handle:
            loaded[celltype] = {
                "ids": list(handle["variant_ids"][:].astype(str)),
                "ref_contrib": handle["ref_contrib"][:],
                "alt_contrib": handle["alt_contrib"][:],
                "ref_onehot": handle["ref_onehot"][:],
                "alt_onehot": handle["alt_onehot"][:],
            }
    return loaded


def chrombpnet_logos(shap, celltype, variant_key):
    """Attribution logos for both alleles: contribution times the base present."""
    entry = shap.get(celltype)
    if entry is None or variant_key not in entry["ids"]:
        return None, None, None
    index = entry["ids"].index(variant_key)
    ref_onehot = entry["ref_onehot"][index]
    alt_onehot = entry["alt_onehot"][index]
    return (
        pd.DataFrame(ref_onehot * entry["ref_contrib"][index], columns=list("ACGT")),
        pd.DataFrame(alt_onehot * entry["alt_contrib"][index], columns=list("ACGT")),
        ref_onehot,
    )


def load_cerberus_ism(path):
    """ISM scores and sequences per variant from a hound_ism_snp output."""
    with h5py.File(path) as handle:
        labels = list(handle["label"][:].astype(str))
        data = {
            label: (
                handle["ref/seqs"][i],
                handle["alt/seqs"][i],
                handle["ref/cov/logSUM"][i],
                handle["alt/cov/logSUM"][i],
            )
            for i, label in enumerate(labels)
        }
    return data


def ism_logo(sequence, scores, track):
    """Per-position contribution: the present base minus the mean alternative."""
    # sequences are stored (4, L); transpose to (L, 4) before taking the base.
    present = np.argmax(sequence.T, axis=1)
    frame = pd.DataFrame(0.0, index=range(ISM_LEN), columns=list("ACGT"))
    for position in range(ISM_LEN):
        base = present[position]
        others = [j for j in range(4) if j != base]
        frame.loc[position, "ACGT"[base]] = float(
            scores[position, base, track] - np.mean(scores[position, others, track])
        )
    return frame


def cerberus_logos(ism, variant_key, celltype):
    if celltype not in ATAC_SUBSET_IDX or variant_key not in ism:
        return None, None
    track = ATAC_SUBSET_IDX[celltype]
    ref_seq, alt_seq, ref_scores, alt_scores = ism[variant_key]
    return ism_logo(ref_seq, ref_scores, track), ism_logo(alt_seq, alt_scores, track)


def style_row(ax, show_yticks=True):
    ax.set_xticks([])
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.tick_params(axis="y", labelsize=5.5, width=0.4, length=1.5, pad=1)
    ax.yaxis.set_major_locator(MaxNLocator(2, prune="both"))
    if not show_yticks:
        ax.set_yticklabels([])
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)


def set_row_label(ax, column, label):
    if column == 0 and label:
        ax.set_ylabel(
            label, fontsize=5.5, rotation=0, ha="right", va="center", labelpad=40
        )


def draw_coverage(ax, x, ref, alt, var_pos, xlim, column, label):
    ax.axvspan(
        var_pos - ISM_WIN_BP, var_pos + ISM_WIN_BP,
        color=ISM_SHADE, alpha=0.55, zorder=0, lw=0,
    )
    ax.axvline(var_pos, color="#888", lw=0.5, ls="--", zorder=2)

    if ref is not None and x is not None:
        inside = (x >= xlim[0]) & (x <= xlim[1])
        xs, ref_in, alt_in = x[inside], ref[inside], alt[inside]
        # Draw the lower-coverage allele first so the difference stays visible.
        if ref_in.sum() <= alt_in.sum():
            ax.plot(xs, alt_in, color=ALT_COLOR, lw=0.4, zorder=3)
            ax.plot(xs, ref_in, color=REF_COLOR, lw=0.4, zorder=4)
        else:
            ax.plot(xs, ref_in, color=REF_COLOR, lw=0.4, zorder=3)
            ax.plot(xs, alt_in, color=ALT_COLOR, lw=0.4, zorder=4)
        ax.set_ylim(bottom=0)

    ax.set_xlim(*xlim)
    style_row(ax)
    set_row_label(ax, column, label)


def draw_attribution(ax, logo, ylim, column, label, var_pos):
    if logo is not None:
        logomaker.Logo(
            logo, color_scheme=NT_COLORS, ax=ax,
            flip_below=True, baseline_width=0, font_name="DejaVu Sans",
        )
        for patch in ax.patches:
            patch.set_zorder(3)
        ax.axvspan(
            var_pos - 0.5, var_pos + 0.5,
            color=HIGHLIGHT, alpha=0.6, zorder=0, lw=0,
        )
        ax.axvline(var_pos, color="#888", lw=0.5, ls="--", zorder=1)
        ax.axhline(0, color="#bbb", lw=0.4, zorder=1)
    else:
        ax.text(
            0.5, 0.5, "n/a", transform=ax.transAxes,
            ha="center", va="center", fontsize=6, color="#aaa",
        )
    ax.set_xlim(-0.5, ISM_LEN - 0.5)
    ax.set_ylim(*ylim)
    style_row(ax)
    set_row_label(ax, column, label)


def shared_ylim(*logos, pad=1.15):
    values = np.concatenate(
        [logo.values.ravel() if logo is not None else [0.0] for logo in logos]
    )
    return (min(0, values.min()) * pad, max(0, values.max()) * pad)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", required=True, help="example definition TSV")
    parser.add_argument("--out_dir", default=None)
    args = parser.parse_args()

    setup_style()
    out_dir = args.out_dir or cfg("FIGURE_DIR")
    os.makedirs(out_dir, exist_ok=True)

    work_dir = cfg("WORK_DIR")
    example_dir = os.path.join(work_dir, "caqtl", "examples")
    chrombpnet_pred_dir = os.path.join(example_dir, "chrombpnet_pred")
    cerberus_pred_dir = os.path.join(example_dir, "cerberus_pred")
    cerberus_ism_h5 = os.path.join(example_dir, "cerberus_ism", "ensemble", "scores.h5")
    shap_dir = os.path.join(work_dir, "chrombpnet", "caqtl", "shap")

    examples = pd.read_csv(args.examples, sep="\t", comment="#")
    print(f"{len(examples)} examples from {args.examples}")

    shap = load_chrombpnet_shap(shap_dir, examples["celltype"].unique())
    cerberus_ism = load_cerberus_ism(cerberus_ism_h5)
    motifs = {mid: parse_pwm(mid) for mid in examples["motif_id"].unique()}

    n_columns = len(examples)
    fig = plt.figure(
        figsize=(n_columns * 1.7, (COV_H * 2 + ISM_H * 4 + PWM_H) * 1.05)
    )
    outer = gridspec.GridSpec(1, n_columns, figure=fig, wspace=0.18)

    for column, example in examples.reset_index(drop=True).iterrows():
        variant_key = example["variant_key"]
        celltype = example["celltype"]
        var_pos = int(example["pos"])

        # ── Predicted coverage ────────────────────────────────────────────────
        chrombpnet_npz = os.path.join(chrombpnet_pred_dir, f"{variant_key}.npz")
        cerberus_npz = os.path.join(cerberus_pred_dir, f"{variant_key}.npz")

        if os.path.exists(chrombpnet_npz):
            data = np.load(chrombpnet_npz)
            out_start = int(data["out_start"][0])
            out_len = int(data["out_len"][0])
            cbp_x = out_start + np.arange(out_len) + 0.5
            cbp_xlim = (out_start, out_start + out_len)
            cbp_ref = data[f"{celltype}_ref"]
            cbp_alt = data[f"{celltype}_alt"]
        else:
            print(f"  missing ChromBPNet prediction: {chrombpnet_npz}")
            cbp_x = cbp_ref = cbp_alt = None
            cbp_xlim = (var_pos - COV_WIN_BP, var_pos + COV_WIN_BP)

        if os.path.exists(cerberus_npz):
            data = np.load(cerberus_npz)
            bins_start = int(data["bins_start"])
            bin_size = int(data["bin_size"])
            ref = data[f"{celltype}_ref"]
            cer_x = bins_start + np.arange(len(ref)) * bin_size + bin_size // 2
            cer_ref = ref
            cer_alt = data[f"{celltype}_alt"]
        else:
            print(f"  missing Cerberus prediction: {cerberus_npz}")
            cer_x = cer_ref = cer_alt = None
        cer_xlim = (var_pos - COV_WIN_BP, var_pos + COV_WIN_BP)

        # ── Attributions ──────────────────────────────────────────────────────
        cbp_ref_logo, cbp_alt_logo, cbp_onehot = chrombpnet_logos(
            shap, celltype, variant_key
        )
        cer_ref_logo, cer_alt_logo = cerberus_logos(
            cerberus_ism, variant_key, celltype
        )

        # Reference and alternate share a y-scale, so the drop at the variant is
        # a real difference rather than an artefact of independent scaling.
        cbp_ylim = shared_ylim(cbp_ref_logo, cbp_alt_logo)
        cer_ylim = shared_ylim(cer_ref_logo, cer_alt_logo)

        # ── Motif ─────────────────────────────────────────────────────────────
        if variant_key in cerberus_ism:
            sequence = cerberus_ism[variant_key][0].T
            motif_var_pos = CERBERUS_VAR_POS
        elif cbp_onehot is not None:
            sequence = cbp_onehot
            motif_var_pos = CHROMBPNET_VAR_POS
        else:
            sequence = None

        if sequence is not None:
            pwm, is_rc, _ = align_pwm(
                sequence, motifs[example["motif_id"]],
                motif_var_pos, ISM_LEN, ic_mode="simple",
            )
            print(f"  {variant_key}: {example['motif_name']} reverse={is_rc}")
        else:
            pwm = pd.DataFrame(np.zeros((ISM_LEN, 4)), columns=list("ACGT"))
            is_rc = False

        peak_ic = pwm.values.max()
        pwm_ylim = (0, peak_ic * 1.15 if peak_ic > 0 else 0.01)

        # ── Layout ────────────────────────────────────────────────────────────
        # The alternate-allele ISM row and the PWM sit flush against each other
        # so the motif reads as an annotation of the attribution above it.
        column_grid = gridspec.GridSpecFromSubplotSpec(
            6, 1, subplot_spec=outer[column],
            height_ratios=[COV_H, ISM_H, ISM_H, COV_H, ISM_H, ISM_H + PWM_H],
            hspace=0.10,
        )
        ism_pwm_grid = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=column_grid[5],
            height_ratios=[ISM_H, PWM_H], hspace=0.0,
        )

        # ChromBPNet coverage
        ax = fig.add_subplot(column_grid[0])
        draw_coverage(ax, cbp_x, cbp_ref, cbp_alt, var_pos, cbp_xlim, column, ROW_LABELS[0])
        ax.set_title(
            f"\n{example['label']}\n{CT_LABELS.get(celltype, celltype)}",
            fontsize=6.5, pad=2, fontweight="bold", color="black",
        )
        ax.text(
            0.5, 1.01, "●", transform=ax.transAxes,
            color=CT_COLORS.get(celltype, "#333"), fontsize=7,
            ha="center", va="bottom", clip_on=False, zorder=10,
        )
        if column == 0:
            ax.legend(
                handles=[
                    mpatches.Patch(color=REF_COLOR, alpha=0.7, label="Ref"),
                    mpatches.Patch(color=ALT_COLOR, alpha=0.7, label="Alt"),
                ],
                fontsize=4, frameon=False, loc="upper left", handlelength=0.8,
                labelspacing=0.15, handletextpad=0.3, borderpad=0,
            )

        draw_attribution(
            fig.add_subplot(column_grid[1]), cbp_ref_logo, cbp_ylim,
            column, ROW_LABELS[1], CHROMBPNET_VAR_POS,
        )
        draw_attribution(
            fig.add_subplot(column_grid[2]), cbp_alt_logo, cbp_ylim,
            column, ROW_LABELS[2], CHROMBPNET_VAR_POS,
        )
        draw_coverage(
            fig.add_subplot(column_grid[3]), cer_x, cer_ref, cer_alt,
            var_pos, cer_xlim, column, ROW_LABELS[3],
        )
        draw_attribution(
            fig.add_subplot(column_grid[4]), cer_ref_logo, cer_ylim,
            column, ROW_LABELS[4], CERBERUS_VAR_POS,
        )
        draw_attribution(
            fig.add_subplot(ism_pwm_grid[0]), cer_alt_logo, cer_ylim,
            column, ROW_LABELS[5], CERBERUS_VAR_POS,
        )

        # PWM row
        ax = fig.add_subplot(ism_pwm_grid[1])
        logomaker.Logo(
            pwm, color_scheme=NT_COLORS, ax=ax,
            flip_below=False, baseline_width=0, font_name="DejaVu Sans",
        )
        for patch in ax.patches:
            patch.set_zorder(3)
        ax.axvspan(
            CERBERUS_VAR_POS - 0.5, CERBERUS_VAR_POS + 0.5,
            color=HIGHLIGHT, alpha=0.6, zorder=0, lw=0,
        )
        ax.axvline(CERBERUS_VAR_POS, color="#888", lw=0.5, ls="--", zorder=1)
        ax.set_xlim(-0.5, ISM_LEN - 0.5)
        ax.set_ylim(*pwm_ylim)
        style_row(ax)
        set_row_label(ax, column, ROW_LABELS[6])
        ax.text(
            0.5, -0.40, example["motif_name"] + (" (−)" if is_rc else ""),
            transform=ax.transAxes, ha="center", va="top",
            fontsize=6, style="italic",
        )

    stem = os.path.join(out_dir, "caqtl_examples")
    fig.savefig(f"{stem}.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {stem}.pdf / .png")


if __name__ == "__main__":
    main()
