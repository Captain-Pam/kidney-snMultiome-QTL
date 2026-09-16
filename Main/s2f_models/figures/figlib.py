"""Shared plotting, motif and indexing helpers for the figure scripts.

Everything here was previously copy-pasted across the individual figure scripts.
Importing it keeps one definition of the things that must agree between panels:
the colour scheme, the track-index conventions, and how a PWM is aligned to an
attribution window.

Configuration comes from `../config.sh`, which this module parses so that both
the shell scripts and the Python scripts read one file. Real environment
variables win over the file, so a path can be overridden for a single run.
"""
import os
import re

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

# ─── Configuration ────────────────────────────────────────────────────────────

_CONFIG_CACHE = None


def load_config(path=None):
    """Parse `config.sh` into a dict, expanding ${VAR} references.

    Real environment variables take precedence, so
    `WORK_DIR=/tmp/x python caqtl_examples.py` works as expected.
    """
    global _CONFIG_CACHE
    if _CONFIG_CACHE is not None and path is None:
        return _CONFIG_CACHE

    if path is None:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(os.path.dirname(here), "config.sh")

    values = {}
    if os.path.exists(path):
        pattern = re.compile(r'^\s*export\s+([A-Za-z_][A-Za-z0-9_]*)=(.*)$')
        with open(path) as handle:
            for line in handle:
                match = pattern.match(line)
                if not match:
                    continue
                key, raw = match.group(1), match.group(2).strip()
                # Strip a trailing comment outside quotes, then the quotes.
                if raw and raw[0] not in "\"'":
                    raw = raw.split("#")[0].strip()
                raw = raw.strip('"').strip("'")
                # Expand ${OTHER} against what we have so far, then the real env.
                raw = re.sub(
                    r"\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?",
                    lambda m: values.get(m.group(1), os.environ.get(m.group(1), "")),
                    raw,
                )
                values[key] = raw

    _CONFIG_CACHE = values
    return values


def cfg(key, default=None):
    """One configuration value, with a clear error if it is missing."""
    value = os.environ.get(key) or load_config().get(key)
    if value in (None, ""):
        if default is not None:
            return default
        raise KeyError(
            f"{key} is not set. Copy config.example.sh to config.sh and fill it in."
        )
    return value


# ─── Style ────────────────────────────────────────────────────────────────────

# Type 42 embeds fonts as editable text rather than outlines, so the PDFs can be
# relabelled in a vector editor.
PDF_RCPARAMS = {"pdf.fonttype": 42, "ps.fonttype": 42}

FONT_DIR = "/usr/share/fonts/truetype/lato"
FONT_FACES = [
    "Lato-Regular.ttf",
    "Lato-Bold.ttf",
    "Lato-Italic.ttf",
    "Lato-BoldItalic.ttf",
]


def setup_style(font_dir=None):
    """Apply the figure style. Falls back to the default font if Lato is absent."""
    matplotlib.rcParams.update(PDF_RCPARAMS)

    font_dir = font_dir or os.environ.get("FIGURE_FONT_DIR", FONT_DIR)
    loaded = False
    for face in FONT_FACES:
        candidate = os.path.join(font_dir, face)
        if os.path.exists(candidate):
            fm.fontManager.addfont(candidate)
            loaded = True
    if loaded:
        matplotlib.rcParams["font.family"] = "Lato"


# Reference and alternate allele, used for every predicted-coverage panel.
REF_COLOR = "#d7191c"
ALT_COLOR = "#2c7bb6"
HIGHLIGHT = "#d0e8ff"

# Nucleotide colours for attribution logos and PWMs.
NT_COLORS = {"A": "#109648", "C": "#255C99", "G": "#F7B32B", "T": "#D62839"}


def style_ax(ax, box_aspect=1):
    """Shared scatter-panel styling: no top/right spines, thin ticks, few labels."""
    ax.spines[["top", "right"]].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(0.4)
    ax.tick_params(labelsize=6, width=0.4, length=2, pad=1.5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=2, prune="both"))
    ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="both"))
    if box_aspect is not None:
        ax.set_box_aspect(box_aspect)


def add_ref_lines(ax, x0, y0):
    ax.axhline(y0, color="lightgrey", lw=0.4, zorder=1)
    ax.axvline(x0, color="lightgrey", lw=0.4, zorder=1)


def add_ols(ax, x, y, color="grey"):
    """Ordinary least-squares fit line across the observed x range."""
    slope, intercept = np.polyfit(x, y, 1)
    ax.plot(
        [x.min(), x.max()],
        [slope * x.min() + intercept, slope * x.max() + intercept],
        ls=":", lw=0.8, color=color, alpha=0.9, zorder=3,
    )


# ─── Cell types ───────────────────────────────────────────────────────────────

# The 12 cell types, in the order used throughout.
CELL_TYPES = [
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PT", "Podocyte", "Stromal", "TAL", "injPT",
]

# Display labels. Note PT vs PTS: the RASQUAL/ChromBPNet side calls proximal
# tubule "PT", the Cerberus targets and eQTL tables call it "PTS", and the
# figures always show "PTS". Use `model_ct()` / `qtl_ct()` to move between them
# rather than hardcoding either spelling.
CT_LABELS = {
    "CNT_CD_PC": "CNT/CD/PC", "DCT": "DCT", "DTL_ATL": "DTL/ATL",
    "EC": "EC", "IC": "IC", "Immune": "Immune", "PEC": "PEC",
    "PT": "PTS", "Podocyte": "Podocyte", "Stromal": "Stromal",
    "TAL": "TAL", "injPT": "injPT",
}

CT_COLORS = {
    "CNT_CD_PC": "#b5bd61", "DCT": "#006fa6", "DTL_ATL": "#ffb3c6",
    "EC": "#d62728", "IC": "#aa40fc", "Immune": "#279e68", "PEC": "#6a3a4c",
    "PT": "#e377c2", "Podocyte": "#ff7f0e", "injPT": "#950046",
    "Stromal": "#8c6d31", "TAL": "#aec7e8",
}

_QTL_TO_MODEL = {"PT": "PTS"}
_MODEL_TO_QTL = {"PTS": "PT"}


def model_ct(celltype):
    """QTL-table cell-type name -> model/targets name (PT -> PTS)."""
    return _QTL_TO_MODEL.get(celltype, celltype)


def qtl_ct(celltype):
    """Model/targets cell-type name -> QTL-table name (PTS -> PT)."""
    return _MODEL_TO_QTL.get(celltype, celltype)


# ─── Track indexing ───────────────────────────────────────────────────────────
#
# Three different track orderings are in play. Picking the wrong one yields
# plausible numbers for the wrong cell type, so always index through these maps.
#
# 1. STRAND_COLLAPSED (124 tracks) — the full human head after RNA strand pairs
#    are averaged. This is what hound_snp and hound_ism_snp emit when given
#    kidney_targets_w5_human{,_local}.txt. ATAC tracks are the even indices
#    0-22, the matching RNA+ tracks the odd indices 1-23.
#
# 2. ATAC_SUBSET (11 tracks) — the small ATAC-only targets file used for the
#    example-locus ISM runs, to keep those outputs small.
#
# 3. The raw 36-row Susztak targets file, stride 3 per cell type
#    (ATAC, RNA+, RNA-). Only relevant if you regenerate targets without the
#    strand collapsing; no current script uses it.

ATAC_IDX = {
    "CNT_CD_PC": 0, "DCT": 2, "DTL_ATL": 4, "EC": 6, "IC": 8, "Immune": 10,
    "PEC": 12, "PT": 14, "Podocyte": 16, "Stromal": 18, "TAL": 20, "injPT": 22,
}

RNA_IDX = {ct: idx + 1 for ct, idx in ATAC_IDX.items()}

ATAC_SUBSET_IDX = {
    "CNT_CD_PC": 0, "EC": 3, "IC": 4, "TAL": 8, "injPT": 9, "PT": 10,
}


# ─── Motifs ───────────────────────────────────────────────────────────────────

_MOTIF_CACHE = {}


def parse_pwm(motif_id, motif_db=None):
    """Position probability matrix for one motif from a MEME-format database."""
    motif_db = motif_db or cfg("MOTIF_DB")
    if motif_db not in _MOTIF_CACHE:
        with open(motif_db) as handle:
            _MOTIF_CACHE[motif_db] = handle.read()
    text = _MOTIF_CACHE[motif_db]

    match = re.search(
        rf"MOTIF {re.escape(motif_id)}.*?\n(.*?)(?=MOTIF|\Z)", text, re.DOTALL
    )
    if match is None:
        raise KeyError(f"{motif_id} not found in {motif_db}")

    rows = [
        line.split()
        for line in match.group(0).splitlines()
        if re.match(r"^\s*[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+", line)
    ]
    return np.array([[float(value) for value in row] for row in rows])


def rc_matrix(matrix):
    """Reverse complement of a position matrix."""
    return matrix[::-1, [3, 2, 1, 0]]


def pwm_to_ic_df(pwm):
    """Scale each position of a PPM by that position's total information content.

    The usual sequence-logo convention: letter heights within a column are
    proportional to their probability, and the column's total height is its
    information content in bits.
    """
    eps = 1e-12
    ic = np.log2(pwm + eps) - np.log2(0.25)
    ic[ic < 0] = 0
    total = (pwm * ic).sum(axis=1, keepdims=True)
    return pd.DataFrame(pwm * total, columns=list("ACGT"))


def pwm_to_ic_simple(pwm, background=0.25):
    """Per-letter information content, IC(i, b) = P(i, b) * log2[P(i, b) / bg].

    Differs from `pwm_to_ic_df` in that each letter is scaled by its own
    contribution rather than by the column total, so letters below background
    frequency vanish instead of being drawn small. Slightly different heights,
    same motif.
    """
    return np.clip(pwm * np.log2(np.clip(pwm, 1e-9, 1) / background), 0, None)


def align_pwm(ref_onehot, prob_pwm, var_pos, n_bases, ic_mode="column"):
    """Place a PWM where it best matches the sequence.

    Scores every offset and both orientations by log-odds against a uniform
    background. When `var_pos` is given, only placements covering it are
    considered — a motif that does not overlap the variant cannot explain its
    effect. Pass `var_pos=None` to search the whole window unconstrained.

    `ic_mode` selects the letter-height convention: "column" scales by the
    column's total information content (`pwm_to_ic_df`), "simple" by each
    letter's own (`pwm_to_ic_simple`).

    Returns (IC dataframe padded to n_bases, is_reverse_complement, offset).
    """
    length = len(prob_pwm)
    background = 0.25
    lo_fwd = np.log2(np.clip(prob_pwm, 1e-9, 1) / background)
    lo_rev = np.log2(np.clip(rc_matrix(prob_pwm), 1e-9, 1) / background)
    seq_idx = ref_onehot.argmax(axis=-1)

    best, best_pos, best_rc = -np.inf, 0, False
    for pos in range(n_bases - length + 1):
        if var_pos is not None and not (pos <= var_pos <= pos + length - 1):
            continue
        window = seq_idx[pos:pos + length]
        fwd = sum(lo_fwd[j, window[j]] for j in range(length))
        rev = sum(lo_rev[j, window[j]] for j in range(length))
        if fwd > best:
            best, best_pos, best_rc = fwd, pos, False
        if rev > best:
            best, best_pos, best_rc = rev, pos, True

    use = rc_matrix(prob_pwm) if best_rc else prob_pwm
    aligned = np.zeros((n_bases, 4))

    if ic_mode == "simple":
        aligned[best_pos:best_pos + length] = pwm_to_ic_simple(use)
        return pd.DataFrame(aligned, columns=list("ACGT")), best_rc, best_pos

    aligned[best_pos:best_pos + length] = np.clip(use, 1e-6, 1)
    frame = pwm_to_ic_df(np.clip(aligned, 1e-6, 1))
    frame.iloc[:best_pos] = 0.0
    frame.iloc[best_pos + length:] = 0.0
    return frame, best_rc, best_pos


def align_pwm_to_cwm(cwm, prob_pwm, n_bases):
    """Place a PWM where its IC profile best matches an attribution matrix.

    TomTom-style: rather than matching the sequence, match the contribution
    scores, so the motif lands where the model actually attributes signal. Used
    when the motif is discovered rather than specified.

    Returns (IC dataframe, is_reverse_complement, offset).
    """
    length = len(prob_pwm)

    def ic_weighted(pwm):
        ic = np.log2(np.clip(pwm, 1e-9, 1)) - np.log2(0.25)
        ic[ic < 0] = 0
        return pwm * ic

    candidates = [(False, ic_weighted(prob_pwm)), (True, ic_weighted(rc_matrix(prob_pwm)))]

    best, best_pos, best_rc = -np.inf, 0, False
    for is_rc, weighted in candidates:
        for pos in range(n_bases - length + 1):
            score = float(np.sum(weighted * cwm[pos:pos + length]))
            if score > best:
                best, best_pos, best_rc = score, pos, is_rc

    use = rc_matrix(prob_pwm) if best_rc else prob_pwm
    aligned = np.zeros((n_bases, 4))
    aligned[best_pos:best_pos + length] = np.clip(use, 1e-6, 1)
    frame = pwm_to_ic_df(np.clip(aligned, 1e-6, 1))
    frame.iloc[:best_pos] = 0.0
    frame.iloc[best_pos + length:] = 0.0
    return frame, best_rc, best_pos


# ─── Attributions ─────────────────────────────────────────────────────────────


def ism_to_logo(ism_scores, onehot):
    """Convert raw ISM scores to the per-position contribution that gets plotted.

    ISM gives the score of every possible substitution. The contribution of the
    base actually present is defined as its score minus the mean of the three
    alternatives; only that base is drawn, so the logo shows what the present
    base contributes relative to changing it.

    `ism_scores` is (L, 4), `onehot` is (L, 4). Returns an (L, 4) dataframe with
    one non-zero entry per row.
    """
    length = ism_scores.shape[0]
    contrib = np.zeros((length, 4))
    for i in range(length):
        present = int(np.argmax(onehot[i]))
        others = [j for j in range(4) if j != present]
        contrib[i, present] = -float(np.mean(ism_scores[i, others]))
    return pd.DataFrame(contrib, columns=list("ACGT"))


def draw_logo(ax, logo_df, ylim, var_pos, n_bases, ylabel=None, labelsize=5.5):
    """Draw an attribution or PWM logo with the variant position marked."""
    import logomaker

    logomaker.Logo(
        logo_df, color_scheme=NT_COLORS, ax=ax,
        baseline_width=0, font_name="DejaVu Sans",
    )
    ax.axvline(var_pos, color="#888", lw=0.5, ls="--", zorder=1)
    ax.set_xlim(-0.5, n_bases - 0.5)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xticks([])
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=labelsize)
    else:
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)


# ─── Annotation ───────────────────────────────────────────────────────────────

_GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')


def parse_gtf_gene(gtf_path, gene_name):
    """Span, exons and strand of one gene.

    Returns (start, end, exons, strand) with 1-based inclusive coordinates.
    Exons are pooled across transcripts, which is what the gene-model track and
    the exonic coverage sums both want.
    """
    exons = []
    start = end = None
    strand = "+"

    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon", "gene"):
                continue
            match = _GENE_NAME_RE.search(fields[8])
            if not match or match.group(1) != gene_name:
                continue

            feature_start, feature_end = int(fields[3]), int(fields[4])
            strand = fields[6]
            if fields[2] in ("transcript", "gene"):
                start = feature_start if start is None else min(start, feature_start)
                end = feature_end if end is None else max(end, feature_end)
            else:
                exons.append((feature_start, feature_end))

    if start is None:
        raise KeyError(f"{gene_name} not found in {gtf_path}")
    return start, end, sorted(set(exons)), strand


# ─── QTL helpers ──────────────────────────────────────────────────────────────


def select_lead(df, score_col, group_col="Feature_ID"):
    """One variant per peak: the one with the largest absolute model score.

    A peak typically has many caQTL variants in tight LD with indistinguishable
    measured effects. Taking the model's own top-scoring variant per peak gives
    one independent point per peak.
    """
    valid = df[df[score_col].notna()].copy()
    idx = valid.groupby(group_col)[score_col].apply(lambda s: s.abs().idxmax())
    return valid.loc[idx].reset_index(drop=True)
