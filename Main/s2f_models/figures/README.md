# Figures

Plotting code for the sequence-to-function panels. Every script reads score
tables and cached prediction files produced by `../chrombpnet` and `../cerberus`
— **no model is loaded here**, so these run on a laptop in the `s2f`
environment, with neither ChromBPNet nor baskerville installed.

Paths come from `../config.sh`. Each script takes `--out_dir`, defaulting to
`${FIGURE_DIR}`, and writes a matching `.pdf` and `.png`.

## Layout

```
figlib.py        shared style, colours, motif alignment, cell-type maps
caqtl_data.py    the joined caQTL score table
eqtl_data.py     fine-mapped eQTL effects joined to Cerberus logSED

chrombpnet_accuracy.py         held-out accuracy, per cell type + one example fold
cerberus_accuracy.py           held-out accuracy, RNA and ATAC per cell type
caqtl_effect_concordance.py    both models vs the measured RASQUAL effect
caqtl_model_concordance.py     Cerberus vs ChromBPNet
caqtl_examples.py              example caQTLs: coverage, attributions, motif
eqtl_effect_concordance.py     Cerberus logSED vs fine-mapped eQTL effect
eqtl_examples.py               example eQTLs: association, PIP, coverage, motif
chr16_caqtl_locus.py           a PT-specific caQTL, stratified by donor genotype
fgf5_gwas_locus.py             an eGFR GWAS locus, credible set to motif

make_caqtl_scored.py           scored caQTL table
make_eqtl_scored.py            scored eQTL table

examples_caqtl.tsv             which variants caqtl_examples.py draws
supp/                          supplementary panels
```

## Inputs

| Script | Reads |
|--------|-------|
| `chrombpnet_accuracy.py` | `${WORK_DIR}/chrombpnet/metrics/{pearson_r,scatter_data}.csv` |
| `cerberus_accuracy.py` | `${CERBERUS_MODEL_DIR}/f<I>c0/eval0/fold<I>/acc.txt` |
| `caqtl_effect_concordance.py`, `caqtl_model_concordance.py`, `make_caqtl_scored.py` | the caQTL ensemble scores + the merged RASQUAL/ChromBPNet table (`caqtl_data.py`) |
| `caqtl_examples.py` | `${WORK_DIR}/caqtl/examples/{chrombpnet_pred,cerberus_pred,cerberus_ism}`, `${WORK_DIR}/chrombpnet/caqtl/shap/<ct>_ism.h5` |
| `eqtl_effect_concordance.py`, `make_eqtl_scored.py` | the logSED credible-set scores + `${EQTL_SUSIE_TSV}` (`eqtl_data.py`) |
| `eqtl_examples.py` | `${WORK_DIR}/eqtl/examples/{cerberus_pred,cerberus_ism}`, `${EQTL_SUSIE_TSV}`, `${EQTL_GTF}` |
| `chr16_caqtl_locus.py` | `${GENOTYPE_VCF}`, `${FRAG_DIR_INDIVIDUAL}`, ChromBPNet SHAP, Cerberus ISM |
| `fgf5_gwas_locus.py` | `${GWAS_SUSIE_DIR}/<locus>`, `${WORK_DIR}/gwas/*` |

`${GENOTYPE_VCF}`, `${CAQTL_SUSIE_TSV}`, `${EQTL_SUSIE_TSV}` and
`${GWAS_SUSIE_DIR}` are controlled-access; see the top-level README.

## `figlib.py`

The single definition of everything that must agree between panels.

**Configuration.** `load_config()` parses `../config.sh` — the same file the
shell scripts source — expanding `${VAR}` references. Real environment variables
win, so `WORK_DIR=/tmp/x python caqtl_examples.py` works. `cfg("KEY")` raises a
readable error when a key is unset rather than failing later on a missing file.

**Style.** `setup_style()` sets Type 42 font embedding, so text in the PDFs stays
editable rather than being converted to outlines, and loads Lato if present,
falling back to the default font otherwise.

**Motifs.** `parse_pwm` reads a MEME database (cached per file).
`align_pwm(onehot, pwm, var_pos, n_bases)` places a motif by log-odds against
the sequence; passing a `var_pos` restricts placements to those covering the
variant, since a motif that does not overlap the variant cannot explain its
effect. `align_pwm_to_cwm` instead matches the attribution profile, TomTom-style,
for the panels where the motif is discovered rather than specified.

`ic_mode` picks the letter-height convention: `"column"` scales each column by
its total information content (the usual sequence-logo convention), `"simple"`
scales each letter by its own, `P(i,b) · log₂[P(i,b)/0.25]`. Both appear in the
published panels.

**Attributions.** `ism_to_logo` converts raw ISM scores to what is plotted: for
each position, the score of the base actually present minus the mean of the
three alternatives.

**Cell types.** `CT_LABELS`, `CT_COLORS`, and the `PT` ↔ `PTS` conversion.

## Cell-type naming

Proximal tubule is spelled **`PT`** by RASQUAL and the ChromBPNet model
directories, **`PTS`** by the Cerberus targets and the eQTL tables, and is always
displayed as **PTS**. `figlib.model_ct()` and `figlib.qtl_ct()` convert between
them; don't hardcode either spelling.

## Track indexing

Three orderings are in play, and using the wrong one silently produces plausible
numbers for the wrong cell type. All three are defined in `figlib.py`:

| Map | Tracks | Where it applies |
|-----|--------|------------------|
| `ATAC_IDX` / `RNA_IDX` | 124 | the full human head after RNA strand pairs are averaged — what `hound_snp` emits for `kidney_targets_w5_human{,_local}.txt`. ATAC is the even indices 0–22, the matching RNA+ the odd 1–23. |
| `ATAC_SUBSET_IDX` | 11 | the ATAC-only targets file used for example-locus ISM, to keep those outputs small |
| (stride 3) | 36 | the raw Susztak targets file before strand collapsing; no current script uses it |

## Allele orientation

Model scores are always **log(ALT/REF)** in hg38 orientation. The eQTL tables
record A1 as the effect allele, which is not always the reference; where
`match_type == "swap"` the slope must be negated before comparison.
`eqtl_data.load_eqtl` does this, and the scored table reports both `eqtl_beta`
and `eqtl_beta_aligned`. Comparing against the raw beta flips the sign for those
variants and washes out the correlation.

## Choosing the example loci

The example variants in `examples_caqtl.tsv` were selected, not picked at
random. caQTL examples needed a large predicted effect in the target cell type,
a specificity ratio of at least 2 against the mean across other cell types,
agreement in direction between the two models, and a motif spanning the variant.
eQTL examples needed PIP ≥ 0.5, agreement in sign between logSED and the eQTL
slope, and a clear motif; the two shown run in opposite directions.

Edit `examples_caqtl.tsv` to draw different variants — you will need their
predictions and attributions cached first.

## Scored tables

`make_caqtl_scored.py` and `make_eqtl_scored.py` write per-variant CSVs carrying
the measured effect alongside both models' predictions, so the concordance
analyses can be reproduced without rerunning a model. Columns are documented in
each script's docstring.

`make_eqtl_scored.py`'s `dist_to_tss` is the distance to the gene's 5' end from
the GTF `gene` records. An earlier hand-built version of this table used a
different, undocumented distance; every other column reproduces exactly.

## Software

Python 3.11 with numpy, pandas, scipy, matplotlib, logomaker, h5py, pysam,
pyBigWig. Lato is optional. See `../environment.yml`.
