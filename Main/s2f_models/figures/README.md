# Figures

Plotting code for the sequence-to-function panels. Every script reads score
tables and cached predictions produced by `../chrombpnet` and `../cerberus` —
**no model is loaded here**, so these run without a GPU.

Paths come from `../config.sh`. Each script takes `--out_dir`, defaulting to
`${FIGURE_DIR}`, and writes a matching `.pdf` and `.png`.

## What each script draws

| Script | Panel |
|--------|-------|
| `chrombpnet_accuracy.py` | ChromBPNet held-out accuracy per cell type, plus one example fold as a density scatter |
| `cerberus_accuracy.py` | Cerberus held-out accuracy per cell type, RNA and ATAC |
| `caqtl_effect_concordance.py` | both models' predicted caQTL effect against the measured one |
| `caqtl_model_concordance.py` | the two models against each other, coloured by measured effect |
| `caqtl_examples.py` | example caQTLs: predicted coverage, attributions on both alleles, disrupted motif |
| `eqtl_effect_concordance.py` | predicted expression effect against fine-mapped eQTL effect size |
| `eqtl_examples.py` | example eQTLs: association, fine-mapping, predicted coverage, motif |
| `chr16_caqtl_locus.py` | one caQTL across cell types, with coverage stratified by donor genotype |
| `fgf5_gwas_locus.py` | a GWAS locus end to end: credible set, model scores, attributions, motif |
| `make_caqtl_scored.py` | per-variant caQTL table: measured effect plus both models' predictions |
| `make_eqtl_scored.py` | the same for fine-mapped eQTLs |
| `supp/` | supplementary variants of the above, plus specificity and effect-vs-PIP panels |

`caqtl_data.py` and `eqtl_data.py` hold the loading logic shared between a
figure and its table, so the two cannot drift apart in which scores they use.

## `figlib.py`

The single definition of everything that must agree between panels.

- **Configuration.** `load_config()` parses `../config.sh` — the same file the
  shell scripts source — expanding `${VAR}` references. Real environment
  variables win, so `WORK_DIR=/tmp/x python caqtl_examples.py` works. `cfg()`
  raises a readable error on an unset key rather than failing later on a
  missing file.
- **Style.** `setup_style()` sets Type 42 font embedding so PDF text stays
  editable rather than becoming outlines, and loads Lato if present, falling
  back to the default font otherwise.
- **Motifs.** `parse_pwm` reads a MEME database. `align_pwm` places a motif by
  log-odds against the sequence; given a variant position it only considers
  placements covering it, since a motif that misses the variant cannot explain
  its effect. `align_pwm_to_cwm` instead matches the attribution profile, for
  the panels where the motif is discovered rather than named.
- **Attributions.** `ism_to_logo` turns raw ISM scores into what gets plotted:
  per position, the base actually present minus the mean of the alternatives.
- **Cell types.** Labels, colours, and the `PT` ↔ `PTS` conversion.

## Things that will bite you

**Cell-type naming.** Proximal tubule is `PT` to RASQUAL and the ChromBPNet
model directories, `PTS` to the Cerberus targets and the eQTL tables, and always
displays as PTS. Use `figlib.model_ct()` / `qtl_ct()`; don't hardcode either.

**Track indexing.** Three orderings coexist, and the wrong one silently gives
plausible numbers for the wrong cell type:

| Map | Tracks | Applies to |
|-----|--------|------------|
| `ATAC_IDX` / `RNA_IDX` | 124 | the full human head after RNA strand pairs are averaged — what `hound_snp` emits for the standard targets file. ATAC is even indices, the matching RNA+ odd. |
| `ATAC_SUBSET_IDX` | 11 | the ATAC-only targets file used for example-locus ISM |
| (stride 3) | 36 | the raw targets file before strand collapsing; unused by current scripts |

**Allele orientation.** Model scores are log(ALT/REF) in hg38 orientation. The
eQTL tables record A1 as the effect allele, which is not always the reference;
where `match_type` is `swap` the slope must be negated. `eqtl_data.load_eqtl`
does this, and the scored table reports both the raw and aligned values.
Comparing against the raw beta flips the sign for those variants.

**Example loci are configured, not hardcoded.** `examples_caqtl.tsv` lists which
variants `caqtl_examples.py` draws. Swapping them requires their predictions and
attributions to be cached first.

**`make_eqtl_scored.py`'s `dist_to_tss`** is the distance to the gene's 5' end
from the GTF `gene` records. An earlier hand-built version of this table used a
different, undocumented distance; every other column reproduces exactly.

## Software

Python with numpy, pandas, scipy, matplotlib, logomaker, h5py, pysam and
pyBigWig — all of which the baskerville environment already provides. Lato is
optional.
