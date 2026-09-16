# Sequence-to-function models

Two deep-learning models that predict chromatin accessibility and gene
expression directly from DNA sequence, trained on the kidney snMultiome atlas
and used to score and interpret molecular QTL variants. Model architectures,
training settings and analysis parameters are described in the paper methods;
this README covers what the code does and how to run it.

- **ChromBPNet** — short-context, one model per cell type, accessibility only.
- **Cerberus** — long-context, a single model covering all cell types and both
  species, predicting accessibility *and* expression.

They are complementary. ChromBPNet sees only local sequence, so any effect it
reports is necessarily local. Cerberus sees the whole locus and can score a
variant's effect on a distant gene — but to compare the two on accessibility,
Cerberus is restricted to a matched local window.

## Folders

| Folder | What it does |
|--------|--------------|
| `data/` | turns aligned multiome reads into per-cell-type coverage, and builds the targets table listing the tracks Cerberus predicts |
| `chrombpnet/` | calls peaks, trains the per-cell-type models, scores variants, and computes per-variant and genome-wide attributions |
| `cerberus/` | builds the training examples, fine-tunes and evaluates the folds, scores variants for accessibility and expression, and runs in-silico mutagenesis |
| `figures/` | draws every panel and writes the scored-variant tables, from cached scores only |

Each folder has its own README with a script-by-script table.

## Getting started

```bash
cp config.example.sh config.sh    # then edit the paths in it
```

`config.sh` is the single place paths live. Shell scripts `source` it; Python
scripts parse the same file through `figlib.load_config()`. Nothing in this
repository is a data file.

**The released dataset starts at per-cell-type fragments and `.w5` coverage.**
`data/` documents how those were produced, but you do not need to rerun it —
point `config.sh` at the downloaded files and start at `chrombpnet/` or
`cerberus/`.

**Redrawing the figures needs neither model.** The plotting scripts read cached
score tables and predictions, so `figures/` runs without a GPU and without
either model framework.

## Environments

Two, because the frameworks conflict:

| Environment | For |
|---|---|
| `baskerville` | Cerberus training, scoring and ISM — and the figures, whose dependencies it already provides |
| `chrombpnet` | ChromBPNet training, attributions, and variant-scorer |

Scripts are invoked through console entry points (`hound_*`, `chrombpnet`), so
no repository paths are hardcoded.

Pin baskerville to a revision providing `--local_window` and the `window`
targets column. Without it, `hound_snp --local_window` is **silently ignored**
and caQTL scores come out whole-sequence rather than local — plausible numbers
that are not comparable to the published ones or to ChromBPNet.

## How the stages chain

```
   released fragments + .w5 coverage
        │
        ├── chrombpnet/  peaks ──> fold models ──> variant scores, attributions
        │                                                  │
        └── cerberus/    training data ──> fold models ──> scores, ISM
                                                           │
                                          figures/  panels + scored tables
```

Both models score **the same caQTL variant set** — the VCF is built once by
`cerberus/score/caqtl/1_prepare_variants.py` and reused — so the two score
columns are directly comparable. All scores are alternate-minus-reference in
hg38 orientation and averaged across folds.

## Data access

Inputs derived from individual-level genotypes are controlled-access. They are
referenced through `config.sh` and never distributed here: donor genotypes,
RASQUAL caQTL results, and the caQTL, eQTL and GWAS fine-mapping outputs.

The public BioProjects contributing tracks to the model are listed in
`data/README.md`.

## Citation

ChromBPNet, TF-MoDISco and variant-scorer are from the Kundaje lab; baskerville
and the Cerberus architecture are from Calico.
