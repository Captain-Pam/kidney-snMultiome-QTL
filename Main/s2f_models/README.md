# Sequence-to-function models

Two deep-learning models that predict chromatin accessibility and gene
expression directly from DNA sequence, trained on the kidney snMultiome atlas
and used to score and interpret molecular QTL variants.

| | **ChromBPNet** | **Cerberus** |
|---|---|---|
| Context | 2,114 bp | 786,432 bp |
| Output | 1 kb base-resolution ATAC profile + total counts | 32 bp bins, 182 human + 373 mouse tracks |
| Scope | one model per cell type | one model, all cell types and both species |
| Predicts | accessibility | accessibility **and** expression |
| Folds | 5 | 8 |
| Framework | TensorFlow | PyTorch (baskerville) |

They are complementary. ChromBPNet sees only the local sequence, so a variant
effect it reports is necessarily local. Cerberus sees the whole locus, so it can
score a variant's effect on a gene hundreds of kilobases away — but to compare
the two on accessibility, Cerberus is restricted to a matched local window.

## Layout

```
data/         pseudobulk coverage and the model targets files
chrombpnet/   peak calling, training, variant scoring, attributions
cerberus/     training data, fine-tuning, evaluation, scoring, ISM
figures/      all panels and the scored-variant tables
```

Each directory has its own README with an input→output table per script.

## Getting started

```bash
cp config.example.sh config.sh    # then edit every path in it
conda env create -f environment.yml
```

`config.sh` is the single place paths live. Shell scripts `source` it; Python
scripts parse the same file through `figlib.load_config()`. Nothing in this
repository is a data file.

**The released dataset starts at per-cell-type fragments and `.w5` coverage.**
`data/` documents how those were produced from raw Cell Ranger output, but you
do not need to rerun it — point `config.sh` at the downloaded files and start at
`chrombpnet/` or `cerberus/`.

**To redraw the figures you need neither model.** The plotting scripts read
cached score tables and predictions, so `figures/` runs in the `s2f` environment
alone. Only the scoring and training stages need the model frameworks.

## Environments

Three, because the two frameworks conflict:

| Environment | For |
|---|---|
| `s2f` (`environment.yml`) | everything in `figures/`, and the ensembling scripts |
| `baskerville` | Cerberus training, `hound_snp`, `hound_ism_snp`, `hound_eval` |
| `chrombpnet` | ChromBPNet training, contributions, `variant-scorer` |

Scripts are invoked through console entry points (`hound_*`, `chrombpnet`), so
no repository paths are hardcoded.

> Two different packages are called **baskerville**. The PyTorch one
> (`baskerville`, formerly `baskerville-torch`) trains and scores Cerberus. The
> TensorFlow one (`calico/baskerville`) supplies `bw_w5.py`, used once in
> `data/` to convert bigWigs. Both are needed; they are not interchangeable.

Pin the PyTorch baskerville to a revision providing `--local_window` and the
`window` targets column. Without it, `hound_snp --local_window` is **silently
ignored** and the caQTL scores are whole-sequence rather than local — plausible
numbers that are not comparable to the published ones or to ChromBPNet.

## Pipeline

```
   released fragments + .w5 coverage
        │
        ├── chrombpnet/  peaks ──> 5 fold models ──> variant scores, SHAP
        │                                                  │
        └── cerberus/    hound_data ──> 8 fold models ──> logSUM / logSED, ISM
                                                           │
                                          figures/  panels + scored tables
```

Both models score **the same caQTL variant set**: the VCF is built once by
`cerberus/score/caqtl/1_prepare_variants.py` and reused, so the two score
columns are directly comparable.

## Statistics

| Statistic | Model | Meaning |
|---|---|---|
| `logfc` | ChromBPNet | log fold-change in predicted total counts between alleles |
| `logSUM` | Cerberus | log ratio of predicted coverage summed over ±512 bp |
| `covgene/logSED` | Cerberus | log ratio of predicted RNA summed over a gene's exons |

All are alternate-minus-reference in hg38 orientation, and all are averaged
across folds. See `figures/README.md` on aligning eQTL effect sizes to that
orientation.

## Data access

Model inputs derived from individual-level genotypes are controlled-access and
are referenced by `config.sh`, never distributed here:

- donor genotypes (`GENOTYPE_VCF`) — the genotype-stratified coverage panel
- RASQUAL caQTL results (`RASQUAL_CAQTL_DIR`)
- caQTL, eQTL and GWAS fine-mapping (`CAQTL_SUSIE_TSV`, `EQTL_SUSIE_TSV`, `GWAS_SUSIE_DIR`)

The public BioProjects contributing tracks to the model are listed in
`data/README.md`.

## Citation

ChromBPNet, TF-MoDISco and variant-scorer are from the Kundaje lab; baskerville
and the Cerberus architecture are from Calico. See each directory's README for
the tools it uses.
