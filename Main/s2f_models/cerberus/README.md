# Fine-tuned Cerberus model

Cerberus is a long-range sequence-to-function model: a large window of DNA in,
binned coverage tracks out, through a convolutional trunk and Hydra
(bi-directional Mamba) blocks with one output head per species. We initialised
from a pretrained Cerberus model and fine-tuned on human and mouse kidney
coverage. Architecture and training hyperparameters are in `params.json` and in
the paper methods.

Everything here uses **baskerville** (the PyTorch port of Calico's baskerville)
through its `hound_*` console entry points. Paths come from `../config.sh`.

## Training and evaluation

| Script | Does |
|--------|------|
| `1_hound_data.sh {hg38\|mm10}` | reads each target's `.w5` coverage over the split intervals and writes the per-fold training examples |
| `2_train_folds.sh` | fine-tunes one model per fold, each from the matching fold of the pretrained checkpoint |
| `3_eval_folds.sh` | per-track accuracy of every model fold against every data fold, for both genomes |
| `params.json` | the as-run architecture and optimiser settings |

Two things that are easy to get wrong:

- **The cross-fold splits are not generated here.** `sequences.bed`,
  `contigs.bed` and friends come with the pretrained dataset and must be in
  place before step 1. Reusing them is what guarantees that a fold held out here
  was also held out during pretraining; generating fresh splits would leak test
  sequences in through the initialisation.
- **Fold indices must stay aligned.** Model `f<I>c0` is initialised from
  pretrained fold `I`. Initialising from a checkpoint that saw fold `I` during
  pretraining defeats the held-out split.

The two data directories become **dataset 0 (human)** and **dataset 1 (mouse)**;
that numbering is the `--head` / `--dataset` argument every downstream script
uses. The accuracy figures read only the diagonal of the evaluation grid — each
model on its own held-out split.

## Provenance of these scripts

`params.json` is the as-run configuration, copied from the trained model.
`3_eval_folds.sh` is a cleaned-up copy of a wrapper that survived.

**`1_hound_data.sh` and `2_train_folds.sh` are reconstructions.** The original
invocations were issued interactively and no job script for them survives. They
were rebuilt from `params.json`, the `statistics.json` the data step wrote, and
templates left by earlier training generations. Two `hound_data` arguments are
inferred rather than recorded — the write-job batch size and worker count —
neither of which changes the contents of the training data, only how it is
produced. Everything that does affect the data comes from the preserved configs.

## `score/` — variant effects

| Stage | Scripts | Statistic |
|-------|---------|-----------|
| caQTL | `caqtl/1_prepare_variants.py`, `2_score_folds.sh` | `logSUM` over a local window |
| eQTL | `eqtl/1_prepare_gtf.py`, `2_prepare_vcf.py`, `3_score_folds.sh` | `covgene/logSED` over the eGene's exons |
| GWAS | `gwas/1_score_credible_set.py`, `2_predict_tracks.py` | both, for one named gene and cell type |

`ensemble_folds.py` then averages any of those outputs across folds — it works
on variant scores, gene scores and ISM alike, so it replaces what used to be
three near-identical scripts. `cerberus_model.py` loads a fold and predicts a
reference/alternate pair, for the places that need the coverage profile itself
rather than a summary score.

**`--local_window` needs the `_local` targets file.** The caQTL score is
restricted to a window around the variant, which is what makes it comparable to
ChromBPNet's much smaller receptive field. That restriction is triggered by the
`window` column in the `_local` targets file; passing `--local_window` with the
plain targets file is **silently ignored** and yields a whole-sequence score
instead. The eQTL scoring deliberately uses the plain file, because there the
gene's exon set defines the region.

**Which variants get scored.** `caqtl/1_prepare_variants.py` keeps only variants
whose scoring window fully contains their peak — otherwise the score measures
something different from the peak-level effect it is compared against.
`eqtl/1_prepare_gtf.py` filters the annotation to eQTL target genes, so each
variant is paired with its eGene rather than with every overlapping gene.

**Allele orientation.** Scores are always log(ALT/REF) in hg38 orientation. The
eQTL table's effect allele is not always the reference; where `match_type` is
`swap` the eQTL slope must be negated before comparison.

## `interpret/` — in-silico mutagenesis

```
bash interpret/run_ism.sh {atac|rna} <variants.vcf> <out_dir>
```

Substitutes every base in a window with all three alternatives and re-scores,
using the accessibility statistic or the expression one. With no fold argument
it submits all folds; average them with `score/ensemble_folds.py`.

The contribution plotted in the figures is derived from these raw scores: the
score of the base actually present minus the mean of the three alternatives.

## Track indexing

Three track orderings are in play and mixing them up produces plausible numbers
for the wrong cell type. They are defined in one place,
`../figures/figlib.py`, and documented in `../figures/README.md`.
`cerberus_model.CerberusFold.check_track` asserts that a track's description
matches what the caller expects, so an index error fails loudly.

## Software

baskerville (PyTorch) with the `hound_*` entry points, PyTorch, pysam.
