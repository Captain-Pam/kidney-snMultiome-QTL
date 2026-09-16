# Fine-tuned Cerberus model

Cerberus is a long-range sequence-to-function model: **786,432 bp** of input
sequence in, binned coverage tracks at **32 bp** resolution out. The trunk is a
DNA convolution stem, a convolutional tower, and Hydra (bi-directional Mamba)
blocks, followed by one output head per species. We initialised from a
pretrained Cerberus model and fine-tuned the whole network on human and mouse
kidney coverage.

Everything here uses **baskerville** (the PyTorch port of Calico's baskerville)
through its `hound_*` console entry points. Paths come from `../config.sh`.

> The package providing `bw_w5.py` in `../data` is the *TensorFlow* baskerville,
> a different codebase that happens to share the name. Both are needed.

## Layout

```
params.json        as-run architecture and optimiser settings
1_hound_data.sh    targets + splits -> per-fold training examples
2_train_folds.sh   fine-tune 8 folds from the pretrained checkpoints
3_eval_folds.sh    per-track accuracy on held-out sequences
score/             variant effect scoring
  ensemble_folds.py    average any hound output across folds
  cerberus_model.py    load a fold, predict a ref/alt pair
  caqtl/  eqtl/  gwas/
interpret/         in-silico mutagenesis
```

## Model

| | |
|---|---|
| Input | 786,432 bp |
| Output | 24,576 bins of 32 bp, no cropping |
| Trunk | ConvDNA (768 ch, k=13) → ConvTower (→1536 ch, ×4) → HydraTower (1536 ch, ×8) |
| Activation | SiLU |
| Human head | 182 targets — 66 ATAC, 58 RNA+, 58 RNA− |
| Mouse head | 373 targets — 127 ATAC, 123 RNA+, 123 RNA− |

Coverage is square-root-summed within bins (`sum_sqrt`) and clipped at 200
(ATAC) or 300 (RNA); see `../data/README.md`.

**Training.** AdamW with schedule-free learning at lr 3 × 10⁻⁴, β = (0.9, 0.99),
weight decay 0.015 with none on biases and 10⁻⁴ on norm and head weights, 1,500
warmup steps, gradient clipping at global norm 1.0, batch size 2, a
multinomial-Poisson (`poisson_mn`) loss, bfloat16 mixed precision, up to 20
epochs. Both species train jointly through the shared trunk.

**Folds.** Eight models, `f0c0 … f7c0`. Model `f<I>c0` holds out data fold I as
test and fold I+1 as validation. Predictions and variant scores are averaged
across all eight, with reverse-complement averaging within each.

**Initialisation.** Each fold starts from the *matching* fold of the pretrained
Cerberus foundation model, `${CERBERUS_PRETRAINED_DIR}/f<I>c0/train/model_best.pth`.
Keeping the fold index aligned is what prevents the initialisation from having
seen this fold's test sequences during pretraining.

## Provenance of these scripts

`params.json` is the as-run configuration, copied from the trained model.

**`1_hound_data.sh` and `2_train_folds.sh` are reconstructions.** The original
invocations were issued interactively and no job script for them survives: there
is no `.sb` or `.sh` under the data or model directories other than the two
evaluation wrappers. They were rebuilt from `params.json`, the `statistics.json`
the data step wrote (sequence length, bin width, fold count, target counts), and
the job-script templates left by earlier, superseded training generations.

Two `hound_data` arguments are inferred rather than recorded: the write-job
batch size `-r 64` (from the on-disk job layout) and the worker count `-p 64`.
Neither changes the contents of the training data, only how it is produced.
Everything that does affect the data — length, bin width, folds, clipping,
summary statistic — comes from the preserved configs.

`3_eval_folds.sh` is a cleaned-up copy of a wrapper that did survive.

## 1. Training data

| Script | In | Out |
|--------|----|-----|
| `1_hound_data.sh {hg38\|mm10}` | targets file, genome FASTA, blacklist, umap, **plus the split files** | `${CERBERUS_DATA_*}/{seqs_cov/*.h5, examples/fold{0..15}.zarr, statistics.json}` |

`hound_data` reads each target's `.w5` coverage over the split intervals, bins
it at 32 bp, and writes the zarr examples the trainer consumes.

**The cross-fold splits are not generated here.** `sequences.bed`,
`contigs.bed`, `contig_components.txt`, `nets*.bed` and `mseqs_unmap.npy` were
taken verbatim from the pretrained Cerberus dataset and must be placed in the
output directory before running this step. This is deliberate: reusing the
foundation model's splits is what guarantees that a fold held out here was also
held out during pretraining. Generating fresh splits with `hound_data_align`
would break that correspondence and leak test sequences through the
initialisation. Obtain them with the pretrained model.

## 2–3. Training and evaluation

| Script | In | Out |
|--------|----|-----|
| `2_train_folds.sh` | `params.json`, both data directories, pretrained checkpoints | `${CERBERUS_MODEL_DIR}/f<I>c0/train/{model_best.pth,params.json,log.txt}` |
| `3_eval_folds.sh` | each fold model + both data directories | `f<I>c0/eval<g>/fold<j>/acc.txt` |

`2_train_folds.sh` writes a per-fold `params.json` with that fold's
initialisation checkpoint substituted in, then hands the folds to
`hound_train_folds`. The two data directories become **dataset 0 (human)** and
**dataset 1 (mouse)**; that numbering is the `--head` / `--dataset` argument
every downstream script uses.

`3_eval_folds.sh` evaluates every model fold against every data fold, for both
genomes. `acc.txt` gives per-track Pearson r and r². The accuracy figures read
only the diagonal — model `f<I>c0` on data fold I, its own held-out split —
but the full grid is cheap and makes fold-specific artefacts visible.

## `score/` — variant effects

| Stage | Script | Statistic |
|-------|--------|-----------|
| caQTL | `caqtl/1_prepare_variants.py`, `2_score_folds.sh` | `logSUM` over ±512 bp |
| eQTL | `eqtl/1_prepare_gtf.py`, `2_prepare_vcf.py`, `3_score_folds.sh` | `covgene/logSED` over the eGene's exons |
| GWAS | `gwas/1_score_credible_set.py`, `2_predict_tracks.py` | both, for one named gene and cell type |

Then, for either bulk stage:

```
python score/ensemble_folds.py \
    --fold_dir '<work>/caqtl/scores/logSUM/local_f{fold}c0' \
    --out_dir  '<work>/caqtl/scores/logSUM/local_ensemble' \
    --keys cov/logSUM --n_chunks 63
```

**`--local_window` needs the `_local` targets file.** The caQTL score is
`logSUM` restricted to a 1,024 bp window centred on the variant, which is what
makes it comparable to ChromBPNet's 2,114 bp receptive field. That restriction
is triggered by the `window` column in `${TARGETS_HUMAN_LOCAL}`; passing
`--local_window` with the plain targets file is **silently ignored** and yields
a whole-sequence score instead. The eQTL scoring deliberately uses the plain
file, because there the gene's exon set defines the region.

**Which variants get scored.** `caqtl/1_prepare_variants.py` keeps only variants
whose ±512 bp window fully contains their RASQUAL peak — otherwise the score
measures something different from the peak-level effect being compared against.
`eqtl/1_prepare_gtf.py` filters the annotation to eQTL target genes, so
`covgene/logSED` pairs each variant with its eGene rather than with every
overlapping gene.

**Allele orientation.** Scores are always log(ALT/REF) in hg38 orientation. The
eQTL table's effect allele is not always the reference; where `match_type ==
"swap"` the eQTL slope must be negated before comparison. See
`eqtl/2_prepare_vcf.py`.

## `interpret/` — in-silico mutagenesis

```
bash interpret/run_ism.sh {atac|rna} <variants.vcf> <out_dir>
```

Substitutes every base in a window with all three alternatives and re-scores:
`logSUM` over ±512 bp for accessibility (50 bp mutagenised), `covgene/logSED`
over the target gene's exons for expression (100 bp mutagenised). With no fold
argument it submits all eight folds; average them with `score/ensemble_folds.py`
using `--keys ref/cov/logSUM alt/cov/logSUM`.

The per-position contribution plotted in the figures is derived from these raw
scores: the score of the base actually present minus the mean of the three
alternatives.

## Track indexing

Three different track orderings appear in this pipeline and mixing them up
produces plausible-looking numbers for the wrong cell type. They are documented
in `../figures/figlib.py`, which is the single place they are defined.
`cerberus_model.CerberusFold.check_track` asserts a track's description matches
what the caller expects, so an index error fails loudly.

## Software

baskerville (PyTorch) with the `hound_*` entry points · PyTorch · pysam ·
numpy, pandas, h5py, zarr
