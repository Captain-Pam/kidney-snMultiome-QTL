# ChromBPNet models

Per-cell-type, local-context models of chromatin accessibility, trained with
[ChromBPNet](https://github.com/kundajelab/chrombpnet) on the pseudobulk snATAC
fragments from `../data`. One model per cell type per fold. Architecture,
hyperparameters and peak-calling settings are in the paper methods.

Cell types: `CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT`.

Paths come from `../config.sh`. Everything here runs in the `chrombpnet`
environment.

## Training and interpretation

| Script | Does |
|--------|------|
| `1_prepare_peaks.sh` | calls peaks with MACS2, filters the blacklist, and samples GC-matched background regions |
| `2_train_folds.sh` | trains a bias-factorized model per cell type × fold; downstream stages use the bias-free `chrombpnet_nobias.h5` |
| `3_contribs_bw.sh` | genome-wide contribution scores over each cell type's peaks |
| `4_modisco.sh` | discovers motifs from those contributions and matches them to a known-motif database |
| `5_export_metrics.py` | collects held-out accuracy into two small CSVs for plotting |
| `h5_to_npz.py` | format shim between steps 3 and 4 |

`5_export_metrics.py` only reads files training already wrote, so it needs no
GPU. Its output feeds `../figures/chrombpnet_accuracy.py`.

## `score/` — variant effects

| Script | Does |
|--------|------|
| `1_prepare_caqtl_variants.sh` | converts the shared caQTL VCF to the flat table variant-scorer reads |
| `2_prepare_eqtl_variants.py` | builds per-cell-type lists of fine-mapped eQTL SNPs, with the eQTL slope aligned to ALT-vs-REF |
| `3_score_folds.sh {caqtl\|eqtl}` | scores a variant list with every fold model |
| `4_ensemble.py` | averages the per-fold scores |

Scoring uses `variant_scoring.py` from a local
[variant-scorer](https://github.com/kundajelab/variant-scorer) clone
(`${VARIANT_SCORER_SRC}`). The reported effect is `logfc`.

Two things worth knowing:

- The caQTL VCF is built **once**, upstream, by
  `../../cerberus/score/caqtl/1_prepare_variants.py`, and reused here. Both
  models therefore score exactly the same variant set and their scores are
  directly comparable.
- `4_ensemble.py` checks that all folds scored the same variants in the same
  order. The folds differ only in their chromosome split, so this is model
  ensembling, not aggregation over different variants.

## `interpret/` — per-variant attributions

```
bash interpret/run_shap.sh <variants.tsv> <out_dir> [celltype ...]
```

Computes DeepLIFT/SHAP contributions on both alleles, averaged over folds, and
saves a window around each variant. Attribution is taken against the **counts**
head: a caQTL changes how many reads fall in a region, which is what that head
predicts, whereas the profile head only explains where within the window they
fall.

The wrapper derives its job count from the cell-type list, so adding a cell type
cannot leave it silently unscored.

## Software

ChromBPNet, TensorFlow, MACS2, bedtools, TF-MoDISco-lite,
[variant-scorer](https://github.com/kundajelab/variant-scorer).
