# ChromBPNet models

Per-cell-type, local-context models of chromatin accessibility, trained with
[ChromBPNet](https://github.com/kundajelab/chrombpnet) on the pseudobulk snATAC
fragments from `../data`. Each model reads **2,114 bp** of sequence and predicts
a **1,000 bp** base-resolution Tn5 insertion profile together with a scalar read
count. One model per cell type per fold.

Cell types: `CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT`.

All paths come from `../config.sh`; see `../README.md`. Everything here runs in
the `chrombpnet` (TensorFlow) environment — not the `s2f` environment used for
plotting.

## Layout

```
1_prepare_peaks.sh    MACS2 peak calling + GC-matched background regions
2_train_folds.sh      chrombpnet pipeline, 12 cell types x 5 folds
3_contribs_bw.sh      genome-wide contribution scores over peaks
4_modisco.sh          de novo motif discovery from those contributions
5_export_metrics.py   held-out accuracy tables for plotting
h5_to_npz.py          format shim used by 4_modisco.sh
score/                variant scoring via variant-scorer
interpret/            per-variant attributions for the example loci
```

## 1–2. Peaks and training

| Script | In | Out |
|--------|----|-----|
| `1_prepare_peaks.sh` | `${FRAG_DIR}/<ct>_merged_fragments.tsv.gz`, blacklist, fold-0 split | `${CHROMBPNET_PEAK_DIR}/<ct>_peaks.bed`, `<ct>_nonpeaks_negatives.bed` |
| `2_train_folds.sh` | peaks, background, fragments, per-fold bias model | `${CHROMBPNET_MODEL_DIR}/<ct>_fold_<k>/models/chrombpnet_nobias.h5` |

**Peak calling.** MACS2 with `--nomodel --shift -100 --extsize 200 --keep-dup all
-p 0.01 -g hs`. The shift/extend pair centres a 200 bp window on each Tn5
insertion; the permissive p-value is deliberate, since ChromBPNet wants a
generous peak set rather than a calibrated one.

**Blacklist.** Peaks overlapping the ENCODE hg38 blacklist **extended by 1,057 bp
on each side** are dropped. The extension is half the model's 2,114 bp receptive
field, so no training window can reach into a blacklisted region.

**Background.** `chrombpnet prep nonpeaks` samples GC-matched non-peak regions.
These are generated once against the fold-0 split and reused for all folds.

**Bias.** Tn5 sequence bias is handled by a pre-trained bias model (ENCODE
**ENCSR291GJU**, HepG2), held fixed while the accessibility component trains.
The bias model is fold-matched so none of its training chromosomes leak into a
given fold's test set. Every downstream stage uses the bias-free model,
`chrombpnet_nobias.h5`.

**Folds.** Five chromosome-held-out splits, e.g. fold 0 tests chr1/chr3/chr6 and
validates on chr8/chr20. **All predictions, variant scores and attributions are
averaged across the five folds.**

## 3–5. Contributions, motifs, metrics

| Script | In | Out |
|--------|----|-----|
| `3_contribs_bw.sh` | fold-0 model, the pipeline's `auxiliary/filtered.peaks.bed` | `${WORK_DIR}/chrombpnet/contribs_bw/<ct>.{counts,profile}_scores.h5`, bigWigs |
| `4_modisco.sh` | `<ct>.profile_scores.h5` | `${WORK_DIR}/chrombpnet/modisco/<ct>/modisco_results.h5` + `report/` |
| `5_export_metrics.py` | each run's `evaluation/chrombpnet_metrics.json`, `chrombpnet_predictions.h5`, `auxiliary/data_unstranded.bw` | `pearson_r.csv`, `scatter_data.csv` |

`5_export_metrics.py` only reads files training already wrote, so it needs no
GPU. Its two CSVs are the input to `../figures/chrombpnet_accuracy.py`.

TF-MoDISco clusters recurring high-attribution subsequences into motifs and then
matches them against a known-motif database (`${MODISCO_REPORT_MEME}`, HOCOMOCO
v12 core). This is genome-wide motif discovery; the per-variant attributions
shown alongside the example loci come from `interpret/` instead.

## `score/` — variant effects

| Script | In | Out |
|--------|----|-----|
| `1_prepare_caqtl_variants.sh` | the caQTL VCF from `../../cerberus/score/caqtl` | `caqtl_variants.tsv` (5-column variant-scorer table) |
| `2_prepare_eqtl_variants.py` | `${EQTL_SUSIE_TSV}` | `variants/<ct>_variants.tsv`, `variants/eqtl_cs_meta.tsv` |
| `3_score_folds.sh {caqtl\|eqtl}` | a variant list + all fold models | `scores_folds/<ct>/fold_<k>/<ct>_fold_<k>.variant_scores.tsv` |
| `4_ensemble.py` | those per-fold tables | `scores_ensemble/<ct>/<ct>.variant_scores.tsv` |

The caQTL VCF is built once, upstream, by the Cerberus side, so both models
score **exactly the same variant set** and their scores are directly comparable.

Scoring runs `variant_scoring.py` from a local clone of
[variant-scorer](https://github.com/kundajelab/variant-scorer)
(`${VARIANT_SCORER_SRC}`) with `-sc original -p <peaks> -n 10 --no_hdf5`. The
reported effect is **`logfc`**, the log fold-change in predicted total counts
between alleles. Passing the cell type's peak set normalises the score against
that cell type's own distribution, which is what makes the quantile columns
comparable across cell types.

`4_ensemble.py` averages every numeric column across folds and checks that all
folds scored the same variants in the same order — the folds differ only in
their chromosome split, so this is model ensembling, not aggregation over
different variants.

**Allele orientation.** Model scores are always `log(ALT/REF)` with REF the hg38
reference base. For eQTL variants whose QTL effect allele is the reference
(`match_type == "swap"`), `2_prepare_eqtl_variants.py` negates the eQTL slope
into `slope_adj` so both measures share that orientation. Comparing against the
raw `slope` silently flips the sign for those variants.

## `interpret/` — per-variant attributions

```
bash interpret/run_shap.sh <variants.tsv> <out_dir> [celltype ...]
```

`run_shap.py` computes DeepLIFT/SHAP contributions on both the reference and the
alternate allele sequence, averaged over folds, and saves a ±25 bp window around
each variant to `<ct>_ism.h5`. The wrapper submits one job per cell type and
derives its job count from the cell-type list, so adding a cell type cannot
leave it silently unscored.

Attribution is taken against the **counts** head, not the profile head: a caQTL
changes how many reads fall in a region, which is what the counts head predicts,
whereas the profile head only explains where within the window they fall.

Stored values are `shap * one_hot` — the attribution of the base actually
present — which is what an attribution logo plots.

## Software

ChromBPNet · TensorFlow · MACS2 2.2.9.1 · bedtools · TF-MoDISco-lite ·
[variant-scorer](https://github.com/kundajelab/variant-scorer) · Python 3 with
numpy, pandas, h5py, pyBigWig, hdf5plugin
