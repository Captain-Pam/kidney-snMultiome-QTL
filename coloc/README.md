# GWAS–eQTL colocalization

Colocalization of the **eGFR GWAS** (kidney function, meta-analysis, N ≈ 1.7M)
against **single-cell eQTLs** across kidney cell types, using the
[`coloc`](https://chr1swallace.github.io/coloc/) package. Two complementary
strategies are provided:

```
coloc_abf/    coloc.abf  — single-causal-variant, LD-free, runs on hg38 sumstats
coloc_susie/  coloc.susie — multiple causal variants, uses SuSiE fits from ../susie
```

A pair is called colocalized when **PP.H4 > 0.8**.

## Paths

`coloc_susie` reads the SuSiE fits produced by the sibling
`../../susie` stage via relative paths; raw sumstats are expected under `data/`
(same layout as `../susie/README.md`, plus `data/reference/gene_model_hg{19,38}.rds`
gene models for the locus plots).

## `coloc_abf/` (single causal variant)

Works directly on summary statistics in **hg38**; no LD panel required.

| File | Inputs | Outputs |
|------|--------|---------|
| `gene_range.R` | tensorQTL eQTL parquets (hg38, `egene_qval < 0.1`) | `gene_gr.rds` — hg38 gene window `GRanges` per cell type |
| `run_coloc.R` | eQTL parquets; expression BED (`sdY`); eGFR per-locus sumstats (`egfr_list.rds`, 100 kb windows); independent loci (`gr_indep.rds`); `gene_gr.rds` | `outputs/<ct>.csv` (all gene×locus tests: PP.H0–H4, top SNP, z-scores), `coloc_res/*.rds` when PP.H4 > 0.8, combined `coloc_0.8_z.csv` |
| `save_coloc.R` | Older variant of `run_coloc.R` on PMBB data paths; a pre-filtered `coloc_0.8.csv` | `coloc_res/*.rds`, `coloc_0.8_z.csv` (annotated with top-SNP REF/ALT and z-scores) |
| `plot.R` | `coloc_res/*.rds`; eGFR + eQTL sumstats; hg38 gene models | `plots/*.pdf` — eGFR vs sc-eQTL locus plots + gene track (`Gviz`) |

Per gene×locus: match SNPs across datasets handling **both allele orientations**
(`id1`/`id2`, flip eQTL `slope`), require ≥ 20 shared SNPs, then
`coloc.abf()` with both datasets as `quant` (eGFR uses `N`, `MAF`; eQTL uses `sdY`).

## `coloc_susie/` (multiple causal variants)

Consumes the fitted SuSiE objects from `../susie` (hg19).

| File | Inputs | Outputs |
|------|--------|---------|
| `run_coloc.R` | GWAS + eQTL SuSiE fits (`rss.rds`) with ≥ 1 credible set; `susie_meta.csv` / `susie_<ct>.csv`; gene/locus `GRanges` for overlap | Per colocalizing pair `outputs/<ct>_<locus>_<gene>/coloc_res.{rds,csv}` (saved when PP.H4 > 0.8); combined `coloc_all.csv` (every CS-pair test) |
| `filter_coloc.R` | `coloc_all.csv`; the saved `coloc_res.rds`; eQTL + eGFR sumstats | `coloc_0.8.csv` — PP.H4 > 0.8 pairs annotated with **top colocalizing SNP** (max `SNP.PP.H4`), REF/ALT, and directional `eQTL_z` / `eGFR_z` (sign-aligned to a common ALT allele) |
| `plot.R` | `coloc_0.8.csv`; eGFR + eQTL `res_susie.csv`; hg19 gene models | `plots/*.pdf` — 5-panel plots: eGFR −log10P, eGFR PIP, eQTL −log10Q, eQTL PIP, gene track; colored by credible set |

Per pair: match GWAS locus to overlapping genes (`findOverlaps`), intersect SNPs
(require ≥ 50 shared), run `coloc.susie(rss1, rss2)` which tests every credible-set
pair.

## coloc.abf vs coloc.susie

| | `coloc_abf` | `coloc_susie` |
|--|-------------|---------------|
| Causal variants | 1 per locus | multiple (per credible set) |
| LD panel | not needed | UKBB EUR hg19 (via `../susie`) |
| Genome build | hg38 (native) | hg19 (eQTLs lifted over) |
| GWAS window | 100 kb | 1 Mb |
| Min shared SNPs | 20 | 50 |
| Coloc call | `coloc.abf(d1, d2)` | `coloc.susie(rss1, rss2)` |

Both threshold at **PP.H4 > 0.8** and report top-SNP directional z-scores so the
effect direction (expression ↔ eGFR) can be interpreted.
