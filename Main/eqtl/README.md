# Single-cell cis-eQTL pipeline

Cell-type-resolved cis-eQTL mapping in kidney snMultiome pseudobulk, using
[tensorQTL](https://github.com/broadinstitute/tensorqtl). The pipeline runs in
three stages: **prepare data → tune hyperparameters → final run**.

Cell types: `CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT`.

All commands are run from this `eqtl/` directory. Absolute paths have been
replaced with placeholders under `data/` — edit them for your environment.

## Layout

```
prepare_data/   pseudobulk, QC, expression normalization, PEER, PCs, covariates
tune/           sweep genotype PCs + PEER factors on the eGene count
run/            final cis-eQTL run and assembly of the output tables
src/            shared tensorQTL wrappers (map_cis, map_nominal)
```

## 1. `prepare_data/`

| Script | In | Out |
|--------|----|-----|
| `1_pseudobulk.py` | `data/rna_clean.h5ad` | `anndata/counts/<ct>.csv` (summed counts), `anndata/cellcounts.csv` |
| `2_sample_qc.R` | pseudobulk counts, `data/genotype/final_QCed.fam` | `sample_qc/<ct>/` (QC'd counts, `sample_lookup.txt`), `sample_qc/sample_qc.csv` |
| `3_make_gct.R` | QC'd counts, `data/reference/gencode.v46...gtf.gz` | `sample_qc/<ct>/{raw,tpm}.gct`, `sample_qc/gencode_v46_rev.gtf` |
| `4_normalize_expression.sh` | GCTs + GTF (GTEx pipeline) | `bed/<ct>.expression.bed.gz` |
| `5_compute_peer.sh` | expression BED (GTEx `run_PEER.R`) | `peer/<ct>_peer<k>.PEER_covariates.txt`, k∈{5,10,15,20} |
| `6_genotype_pcs.py` | `data/genotype/genotype_pca.eigenvec` | `pca/PC_<n>.txt`, n∈1..20 |
| `7_make_covariates.py` | `data/metadata/multiome_metadata.txt`, `data/rna_clean.h5ad` | `covs/covs_<ct>.txt` |

**Sample/gene QC** (`2_sample_qc.R`): drop low-depth samples (log10 counts <
mean − 1.5 SD), keep genes expressed in > 20% of samples, drop expression PCA
outliers (|z(PC1/PC2)| > 2.5), keep samples with WGS genotypes.

**Expression** (`4_...sh`): GTEx TMM normalization, thresholds
`tpm > 0.1`, `count > 6`, expressed in ≥ 20% of samples.

**Covariates** (`7_...py`): Gender, Age, log2(nCount_RNA+1), log2(cell_count+1),
percent_mito, percent_ribo — combined with PEER factors at mapping time.
Genotype PCs are prepared for tuning only.

## 2. `tune/`

`run_tune.sh` sweeps the two nuisance-covariate counts, scoring each by the
**number of eGenes at permutation p < 0.05**:

- **Genotype PCs** — swept on one cell type (PTS). PC = 0 was optimal, so **no
  genotype PCs are used**.
- **PEER factors** — swept per cell type over k ∈ {5, 10, 15, 20}
  (`select_hyperparams.py` → `tune/peer_optimal.csv`).

## 3. `run/`

`run_eqtl.sh` runs, per cell type, two tensorQTL passes with **5 PEER factors**
(fixed across all cell types) and no genotype PCs:

- `src/map_cis.py` — permutation pass, window **±250 kb**, **MAF ≥ 0.1**;
  genome-wide FDR via `calculate_qvalues` (λ = 0.85) → per-gene `qval`,
  `pval_nominal_threshold`.
- `src/map_nominal.py` — nominal pass, window **±500 kb**, all variant-gene
  pairs (per chr).

`assemble_parquet.py` merges the two passes (filtering `af ≥ 0.1`,
`|start_distance| < 250 kb`, adding REF/ALT alleles) into:

- `run/all_parquets/<ct>.parquet` — all cis pairs with gene-level `egene_qval`.
- `run/fdr/<ct>.parquet` — significant eGenes only (`egene_qval < 0.1`).

These parquet tables are the inputs to the fine-mapping (`../susie`) and
colocalization (`../coloc`) pipelines.

## Expected `data/` inputs

| Path | Contents |
|------|----------|
| `data/rna_clean.h5ad` | QC'd RNA object (`obs.multivi_cell_type`, `obs.demultiplex_sample`, `nCount_RNA`, `percent_mito`, `percent_ribo`) |
| `data/metadata/multiome_metadata.txt` | `Gender`, `Age` per sample |
| `data/genotype/final_QCed.{bed,bim,fam}` | WGS genotypes (plink, hg38) |
| `data/genotype/genotype_pca.eigenvec` | genotype PCA (e.g. `plink2 --pca 20`) |
| `data/reference/gencode.v46.basic.annotation.gtf.gz` | gene annotation |
| `gtex-pipeline/` | external [GTEx qtl pipeline](https://github.com/broadinstitute/gtex-pipeline) (`eqtl_prepare_expression.py`, `run_PEER.R`) |

## Key parameters

cis window ±250 kb · MAF ≥ 0.1 · 5 PEER factors (all cell types) ·
genotype PCs = 0 · gene-level FDR `qval < 0.1` (Storey, λ = 0.85).
