# scTWAS and FOCUS workflow

This directory contains the final 250 kb cis-window workflow used to train kidney cell-type scTWAS models, test GWAS associations, and perform FOCUS fine-mapping.

## Input contract

All genotype, GWAS, gene-annotation, and LD-reference coordinates must use the same hg38 build.

The Stage 1 inputs are:

- Per-cell-type raw pseudobulk count matrices in CSV format, with genes as rows and donors as columns.
- A QCed donor genotype in PLINK `bed/bim/fam` format.
- eQTL covariate files containing `Gender` and `Age` rows.
- PEER covariate files containing at least five factors.
- A gene-annotation table containing `seqnames`, `gene_type`, `gene_symbol`, `TSS`, `start`, and `end`.

Before running Stage 1, create `1_stage1_weights/2_data/celltypes.tsv` with the tab-separated columns `celltype`, `counts_file`, `covariate_file`, and `peer5_file`. The manifest maps each of the 12 cell types to its pseudobulk count, covariate, and PEER5 file names. The covariate files are read from `<eqtl_dir>/3_covariates`, and the PEER files are read from `<eqtl_dir>/5_peer`.

Before running Stage 2, create `2_stage2_focus/2_data/gwas_manifest.tsv` with the tab-separated columns `trait`, `input_file`, `has_header`, `chr_col`, `pos_col`, `a1_col`, `a2_col`, `beta_col`, `se_col`, `z_col`, `p_col`, `n_col`, and `notes`.

GWAS column indices are 1-based. Supply `z_col` directly when available; otherwise leave it blank and provide `beta_col` and `se_col`. The GWAS chromosome and position columns must contain hg38 coordinates. `N` is required for the FOCUS munging step.

## Stage 1: train cell-type weights

Set the workflow directory to the absolute path of this directory:

```bash
export SCTWAS_DIR=/path/to/kidney-snMultiome-QTL/Main/sctwas
```

Build donor-level SCT objects:

```bash
Rscript 1_stage1_weights/1_build_sct_objects.R \
  --workflow-dir "$SCTWAS_DIR" \
  --counts-dir /path/to/pseudobulk_counts \
  --fam /path/to/qced_genotype.fam \
  --celltype all
```

Build the Gender, Age, and PEER1-PEER5 quantitative covariates:

```bash
Rscript 1_stage1_weights/2_build_qcovar_all_celltypes.R \
  --workflow-dir "$SCTWAS_DIR" \
  --eqtl-dir /path/to/eqtl_workflow \
  --celltype all
```

Build autosomal gene manifests and TSS +/- 250 kb cis ranges:

```bash
Rscript 1_stage1_weights/3_build_gene_manifest_250kb.R \
  --workflow-dir "$SCTWAS_DIR" \
  --gencode /path/to/gene_annotation.tsv.gz \
  --celltype all
```

Prepare chromosome-specific genotype files with MAF >= 0.05:

```bash
PLINK_BIN=/path/to/plink \
  bash 1_stage1_weights/4_prepare_genotype_by_chr.sh \
  "$SCTWAS_DIR" /path/to/qced_genotype
```

Submit Stage 1 from its code directory so that Slurm logs are written there:

```bash
cd "$SCTWAS_DIR/1_stage1_weights"
export RSCRIPT_BIN=/path/to/Rscript
export PLINK_BIN=/path/to/plink
sbatch 6_submit_stage1_array.sh
```

The fixed training configuration is:

| Parameter | Value |
| --- | --- |
| Cell types | 12 |
| Chromosomes | 1-22 |
| Cis half-window | 250 kb |
| Genotype MAF | 0.05 |
| Cross-validation folds | 5 |
| SCT iterations | 5 |
| Elastic-net alpha | 0.5 |
| Slurm array | 264 tasks |

After all chromosome jobs finish, combine the summaries and build one `weights_index.txt` per cell type. Models with `cv_R2_max > 0` are retained for Stage 2.

```bash
Rscript 7_combine_stage1_weights.R \
  --workflow-dir "$SCTWAS_DIR" \
  --celltype all
```

## Stage 2: test GWAS associations

Return to the workflow directory, standardize the GWAS inputs, and build the hg38 1KGP EUR LD reference over the union of the 250 kb cis ranges:

```bash
cd "$SCTWAS_DIR"

Rscript 2_stage2_focus/1_prepare_gwas_sumstats.R \
  --workflow-dir "$SCTWAS_DIR" \
  --trait all

PLINK_BIN=/path/to/plink \
  bash 2_stage2_focus/2_build_ldref_250kb_union.sh \
  "$SCTWAS_DIR" /path/to/1kgp_eur_hg38 /path/to/eur.keep

Rscript 2_stage2_focus/3_make_stage2_task_table.R \
  --workflow-dir "$SCTWAS_DIR" \
  --trait all
```

Submit one Slurm task per row of `2_stage2_focus/2_data/stage2_tasks.tsv`. Replace `N` with the number of data rows in that file, and use the same array range for FOCUS:

```bash
cd "$SCTWAS_DIR/2_stage2_focus"
export RSCRIPT_BIN=/path/to/Rscript
sbatch --array=1-N 5_submit_stage2_array.sh
```

Each task runs `twas_stage2_sumstats()` across chromosomes and combines the available chromosome-level results for one GWAS-cell-type pair.

## FOCUS fine-mapping

For each GWAS trait, munge the standardized summary statistics:

```bash
FOCUS_BIN=/path/to/focus \
  bash 6_munge_focus_sumstats.sh "$SCTWAS_DIR" example_trait
```

Build one FOCUS SQLite weight database per cell type:

```bash
Rscript 7_build_focus_db_celltype.R \
  --workflow-dir "$SCTWAS_DIR" \
  --celltype PTS
```

Repeat the database command for the required cell types. Then submit FOCUS from the Stage 2 code directory:

```bash
export FOCUS_BIN=/path/to/focus
export LD_BLOCK_BED=/path/to/grch38.eur.loci.bed
sbatch --array=1-N 9_submit_focus_array.sh
```

FOCUS runs on chromosomes 1-22 with `--p-threshold 5e-8`. Each array task combines the chromosome-level files into `<celltype>.focus_all.tsv`.

## Basic result summary

Summarize the TWAS and FOCUS results after all jobs finish:

```bash
Rscript 10_summarize_trait_celltype_results.R \
  --workflow-dir "$SCTWAS_DIR" \
  --trait all \
  --celltype all
```

The summary annotates TWAS results with BH-adjusted P values after requiring `MODELCV.R2 > 0.01`, and reports FOCUS PIP counts above 0.5, 0.8, and 0.9.

## Outputs

| Output | Description |
| --- | --- |
| `1_stage1_weights/2_data/sct_objects/` | Donor-level SCT objects |
| `1_stage1_weights/2_data/qcovar/` | Gender, Age, and PEER1-PEER5 covariates |
| `1_stage1_weights/2_data/gene_manifests/` | Gene tables and 250 kb cis ranges |
| `1_stage1_weights/3_outputs/weights/` | Cell-type scTWAS weights and indices |
| `2_stage2_focus/2_data/sumstats_by_trait/` | Standardized scTWAS and FOCUS GWAS inputs |
| `2_stage2_focus/2_data/ldref/` | Chromosome-specific hg38 EUR LD references |
| `2_stage2_focus/3_outputs/<trait>/stage2/` | scTWAS Stage 2 results |
| `2_stage2_focus/3_outputs/focus_db/` | Cell-type FOCUS weight databases |
| `2_stage2_focus/3_outputs/<trait>/focus/` | FOCUS fine-mapping results |
| `2_stage2_focus/3_outputs/summaries/` | Combined trait-by-cell-type summaries |

## Dependencies

- R 4.4
- R packages: `argparser`, `Seurat`, `Matrix`, `data.table`, `scTWAS`, `DBI`, and `RSQLite`
- PLINK with `--extract range` support
- FOCUS with the hg38 EUR LD-block BED file
