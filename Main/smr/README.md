# SMR/HEIDI workflow

This folder contains the core pairwise GWAS-eQTL workflow used to prepare hg38 inputs and run SMR/HEIDI.

All variant identifiers in the QTL, GWAS, and LD-reference files must represent the same hg38 coordinates.

## 1. Prepare the eQTL BESD files

`1_prepare_eqtl/1_build_eqtl_qfile.R` converts one tensorQTL Parquet file into an SMR qfile. The required tensorQTL columns are `phenotype_id`, `variant_id`, `af`, `A1`, `A2`, `slope`, `slope_se`, and `pval_nominal`. The gene annotation must provide `seqnames`, `TSS`, `strand`, `gene_id`, and `gene_symbol`.

The output columns are `SNP`, `Chr`, `BP`, `A1`, `A2`, `Freq`, `Probe`, `Probe_Chr`, `Probe_bp`, `Gene`, `Orientation`, `b`, `se`, and `p`. Probes with invalid annotation coordinates are removed within this script.

```bash
Rscript 1_prepare_eqtl/1_build_eqtl_qfile.R \
  --input cell_type.tensorqtl.parquet \
  --annotation gene_annotation.tsv \
  --output results/cell_type.qfile.tsv
```

Convert the qfile to dense BESD format:

```bash
SMR_BIN=/path/to/smr bash 1_prepare_eqtl/2_make_eqtl_besd.sh \
  results/cell_type.qfile.tsv \
  results/cell_type_eqtl
```

## 2. Prepare the GWAS input

`2_prepare_gwas/1_build_gwas_input.R` expects a unified hg38 GWAS table containing `SNP`, `A1`, `A2`, `FREQ_A1`, `BETA`, `SE`, `P`, `N`, `hg38_chr`, and `hg38_bp`. It writes the eight SMR columns `SNP`, `A1`, `A2`, `freq`, `b`, `se`, `p`, and `n`, with `SNP` rebuilt as hg38 `chr:bp`.

For duplicated `chr:bp` identifiers, the script retains one row deterministically by preferring single-nucleotide alleles, smaller P values, and larger sample sizes. `--freq_impute` is optional and replaces only missing allele frequencies.

```bash
Rscript 2_prepare_gwas/1_build_gwas_input.R \
  --input trait.unified.with_hg38.tsv.gz \
  --output results/trait.smr.tsv

# Optional, when a justified frequency value must replace missing FREQ_A1 values:
Rscript 2_prepare_gwas/1_build_gwas_input.R \
  --input trait.unified.with_hg38.tsv.gz \
  --output results/trait.smr.tsv \
  --freq_impute 0.5
```

## 3. Run SMR/HEIDI

The runner performs one GWAS-eQTL comparison. It uses `--peqtl-smr 5e-8`, requires `diff_freq` as the fifth argument, and otherwise uses the SMR 1.4.0 default HEIDI settings.

```bash
SMR_BIN=/path/to/smr bash 3_run_smr/1_run_smr.sh \
  results/trait.smr.tsv \
  results/cell_type_eqtl \
  /path/to/hg38_ld_reference \
  results/trait_cell_type \
  0.2
```

## Dependencies

- R 4.4 with `argparser`, `arrow`, and `data.table`
- SMR 1.4.0
- An hg38 PLINK LD reference (`.bed`, `.bim`, and `.fam`)
