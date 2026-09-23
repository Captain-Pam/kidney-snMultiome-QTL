# Cell-type-resolved cis-caQTL pipeline

This directory documents the cis-caQTL analysis of kidney single-nucleus
multiome data. The workflow starts from Cell Ranger ARC ATAC BAMs and phased
whole-genome sequencing variants and proceeds through WASP correction,
cell-type pseudobulk preparation, RASQUAL mapping and empirical-FDR caPeak
calling.

All commands are run from this `caqtl/` directory. Large data and reference
files are represented by repository-relative paths under `data/`; replace
these paths as appropriate for the computing environment.

Cell types: `CNT_CD_PC DCT DTL_ATL EC IC Immune PEC PTS Podocyte Stromal TAL injPT`.

## Layout

```text
prepare_data/   WASP correction, pseudobulk QC, PEER and RASQUAL inputs
tune/           PEER-factor and genotype-PC tuning record
run/            final RASQUAL tasks, result merging and caPeak calling
src/            shared empirical-FDR function
```

## Analysis workflow

### 1. WASP correction

`prepare_data/1_prepare_wasp_variants.sh` converts phased, non-imputed hg38
VCFs to the SNP and haplotype HDF5 files used by WASP. The next three scripts
extract each donor's ATAC reads from a library-level Cell Ranger ARC BAM, run
the complete WASP allele-swap and BWA-MEM remapping workflow, remove reads
marked with the duplicate flag, and split each corrected donor BAM into the 12
cell types.

The donor barcode lists used before WASP and the donor-by-cell-type barcode
lists used after WASP are separate inputs.

### 2. Pseudobulk, sample QC and covariates

`prepare_data/5_make_pseudobulk.py` aggregates the peak matrix by donor and
cell type. `6_sample_qc.R` applies the study filters:

- remove samples below mean log10 pseudobulk depth minus two standard
  deviations;
- calculate PCA from chr1 peaks detected in more than 20% of donors after
  library-size scaling and `log2(value + 1)` transformation;
- remove samples with `|z(PC1)| > 2.5` or `|z(PC2)| > 2.5`;
- retain donors with matched WGS data; and
- retain peaks with at least three fragments in at least 5% of donors, with a
  minimum of five donors.

PEER factors are calculated from chr1 peaks with TPM at least 0.1 and raw
count at least 6 in at least 20% of donors. The PEER matrix is quantile
normalized and each peak is inverse-normal transformed across donors.

Final covariates are age, numeric sex, z-scored `log10(nucleus count + 1)`,
z-scored median TSS enrichment, z-scored median nucleosome ratio and two PEER
factors. Genotype PCs were evaluated during tuning but were not included in
the final model.

### 3. RASQUAL mapping

RASQUAL is run across the 22 autosomes using a cis window of plus or minus
10 kb around each peak midpoint. Final tasks comprise:

- observed lead mode;
- two independent phenotype-permutation runs using `-r`; and
- observed full mode for all filtered peaks.

The RASQUAL version used in the study seeds its random-number generator from
the current time and process ID; fixed permutation seeds are therefore not
specified. Only fits with convergence status 0 are used downstream.

### 4. Empirical FDR

Observed lead results are compared with the equally pooled results from the
two permutation runs, separately for each cell type. Peaks passing empirical
FDR < 0.1 are designated caPeaks, and all converged full-mode associations for
those peaks are retained.

## Expected inputs

| Path | Contents |
|---|---|
| `data/atac_clean.h5ad` | Peak-by-nucleus count matrix; `obs` includes `demultiplex_sample`, `multivi_cell_type2`, `TSSEnrichment` and `NucleosomeRatio` |
| `data/manifests/donor_bams.tsv` | Header plus `donor_id` and library-level `bam_path` columns |
| `data/barcodes/donor/<donor>.txt` | All ATAC barcodes assigned to one donor |
| `data/barcodes/cell_type/<donor>/<cell_type>.txt` | Donor-by-cell-type barcode lists |
| `data/metadata/donor_metadata.tsv` | `sample_id`, numeric `Age` and numeric `Gender` |
| `data/genotype/phased/chr<1-22>.vcf.gz` | Phased, biallelic, non-imputed hg38 VCFs |
| `data/genotype/final_QCed.fam` | Donors retained after WGS QC |
| `data/genotype/genotype_pca.eigenvec` | Precomputed PLINK genotype PCA used only for tuning |
| `data/genotype/qc_variants.txt` | One hg38 `chr:position:allele1:allele2` ID per line; both allele orders are included |
| `data/reference/genome.fa` | hg38 FASTA with BWA index |
| `data/reference/chromInfo.hg38.txt.gz` | hg38 chromosome metadata for WASP `snp2h5` |

The scripts also use external installations of WASP 0.3.4, subset-bam 1.1.0,
BWA-MEM, samtools, bcftools/tabix, PEER, RASQUAL and the R package
`rasqualTools`. Software locations are grouped near the beginning of the
scripts that use them.

## Canonical ordering

For each cell type, `samples.qc.txt` defines the donor order used by the BAM
list, count matrix, covariate matrix and VCF. `master_peaks.tsv` is sorted by
chr1-chr22, start and end and contains the stable one-based `rasqual_index`
used by all RASQUAL tasks.
