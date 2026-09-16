# CellBender preprocessing

This directory contains the code used to extract the Gene Expression modality
from single-nucleus multiome data and remove ambient RNA with CellBender.

## Workflow

1. `1_extract_gex_from_multiome.R` extracts Gene Expression features from a
   Cell Ranger ARC multiome matrix and writes a gene-expression-only 10x matrix.
2. `2_run_cellbender.sh` applies CellBender `remove-background` to the extracted
   gene-expression matrix.
3. `3_compare_raw_and_cellbender.R` compares the extracted pre-CellBender
   counts with the CellBender-filtered counts across matched samples.

## Input and output

`1_extract_gex_from_multiome.R` reads a Cell Ranger ARC
`raw_feature_bc_matrix` directory containing `barcodes.tsv.gz`,
`features.tsv.gz`, and `matrix.mtx.gz`. It writes the same three files after
retaining only features labelled as `Gene Expression`.

`2_run_cellbender.sh` reads the gene-expression-only matrix directory and
writes a CellBender-filtered HDF5 matrix and a checkpoint archive.

`3_compare_raw_and_cellbender.R` matches sample-named input directories,
creates a dual-assay Seurat object containing pre- and post-CellBender counts,
and writes the following QC outputs:

- `dual_seurat.qs`: dual-assay Seurat object used for the comparison.
- `Raw_counts_vs_cellbender_median_stats.csv`: per-sample median changes in
  detected counts and features.
- `Raw_counts_vs_cellbender_most_changed_features.csv`: feature-level count
  changes across matched samples.
- `Raw_counts_vs_cellbender.jpg`: feature-level comparisons using 25% and 50%
  count-difference thresholds.

## Software

- R 4.4.2
- data.table 1.15.4
- gzip 1.12
- CellBender 0.3.2
- Seurat
- scCustomize
- qs
- ggplot2
- patchwork
- Slurm

## Initial CellBender parameters

The values below were used as the initial parameters. Parameters were adjusted
for individual samples based on CellBender diagnostic results.

| Parameter | Value |
| --- | --- |
| Learning rate | 0.0001 |
| Epochs | 150 |
| False positive rate | 0.01 |
| Low-count threshold | 100 |

## Notes

CellBender was applied only to the Gene Expression modality. ATAC peak counts
are not modified by the scripts in this directory.
