# Cell Ranger ARC preprocessing

This directory contains the code used to align and quantify paired Gene
Expression and Chromatin Accessibility reads from single-nucleus multiome
libraries with Cell Ranger ARC.

## Workflow

`1_run_cellranger_arc.sh` runs `cellranger-arc count` for one multiome library
within a single Slurm job. Separate libraries are processed with independent
job submissions.

## Input and output

The script accepts a run identifier, a libraries CSV file, a compatible
reference directory and an output directory. The libraries CSV header must be
`fastqs,sample,library_type` and must describe both the `Gene Expression` and
`Chromatin Accessibility` FASTQ inputs for the library.

Cell Ranger ARC writes the analysis results to
`<output_directory>/<run_id>/outs`. These outputs include the feature-barcode
matrices and ATAC fragment files used in downstream analyses.

## Software and reference

- Cell Ranger ARC 2.0.2
- `refdata-cellranger-arc-GRCh38-2024-A`
- Slurm

## Resource parameters

Cell Ranger ARC was run in local mode within each Slurm job using 30 local
cores and 70 GB of local memory. The Slurm job requests 80 GB of memory to
provide additional memory for the parent process and system overhead. Default
Cell Ranger ARC alignment and counting parameters were otherwise retained.

## Notes

The Gene Expression matrix is used for CellBender background removal, whereas
the ATAC fragment output is used for downstream chromatin-accessibility
processing.
