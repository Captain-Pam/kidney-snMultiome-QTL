# Open4Gene peak–gene linkage

This directory contains the peak–gene linkage analysis workflow using Open4Gene, a method previously developed by our group.

Peak–gene linkage analysis is computationally intensive. Although newer versions of Open4Gene use `fasthurdle` to improve computational efficiency, substantial computing resources and runtime may still be required.

## Input contract

- RNA and ATAC chunks follow `split_chunks_add_injPT/chunk_###/snRNA_chunk.h5ad` and `split_chunks_add_injPT/chunk_###/snATAC_chunk.h5ad`. Cells are observations, genes or peaks are variables, and the files share cell barcodes.
- Cell metadata is in `cell_rna_atac_metadata.tsv.gz`; its first column contains cell barcodes. Required columns are `cell_state`, `nCount_RNA`, `percent_mito`, and `nFrags`.
- Candidate pairs are in `peak_gene_pairs_tss_250kb_with_chunk.tsv`, with required columns `Peak`, `Gene`, and `chunk_name`.
- All genomic coordinates use hg38. Pairs are limited before Open4Gene to peaks within TSS ±250 kb; this workflow does not create or re-filter pairs.

## Analysis workflow

To make this computationally intensive analysis tractable, the input is divided into chunks that are processed independently as a Slurm array. Within each chunk, candidate peak–gene pairs are further divided into cell-type-specific batches and analyzed in parallel; the chunk-level results are then merged into the final output.

`1_run_open4gene_chunk.R` analyzes one chunk, `2_submit_open4gene_array.sh` submits the chunk jobs as a Slurm array, and `3_merge_open4gene_results.R` merges chunk results and computes the final FDR values.

Submit the array from this directory:

```bash
cd /path/to/kidney-snMultiome-QTL/Main/peak_gene_linkage
export INPUT_DIR=/path/to/input
export OUTPUT_DIR=/path/to/output
# Optional executable overrides:
export RSCRIPT_BIN=Rscript
export PYTHON_BIN=python
sbatch 2_submit_open4gene_array.sh
```

Run one chunk directly:

```bash
Rscript 1_run_open4gene_chunk.R --chunk-name chunk_001 --input-dir "$INPUT_DIR" --output-dir "$OUTPUT_DIR" --batches 20 --cores 20
```

Merge completed chunk results:

```bash
Rscript 3_merge_open4gene_results.R --result-dir "$OUTPUT_DIR"
```

## Fixed Open4Gene configuration

| Setting | Value |
|---|---|
| Covariates | `log_nCount_RNA`, `percent_mito`, `log_nFrags` |
| Gene filter | Fraction >0.01; at least 50 cells |
| Peak filter | Fraction >0.005; at least 50 cells |
| Binary | `FALSE` |
| Method | `fasthurdle` |
| MinNum.Cells | 5 |
| Batches per cell type | 20 |
| Slurm array | 1–60, concurrency limit 10 (`%10`) |
| CPUs | 20 |
| Memory | 200G |

## Outputs

| File | Description |
|---|---|
| `res_chunk_###.tsv` | Open4Gene results for an individual chunk. |
| `status_chunk_###.tsv` | Per-chunk run status. |
| `combine_peak_gene.txt.gz` | Combined results with `zero_FDR` and `count_FDR`, calculated by Benjamini–Hochberg correction separately within each `Celltype`. |
| `hit_fdr005_peak_gene.txt.gz` | Combined rows where either `zero_FDR` or `count_FDR` is <0.05. |

## Dependencies

R 4.4; R packages `Open4Gene`, `fasthurdle`, `argparser`, `data.table`, `Matrix`, and `reticulate`; Python package `scanpy`; and Slurm. The historical Open4Gene version is 1.1.0 from [hbliu/Open4Gene](https://github.com/hbliu/Open4Gene).
