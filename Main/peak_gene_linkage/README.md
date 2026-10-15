# Open4Gene peak–gene linkage

Public workflow for reproducing `3_open4gene_results_2` with the original Open4Gene implementation.

## Scripts

| Script | Purpose |
| --- | --- |
| `1_run_open4gene_chunk.R` | Run Open4Gene for one input chunk. |
| `2_submit_open4gene_array.sh` | Submit the 60-chunk Slurm array. |
| `3_merge_open4gene_results.R` | Combine chunk results and calculate FDR. |

## Inputs

Arrange the inputs under `INPUT_DIR`:

```text
INPUT_DIR/
├── split_chunks_add_injPT/
│   ├── chunk_001/
│   │   ├── snRNA_chunk.h5ad
│   │   └── snATAC_chunk.h5ad
│   ├── ...
│   └── chunk_060/
│       ├── snRNA_chunk.h5ad
│       └── snATAC_chunk.h5ad
├── cell_rna_atac_metadata.tsv.gz
└── peak_gene_pairs_tss_250kb_with_chunk.tsv
```

Coordinates use hg38. H5AD observations are cells and variables are features; RNA and ATAC inputs share cell barcodes. The metadata first column contains barcodes, with additional columns `cell_state`, `nCount_RNA`, `percent_mito`, and `nFrags`. The candidate-pair table contains `Peak`, `Gene`, and `chunk_name`.

Candidate pairs are **already filtered to TSS ±250 kb before Open4Gene**. This workflow does not generate candidate pairs or apply a post hoc distance filter.

## Usage

Set input and output directories, optionally select executable paths, and submit from this directory:

```bash
export INPUT_DIR=/path/to/inputs
export OUTPUT_DIR=/path/to/outputs
# Optional executable overrides:
export RSCRIPT_BIN=/path/to/Rscript
export PYTHON_BIN=/path/to/python
sbatch 2_submit_open4gene_array.sh
```

Alternatively, run one chunk directly:

```bash
Rscript 1_run_open4gene_chunk.R --chunk-name chunk_001 --input-dir "$INPUT_DIR" --output-dir "$OUTPUT_DIR" --batches 20 --cores 20
```

After all chunks finish, merge the results:

```bash
Rscript 3_merge_open4gene_results.R --result-dir "$OUTPUT_DIR"
```

## Outputs

The workflow writes per-chunk `res_chunk_###.tsv` and `status_chunk_###.tsv` files under `OUTPUT_DIR`. Slurm writes `open4gene_<job>_<task>.out` and `open4gene_<job>_<task>.err` logs in the submission directory. The merge step writes `combine_peak_gene.txt.gz` and `hit_fdr005_peak_gene.txt.gz` in the result directory. Benjamini–Hochberg correction is performed separately for `zero_FDR` and `count_FDR` within each `Celltype`; the hit table retains pairs with either value < 0.05.

## Fixed parameters

- Covariates: `log_nCount_RNA`, `percent_mito`, and `log_nFrags`.
- Gene detection: proportion > 0.01 and at least 50 cells.
- Peak detection: proportion > 0.005 and at least 50 cells.
- Open4Gene: `Binary = FALSE`, `Method = "fasthurdle"`, `MinNum.Cells = 5`, and 20 batches per cell type.
- Slurm: array `1-60%10`, 20 CPUs, 200 GB memory, and five days per task.

## Scope

This public workflow excludes input data, generated results, logs, plots, `7_power`, and patched Open4Gene implementations.

## Dependencies

- R 4.4 with Open4Gene, fasthurdle, argparser, data.table, Matrix, and reticulate.
- Python with scanpy.
- Slurm for array submission.

The historical analysis used Open4Gene 1.1.0 from `hbliu/Open4Gene`, commit `6e7f36aa80e6ff778737fc25978337342394cc68`.
