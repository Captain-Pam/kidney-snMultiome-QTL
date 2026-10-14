# scVI RNA integration

This directory contains the workflow for integrating post-QC
single-cell or single-nucleus RNA data with scVI.

## Input contract

The input must be a post-QC RNA `.h5ad` file with raw, non-negative integer
counts in `adata.X`. The workflow validates this requirement before it copies
the matrix to `adata.layers["counts"]`.

The following `adata.obs` columns are required:

| Column | Role |
| --- | --- |
| `orig_ident` | Primary library-level technical batch, passed as `batch_key` |
| `percent_mito` | Continuous mitochondrial-content covariate |
| `nCount_RNA` | Continuous total-RNA-count covariate |

`percent_mito` and `nCount_RNA` must contain numeric, finite values.

## Preprocessing

The workflow applies these steps in order:

1. Validate the required metadata and raw counts.
2. Copy raw counts from `.X` to `layers["counts"]`.
3. Retain cells with at least 200 detected genes.
4. Retain genes detected in at least 3 cells.
5. Normalize each cell to a total of 10,000 counts.
6. Apply `log1p` to `.X`.
7. Select 5,000 highly variable genes (`n_top_genes=5000`) with
   `flavor="seurat_v3"`, using the raw `counts` layer and `orig_ident` as the
   batch key.

The normalized and log-transformed matrix remains in `.X`; scVI is trained from
the raw `counts` layer.

## Fixed scVI configuration

AnnData registration uses:

```python
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="orig_ident",
    continuous_covariate_keys=["percent_mito", "nCount_RNA"],
)
```

The model and training parameters are fixed as follows:

| Parameter | Value |
| --- | --- |
| Hidden layers | 3 |
| Latent dimensions | 30 |
| Gene likelihood | `zinb` |
| Dropout rate | 0.1 |
| Maximum epochs | 150 |
| Learning rate | `1e-4` |
| Batch size | 2,048 |
| Early stopping | Enabled |
| Early-stopping patience | 25 epochs |
| Validation interval | Every epoch |
| Random seed | 0 |

The latent representation is stored in `obsm["X_scVI"]`. Neighbors are
computed from `X_scVI`; Leiden clustering uses resolution 1, the `igraph`
flavor, 2 iterations, and `directed=False`. UMAP uses `min_dist=0.1`,
`spread=1`, and random seed 0. Its coordinates are retained in both
`obsm["X_umap"]` and `obsm["X_umap_scvi"]`.

## Optional barcode exclusion

The script includes a commented example immediately after data loading. To use
it, supply a headerless one-column text file containing barcodes. The example
reads them as strings and excludes entries that match `adata.obs_names` before
preprocessing. Project-specific barcode files and barcode values are not
included here.

## Run

```bash
python 1_run_scvi_integration.py \
  --input-h5ad post_qc_rna.h5ad \
  --output-dir scvi_results
```

The output directory contains deterministic names:

| Output | Description |
| --- | --- |
| `scvi_model/` | Trained scVI model |
| `scvi_integrated.h5ad` | Processed AnnData with `X_scVI`, Leiden, and UMAP results |
| `training_history.png` | Training and validation reconstruction-loss history |
| `scvi_workflow.log` | Workflow stages and software versions |

## Reproducibility notes

Reference environment metadata from the analysis environment:

- Python 3.12.2
- NumPy 1.26.4
- scvi-tools 1.2.2.post2
- Scanpy 1.10.4
- AnnData 0.11.3
