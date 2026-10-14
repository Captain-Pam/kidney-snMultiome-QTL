# PeakVI ATAC integration

This directory contains the workflow for integrating post-QC
single-cell or single-nucleus ATAC data with PeakVI.

## Input contract

The input must be one post-QC `.h5ad` file with cells as observations, called
ATAC peaks as variables, raw finite non-negative integer peak counts in `adata.X`, and
`Sample` in `adata.obs` as the library-level technical batch.

Normalized or transformed data are invalid. A truly raw binary matrix remains
mathematically valid because its entries are non-negative integer peak counts.

## Preprocessing

The workflow applies these steps in order:

1. Validate the entire input, including metadata and raw counts.
2. Optionally exclude barcodes with the disabled, user-enabled example that
   reads a headerless one-column file.
3. Require at least one remaining cell.
4. Compute `min_cells=max(1, int(n_obs * 0.01))`.
5. Filter peaks using `min_cells`.
6. Binarize positive counts.
7. Require at least one remaining peak.

Validating the complete input before removal and binarization prevents excluded
cells or the binary conversion from masking invalid values.

## Fixed PeakVI configuration

AnnData registration and model construction use:

```python
scvi.model.PEAKVI.setup_anndata(adata, batch_key="Sample")
model = scvi.model.PEAKVI(
    adata,
    model_depth=True,
)
```

`n_hidden` and `n_latent` are not fixed in the original workflow. They are
automatically derived from the number of retained peaks by PeakVI.
`model_depth=True` is shown explicitly because PeakVI enables it by default to
model cell-specific accessibility depth.

No additional categorical or continuous covariates are registered. RNA
metrics such as mitochondrial percentage and RNA UMI counts are inappropriate
for this ATAC workflow. PeakVI models cell-specific accessibility depth
internally with `model_depth=True`, so `nFrags` is not added as a continuous
covariate.

| Parameter | Value |
| --- | --- |
| Maximum epochs | 350 |
| Learning rate | `1e-4` |
| Batch size | 1,024 |
| Early stopping | Enabled |
| Early-stopping patience | 50 epochs |
| Validation interval | Every epoch |
| Random seed | 0 |

The latent representation is stored in `obsm["X_peakVI"]`, and neighbors are
computed from it. Leiden clustering uses resolution 1, the `igraph` flavor, 2
iterations, `directed=False`, and random seed 0. UMAP uses `min_dist=0.1`,
`spread=1`, and random seed 0. Its coordinates are stored in `obsm["X_umap"]`
and copied to `obsm["X_umap_peakvi"]`.

## Optional barcode exclusion

The script includes a disabled example immediately after complete-input
validation. When enabled by the user, it reads a headerless one-column string
file and excludes entries that match `adata.obs_names`. No project barcode
values or historical MultiVI cluster list is included or applied.

## Run

```bash
python 1_run_peakvi_integration.py \
  --input-h5ad post_qc_atac.h5ad \
  --output-dir peakvi_results
```

The output directory contains deterministic names:

| Output | Description |
| --- | --- |
| `peakvi_model/` | Trained PeakVI model |
| `peakvi_integrated.h5ad` | Integrated object with `X_peakVI`, the neighbor graph, Leiden, `X_umap`, and `X_umap_peakvi` |
| `latent.npy` | PeakVI latent representation, also stored in `obsm["X_peakVI"]` |
| `training_history.png` | Training and validation reconstruction-loss history |
| `peakvi_leiden.png` | Leiden clusters displayed on the PeakVI UMAP |
| `peakvi_workflow.log` | Workflow stages and software versions |

## Reproducibility notes

Reference environment metadata from the analysis environment:

- Python 3.12.2
- NumPy 1.26.4
- scvi-tools 1.2.2.post2
- Scanpy 1.10.4
- AnnData 0.11.3
