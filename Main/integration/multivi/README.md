# MultiVI RNA--ATAC integration

This directory contains the workflow for integrating post-QC RNA and ATAC data
with MultiVI. Both AnnData and MuData input formats are supported.

## Input contract

`1_run_multivi_integration.py` accepts one of the following:

- `.h5ad`: a joint AnnData object with genes followed by peaks.
- `.h5mu`: a MuData object with modalities named `rna` and `atac`.

For both formats, count matrices must be stored in `X`, and `obs` must contain:

- `modality`: the observation-level modality/batch label used by MultiVI.
- `Sample`: the sample-level categorical technical covariate.

AnnData inputs must also contain `var["modality"]`, using
`Gene_Expression` for genes and `Peaks` for accessible regions. Input data are
expected to have completed all project-specific filtering before this script is
run.

## Preprocessing

The workflow applies these steps in order:

1. Select the AnnData or MuData branch from the input file extension.
2. Read the prepared joint object.
3. For MuData, update the joint metadata and retain observations shared by the
   `rna` and `atac` modalities.
4. Validate the required observation metadata and modality definitions.
5. Remove unused layers while retaining raw count matrices in `X`.

The script does not repeat feature selection or the earlier 0.5% feature
filtering step. Project-specific barcode and cluster removal must be completed
before running this workflow.

## Fixed MultiVI configuration

AnnData registration uses:

```python
scvi.model.MULTIVI.setup_anndata(
    data,
    batch_key="modality",
    categorical_covariate_keys=["Sample"],
)
```

MuData registration uses:

```python
scvi.model.MULTIVI.setup_mudata(
    data,
    modalities={"rna_layer": "rna", "atac_layer": "atac"},
    batch_key="modality",
    categorical_covariate_keys=["Sample"],
)
```

Training follows the final analysis script:

| Parameter | Value |
| --- | --- |
| Maximum epochs | 60 |
| Batch size | 4,096 |
| Validation interval | Every epoch |
| Early stopping | Enabled |
| Random seed | 0 |

## Run

AnnData example:

```bash
python 1_run_multivi_integration.py \
  --input-data prepared_joint_data.h5ad \
  --output-dir multivi_results
```

MuData example:

```bash
python 1_run_multivi_integration.py \
  --input-data prepared_joint_data.h5mu \
  --output-dir multivi_results
```

## Outputs

The output directory contains the following files:

| Output | Description |
| --- | --- |
| `trained_model/` | Saved MultiVI model |
| `mdata.h5ad` or `mdata.h5mu` | Integrated object with latent space, neighbors, UMAP, and Leiden clusters |
| `latent.npy` | MultiVI latent representation, also stored in `obsm["multivi"]` |
| `train.jpg` | Training and validation reconstruction-loss curves |
| `multivi_leiden.png` | Leiden clusters displayed on the MultiVI UMAP |
| `multivi_YYYYMMDD_HHMMSS.log` | Timestamped progress log |

Neighbors are calculated from `obsm["multivi"]`. UMAP uses `min_dist=0.1`,
and Leiden clustering uses the Scanpy defaults, matching the final analysis
script.

## Reproducibility notes

Reference environment metadata from the analysis environment:

- Python 3.12.2
- NumPy 1.26.4
- scvi-tools 1.2.2.post2
- Scanpy 1.10.4
- AnnData 0.11.3
