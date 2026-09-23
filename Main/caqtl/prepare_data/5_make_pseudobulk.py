#!/usr/bin/env python
"""Step 5 - Aggregate peak counts and ATAC QC by donor and cell type.

Input:
  data/atac_clean.h5ad
Outputs, per cell type:
  prepare_data/output/pseudobulk/<cell_type>/counts.mtx.gz
  prepare_data/output/pseudobulk/<cell_type>/peaks.tsv
  prepare_data/output/pseudobulk/<cell_type>/samples.tsv
  prepare_data/output/pseudobulk/<cell_type>/qc_summary.tsv

Run from the caqtl/ directory:
  python prepare_data/5_make_pseudobulk.py
"""

from pathlib import Path
import gzip

import anndata
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


INPUT_H5AD = Path("data/atac_clean.h5ad")
OUTPUT_ROOT = Path("prepare_data/output/pseudobulk")
SAMPLE_COLUMN = "demultiplex_sample"
CELL_TYPE_COLUMN = "multivi_cell_type2"
CELL_TYPES = (
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT",
)


def write_sparse_matrix(path: Path, matrix: sparse.spmatrix) -> None:
    """Write a sparse Matrix Market file through gzip."""
    with gzip.open(path, "wb") as handle:
        mmwrite(handle, matrix)


adata = anndata.read_h5ad(INPUT_H5AD)
missing_obs = {
    SAMPLE_COLUMN, CELL_TYPE_COLUMN, "TSSEnrichment", "NucleosomeRatio"
} - set(adata.obs.columns)
if missing_obs:
    raise ValueError("Missing required obs columns: " + ", ".join(sorted(missing_obs)))

matrix = adata.X
if not sparse.issparse(matrix):
    matrix = sparse.csr_matrix(matrix)
else:
    matrix = matrix.tocsr()

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
peak_ids = pd.Index(adata.var_names.astype(str), name="peak_id")

for cell_type in CELL_TYPES:
    cell_mask = adata.obs[CELL_TYPE_COLUMN].astype(str).to_numpy() == cell_type
    if not cell_mask.any():
        raise ValueError(f"No nuclei found for cell type {cell_type}")

    obs = adata.obs.loc[cell_mask].copy()
    donors = sorted(obs[SAMPLE_COLUMN].astype(str).unique())
    donor_codes = pd.Categorical(
        obs[SAMPLE_COLUMN].astype(str), categories=donors, ordered=True
    ).codes
    grouping = sparse.csr_matrix(
        (np.ones(obs.shape[0], dtype=np.int32),
         (donor_codes, np.arange(obs.shape[0]))),
        shape=(len(donors), obs.shape[0]),
    )
    donor_by_peak = grouping @ matrix[cell_mask, :]

    out_dir = OUTPUT_ROOT / cell_type
    out_dir.mkdir(parents=True, exist_ok=True)
    write_sparse_matrix(out_dir / "counts.mtx.gz", donor_by_peak.T.tocoo())
    pd.DataFrame({"peak_id": peak_ids}).to_csv(
        out_dir / "peaks.tsv", sep="\t", index=False
    )
    pd.Series(donors, name="sample_id").to_csv(
        out_dir / "samples.tsv", sep="\t", index=False, header=False
    )

    obs[SAMPLE_COLUMN] = obs[SAMPLE_COLUMN].astype(str)
    qc = obs.groupby(SAMPLE_COLUMN, sort=True).agg(
        cell_count=(SAMPLE_COLUMN, "size"),
        median_TSSEnrichment=("TSSEnrichment", "median"),
        median_NucleosomeRatio=("NucleosomeRatio", "median"),
    )
    qc.index.name = "sample_id"
    qc.to_csv(out_dir / "qc_summary.tsv", sep="\t")
    print(f"{cell_type}: {len(donors)} donors and {adata.n_vars} peaks")
