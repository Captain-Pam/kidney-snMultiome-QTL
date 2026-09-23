#!/usr/bin/env python
"""Step 8 - Assemble ordered base, PEER and genotype-PC covariate tables.

Inputs:
  data/metadata/donor_metadata.tsv
  data/genotype/genotype_pca.eigenvec
  prepare_data/output/sample_qc/<cell_type>/{samples.qc.txt,qc_summary.tsv}
  prepare_data/output/peer/<cell_type>_peer<k>.PEER_covariates.txt
Outputs:
  prepare_data/output/covariates/<cell_type>/*.tsv

Run from the caqtl/ directory:
  python prepare_data/8_make_covariates.py
"""

from pathlib import Path
import shutil

import numpy as np
import pandas as pd


CELL_TYPES = (
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT",
)
PEER_COUNTS = (0, 1, 2, 3, 4, 5, 10)
PC_COUNTS = (1, 2, 3, 4, 5, 10, 15)
METADATA_FILE = Path("data/metadata/donor_metadata.tsv")
PCA_FILE = Path("data/genotype/genotype_pca.eigenvec")
SAMPLE_QC_ROOT = Path("prepare_data/output/sample_qc")
PEER_ROOT = Path("prepare_data/output/peer")
OUTPUT_ROOT = Path("prepare_data/output/covariates")


def zscore(values: pd.Series) -> pd.Series:
    """Return the population-standardized values used in the original code."""
    numeric = pd.to_numeric(values)
    return (numeric - numeric.mean()) / numeric.std(ddof=0)


metadata = pd.read_csv(METADATA_FILE, sep="\t", dtype={"sample_id": str})
metadata = metadata.set_index("sample_id")
required_metadata = {"Age", "Gender"}
if not required_metadata.issubset(metadata.columns):
    raise ValueError("donor_metadata.tsv must contain sample_id, Age and Gender")

# PLINK eigenvec format: FID IID PC1 ... PCn, with or without a commented header.
pca_raw = pd.read_csv(PCA_FILE, sep=r"\s+", header=None, comment="#")
pca = pca_raw.set_index(1).drop(columns=0).apply(pd.to_numeric)
pca.index = pca.index.astype(str)
pca.columns = [f"geno_PC{i}" for i in range(1, pca.shape[1] + 1)]
if pca.shape[1] < max(PC_COUNTS):
    raise ValueError("genotype_pca.eigenvec must contain at least 15 PCs")

for cell_type in CELL_TYPES:
    qc_dir = SAMPLE_QC_ROOT / cell_type
    samples = pd.read_csv(
        qc_dir / "samples.qc.txt", header=None, names=["sample_id"], dtype=str
    )["sample_id"].tolist()
    qc = pd.read_csv(qc_dir / "qc_summary.tsv", sep="\t", dtype={"sample_id": str})
    qc = qc.set_index("sample_id").loc[samples]

    base = metadata.loc[samples, ["Gender", "Age"]].copy()
    base["log_count_zscore"] = zscore(np.log10(qc["cell_count"] + 1))
    base["median_TSSEnrichment"] = zscore(qc["median_TSSEnrichment"])
    base["median_NucleosomeRatio"] = zscore(qc["median_NucleosomeRatio"])
    base.index.name = "sample_id"

    out_dir = OUTPUT_ROOT / cell_type
    out_dir.mkdir(parents=True, exist_ok=True)
    base.to_csv(out_dir / "base_peer0_pc0.tsv", sep="\t")

    peer2 = None
    for n_peer in PEER_COUNTS[1:]:
        peer_file = PEER_ROOT / f"{cell_type}_peer{n_peer}.PEER_covariates.txt"
        peer = pd.read_csv(peer_file, sep="\t", index_col=0).T
        peer.index = peer.index.astype(str)
        peer = peer.loc[samples].iloc[:, :n_peer]
        peer.columns = [f"InferredCov{i}" for i in range(1, n_peer + 1)]
        combined = pd.concat([base, peer], axis=1)
        combined.to_csv(out_dir / f"base_peer{n_peer}_pc0.tsv", sep="\t")
        if n_peer == 2:
            peer2 = peer

    if peer2 is None:
        raise ValueError(f"Two-factor PEER table was not created for {cell_type}")
    for n_pc in PC_COUNTS:
        pc_table = pca.loc[samples].iloc[:, :n_pc]
        combined = pd.concat([base, peer2, pc_table], axis=1)
        combined.to_csv(out_dir / f"base_peer2_pc{n_pc}.tsv", sep="\t")

    shutil.copyfile(out_dir / "base_peer2_pc0.tsv", out_dir / "final_peer2_pc0.tsv")
    print(f"{cell_type}: wrote covariates for {len(samples)} donors")
