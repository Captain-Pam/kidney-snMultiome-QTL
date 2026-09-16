"""Loading fine-mapped eQTL effects joined to Cerberus logSED scores.

Shared by the eQTL concordance figure and the scored-variant table.

Note the cell-type naming: the eQTL tables and the Cerberus targets both spell
proximal tubule "PTS", whereas the RASQUAL/ChromBPNet side spells it "PT".
`figlib.qtl_ct` converts to the key used by the shared colour and label maps.
"""
import glob
import os

import h5py
import numpy as np
import pandas as pd

from figlib import RNA_IDX, cfg, qtl_ct

# Track index per cell type, in the model's own naming.
CT_RNA_IDX = {
    "CNT_CD_PC": RNA_IDX["CNT_CD_PC"], "DCT": RNA_IDX["DCT"],
    "DTL_ATL": RNA_IDX["DTL_ATL"], "EC": RNA_IDX["EC"], "IC": RNA_IDX["IC"],
    "Immune": RNA_IDX["Immune"], "PEC": RNA_IDX["PEC"], "PTS": RNA_IDX["PT"],
    "Podocyte": RNA_IDX["Podocyte"], "Stromal": RNA_IDX["Stromal"],
    "TAL": RNA_IDX["TAL"], "injPT": RNA_IDX["injPT"],
}

CELL_TYPES = list(CT_RNA_IDX)


def scores_dir():
    return os.path.join(cfg("WORK_DIR"), "eqtl", "scores", "logSED", "credible_set")


def load_logsed(directory=None, gene_map_path=None):
    """One row per scored SNP-gene pair, with the full per-track score vector."""
    directory = directory or scores_dir()
    gene_map = pd.read_csv(gene_map_path or cfg("EQTL_GENE_ID_MAP"), sep="\t")
    id_to_name = dict(zip(gene_map["gene_id"], gene_map["gene_name"]))

    chunks = sorted(
        glob.glob(os.path.join(directory, "chunk_*")),
        key=lambda p: int(p.rsplit("_", 1)[1]),
    )
    print(f"Loading {len(chunks)} score chunks from {directory}")

    rows = []
    for chunk in chunks:
        path = os.path.join(chunk, "scores.h5")
        if not os.path.exists(path):
            continue
        try:
            handle = h5py.File(path)
        except OSError:
            print(f"  unreadable, skipping: {path}")
            continue
        with handle as scores_file:
            if scores_file["progress_status"][()].decode() != "completed":
                print(f"  incomplete, skipping: {path}")
                continue
            snps = scores_file["snp"][:].astype(str)
            gene_ids = scores_file["gene_ids"][:].astype(str)
            snp_idx = scores_file["snp_idx"][:]
            gene_idx = scores_file["gene_idx"][:]
            scores = scores_file["covgene/logSED"][:]

        # snp_idx / gene_idx index into the snp and gene_ids tables; each
        # element of `scores` is one SNP-gene pair.
        for pair in range(len(snp_idx)):
            gene_id = gene_ids[gene_idx[pair]]
            rows.append(
                (
                    snps[snp_idx[pair]],
                    id_to_name.get(gene_id, gene_id.split(".")[0]),
                    scores[pair],
                )
            )

    frame = pd.DataFrame(rows, columns=["variant_id", "gene_name", "scores"])
    print(f"  {len(frame):,} SNP-gene pairs")
    return frame


def load_eqtl(eqtl_tsv=None):
    """Credible-set eQTL rows, with the slope aligned to hg38 ALT-vs-REF.

    Model scores are log(ALT/REF). Where the eQTL effect allele is the hg38
    reference (`match_type == "swap"`), the slope is negated so both quantities
    describe the same allele; comparing the raw slope would flip the sign for
    those variants and wash out the correlation.
    """
    eqtl = pd.read_csv(
        eqtl_tsv or cfg("EQTL_SUSIE_TSV"),
        sep="\t",
        compression="gzip",
        usecols=[
            "celltype", "variant_id", "phenotype_id", "A1", "A2",
            "slope", "match_type", "cs", "variable_prob",
        ],
    )
    eqtl = eqtl[eqtl["cs"] != -1].copy()
    eqtl.loc[eqtl["match_type"] == "swap", "slope"] *= -1
    return eqtl


def load_matched(pip_threshold=0.5):
    """eQTL effects joined to the cell-type-matched Cerberus logSED score."""
    pairs = load_logsed()
    eqtl = load_eqtl()

    merged = eqtl.merge(
        pairs,
        left_on=["variant_id", "phenotype_id"],
        right_on=["variant_id", "gene_name"],
        how="inner",
    )
    print(f"  {len(merged):,} matched variant-gene pairs")

    scores = np.vstack(merged["scores"].values)
    track = merged["celltype"].map(CT_RNA_IDX)
    merged["logSED"] = scores[np.arange(len(merged)), track.values].astype(float)
    merged["ct_key"] = merged["celltype"].map(qtl_ct)

    subset = merged[merged["variable_prob"] >= pip_threshold].dropna(
        subset=["logSED", "slope"]
    )
    print(f"  PIP >= {pip_threshold}: {len(subset):,} rows")
    return subset.copy()
