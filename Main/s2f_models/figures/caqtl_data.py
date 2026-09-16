"""Loading the joined caQTL score table.

Shared by the two caQTL concordance figures and the scored-variant table, so
they cannot drift apart in which scores they use.
"""
import glob
import os

import h5py
import numpy as np
import pandas as pd

from figlib import ATAC_IDX, cfg

CELL_TYPES = list(ATAC_IDX)


def cerberus_scores_dir():
    return os.path.join(cfg("WORK_DIR"), "caqtl", "scores", "logSUM", "local_ensemble")


def merged_table_path():
    return os.path.join(
        cfg("WORK_DIR"), "chrombpnet", "caqtl", "caqtl_scores_merged.tsv.gz"
    )


def variant_key(chrom, pos, ref, alt):
    return chrom + "_" + pos.astype(str) + "_" + ref + "_" + alt


def load_cerberus_scores(scores_dir=None):
    """Per-variant ATAC logSUM from the fold ensemble, one column per cell type.

    Returned indexed by `<chrom>_<pos>_<ref>_<alt>`.
    """
    scores_dir = scores_dir or cerberus_scores_dir()
    frames, arrays = [], []

    chunks = sorted(
        glob.glob(os.path.join(scores_dir, "chunk_*")),
        key=lambda p: int(p.split("_")[-1]),
    )
    for chunk in chunks:
        path = os.path.join(chunk, "scores.h5")
        if not os.path.exists(path):
            continue
        with h5py.File(path) as handle:
            frames.append(
                pd.DataFrame(
                    {
                        "chrom": handle["chr"][:].astype(str),
                        "pos": handle["pos"][:],
                        "ref": handle["ref_allele"][:].astype(str),
                        "alt": handle["alt_allele"][:].astype(str),
                    }
                )
            )
            arrays.append(handle["cov/logSUM"][:].astype(np.float32))

    if not frames:
        raise SystemExit(
            f"No scored chunks under {scores_dir}. Run the Cerberus caQTL "
            "scoring and fold ensembling first."
        )

    scores = pd.concat(frames, ignore_index=True)
    values = np.concatenate(arrays, axis=0)
    scores["_key"] = variant_key(
        scores["chrom"], scores["pos"], scores["ref"], scores["alt"]
    )
    for celltype, index in ATAC_IDX.items():
        scores[f"logSUM_{celltype}"] = values[:, index]
    return scores.set_index("_key")


def load_merged(scores_dir=None, merged_tsv=None):
    """RASQUAL + ChromBPNet table with current Cerberus logSUM columns joined on.

    The merged table may already carry `logSUM_*` columns from an earlier
    scoring run; they are overwritten here so that both models' scores always
    come from the current fold ensemble.
    """
    merged = pd.read_csv(merged_tsv or merged_table_path(), sep="\t")
    merged["_key"] = variant_key(
        merged["Chromosome"],
        merged["SNP_position"],
        merged["Ref_allele"],
        merged["Alt_allele"],
    )

    cerberus = load_cerberus_scores(scores_dir)
    for celltype in CELL_TYPES:
        merged[f"logSUM_{celltype}"] = merged["_key"].map(
            cerberus[f"logSUM_{celltype}"]
        )
    return merged
