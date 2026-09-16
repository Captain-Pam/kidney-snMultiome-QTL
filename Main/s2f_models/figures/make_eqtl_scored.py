#!/usr/bin/env python
"""Build the scored eQTL table: one row per fine-mapped variant-gene-cell-type triple.

The counterpart of `make_caqtl_scored.py` for expression. Restricted to
confidently fine-mapped credible-set variants (PIP >= 0.5), which is the set the
concordance analysis uses.

| Column | Meaning |
|---|---|
| `variant_id`, `chrom`, `pos`, `ref`, `alt`, `is_snp` | variant identity |
| `gene`, `celltype`, `dist_to_tss` | the eQTL association |

`dist_to_tss` is `pos` minus the gene's 5' end taken from the GTF `gene`
records. The version of this table circulated before the pipeline was cleaned up
used a different, undocumented distance definition, so those values differ for
minus-strand genes; every other column reproduces exactly.
| `match_type` | whether the eQTL effect allele is the hg38 reference |
| `pip` | SuSiE posterior inclusion probability |
| `eqtl_beta` | effect size as reported |
| `eqtl_beta_aligned` | the same, negated for `swap` rows so it is ALT-vs-REF |
| `chrombpnet_logfc` | ChromBPNet accessibility effect, fold ensemble |
| `borzoi_atac_logSUM` | Cerberus accessibility effect, fold ensemble |
| `borzoi_rna_logSED` | Cerberus expression effect over the eGene's exons |

Compare model scores against `eqtl_beta_aligned`, never the raw `eqtl_beta`:
model scores are log(ALT/REF), and for `swap` rows the raw beta describes the
opposite allele.

The accessibility columns are optional — they come from scoring the same eQTL
variants with the accessibility statistics, and are left empty if those score
directories are absent.

Usage:
    python make_eqtl_scored.py [--out <path>] [--pip 0.5]
"""
import argparse
import glob
import os
import re

import h5py
import numpy as np
import pandas as pd

from eqtl_data import load_matched
from figlib import ATAC_IDX, cfg, qtl_ct

_GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def load_tss(gtf_path, genes):
    """Transcription start site per gene, from the GTF `gene` records.

    The TSS is the 5' end: `start` on the plus strand, `end` on the minus
    strand. Genes are named by symbol or by unversioned Ensembl ID depending on
    the eQTL table, so both are matched.

    Scans the GTF once.
    """
    wanted = set(genes)
    spans = {}

    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            name_match = _GENE_NAME_RE.search(fields[8])
            id_match = _GENE_ID_RE.search(fields[8])
            candidates = set()
            if name_match:
                candidates.add(name_match.group(1))
            if id_match:
                candidates.add(id_match.group(1).split(".")[0])

            for key in candidates & wanted:
                spans[key] = (int(fields[3]), int(fields[4]), fields[6])

    missing = wanted - set(spans)
    if missing:
        print(f"  {len(missing)} genes absent from the annotation, dist_to_tss empty")

    return {
        gene: (start if strand == "+" else end)
        for gene, (start, end, strand) in spans.items()
    }


def load_chrombpnet_scores(scores_dir):
    """ChromBPNet logfc per (cell type, variant), from the fold ensemble."""
    scores = {}
    if not os.path.isdir(scores_dir):
        print(f"  ChromBPNet scores not found at {scores_dir}, column left empty")
        return scores

    for path in sorted(glob.glob(os.path.join(scores_dir, "*", "*.variant_scores.tsv"))):
        if path.endswith(".shuffled.tsv"):
            continue
        celltype = os.path.basename(os.path.dirname(path))
        frame = pd.read_csv(path, sep="\t")
        for variant_id, logfc in zip(frame["variant_id"], frame["logfc"]):
            scores[(celltype, variant_id)] = logfc

    print(f"  {len(scores):,} ChromBPNet cell-type x variant scores")
    return scores


def load_cerberus_atac(scores_dir):
    """Cerberus ATAC logSUM per variant, one column per cell type."""
    if not os.path.isdir(scores_dir):
        print(f"  Cerberus ATAC scores not found at {scores_dir}, column left empty")
        return None

    frames, arrays = [], []
    for chunk in sorted(glob.glob(os.path.join(scores_dir, "chunk_*"))):
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
        return None

    scores = pd.concat(frames, ignore_index=True)
    values = np.concatenate(arrays, axis=0)
    scores["_key"] = (
        scores["chrom"] + ":" + scores["pos"].astype(str)
        + ":" + scores["ref"] + ":" + scores["alt"]
    )
    for celltype, index in ATAC_IDX.items():
        scores[celltype] = values[:, index]
    print(f"  {len(scores):,} Cerberus ATAC-scored variants")
    return scores.set_index("_key")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=None)
    parser.add_argument("--pip", type=float, default=0.5)
    args = parser.parse_args()

    out_path = args.out or os.path.join(cfg("FIGURE_DIR"), "eqtl_finemapped_scored.csv")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    work_dir = cfg("WORK_DIR")
    data = load_matched(pip_threshold=args.pip)

    # variant_id is chrom:pos; the alleles come from the eQTL table's own
    # match_type convention, matching cerberus/score/eqtl/2_prepare_vcf.py.
    positions = data["variant_id"].str.split(":", expand=True)
    data["chrom"] = positions[0]
    data["pos"] = positions[1].astype(int)

    print("Loading annotation...")
    tss = load_tss(cfg("EQTL_GTF"), data["phenotype_id"].unique())

    print("Loading model scores...")
    chrombpnet = load_chrombpnet_scores(
        os.path.join(work_dir, "chrombpnet", "eqtl", "scores_ensemble")
    )
    cerberus_atac = load_cerberus_atac(
        os.path.join(work_dir, "eqtl", "scores", "atac_logSUM", "ensemble")
    )

    # The eQTL table records A1 as the effect allele; match_type says whether
    # that is the hg38 reference. See cerberus/score/eqtl/2_prepare_vcf.py.
    is_swap = data["match_type"] == "swap"
    ref = np.where(is_swap, data["A1"], data["A2"])
    alt = np.where(is_swap, data["A2"], data["A1"])

    table = pd.DataFrame(
        {
            "variant_id": (
                data["chrom"] + ":" + data["pos"].astype(str)
                + ":" + ref + ":" + alt
            ),
            "chrom": data["chrom"].values,
            "pos": data["pos"].values,
            "ref": ref,
            "alt": alt,
            "is_snp": (pd.Series(ref).str.len() == 1).values
            & (pd.Series(alt).str.len() == 1).values,
            "gene": data["phenotype_id"].values,
            "celltype": data["celltype"].values,
            "dist_to_tss": [
                (position - tss[gene]) if gene in tss else pd.NA
                for position, gene in zip(data["pos"], data["phenotype_id"])
            ],
            "match_type": data["match_type"].values,
            "pip": data["variable_prob"].values,
            # load_eqtl already negated the slope for swap rows, so what it
            # holds is the aligned quantity; undo it to report both.
            "eqtl_beta_aligned": data["slope"].values,
            "borzoi_rna_logSED": data["logSED"].values,
        }
    )
    table["eqtl_beta"] = np.where(
        is_swap.values, -table["eqtl_beta_aligned"], table["eqtl_beta_aligned"]
    )

    table["chrombpnet_logfc"] = [
        chrombpnet.get((celltype, variant_id), np.nan)
        for celltype, variant_id in zip(table["celltype"], table["variant_id"])
    ]

    if cerberus_atac is not None:
        table["borzoi_atac_logSUM"] = [
            cerberus_atac[qtl_ct(celltype)].get(variant_id, np.nan)
            if qtl_ct(celltype) in cerberus_atac.columns
            else np.nan
            for celltype, variant_id in zip(table["celltype"], table["variant_id"])
        ]
    else:
        table["borzoi_atac_logSUM"] = np.nan

    columns = [
        "variant_id", "chrom", "pos", "ref", "alt", "is_snp", "gene", "celltype",
        "dist_to_tss", "match_type", "pip", "eqtl_beta", "eqtl_beta_aligned",
        "chrombpnet_logfc", "borzoi_atac_logSUM", "borzoi_rna_logSED",
    ]
    table = table[columns]

    table.to_csv(out_path, index=False)
    print(
        f"\n{len(table):,} rows, {table['variant_id'].nunique():,} variants, "
        f"{table['gene'].nunique():,} genes -> {out_path}"
    )


if __name__ == "__main__":
    main()
