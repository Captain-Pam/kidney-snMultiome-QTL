#!/usr/bin/env python
"""Restrict the gene annotation to genes that are eQTL targets.

`hound_snp --stats covgene/logSED` pairs each variant with every gene in the
GTF whose span it falls inside, and scores the variant's effect on that gene's
exonic coverage. Handing it the full GENCODE annotation would produce a score
for every overlapping gene, most of which are not the eGene being tested — so
the GTF is filtered to the eQTL target genes first.

Target genes are named either by symbol or by Ensembl ID in the eQTL table, so
both are matched. Ensembl IDs are matched on the unversioned prefix.

Also writes a gene_id -> gene_name map, because the scores come back keyed by
versioned gene_id while the eQTL table uses symbols.

Usage:
    python 1_prepare_gtf.py --out_gtf <...>.gtf --out_map <...>.tsv
"""
import argparse
import os
import re

import pandas as pd

RE_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
RE_GENE_ID = re.compile(r'gene_id "([^"]+)"')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_gtf", required=True)
    parser.add_argument("--out_map", required=True)
    parser.add_argument("--eqtl_tsv", default=os.environ.get("EQTL_SUSIE_TSV"))
    parser.add_argument("--gencode_gtf", default=os.environ.get("GENCODE_GTF"))
    args = parser.parse_args()

    print(f"Reading target genes from {args.eqtl_tsv}")
    eqtl = pd.read_csv(
        args.eqtl_tsv, sep="\t", compression="gzip", usecols=["phenotype_id"]
    )
    targets = set(eqtl["phenotype_id"].unique())
    print(f"  {len(targets):,} unique target genes")

    ensembl_targets = {g for g in targets if g.startswith("ENSG")}
    symbol_targets = targets - ensembl_targets

    os.makedirs(os.path.dirname(os.path.abspath(args.out_gtf)), exist_ok=True)

    print(f"Filtering {args.gencode_gtf}")
    kept = 0
    gene_id_to_name = {}

    with open(args.gencode_gtf) as fin, open(args.out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            gene_id_match = RE_GENE_ID.search(line)
            gene_name_match = RE_GENE_NAME.search(line)
            gene_id = gene_id_match.group(1) if gene_id_match else None
            gene_name = gene_name_match.group(1) if gene_name_match else None
            gene_id_base = gene_id.split(".")[0] if gene_id else None

            if gene_name in symbol_targets or gene_id_base in ensembl_targets:
                fout.write(line)
                kept += 1
                if gene_id:
                    gene_id_to_name[gene_id] = gene_name or gene_id_base

    print(f"  {kept:,} lines -> {args.out_gtf}")

    mapping = pd.DataFrame(
        sorted(gene_id_to_name.items()), columns=["gene_id", "gene_name"]
    )
    mapping.to_csv(args.out_map, sep="\t", index=False)
    print(f"  {len(mapping):,} gene_id -> gene_name entries -> {args.out_map}")

    missing = targets - set(mapping["gene_name"])
    if missing:
        example = sorted(missing)[:10]
        print(
            f"\n{len(missing):,} eQTL target genes absent from the annotation "
            f"(these cannot be scored): {example}"
            f"{'...' if len(missing) > 10 else ''}"
        )


if __name__ == "__main__":
    main()
