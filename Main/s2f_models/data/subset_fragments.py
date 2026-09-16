#!/usr/bin/env python
"""Subset a 10x ATAC fragment file to a set of cell barcodes.

Usage: python subset_fragments.py <fragments.tsv.gz> <barcodes.csv> <out.tsv>

<barcodes.csv> is one barcode per line, no header. Comment lines in the
fragment file are preserved so the output stays a valid fragment file.
"""
import gzip
import sys


def subset_fragments(fragments_file, barcodes_file, output_file):
    with open(barcodes_file) as f:
        barcodes = {line.strip() for line in f if line.strip()}
    print(f"Loaded {len(barcodes)} barcodes")

    count = 0
    with gzip.open(fragments_file, "rt") as fin, open(output_file, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] in barcodes:
                fout.write(line)
                count += 1

    print(f"Wrote {count} fragments to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    subset_fragments(sys.argv[1], sys.argv[2], sys.argv[3])
