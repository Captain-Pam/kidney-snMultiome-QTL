#!/usr/bin/env python
"""Build baskerville targets files from a set of per-cell-type .w5 coverage files.

A targets file is the table that tells `hound_data` which coverage tracks the
model predicts, how to summarise them into bins, and which pairs of tracks are
the two strands of the same RNA assay.

Each dataset directory is expected to hold files named

    <cell_type>_atac.w5      chromatin accessibility
    <cell_type>_rna+.w5      forward-strand RNA
    <cell_type>_rna-.w5      reverse-strand RNA

Usage
-----
    python 3_make_targets.py --manifest manifest.tsv --out targets/kidney_targets_w5_human.txt

`manifest.tsv` is a headerless, tab-separated list of the datasets to include,
in the order their tracks should appear:

    prefix       dataset_label             w5_directory
    Susztak      Susztak_CKD_multiome      /path/to/susztak/bw_w5
    PRJNA909297  Ledru_Control_multiome    /path/to/PRJNA909297/bw_w5

Pass --local to additionally write the `_local` variant, which carries an extra
`window` column set to "local". That column is what makes `hound_snp
--local_window` restrict its statistic to a window around the variant; without
it the flag has no effect. See cerberus/score/caqtl/.
"""
import argparse
import os

import pandas as pd

# Coverage is square-root-summed within each bin and clipped, separately per
# assay, to keep a few very high-coverage loci from dominating the loss.
CLIP = {"atac": 200, "rna+": 300, "rna-": 300}
SUM_STAT = "sum_sqrt"
SCALE = 1

ASSAY_LABEL = {"atac": "ATAC", "rna+": "RNA+", "rna-": "RNA-"}


def track_assay(filename):
    """Return 'atac', 'rna+', 'rna-', or None for a .w5 filename."""
    stem = filename[: -len(".w5")]
    for assay in ("atac", "rna+", "rna-"):
        if stem.endswith("_" + assay):
            return assay
    return None


def cell_type_of(filename, assay):
    """Strip the trailing _<assay> and any _merged infix from a .w5 filename."""
    stem = filename[: -len(".w5")]
    stem = stem[: -(len(assay) + 1)]
    return stem.replace("_merged", "")


def build_rows(manifest):
    rows = []
    for prefix, label, w5_dir in manifest:
        if not os.path.isdir(w5_dir):
            raise FileNotFoundError(f"{label}: no such directory: {w5_dir}")

        # Sorted order puts <ct>_atac, <ct>_rna+, <ct>_rna- adjacent and in that
        # order, which the strand_pair arithmetic below relies on.
        files = sorted(f for f in os.listdir(w5_dir) if f.endswith(".w5"))
        if not files:
            raise FileNotFoundError(f"{label}: no .w5 files in {w5_dir}")

        for filename in files:
            assay = track_assay(filename)
            if assay is None:
                continue

            index = len(rows)
            if assay == "rna+":
                strand_pair = index + 1  # the matching rna- track
            elif assay == "rna-":
                strand_pair = index - 1  # the matching rna+ track
            else:
                strand_pair = index  # ATAC is unstranded, so it pairs with itself

            cell_type = cell_type_of(filename, assay)
            identifier = f"{prefix}_{filename[:-len('.w5')].upper().replace('_', '')}"

            rows.append(
                {
                    "identifier": identifier,
                    "file": os.path.join(w5_dir, filename),
                    "clip": CLIP[assay],
                    "scale": SCALE,
                    "sum_stat": SUM_STAT,
                    "strand_pair": strand_pair,
                    "description": f"{ASSAY_LABEL[assay]}:{cell_type}:{label}",
                }
            )

        print(f"{label}: {len(files)} tracks from {w5_dir}")

    return pd.DataFrame(rows)


def read_manifest(path):
    manifest = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            prefix, label, w5_dir = line.split("\t")
            manifest.append((prefix, label, w5_dir))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="datasets to include")
    parser.add_argument("--out", required=True, help="output targets file")
    parser.add_argument(
        "--local",
        action="store_true",
        help="also write a _local variant carrying window=local",
    )
    args = parser.parse_args()

    targets = build_rows(read_manifest(args.manifest))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    targets.to_csv(args.out, sep="\t", index=True)
    print(f"\n{len(targets)} targets -> {args.out}")

    if args.local:
        local = targets.copy()
        local["window"] = "local"
        base, ext = os.path.splitext(args.out)
        local_path = f"{base}_local{ext}"
        local.to_csv(local_path, sep="\t", index=True)
        print(f"{len(local)} targets -> {local_path}")


if __name__ == "__main__":
    main()
