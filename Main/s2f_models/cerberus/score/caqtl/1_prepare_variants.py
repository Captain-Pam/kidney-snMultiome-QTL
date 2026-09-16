#!/usr/bin/env python
"""Collect RASQUAL caQTL variants into a single VCF for scoring.

Pools the per-cell-type FDR 0.1 caQTL tables, keeps variants whose scoring
window fully contains their peak, deduplicates and sorts. The resulting VCF is
scored by both models — Cerberus here and ChromBPNet via
../../../chrombpnet/score — so the two score sets are directly comparable.

The window filter matters. Cerberus scores a caQTL with `logSUM` over a 1,024 bp
window centred on the variant; if that window does not cover the whole peak, the
score measures something different from RASQUAL's peak-level effect and is not
comparable across variants. Variants failing this are dropped rather than
scored and filtered later.

Usage:
    python 1_prepare_variants.py --out <work>/caqtl/caqtl_full_fdr0.1_inPeaks.vcf
"""
import argparse
import os

import pandas as pd

CELL_TYPES = [
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT",
]

CHR_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHR_ORDER["chrX"] = 23
CHR_ORDER["chrY"] = 24

# Half of hound_snp's --local_window 1024.
LOCAL_HALF = 512


def peak_within_window(row):
    """True if the +/-512 bp scoring window fully contains the variant's peak.

    Feature_ID encodes the peak as chr-start-end, e.g. chr1-100461142-100461644.
    """
    parts = row["Feature_ID"].split("-")
    peak_start = int(parts[-2])
    peak_end = int(parts[-1])
    pos = int(row["SNP_position"]) - 1  # VCF is 1-based, peaks are 0-based
    return (pos - LOCAL_HALF <= peak_start) and (peak_end <= pos + LOCAL_HALF)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="output VCF")
    parser.add_argument(
        "--caqtl_dir",
        default=os.environ.get("RASQUAL_CAQTL_DIR"),
        help="directory of <ct>_full_caPeaks_FDR0.1.tsv files",
    )
    args = parser.parse_args()

    frames = []
    for celltype in CELL_TYPES:
        path = os.path.join(args.caqtl_dir, f"{celltype}_full_caPeaks_FDR0.1.tsv")
        if not os.path.exists(path):
            print(f"  warning: missing {path}")
            continue
        frame = pd.read_csv(path, sep="\t")
        frame["cell_type"] = celltype
        frames.append(frame)
        print(f"  {celltype}: {len(frame):,} rows")

    if not frames:
        raise SystemExit(f"No caQTL tables found under {args.caqtl_dir}")

    variants = pd.concat(frames, ignore_index=True)
    print(f"total: {len(variants):,} rows")

    variants["Chromosome"] = variants["Chromosome"].apply(
        lambda value: value if str(value).startswith("chr") else f"chr{value}"
    )
    variants = variants[variants["Chromosome"].isin(CHR_ORDER)]
    print(f"  {len(variants):,} on canonical chromosomes")

    keep = variants.apply(peak_within_window, axis=1)
    print(f"  {(~keep).sum():,} dropped: scoring window does not cover the peak")
    variants = variants[keep].copy()

    unique = variants.drop_duplicates(
        subset=["Chromosome", "SNP_position", "Ref_allele", "Alt_allele"]
    ).copy()
    print(f"  {len(unique):,} unique variants")

    unique["_order"] = unique["Chromosome"].map(CHR_ORDER)
    unique = unique.sort_values(["_order", "SNP_position"])

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=1_prepare_variants.py\n")
        handle.write(
            "##INFO=<ID=PEAK_FEATURE,Number=1,Type=String,"
            'Description="RASQUAL peak feature ID (chr-start-end)">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in unique.iterrows():
            handle.write(
                f"{row['Chromosome']}\t{int(row['SNP_position'])}\t{row['rs_ID']}\t"
                f"{row['Ref_allele']}\t{row['Alt_allele']}\t.\tPASS\t"
                f"PEAK_FEATURE={row['Feature_ID']}\n"
            )

    print(f"\n{len(unique):,} variants -> {args.out}")
    print(f"For 1000-variant chunks, score with --array=0-{len(unique) // 1000}")


if __name__ == "__main__":
    main()
