#!/usr/bin/env python3
"""
Purpose:
    Merge cell-type-specific MACS3 peaks, remove excluded regions, construct a
    nucleus-by-peak matrix and apply the final peak-coverage filters.

Input:
    A SnapATAC2 AnnDataSet containing MACS3 peak calls, final cell-type/state
    annotations and sample identifiers, plus the ENCODE hg38 blacklist BED.

Output:
    A blacklist-filtered consensus peak BED, peak-filtering summaries, the
    final retained peak list and an LZF-compressed filtered H5AD peak matrix.

Example:
    python 2_build_and_filter_peak_matrix.py \
        --dataset data/atac_qc.h5ads \
        --blacklist data/reference/hg38-blacklist.v2.bed \
        --output-dir results/peak_calling
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr
import scipy.sparse as sp
import snapatac2 as snap


CELL_TYPES = (
    "CNT_CD_PC",
    "DCT",
    "DTL_ATL",
    "EC",
    "IC",
    "Immune",
    "PEC",
    "PTS",
    "Podocyte",
    "Stromal",
    "TAL",
    "injPT",
)

PEAK_HALF_WIDTH = 250
CELL_FRACTION = 0.005
MIN_CELLS = 50
SAMPLE_FRACTION = 0.05
MIN_SAMPLES = 5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build and filter the final single-nucleus ATAC peak matrix."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="SnapATAC2 AnnDataSet containing the MACS3 results.",
    )
    parser.add_argument(
        "--blacklist",
        required=True,
        help="ENCODE hg38 blacklist BED file.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for the consensus peaks, summaries and filtered matrix.",
    )
    parser.add_argument(
        "--cell-type-column",
        default="multivi_cell_type2",
        help="obs column containing the final cell-type/state labels.",
    )
    parser.add_argument(
        "--sample-column",
        default="demultiplex_sample",
        help="obs column containing sample identifiers.",
    )
    parser.add_argument(
        "--peak-key",
        default="macs3_cell_types",
        help="uns key containing the cell-type-specific MACS3 results.",
    )
    return parser.parse_args()


def prepare_consensus_peaks(dataset, peak_key, blacklist_path):
    """Merge MACS3 peaks and remove non-autosomal and blacklisted regions."""
    merged = snap.tl.merge_peaks(
        dataset.uns[peak_key],
        chrom_sizes=snap.genome.hg38,
        half_width=PEAK_HALF_WIDTH,
    )

    peaks = merged.with_columns(
        pl.col("Peaks").str.replace(r":.*$", "").alias("Chromosome"),
        pl.col("Peaks").str.extract(r":(\d+)-", 1).cast(pl.Int64).alias("Start"),
        pl.col("Peaks").str.extract(r"-(\d+)$", 1).cast(pl.Int64).alias("End"),
    )
    autosomes = [f"chr{chromosome}" for chromosome in range(1, 23)]
    peaks = peaks.filter(pl.col("Chromosome").is_in(autosomes))

    blacklist = pd.read_csv(
        blacklist_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],
    )
    blacklist.columns = ["Chromosome", "Start", "End"]

    peak_ranges = pr.PyRanges(
        peaks.select(["Chromosome", "Start", "End"]).to_pandas()
    )
    blacklist_ranges = pr.PyRanges(blacklist)
    clean_ranges = pl.from_pandas(
        peak_ranges.overlap(blacklist_ranges, invert=True).as_df()[
            ["Chromosome", "Start", "End"]
        ]
    )

    return peaks.join(
        clean_ranges,
        on=["Chromosome", "Start", "End"],
        how="semi",
    )


def filter_peaks_by_coverage(peak_matrix, sample_column, cell_type_column):
    """Apply the final per-cell-type nucleus and sample coverage criteria."""
    matrix = peak_matrix.X
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()

    # Convert counts to binary detection without densifying the peak matrix.
    detected = matrix.astype(np.int32, copy=True)
    detected.data.fill(1)

    samples = np.asarray(peak_matrix.obs[sample_column]).astype(str)
    cell_types = np.asarray(peak_matrix.obs[cell_type_column]).astype(str)

    # Aggregate non-zero nuclei for every sample-by-cell-type group.
    group_index = pd.MultiIndex.from_arrays([samples, cell_types])
    group_codes, unique_groups = pd.factorize(group_index, sort=True)
    group_selector = sp.csr_matrix(
        (
            np.ones(peak_matrix.n_obs, dtype=np.int32),
            (group_codes, np.arange(peak_matrix.n_obs)),
        ),
        shape=(len(unique_groups), peak_matrix.n_obs),
    )
    cells_per_group = np.bincount(group_codes)
    cells_by_group_peak = group_selector @ detected

    sample_support = cells_by_group_peak.copy()
    sample_support.data.fill(1)
    group_cell_types = unique_groups.get_level_values(1).to_numpy(dtype=str)

    keep_any_cell_type = np.zeros(peak_matrix.n_vars, dtype=bool)
    summary_rows = []

    for cell_type in CELL_TYPES:
        group_rows = np.flatnonzero(group_cell_types == cell_type)
        n_samples = int(group_rows.size)
        n_cells = int(cells_per_group[group_rows].sum())

        cells_with_peak = np.asarray(
            cells_by_group_peak[group_rows, :].sum(axis=0)
        ).ravel()
        samples_with_peak = np.asarray(
            sample_support[group_rows, :].sum(axis=0)
        ).ravel()

        cell_threshold = max(MIN_CELLS, int(np.ceil(CELL_FRACTION * n_cells)))
        sample_threshold = max(
            MIN_SAMPLES,
            int(np.ceil(SAMPLE_FRACTION * n_samples)),
        )

        pass_cells = cells_with_peak >= cell_threshold
        pass_samples = samples_with_peak >= sample_threshold
        keep_cell_type = pass_cells & pass_samples
        keep_any_cell_type |= keep_cell_type

        summary_rows.append(
            {
                "cell_type": cell_type,
                "n_cells": n_cells,
                "n_samples": n_samples,
                "cell_threshold": cell_threshold,
                "sample_threshold": sample_threshold,
                "retained_peaks": int(keep_cell_type.sum()),
            }
        )

    return keep_any_cell_type, pd.DataFrame(summary_rows)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("[1/5] Opening the SnapATAC2 dataset.")
    dataset = snap.read_dataset(args.dataset, mode="r")

    try:
        print("[2/5] Merging peaks and removing excluded regions.")
        consensus_peaks = prepare_consensus_peaks(
            dataset,
            peak_key=args.peak_key,
            blacklist_path=args.blacklist,
        )
        consensus_peaks.select(["Chromosome", "Start", "End"]).write_csv(
            output_dir / "consensus_peaks.hg38.autosomes.blacklist_filtered.bed",
            separator="\t",
            include_header=False,
        )

        print("[3/5] Building the paired-insertion peak matrix.")
        peak_matrix = snap.pp.make_peak_matrix(
            dataset,
            use_rep=consensus_peaks["Peaks"],
            counting_strategy="paired-insertion",
        )
    finally:
        dataset.close()

    print("[4/5] Applying cell-type-specific peak coverage filters.")
    keep_mask, summary = filter_peaks_by_coverage(
        peak_matrix,
        sample_column=args.sample_column,
        cell_type_column=args.cell_type_column,
    )
    summary.to_csv(output_dir / "peak_filter_summary.tsv", sep="\t", index=False)
    pd.Series(
        np.asarray(peak_matrix.var_names)[keep_mask],
        name="peak",
    ).to_csv(output_dir / "retained_peaks.tsv", sep="\t", index=False)

    print("[5/5] Writing the filtered peak matrix.")
    filtered_matrix = peak_matrix[:, keep_mask].copy()
    filtered_matrix.write_h5ad(
        output_dir / "peak_matrix.filtered.h5ad",
        compression="lzf",
    )

    print(
        f"Retained {int(keep_mask.sum()):,} of {peak_matrix.n_vars:,} consensus peaks."
    )


if __name__ == "__main__":
    main()
