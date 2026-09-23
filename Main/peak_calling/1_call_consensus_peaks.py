#!/usr/bin/env python3
"""
Purpose:
    Call MACS3 peaks for the 12 final kidney cell types and states.

Input:
    A writable SnapATAC2 AnnDataSet containing fragment data and final
    cell-type/state annotations in obs.

Output:
    Cell-type-specific MACS3 results stored in the input AnnDataSet under the
    specified uns key.

Example:
    python 1_call_consensus_peaks.py \
        --dataset data/atac_qc.h5ads \
        --cell-type-column multivi_cell_type2 \
        --n-jobs 12
"""

import argparse

import numpy as np
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Call MACS3 peaks for the 12 final kidney cell types and states."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Writable SnapATAC2 AnnDataSet (.h5ads).",
    )
    parser.add_argument(
        "--cell-type-column",
        default="multivi_cell_type2",
        help="obs column containing the final cell-type/state labels.",
    )
    parser.add_argument(
        "--key-added",
        default="macs3_cell_types",
        help="uns key used to store the MACS3 results.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=12,
        help="Number of parallel MACS3 jobs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("[1/2] Opening the annotated SnapATAC2 dataset.")
    dataset = snap.read_dataset(args.dataset, mode="r+")

    try:
        # Confirm that all final analysis groups are represented.
        observed = set(np.asarray(dataset.obs[args.cell_type_column]).astype(str))
        missing = [cell_type for cell_type in CELL_TYPES if cell_type not in observed]
        if missing:
            raise ValueError(
                "Missing final cell-type/state labels: " + ", ".join(missing)
            )

        # Call peaks directly from the final 12 cell-type/state annotations.
        print("[2/2] Calling MACS3 peaks for the 12 cell types and states.")
        snap.tl.macs3(
            dataset,
            groupby=args.cell_type_column,
            selections=set(CELL_TYPES),
            n_jobs=args.n_jobs,
            key_added=args.key_added,
        )
    finally:
        dataset.close()

    print(f"MACS3 peak calls were stored in dataset.uns['{args.key_added}'].")


if __name__ == "__main__":
    main()
