#!/usr/bin/env python
"""Convert ChromBPNet contribution scores to the npz pair TF-MoDISco expects.

`chrombpnet contribs_bw` writes an HDF5 file holding, for each scored region:

    raw/seq             (n, 4, L) int8      one-hot encoded sequence
    projected_shap/seq  (n, 4, L) float16   hypothetical contribution scores

Both are already length-last, which is the layout `modisco motifs` wants; this
just splits them into two float32 npz files.

Usage:
    python h5_to_npz.py --h5 <ct>.profile_scores.h5 --ohe ohe.npz --shap shap.npz
"""
import argparse

import h5py
import hdf5plugin  # noqa: F401  — registers the BLOSC filter chrombpnet writes with
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", required=True, help="chrombpnet contribution scores")
    parser.add_argument("--ohe", required=True, help="output npz, one-hot sequences")
    parser.add_argument("--shap", required=True, help="output npz, attributions")
    args = parser.parse_args()

    with h5py.File(args.h5, "r") as handle:
        ohe = handle["raw/seq"][:]
        shap = handle["projected_shap/seq"][:]

    print(f"sequences:    {ohe.shape} {ohe.dtype}")
    print(f"attributions: {shap.shape} {shap.dtype}")

    np.savez_compressed(args.ohe, ohe.astype(np.float32))
    np.savez_compressed(args.shap, shap.astype(np.float32))
    print(f"wrote {args.ohe} and {args.shap}")


if __name__ == "__main__":
    main()
