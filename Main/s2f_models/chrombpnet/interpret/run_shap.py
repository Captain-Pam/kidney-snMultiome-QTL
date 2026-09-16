#!/usr/bin/env python
"""Per-variant DeepLIFT/SHAP contribution scores from ChromBPNet.

For one cell type, computes contribution scores on both the reference and the
alternate allele sequence of every variant in a list, averaged over folds, and
saves a window around the variant.

Contributions are `shap * one_hot`, i.e. the attribution of the base actually
present, which is what an attribution logo plots. Attribution is taken against
the **counts** head by default: the profile head explains where reads fall
within the window, whereas a caQTL is a change in how many reads fall there.

Output: <out_dir>/<cell_type>_ism.h5
    variant_ids   (N,)         variant identifiers, as byte strings
    ref_contrib   (N, W, 4)    contribution scores, reference allele
    alt_contrib   (N, W, 4)    contribution scores, alternate allele
    ref_onehot    (N, W, 4)    one-hot reference sequence
    alt_onehot    (N, W, 4)    one-hot alternate sequence
  with W = 2 * half_win, and the variant at index half_win.

Usage:
    python run_shap.py --celltype PTS --variants variants.tsv --out_dir <dir>

`--variants` is the five-column variant-scorer table written by
../score/1_prepare_caqtl_variants.sh: chr, pos, variant_id, ref, alt.
"""
import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd

# ChromBPNet models are TensorFlow 1.x graphs; SHAP needs v1 behaviour.
sys.path.insert(0, os.environ["VARIANT_SCORER_SRC"])
import tensorflow as tf  # noqa: E402

tf.compat.v1.disable_v2_behavior()

from utils.helpers import load_model_wrapper  # noqa: E402
from utils.shap_utils import fetch_shap  # noqa: E402

# ChromBPNet's receptive field. The variant sits at the centre of the input.
INPUT_LEN = 2114
VAR_IDX = INPUT_LEN // 2

# Cell-type naming differs between the QTL tables and the model directories.
QTL_TO_MODEL = {"PT": "PTS"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--celltype", required=True, help="cell type, QTL naming")
    parser.add_argument("--variants", required=True, help="variant table (5 columns)")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument(
        "--half_win", type=int, default=25, help="bases saved either side of the variant"
    )
    parser.add_argument("--shap_type", default="counts", choices=["counts", "profile"])
    parser.add_argument(
        "--n_folds", type=int, default=int(os.environ.get("CHROMBPNET_N_FOLDS", 5))
    )
    args = parser.parse_args()

    model_ct = QTL_TO_MODEL.get(args.celltype, args.celltype)
    models_base = os.environ["CHROMBPNET_MODEL_DIR"]
    genome = os.environ["CHROMBPNET_FASTA"]

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"{args.celltype}_ism.h5")
    if os.path.exists(out_path):
        print(f"exists, skipping: {out_path}")
        return

    variants = pd.read_csv(
        args.variants,
        sep="\t",
        header=None,
        names=["chr", "pos", "variant_id", "allele1", "allele2"],
    )
    variants = variants.drop_duplicates(subset=["chr", "pos", "allele1", "allele2"])
    variants = variants.reset_index(drop=True)
    print(f"{args.celltype}: {len(variants)} variants")

    ref_sum = alt_sum = None
    ref_onehot = alt_onehot = None
    variant_ids = None
    folds_done = 0

    for fold in range(args.n_folds):
        model_path = os.path.join(
            models_base, f"{model_ct}_fold_{fold}", "models", "chrombpnet_nobias.h5"
        )
        if not os.path.exists(model_path):
            print(f"  fold {fold}: no model, skipping")
            continue

        print(f"  fold {fold}: {model_path}")
        model = load_model_wrapper(model_path)

        variant_ids, ref_inp, alt_inp, ref_shap, alt_shap = fetch_shap(
            model, variants, INPUT_LEN, genome, batch_size=64, shap_type=args.shap_type
        )

        ref_contrib = ref_shap * ref_inp
        alt_contrib = alt_shap * alt_inp

        if ref_sum is None:
            ref_sum, alt_sum = ref_contrib.copy(), alt_contrib.copy()
            # The one-hot sequences are the same for every fold.
            ref_onehot, alt_onehot = ref_inp.copy(), alt_inp.copy()
        else:
            ref_sum += ref_contrib
            alt_sum += alt_contrib
        folds_done += 1

        del model
        tf.keras.backend.clear_session()

    if folds_done == 0:
        sys.exit(f"No fold models found under {models_base} for {model_ct}")
    if folds_done < args.n_folds:
        print(f"  warning: averaging over {folds_done}/{args.n_folds} folds")

    ref_avg = ref_sum / folds_done
    alt_avg = alt_sum / folds_done

    start = VAR_IDX - args.half_win
    end = VAR_IDX + args.half_win

    with h5py.File(out_path, "w") as handle:
        handle.create_dataset("variant_ids", data=np.array(variant_ids, dtype="S50"))
        handle.create_dataset("ref_contrib", data=ref_avg[:, start:end, :].astype(np.float32))
        handle.create_dataset("alt_contrib", data=alt_avg[:, start:end, :].astype(np.float32))
        handle.create_dataset("ref_onehot", data=ref_onehot[:, start:end, :].astype(np.float32))
        handle.create_dataset("alt_onehot", data=alt_onehot[:, start:end, :].astype(np.float32))
        handle.attrs["celltype"] = args.celltype
        handle.attrs["shap_type"] = args.shap_type
        handle.attrs["half_win"] = args.half_win
        handle.attrs["n_folds"] = folds_done
        handle.attrs["n_variants"] = len(variant_ids)

    print(f"wrote {out_path} ({len(variant_ids)} variants, +/-{args.half_win} bp)")


if __name__ == "__main__":
    main()
