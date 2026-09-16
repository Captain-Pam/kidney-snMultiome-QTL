#!/usr/bin/env python
"""Average Cerberus scores across model folds.

The eight folds differ only in which sequences were held out during training.
Averaging their scores is the ensembling step, and every published score is a
fold average — single-fold scores are noticeably noisier.

Works on any hound_snp / hound_ism_snp output: name the datasets to average
with --keys and everything else is copied through from fold 0.

    variant scores    --keys cov/logSUM
    gene scores       --keys covgene/logSED
    ISM               --keys ref/cov/logSUM alt/cov/logSUM

Usage
-----
    # 63 chunks of caQTL scores
    python ensemble_folds.py \
        --fold_dir '<work>/caqtl/scores/logSUM/local_f{fold}c0' \
        --out_dir  '<work>/caqtl/scores/logSUM/local_ensemble' \
        --keys cov/logSUM --n_chunks 63

    # a single unchunked ISM run
    python ensemble_folds.py \
        --fold_dir '<work>/ism/f{fold}c0' --out_dir '<work>/ism/ensemble' \
        --keys ref/cov/logSUM alt/cov/logSUM

`--fold_dir` must contain the literal placeholder `{fold}`.
"""
import argparse
import os

import h5py
import numpy as np


def is_complete(path):
    """hound_snp marks finished files; a crashed job leaves a partial one."""
    if not os.path.exists(path):
        return False
    try:
        with h5py.File(path, "r") as handle:
            if "progress_status" not in handle:
                return True  # not all writers record status
            return handle["progress_status"][()].decode() == "completed"
    except OSError:
        return False


def fold_paths(fold_dir, n_folds, chunk):
    for fold in range(n_folds):
        directory = fold_dir.format(fold=fold)
        if chunk is not None:
            directory = os.path.join(directory, f"chunk_{chunk}")
        yield os.path.join(directory, "scores.h5")


def ensemble_one(fold_dir, out_dir, keys, n_folds, chunk):
    if chunk is not None:
        out_dir = os.path.join(out_dir, f"chunk_{chunk}")
    out_path = os.path.join(out_dir, "scores.h5")

    if is_complete(out_path):
        return "skip"

    paths = list(fold_paths(fold_dir, n_folds, chunk))
    for path in paths:
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        if not is_complete(path):
            raise RuntimeError(f"incomplete, rerun this fold: {path}")

    # Accumulate in float32; the files themselves are float16, and summing eight
    # of those loses precision that matters for the smaller effect sizes.
    totals = {key: None for key in keys}

    for path in paths:
        with h5py.File(path, "r") as handle:
            for key in keys:
                values = handle[key][:].astype(np.float32)
                totals[key] = values if totals[key] is None else totals[key] + values

    os.makedirs(out_dir, exist_ok=True)
    keys = set(keys)
    with h5py.File(paths[0], "r") as src, h5py.File(out_path, "w") as dst:
        # Copy every dataset except the ones being averaged, then write the
        # averages in their place. visititems walks nested groups, so this keeps
        # siblings like cov/logSUM_quantiles alongside an averaged cov/logSUM.
        def copy_unless_averaged(name, obj):
            if isinstance(obj, h5py.Dataset) and name not in keys:
                dst.create_dataset(name, data=obj[:])

        src.visititems(copy_unless_averaged)

        for key in keys:
            dst.create_dataset(
                key,
                data=(totals[key] / n_folds).astype(np.float16),
                compression="gzip",
                compression_opts=4,
            )

    return "done"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold_dir", required=True, help="path containing {fold}")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--keys", nargs="+", required=True, help="datasets to average")
    parser.add_argument(
        "--n_folds", type=int, default=int(os.environ.get("CERBERUS_N_FOLDS", 8))
    )
    parser.add_argument(
        "--n_chunks", type=int, default=None, help="omit for an unchunked run"
    )
    args = parser.parse_args()

    if "{fold}" not in args.fold_dir:
        raise SystemExit("--fold_dir must contain the placeholder {fold}")

    chunks = [None] if args.n_chunks is None else list(range(args.n_chunks))
    counts = {"done": 0, "skip": 0, "error": 0}

    for chunk in chunks:
        label = "single" if chunk is None else f"chunk_{chunk}"
        try:
            status = ensemble_one(
                args.fold_dir, args.out_dir, args.keys, args.n_folds, chunk
            )
            counts[status] += 1
            if status == "done":
                print(f"  {label}: ensembled")
        except (FileNotFoundError, RuntimeError) as exc:
            counts["error"] += 1
            print(f"  {label}: ERROR {exc}")

    print(
        f"\n{counts['done']} ensembled, {counts['skip']} already done, "
        f"{counts['error']} errors -> {args.out_dir}"
    )


if __name__ == "__main__":
    main()
