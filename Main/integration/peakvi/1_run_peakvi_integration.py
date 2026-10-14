#!/usr/bin/env python3
"""Prepare single-cell chromatin accessibility data for PeakVI integration.

Purpose:
    Validate and preprocess a post-QC chromatin accessibility AnnData object
    for a reproducible PeakVI integration workflow.

Input:
    An H5AD file with raw non-negative integer peak counts in ``adata.X`` and
    a ``Sample`` column in ``adata.obs``.

Preprocessing:
    Validate the complete raw matrix before optional barcode exclusion,
    retain peaks detected in at least ``max(1, int(n_obs * 0.01))`` cells,
    and then binarize positive counts.

Fixed configuration:
    PeakVI uses ``batch_key="Sample"`` and ``model_depth=True``. Hidden and
    latent dimensions are automatically derived by PeakVI. Training uses at
    most 350 epochs, learning rate ``1e-4``, batch size 1024,
    ``early_stopping=True``, patience 50, ``check_val_every_n_epoch=1``, and
    seed 0.

Output:
    A ``peakvi_model`` directory, ``peakvi_integrated.h5ad``,
    ``latent.npy``, ``training_history.png``, ``peakvi_leiden.png``, and
    ``peakvi_workflow.log`` in the requested output directory.

Dependencies:
    Python 3.12, NumPy, AnnData, Scanpy, scvi-tools, PyTorch, and Matplotlib.
    SciPy provides the sparse matrix backend.

"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys


LOGGER_NAME = "peakvi_workflow"
REQUIRED_OBS_COLUMNS = ("Sample",)
OUTPUT_NAMES = (
    "peakvi_model",
    "peakvi_integrated.h5ad",
    "latent.npy",
    "training_history.png",
    "peakvi_leiden.png",
    "peakvi_workflow.log",
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the two required workflow paths."""
    parser = argparse.ArgumentParser(
        description="Integrate post-QC chromatin accessibility data with PeakVI."
    )
    parser.add_argument(
        "--input-h5ad",
        required=True,
        type=Path,
        help="Post-QC chromatin accessibility H5AD with raw counts in .X.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for the model, integrated H5AD, plot, and log.",
    )
    return parser.parse_args(argv)


def prepare_output_directory(output_dir: Path) -> Path:
    """Create the output directory and refuse to overwrite workflow outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    existing_targets = [
        output_dir / output_name
        for output_name in OUTPUT_NAMES
        if (output_dir / output_name).exists()
    ]
    if existing_targets:
        existing_text = ", ".join(str(target) for target in existing_targets)
        raise FileExistsError(
            f"Refusing to overwrite existing PeakVI outputs: {existing_text}"
        )
    return output_dir


def configure_logging(output_dir: Path) -> logging.Logger:
    """Log each workflow message once to standard output and a file."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    # Reconfiguration keeps repeated calls from duplicating handlers.
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    file_handler = logging.FileHandler(output_dir / "peakvi_workflow.log")
    file_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    return logger


def log_software_versions(logger: logging.Logger) -> None:
    """Record the runtime versions needed to interpret the analysis."""
    import anndata
    import scanpy
    import scvi

    logger.info("Python version: %s", sys.version.split()[0])
    logger.info("scvi-tools version: %s", scvi.__version__)
    logger.info("Scanpy version: %s", scanpy.__version__)
    logger.info("AnnData version: %s", anndata.__version__)


def validate_input(adata) -> None:
    """Validate matrix dimensions, required metadata, and raw peak counts."""
    import numpy as np
    from scipy import sparse

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("adata must contain at least one cell and one peak.")

    missing_columns = [
        column for column in REQUIRED_OBS_COLUMNS if column not in adata.obs
    ]
    if missing_columns:
        raise ValueError(
            "Missing required adata.obs columns: " + ", ".join(missing_columns)
        )

    # Sparse matrices store only explicit values; implicit entries are zero.
    count_values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    try:
        valid_counts = (
            np.all(np.isfinite(count_values))
            and np.all(count_values >= 0)
            and np.all(np.equal(count_values, np.floor(count_values)))
        )
    except (TypeError, ValueError):
        valid_counts = False
    if not valid_counts:
        raise ValueError(
            "adata.X must contain raw non-negative integer peak counts."
        )


def prepare_anndata(adata):
    """Filter peaks detected in at least 1% of cells and binarize counts."""
    import numpy as np
    import scanpy as sc

    if adata.n_obs == 0:
        raise ValueError(
            "adata must contain at least one cell after optional barcode exclusion."
        )

    min_cells = max(1, int(adata.n_obs * 0.01))
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if adata.n_vars == 0:
        raise ValueError("No peaks remain after prevalence filtering.")

    # Preserve sparse storage while converting every positive count to one.
    adata.X = (adata.X > 0).astype(np.int8)
    if hasattr(adata.X, "eliminate_zeros"):
        adata.X.eliminate_zeros()
    return adata


def plot_training_history(model, output_path: Path) -> None:
    """Plot training and validation reconstruction losses."""
    import matplotlib.pyplot as plt

    train_key = "reconstruction_loss_train"
    validation_key = "reconstruction_loss_validation"
    train_loss = model.history[train_key][train_key]
    validation_loss = model.history[validation_key][validation_key]

    figure, axis = plt.subplots(figsize=(6, 4))
    axis.plot(train_loss, label="train")
    axis.plot(validation_loss, label="validation")
    axis.axhline(validation_loss.min(), color="black")
    axis.set_xlabel("Epoch")
    axis.set_ylabel("Reconstruction loss")
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_peakvi(adata, output_dir: Path):
    """Train PeakVI and write the fixed latent-space workflow outputs."""
    import matplotlib.pyplot as plt
    import numpy as np
    import scanpy as sc
    import scvi
    import torch

    output_dir = Path(output_dir)
    logger = logging.getLogger(LOGGER_NAME)

    scvi.settings.seed = 0
    torch.set_float32_matmul_precision("high")

    logger.info("Registering AnnData with Sample as the PeakVI batch key.")
    scvi.model.PEAKVI.setup_anndata(adata, batch_key="Sample")
    model = scvi.model.PEAKVI(
        adata,
        model_depth=True,  # Peak depth is modeled internally by PeakVI.
    )

    logger.info("Training PeakVI model.")
    model.train(
        max_epochs=350,
        plan_kwargs={"lr": 1e-4},
        early_stopping=True,
        batch_size=1024,
        early_stopping_patience=50,
        check_val_every_n_epoch=1,
    )

    latent = model.get_latent_representation()
    adata.obsm["X_peakVI"] = latent
    np.save(output_dir / "latent.npy", latent)
    plot_training_history(model, output_dir / "training_history.png")
    model.save(output_dir / "peakvi_model", overwrite=True)

    logger.info("Computing neighbors, Leiden clusters, and UMAP coordinates.")
    sc.pp.neighbors(adata, use_rep="X_peakVI", random_state=0)
    sc.tl.leiden(
        adata,
        resolution=1,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=0,
    )
    sc.tl.umap(adata, random_state=0, min_dist=0.1, spread=1)
    adata.obsm["X_umap_peakvi"] = adata.obsm["X_umap"].copy()
    adata.write_h5ad(output_dir / "peakvi_integrated.h5ad", compression="lzf")

    sc.pl.umap(
        adata,
        color=["leiden"],
        title="PeakVI clustering - Leiden",
        frameon=False,
        legend_loc="on data",
        show=False,
    )
    plt.savefig(output_dir / "peakvi_leiden.png", dpi=100, bbox_inches="tight")
    plt.close()

    logger.info("PeakVI integration completed successfully.")
    return model


def main(argv: list[str] | None = None) -> None:
    """Run the complete PeakVI integration workflow."""
    args = parse_args(argv)
    output_dir = prepare_output_directory(args.output_dir)
    logger = configure_logging(output_dir)

    try:
        log_software_versions(logger)

        import scanpy as sc

        logger.info("Reading post-QC ATAC H5AD: %s", args.input_h5ad)
        adata = sc.read_h5ad(args.input_h5ad)

        logger.info("Validating complete input AnnData.")
        validate_input(adata)

        # Optional: remove barcodes listed in a headerless one-column text file.
        # import pandas as pd
        # excluded_barcodes = pd.read_csv(
        #     "barcodes_to_exclude.txt", header=None, dtype=str
        # )[0]
        # adata = adata[~adata.obs_names.isin(excluded_barcodes)].copy()

        logger.info("Filtering peaks by prevalence and binarizing counts.")
        adata = prepare_anndata(adata)
        run_peakvi(adata, output_dir)
    except Exception:
        logger.exception("PeakVI workflow failed.")
        raise


if __name__ == "__main__":
    main()
