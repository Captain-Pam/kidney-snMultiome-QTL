#!/usr/bin/env python3
"""Run scVI integration for a post-QC single-cell RNA AnnData object.

Purpose:
    Validate and preprocess post-QC RNA counts, train scVI, and construct a
    reproducible latent-space neighborhood graph, Leiden clustering, and UMAP.

Input:
    An H5AD file with raw non-negative integer RNA counts in ``.X`` and the
    ``orig_ident``, ``percent_mito``, and ``nCount_RNA`` columns in ``.obs``.

Output:
    A trained model directory, integrated H5AD file, training-history plot,
    and workflow log in the requested output directory.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


LOGGER_NAME = "scvi_workflow"
REQUIRED_OBS_COLUMNS = ("orig_ident", "percent_mito", "nCount_RNA")
CONTINUOUS_COVARIATES = ("percent_mito", "nCount_RNA")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the two required workflow paths."""
    parser = argparse.ArgumentParser(
        description="Integrate post-QC single-cell RNA data with scVI."
    )
    parser.add_argument(
        "--input-h5ad",
        required=True,
        type=Path,
        help="Post-QC RNA H5AD with raw counts in .X.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for the model, integrated H5AD, plot, and log.",
    )
    return parser.parse_args(argv)


def configure_logging(output_dir: Path) -> logging.Logger:
    """Log each workflow message once to standard output and a file."""
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
    file_handler = logging.FileHandler(output_dir / "scvi_workflow.log")
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
    """Validate required metadata and raw count values."""
    import numpy as np

    missing_columns = [
        column for column in REQUIRED_OBS_COLUMNS if column not in adata.obs
    ]
    if missing_columns:
        raise ValueError(
            "Missing required adata.obs columns: " + ", ".join(missing_columns)
        )

    # Continuous covariates must be usable as finite numeric model inputs.
    for column in CONTINUOUS_COVARIATES:
        try:
            values = np.asarray(adata.obs[column], dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"adata.obs['{column}'] must contain numeric finite values."
            ) from error
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"adata.obs['{column}'] must contain numeric finite values."
            )

    # Sparse matrices expose stored values through .data; dense matrices do not.
    count_values = adata.X.data if hasattr(adata.X, "format") else np.asarray(adata.X)
    try:
        valid_counts = (
            np.all(np.isfinite(count_values))
            and np.all(count_values >= 0)
            and np.all(np.equal(count_values, np.floor(count_values)))
        )
    except (TypeError, ValueError):
        valid_counts = False
    if not valid_counts:
        raise ValueError("adata.X must contain raw non-negative integer counts.")


def prepare_anndata(adata, n_top_genes: int = 5000):
    """Filter, normalize, log-transform, and select batch-aware HVGs."""
    import scanpy as sc

    validate_input(adata)
    adata.layers["counts"] = adata.X.copy()

    # Apply the same ordered preprocessing steps used by the public workflow.
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=n_top_genes,
        layer="counts",
        batch_key="orig_ident",
        subset=True,
    )
    return adata


def plot_training_history(model, output_path: Path) -> None:
    """Save scVI training and validation reconstruction losses."""
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
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_scvi(adata, output_dir: Path):
    """Train scVI and write the model, embedding, clustering, and H5AD output."""
    import scanpy as sc
    import scvi
    import torch

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(LOGGER_NAME)

    # Fix stochastic settings before constructing the model.
    scvi.settings.seed = 0
    torch.set_float32_matmul_precision("high")

    logger.info("Registering AnnData fields for scVI.")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="orig_ident",
        continuous_covariate_keys=["percent_mito", "nCount_RNA"],
    )
    model = scvi.model.SCVI(
        adata,
        n_layers=3,
        n_latent=30,
        gene_likelihood="zinb",
        dropout_rate=0.1,
    )

    logger.info("Training scVI for up to 150 epochs.")
    model.train(
        max_epochs=150,
        plan_kwargs={"lr": 1e-4},
        early_stopping=True,
        batch_size=2048,
        early_stopping_patience=25,
        check_val_every_n_epoch=1,
    )

    logger.info("Saving the latent representation, history, and trained model.")
    adata.obsm["X_scVI"] = model.get_latent_representation()
    plot_training_history(model, output_dir / "training_history.png")
    model.save(output_dir / "scvi_model", overwrite=True)

    # Build deterministic downstream embeddings from the scVI latent space.
    logger.info("Computing neighbors, Leiden clusters, and UMAP.")
    sc.pp.neighbors(adata, use_rep="X_scVI", random_state=0)
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
    adata.obsm["X_umap_scvi"] = adata.obsm["X_umap"].copy()

    output_h5ad = output_dir / "scvi_integrated.h5ad"
    adata.write_h5ad(output_h5ad, compression="lzf")
    logger.info("Workflow complete: %s", output_h5ad)
    return model


def main(argv: list[str] | None = None) -> None:
    """Load the input, run preprocessing and integration, and write outputs."""
    args = parse_args(argv)
    logger = configure_logging(args.output_dir)

    try:
        log_software_versions(logger)
        logger.info("Reading post-QC RNA data: %s", args.input_h5ad)
        import scanpy as sc

        adata = sc.read_h5ad(args.input_h5ad)

        # Optional: remove barcodes listed in a headerless one-column text file.
        # import pandas as pd
        # excluded_barcodes = pd.read_csv(
        #     "barcodes_to_exclude.txt", header=None, dtype=str
        # )[0]
        # adata = adata[~adata.obs_names.isin(excluded_barcodes)].copy()

        logger.info("Validating and preprocessing the AnnData object.")
        adata = prepare_anndata(adata)
        run_scvi(adata, args.output_dir)
    except Exception:
        logger.exception("scVI workflow failed.")
        raise


if __name__ == "__main__":
    main()
