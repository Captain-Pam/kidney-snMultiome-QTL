#!/usr/bin/env python3
"""Run MultiVI integration from a prepared AnnData or MuData object.

Purpose:
    Train MultiVI on a joint RNA--ATAC object after project-specific quality
    control, cell removal, and feature preparation have been completed.

Input:
    One ``.h5ad`` AnnData object or one ``.h5mu`` MuData object. Count matrices
    must be stored in ``X``. The object must contain ``Sample`` and ``modality``
    in ``obs``. AnnData features must additionally contain ``var['modality']``
    with ``Gene_Expression`` and ``Peaks`` labels. MuData inputs must contain
    modalities named ``rna`` and ``atac``.

Output:
    A trained model, latent representation, training-history plot, Leiden UMAP,
    timestamped log, and an integrated ``mdata.h5ad`` or ``mdata.h5mu`` file.

Example:
    python 1_run_multivi_integration.py \
        --input-data /path/to/prepared_joint_data.h5ad \
        --output-dir /path/to/multivi_results
"""

from __future__ import annotations

import argparse
from datetime import datetime
import gc
import logging
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse input and output paths."""
    parser = argparse.ArgumentParser(
        description="Run MultiVI from a prepared AnnData or MuData object."
    )
    parser.add_argument(
        "--input-data",
        required=True,
        type=Path,
        help="Prepared joint RNA--ATAC data in .h5ad or .h5mu format.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for the trained model and integration results.",
    )
    return parser.parse_args()


def configure_logging(output_dir: Path) -> logging.Logger:
    """Write progress messages to the terminal and a timestamped log file."""
    log_path = output_dir / f"multivi_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def read_input(input_path: Path, logger: logging.Logger):
    """Read an AnnData or MuData object based on its file extension."""
    suffix = input_path.suffix.lower()
    if suffix == ".h5ad":
        import scanpy as sc

        logger.info("Reading AnnData: %s", input_path)
        return sc.read_h5ad(input_path), "anndata"

    if suffix == ".h5mu":
        import muon as mu

        logger.info("Reading MuData: %s", input_path)
        mdata = mu.read_h5mu(input_path)
        mdata.update()
        mu.pp.intersect_obs(mdata)
        return mdata, "mudata"

    raise ValueError("Input must use the .h5ad or .h5mu extension.")


def validate_input(data, data_format: str) -> None:
    """Check metadata and modality labels required by the final workflow."""
    missing_obs = [key for key in ("Sample", "modality") if key not in data.obs]
    if missing_obs:
        raise ValueError(
            "Input obs is missing required columns: " + ", ".join(missing_obs)
        )

    if data_format == "anndata":
        if "modality" not in data.var:
            raise ValueError("AnnData var must contain the modality column.")
        feature_modalities = set(data.var["modality"].astype(str))
        required = {"Gene_Expression", "Peaks"}
        if not required.issubset(feature_modalities):
            raise ValueError(
                "AnnData var['modality'] must contain Gene_Expression and Peaks."
            )
    else:
        missing_modalities = [key for key in ("rna", "atac") if key not in data.mod]
        if missing_modalities:
            raise ValueError(
                "MuData is missing required modalities: "
                + ", ".join(missing_modalities)
            )


def clear_unused_layers(data, data_format: str, logger: logging.Logger) -> None:
    """Remove unused layers while retaining count matrices stored in X."""
    logger.info("Removing unused layers; count matrices in X are retained.")
    if data_format == "anndata":
        data.layers.clear()
    else:
        for modality in data.mod.values():
            modality.layers.clear()
    gc.collect()


def setup_model(data, data_format: str, logger: logging.Logger):
    """Register the selected container and initialize MultiVI."""
    import scvi

    if data_format == "mudata":
        logger.info("Registering MuData for MultiVI.")
        scvi.model.MULTIVI.setup_mudata(
            data,
            modalities={"rna_layer": "rna", "atac_layer": "atac"},
            batch_key="modality",
            categorical_covariate_keys=["Sample"],
        )
        return scvi.model.MULTIVI(
            data,
            n_genes=data.mod["rna"].n_vars,
            n_regions=data.mod["atac"].n_vars,
        )

    logger.info("Registering AnnData for MultiVI.")
    scvi.model.MULTIVI.setup_anndata(
        data,
        batch_key="modality",
        categorical_covariate_keys=["Sample"],
    )
    return scvi.model.MULTIVI(
        data,
        n_genes=int((data.var["modality"] == "Gene_Expression").sum()),
        n_regions=int((data.var["modality"] == "Peaks").sum()),
    )


def plot_training_history(model, output_path: Path) -> None:
    """Plot training and validation reconstruction loss."""
    import matplotlib.pyplot as plt

    train_key = "reconstruction_loss_train"
    validation_key = "reconstruction_loss_validation"
    train_loss = model.history[train_key][train_key]
    validation_loss = model.history[validation_key][validation_key]

    figure, axis = plt.subplots(figsize=(8, 6))
    axis.plot(train_loss, label="train")
    axis.plot(validation_loss, label="validation")
    axis.axhline(validation_loss.min(), color="black")
    axis.set_xlabel("Epoch")
    axis.set_ylabel("Reconstruction loss")
    axis.legend()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def save_embedding_and_clusters(data, data_format: str, model, output_dir: Path) -> None:
    """Save the latent space and calculate neighbors, UMAP, and Leiden clusters."""
    import matplotlib.pyplot as plt
    import numpy as np
    import scanpy as sc

    latent = model.get_latent_representation()
    data.obsm["multivi"] = latent
    np.save(output_dir / "latent.npy", latent)

    sc.pp.neighbors(data, use_rep="multivi")
    sc.tl.umap(data, min_dist=0.1)
    sc.tl.leiden(data)

    if data_format == "mudata":
        data.write_h5mu(output_dir / "mdata.h5mu")
    else:
        data.write_h5ad(output_dir / "mdata.h5ad", compression="lzf")

    sc.pl.umap(
        data,
        color=["leiden"],
        title="MultiVI clustering - Leiden",
        frameon=False,
        legend_loc="on data",
        show=False,
    )
    plt.savefig(output_dir / "multivi_leiden.png", dpi=100, bbox_inches="tight")
    plt.close()


def main() -> None:
    """Run the final MultiVI integration workflow."""
    args = parse_args()

    import scvi
    import torch

    args.output_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(args.output_dir)
    logger.info("Script started")
    logger.info("Output directory: %s", args.output_dir)

    scvi.settings.seed = 0
    torch.set_float32_matmul_precision("high")

    data, data_format = read_input(args.input_data, logger)
    validate_input(data, data_format)
    clear_unused_layers(data, data_format, logger)
    model = setup_model(data, data_format, logger)

    logger.info("Starting training")
    model.train(
        max_epochs=60,
        batch_size=4096,
        check_val_every_n_epoch=1,
        early_stopping=True,
    )
    logger.info("Training completed")

    plot_training_history(model, args.output_dir / "train.jpg")
    model.save(args.output_dir / "trained_model", overwrite=True)
    save_embedding_and_clusters(data, data_format, model, args.output_dir)
    logger.info("All done")


if __name__ == "__main__":
    main()
