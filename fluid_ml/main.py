"""Entrypoint to train the fluid model and generate submission."""
import argparse
import datetime
import logging
import json
import shutil
from pathlib import Path
import torch

import visualize
from constants import *  # noqa: F403, WPS347
from dataloader import FluidDataset
import submit
from train import Trainer


def main(config_path: str | Path | None):
    """Entrypoint to train the model.

    Args:
        config_path: Path to configuration file.
    """
    config = get_config(config_path)
    save_config(config)
    save_path = Path(config[PATHS][SAVE_PATH])

    torch.manual_seed(config.get(RANDOM_SEED, DEFAULT_RANDOM_SEED))

    trainer = _train_model(config, save_path)
    _write_submission_data(save_path, trainer)


def _write_submission_data(save_path, trainer):
    inputs_path = (Path(__file__).resolve().parent.parent /
                   RESOURCES / INPUTS_FILE_NAME)

    visualize.loss_plot(save_path / "losses.png", [
        (trainer.train_losses, "Train Loss"),
        (trainer.val_losses, "Validation Loss"),
        (trainer.test_losses, "Test Loss"),
    ], title="Training, Test and Validation Loss Over Epochs", )

    model = submit.init_model(save_path)
    submit.generate_kaggle_submission(
        model, save_path, inputs_path
    )
    submit.generate_interpolation_plots(model, trainer, save_path)
    submit.generate_extrapolation_plots(
        model, save_path, trainer.dataset.stats
    )

    submit.run_notebook(
        submit.copy_notebook(save_path)
    )


def _train_model(config, save_path):
    dataset = FluidDataset(
        Path(config[PATHS][TRAIN_FILES_PATH])
    )
    dataset.normalize()
    dataset.save(save_path)
    save_model_init(config)

    trainer = Trainer(dataset, config=config)
    try:
        trainer.train()
    except KeyboardInterrupt:
        trainer.log.info("Training interrupted by user. Saving...")
    trainer.save_stats()
    return trainer


def get_config(config_path: str | Path | None) -> dict:
    """Reads configuration for the training run from file.

    Args:
        config_path: Path to configuration file.

    Returns:
        A dictionary containing configuration parameters.
    """
    config = {}
    if config_path and Path(config_path).exists():
        with open(config_path, "r", encoding="utf-8") as config_file:
            config = json.load(config_file)

    current_file_path = Path(__file__).resolve()
    train_files_path = (current_file_path.parent.parent /
                        BUILD / TRAIN)
    model_path = (current_file_path.parent.parent /
                  MODELS / get_unique_folder_name())
    model_path.mkdir(parents=True, exist_ok=True)

    config[PATHS] = {
        TRAIN_FILES_PATH: str(train_files_path),
        SAVE_PATH: str(model_path),
        MODEL: str(model_path / DEFAULT_MODEL_NAME),
    }

    return config


def get_unique_folder_name() -> str:
    """
    Generate a unique folder name based on the current timestamp.

    Returns:
        A string representing the unique folder name.
    """
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def save_config(config: dict) -> None:
    """
    Save the configuration dictionary to a JSON file.

    Args:
        config: The configuration dictionary to save.
    """
    save_path = Path(config[PATHS][SAVE_PATH]) / CONFIG_FILE_NAME
    save_data = config.copy()
    if PATHS in save_data:
        save_data.pop(PATHS)
    with open(save_path, "w", encoding="utf-8") as config_file:
        json.dump(save_data, config_file, indent=4, default=str)


def save_model_init(config: dict) -> None:
    """
    Copy the model definition file to the submission directory.

    Args:
        config: The configuration dictionary containing the save path.
    """
    save_path = Path(config[PATHS][SAVE_PATH]) / INIT_MODEL_PY
    model_path = Path(__file__).resolve().parent / "model.py"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(model_path, save_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a fluid model.")
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to the configuration JSON file.",
    )
    args = parser.parse_args()

    cfg_path = args.config
    if not cfg_path:
        cfg_path = (Path(__file__).resolve().parent.parent
                    / "cfg" / "ml" / "config.json")

    logging.basicConfig(level=logging.INFO)
    main(config_path=cfg_path)
