import argparse
import datetime
import logging
import json
import shutil
from pathlib import Path

from dataloader import FluidDataset
from train import Trainer
from submit import generate_submission
import constants


def main(config_path: str | Path | None):
    """Entrypoint to train the model.

    Args:
        config_path: Path to configuration file.
    """
    config = get_config(config_path)
    save_config(config)
    save_path = Path(config[constants.PATHS_KEY][constants.BASE_SAVE_PATH_KEY])

    dataset = FluidDataset(
        Path(config[constants.PATHS_KEY][constants.TRAIN_FILES_PATH_KEY])
    )
    dataset.normalize()
    dataset.save(save_path)

    save_model_init(config)

    trainer = Trainer(dataset, config=config)
    trainer.train()
    trainer.save_stats()

    inputs_path = (Path(__file__).resolve().parent.parent /
                   constants.RESOURCE_FOLDER_NAME /
                   constants.INPUTS_FILE_NAME)
    generate_submission(save_path, inputs_path)


def get_config(config_path: str | Path | None) -> dict:
    """Reads configuration for the training run from file.

    Args:
        config_path: Path to configuration file.

    Returns:
        A dictionary containing configuration parameters.
    """
    config = {}
    if config_path and Path(config_path).exists():
        with open(config_path, 'r') as config_file:
            config = json.load(config_file)

    current_file_path = Path(__file__).resolve()
    train_files_path = (current_file_path.parent.parent /
                        constants.BUILD_PATH /
                        constants.TRAIN_FOLDER_NAME)
    model_path = (current_file_path.parent.parent /
                  constants.MODELS_SAVE_PATH /
                  get_unique_folder_name())
    model_path.mkdir(parents=True, exist_ok=True)

    config[constants.PATHS_KEY] = {
        constants.TRAIN_FILES_PATH_KEY: str(train_files_path),
        constants.BASE_SAVE_PATH_KEY: str(model_path),
        constants.MODEL_SAVE_PATH_KEY: str(model_path / constants.DEFAULT_MODEL_NAME),
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
    save_path = Path(
        config[constants.PATHS_KEY][constants.BASE_SAVE_PATH_KEY]
    ) / constants.CONFIG_FILE_NAME
    save_data = config.copy()
    if constants.PATHS_KEY in save_data:
        save_data.pop(constants.PATHS_KEY)
    with open(save_path, "w") as config_file:
        json.dump(save_data, config_file, indent=4, default=str)


def save_model_init(config: dict) -> None:
    """
    Copy the model definition file to the submission directory.

    Args:
        config: The configuration dictionary containing the save path.
    """
    save_path = Path(
        config[constants.PATHS_KEY][constants.BASE_SAVE_PATH_KEY]
    ) / constants.INIT_MODEL_PY
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

    config = args.config
    if not config:
        config = (Path(__file__).resolve().parent.parent
                  / "cfg" / "ml" / "config.json"
                  )

    logging.basicConfig(level=logging.INFO)
    main(config_path=config)
