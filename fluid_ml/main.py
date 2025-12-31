import datetime
import logging
import json
import shutil
from pathlib import Path

from dataloader import FluidDataset
from train import Trainer
from submit import generate_submission
import constants


def main():
    """Set up configuration, data, and training for the model."""
    config = get_config()
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
    trainer.save_stats(save_plot=True)

    inputs_path = (Path(__file__).resolve().parent.parent /
                   constants.RESOURCE_FOLDER_NAME /
                   constants.INPUTS_FILE_NAME)
    generate_submission(save_path, inputs_path)


# this could either read config from file, get from CLI, ...
def get_config() -> dict:
    """
    Get the base configuration for the training run.

    Returns:
        A dictionary containing configuration parameters.
    """
    config = {
        constants.EPOCHS_KEY: constants.DEFAULT_EPOCHS,
        constants.BATCH_SIZE_KEY: constants.DEFAULT_BATCH_SIZE,
        constants.LEARNING_RATE_KEY: constants.DEFAULT_LR,
        constants.NUM_HIDDEN_LAYERS_KEY: constants.DEFAULT_NUM_HIDDEN_LAYERS,
        constants.KERNEL_SIZE_KEY: constants.DEFAULT_KERNEL_SIZE,
        constants.IN_CHANNELS_KEY: constants.DEFAULT_IN_CHANNELS,
        constants.OUT_CHANNELS_KEY: constants.DEFAULT_OUT_CHANNELS,
        constants.HIDDEN_CHANNELS_KEY: constants.DEFAULT_HIDDEN_CHANNELS,
        constants.OUTPUT_ACTIVATION_KEY: constants.DEFAULT_OUTPUT_ACTIVATION,
        constants.USE_BIAS_KEY: constants.DEFAULT_USE_BIAS,
        constants.PADDING_MODE_KEY: constants.DEFAULT_PADDING_MODE,
    }

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


# this is nasty, would still like to have something like the current date/time
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
    logging.basicConfig(level=logging.INFO)
    main()
