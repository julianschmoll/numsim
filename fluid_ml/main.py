import datetime
import logging
import json
import shutil
from pathlib import Path

from dataloader import FluidDataset
from train import Trainer
from submit import generate_submission

SAVE_PATH_KEY = "save_path"
MODEL_SAVE_PATH_KEY = "model_save_path"


def main():
    """Set up configuration, data, and training for the model."""
    config = get_config()
    save_config(config)

    dataset = FluidDataset()
    dataset.create(config["train_files_path"])
    dataset.normalize()
    dataset.save(config[SAVE_PATH_KEY])

    save_model_init(config)

    trainer = Trainer(dataset, config=config)
    trainer.train()
    trainer.save_stats(save_plot=True)

    inputs_path = Path(__file__).resolve().parent.parent / "resources" / "inputs.pt"
    generate_submission(config[SAVE_PATH_KEY], inputs_path)


# this could either read config from file, get from CLI, ...
def get_config() -> dict:
    """
    Get the base configuration for the training run.

    Returns:
        A dictionary containing configuration parameters.
    """
    config = {
        "epochs": 5000,
        "batch_size": 32,
        "lr": 5e-5,
        "num_hidden_layers": 5,
        "kernel_size": 7,
        "in_channels": 1,
        "out_channels": 2,
        "hidden_channels": 16,
        "output_activation": None,
        "use_bias": True,
        "padding_mode": "zeros",
    }

    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"
    model_path = current_file_path.parent.parent / "models" / get_unique_folder_name()
    model_path.mkdir(parents=True, exist_ok=True)

    config["train_files_path"] = train_files_path
    config[SAVE_PATH_KEY] = model_path
    config[MODEL_SAVE_PATH_KEY] = model_path / "model.pt"

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
    save_path = Path(config[SAVE_PATH_KEY]) / "config.json"
    with open(save_path, "w") as config_file:
        json.dump(config, config_file, indent=4, default=str)


def save_model_init(config: dict) -> None:
    """
    Copy the model definition file to the submission directory.

    Args:
        config: The configuration dictionary containing the save path.
    """
    save_path = config[SAVE_PATH_KEY] / "mymodel.py"
    model_path = Path(__file__).resolve().parent / "model.py"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(model_path, save_path)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
