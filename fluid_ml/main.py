import torch
from pathlib import Path
import datetime

from dataloader import FluidDataset
from train import Trainer
import shutil

import json


def main():
    config = get_config()
    save_config(config)

    dataset = FluidDataset()
    dataset.create(config["train_files_path"])
    dataset.normalize()
    dataset.save(config["save_path"])

    save_model_init(config)

    trainer = Trainer(dataset, config=config)
    trainer.train()
    trainer.save_stats(save_plot=True)

# this could either read config from file, get from CLI, ...
def get_config() -> dict[str, int | float | Path]:
    config = {"epochs": 5000, "batch_size": 32, "lr": 5e-5}

    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"
    model_path = current_file_path.parent.parent / "models" / get_unique_folder_name()
    model_path.mkdir(parents=True, exist_ok=True)

    config["train_files_path"] = train_files_path
    config["save_path"] = model_path
    config["model_save_path"] = model_path / "model.pt"

    return config


# this is nasty, would still like to have something like the current date/time so we can remember
def get_unique_folder_name():
    return str(datetime.datetime.now()).replace("-", "").replace(" ", "_").replace(":", "").split(".")[0]


def save_config(config: dict) -> None:
    save_path = Path(config["save_path"]) / "config.json"
    with open(save_path, "w") as f:
        json.dump(config, f, indent=4, default=str)


def save_model_init(config: dict[str, int | float | Path]) -> None:
    save_path = config["save_path"] / "mymodel.py"
    model_path = Path(__file__).resolve().parent / "model.py"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(model_path, save_path)

if __name__ == "__main__":
    main()
