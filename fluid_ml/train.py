"""Trainer module for FluidCNN model training."""
from pathlib import Path
import logging

import torch
import yaml
from torch.utils.data import DataLoader, Subset, random_split

from model import FluidCNN
import evaluate
from constants import *  # noqa: F403, WPS347


def _print_epoch(epoch: int) -> bool:
    """Checks if the current epoch should be printed.

    Args:
        epoch: The current epoch number.

    Returns:
        If the epoch should be printed.
    """
    return epoch % 100 == 0 or epoch == 1


def _get_subsets(cfg, dataset):
    """Splits the dataset into training, validation and test subsets.

    Args:
        cfg: Configuration dict.
        dataset: The dataset to split.

    Returns:
        The test, train and validation subsets.
    """
    n_total = len(dataset)
    n_train = int(cfg.get(TRAIN_RATIO, DEFAULT_TRAIN_RATIO) * n_total)
    n_test = int((n_total - n_train) / 2)
    n_val = n_total - n_train - n_test

    split_fn = random_split if cfg.get(RANDOM_SPLIT) else _sorted_split
    return split_fn(dataset, [n_train, n_test, n_val])


def _sorted_split(dataset, lengths):
    """Splits the dataset into subsets of given lengths without shuffling.

    This is implemented here because we hope to archive better extrapolation
    performance, when the model trains on the earlier part of the dataset and
    is tested/validated on the later part.

    Args:
        dataset: The dataset to split.
        lengths: A list of lengths for each subset.

    Returns:
        A list of Subset objects corresponding to the splits.
    """
    if sum(lengths) != len(dataset):
        raise ValueError(
            "Sum of input lengths does not equal the length of the input dataset!"
        )

    sets = []
    start_idx = 0

    for length in lengths:
        end_idx = start_idx + length
        indices = range(start_idx, end_idx)
        sets.append(Subset(dataset, indices))
        start_idx = end_idx

    return sets


class Trainer:  # noqa: WPS230, pylint: disable=too-many-instance-attributes
    """Trainer class for training a FluidCNN model on a FluidDataset."""
    def __init__(self, dataset, config=None):
        """Initialize the Trainer with dataset and configuration.

        Parameters:
            dataset (FluidDataset): The dataset to train on.
            config (dict, optional): Configuration parameters for training.
        """
        cfg = config or {}
        train_set, test_set, val_set = _get_subsets(cfg, dataset)
        self.dataset = dataset

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = FluidCNN(config=cfg).to(self.device)
        self._optimizer = torch.optim.Adam(
            self.model.parameters(),
            lr=cfg.get(LEARNING_RATE, DEFAULT_LR),
            weight_decay=cfg.get(WEIGHT_DECAY, 0)
        )
        self._scheduler = None
        if cfg.get(USE_LR_SCHEDULER):
            self._scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self._optimizer, mode="min", factor=0.5, patience=10
            )

        self._use_early_stopping = cfg.get(USE_EARLY_STOPPING, False)
        if self._use_early_stopping:
            self._early_stopping_patience = cfg.get(
                EARLY_STOPPING_PATIENCE, DEFAULT_EARLY_STOPPING_PATIENCE
            )
            self._early_stopping_counter = 0

        self._criterion = cfg.get(CRITERION, torch.nn.MSELoss)()

        batch_size = cfg.get(BATCH_SIZE, DEFAULT_BATCH_SIZE)
        self.train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True)
        self.val_loader = DataLoader(val_set, batch_size=batch_size, shuffle=False)
        self.test_loader = DataLoader(test_set, batch_size=batch_size, shuffle=False)

        self.train_losses = []
        self.val_losses = []
        self.test_losses = []

        self._best_val_loss = float("inf")
        self._save_path = Path(cfg.get(PATHS, {}).get(MODEL, "model.pt"))
        self._epochs = cfg.get(EPOCHS, DEFAULT_EPOCHS)

        self._early_stopping_counter = 0

        self.log = logging.getLogger("Trainer")
        self.log.info(f"Training on {self.device}")

    def train(self):
        """Trains the model according to the specified configuration."""
        last_saved_text = ""
        for epoch in range(1, self._epochs + 1):
            epoch_loss = self._train_epoch()

            self.train_losses.append(epoch_loss / len(self.train_loader))
            self.test_losses.append(self._get_loss(self.test_loader))
            self.val_losses.append(self._get_loss(self.val_loader))

            if self._scheduler:
                self._scheduler.step(self.val_losses[-1])

            if self.val_losses[-1] < self._best_val_loss:
                self._best_val_loss = self.val_losses[-1]
                torch.save(self.model.state_dict(), self._save_path)
                last_saved_text = f" (last saved model: {epoch})"
                self._early_stopping_counter = 0

            if self._early_stop(epoch):
                break

            if _print_epoch(epoch):
                self.log.info(
                    f"Epoch {epoch}/{self._epochs} | "  # noqa: WPS237
                    f"Train: {self.train_losses[-1]:.4e} | "
                    f"Val: {self.val_losses[-1]:.4e} | "
                    f"Test: {self.test_losses[-1]:.4e}{last_saved_text}"
                )
                last_saved_text = ""

    def save_stats(self):
        """Saves model statistics and evaluation metrics to a YAML file.

        Returns:
            Path: Path to the saved YAML file.
        """
        save_path = self._save_path.parent / "stats.yaml"
        stats = {
            "best_train_loss": min(self.train_losses),
            "best_val_loss": self._best_val_loss,
            "best_test_loss": min(self.test_losses),
            "epochs": len(self.train_losses),
            "metrics": {
                "train": evaluate.evaluation_metrics(
                    self.train_loader, self.model, self.dataset, self.device
                ),
                "validation": evaluate.evaluation_metrics(
                    self.val_loader, self.model, self.dataset, self.device
                ),
                "test": evaluate.evaluation_metrics(
                    self.test_loader, self.model, self.dataset, self.device
                )
            },
        }
        with open(save_path, "w", encoding="utf-8") as stats_file:
            yaml.dump(stats, stats_file)
        return save_path

    def _early_stop(self, epoch):
        """Checks if early stopping criteria are met.

        Parameters:
            epoch (int): The current epoch number.

        Returns:
            bool: True if training should stop early, False otherwise.
        """
        if self._use_early_stopping:
            self._early_stopping_counter += 1
            if self._early_stopping_counter >= self._early_stopping_patience:
                self.log.info(f"Early stopping at epoch {epoch}")
                return True
        return False

    def _train_epoch(self):
        """Trains the model for one epoch.

        Returns:
            float: The total loss for the epoch.
        """
        self.model.train()
        epoch_loss = .0

        for inputs, labels in self.train_loader:
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            self._optimizer.zero_grad()

            outputs = self.model(inputs)
            loss = self._criterion(outputs, labels)

            loss.backward()
            self._optimizer.step()

            epoch_loss += loss.item()
        return epoch_loss

    def _get_loss(self, loader):
        """Calculates average test loss over dataset specified by loader.

        Args:
            loader (DataLoader): DataLoader for the dataset to evaluate.

        Returns:
            float: The average test loss.
        """
        self.model.eval()
        total_loss = .0

        with torch.no_grad():
            for inputs, labels in loader:
                inputs = inputs.to(self.device)
                labels = labels.to(self.device)
                outputs = self.model(inputs)
                loss = self._criterion(outputs, labels)
                total_loss += loss.item()

        return total_loss / len(loader)
