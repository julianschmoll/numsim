from pathlib import Path
import logging

from matplotlib import pyplot as plt
import torch
import yaml
from torch.utils.data import DataLoader, random_split

from dataloader import FluidDataset
from model import FluidCNN
import constants


class Trainer:
    """Trainer class for training a FluidCNN model on a FluidDataset."""
    def __init__(self, dataset, config=None):
        """Initialize the Trainer with dataset and configuration.

        Parameters:
            dataset (FluidDataset): The dataset to train on.
            config (dict, optional): Configuration parameters for training.
        """
        cfg = config or {}
        n_train = int(
            cfg.get(
                constants.TRAIN_RATIO_KEY,
                constants.DEFAULT_TRAIN_RATIO) * len(dataset)
        )
        n_test = len(dataset) - n_train
        train_set, test_set = random_split(dataset, [n_train, n_test])

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = FluidCNN(config=cfg).to(self.device)
        self._optimizer = torch.optim.Adam(
            self.model.parameters(),
            lr=cfg.get(constants.LEARNING_RATE_KEY, constants.DEFAULT_LR),
            weight_decay=cfg.get(constants.WEIGHT_DECAY_KEY, 0)
        )
        self._scheduler = None
        if cfg.get(constants.USE_LR_SCHEDULER_KEY):
            self._scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self._optimizer, mode="min", factor=0.5, patience=10
            )

        self._use_early_stopping = cfg.get(constants.USE_EARLY_STOPPING_KEY, False)
        if self._use_early_stopping:
            self._early_stopping_patience = cfg.get(
                constants.EARLY_STOPPING_PATIENCE_KEY,
                constants.DEFAULT_EARLY_STOPPING_PATIENCE
            )
            self._early_stopping_counter = 0

        self._criterion = cfg.get(constants.CRITERION_KEY, torch.nn.MSELoss)()

        self._train_loader = DataLoader(
            train_set, batch_size=cfg.get(
                constants.BATCH_SIZE_KEY, constants.DEFAULT_BATCH_SIZE
            ), shuffle=True
        )
        self._test_loader = DataLoader(
            test_set, batch_size=cfg.get(
                constants.BATCH_SIZE_KEY, constants.DEFAULT_BATCH_SIZE
            ), shuffle=False
        )

        self._train_losses = []
        self._test_losses = []
        self._best_test_loss = float("inf")
        self._save_path = Path(
            cfg.get(constants.PATHS_KEY, {}).get(
                constants.MODEL_SAVE_PATH_KEY, "model.pt"
            )
        )
        self._epochs = cfg.get(constants.EPOCHS_KEY, constants.DEFAULT_EPOCHS)

        self._early_stopping_counter = 0

        self.log = logging.getLogger("Trainer")
        self.log.info(f"Training on {self.device}")

    def train(self):
        """Trains the model according to the specified configuration."""
        for epoch in range(self._epochs):
            last_saved_text = ""
            epoch_loss = self._train_epoch()
            avg_train_loss = epoch_loss / len(self._train_loader)
            avg_test_loss = self._get_test_loss()

            if self._scheduler:
                self._scheduler.step(avg_test_loss)

            self._train_losses.append(avg_train_loss)
            self._test_losses.append(avg_test_loss)

            if avg_test_loss < self._best_test_loss:
                self._best_test_loss = avg_test_loss
                torch.save(self.model.state_dict(), self._save_path)
                last_saved_text = f" (last saved model: {epoch})"
                self._early_stopping_counter = 0

            if self._early_stop(epoch):
                break

            if epoch % 100 == 0:
                self.log.info(
                    f"Epoch {epoch}/{self._epochs} | "
                    f"Train: {avg_train_loss:.4e} | "
                    f"Test: {avg_test_loss:.4e}{last_saved_text}"
                )

    def save_stats(self):
        """Saves model statistics to a YAML file."""
        save_path = self._save_path.parent / "stats.yaml"
        stats = {
            "best_train_loss": min(self._train_losses),
            "best_test_loss": self._best_test_loss,
            "epochs": len(self._train_losses),
        }
        with open(save_path, "w") as stats_file:
            yaml.dump(stats, stats_file)
        self._save_plot()

    def _save_plot(self):
        """Plot training and test losses over epochs and save the figure."""
        plt.figure()
        plt.semilogy(self._train_losses, label="Train loss")
        plt.semilogy(self._test_losses, label="Test loss")
        plt.xlabel("Epoch")
        plt.ylabel("MSE (log scale)")
        plt.legend()
        plt.title("Training and test loss")
        plt.savefig(Path(self._save_path.parent / "losses.png"))

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

        for inputs, labels in self._train_loader:
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            self._optimizer.zero_grad()

            outputs = self.model(inputs)
            loss = self._criterion(outputs, labels)

            loss.backward()
            self._optimizer.step()

            epoch_loss += loss.item()
        return epoch_loss

    def _get_test_loss(self):
        """Calculates the average test loss over the test dataset.

        Returns:
            float: The average test loss.
        """
        self.model.eval()
        total_test_loss = .0

        with torch.no_grad():
            for inputs, labels in self._test_loader:
                inputs = inputs.to(self.device)
                labels = labels.to(self.device)
                outputs = self.model(inputs)
                loss = self._criterion(outputs, labels)
                total_test_loss += loss.item()

        return total_test_loss / len(self._test_loader)


if __name__ == "__main__":
    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"

    dataset = FluidDataset(train_files_path)

    config = {
        constants.EPOCHS_KEY: constants.DEFAULT_EPOCHS,
        constants.BATCH_SIZE_KEY: constants.DEFAULT_BATCH_SIZE,
        constants.LEARNING_RATE_KEY: constants.DEFAULT_LR
    }
    trainer = Trainer(dataset, config=config)
    trainer.train()
