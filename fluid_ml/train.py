import torch
from torch import optim, nn
from torch.utils.data import DataLoader, random_split
from pathlib import Path
import matplotlib.pyplot as plt
from dataloader import FluidDataset
from model import FluidCNN
import yaml


class Trainer:
    def __init__(self, dataset, config=None):
        cfg = config or {}
        n_samples = len(dataset)
        batch_size = cfg.get("batch_size", 32)

        train_ratio = cfg.get("train_ratio", 0.8)
        n_train = int(train_ratio * n_samples)
        n_test = n_samples - n_train
        train_set, test_set = random_split(dataset, [n_train, n_test])

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = FluidCNN(config=cfg).to(self.device)

        lr = cfg.get("lr", 5e-5)
        weight_decay = cfg.get("weight_decay", 0)
        self.optimizer = optim.Adam(self.model.parameters(), lr=lr, weight_decay=weight_decay)

        self.scheduler = None
        if cfg.get("use_lr_scheduler", False):
            self.scheduler = optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer, mode='min', factor=0.5, patience=10
            )

        self.use_early_stopping = cfg.get("use_early_stopping", False)
        if self.use_early_stopping:
            self.early_stopping_patience = cfg.get("early_stopping_patience", 20)
            self.early_stopping_counter = 0

        criterion_class = cfg.get("criterion", nn.MSELoss)
        self.criterion = criterion_class()

        self.train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True)
        self.test_loader = DataLoader(test_set, batch_size=batch_size, shuffle=False)

        self.train_losses = []
        self.test_losses = []
        self.best_test_loss = float('inf')
        self.save_path = cfg.get("model_save_path", "model.pt")
        self.epochs = cfg.get("epochs", 5000)

        print(f"Training on {self.device}")

    def _get_test_loss(self, epoch):
        self.model.eval()
        total_test_loss = 0.0

        with torch.no_grad():
            for inputs, labels in self.test_loader:
                inputs, labels = inputs.to(self.device), labels.to(self.device)
                outputs = self.model(inputs)
                loss = self.criterion(outputs, labels)
                total_test_loss += loss.item()

        avg_test_loss = total_test_loss / len(self.test_loader)
        return avg_test_loss

    def train(self):
        status = ""

        for epoch in range(self.epochs):
            self.model.train()
            epoch_loss = 0.0

            for inputs, labels in self.train_loader:
                inputs, labels = inputs.to(self.device), labels.to(self.device)

                self.optimizer.zero_grad()

                outputs = self.model(inputs)
                loss = self.criterion(outputs, labels)

                loss.backward()
                self.optimizer.step()

                epoch_loss += loss.item()

            avg_train_loss = epoch_loss / len(self.train_loader)
            avg_test_loss = self._get_test_loss(epoch)

            if self.scheduler:
                self.scheduler.step(avg_test_loss)

            self.train_losses.append(avg_train_loss)
            self.test_losses.append(avg_test_loss)

            if avg_test_loss < self.best_test_loss:
                self.best_test_loss = avg_test_loss
                torch.save(self.model.state_dict(), self.save_path)
                status = f" (last saved model: {epoch})"
                if self.use_early_stopping:
                    self.early_stopping_counter = 0
            elif self.use_early_stopping:
                self.early_stopping_counter += 1
                if self.early_stopping_counter >= self.early_stopping_patience:
                    print(f"Early stopping at epoch {epoch}")
                    break

            if epoch % 100 == 0 or epoch == 0:
                print(f"Epoch {epoch}/{self.epochs} | "
                      f"Train: {avg_train_loss:.4e} | "
                      f"Test: {avg_test_loss:.4e}{status}")
                status = ""

    def plot_losses(self, save=False):
        plt.figure()
        plt.semilogy(self.train_losses, label="Train loss")
        plt.semilogy(self.test_losses, label="Test loss")
        plt.xlabel("Epoch")
        plt.ylabel("MSE (log scale)")
        plt.legend()
        plt.title("Training and test loss")
        if not save:
            plt.show()
        else:
            plt.savefig(Path(self.save_path.parent / "losses.png"))

    def save_stats(self, save_plot = False):
        save_path = self.save_path.parent / "stats.yaml"
        stats = {
            "best_train_loss": min(self.train_losses),
            "best_test_loss": self.best_test_loss,
            "epochs": len(self.train_losses),
        }
        with open(save_path, "w") as f:
            yaml.dump(stats, f)
        if save_plot:
            self.plot_losses(save=True)


if __name__ == "__main__":
    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"

    dataset = FluidDataset()
    dataset.create(train_files_path)

    config = {"epochs": 5000, "batch_size": 32, "lr": 5e-5}
    trainer = Trainer(dataset, config=config)
    trainer.train()
