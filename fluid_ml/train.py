from torch import optim, nn
from torch.utils.data import DataLoader, Subset

from pathlib import Path
from mathplotlib import plt
from fluid_ml.dataloader import FluidDataset
from fluid_ml.model import FluidCNN
from fluid_ml import evaluate

import torch

class Trainer:
    def __init__(self, dataset, config=None):
        config = config or {}
        n_samples = len(dataset)
        batch_size = config.get("batch_size", 32)

        train_ratio = config.get("train_ratio", 0.8)
        n_train = int(train_ratio * n_samples)

        indices = torch.randperm(n_samples)
        train_indices = indices[:n_train]
        test_indices = indices[n_train:]

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = FluidCNN(config=config).to(self.device)

        lr = config.get("lr", 5e-5)
        self.optimizer = optim.Adam(self.model.parameters(), lr=lr)

        criterion_class = config.get("criterion", nn.MSELoss)
        self.criterion = criterion_class()

        train_subset = Subset(dataset, train_indices)
        test_subset = Subset(dataset, test_indices)

        self.train_loader = DataLoader(train_subset, batch_size=batch_size, shuffle=True)
        self.test_loader = DataLoader(test_subset, batch_size=batch_size, shuffle=False)

        self.train_losses = []
        self.test_losses = []
        self.improved = False
        self.save_path = config.get("save_path", "model.pt")

        print(f"Training on {self.device}")

    def _evaluate(self):
        return
        self.model.eval()

        with torch.no_grad():
            prediction = self.model(self.test_loader)
            test_loss = self.criterion(prediction, y_test_d)

        self.train_losses.append(loss.item())
        self.test_losses.append(test_loss.item())

        if (epoch + 1) % 200 == 0:
            print(f"Epoch {epoch+1}/{n_epochs:4d}")
            print(f"    Train Loss: {loss.item():.4e})")
            print(f"    Test Loss: {test_loss.item():.4e}\n")

    def plot_losses(self):
        plt.figure()
        plt.semilogy(self.train_losses, label="Train loss")
        plt.semilogy(self.test_losses, label="Test loss")
        plt.xlabel("Epoch")
        plt.ylabel("MSE (log scale)")
        plt.legend()
        plt.title("Training and test loss")
        plt.show()

    def train(self, epochs=10):
        self.model.train()

        for epoch in range(epochs):
            epoch_loss = 0.0
            for batch_idx, (inputs, labels) in enumerate(self.train_loader):
                inputs, labels = inputs.to(self.device), labels.to(self.device)

                outputs = self.model(inputs)
                loss = self.criterion(outputs, labels)

                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

                epoch_loss += loss.item()

            avg_loss = epoch_loss / len(self.train_loader)
            self._evaluate()

            if self.improved:
                torch.save(self.model.state_dict(), "model.pt")


if __name__ == "__main__":
    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"

    dataset = FluidDataset()
    dataset.create(train_files_path)
    trainer = Trainer(dataset)
    trainer.train()
    trainer.plot_losses()