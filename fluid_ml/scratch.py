from torch import optim, nn
from torch.utils.data import DataLoader

from pathlib import Path

from fluid_ml.dataloader import FluidDataset
from fluid_ml.model import FluidCNN
from fluid_ml import evaluate

import torch

def train(dataset, epochs=5000, batch_size=32):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Training on {device}")

    model = FluidCNN(config=None).to(device)
    optimizer = optim.Adam(model.parameters(), lr=5e-5)
    criterion = nn.MSELoss()

    train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model.train()
    for epoch in range(epochs):
        epoch_loss = 0.0
        for batch_idx, (inputs, labels) in enumerate(train_loader):
            inputs, labels = inputs.to(device), labels.to(device)

            outputs = model(inputs)
            loss = criterion(outputs, labels)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()

        if (epoch + 1) % 100 == 0:
            avg_loss = epoch_loss / len(train_loader)
            print(f"Epoch [{epoch + 1}/{epochs}], Loss: {avg_loss:.8f}")

    return model


if __name__ == "__main__":
    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"

    dataset = FluidDataset()
    dataset.create(train_files_path)

    trained_model = train(dataset, epochs=500)
    trained_model.eval()

    sample_input, sample_label = dataset[0]
    sample_input = sample_input.unsqueeze(0).to(next(trained_model.parameters()).device)

    with torch.no_grad():
        prediction = trained_model(sample_input).cpu()

    dataset.denormalize()

    evaluate.visualize(sample_input.cpu(), sample_label.unsqueeze(0), title="Ground Truth", quiver=True)
    evaluate.visualize(sample_input.cpu(), prediction, title="Prediction", quiver=True)
