import torch
import numpy as np
import pyvista as pv
import os
import yaml
from pathlib import Path
from torch.utils.data import Dataset, DataLoader

from fluid_ml import evaluate


class FluidDataset(Dataset):
    def __init__(self):
        self.inputs = None
        self.labels = None

    def __len__(self):
        return len(self.inputs) if self.inputs is not None else 0

    def __getitem__(self, idx):
        x = torch.from_numpy(self.inputs[idx]).float()
        y = torch.from_numpy(self.labels[idx]).float()
        return x, y

    def create(self, base_dir, n_samples):
        all_inputs = []
        all_labels = []
        base_path = Path(base_dir)

        for i in range(n_samples):
            folder = base_path / f"out_{i:04d}"
            if not folder.exists():
                continue

            vti_files = sorted([f for f in os.listdir(folder) if f.startswith("output_") and f.endswith(".vti")])

            for vti_file in vti_files:
                mesh = pv.read(folder / vti_file)
                nx, ny, nz = mesh.dimensions

                vel_data = mesh.point_data["velocity"].reshape((ny, nx, 3), order="F")
                vel_data = np.transpose(vel_data, (1, 0, 2))

                u_field = vel_data[:, :, 0]  # Horizontal
                v_field = vel_data[:, :, 1]  # Vertical

                # Handcrafted input: A zero-filled grid with the top lid velocity
                # Note: Index -1 is the "Top" in (ny, nx) layout
                ux_val = u_field[-1, 1:-1].mean()
                input_channel = np.zeros((1, ny, nx))
                input_channel[0, -1, 1:-1] = ux_val

                # Label: Stack u and v into (2, ny, nx)
                label_channels = np.stack([u_field, v_field], axis=0)

                all_inputs.append(input_channel)
                all_labels.append(label_channels)

        self.inputs = np.array(all_inputs)
        self.labels = np.array(all_labels)

    def normalize(self):
        """Normalizes each channel and returns min/max stats."""

        def get_stats_and_norm(data):
            stats = []
            for c in range(data.shape[1]):
                c_min, c_max = data[:, c].min(), data[:, c].max()
                if c_max - c_min > 1e-9:
                    data[:, c] = (data[:, c] - c_min) / (c_max - c_min)
                stats.append({"min": float(c_min), "max": float(c_max)})
            return stats

        in_stats = get_stats_and_norm(self.inputs)
        lb_stats = get_stats_and_norm(self.labels)
        return in_stats, lb_stats

    def save(self, folder_path):
        os.makedirs(folder_path, exist_ok=True)
        np.save(f"{folder_path}/inputs.npy", self.inputs)
        np.save(f"{folder_path}/labels.npy", self.labels)

        # I think the normalization here is wrong, we should already
        # normalize the data in set probably
        in_stats, lb_stats = self.normalize()
        config = {
            "inputs": {"u": in_stats[0]},
            "labels": {"u": lb_stats[0], "v": lb_stats[1]}
        }
        with open(f"{folder_path}/min_max.yaml", "w") as f:
            yaml.dump(config, f)

    def load(self, folder_path):
        self.inputs = np.load(f"{folder_path}/inputs.npy")
        self.labels = np.load(f"{folder_path}/labels.npy")


# this is just an entrypoint to test the dataloader, we will not use this in the actual pipeline
if __name__ == "__main__":
    dataset = FluidDataset()
    current_file_path = Path(__file__).resolve()
    train_dir = current_file_path.parent.parent / "build" / "train"
    dataset.create(train_dir, n_samples=101)
    train_loader = DataLoader(dataset, batch_size=32, shuffle=True)
    sample_inputs, sample_labels = next(iter(train_loader))
    evaluate.visualize(sample_inputs, sample_labels)
