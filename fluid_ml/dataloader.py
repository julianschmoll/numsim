import torch
import numpy as np
import pyvista as pv
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
        x = torch.from_numpy(self.inputs[idx])
        y = torch.from_numpy(self.labels[idx])
        return x, y

    def __str__(self):
        return f"FluidDataset with {self.__len__()} inputs/labels."

    def create(self, base_dir):
        inputs = []
        labels = []
        base_path = Path(base_dir)

        folders = sorted(
            path for path in base_path.iterdir()
            if path.is_dir() and path.name.startswith("out_")
        )

        for folder in folders:
            vti_files = sorted(
                file for file in folder.iterdir()
                if file.is_file() and file.name.startswith("output_") and file.suffix == ".vti"
            )

            for vti_file in vti_files:
                mesh = pv.read(folder / vti_file)
                hx, hy, _ = mesh.dimensions
                vel_data = mesh.point_data["velocity"].reshape((hx, hy, 3), order="C")

                u_field = vel_data[:, :, 0]
                v_field = vel_data[:, :, 1]

                # Handcrafted input as described in readthedocs
                ux_val = u_field[-1, 1:-1].mean()
                input_channel = np.zeros((1, hy, hx))
                input_channel[0, -1, 1:-1] = ux_val
                label_channels = np.stack([u_field, v_field], axis=0)

                inputs.append(input_channel)
                labels.append(label_channels)

        self.inputs = np.array(inputs)
        self.labels = np.array(labels)

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
        folder = Path(folder_path)
        folder.mkdir(parents=True, exist_ok=True)

        np.save(folder / "inputs.npy", self.inputs)
        np.save(folder / "labels.npy", self.labels)

        # I think the normalization here is wrong, we should already
        # normalize the data in set probably
        in_stats, lb_stats = self.normalize()
        config = {
            "inputs": {"u": in_stats[0]},
            "labels": {"u": lb_stats[0], "v": lb_stats[1]},
        }

        with open(folder / "min_max.yaml", "w") as f:
            yaml.dump(config, f)

    def load(self, folder_path):
        folder = Path(folder_path)
        self.inputs = np.load(folder / "inputs.npy")
        self.labels = np.load(folder / "labels.npy")


# this is just an entrypoint to test the dataloader, we will not use this in the actual pipeline
if __name__ == "__main__":
    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"

    dataset = FluidDataset()
    dataset.create(train_files_path)
    print(dataset)

    train_loader = DataLoader(dataset, batch_size=32, shuffle=False)
    evaluate.visualize(*next(iter(train_loader)))
