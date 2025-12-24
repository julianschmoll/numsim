from pathlib import Path

import numpy as np
import pyvista as pv
import torch
import yaml
from torch.utils.data import DataLoader, Dataset

import evaluate


class FluidDataset(Dataset):
    inputs: np.ndarray # shape: (#samples, #fields, #cells x, #cells y)
    labels: np.ndarray
    stats:  dict[str, dict[str, dict[str, float]]]
    normalized: bool


    def __init__(self):
        self.inputs = np.array([], dtype=np.float32)
        self.labels = np.array([], dtype=np.float32)
        self.stats = {}
        self.normalized = False


    def __len__(self):
        return len(self.inputs)


    def __getitem__(self, idx) -> tuple[torch.Tensor, torch.Tensor]:
        x = torch.from_numpy(self.inputs[idx])
        y = torch.from_numpy(self.labels[idx])
        return x, y


    def __str__(self):
        string = f"FluidDataset with {self.__len__()} inputs/labels.\n"
        string += "Dataset Statistics (Normalized 0-1):\n"
        for category, channels in self.stats.items():
            string += f"  {category}:\n"
            for channel, values in channels.items():
                string += f"    - {channel}: min={values['min']:.4f}, max={values['max']:.4f}\n"
        return string


    def create(self, base_dir: Path | str):
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

        self.inputs = np.array(inputs, dtype=np.float32)
        self.labels = np.array(labels, dtype=np.float32)
        self.normalize()


    def normalize(self) -> None:
        """Normalize the data if it's not already normalized"""

        def process(data: np.ndarray):
            c_min, c_max = float(data.min()), float(data.max())
            if (c_max - c_min) > 1e-9:
                data[:] = (data - c_min) / (c_max - c_min)
            return {"min": c_min, "max": c_max}
        
        if self.normalized:
            return

        self.stats = {
            "inputs": {
                "u": process(self.inputs[:, 0])
            },
            "labels": {
                "u": process(self.labels[:, 0]),
                "v": process(self.labels[:, 1])
            }
        }
        self.normalized = True


    def denormalize(self) -> None:
        """Denormalize the dataset if it is normalized."""
        
        def process(data: np.ndarray, c_stats: dict[str, float]):
            c_min, c_max = c_stats["min"], c_stats["max"]
            if (c_max - c_min) > 1e-9:
                data[:] = data * (c_max - c_min) + c_min

        if not self.normalized:
            return

        stats_in_u = self.stats["inputs"]["u"]
        stats_label_u = self.stats["labels"]["u"]
        stats_label_v = self.stats["labels"]["v"]
        process(self.inputs[:, 0], stats_in_u)
        process(self.labels[:, 0], stats_label_u)
        process(self.labels[:, 1], stats_label_v)
        self.normalized = False


    def save(self, dataset_path: str | Path):
        """Save the dataset and normalization info to `dataset_path`."""

        folder = Path(dataset_path)
        folder.mkdir(parents=True, exist_ok=True)

        stats = { **self.stats, "normalized": self.normalized }

        np.save(folder / "inputs.npy", self.inputs)
        np.save(folder / "labels.npy", self.labels)

        with open(folder / "min_max.yaml", "w") as f:
            yaml.dump(stats, f)


    def load(self, dataset_path: str | Path):
        """Load a dataset from `dataset_path` without applying any normalization or denormalization."""

        folder = Path(dataset_path)

        self.inputs = np.load(folder / "inputs.npy")
        self.labels = np.load(folder / "labels.npy")

        with open(folder / "min_max.yaml", "r") as f:
            try:
                stats = yaml.safe_load(f)
                self.normalized = stats["normalized"]
                self.stats = {
                    "inputs": stats["inputs"],
                    "labels": stats["labels"]
                }
            except yaml.YAMLError as exception:
                print(f"Error while loading min_max.yaml:\n{exception}")


# this is just an entrypoint to test the dataloader, we will not use this in the actual pipeline
if __name__ == "__main__":
    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"

    dataset = FluidDataset()
    #dataset.create(train_files_path)
    dataset.load(current_file_path.parent)
    print(dataset)

    train_loader = DataLoader(dataset, batch_size=32, shuffle=False)
    evaluate.visualize(*next(iter(train_loader)))
