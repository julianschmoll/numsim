from pathlib import Path

import numpy as np
import pyvista as pv
import torch
import yaml
from torch.utils.data import DataLoader, Dataset

import evaluate
import constants


def _normalize_channel(channel_data):
    """Normalize a single data channel in-place and return its stats."""
    c_min = float(channel_data.min())
    c_max = float(channel_data.max())
    if (c_max - c_min) > constants.EPSILON:
        channel_data[...] = (channel_data - c_min) / (c_max - c_min)
    return {"min": c_min, "max": c_max}


def _denormalize_channel(channel_data, channel_stats):
    """Denormalize a single data channel in-place using provided stats."""
    c_min, c_max = channel_stats["min"], channel_stats["max"]
    if (c_max - c_min) > constants.EPSILON:
        channel_data[...] = channel_data * (c_max - c_min) + c_min


def _get_input_channel(hx, hy, u_field):
    ux_val = u_field[-1, 1:-1].mean()
    input_channel = np.zeros((1, hy, hx))
    input_channel[0, -1, 1:-1] = ux_val
    return input_channel


def _read_vti_file(path):
    mesh = pv.read(path)
    dims = mesh.dimensions[0], mesh.dimensions[1]
    velocity_data = mesh.point_data["velocity"].reshape(
        (*dims, 3), order="C"
    )
    return {
        "dimensions": dims,
        "u_field": velocity_data[..., 0],
        "v_field": velocity_data[..., 1],
    }


def _load_data(base_dir):
    """Helper generator to keep __init__ clean."""
    output_folders = sorted(
        path
        for path in Path(base_dir).iterdir()
        if path.is_dir() and path.name.startswith("out_")
    )
    for folder in output_folders:
        vti_files = sorted(
            vti_file
            for vti_file in folder.iterdir()
            if vti_file.is_file()
            and vti_file.name.startswith("output_")
            and vti_file.suffix == ".vti"
        )
        for vti_file in vti_files:
            vti_data = _read_vti_file(folder / vti_file)
            yield (
                _get_input_channel(*vti_data["dimensions"], vti_data["u_field"]),
                np.stack([vti_data["u_field"], vti_data["v_field"]], axis=0)
            )


class FluidDataset(Dataset):
    def __init__(self, base_dir):
        self.stats = {}
        self.normalized = False

        fluid_data = list(_load_data(base_dir))
        inputs, labels = zip(*fluid_data)

        self.inputs = np.array(inputs, dtype=np.float32)
        self.labels = np.array(labels, dtype=np.float32)

    def __len__(self):
        return len(self.inputs)

    def __str__(self):
        """Return a string representation of the dataset including statistics."""
        lines = [
            f"FluidDataset with {len(self)} inputs/labels.",
            "Dataset Statistics (Normalized 0-1):",
        ]
        for category, channels in self.stats.items():
            lines.append(f"  {category}:")
            for channel, channel_stats in channels.items():
                lines.append(f"    - {channel}: "
                             f"min={channel_stats["min"]:.4f}, "
                             f"max={channel_stats["max"]:.4f}")
        return "\n".join(lines)

    def __getitem__(self, idx) -> tuple[torch.Tensor, torch.Tensor]:
        input_ = torch.from_numpy(self.inputs[idx])
        label = torch.from_numpy(self.labels[idx])
        return input_, label

    def normalize(self) -> None:
        """Normalize the data if it's not already normalized."""
        if self.normalized:
            return

        self.stats = {
            constants.INPUTS_KEY: {
                constants.U_CHANNEL_KEY: _normalize_channel(self.inputs[..., 0])
            },
            constants.LABELS_KEY: {
                constants.U_CHANNEL_KEY: _normalize_channel(self.labels[..., 0]),
                constants.V_CHANNEL_KEY: _normalize_channel(self.labels[..., 1]),
            },
        }
        self.normalized = True

    def denormalize(self) -> None:
        """Denormalize the dataset if it is normalized."""
        if not self.normalized:
            return

        stats_in_u = self.stats[constants.INPUTS_KEY][constants.U_CHANNEL_KEY]
        stats_label_u = self.stats[constants.LABELS_KEY][constants.U_CHANNEL_KEY]
        stats_label_v = self.stats[constants.LABELS_KEY][constants.V_CHANNEL_KEY]

        _denormalize_channel(self.inputs[..., 0], stats_in_u)
        _denormalize_channel(self.labels[..., 0], stats_label_u)
        _denormalize_channel(self.labels[..., 1], stats_label_v)

        self.normalized = False

    def save(self, dataset_path: str | Path):
        """Save the dataset and normalization info to `dataset_path`."""

        folder = Path(dataset_path)
        folder.mkdir(parents=True, exist_ok=True)

        with open(folder / "min_max.yaml", "w") as min_max_yaml:
            yaml.dump(self.stats, min_max_yaml)


if __name__ == "__main__":
    current_file_path = Path(__file__).resolve()
    train_files_path = current_file_path.parent.parent / "build" / "train"

    dataset = FluidDataset(train_files_path)
    dataset.save(current_file_path.parent)

    train_loader = DataLoader(
        dataset, batch_size=constants.DEFAULT_BATCH_SIZE, shuffle=False
    )
    evaluate.visualize(*next(iter(train_loader)))
