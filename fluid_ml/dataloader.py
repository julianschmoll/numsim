"""DataLoader for fluid simulation data stored in VTI files."""
from pathlib import Path

import numpy as np
import pyvista as pv
import torch
import yaml
from torch.utils.data import Dataset

import constants
import normalization


def _get_input_channel(hx, hy, u_field):
    """Crafts input channel based on mean value in u field.

    Args:
        hx: Width of the Data.
        hy: Height of the Data.
        u_field: U Velocity field.

    Returns:
        input_channel: Input channel corresponding to u_field.
    """
    ux_val = u_field[-1, 1:-1].mean()
    input_channel = np.zeros((1, hy, hx), dtype=np.float32)
    input_channel[0, -1, 1:-1] = ux_val
    return input_channel


def _read_vti_file(path):
    """Reads vti files with pyvista.

    Args:
        path: Path to the vti file.

    Returns:
        dict: VTI Data.
    """
    mesh = pv.read(path)
    hx, hy, _ = mesh.dimensions
    velocity_data = mesh.point_data["velocity"].reshape((hx, hy, 3), order="C")

    return {
        "hx": hx,
        "hy": hy,
        "u_field": velocity_data[..., 0],
        "v_field": velocity_data[..., 1],
    }


def _load_data(base_dir):
    """Generator for loading data from vti files in directory.

    Args:
        base_dir: Base directory for data files.

    Yields:
        (input_channel, output_channel): Data from the current file.
    """
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
                _get_input_channel(vti_data["hx"], vti_data["hy"], vti_data["u_field"]),
                np.stack([vti_data["u_field"], vti_data["v_field"]], axis=0)
            )


class FluidDataset(Dataset):
    """Dataset class to load fluid data from vti files."""
    def __init__(self, base_dir):
        """Initializes the dataset.

        Args:
            base_dir (Path): Path to the data directory.
        """
        self.stats = {}
        self.normalized = False

        fluid_data = list(_load_data(base_dir))
        inputs, labels = zip(*fluid_data)

        self.inputs = np.array(inputs, dtype=np.float32)
        self.labels = np.array(labels, dtype=np.float32)

    def __len__(self):
        """Returns the length of the dataset.

        Returns:
            int: Length of the dataset.
        """
        return len(self.inputs)

    def __str__(self):
        """Return a string representation of the dataset including statistics.

        Returns:
            str: String representation of the dataset.
        """
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
        """Returns the input and label tensors.

        Args:
            idx (int): Index of the data sample.

        Returns:
            tuple[torch.Tensor, torch.Tensor]: Input and label tensors.
        """
        input_ = torch.from_numpy(self.inputs[idx])
        label = torch.from_numpy(self.labels[idx])
        return input_, label

    def normalize(self) -> None:
        """Normalize the data if it's not already normalized."""
        if self.normalized:
            return

        self.stats = {
            constants.INPUTS_KEY: {
                constants.U_CHANNEL_KEY:
                    normalization.normalize_channel(self.inputs[:, 0])  # noqa: WPS478
            },
            constants.LABELS_KEY: {
                constants.U_CHANNEL_KEY:
                    normalization.normalize_channel(self.labels[:, 0]),  # noqa: WPS478
                constants.V_CHANNEL_KEY:
                    normalization.normalize_channel(self.labels[:, 1]),  # noqa: WPS478
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

        normalization.denormalize_channel(
            self.inputs[:, 0], stats_in_u  # noqa: WPS478
        )
        normalization.denormalize_channel(
            self.labels[:, 0], stats_label_u  # noqa: WPS478
        )
        normalization.denormalize_channel(
            self.labels[:, 1], stats_label_v  # noqa: WPS478
        )

        self.normalized = False

    def save(self, dataset_path: str | Path):
        """Save the dataset and normalization info to `dataset_path`.

        Args:
            dataset_path (Path): Folder to save the dataset to.

        Returns:
            Path: Path to the saved dataset statistics.
        """
        folder = Path(dataset_path)
        folder.mkdir(parents=True, exist_ok=True)

        save_path = folder / "min_max.yaml"

        with open(save_path, "w", encoding="utf-8") as min_max_yaml:
            yaml.dump(self.stats, min_max_yaml)

        return save_path
