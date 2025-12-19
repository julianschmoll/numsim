# some tests on loading the data

import numpy as np
import pyvista as pv
import os
from pathlib import Path

def prepare_dataset(base_dir, n_samples):
    all_inputs = []
    all_labels = []

    base_path = Path(base_dir)

    for i in range(n_samples):
        folder = base_path / f"out_{i:03d}"
        if not folder.exists():
            print(f"Skipping: {folder} does not exist.")
            continue
        vti_files = sorted([f for f in os.listdir(folder) if f.startswith('output_') and f.endswith('.vti')])
        if not vti_files:
            print(f"Skipping: {folder} does not contain any VTI files.")
            continue

        print(f"Processing sample {i:03d}...")
        mesh = pv.read(folder / vti_files[-1])
        nx, ny, _ = mesh.dimensions

        vel_data = mesh.point_data['velocity'].reshape((ny, nx, 3), order='F')

        u_field = vel_data[:, :, 0] # Index 0 is u
        v_field = vel_data[:, :, 1] # Index 1 is v

        ux_val = u_field[-1, 1:-1].mean()
        input_channel = np.zeros((1, ny, nx))
        input_channel[0, -1, 1:-1] = ux_val
        label_channels = np.stack([u_field, v_field], axis=0)

        all_inputs.append(input_channel)
        all_labels.append(label_channels)

    return np.array(all_inputs), np.array(all_labels)

def normalize_and_save_stats(data):
    stats = []
    for c in range(data.shape[1]):
        c_min = data[:, c].min()
        c_max = data[:, c].max()
        if c_max - c_min > 1e-9:
            data[:, c] = (data[:, c] - c_min) / (c_max - c_min)
        stats.append({'min': float(c_min), 'max': float(c_max)})
    return data, stats


train_dir = "../build/train"
inputs, labels = prepare_dataset(train_dir, 101)

inputs, in_stats = normalize_and_save_stats(inputs)
labels, lb_stats = normalize_and_save_stats(labels)

config_content = f"""
inputs:
  u:
    max: {in_stats[0]['max']}
    min: {in_stats[0]['min']}
labels:
  u:
    max: {lb_stats[0]['max']}
    min: {lb_stats[0]['min']}
  v:
    max: {lb_stats[1]['max']}
    min: {lb_stats[1]['min']}
"""

with open('min_max.yaml', 'w') as f:
    f.write(config_content.strip())

np.save('inputs.npy', inputs.astype(np.float32))
np.save('labels.npy', labels.astype(np.float32))

print("\nDone!")
print(f"Final Input Shape:  {inputs.shape}")  # (Samples, 1, H, W)
print(f"Final Label Shape:  {labels.shape}")  # (Samples, 2, H, W)
print("Files saved: inputs.npy, labels.npy, min_max.yaml")
