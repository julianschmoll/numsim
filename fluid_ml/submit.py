import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import yaml

import constants
import normalization
from visualize import visualize


def generate_submission(submission_dir, inputs_file):
    """Generate a submission file from a trained model and inputs.

    Args:
        submission_dir: The directory containing the trained model and statistics.
        inputs_file: The file containing the input tensors.

    """
    sys.path.append(str(submission_dir))
    from mymodel import init_my_model

    submission_dir = Path(submission_dir)
    stats_path = submission_dir / "min_max.yaml"

    with open(stats_path, "r") as min_max_yaml:
        min_max_stats = yaml.safe_load(min_max_yaml)

    model = init_my_model()
    model.load_state_dict(torch.load(submission_dir / "model.pt", map_location="cpu"))
    model.eval()

    write_csv(
        normalization.rescale_inputs(
            torch.load(inputs_file), min_max_stats
        ).to(torch.float32),
        model,
        submission_dir,
        stats_path
    )


def write_csv(inputs_scaled, model, submission_dir, stats_path):
    """Write model predictions to a CSV file.

    Args:
        inputs_scaled: The scaled input tensors.
        model: The trained model for making predictions.
        submission_dir: The directory to save the submission CSV file.
        stats_path: Path to the statistics file for visualization.
    """
    rows = []
    with torch.no_grad():
        for index, input_tensor in enumerate(inputs_scaled):
            prediction = model(input_tensor.unsqueeze(0)).detach()
            rows.append(
                np.concatenate(([index], prediction.numpy().reshape(-1)))
            )

            if index == len(inputs_scaled) - 1:
                visualize_prediction(input_tensor, prediction, stats_path)

    cols = ["id"] + [
        f"val{value_index}"
        for value_index in range(
            constants.N_CHANNELS * constants.IMG_SIZE * constants.IMG_SIZE
        )
    ]
    pd.DataFrame(rows, columns=cols).to_csv(
        submission_dir / "submission.csv", index=False
    )


def visualize_prediction(input_tensor, prediction, stats_path):
    """Visualizes the prediction.

    Args:
        input_tensor: The input tensor.
        prediction: The prediction.
        stats_path: Path to the statistics file for visualization.
    """
    visualize(
        input_tensor.unsqueeze(0),
        prediction,
        title="Prediction Visualization",
        stats_path=stats_path,
    )
    visualize(
        input_tensor.unsqueeze(0),
        prediction,
        title="Prediction Visualization (Quiver)",
        quiver=True,
        stats_path=stats_path,
    )


if __name__ == "__main__":
    models_dir = Path(__file__).resolve().parent.parent / "models"
    inputs_path = Path(__file__).resolve().parent.parent / "resources" / "inputs.pt"
    model_dirs = [directory for directory in models_dir.iterdir() if directory.is_dir()]

    if not model_dirs:
        raise RuntimeError("No models found in models directory.")

    latest_submission = max(model_dirs, key=lambda path: path.stat().st_mtime)
    generate_submission(latest_submission, inputs_path)
