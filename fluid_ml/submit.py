import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import yaml

import constants

DEFAULT_OUT_RANGE = (torch.tensor([.0]), torch.tensor([1.0]))


def rescale(tensor, mins_in, maxs_in, out_range=DEFAULT_OUT_RANGE):
    """Rescale a tensor to a new range.

    Args:
        tensor: The input tensor to be rescaled.
        mins_in: The minimum values of the input tensor.
        maxs_in: The maximum values of the input tensor.
        out_range: A tuple containing the desired output range (min, max).

    Returns:
        The rescaled tensor.
    """
    min_out, max_out = out_range
    mins_in = mins_in.unsqueeze(-1).unsqueeze(-1)
    maxs_in = maxs_in.unsqueeze(-1).unsqueeze(-1)
    min_out = min_out.unsqueeze(-1).unsqueeze(-1)
    max_out = max_out.unsqueeze(-1).unsqueeze(-1)

    return (tensor - mins_in) / (maxs_in - mins_in) * (max_out - min_out) + min_out


def rescale_inputs(inputs, min_max_stats):
    """Rescale input tensors using provided statistics.

    Args:
        inputs: The input tensor to be rescaled.
        min_max_stats: A dictionary containing min and max statistics for rescaling.

    Returns:
        The rescaled input tensor.
    """
    inputs_mins = torch.tensor(
        [min_max_stats[constants.INPUTS_KEY][constants.U_CHANNEL_KEY]["min"]]
    )
    inputs_maxs = torch.tensor(
        [min_max_stats[constants.INPUTS_KEY][constants.U_CHANNEL_KEY]["max"]]
    )
    return rescale(inputs, inputs_mins, inputs_maxs)


def generate_submission(submission_dir, inputs_file):
    """Generate a submission file from a trained model and inputs.

    Args:
        submission_dir: The directory containing the trained model and statistics.
        inputs_file: The file containing the input tensors.

    """
    sys.path.append(str(submission_dir))
    from mymodel import init_my_model

    submission_dir = Path(submission_dir)

    with open(submission_dir / "min_max.yaml", "r") as min_max_yaml:
        min_max_stats = yaml.safe_load(min_max_yaml)

    model = init_my_model()
    model.load_state_dict(torch.load(submission_dir / "model.pt", map_location="cpu"))
    model.eval()

    inputs_scaled = rescale_inputs(torch.load(inputs_file), min_max_stats)

    write_csv(inputs_scaled.to(torch.float32), model, submission_dir)


def write_csv(inputs_scaled, model, submission_dir):
    """Write model predictions to a CSV file.

    Args:
        inputs_scaled: The scaled input tensors.
        model: The trained model for making predictions.
        submission_dir: The directory to save the submission CSV file.
    """
    rows = []
    with torch.no_grad():
        for index, input_tensor in enumerate(inputs_scaled):
            prediction = model(input_tensor.unsqueeze(0)).detach().numpy()
            rows.append(
                np.concatenate(([index], prediction.reshape(-1)))
            )

    cols = ["id"] + [
        f"val{value_index}"
        for value_index in range(
            constants.N_CHANNELS * constants.IMG_SIZE * constants.IMG_SIZE
        )
    ]
    pd.DataFrame(rows, columns=cols).to_csv(
        submission_dir / "submission.csv", index=False
    )


if __name__ == "__main__":
    models_dir = Path(__file__).resolve().parent.parent / "models"
    inputs_path = Path(__file__).resolve().parent.parent / "resources" / "inputs.pt"
    model_dirs = [directory for directory in models_dir.iterdir() if directory.is_dir()]

    if not model_dirs:
        raise RuntimeError("No models found in models directory.")

    latest_submission = max(model_dirs, key=lambda path: path.stat().st_mtime)
    generate_submission(latest_submission, inputs_path)
