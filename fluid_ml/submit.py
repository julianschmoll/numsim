"""Generates submission files from trained models."""
import shutil
import sys
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import yaml

from constants import *  # noqa: F403, WPS347
import normalization
from visualize import visualize


def generate_submission(submission_dir, inputs_file):
    """Generate a submission file from a trained model and inputs.

    Args:
        inputs_file: The file containing the input tensors.

    """
    submission_dir = Path(submission_dir)
    with open(submission_dir / MIN_MAX_YAML, "r", encoding="utf-8") as min_max_yaml:
        min_max_stats = yaml.safe_load(min_max_yaml)

    sys.path.append(str(submission_dir))
    from mymodel import init_my_model  # pylint: disable=import-outside-toplevel

    model = init_my_model()
    model.load_state_dict(
        torch.load(submission_dir / "model.pt", map_location="cpu")
    )
    model.eval()

    write_csv(
        normalization.rescale_inputs(
            torch.load(inputs_file), min_max_stats
        ).to(torch.float32),
        model,
        submission_dir,
    )

    save_plots(model, submission_dir, min_max_stats)

    run_notebook(copy_notebook(submission_dir))


def save_plots(model, submission_dir: Path, min_max_stats):
    """Saves visualizations of model predictions at different flow speeds.

    Args:
        model: The trained model for making predictions.
        submission_dir: The directory to save the visualizations.
        min_max_stats: Stats for rescaling
    """
    ground_truth_inputs = torch.from_numpy(
        np.load(submission_dir / "np_arrays" / "inputs.npy")
    )
    ground_truth_labels = torch.from_numpy(
        np.load(submission_dir / "np_arrays" / "labels.npy")
    )
    print("GT inputs:", ground_truth_inputs.shape)
    print("\t", ground_truth_inputs.min(), ground_truth_inputs.max())
    print("GT labels:", ground_truth_labels.shape)
    print("\t", ground_truth_labels.min(), ground_truth_labels.max())
    den_gt_inputs, den_gt_labels = normalization.denormalize(
        ground_truth_inputs, ground_truth_labels,
        submission_dir / MIN_MAX_YAML
    )
    print("Den GT inputs:", den_gt_inputs.shape)
    print("\t", den_gt_inputs.min(), den_gt_inputs.max())
    print("Den GT labels:", den_gt_labels.shape)
    print("\t", den_gt_labels.min(), den_gt_labels.max())

    flow_speeds = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0]
    for flow_speed in flow_speeds:
        inputs, labels = normalization.denormalize(
            *model.predict(flow_speed, IMG_SIZE, IMG_SIZE, min_max_stats),
            submission_dir / MIN_MAX_YAML
        )
        save_visualization(
            inputs, labels,
            submission_dir,
            title=f"Prediction with Flow Speed: {flow_speed}"
        )
        if 0.5 <= flow_speed <= 1.5:
            sample_index = int((flow_speed - 0.5) * 100)
            print(flow_speed, sample_index)
            save_error_visualization(
                labels,
                den_gt_labels[sample_index:sample_index + 1],
                submission_dir,
                title=f"Error with Flow Speed: {flow_speed}"
            )
        else:
            pass
            # ToDo: Plot ground truth error plots from outside training data
            #  (generate it first)


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
            prediction = model(input_tensor.unsqueeze(0)).detach()
            rows.append(
                np.concatenate(([index], prediction.numpy().reshape(-1)))
            )

    cols = ["id"] + [
        f"val{value_index}"
        for value_index in range(N_CHANNELS * IMG_SIZE * IMG_SIZE)
    ]
    pd.DataFrame(rows, columns=cols).to_csv(
        submission_dir / "submission.csv", index=False
    )


def save_visualization(input_tensor, prediction, submission_dir, title=""):
    """Visualizes the prediction.

    Args:
        input_tensor: The input tensor.
        prediction: The prediction.
        submission_dir: Path to the submission directory.
        title: Title for the visualization.
    """
    sanitized_title = title.replace(
        " ", "_"
    ).replace(
        ":", ""
    ).lower() if title else "prediction_visualization"

    fig = visualize(
        [
            (input_tensor[0, 0], "Input"),
            (prediction[0, 0], "Prediction U"),
            (prediction[0, 1], "Prediction V")
        ],
        title=title if title else "Prediction Visualization",
    )
    fig.savefig(submission_dir / f"{sanitized_title}.png")

    quiver_fig = visualize(
        [(input_tensor[0], "Input"), (prediction[0], "Prediction")],
        title=title if title else "Prediction Visualization (Quiver)",
        plt_fn="quiver",
    )
    quiver_fig.savefig(submission_dir / f"{sanitized_title}_quiver.png")


def save_error_visualization(prediction, ground_truth, submission_dir, title=""):
    """Visualizes the error of the prediction.

    Args:
        prediction: The prediction of the flow.
        ground_truth: The simulated result of the flow.
        submission_dir: Path to the submission directory.
        title: Title for the visualization.
    """
    sanitized_title = title.replace(
        " ", "_"
    ).replace(
        ":", ""
    ).lower() if title else "error_visualization"

    prediction_u, prediction_v = prediction[0]
    prediction_magnitude = (prediction_u ** 2 + prediction_v ** 2) ** 0.5
    ground_truth_u, ground_truth_v = ground_truth[0]
    ground_truth_magnitude = (ground_truth_u ** 2 + ground_truth_v ** 2) ** 0.5

    error_u = abs(ground_truth_u - prediction_u)
    error_v = abs(ground_truth_v - prediction_v)
    error_magnitude = (error_u ** 2 + error_v ** 2) ** 0.5

    fig = visualize(
        [
            (ground_truth_u, "Ground Truth U"),
            (ground_truth_v, "Ground Truth V"),
            (ground_truth_magnitude, "Ground Truth Magnitude"),
            (prediction_u, "Prediction U"),
            (prediction_v, "Prediction V"),
            (prediction_magnitude, "Prediction Magnitude"),
            (error_u, "Absolute Error U"),
            (error_v, "Absolute Error V"),
            (error_magnitude, "Root Squared Error Magnitude"),
        ],
        title=title if title else "Error Visualization",
        num_rows=3
    )
    fig.savefig(submission_dir / f"{sanitized_title}.png")


def copy_notebook(destination_dir):
    """Copies result jupyter notebook to the destination directory.

    Args:
        destination_dir: The directory to copy the notebook to.

    Returns:
        Path to the copied notebook.
    """
    source_notebook = (Path(__file__).resolve().parent.parent
                       / "resources" / "results.ipynb")
    destination_notebook = Path(destination_dir) / "results.ipynb"
    shutil.copy(source_notebook, destination_notebook)
    return destination_notebook


def run_notebook(notebook_path, save_pdf=True):
    """Runs a Jupyter notebook and optionally saves it as a PDF.

    Args:
        notebook_path: Path to the Jupyter notebook.
        save_pdf: Whether to save the executed notebook as a PDF.
    """
    run_notebook_cmd = [
        "jupyter", "nbconvert",
        "--to", "notebook",
        "--execute",
        "--inplace",
        notebook_path
    ]
    try:
        subprocess.run(
            run_notebook_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Error executing notebook {notebook_path}: {error}"
        ) from error

    if not save_pdf:
        return

    pdf_cmd = ["jupyter", "nbconvert", "--to", "pdf", "--no-input", notebook_path]

    # this requires a LaTeX installation
    try:
        subprocess.run(
            pdf_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Error converting notebook {notebook_path} to PDF: {error}"
        ) from error


if __name__ == "__main__":
    models_dir = Path(__file__).resolve().parent.parent / "models"
    inputs_path = Path(__file__).resolve().parent.parent / "resources" / "inputs.pt"
    model_dirs = [directory for directory in models_dir.iterdir() if directory.is_dir()]

    if not model_dirs:
        raise RuntimeError("No models found in models directory.")

    latest_submission = max(model_dirs, key=lambda path: path.stat().st_mtime)

    generate_submission(latest_submission, inputs_path)
