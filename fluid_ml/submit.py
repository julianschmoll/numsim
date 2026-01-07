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
import visualize


def init_model(submission_dir):
    sys.path.append(str(submission_dir))
    from mymodel import init_my_model  # pylint: disable=import-outside-toplevel

    submission_dir = Path(submission_dir)

    model = init_my_model()
    model.load_state_dict(torch.load(submission_dir / "model.pt", map_location="cpu"))
    model.eval()

    return model

def generate_kaggle_submission(model, submission_dir, inputs_file):
    """Generate a submission file from a trained model and inputs.

    Args:
        model (torch.nn.Module): The trained model.
        submission_dir: The directory containing the trained model and statistics.
        inputs_file: The file containing the input tensors.

    """
    with open(submission_dir / MIN_MAX_YAML, "r", encoding="utf-8") as min_max_yaml:
        min_max_stats = yaml.safe_load(min_max_yaml)

    _write_csv(
        normalization.rescale_inputs(
            torch.load(inputs_file), min_max_stats
        ).to(torch.float32),
        model,
        submission_dir,
    )


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


def generate_interpolation_plots(model, trainer, submission_dir):
    """Saves visualizations of model predictions at different flow speeds.

    Args:
        model: The trained model for making predictions.
        trainer: The trainer containing the trained model and data loaders.
        submission_dir: The directory to save the submission to.
    """
    plot_dir = submission_dir / "plots" / "interpolation"
    plot_dir.mkdir(parents=True, exist_ok=True)

    train_loader = trainer.train_loader
    first_train_sample = train_loader.dataset[0:1]
    last_train_sample = train_loader.dataset[-1:]
    training_samples = [first_train_sample, last_train_sample]

    for training_sample in training_samples:
        training_input = training_sample[0]
        training_label = training_sample[1]
        den_training_input, den_training_label = normalization.denormalize(
            training_input, training_label,
            submission_dir / MIN_MAX_YAML
        )

        flow_speed = round(float(den_training_input[0, 0, -1, 1:-1].mean()), 2)

        _, den_pred_label = normalization.denormalize(
            training_input,
            model(training_input).detach(),
            submission_dir / MIN_MAX_YAML
        )
        visualize.prediction_visualization(
            den_training_input,
            den_pred_label,
            plot_dir,
            title=f"Prediction with Flow Speed: {flow_speed}"
        )
        visualize.error_visualization(
            den_pred_label,
            den_training_label,
            plot_dir,
            title=f"Error with Flow Speed: {flow_speed}"
        )


def generate_extrapolation_plots(model, plot_dir, resource_dir, stats):
    _extrapolate_flow_seeds(model, plot_dir, resource_dir, stats, [0.25, 3.0])
    _extrapolate_resolution(model, plot_dir, resource_dir, stats, [(41, 41)])
    _extrapolate_boundary(model, plot_dir, stats)


def _extrapolate_boundary(model, plot_dir, stats):
    # from here we don't have ground truth data,
    # as we could simply flip/rotate the tensor
    # to make those interpolation problems

    # extrapolation with negative flow speed
    input_tensor = _tensor_from_flow_speed(-1.0, IMG_SIZE, IMG_SIZE, stats)
    neg_speed_prediction = model(input_tensor).detach()
    visualize.prediction_visualization(
        input_tensor, neg_speed_prediction, plot_dir, f"Negative Flow Speed"
    )

    # extrapolation with boundary at different position
    flow_speed = 1.0
    input_channel = np.zeros((1, IMG_SIZE, IMG_SIZE), dtype=np.float32)
    input_channel[0, 1:-1, -1] = flow_speed
    input_tensor = torch.from_numpy(input_channel).unsqueeze(0)
    out_min = torch.tensor(stats[INPUTS][U][MIN])
    out_max = torch.tensor(stats[INPUTS][U][MAX])
    n_input_tensor = normalization.rescale(input_tensor, out_min, out_max)
    other_boundary_prediction = normalization.rescale(
        model(n_input_tensor).detach(), torch.tensor(.0), torch.tensor(.1),
        out_range=(out_min, out_max)
    )

    visualize.prediction_visualization(
        input_tensor, other_boundary_prediction, plot_dir, f"Right Boundary"
    )


def _extrapolate_resolution(model, plot_dir, resource_dir, stats, resolutions):
    for hx, hy in resolutions:
        input_tensor = _tensor_from_flow_speed(1.0, hx, hy, stats)
        truth = torch.from_numpy(
            np.load(Path(resource_dir) / f"{hx}_cells_labels.npy")
        )
        _, prediction = normalization.denormalize(
            input_tensor,
            model(input_tensor).detach(),
            Path(plot_dir).parent.parent / MIN_MAX_YAML
        )

        visualize.error_visualization(
            prediction, truth, plot_dir, f"Extrapolation at {hx}x{hy} Resolution"
        )


def _extrapolate_flow_seeds(model, plot_dir, resource_dir, stats, extrapolation_speeds):
    for flow_speed in extrapolation_speeds:
        input_tensor = _tensor_from_flow_speed(flow_speed, IMG_SIZE, IMG_SIZE, stats)
        truth = torch.from_numpy(
            np.load(Path(resource_dir) / f"flow_speed_{flow_speed}_labels.npy")
        )
        _, prediction = normalization.denormalize(
            input_tensor,
            model(input_tensor).detach(),
            Path(plot_dir).parent.parent / MIN_MAX_YAML
        )

        visualize.error_visualization(
            prediction, truth, plot_dir, f"Extrapolation at Flow Speed: {flow_speed}"
        )


def _write_csv(inputs_scaled, model, submission_dir):
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


def _tensor_from_flow_speed(flow_speed, hx, hy, min_max_stats):
    """Creates an input tensor from a flow speed value.

    Args:
            flow_speed: Speed of flow at lid.
            hx: Width of the domain.
            hy: Height of the domain.

    Returns:
        Normalized input tensor.
    """
    input_channel = np.zeros((1, hx, hy), dtype=np.float32)
    input_channel[0, -1, 1:-1] = flow_speed
    input_tensor = torch.from_numpy(input_channel).unsqueeze(0)
    out_min = torch.tensor(min_max_stats[INPUTS][U][MIN])
    out_max = torch.tensor(min_max_stats[INPUTS][U][MAX])
    n_input_tensor = normalization.rescale(input_tensor, out_min, out_max) # scale to 0-1
    return n_input_tensor


if __name__ == "__main__":
    models_dir = Path(__file__).resolve().parent.parent / "models"
    inputs_path = Path(__file__).resolve().parent.parent / "resources" / "inputs.pt"
    model_dirs = [directory for directory in models_dir.iterdir() if directory.is_dir()]

    if not model_dirs:
        raise RuntimeError("No models found in models directory.")

    latest_submission = max(model_dirs, key=lambda path: path.stat().st_mtime)
    generate_submission(latest_submission, inputs_path)
