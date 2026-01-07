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
import prediction_plot


def init_model(submission_dir):
    """Inits the model.

    Args:
        submission_dir: Path to the submission directory.
    """
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

    # first and last samples in the training set
    training_samples = [
        trainer.train_loader.dataset[:1], trainer.train_loader.dataset[-1:]
    ]

    for training_sample in training_samples:
        prediction_plot.interpolate_flow_speeds(
            model,
            plot_dir,
            submission_dir,
            training_sample
        )


def generate_extrapolation_plots(model, save_path, stats):
    """Saves visualizations of model predictions for extrapolation scenarios.

    Args:
        model: The trained model for making predictions.
        save_path: The directory to save the submission to.
        stats: Normalization statistics for inputs and labels.
    """
    plot_dir = Path(save_path) / "plots" / "extrapolation"
    plot_dir.mkdir(parents=True, exist_ok=True)
    resource_dir = Path(__file__).resolve().parent.parent / RESOURCES

    prediction_plot.extrapolate_flow_speeds(
        model,
        plot_dir,
        resource_dir,
        stats,
        [0.25, 3.0]
    )
    prediction_plot.extrapolate_resolution(
        model,
        plot_dir,
        resource_dir,
        stats,
        [(41, 41)]
    )
    prediction_plot.extrapolate_boundary(model, plot_dir, stats)


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
