"""Plotting utilities for model predictions."""
from pathlib import Path

import numpy as np
import torch

import normalization
import visualize
from constants import *  # noqa: F403, WPS347


def interpolate_flow_speeds(model, plot_dir, submission_dir, training_sample):
    """Plots prediction and error viz of a training sample input.

    Args:
        model: The trained model.
        plot_dir: The directory to save the figure to.
        submission_dir: The directory of the submission.
        training_sample: The training sample.
    """
    den_training_input, den_training_label = normalization.denormalize(
        training_sample[0], training_sample[1],
        submission_dir / MIN_MAX_YAML
    )
    flow_speed = round(
        float(
            den_training_input[0, 0, -1, 1:-1].mean()
        ), 2
    )

    _, den_pred_label = normalization.denormalize(
        training_sample[0],
        model(training_sample[0]).detach(),
        submission_dir / MIN_MAX_YAML
    )
    plot_subdir = plot_dir / f"flow_speed_{flow_speed}"
    plot_subdir.mkdir(parents=True, exist_ok=True)

    _prediction_visualization(
        den_training_input,
        den_pred_label,
        plot_subdir,
        title=f"Prediction with Flow Speed: {flow_speed}"
    )
    _error_visualization(
        den_pred_label,
        den_training_label,
        plot_subdir,
        title=f"Error with Flow Speed: {flow_speed}"
    )


def extrapolate_flow_speeds(model, plot_dir, resource_dir, stats, speeds):
    """Plots prediction and error viz of a tensor generated from a list of
    extrapolated flow speeds.

    Args:
        model: The trained model.
        plot_dir: The directory to save the figure to.
        resource_dir: The directory of the ground truth flow data.
        stats: Normalization statistics for inputs and labels of the model.
        speeds: The list of extrapolated flow speeds.
    """

    for flow_speed in speeds:
        plt_subdir = plot_dir / f"flow_speed_{flow_speed}"
        plt_subdir.mkdir(parents=True, exist_ok=True)
        input_tensor = _tensor_from_flow_speed(flow_speed, IMG_SIZE, IMG_SIZE, stats)
        truth = torch.from_numpy(
            np.load(Path(resource_dir) / f"flow_speed_{flow_speed}_labels.npy")
        )
        _, prediction = normalization.denormalize(
            input_tensor,
            model(input_tensor).detach(),
            Path(plot_dir).parent.parent / MIN_MAX_YAML
        )

        _error_visualization(
            prediction, truth, plt_subdir, f"Extrapolation at Flow Speed: {flow_speed}"
        )


def extrapolate_resolution(model, plot_dir, resource_dir, stats, resolutions):
    """Plots prediction and error viz of a tensor generated from a list of resolutions.

    Args:
        model: The trained model.
        plot_dir: The directory to save the figure to.
        resource_dir: The directory of the ground truth flow data.
        stats: Normalization statistics for inputs and labels of the model.
        resolutions: The list of resolutions.
    """

    for resolution in resolutions:
        plt_subdir = plot_dir / f"{resolution[0]}x{resolution[1]}_resolution"
        plt_subdir.mkdir(parents=True, exist_ok=True)
        input_tensor = _tensor_from_flow_speed(1.0, *resolution, stats)

        truth = torch.from_numpy(
            np.load(
                Path(resource_dir) / f"{resolution[0]}x{resolution[1]}_cells_labels.npy"
            )
        )
        _, prediction = normalization.denormalize(
            input_tensor,
            model(input_tensor).detach(),
            Path(plot_dir).parent.parent / MIN_MAX_YAML
        )

        _error_visualization(
            prediction, truth, plt_subdir,
            f"Extrapolation at {resolution[0]}x{resolution[1]} Resolution"
        )


def extrapolate_boundary(model, plot_dir, stats):
    """Plots prediction and error viz of a tensor generated from different boundaries.

    Args:
        model: The trained model.
        plot_dir: The directory to save the figure to.
        stats: Normalization statistics for inputs and labels of the model.
    """

    # extrapolation with negative flow speed
    plot_subdir = plot_dir / "negative_flow_speed"
    plot_subdir.mkdir(parents=True, exist_ok=True)
    input_tensor = _tensor_from_flow_speed(-1.0, IMG_SIZE, IMG_SIZE, stats)
    _prediction_visualization(
        input_tensor, model(input_tensor).detach(), plot_subdir, "Negative Flow Speed"
    )

    # extrapolation with boundary at different position
    plot_subdir = plot_dir / "bottom_boundary"
    plot_subdir.mkdir(parents=True, exist_ok=True)
    input_channel = np.zeros((1, IMG_SIZE, IMG_SIZE), dtype=np.float32)
    input_channel[-1, 0:1, 1:-1] = 1.0
    input_tensor = torch.from_numpy(input_channel).unsqueeze(0)
    out_range = [
        torch.tensor(stats[INPUTS][U][MIN]),
        torch.tensor(stats[INPUTS][U][MAX])
    ]
    n_input_tensor = normalization.rescale(input_tensor, *out_range)
    _prediction_visualization(
        input_tensor,
        normalization.rescale(
            model(n_input_tensor).detach(), torch.tensor(.0), torch.tensor(.1),
            out_range=out_range
        ),
        plot_subdir,
        "Bottom Boundary"
    )


def _error_visualization(  # noqa: WPS210, as it makes the method more readable
        prediction, ground_truth, submission_dir, title=""
):
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

    fig = visualize.visualize(
        [
            (ground_truth_u, "Ground Truth U"),
            (ground_truth_v, "Ground Truth V"),
            (ground_truth_magnitude, "Ground Truth Magnitude"),
            (prediction_u, "Prediction U"),
            (prediction_v, "Prediction V"),
            (prediction_magnitude, "Prediction Magnitude"),
            (error_u, "Absolute Error U"),
            (error_v, "Absolute Error V"),
            (error_magnitude, "Combined Error Magnitude"),
        ],
        title=title if title else "Error Visualization",
        num_rows=3
    )
    fig.savefig(submission_dir / f"{sanitized_title}.png")


def _prediction_visualization(input_tensor, prediction, submission_dir, title=""):
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

    fig = visualize.visualize(
        [
            (input_tensor[0, 0], "Input"),
            (prediction[0, 0], "Prediction U"),
            (prediction[0, 1], "Prediction V")
        ],
        title=title if title else "Prediction Visualization",
    )
    fig.savefig(submission_dir / f"{sanitized_title}.png")

    quiver_fig = visualize.visualize(
        [(input_tensor[0], "Input"), (prediction[0], "Prediction")],
        title=title if title else "Prediction Visualization (Quiver)",
        plt_fn="quiver",
    )
    quiver_fig.savefig(submission_dir / f"{sanitized_title}_quiver.png")


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
    # scale to 0-1
    return normalization.rescale(input_tensor, out_min, out_max)
