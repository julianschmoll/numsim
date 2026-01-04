"""Evaluation utilities for fluid dynamics models."""
import numpy as np
import torch

import normalization
from constants import *  # noqa: F403, WPS347


def rmse(pred, truth):
    """Root Mean Square Error calculation.

    Args:
        pred (np.ndarray): Predicted values.
        truth (np.ndarray): True values.

    Returns:
        float: Root Mean Square Error between predictions and true values.
    """
    return float(np.sqrt(np.mean((pred - truth) ** 2)))


def mae(pred, truth):
    """Mean Absolute Error calculation.

    Args:
        pred (np.ndarray): Predicted values.
        truth (np.ndarray): True values.

    Returns:
        float: Mean Absolute Error between predictions and true values.
    """
    return float(np.mean(np.abs(pred - truth)))


# linter says we have too many variables but splitting function makes code less readable
def evaluation_metrics(loader, model, dataset, device):  # noqa: WPS210
    """Calculates RMSE and MAE on un-normalized data for a given loader.

    Args:
        loader (DataLoader): The data loader to evaluate on.
        model (torch.nn.Module): The trained model.
        dataset (FluidDataset): The dataset containing normalization stats.
        device (torch.device): The device to run the evaluation on.

    Returns:
        dict: A dictionary containing the calculated metrics.
    """
    labels, predictions = _get_predictions(device, loader, model)

    # we ignore slicing warnings as I did not find another working way to slice
    predicted_u_field = predictions[:, 0].copy()  # noqa: WPS478
    predicted_v_field = predictions[:, 1].copy()  # noqa: WPS478
    true_u_field = labels[:, 0].copy()  # noqa: WPS478
    true_v_field = labels[:, 1].copy()  # noqa: WPS478

    normalization.denormalize_channel(predicted_u_field, dataset.stats[LABELS][U])
    normalization.denormalize_channel(predicted_v_field, dataset.stats[LABELS][V])
    normalization.denormalize_channel(true_u_field, dataset.stats[LABELS][U])
    normalization.denormalize_channel(true_v_field, dataset.stats[LABELS][V])

    return {
        "u": {
            "rmse": rmse(predicted_u_field, true_u_field),
            "mae": mae(predicted_u_field, true_u_field)
        },
        "v": {
            "rmse": rmse(predicted_v_field, true_v_field),
            "mae": mae(predicted_v_field, true_v_field)
        },
    }


def _get_predictions(device, loader, model):
    """Get model predictions and true labels from a data loader.

    Args:
        device: The device to run the model on.
        loader: The data loader to get predictions from.
        model: The trained model.

    Returns:
        A tuple of numpy arrays: (true labels, model predictions).
    """
    model.eval()
    all_preds = []
    all_labels = []

    with torch.no_grad():
        for inputs, labels in loader:
            inputs = inputs.to(device)
            outputs = model(inputs)
            all_preds.append(outputs.cpu().numpy())
            all_labels.append(labels.cpu().numpy())

    return np.concatenate(all_labels), np.concatenate(all_preds)
