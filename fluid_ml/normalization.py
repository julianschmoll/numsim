"""Normalization and denormalization utilities."""
import torch
import yaml

from constants import *  # noqa: F403, WPS347


DEFAULT_OUT_RANGE = (torch.tensor([.0]), torch.tensor([1.0]))


def denormalize(inputs, labels, stats_path):
    """Denormalizes the inputs and labels.

    Args:
        inputs: Input tensor
        labels: Labels tensor.
        stats_path: Path to statistics file.

    Returns:
        Denormalized tensor
    """
    with open(stats_path, "r", encoding="utf-8") as min_max_yaml:
        stats = yaml.safe_load(min_max_yaml)

    if not stats.get("inputs") and not stats.get("labels"):
        return inputs, labels

    return rescale_inputs(inputs, stats), rescale_labels(labels, stats)


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
    inputs_mins = torch.tensor([min_max_stats[INPUTS][U][MIN]])
    inputs_maxs = torch.tensor([min_max_stats[INPUTS][U][MAX]])
    return rescale(inputs, inputs_mins, inputs_maxs)


def rescale_labels(labels, stats):
    """Rescale label tensors using provided statistics.

    Args:
        labels: The label tensor to be rescaled.
        stats: Dictionary containing min and max statistics for rescaling.

    Returns:
        The rescaled label tensor.
    """
    labels_min = torch.tensor([
        stats[LABELS][U][MIN],
        stats[LABELS][V][MIN]
    ])
    labels_max = torch.tensor([
        stats[LABELS][U][MAX],
        stats[LABELS][V][MAX]
    ])
    return rescale(labels, *DEFAULT_OUT_RANGE, out_range=(labels_min, labels_max))


def normalize_channel(channel_data):
    """Normalize a single data channel in-place and return its stats.

    Args:
        channel_data: The data channel to be normalized.

    Returns:
        dict: Statistics (old min and max values) about normalization.
    """
    c_min = float(channel_data.min())
    c_max = float(channel_data.max())
    if (c_max - c_min) > EPSILON:
        channel_data[...] = (channel_data - c_min) / (c_max - c_min)
    return {MIN: c_min, MAX: c_max}


def denormalize_channel(channel_data, channel_stats):
    """Denormalize a single data channel in-place using provided stats.

    Args:
        channel_data: The data channel to be denormalized.
        channel_stats: Dictionary containing min, max statistics for denormalization.
    """
    c_min, c_max = channel_stats[MIN], channel_stats[MAX]
    if (c_max - c_min) > EPSILON:
        channel_data[...] = channel_data * (c_max - c_min) + c_min
