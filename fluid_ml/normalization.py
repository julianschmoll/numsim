import torch
import yaml

import constants


DEFAULT_OUT_RANGE = (torch.tensor([.0]), torch.tensor([1.0]))


def denormalize(inputs, labels, stats_path):
    with open(stats_path, "r") as min_max_yaml:
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
    inputs_mins = torch.tensor(
        [min_max_stats[constants.INPUTS_KEY][constants.U_CHANNEL_KEY][constants.MIN]]
    )
    inputs_maxs = torch.tensor(
        [min_max_stats[constants.INPUTS_KEY][constants.U_CHANNEL_KEY][constants.MAX]]
    )
    return rescale(inputs, inputs_mins, inputs_maxs)


def rescale_labels(labels, stats):
    labels_min = torch.tensor([
        stats[constants.LABELS_KEY][constants.U_CHANNEL_KEY][constants.MIN],
        stats[constants.LABELS_KEY][constants.V_CHANNEL_KEY][constants.MIN]
    ])
    labels_max = torch.tensor([
        stats[constants.LABELS_KEY][constants.U_CHANNEL_KEY][constants.MAX],
        stats[constants.LABELS_KEY][constants.V_CHANNEL_KEY][constants.MAX]
    ])
    return rescale(labels, *DEFAULT_OUT_RANGE, out_range=(labels_min, labels_max))


def normalize_channel(channel_data):
    """Normalize a single data channel in-place and return its stats."""
    c_min = float(channel_data.min())
    c_max = float(channel_data.max())
    if (c_max - c_min) > constants.EPSILON:
        channel_data[...] = (channel_data - c_min) / (c_max - c_min)
    return {constants.MIN: c_min, constants.MAX: c_max}


def denormalize_channel(channel_data, channel_stats):
    """Denormalize a single data channel in-place using provided stats."""
    c_min, c_max = channel_stats[constants.MIN], channel_stats[constants.MAX]
    if (c_max - c_min) > constants.EPSILON:
        channel_data[...] = channel_data * (c_max - c_min) + c_min
