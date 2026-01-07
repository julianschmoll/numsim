"""FluidCNN model for fluid flow predictions."""
import json
from pathlib import Path

import numpy as np
import torch

# We bake this in here instead of importing constants
ACTIVATION = "activation"
IN_CHANNELS = "in_channels"
HIDDEN_CHANNELS = "hidden_channels"
OUT_CHANNELS = "out_channels"
KERNEL_SIZE = "kernel_size"
PADDING_MODE = "padding_mode"
USE_BIAS = "use_bias"
NUM_HIDDEN_LAYERS = "num_hidden_layers"
DROPOUT_RATE = "dropout_rate"
OUTPUT_ACTIVATION = "output_activation"


class FluidCNN(torch.nn.Module):
    """Configurable CNN for fluid flow predictions."""
    def __init__(self, config=None):
        """Initialize the FluidCNN model.

        Args:
            config (dict, optional): Configuration parameters for the model.
        """
        super().__init__()

        if config is None:
            config = {}

        activation_layer = getattr(torch.nn, config.get(ACTIVATION, "ReLU"))

        layers = [
            torch.nn.Conv2d(
                in_channels=config.get(IN_CHANNELS, 1),
                out_channels=config.get(HIDDEN_CHANNELS, 8),
                kernel_size=config.get(KERNEL_SIZE, 7),
                stride=1,
                padding="same",
                padding_mode=config.get(PADDING_MODE, "zeros"),
                bias=config.get(USE_BIAS, True),
            ),
            activation_layer(),
        ]

        for _ in range(config.get(NUM_HIDDEN_LAYERS, 5)):
            layers.append(
                torch.nn.Conv2d(
                    in_channels=config.get(HIDDEN_CHANNELS, 8),
                    out_channels=config.get(HIDDEN_CHANNELS, 8),
                    kernel_size=config.get(KERNEL_SIZE, 7),
                    stride=1,
                    padding="same",
                    padding_mode=config.get(PADDING_MODE, "zeros"),
                    bias=config.get(USE_BIAS, True),
                )
            )
            layers.append(activation_layer())
            if config.get(DROPOUT_RATE):
                layers.append(torch.nn.Dropout2d(p=config[DROPOUT_RATE]))

        layers.append(
            torch.nn.Conv2d(
                in_channels=config.get(HIDDEN_CHANNELS, 8),
                out_channels=config.get(OUT_CHANNELS, 2),
                kernel_size=config.get(KERNEL_SIZE, 7),
                stride=1,
                padding="same",
                padding_mode=config.get(PADDING_MODE, "zeros"),
                bias=config.get(USE_BIAS, True),
            )
        )

        if config.get(OUTPUT_ACTIVATION):
            layers.append(config[OUTPUT_ACTIVATION]())

        self.net = torch.nn.Sequential(*layers)

    def forward(self, input_tensor):
        """Forward pass of the FluidCNN model.

        Args:
            input_tensor (torch.Tensor): Input tensor to the model.

        Returns:
            torch.Tensor: Output tensor from the model.
        """
        return self.net(input_tensor)

    def predict(self, flow_speed, hx, hy):
        """Predicts lid driven cavity scenario.

        Args:
            flow_speed: Speed of flow at lid.
            hx: Width of the domain.
            hy: Height of the domain.

        Returns:
            Predicted tensor.
        """
        input_channel = np.zeros((1, hx, hy), dtype=np.float32)
        input_channel[0, -1, 1:-1] = flow_speed
        input_tensor = torch.from_numpy(input_channel).unsqueeze(0)
        device = next(self.parameters()).device
        self.eval()
        with torch.no_grad():
            output = self.forward(input_tensor.to(device))
        return input_tensor.detach(), output.detach()


class Sine(torch.nn.Module):
    """Sine activation function."""
    def forward(self, input_tensor):
        """Forward pass of the Sine activation function.

        Args:
            input_tensor (torch.Tensor): Input tensor to the activation function.

        Returns:
            torch.Tensor: Output tensor after applying sine activation.
        """
        return torch.sin(input_tensor)


# this has to be named like this so the submission system can read the model
def init_my_model():
    """Initialize the FluidCNN model from a configuration file.

    Returns:
        FluidCNN: An instance of the FluidCNN model.
    """
    config_path = Path(__file__).resolve().parent / "config.json"

    if not config_path.is_file():
        raise RuntimeError(f"Config file not found at {config_path}")

    with open(config_path, "r", encoding="utf-8") as config_file:
        config = json.load(config_file)

    return FluidCNN(config)
