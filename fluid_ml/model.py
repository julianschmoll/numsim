import json
from pathlib import Path

from torch import nn

# This could be in a seperate constants.py file, but we keep it here for simplicity
ACTIVATION_KEY = "activation"
IN_CHANNELS_KEY = "in_channels"
HIDDEN_CHANNELS_KEY = "hidden_channels"
OUT_CHANNELS_KEY = "out_channels"
KERNEL_SIZE_KEY = "kernel_size"
PADDING_MODE_KEY = "padding_mode"
USE_BIAS_KEY = "use_bias"
NUM_HIDDEN_LAYERS_KEY = "num_hidden_layers"
DROPOUT_RATE_KEY = "dropout_rate"
OUTPUT_ACTIVATION_KEY = "output_activation"


class FluidCNN(nn.Module):
    """Configurable CNN for fluid flow predictions."""
    def __init__(self, config=None):
        """Initialize the FluidCNN model.

        Args:
            config (dict, optional): Configuration parameters for the model.
        """
        super().__init__()

        if config is None:
            config = {}

        activation_layer = getattr(nn, config.get(ACTIVATION_KEY, "ReLU"))

        layers = [
            nn.Conv2d(
                in_channels=config.get(IN_CHANNELS_KEY, 1),
                out_channels=config.get(HIDDEN_CHANNELS_KEY, 8),
                kernel_size=config.get(KERNEL_SIZE_KEY, 7),
                stride=1,
                padding="same",
                padding_mode=config.get(PADDING_MODE_KEY, "zeros"),
                bias=config.get(USE_BIAS_KEY, True),
            ),
            activation_layer(),
        ]

        for _ in range(config.get(NUM_HIDDEN_LAYERS_KEY, 5)):
            layers.append(
                nn.Conv2d(
                    in_channels=config.get(HIDDEN_CHANNELS_KEY, 8),
                    out_channels=config.get(HIDDEN_CHANNELS_KEY, 8),
                    kernel_size=config.get(KERNEL_SIZE_KEY, 7),
                    stride=1,
                    padding="same",
                    padding_mode=config.get(PADDING_MODE_KEY, "zeros"),
                    bias=config.get(USE_BIAS_KEY, True),
                )
            )
            layers.append(activation_layer())
            if config.get(DROPOUT_RATE_KEY):
                layers.append(nn.Dropout2d(p=config[DROPOUT_RATE_KEY]))

        layers.append(
            nn.Conv2d(
                in_channels=config.get(HIDDEN_CHANNELS_KEY, 8),
                out_channels=config.get(OUT_CHANNELS_KEY, 2),
                kernel_size=config.get(KERNEL_SIZE_KEY, 7),
                stride=1,
                padding="same",
                padding_mode=config.get(PADDING_MODE_KEY, "zeros"),
                bias=config.get(USE_BIAS_KEY, True),
            )
        )

        if config.get(OUTPUT_ACTIVATION_KEY):
            layers.append(config[OUTPUT_ACTIVATION_KEY]())

        self.net = nn.Sequential(*layers)

    def forward(self, input_tensor):
        """Forward pass of the FluidCNN model.

        Args:
            input_tensor (torch.Tensor): Input tensor to the model.

        Returns:
            torch.Tensor: Output tensor from the model.
        """
        return self.net(input_tensor)


# this has to be named like this so the submission system can read the model
def init_my_model():
    """Initialize the FluidCNN model from a configuration file.

    Returns:
        FluidCNN: An instance of the FluidCNN model.
    """
    config_path = Path(__file__).resolve().parent / "config.json"

    if not config_path.is_file():
        raise RuntimeError(f"Config file not found at {config_path}")

    with open(config_path, "r") as config_file:
        config = json.load(config_file)

    activation_str = config.get(ACTIVATION_KEY)
    if activation_str and isinstance(activation_str, str):
        config[ACTIVATION_KEY] = getattr(nn, activation_str)

    return FluidCNN(config)
