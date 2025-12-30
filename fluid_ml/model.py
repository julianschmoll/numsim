from torch import nn
from pathlib import Path
import json


class FluidCNN(nn.Module):
    def __init__(self, config=None):
        super().__init__()

        if config is None:
            config = {}

        num_hidden_layers = config.get("num_hidden_layers", 5)
        kernel_size = config.get("kernel_size", 7)
        in_channels = config.get("in_channels", 1)
        hidden_channels = config.get("hidden_channels", 16)
        out_channels = config.get("out_channels", 2)
        activation_str = config.get("activation", "ReLU")
        activation_layer = getattr(nn, activation_str)
        output_activation = config.get("output_activation")
        use_bias = config.get("use_bias", True)
        padding_mode = config.get("padding_mode", "zeros")
        dropout_rate = config.get("dropout_rate", 0.0)

        layers = [
            nn.Conv2d(
                in_channels=in_channels,
                out_channels=hidden_channels,
                kernel_size=kernel_size,
                stride=1,
                padding="same",
                padding_mode=padding_mode,
                bias=use_bias,
            ),
            activation_layer()
        ]

        for _ in range(num_hidden_layers):
            layers.append(nn.Conv2d(
                in_channels=hidden_channels,
                out_channels=hidden_channels,
                kernel_size=kernel_size,
                stride=1,
                padding="same",
                padding_mode=padding_mode,
                bias=use_bias,
            ))
            layers.append(activation_layer())
            if dropout_rate > 0:
                layers.append(nn.Dropout2d(p=dropout_rate))

        layers.append(nn.Conv2d(
            in_channels=hidden_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            stride=1,
            padding="same",
            padding_mode=padding_mode,
            bias=use_bias,
        ))

        if output_activation:
            layers.append(output_activation())

        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


# this has to be named like this so the submission system can read the model
def init_my_model():
    config_path = Path(__file__).resolve().parent / "config.json"

    if not config_path.is_file():
        print(f"Warning: Config not found at {config_path}. Using default architecture.")
        return FluidCNN()

    with open(config_path, "r") as f:
        config = json.load(f)

    activation_str = config.get("activation")
    if activation_str and isinstance(activation_str, str):
        config["activation"] = getattr(nn, activation_str)

    return FluidCNN(config)
