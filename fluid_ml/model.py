from torch import nn


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
        activation_layer = config.get("activation", nn.ReLU)
        output_activation = config.get("output_activation")
        use_bias = config.get("use_bias", True)
        padding_mode = config.get("padding_mode", "zeros")

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
