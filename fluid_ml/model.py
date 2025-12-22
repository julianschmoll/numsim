from torch import nn


class FluidCNN(nn.Module):
    def __init__(self, config):
        super().__init__()
        layers = [
            nn.Conv2d(
                in_channels=1,
                out_channels=16,
                kernel_size=7,
                stride=1,
                padding="same"
            ),
            nn.ReLU()
        ]

        for _ in range(5):
            layers.append(nn.Conv2d(
                in_channels=16,
                out_channels=16,
                kernel_size=7,
                stride=1,
                padding="same"
            ))
            layers.append(nn.ReLU())

        layers.append(nn.Conv2d(
            in_channels=16,
            out_channels=2,
            kernel_size=7,
            stride=1,
            padding="same"
        ))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)
