from matplotlib import pyplot as plt

import constants
import normalization


def visualize(inputs, labels, title="Visualization", quiver=False, stats_path=None):
    """
    Visualizes the input and label channels.

    Args:
        inputs: Input tensor.
        labels: Label tensor.
        title: Title of the figure.
        quiver: Whether to display as quiver plot.
        stats_path: Path to the directory containing 'min_max.yaml'.
    """
    fig_width = constants.FIG_WIDTH * 2 if quiver else constants.FIG_WIDTH * 3
    fig = plt.figure(figsize=(fig_width, constants.FIG_HEIGHT))
    plt.suptitle(
        title,
        fontsize=constants.TITLE_FONT_SIZE,
        fontweight="bold",
        y=constants.TITLE_Y_POSITION,
    )

    if stats_path:
        inputs, labels = normalization.denormalize(inputs, labels, stats_path)

    if quiver:
        _quiver(fig, inputs[0], "Flow input", (1, 2, 1))
        _quiver(fig, labels[0], "Flow output", (1, 2, 2))
    else:
        _implot(fig, inputs[0, 0], "Input $u$", (1, 3, 1))
        _implot(fig, labels[0, 0], "Output $u$", (1, 3, 2))
        _implot(fig, labels[0, 1], "Output $v$", (1, 3, 3))

    plt.tight_layout()
    plt.show()


def _implot(fig, image_data, label, position):
    """
    Adds a subplot to the figure with given data.

    Args:
        fig: Figure to add the subplot to.
        image_data: Data to be plotted.
        label: Label of the subplot (legend).
        position: Position of the subplot.
    """
    axis = fig.add_subplot(*position)
    image = axis.imshow(image_data, origin="lower")
    axis.set_title(label)
    plt.colorbar(image, ax=axis)


def _quiver(fig, vector_data, label, position):
    """
    Adds a quiver subplot to the figure with given data.

    Args:
        fig: Figure to add the subplot to.
        vector_data: Data to be plotted.
        label: Label of the subplot (legend).
        position: Position of the subplot.
    """
    axis = fig.add_subplot(*position)
    axis.set_title(label)

    u_component = vector_data[0]
    v_component = 0
    if vector_data.shape[0] > 1:
        v_component = vector_data[1]

    magnitude = (u_component**2 + v_component**2) ** 0.5
    plt.quiver(u_component, v_component, magnitude, cmap=constants.COLORMAP)
