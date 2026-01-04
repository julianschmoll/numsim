"""Visualization utilities for fluid data."""
from types import MappingProxyType

from matplotlib import pyplot as plt

from constants import *  # noqa: F403, WPS347


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
    plt.quiver(u_component, v_component, magnitude, cmap=COLORMAP)


PLOT_FUNCTIONS = MappingProxyType(
    {
        "implot": _implot,
        "quiver": _quiver,
    }
)


def visualize(input_data, title="Visualization", plt_fn="implot", num_rows=1):
    """
    Visualizes the input and label channels.

    Args:
        input_data: List of tuples containing data and labels to be plotted.
        title: Title of the figure.
        plt_fn: Name of the plotting function to use (e.g. 'implot' or 'quiver').
        num_rows: Number of rows in the figure.
    """
    num_cols = (len(input_data) + num_rows - 1) // num_rows

    fig = plt.figure(figsize=(FIG_WIDTH * num_cols, FIG_HEIGHT * num_rows))

    plt.suptitle(
        title, fontsize=TITLE_FONT_SIZE, fontweight="bold", y=TITLE_Y_POSITION,
    )

    plot_function = PLOT_FUNCTIONS.get(plt_fn)
    if not plot_function:
        raise ValueError(f"Unknown plot function: {plt_fn}")

    for idx, plt_data in enumerate(input_data):
        plot_function(fig, *plt_data, (num_rows, num_cols, idx + 1))

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def loss_plot(save_path, loss_data, title="Loss Plot"):
    """Saves Semilogy plot of MSE loss.

    Args:
        save_path: Path to save the plot.
        loss_data: List of tuples containing loss values and their labels.
        title: Title of the plot.
    """
    plt.figure()
    for loss, label in loss_data:
        plt.semilogy(loss, label=label)
    plt.xlabel("Epoch")
    plt.ylabel("MSE (log scale)")
    plt.legend()
    plt.title(title)
    plt.savefig(save_path)
