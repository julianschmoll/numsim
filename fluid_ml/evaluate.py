import matplotlib.pyplot as plt


def visualize(inputs, labels, title="Visualization"):
    """Visualizes the input and label channels

    Args:
        inputs: input tensor
        labels: label tensor
        title: title of the figure
    """
    fig = plt.figure(figsize=(15, 5))
    plt.suptitle(title, fontsize=20, fontweight="bold", y=0.95)

    _plot(fig, inputs[0, 0], "Input $u$", (1, 3, 1))
    _plot(fig, labels[0, 0], "Output $u$", (1, 3, 2))
    _plot(fig, labels[0, 1], "Output $v$", (1, 3, 3))

    plt.tight_layout()
    plt.show()


def _plot(fig, data, label, position):
    """Adds a subplot to the figure with given data.

    Args:
        fig: figure to add the subplot to
        data: data to be plotted
        label: label of the subplot (legend)
        position: position of the subplot
    """
    ax = fig.add_subplot(*position)
    im = ax.imshow(data, origin='lower')
    ax.set_title(label)
    plt.colorbar(im, ax=ax)
