"""Colormap utilities for reading, registering, and rendering colormaps."""

import os
from typing import List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_color_maps(path_name: str) -> List[str]:
    """Read all colormaps from a folder and register them with matplotlib.

    Each ``.txt`` file in *path_name* is expected to contain whitespace-separated
    RGB values (one row per colour, values in 0-1 range).

    Parameters
    ----------
    path_name : str
        Directory containing colormap text files.

    Returns
    -------
    List[str]
        Full list of registered matplotlib colormap names.
    """
    for file in os.listdir(path_name):
        if file.endswith(".txt"):
            name = os.path.splitext(file)[0]
            rgb = read_colormap(os.path.join(path_name, file))
            cmap = mpl.colors.ListedColormap(rgb, name=name)
            mpl.colormaps.register(cmap=cmap)
    return plt.colormaps()


def cm2png(
    cmap: mpl.colors.Colormap,
    file_name: str = "colorbar.png",
    orientation: str = "horizontal",
    vmin: float = 0.0,
    vmax: float = 1.0,
    legend_title: str = "",
    legend_label: str = "",
    units: str = "",
    unit_string: str = "",
    decimals: int = -1,
) -> None:
    """Render a colormap to a PNG image.

    Parameters
    ----------
    cmap : matplotlib.colors.Colormap
        Colormap to render.
    file_name : str
        Output PNG path.
    orientation : str
        ``"horizontal"`` or ``"vertical"``.
    vmin, vmax : float
        Value range for the colour scale.
    legend_label : str
        Label shown alongside the colorbar.
    """
    if orientation == "horizontal":
        fig = plt.figure(figsize=(2.5, 1))
        ax = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    else:
        fig = plt.figure(figsize=(1, 2.5))
        ax = fig.add_axes([0.80, 0.05, 0.15, 0.90])

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = mpl.colorbar.ColorbarBase(
        ax, cmap=cmap, norm=norm, orientation=orientation, label=legend_label
    )
    cb.ax.tick_params(labelsize=6)

    fig.savefig(file_name, dpi=150, bbox_inches="tight")
    plt.close(fig)


def read_colormap(file_name: str) -> np.ndarray:
    """Read a colormap from a whitespace-separated RGB text file.

    Parameters
    ----------
    file_name : str
        Path to the colormap file.

    Returns
    -------
    np.ndarray
        Array of shape ``(N, 3)`` with RGB values.
    """
    df = pd.read_csv(
        file_name,
        index_col=False,
        header=None,
        sep=r"\s+",
        names=["r", "g", "b"],
    )
    return df.to_numpy()


def rgb2hex(rgb: Tuple[int, int, int]) -> str:
    """Convert an RGB tuple to a hex colour string.

    Parameters
    ----------
    rgb : Tuple[int, int, int]
        Red, green, blue values (0-255).

    Returns
    -------
    str
        Six-character hex string (no leading ``#``).
    """
    return "%02x%02x%02x" % rgb
