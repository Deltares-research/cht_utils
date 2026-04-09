"""Spatial interpolation utilities for regular and unstructured grids."""

from typing import Union

import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata


def interp2(
    x0: np.ndarray,
    y0: np.ndarray,
    z0: np.ndarray,
    x1: np.ndarray,
    y1: np.ndarray,
    method: str = "linear",
) -> np.ndarray:
    """Interpolate from a regular grid to arbitrary points.

    Parameters
    ----------
    x0 : np.ndarray
        1-D array of x-coordinates of the source grid.
    y0 : np.ndarray
        1-D array of y-coordinates of the source grid.
    z0 : np.ndarray
        2-D array of values on the source grid, shape ``(len(y0), len(x0))``.
    x1 : np.ndarray
        Target x-coordinates (1-D or 2-D).
    y1 : np.ndarray
        Target y-coordinates (1-D or 2-D).
    method : str
        Interpolation method (``"linear"``, ``"nearest"``).

    Returns
    -------
    np.ndarray
        Interpolated values at the target locations.
    """
    f = RegularGridInterpolator(
        (y0, x0), z0, bounds_error=False, fill_value=np.nan, method=method
    )
    if x1.ndim > 1:
        sz = x1.shape
        z1 = f((y1.ravel(), x1.ravel())).reshape(sz)
    else:
        z1 = f((y1, x1))
    return z1


def interp2_bilinear(
    xp: np.ndarray,
    yp: np.ndarray,
    zp: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Bilinear interpolation on a regular grid.

    Parameters
    ----------
    xp : np.ndarray
        1-D array of source x-coordinates (monotonically increasing).
    yp : np.ndarray
        1-D array of source y-coordinates (monotonically increasing).
    zp : np.ndarray
        2-D source values, shape ``(len(yp), len(xp))``.
    x : np.ndarray
        Target x-coordinates.
    y : np.ndarray
        Target y-coordinates.

    Returns
    -------
    np.ndarray
        Interpolated values at the target locations.
    """
    dx = xp[1] - xp[0]
    dy = yp[1] - yp[0]

    col = (x - xp[0]) / dx
    row = (y - yp[0]) / dy

    c0 = np.floor(col).astype(int)
    r0 = np.floor(row).astype(int)
    c1 = c0 + 1
    r1 = r0 + 1

    # Clip to valid range
    c0 = np.clip(c0, 0, len(xp) - 1)
    c1 = np.clip(c1, 0, len(xp) - 1)
    r0 = np.clip(r0, 0, len(yp) - 1)
    r1 = np.clip(r1, 0, len(yp) - 1)

    # Fractional parts
    dc = col - np.floor(col)
    dr = row - np.floor(row)

    z = (
        zp[r0, c0] * (1 - dc) * (1 - dr)
        + zp[r0, c1] * dc * (1 - dr)
        + zp[r1, c0] * (1 - dc) * dr
        + zp[r1, c1] * dc * dr
    )
    return z


def interp3(
    x0: np.ndarray,
    y0: np.ndarray,
    z0: np.ndarray,
    x1: np.ndarray,
    y1: np.ndarray,
    method: str = "linear",
) -> np.ndarray:
    """Interpolate from unstructured or curvilinear points.

    Uses :func:`scipy.interpolate.griddata` for scattered-data interpolation.

    Parameters
    ----------
    x0 : np.ndarray
        Source x-coordinates (1-D or 2-D).
    y0 : np.ndarray
        Source y-coordinates (1-D or 2-D).
    z0 : np.ndarray
        Source values.
    x1 : np.ndarray
        Target x-coordinates (1-D or 2-D).
    y1 : np.ndarray
        Target y-coordinates (1-D or 2-D).
    method : str
        Interpolation method (``"linear"``, ``"nearest"``, ``"cubic"``).

    Returns
    -------
    np.ndarray
        Interpolated values at the target locations.
    """
    if x1.ndim > 1:
        sz = x1.shape
        z1 = griddata(
            (x0.ravel(), y0.ravel()),
            z0.ravel(),
            (x1.ravel(), y1.ravel()),
            method=method,
        ).reshape(sz)
    else:
        z1 = griddata(
            (x0.ravel(), y0.ravel()),
            z0.ravel(),
            (x1, y1),
            method=method,
        )
    return z1
