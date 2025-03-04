# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 16:28:15 2021

@author: ormondt
"""

from scipy.interpolate import RegularGridInterpolator
import numpy as np


def interp2(x0, y0, z0, x1, y1, method="linear"):
    """Interpolate from regular grid (numpy arrays x0, y0, z0) to points x1, y1 (can be vectors or 2d matrices)"""

    f = RegularGridInterpolator(
        (y0, x0), z0, bounds_error=False, fill_value=np.nan, method=method
    )
    # reshape x1 and y1
    if x1.ndim > 1:
        sz = x1.shape
        x1 = x1.reshape(sz[0] * sz[1])
        y1 = y1.reshape(sz[0] * sz[1])
        # interpolate
        z1 = f((y1, x1)).reshape(sz)
    else:
        z1 = f((y1, x1))

    return z1
