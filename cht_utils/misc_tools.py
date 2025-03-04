# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 16:28:15 2021

@author: ormondt
"""

from scipy.interpolate import RegularGridInterpolator, griddata
import numpy as np
import json
import yaml


def interp2(x0, y0, z0, x1, y1, method="linear"):

    # meanx = np.mean(x0)
    # meany = np.mean(y0)
    # x0 -= meanx
    # y0 -= meany
    # x1 -= meanx
    # y1 -= meany

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


def interp2_bilinear(xp, yp, zp, x, y):
    """
    Bilinearly interpolate over regular 2D grid.

    `xp` and `yp` are 1D arrays defining grid coordinates of sizes :math:`N_x`
    and :math:`N_y` respectively, and `zp` is the 2D array, shape
    :math:`(N_x, N_y)`, containing the gridded data points which are being
    interpolated from. Note that the coordinate grid should be regular, i.e.
    uniform grid spacing. `x` and `y` are either scalars or 1D arrays giving
    the coordinates of the points at which to interpolate. If these are outside
    the boundaries of the coordinate grid, the resulting interpolated values
    are evaluated at the boundary.

    Parameters
    ----------
    x : 1D array or scalar
        x-coordinates of interpolating point(s).
    y : 1D array or scalar
        y-coordinates of interpolating point(s).
    xp : 1D array, shape M
        x-coordinates of data points zp. Note that this should be a *regular*
        grid, i.e. uniform spacing.
    yp : 1D array, shape N
        y-coordinates of data points zp. Note that this should be a *regular*
        grid, i.e. uniform spacing.
    zp : 2D array, shape (M, N)
        Data points on grid from which to interpolate.

    Returns
    -------
    z : 1D array or scalar
        Interpolated values at given point(s).

    """
    # Transpose zp
    zp = zp.T

    sz = None
    if x.ndim > 1:
        # Flatten x and y
        sz = x.shape
        x = x.reshape(sz[0] * sz[1])
        y = y.reshape(sz[0] * sz[1])

    # if scalar, turn into array
    scalar = False
    if not isinstance(x, (list, np.ndarray)):
        scalar = True
        x = np.array([x])
        y = np.array([y])

    # grid spacings and sizes
    hx = xp[1] - xp[0]
    hy = yp[1] - yp[0]
    Nx = xp.size
    Ny = yp.size

    # snap beyond-boundary points to boundary
    x[x < xp[0]] = xp[0]
    y[y < yp[0]] = yp[0]
    x[x > xp[-1]] = xp[-1]
    y[y > yp[-1]] = yp[-1]

    # find indices of surrounding points
    i1 = np.floor((x - xp[0]) / hx).astype(int)
    i1[i1 == Nx - 1] = Nx - 2
    j1 = np.floor((y - yp[0]) / hy).astype(int)
    j1[j1 == Ny - 1] = Ny - 2
    i2 = i1 + 1
    j2 = j1 + 1

    # get coords and func at surrounding points
    x1 = xp[i1]
    x2 = xp[i2]
    y1 = yp[j1]
    y2 = yp[j2]
    z11 = zp[i1, j1]
    z21 = zp[i2, j1]
    z12 = zp[i1, j2]
    z22 = zp[i2, j2]

    # interpolate
    t11 = z11 * (x2 - x) * (y2 - y)
    t21 = z21 * (x - x1) * (y2 - y)
    t12 = z12 * (x2 - x) * (y - y1)
    t22 = z22 * (x - x1) * (y - y1)
    z = (t11 + t21 + t12 + t22) / (hx * hy)
    if scalar:
        z = z[0]

    if sz is not None:
        z = z.reshape(sz)

    return z


def interp3(x0, y0, z0, x1, y1, method="linear"):
    """interpolation from unstructured or curvilinear source data

    Args:
        x0 (_type_): x coords source grid
        y0 (_type_): y coords source grid
        z0 (_type_): z coords source grid
        x1 (_type_): x coords target grid
        y1 (_type_): y coords target grid
        method (str, optional): interpolator method. Defaults to "linear".

    Returns:
        _type_: z coords target grid
    """
    sz0 = x0.shape
    x0 = x0.reshape(sz0[0] * sz0[1])
    y0 = y0.reshape(sz0[0] * sz0[1])
    z0 = z0.reshape(sz0[0] * sz0[1])

    # reshape x1 and y1
    if np.atleast_1d(x1).ndim > 1:
        sz = x1.shape
        x1 = x1.reshape(sz[0] * sz[1])
        y1 = y1.reshape(sz[0] * sz[1])
        # interpolate
        z1 = griddata(
            np.array([x0, y0]).T,
            z0.T,
            np.array([x1, y1]).T,
            fill_value=np.nan,
            method=method,
        )

        return z1.reshape(sz)
    else:
        z1 = griddata(
            np.array([x0, y0]).T,
            z0.T,
            np.array([x1, y1]).T,
            fill_value=np.nan,
            method=method,
        )
        return z1


def findreplace(file_name, str1, str2):

    # read input file
    fin = open(file_name, "rt")
    # read file contents to string
    data = fin.read()
    # replace all occurrences of the required string
    data = data.replace(str1, str2)
    # close the input file
    fin.close()
    # open the input file in write mode
    fin = open(file_name, "wt")
    # overrite the input file with the resulting data
    fin.write(data)
    # close the file
    fin.close()


def read_json_js(file_name):
    # Read json javascript file (skipping the first line), and return json object
    fid = open(file_name, "r")
    lines = fid.readlines()
    fid.close()
    lines = lines[1:]
    jsn_string = ""
    for line in lines:
        # Strips the newline character
        jsn_string += line.strip()
    jsn = json.loads(jsn_string)
    return jsn


def write_json_js(file_name, jsn, first_line):
    # Writes json javascript file
    if type(jsn) == list:
        f = open(file_name, "w")
        f.write(first_line + "\n")
        f.write("[\n")
        for ix, x in enumerate(jsn):
            json_string = json.dumps(x)
            if ix < len(jsn) - 1:
                f.write(json_string + ",")
            else:
                f.write(json_string)
            f.write("\n")
        f.write("]\n")
        f.close()
    else:
        f = open(file_name, "w")
        f.write(first_line + "\n")
        json_string = json.dumps(jsn)
        f.write(json_string + "\n")
        f.close()


def write_csv_js(file_name, csv_string, first_line):
    # If folder does not exist: make folder
    import os

    os.makedirs(os.path.dirname(file_name), exist_ok=True)

    # Writes json javascript file
    csv_string = csv_string.replace(chr(13), "")
    f = open(file_name, "w")
    f.write(first_line + "\n")
    f.write(csv_string)
    f.write("`;")
    f.close()


def rgb2hex(rgb):
    return "%02x%02x%02x" % rgb


def dict2yaml(file_name, dct, sort_keys=False):
    yaml_string = yaml.dump(dct, sort_keys=sort_keys)
    file = open(file_name, "w")
    file.write(yaml_string)
    file.close()


def yaml2dict(file_name):
    file = open(file_name, "r")
    dct = yaml.load(file, Loader=yaml.FullLoader)
    return dct
