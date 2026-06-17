"""Utility functions for tile coordinate conversions, elevation/PNG encoding, and interpolation."""

from __future__ import annotations

import glob
import math
import os

import numpy as np
from numpy.typing import NDArray
from PIL import Image
from scipy.interpolate import RegularGridInterpolator


def get_zoom_level_for_resolution(dx: float) -> int:
    """Determine the tile zoom level that matches a given spatial resolution.

    Parameters
    ----------
    dx : float
        Desired pixel size in metres.

    Returns
    -------
    int
        Zoom level (0-23) whose native pixel size is just below *dx*.
    """
    dxy = 156543.03 / 2 ** np.arange(24)
    izoom = np.where(dxy < dx)[0]
    if len(izoom) == 0:
        izoom = 23
    else:
        izoom = int(izoom[0])
    return izoom


def get_zoom_level(npixels: int, lat_range: list[float], max_zoom: int) -> int:
    """Determine the zoom level needed to cover a latitude range with a given pixel count.

    Parameters
    ----------
    npixels : int
        Number of pixels available in the latitude direction.
    lat_range : list[float]
        Two-element list ``[lat_min, lat_max]`` in degrees.
    max_zoom : int
        Maximum allowed zoom level.

    Returns
    -------
    int
        Appropriate zoom level.
    """
    dxr = (lat_range[1] - lat_range[0]) * 111111 / npixels
    dxy = 156543.03 / 2 ** np.arange(max_zoom + 1)
    izoom = np.where(dxy < dxr)[0]
    if len(izoom) == 0:
        izoom = max_zoom
    else:
        izoom = izoom[0]
    return izoom


def webmercator_to_lat_lon(easting: float, northing: float) -> tuple[float, float]:
    """Convert Web Mercator coordinates to latitude and longitude.

    Parameters
    ----------
    easting : float
        Web Mercator easting in metres.
    northing : float
        Web Mercator northing in metres.

    Returns
    -------
    tuple[float, float]
        ``(latitude, longitude)`` in degrees.
    """
    lon = (easting / 20037508.34) * 180
    lat = (180 / math.pi) * (
        2 * math.atan(math.exp(northing / 20037508.34 * math.pi)) - (math.pi / 2)
    )
    return lat, lon


def lat_lon_to_webmercator(lat: float, lon: float) -> tuple[float, float]:
    """Convert latitude and longitude to Web Mercator coordinates.

    Parameters
    ----------
    lat : float
        Latitude in degrees.
    lon : float
        Longitude in degrees.

    Returns
    -------
    tuple[float, float]
        ``(x, y)`` in Web Mercator metres.
    """
    x = lon * 20037508.34 / 180
    y = (math.log(math.tan((90 + lat) * math.pi / 360)) / math.pi) * 20037508.34
    return x, y


def lat_lon_to_tile_indices(lat: float, lon: float, zoom: int) -> tuple[int, int]:
    """Convert latitude/longitude to slippy-map tile indices.

    Parameters
    ----------
    lat : float
        Latitude in degrees.
    lon : float
        Longitude in degrees.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[int, int]
        ``(tile_x, tile_y)`` column and row indices.
    """
    tile_x = int((lon + 180) / 360 * (2**zoom))
    tile_y = int(
        (
            1
            - (
                math.log(math.tan(math.radians(lat)) + 1 / math.cos(math.radians(lat)))
                / math.pi
            )
        )
        / 2
        * (2**zoom)
    )
    return tile_x, tile_y


def xy2num(easting: float, northing: float, zoom: int) -> tuple[int, int]:
    """Convert Web Mercator coordinates to tile indices.

    Parameters
    ----------
    easting : float
        Web Mercator easting in metres.
    northing : float
        Web Mercator northing in metres.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[int, int]
        ``(tile_x, tile_y)`` column and row indices.
    """
    lat, lon = webmercator_to_lat_lon(easting, northing)
    ix, it = lat_lon_to_tile_indices(lat, lon, zoom)
    return ix, it


def deg2num(lat_deg: float, lon_deg: float, zoom: int) -> tuple[int, int]:
    """Return the column and row index of a slippy tile for a given lat/lon.

    Parameters
    ----------
    lat_deg : float
        Latitude in degrees.
    lon_deg : float
        Longitude in degrees.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[int, int]
        ``(xtile, ytile)`` column and row indices.
    """
    lat_rad = math.radians(lat_deg)
    n = 2**zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
    return (xtile, ytile)


def num2deg(xtile: int, ytile: int, zoom: int) -> tuple[float, float]:
    """Return the upper-left latitude and longitude of a slippy tile.

    Parameters
    ----------
    xtile : int
        Tile column index.
    ytile : int
        Tile row index.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[float, float]
        ``(latitude, longitude)`` of the upper-left corner in degrees.
    """
    n = 2**zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return (lat_deg, lon_deg)


def num2xy(xtile: int, ytile: int, zoom: int) -> tuple[float, float]:
    """Return the upper-left Web Mercator x/y of a slippy tile.

    Parameters
    ----------
    xtile : int
        Tile column index.
    ytile : int
        Tile row index.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[float, float]
        ``(x, y)`` in Web Mercator metres.
    """
    n = 2**zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    x, y = lat_lon_to_webmercator(lat_deg, lon_deg)
    return x, y


def num2deg_ll(xtile: int, ytile: int, zoom: int) -> tuple[float, float]:
    """Return the lower-left latitude and longitude of a slippy tile (old format).

    Parameters
    ----------
    xtile : int
        Tile column index.
    ytile : int
        Tile row index.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[float, float]
        ``(latitude, longitude)`` of the lower-left corner in degrees.
    """
    n = 2**zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(-lat_rad)
    return (lat_deg, lon_deg)


def num2deg_ur(xtile: int, ytile: int, zoom: int) -> tuple[float, float]:
    """Return the upper-right latitude and longitude of a slippy tile (old format).

    Parameters
    ----------
    xtile : int
        Tile column index.
    ytile : int
        Tile row index.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[float, float]
        ``(latitude, longitude)`` of the upper-right corner in degrees.
    """
    n = 2**zoom
    lon_deg = (xtile + 1) / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * (ytile + 1) / n)))
    lat_deg = math.degrees(-lat_rad)
    return (lat_deg, lon_deg)


def get_lower_left_corner(tile_x: int, tile_y: int, zoom: int) -> tuple[float, float]:
    """Return the lower-left Web Mercator coordinates of a tile.

    Parameters
    ----------
    tile_x : int
        Tile column index.
    tile_y : int
        Tile row index.
    zoom : int
        Zoom level.

    Returns
    -------
    tuple[float, float]
        ``(ll_x, ll_y)`` lower-left corner in Web Mercator metres.
    """
    total_size = 20037508.34 * 2
    tile_size = total_size / (2**zoom)
    ll_x = tile_x * tile_size - 20037508.34
    ll_y = 20037508.34 - (tile_y + 1) * tile_size
    return ll_x, ll_y


def elevation2png(
    val: NDArray[np.floating],
    png_file: str,
    encoder: str = "terrarium",
    encoder_vmin: float = 0.0,
    encoder_vmax: float = 1.0,
    compress_level: int = 6,
) -> None:
    """Convert a 256x256 NumPy array to a PNG tile using a specified encoder.

    Parameters
    ----------
    val : NDArray[np.floating]
        256x256 array of elevation or data values.
    png_file : str
        Output PNG file path.
    encoder : str
        Encoding scheme. One of ``"terrarium"``, ``"terrarium16"``, ``"uint8"``,
        ``"uint16"``, ``"uint24"``, ``"uint32"``, ``"float8"``, ``"float16"``,
        ``"float24"``, ``"float32"``.
    encoder_vmin : float
        Minimum value for float encoders.
    encoder_vmax : float
        Maximum value for float encoders.
    compress_level : int
        PNG compression level (0-9).

    Raises
    ------
    ValueError
        If values exceed the range supported by the chosen encoder.
    """
    if encoder == "terrarium":
        rgb = np.zeros((256, 256, 3), "uint8")
        val += 32768.0
        rgb[:, :, 0] = np.floor(val / 256).astype(int)
        rgb[:, :, 1] = np.floor(val % 256)
        rgb[:, :, 2] = np.floor((val - np.floor(val)) * 256).astype(int)
    elif encoder == "terrarium16":
        rgb = np.zeros((256, 256, 3), "uint8")
        val += 32768.0
        rgb[:, :, 0] = np.floor(val / 256).astype(int)
        rgb[:, :, 1] = np.floor(val % 256).astype(int)
    elif encoder == "uint8":
        if np.any(val >= 255):
            raise ValueError(
                "Some values in are equal to or larger than 255. This is not allowed for encoder 'uint8'."
            )
        rgb = np.zeros((256, 256, 3), "uint8") + 255
        r = val + 0
        r[np.where(val < 0)] = 255
        rgb[:, :, 0] = r
    elif encoder == "uint16":
        if np.any(val >= 65535):
            raise ValueError(
                "Some values are equal to or larger than 65535. This is not allowed for encoder 'uint16'."
            )
        rgb = np.zeros((256, 256, 3), "uint8") + 255
        r = (val // 256) % 256
        g = val % 256
        r[np.where(val < 0)] = 255
        g[np.where(val < 0)] = 255
        rgb[:, :, 0] = r
        rgb[:, :, 1] = g
    elif encoder == "uint24":
        if np.any(val >= 16777215):
            raise ValueError(
                "Some values are equal to or larger than 16777215. This is not allowed for encoder 'uint24'."
            )
        rgb = np.zeros((256, 256, 3), "uint8") + 255
        r = (val // 256**2) % 256
        g = (val // 256) % 256
        b = val % 256
        r[np.where(val < 0)] = 255
        g[np.where(val < 0)] = 255
        b[np.where(val < 0)] = 255
        rgb[:, :, 0] = r
        rgb[:, :, 1] = g
        rgb[:, :, 2] = b
    elif encoder == "uint32":
        if np.any(val >= 4294967295):
            raise ValueError(
                "Some values are equal to or larger than 4294967295. This is not allowed for encoder 'uint32'."
            )
        rgb = np.zeros((256, 256, 4), "uint8") + 255
        r = (val // 256**3) % 256
        g = (val // 256**2) % 256
        b = (val // 256) % 256
        a = val % 256
        r[np.where(val < 0)] = 255
        g[np.where(val < 0)] = 255
        b[np.where(val < 0)] = 255
        a[np.where(val < 0)] = 255
        rgb[:, :, 0] = r
        rgb[:, :, 1] = g
        rgb[:, :, 2] = b
        rgb[:, :, 3] = a
    elif encoder == "float8":
        val = np.maximum(val, encoder_vmin)
        val = np.minimum(val, encoder_vmax)
        val = val - encoder_vmin
        i = np.floor(val * 254 / (encoder_vmax - encoder_vmin)).astype(int) + 1
        i[np.isnan(val)] = 0
        rgb = np.zeros((256, 256, 3), "uint8")
        rgb[:, :, 0] = i
    elif encoder == "float16":
        val = np.maximum(val, encoder_vmin)
        val = np.minimum(val, encoder_vmax)
        val = val - encoder_vmin
        i = np.floor(val * 65534 / (encoder_vmax - encoder_vmin)).astype(int) + 1
        i[np.isnan(val)] = 0
        rgb = np.zeros((256, 256, 3), "uint8")
        rgb[:, :, 0] = (i // 256) % 256
        rgb[:, :, 1] = i % 256
    elif encoder == "float24":
        val = np.maximum(val, encoder_vmin)
        val = np.minimum(val, encoder_vmax)
        val = val - encoder_vmin
        i = np.floor(val * 16777214 / (encoder_vmax - encoder_vmin)).astype(int) + 1
        i[np.isnan(val)] = 0
        rgb = np.zeros((256, 256, 3), "uint8")
        rgb[:, :, 0] = (i // 256**2) % 256
        rgb[:, :, 1] = (i // 256) % 256
        rgb[:, :, 2] = i % 256
    elif encoder == "float32":
        val = np.maximum(val, encoder_vmin)
        val = np.minimum(val, encoder_vmax)
        val = val - encoder_vmin
        i = np.floor(val * 4294967294 / (encoder_vmax - encoder_vmin)).astype(int) + 1
        i[np.isnan(val)] = 0
        rgb = np.zeros((256, 256, 4), "uint8")
        rgb[:, :, 0] = (i // 256**3) % 256
        rgb[:, :, 1] = (i // 256**2) % 256
        rgb[:, :, 2] = (i // 256) % 256
        rgb[:, :, 3] = i % 256

    img = Image.fromarray(rgb)
    img.save(png_file, compress_level=compress_level)


def png2elevation(
    png_file: str,
    encoder: str = "terrarium",
    encoder_vmin: float = 0.0,
    encoder_vmax: float = 1.0,
) -> NDArray[np.floating]:
    """Convert a PNG tile back to an elevation array using a specified encoder.

    Parameters
    ----------
    png_file : str
        Input PNG file path.
    encoder : str
        Encoding scheme used when the tile was created.
    encoder_vmin : float
        Minimum value for float encoders.
    encoder_vmax : float
        Maximum value for float encoders.

    Returns
    -------
    NDArray[np.floating]
        256x256 array of decoded elevation or data values.
    """
    img = Image.open(png_file)
    if encoder == "terrarium":
        rgb = np.array(img.convert("RGB")).astype(float)
        elevation = (rgb[:, :, 0] * 256 + rgb[:, :, 1] + rgb[:, :, 2] / 256) - 32768.0
        elevation[np.where(elevation < -32767.0)] = np.nan
    elif encoder == "terrarium16":
        rgb = np.array(img.convert("RGB")).astype(float)
        elevation = (rgb[:, :, 0] * 256 + rgb[:, :, 1]) - 32768.0
        elevation[np.where(elevation < -32767.0)] = np.nan
    elif encoder == "uint8":
        rgb = np.array(img.convert("RGB")).astype(int)
        elevation = rgb[:, :, 0]
        elevation[np.where(elevation == 255)] = -1
    elif encoder == "uint16":
        rgb = np.array(img.convert("RGB")).astype(int)
        elevation = rgb[:, :, 0] * 256 + rgb[:, :, 1]
        elevation[np.where(elevation == 65535)] = -1
    elif encoder == "uint24":
        rgb = np.array(img.convert("RGB")).astype(int)
        elevation = rgb[:, :, 0] * 65536 + rgb[:, :, 1] * 256 + rgb[:, :, 2]
        elevation[np.where(elevation == 16777215)] = -1
    elif encoder == "uint32":
        rgb = np.array(img.convert("RGBA")).astype(int)
        elevation = (
            rgb[:, :, 0] * 16777216
            + rgb[:, :, 1] * 65536
            + rgb[:, :, 2] * 256
            + rgb[:, :, 3]
        )
        elevation[np.where(elevation == 4294967295)] = -1
    elif encoder == "float8":
        rgb = np.array(img.convert("RGB")).astype(float)
        i = rgb[:, :, 0]
        elevation = encoder_vmin + (encoder_vmax - encoder_vmin) * i / 254
        elevation[np.where(i == 0)] = np.nan
    elif encoder == "float16":
        rgb = np.array(img.convert("RGB")).astype(float)
        i = rgb[:, :, 0] * 256 + rgb[:, :, 1]
        elevation = encoder_vmin + (encoder_vmax - encoder_vmin) * i / 65534
        elevation[np.where(i == 0)] = np.nan
    elif encoder == "float24":
        rgb = np.array(img.convert("RGB")).astype(float)
        i = rgb[:, :, 0] * 65536 + rgb[:, :, 1] * 256 + rgb[:, :, 2]
        elevation = encoder_vmin + (encoder_vmax - encoder_vmin) * i / 16777214
        elevation[np.where(i == 0)] = np.nan
    elif encoder == "float32":
        rgb = np.array(img.convert("RGBA")).astype(float)
        i = (
            rgb[:, :, 0] * 16777216
            + rgb[:, :, 1] * 65536
            + rgb[:, :, 2] * 256
            + rgb[:, :, 3]
        )
        elevation = encoder_vmin + (encoder_vmax - encoder_vmin) * i / 4294967294
        elevation[np.where(i == 0)] = np.nan
    return elevation


def png2int(png_file: str, idummy: int) -> NDArray[np.signedinteger]:
    """Convert a PNG tile to an integer index array.

    Parameters
    ----------
    png_file : str
        Input PNG file path.
    idummy : int
        Value to assign to pixels that decode as the maximum RGBA integer (no-data).

    Returns
    -------
    NDArray[np.signedinteger]
        256x256 integer index array.
    """
    image = Image.open(png_file)
    rgba = np.array(image.convert("RGBA")).astype(int)
    ind = (
        (rgba[:, :, 0] * 256**3)
        + (rgba[:, :, 1] * 256**2)
        + (rgba[:, :, 2] * 256)
        + rgba[:, :, 3]
    )
    ind[np.where(ind == 4294967295)] = idummy
    return ind


def int2png(val: NDArray[np.signedinteger], png_file: str) -> None:
    """Convert an integer index array to a PNG tile.

    Parameters
    ----------
    val : NDArray[np.signedinteger]
        256x256 integer array. Negative values are encoded as no-data (all 255).
    png_file : str
        Output PNG file path.
    """
    rgba = np.zeros((256, 256, 4), "uint8") + 255
    r = (val // 256**3) % 256
    g = (val // 256**2) % 256
    b = (val // 256) % 256
    a = val % 256
    r[np.where(val < 0)] = 255
    g[np.where(val < 0)] = 255
    b[np.where(val < 0)] = 255
    a[np.where(val < 0)] = 255
    rgba[:, :, 0] = r
    rgba[:, :, 1] = g
    rgba[:, :, 2] = b
    rgba[:, :, 3] = a
    img = Image.fromarray(rgba)
    img.save(png_file)


def makedir(path: str) -> None:
    """Create a directory (and parents) if it does not already exist.

    Parameters
    ----------
    path : str
        Directory path to create.
    """
    if not os.path.exists(path):
        os.makedirs(path)


def list_files(src: str) -> list[str]:
    """List files matching a glob pattern.

    Parameters
    ----------
    src : str
        Glob pattern (e.g. ``"/data/*.png"``).

    Returns
    -------
    list[str]
        List of matching file paths.
    """
    file_list = []
    full_list = glob.glob(src)
    for item in full_list:
        if os.path.isfile(item):
            file_list.append(item)
    return file_list


def list_folders(src: str, basename: bool = False) -> list[str]:
    """List directories matching a glob pattern.

    Parameters
    ----------
    src : str
        Glob pattern.
    basename : bool
        If True, return only the directory base names instead of full paths.

    Returns
    -------
    list[str]
        List of matching directory paths (or base names).
    """
    folder_list = []
    full_list = glob.glob(src)
    for item in full_list:
        if os.path.isdir(item):
            if basename:
                folder_list.append(os.path.basename(item))
            else:
                folder_list.append(item)

    return folder_list


def interp2(
    x0: NDArray[np.floating],
    y0: NDArray[np.floating],
    z0: NDArray[np.floating],
    x1: NDArray[np.floating],
    y1: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Bilinear interpolation from a regular grid onto scattered target points.

    Parameters
    ----------
    x0 : NDArray[np.floating]
        1-D x-coordinates of the source grid.
    y0 : NDArray[np.floating]
        1-D y-coordinates of the source grid.
    z0 : NDArray[np.floating]
        2-D source values, shape ``(len(y0), len(x0))``.
    x1 : NDArray[np.floating]
        2-D target x-coordinates.
    y1 : NDArray[np.floating]
        2-D target y-coordinates, same shape as *x1*.

    Returns
    -------
    NDArray[np.floating]
        Interpolated values at target locations, same shape as *x1*.
    """
    f = RegularGridInterpolator((y0, x0), z0, bounds_error=False, fill_value=np.nan)
    sz = x1.shape
    x1 = x1.reshape(sz[0] * sz[1])
    y1 = y1.reshape(sz[0] * sz[1])
    z1 = f((y1, x1)).reshape(sz)
    return z1


def binary_search(
    val_array: NDArray[np.number], vals: NDArray[np.number]
) -> NDArray[np.signedinteger]:
    """Find indices of *vals* within a sorted *val_array* using binary search.

    Parameters
    ----------
    val_array : NDArray[np.number]
        Sorted 1-D array of reference values.
    vals : NDArray[np.number]
        Values to look up.

    Returns
    -------
    NDArray[np.signedinteger]
        Index into *val_array* for each element of *vals*, or -1 if not found.
    """
    indx = np.searchsorted(val_array, vals)
    not_ok = np.where(indx == len(val_array))[0]
    indx[np.where(indx == len(val_array))[0]] = 0
    is_ok = np.where(val_array[indx] == vals)[0]
    indices = np.zeros(len(vals), dtype=int) - 1
    indices[is_ok] = indx[is_ok]
    indices[not_ok] = -1
    return indices
