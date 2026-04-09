"""Convert XYZ point-cloud data to Cloud-Optimized GeoTIFF (COG) format."""

import os
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS
from rasterio.enums import Resampling
from rasterio.transform import from_origin


def xyz_to_cog(
    xyz_path: str,
    output_cog_path: str,
    crs: Union[int, CRS] = 4326,
) -> bool:
    """Grid scattered XYZ data and write as a COG.

    Parameters
    ----------
    xyz_path : str
        Path to a whitespace-separated XYZ file (columns: x, y, z).
    output_cog_path : str
        Output COG file path.
    crs : int or pyproj.CRS
        Coordinate reference system (default: EPSG:4326).

    Returns
    -------
    bool
        ``True`` on success, ``False`` on error.
    """
    try:
        df = pd.read_csv(xyz_path, sep=r"\s+", names=["x", "y", "z"])

        x_unique = np.sort(df["x"].unique())
        y_unique = np.sort(df["y"].unique())

        dx = np.round(np.min(np.diff(x_unique)), 6)
        dy = np.round(np.min(np.diff(y_unique)), 6)

        z_grid = np.full((len(y_unique), len(x_unique)), np.nan)
        x_to_idx = {x: i for i, x in enumerate(x_unique)}
        y_to_idx = {y: i for i, y in enumerate(y_unique)}

        for row in df.itertuples(index=False):
            z_grid[y_to_idx[row.y], x_to_idx[row.x]] = row.z

        da = xr.DataArray(
            data=z_grid,
            dims=["y", "x"],
            coords={"x": x_unique, "y": y_unique},
            name="topography",
        )

        if y_unique[0] > y_unique[-1]:
            da = da.sortby("y")

        transform = from_origin(
            x_unique[0] - dx / 2, y_unique[-1] + dy / 2, dx, dy
        )
        da.rio.write_transform(transform, inplace=True)

        if isinstance(crs, int):
            crs_string = f"EPSG:{crs}"
        elif isinstance(crs, CRS):
            crs_string = f"EPSG:{crs.to_epsg()}"
        else:
            crs_string = str(crs)

        da.rio.write_crs(crs_string, inplace=True)
        da.rio.to_raster(
            output_cog_path,
            driver="COG",
            compress="deflate",
            dtype="float32",
            nodata=np.nan,
            resampling=Resampling.average,
        )

        print(f"Saved to: {os.path.abspath(output_cog_path)}")
        return True

    except Exception as e:
        print(f"Error: {e}")
        return False
