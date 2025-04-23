import numpy as np
import pandas as pd
import xarray as xr
# import rioxarray
from rasterio.transform import from_origin
from rasterio.enums import Resampling
import os
from pyproj import CRS
# import matplotlib.pyplot as plt
# import contextily as ctx

def xyz_to_cog(xyz_path, output_cog_path, crs=4326):

    try:

        # === Step 1: Load the XYZ data ===
        df = pd.read_csv(xyz_path, sep=r'\s+', names=["x", "y", "z"])

        # === Step 2: Create grid ===
        x_unique = np.sort(df["x"].unique())
        y_unique = np.sort(df["y"].unique())

        dx = np.round(np.min(np.diff(x_unique)), 6)
        dy = np.round(np.min(np.diff(y_unique)), 6)

        # Create 2D grid filled with NaNs
        z_grid = np.full((len(y_unique), len(x_unique)), np.nan)

        x_to_idx = {x: i for i, x in enumerate(x_unique)}
        y_to_idx = {y: i for i, y in enumerate(y_unique)}

        for row in df.itertuples(index=False):
            x_idx = x_to_idx[row.x]
            y_idx = y_to_idx[row.y]
            z_grid[y_idx, x_idx] = row.z

        # === Step 3: Create xarray DataArray ===
        da = xr.DataArray(
            data=z_grid,
            dims=["y", "x"],
            coords={"x": x_unique, "y": y_unique},
            name="topography"
        )

        # Optional: flip y if needed to make north up
        if y_unique[0] > y_unique[-1]:
            da = da.sortby("y")

        # === Step 4: Set spatial attributes with rioxarray ===
        # Use the upper-left corner for the origin
        transform = from_origin(x_unique[0] - dx / 2, y_unique[-1] + dy / 2, dx, dy)

        da.rio.write_transform(transform, inplace=True)

        # Assign CRS (replace this string with your actual CRS, e.g., "EPSG:32633")
        if isinstance(crs, int):
            crs_string = f"EPSG:{crs}"
        elif isinstance(crs, CRS):
            # if crs is a pyproj CRS object, convert to string
            crs_string = f"EPSG:{crs.to_epsg()}"

        da.rio.write_crs(crs_string, inplace=True)

        # === Step 5: Export as Cloud Optimized GeoTIFF ===

        da.rio.to_raster(
            output_cog_path,
            driver="COG",
            compress="deflate",
            dtype="float32",
            nodata=np.nan,
            resampling=Resampling.average
        )

        print(f"Saved to: {os.path.abspath(output_cog_path)}")

        return True

        # # === Step 1: Read the COG ===
        # da = rioxarray.open_rasterio("topography_cog.tif", masked=True).squeeze()

        # # === Step 2: Reproject to Web Mercator (EPSG:3857) for OSM tiles ===
        # da_web_mercator = da.rio.reproject("EPSG:3857")

        # # === Step 3: Plot on top of OSM background ===
        # fig, ax = plt.subplots(figsize=(10, 8))

        # img = da_web_mercator.plot(
        #     ax=ax,
        #     cmap="terrain",
        #     alpha=0.7,
        #     cbar_kwargs={"label": "Elevation (m)"}
        # )

        # # Add OpenStreetMap background
        # ctx.add_basemap(ax, crs=da_web_mercator.rio.crs, source=ctx.providers.OpenStreetMap.Mapnik)

        # # Set axis labels (optional)
        # ax.set_xlabel("Easting (m)")
        # ax.set_ylabel("Northing (m)")

        # plt.title("Topography over OpenStreetMap")
        # plt.tight_layout()
        # plt.show()

    except Exception as e:
        print(f"Error: {e}")
        return False
       
