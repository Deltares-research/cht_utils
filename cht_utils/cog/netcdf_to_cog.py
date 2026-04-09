"""Convert NetCDF variables to Cloud-Optimized GeoTIFF (COG) format."""

from typing import Optional

import xarray as xr
from pyproj import CRS


def netcdf_to_cog(
    netcdf_path: str,
    variable_name: str,
    output_cog_path: str,
    time_index: Optional[int] = None,
) -> bool:
    """Convert a NetCDF variable to a Cloud-Optimized GeoTIFF.

    Assumes EPSG:4326 if no CRS information is found in the dataset.

    Parameters
    ----------
    netcdf_path : str
        Path to the NetCDF file.
    variable_name : str
        Name of the variable to extract.
    output_cog_path : str
        Output path for the COG.
    time_index : int or None
        Optional time dimension index to select.

    Returns
    -------
    bool
        ``True`` on success, ``False`` on error.
    """
    try:
        ds = xr.open_dataset(netcdf_path)

        if variable_name not in ds:
            raise ValueError(
                f"Variable '{variable_name}' not found in the NetCDF file."
            )
        da = ds[variable_name]

        if "time" in da.dims and time_index is not None:
            da = da.isel(time=time_index)

        crs = _get_crs(ds)
        crs_str = f"EPSG:{crs.to_epsg()}"

        if crs.is_projected:
            da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
        else:
            da.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)

        da.rio.write_crs(crs_str, inplace=True)
        da.rio.to_raster(
            output_cog_path,
            driver="COG",
            compress="deflate",
            blocksize=512,
            overview_resampling="average",
            dtype=da.dtype,
        )
        print(f"COG saved to: {output_cog_path}")
        return True

    except Exception as e:
        print(f"Error during conversion: {e}")
        return False


def _get_crs(ds: xr.Dataset) -> CRS:
    """Extract CRS from a dataset, defaulting to EPSG:4326."""
    if "crs" not in ds:
        return CRS.from_epsg(4326)

    for attr_name in ("crs_wkt", "spatial_ref"):
        wkt = ds["crs"].attrs.get(attr_name)
        if wkt:
            crs = CRS.from_wkt(wkt)
            if crs.to_epsg() is None and "NAD83" in crs.name:
                return CRS.from_epsg(4269)
            return crs

    return CRS.from_epsg(4326)
