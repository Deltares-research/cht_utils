import xarray as xr
import rioxarray
from rio_cogeo.cogeo import cog_translate
from rio_cogeo.profiles import cog_profiles
from rasterio import Env

def netcdf_to_cog(netcdf_path, variable_name, output_cog_path, time_index=None):
    """
    Convert a NetCDF variable to a Cloud-Optimized GeoTIFF (COG).

    Parameters:
        netcdf_path (str): Path to the NetCDF file.
        variable_name (str): Name of the variable to extract.
        output_cog_path (str): Output path for the COG.
        time_index (int or None): Optional index for time dimension if present.
    """
    # Load dataset
    ds = xr.open_dataset(netcdf_path)

    # Extract variable
    if variable_name not in ds:
        raise ValueError(f"Variable '{variable_name}' not found in the NetCDF file.")
    da = ds[variable_name]

    # Handle time slicing if needed
    if "time" in da.dims and time_index is not None:
        da = da.isel(time=time_index)

    # Set spatial dimensions
    da.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)

    # Set CRS (modify as needed)
    da.rio.write_crs("EPSG:4326", inplace=True)

    # Save to temporary GeoTIFF
    temp_tif = "temp_output.tif"
    da.rio.to_raster(temp_tif)

    # Convert to Cloud-Optimized GeoTIFF
    profile = cog_profiles.get("deflate")
    with Env():
        cog_translate(
            temp_tif,
            output_cog_path,
            profile,
            in_memory=False,
        )

    print(f"COG saved to: {output_cog_path}")
