import xarray as xr
from pyproj import CRS

def netcdf_to_cog(netcdf_path, variable_name, output_cog_path, time_index=None):
    """
    Convert a NetCDF variable to a Cloud-Optimized GeoTIFF (COG).

    Parameters:
        netcdf_path (str): Path to the NetCDF file.
        variable_name (str): Name of the variable to extract.
        output_cog_path (str): Output path for the COG.
        time_index (int or None): Optional index for time dimension if present.
    """

    try:

        # Load dataset
        ds = xr.open_dataset(netcdf_path)

        # Extract variable
        if variable_name not in ds:
            raise ValueError(f"Variable '{variable_name}' not found in the NetCDF file.")
        da = ds[variable_name]

        # Handle time slicing if needed
        if "time" in da.dims and time_index is not None:
            da = da.isel(time=time_index)

        # Check for CRS
        if "crs" not in ds:
            # Assume 4326 if no CRS is found
            crs = CRS.from_epsg(4326)
        else:
            # Get the CRS from the dataset
            wkt_options = ["crs_wkt", "spatial_ref"]
            wkt = None
            for option in wkt_options:
                if option in ds["crs"].attrs:
                    wkt = ds["crs"].attrs[option]
                    break

            if wkt is None:
                # Assume 4326 if no WKT is found
                crs = CRS.from_epsg(4326)    
            else:
                crs = CRS.from_wkt(wkt)
                if crs.to_epsg() is None:
                    if "NAD83" in crs.name:
                        # Handle NAD83 case
                        crs = CRS.from_epsg(4269)

        crs_str = f"EPSG:{crs.to_epsg()}"

        # Set spatial dimensions
        # Check if crs is projected or geographic
        if crs.is_projected:
            da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
        else:
            # For geographic coordinates, use lon/lat
            da.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)

        # Set CRS (modify as needed)
        da.rio.write_crs(crs_str, inplace=True)

        # And write to COG 
        da.rio.to_raster(
            output_cog_path,
            driver="COG",
            compress="deflate",       # or "lzw" (optional compression)
            blocksize=512,            # typical block size (optional)
            overview_resampling="average",  # build internal overviews
            dtype=da.dtype            # preserve original data type
        )
        print(f"COG saved to: {output_cog_path}")

        return True
    
    except Exception as e:
        print(f"Error during conversion: {e}")
        return False
