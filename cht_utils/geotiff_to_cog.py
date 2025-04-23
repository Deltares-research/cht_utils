import os

import rioxarray
import rasterio
from pyproj import CRS

def is_cog(tif_path):
    """
    Check if a GeoTIFF is a Cloud-Optimized GeoTIFF (COG).

    Parameters:
        tif_path (str): Path to the GeoTIFF file.

    Returns:
        bool: True if the file is a COG, False otherwise.
    """
    try:
        with rasterio.open(tif_path) as src:
            # Check if the file is tiled
            if not src.is_tiled:
                return False

            # Check for overviews (COGs should have overviews)
            if not src.overviews(1):
                return False

            # Check for proper metadata (COGs have 'COG' in their metadata)
            metadata = src.tags()
            if "COG" not in metadata.get("TIFFTAG_IMAGELENGTH", ""):
                return False

            # If all checks pass, the file is a COG
            return True

    except Exception as e:
        print(f"Error checking COG: {e}")
        return False



def geotiff_to_cog(geotiff_path, output_cog_path, resampling="average"):
    """
    Convert a GeoTIFF to a Cloud-Optimized GeoTIFF (COG) using rioxarray.

    Parameters:
        geotiff_path (str): Path to the input GeoTIFF file.
        output_cog_path (str): Path to the output COG file.
        resampling (str): Resampling method for overviews (default: 'average').
    """
    try:

        # First check if the file is already a COG
        if is_cog(geotiff_path):
            print(f"{geotiff_path} is already a COG ! Copying to {output_cog_path}")
            # Copy the file to the new location
            os.path.copy(geotiff_path, output_cog_path)
            return True

        # Load GeoTIFF as a DataArray with spatial referencing
        da = rioxarray.open_rasterio(geotiff_path, masked=True)

        # Get CRS from file (default to EPSG:4326 if not found)
        crs = da.rio.crs
        if crs is None:
            print("No CRS found in source; defaulting to EPSG:4326")
            crs = CRS.from_epsg(4326)
            da = da.rio.write_crs(crs)

        # Write to COG
        da.rio.to_raster(
            output_cog_path,
            driver="COG",
            compress="deflate",
            blocksize=512,
            overview_resampling=resampling,
            dtype=str(da.dtype)
        )

        print(f"COG saved to: {output_cog_path}")
        return True

    except Exception as e:
        print(f"Error during conversion: {e}")
        return False
