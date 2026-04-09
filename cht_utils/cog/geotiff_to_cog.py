"""Convert GeoTIFF files to Cloud-Optimized GeoTIFF (COG) format."""

import shutil

import rasterio
from rasterio.shutil import copy as rio_copy


def is_cog(tif_path: str) -> bool:
    """Check if a GeoTIFF is Cloud-Optimized.

    Verifies that the file is tiled and has overviews.

    Parameters
    ----------
    tif_path : str
        Path to the GeoTIFF file.

    Returns
    -------
    bool
        ``True`` if the file is a COG.
    """
    try:
        with rasterio.open(tif_path) as src:
            if not src.is_tiled:
                return False
            overviews = src.overviews(1)
            return len(overviews) > 0
    except Exception as e:
        print(f"Error checking COG: {e}")
        return False


def geotiff_to_cog(
    geotiff_path: str,
    output_cog_path: str,
    resampling: str = "average",
) -> bool:
    """Convert a GeoTIFF to Cloud-Optimized GeoTIFF.

    If the input is already a COG, it is simply copied.

    Parameters
    ----------
    geotiff_path : str
        Path to the input GeoTIFF.
    output_cog_path : str
        Path for the output COG.
    resampling : str
        Resampling method for overviews (default: ``"average"``).

    Returns
    -------
    bool
        ``True`` on success, ``False`` on error.
    """
    try:
        if is_cog(geotiff_path):
            print(f"{geotiff_path} is already a COG! Copying to {output_cog_path}")
            shutil.copy(geotiff_path, output_cog_path)
            return True

        with rasterio.open(geotiff_path) as src:
            rio_copy(
                src,
                output_cog_path,
                driver="COG",
                compress="deflate",
                blocksize=512,
                overview_resampling=resampling,
            )

        print(f"COG saved to: {output_cog_path}")
        return True

    except Exception as e:
        print(f"Error during conversion: {e}")
        return False
