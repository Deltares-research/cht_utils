import xarray as xr
import rioxarray


# Load the NetCDF file
ds = xr.open_dataset("D:/metocean_data/open/NGDC/crm_vol2_2023.nc")

# Select the variable you want to export
# Replace 'your_variable' with the actual variable name
da = ds['z']

# If necessary, select a specific time slice or level
# da = da.isel(time=0)  # if there's a time dimension

# Ensure the data has a CRS and spatial coordinates set
da.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
da.rio.write_crs("EPSG:4326", inplace=True)  # or your appropriate CRS

# Export to a GeoTIFF (initial step before creating a COG)
da.rio.to_raster("D:/metocean_data/open/NGDC/crm_vol2_2023.tif")

# Convert the GeoTIFF to a COG
from rio_cogeo.cogeo import cog_translate
from rio_cogeo.profiles import cog_profiles
from rasterio import Env
from rasterio.shutil import copy as rio_copy

# Define profile
profile = cog_profiles.get("deflate")

# Translate to COG
with Env():
    cog_translate(
        "D:/metocean_data/open/NGDC/crm_vol2_2023.tif",
        "D:/metocean_data/open/NGDC/crm_vol2_2023_cog.tif",
        profile,
        in_memory=False,
    )
