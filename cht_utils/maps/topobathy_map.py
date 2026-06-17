"""Topography/bathymetry map generation from Cloud Optimized GeoTIFF (COG) files.

Provides the TopoBathyMap class for reading elevation data, applying color
maps, writing output files (GeoTIFF/NetCDF), creating map overlays, and
plotting with matplotlib.
"""

import logging
from pathlib import Path

import contextily as ctx
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import rioxarray
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from pyproj import Transformer
from rasterio.warp import Resampling

from cht_utils.maps.flood_map import (
    get_appropriate_overview_level,
    get_rgb_data_array,
    reproject_bbox,
)

__all__ = ["TopoBathyMap"]

logger = logging.getLogger(__name__)


class TopoBathyMap:
    """Manages reading, processing, and visualising topography/bathymetry data.

    Wraps a Cloud Optimized GeoTIFF elevation file and provides methods for
    reading at various overview levels, applying elevation masks, writing
    output as GeoTIFF or NetCDF, creating PNG map overlays, and producing
    matplotlib plots.

    Parameters
    ----------
    topobathy_file : str | Path | None
        Path to the topobathy COG file.
    zbmin : float
        Minimum allowable elevation value; values below are masked to NaN.
    zbmax : float
        Maximum allowable elevation value; values above are masked to NaN.
    max_pixel_size : float
        Maximum pixel size for selecting the appropriate overview level.
    data_array_name : str
        Name of the elevation variable in the output dataset.
    cmap : str | None
        Matplotlib colormap name for coloring.
    cmin : float | None
        Minimum value for colormap normalization (scaled by ``scale_factor``).
    cmax : float | None
        Maximum value for colormap normalization (scaled by ``scale_factor``).
    color_values : list[dict] | None
        Discrete color definitions with ``lower_value``, ``upper_value``,
        and ``color`` keys.
    scale_factor : float
        Multiplicative factor applied to elevation values on read.
    """

    def __init__(
        self,
        topobathy_file: str | Path | None = None,
        zbmin: float = -999999.9,
        zbmax: float = 99999.9,
        max_pixel_size: float = 0.0,
        data_array_name: str = "elevation",
        cmap: str | None = None,
        cmin: float | None = None,
        cmax: float | None = None,
        color_values: list[dict] | None = None,
        scale_factor: float = 1.0,
    ) -> None:
        self.topobathy_file = topobathy_file
        self.zb = rasterio.open(topobathy_file) if topobathy_file else None
        self.zbmin = zbmin
        self.zbmax = zbmax
        self.max_pixel_size = max_pixel_size
        self.data_array_name = data_array_name
        self.cmap = cmap
        self.cmin = cmin * scale_factor
        self.cmax = cmax * scale_factor
        self.color_values = color_values
        self.scale_factor = scale_factor
        self.ds = xr.Dataset()

    def set_topobathy_file(self, topobathy_file: str | Path) -> None:
        """Set the topobathy file and open it with rasterio.

        Parameters
        ----------
        topobathy_file : str | Path
            Path to the topobathy COG file.
        """
        self.topobathy_file = topobathy_file
        self.zb = rasterio.open(self.topobathy_file)

    def close(self) -> None:
        """Close the rasterio dataset and the xarray dataset."""
        if self.zb is not None:
            self.zb.close()
        self.ds.close()

    def read(self, tiffile: str | Path) -> None:
        """Read a GeoTIFF file with elevation data into the dataset.

        Parameters
        ----------
        tiffile : str | Path
            Path to the GeoTIFF file.
        """
        self.ds = xr.Dataset()
        self.ds[self.data_array_name] = rioxarray.open_rasterio(
            tiffile, masked=True
        ).squeeze()

    def make(
        self,
        max_pixel_size: float = 0.0,
        bbox: tuple[float, float, float, float] | None = None,
    ) -> xr.Dataset:
        """Generate a topobathy dataset from the elevation COG.

        Reads the data at an appropriate overview level, applies the scale
        factor, clips to the bounding box, and masks values outside
        ``[zbmin, zbmax]``.

        Parameters
        ----------
        max_pixel_size : float
            Maximum pixel size in metres for overview selection. If 0.0,
            the native resolution is used.
        bbox : tuple[float, float, float, float] | None
            Bounding box ``(minx, miny, maxx, maxy)`` to clip the data.

        Returns
        -------
        xr.Dataset
            Dataset containing the masked elevation data array.
        """
        overview_level = 0

        if max_pixel_size > 0.0:
            overview_level = get_appropriate_overview_level(self.zb, max_pixel_size)

        if overview_level == 0:
            zb = rioxarray.open_rasterio(self.zb)
        else:
            zb = rioxarray.open_rasterio(self.zb, overview_level=overview_level)

        zb = zb * self.scale_factor

        if "band" in zb.dims and zb.sizes["band"] == 1:
            zb = zb.squeeze(dim="band", drop=True)

        if bbox is not None:
            zb = zb.rio.clip_box(minx=bbox[0], miny=bbox[1], maxx=bbox[2], maxy=bbox[3])

        elevation = zb.to_numpy()[:]
        elevation[elevation < self.zbmin] = np.nan
        elevation[elevation > self.zbmax] = np.nan

        self.ds = xr.Dataset()
        self.ds[self.data_array_name] = xr.DataArray(
            elevation, dims=["y", "x"], coords={"y": zb.y, "x": zb.x}
        )
        self.ds[self.data_array_name] = self.ds[self.data_array_name].rio.write_crs(
            zb.rio.crs, inplace=True
        )

        return self.ds

    def write(self, output_file: str | Path = "") -> None:
        """Write the topobathy data to a GeoTIFF or NetCDF file.

        Parameters
        ----------
        output_file : str | Path
            Output file path. Extension determines the format:
            ``".tif"`` for COG GeoTIFF, ``".nc"`` for NetCDF.
        """
        if output_file.endswith(".nc"):
            self.ds.to_netcdf(output_file)

        elif output_file.endswith(".tif"):
            if self.cmap is not None:
                rgb_da = get_rgb_data_array(
                    self.ds[self.data_array_name],
                    cmap=self.cmap,
                    cmin=self.cmin,
                    cmax=self.cmax,
                    color_values=self.color_values,
                )

                rgb_da.rio.to_raster(
                    output_file,
                    driver="COG",
                    compress="deflate",
                    blocksize=512,
                    overview_resampling="nearest",
                )

            else:
                self.ds[self.data_array_name].rio.to_raster(
                    output_file,
                    driver="COG",
                    compress="deflate",
                    blocksize=512,
                    overview_resampling="nearest",
                )

    def map_overlay(
        self,
        file_name: str,
        xlim: list[float] | None = None,
        ylim: list[float] | None = None,
        width: int = 800,
    ) -> bool:
        """Create a PNG map overlay of the topobathy data in EPSG:3857.

        Parameters
        ----------
        file_name : str
            Output PNG file path.
        xlim : list[float] | None
            Longitude extent ``[lon_min, lon_max]``.
        ylim : list[float] | None
            Latitude extent ``[lat_min, lat_max]``.
        width : int
            Width in pixels for resolution calculation.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        if self.ds is None:
            return False

        try:
            lon_min = xlim[0]
            lat_min = ylim[0]
            lon_max = xlim[1]
            lat_max = ylim[1]

            transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
            x_min, y_min = transformer.transform(lon_min, lat_min)
            x_max, y_max = transformer.transform(lon_max, lat_max)

            dxy = (x_max - x_min) / width

            bbox = reproject_bbox(
                lon_min,
                lat_min,
                lon_max,
                lat_max,
                crs_src="EPSG:4326",
                crs_dst=self.zb.crs,
                buffer=0.05,
            )

            self.make(max_pixel_size=dxy, bbox=bbox)

            rgb_da = get_rgb_data_array(
                self.ds[self.data_array_name],
                cmap=self.cmap,
                cmin=self.cmin,
                cmax=self.cmax,
                color_values=self.color_values,
            )

            rgb_3857 = rgb_da.rio.reproject(
                "EPSG:3857", resampling=Resampling.bilinear, nodata=0
            )

            rgb_3857 = rgb_3857.rio.pad_box(
                minx=x_min, miny=y_min, maxx=x_max, maxy=y_max, constant_values=0
            )

            rgb_crop = rgb_3857.rio.clip_box(
                minx=x_min, miny=y_min, maxx=x_max, maxy=y_max
            )

            rgba = rgb_crop.transpose("y", "x", "band").to_numpy().astype("uint8")

            plt.imsave(file_name, rgba)

            return True

        except Exception as e:
            logger.exception(f"Error in map_overlay: {e}")
            return False

    def plot(
        self,
        pngfile: str,
        zoom: int | None = None,
        title: str = "Elevation (m)",
        color_values: list[dict] | None = None,
        cmap: str = "terrain",
        vmin: float | None = None,
        vmax: float | None = None,
        lon_lim: list[float] | None = None,
        lat_lim: list[float] | None = None,
        width: float = 10.0,
        background: str = "EsriWorldImagery",
    ) -> None:
        """Plot the topobathy data with a basemap and save to PNG.

        Parameters
        ----------
        pngfile : str
            Output PNG file path.
        zoom : int | None
            Basemap zoom level. If None, auto-detected.
        title : str
            Plot title.
        color_values : list[dict] | None
            Discrete color definitions. If a string is passed, a default
            elevation color scheme is used.
        cmap : str
            Matplotlib colormap for continuous coloring.
        vmin : float | None
            Minimum value for color mapping. Auto-detected if None.
        vmax : float | None
            Maximum value for color mapping. Auto-detected if None.
        lon_lim : list[float] | None
            Longitude limits ``[lon_min, lon_max]`` for the plot extent.
        lat_lim : list[float] | None
            Latitude limits ``[lat_min, lat_max]`` for the plot extent.
        width : float
            Figure width in inches.
        background : str
            Basemap provider: ``"osm"`` for OpenStreetMap or
            ``"EsriWorldImagery"`` (default).
        """
        if lon_lim is None or lat_lim is None:
            lon_min = self.ds.x.min().to_numpy()
            lat_min = self.ds.y.min().to_numpy()
            lon_max = self.ds.x.max().to_numpy()
            lat_max = self.ds.y.max().to_numpy()
            crs = self.ds[self.data_array_name].rio.crs
            transformer = Transformer.from_crs(crs, "EPSG:3857", always_xy=True)
            x_min, y_min = transformer.transform(lon_min, lat_min)
            x_max, y_max = transformer.transform(lon_max, lat_max)
        else:
            transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
            x_min, y_min = transformer.transform(lon_lim[0], lat_lim[0])
            x_max, y_max = transformer.transform(lon_lim[1], lat_lim[1])

        da_3857 = self.ds[self.data_array_name].rio.reproject("EPSG:3857")

        if vmin is None:
            vmin = float(da_3857.min())
        if vmax is None:
            vmax = float(da_3857.max())

        if color_values is None:
            discrete_colors = False
        else:
            discrete_colors = True
            if isinstance(color_values, str):
                color_values = []
                color_values.append(
                    {"color": "blue", "lower_value": -100, "upper_value": 0}
                )
                color_values.append(
                    {"color": "lightblue", "lower_value": 0, "upper_value": 10}
                )
                color_values.append(
                    {"color": "green", "lower_value": 10, "upper_value": 50}
                )
                color_values.append(
                    {"color": "brown", "lower_value": 50, "upper_value": 100}
                )
                color_values.append({"color": "white", "lower_value": 100})

        aspect_ratio = (y_max - y_min) / (x_max - x_min)
        fig, ax = plt.subplots(figsize=(width, aspect_ratio * width))

        if discrete_colors:
            masked = da_3857.where(da_3857 >= color_values[0]["lower_value"])

            classified = xr.full_like(masked, np.nan)
            colors = []
            labels = []
            for icolor, color_value in enumerate(color_values):
                if "upper_value" in color_value:
                    lv = color_value["lower_value"]
                    uv = color_value["upper_value"]
                    classified = classified.where(
                        ~((masked > lv) & (masked <= uv)), icolor + 1
                    )
                    labels.append(f"{lv}--{uv} m")
                else:
                    lv = color_value["lower_value"]
                    classified = classified.where(~(masked > lv), icolor + 1)
                    labels.append(f">{lv} m")
                colors.append(color_value["color"])

            cmap_obj = ListedColormap(colors)
            bounds = list(range(1, len(colors) + 2))
            norm = BoundaryNorm(bounds, cmap_obj.N)

            classified.plot(ax=ax, cmap=cmap_obj, norm=norm, add_colorbar=False)

            legend_elements = []
            for i, color_value in enumerate(color_values):
                legend_elements.append(
                    Patch(facecolor=color_value["color"], label=labels[i])
                )
            plt.legend(handles=legend_elements, title="Elevation", loc="lower right")

        else:
            da_3857.plot(
                ax=ax,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                add_colorbar=True,
                cbar_kwargs={"label": "Elevation (m)"},
                alpha=0.75,
            )

        if background.lower() == "osm":
            if zoom is None:
                ctx.add_basemap(
                    ax, crs=da_3857.rio.crs, source=ctx.providers.OpenStreetMap.Mapnik
                )
            else:
                ctx.add_basemap(
                    ax,
                    crs=da_3857.rio.crs,
                    source=ctx.providers.OpenStreetMap.Mapnik,
                    zoom=zoom,
                )
        else:
            if zoom is None:
                ctx.add_basemap(
                    ax, crs=da_3857.rio.crs, source=ctx.providers.Esri.WorldImagery
                )
            else:
                ctx.add_basemap(
                    ax,
                    crs=da_3857.rio.crs,
                    source=ctx.providers.Esri.WorldImagery,
                    zoom=zoom,
                )

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        ax.set_axis_off()
        plt.title(title)

        plt.tight_layout()
        plt.savefig(pngfile, dpi=300, bbox_inches="tight", pad_inches=0.1)
        plt.close()
