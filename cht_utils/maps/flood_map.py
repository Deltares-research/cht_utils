"""Flood map generation from water level data, topobathy, and index rasters.

Provides the FloodMap class for computing flood depth from water levels and
topography, writing output as GeoTIFF or NetCDF, creating map overlays, and
plotting. Also includes legacy tile-based flood map generation functions and
shared utilities for overview level selection, RGB conversion, and bounding
box reprojection.
"""

import logging
import os
from pathlib import Path

import contextily as ctx
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import rioxarray
import xarray as xr
from matplotlib import cm
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from PIL import Image
from pyproj import Transformer
from rasterio.warp import Resampling

import cht_utils.maps.fileops as fo
from cht_utils.maps.utils import deg2num, num2deg, png2elevation, png2int

logger = logging.getLogger(__name__)

__all__ = [
    "FloodMap",
    "get_appropriate_overview_level",
    "get_rgb_data_array",
    "reproject_bbox",
    "make_flood_map_tiles",
    "make_flood_map_overlay_v2",
]


class FloodMap:
    """Compute and visualise flood depth maps from water level and topobathy data.

    Uses Cloud Optimized GeoTIFF (COG) files for topography and cell
    indices, and combines them with water level arrays to produce flood
    depth grids. Supports writing output as GeoTIFF/NetCDF, creating
    PNG map overlays, and matplotlib plotting.

    Parameters
    ----------
    topobathy_file : str | Path | None
        Path to the topobathy COG file.
    index_file : str | Path | None
        Path to the cell-index COG file.
    zbmin : float
        Minimum allowable topobathy value; below this is masked.
    zbmax : float
        Maximum allowable topobathy value; above this is masked.
    hmin : float
        Minimum water depth threshold; shallower areas are masked.
    max_pixel_size : float
        Maximum pixel size for overview level selection.
    data_array_name : str
        Name of the depth variable in the output dataset.
    cmap : str | None
        Matplotlib colormap name.
    cmin : float | None
        Minimum value for colormap normalization.
    cmax : float | None
        Maximum value for colormap normalization.
    color_values : list[dict] | None
        Discrete color definitions with ``lower_value``, ``upper_value``,
        and ``color``/``rgb`` keys.
    """

    def __init__(
        self,
        topobathy_file: str | Path | None = None,
        index_file: str | Path | None = None,
        zbmin: float = 0.0,
        zbmax: float = 99999.9,
        hmin: float = 0.1,
        max_pixel_size: float = 0.0,
        data_array_name: str = "water_depth",
        cmap: str | None = None,
        cmin: float | None = None,
        cmax: float | None = None,
        color_values: list[dict] | None = None,
    ) -> None:
        self.topobathy_file = None
        self.index_file = None
        self.zb = None
        self.indices = None
        self.zbmin = zbmin
        self.zbmax = zbmax
        self.hmin = hmin
        self.max_pixel_size = max_pixel_size
        self.data_array_name = data_array_name
        self.color_values = color_values if color_values is not None else "default"
        self.cmap = cmap if cmap is not None else "jet"
        self.cmin = cmin if cmin is not None else 0.0
        self.cmax = cmax if cmax is not None else 1.0
        self.discrete_colors = color_values is not None
        self.ds = xr.Dataset()

        if topobathy_file is not None:
            self.set_topobathy_file(topobathy_file)
        if index_file is not None:
            self.set_index_file(index_file)

        self.legend = {}
        self.legend["title"] = "Flood Depth (m)"
        self.legend["contour"] = []
        self.legend["contour"].append(
            {"color": "#FF0000", "lower_value": 2.0, "text": "2.0+ m"}
        )
        self.legend["contour"].append(
            {
                "color": "#FFA500",
                "lower_value": 1.0,
                "upper_value": 2.0,
                "text": "1.0--2.0 m",
            }
        )
        self.legend["contour"].append(
            {
                "color": "#FFFF00",
                "lower_value": 0.3,
                "upper_value": 1.0,
                "text": "0.3--1.0 m",
            }
        )
        self.legend["contour"].append(
            {
                "color": "#00FF00",
                "lower_value": 0.1,
                "upper_value": 0.3,
                "text": "0.1--0.3 m",
            }
        )

    def set_topobathy_file(self, topobathy_file: str | Path) -> None:
        """Set the topobathy file and open it with rasterio.

        Parameters
        ----------
        topobathy_file : str | Path
            Path to the topobathy COG file.
        """
        self.topobathy_file = topobathy_file
        self.zb = rasterio.open(self.topobathy_file)

    def set_index_file(self, index_file: str | Path) -> None:
        """Set the index file and open it with rasterio.

        Parameters
        ----------
        index_file : str | Path
            Path to the cell-index COG file.
        """
        self.index_file = index_file
        self.indices = rasterio.open(self.index_file)

    def close(self) -> None:
        """Close the topobathy, index, and dataset file handles."""
        if self.zb is not None:
            self.zb.close()
        if self.indices is not None:
            self.indices.close()
        self.ds.close()

    def read(self, tiffile: str | Path) -> None:
        """Read a GeoTIFF file with pre-computed flood depth data.

        Parameters
        ----------
        tiffile : str | Path
            Path to the GeoTIFF file.
        """
        self.ds = xr.Dataset()
        self.ds["water_depth"] = rioxarray.open_rasterio(tiffile, masked=True).squeeze()

    def set_water_level(self, zs: float | np.ndarray) -> None:
        """Set the water level data used for flood depth computation.

        Parameters
        ----------
        zs : float | np.ndarray
            A scalar or 1-D array of water levels indexed by cell index.
        """
        self.zs = zs

    def make(
        self,
        max_pixel_size: float = 0.0,
        bbox: tuple[float, float, float, float] | None = None,
    ) -> xr.Dataset:
        """Compute flood depth from water levels, topobathy, and cell indices.

        Reads topobathy and index COGs at the appropriate overview level,
        computes ``h = zs[index] - zb``, and masks areas that are too
        shallow or outside elevation bounds.

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
            Dataset containing the computed flood depth array.
        """
        overview_level = 0

        if max_pixel_size > 0.0:
            overview_level = get_appropriate_overview_level(self.zb, max_pixel_size)

        if overview_level == 0:
            zb = rioxarray.open_rasterio(self.zb)
        else:
            zb = rioxarray.open_rasterio(self.zb, overview_level=overview_level)
        if "band" in zb.dims and zb.sizes["band"] == 1:
            zb = zb.squeeze(dim="band", drop=True)
        if bbox is not None:
            zb = zb.rio.clip_box(minx=bbox[0], miny=bbox[1], maxx=bbox[2], maxy=bbox[3])

        if overview_level == 0:
            indices = rioxarray.open_rasterio(self.indices)
        else:
            indices = rioxarray.open_rasterio(
                self.indices, overview_level=overview_level
            )
        if "band" in indices.dims and indices.sizes["band"] == 1:
            indices = indices.squeeze(dim="band", drop=True)
        if bbox is not None:
            indices = indices.rio.clip_box(
                minx=bbox[0], miny=bbox[1], maxx=bbox[2], maxy=bbox[3]
            )

        nan_val_indices = indices.attrs["_FillValue"]
        no_data_mask = indices == nan_val_indices
        indices = np.squeeze(indices.to_numpy()[:])
        indices[np.where(indices == nan_val_indices)] = 0

        if isinstance(self.zs, float):
            h = np.full(zb.shape, self.zs) - zb.to_numpy()[:]
        else:
            h = self.zs[indices] - zb.to_numpy()[:]
        h[no_data_mask] = np.nan
        h[h < self.hmin] = np.nan
        h[zb.to_numpy()[:] < self.zbmin] = np.nan
        h[zb.to_numpy()[:] > self.zbmax] = np.nan

        self.ds = xr.Dataset()
        self.ds[self.data_array_name] = xr.DataArray(
            h, dims=["y", "x"], coords={"y": zb.y, "x": zb.x}
        )
        self.ds[self.data_array_name] = self.ds[self.data_array_name].rio.write_crs(
            zb.rio.crs, inplace=True
        )

    def write(self, output_file: str | Path = "") -> None:
        """Write the flood map to a GeoTIFF or NetCDF file.

        Parameters
        ----------
        output_file : str | Path
            Output file path. Extension determines format: ``".tif"`` for
            COG GeoTIFF, ``".nc"`` for NetCDF.
        """
        if output_file.endswith(".nc"):
            self.ds.to_netcdf(output_file)

        elif output_file.endswith(".tif"):
            if self.cmap is not None:
                rgb_da = get_rgb_data_array(
                    self.ds[self.data_array_name],
                    color_values=self.color_values,
                    cmap=self.cmap,
                    cmin=self.cmin,
                    cmax=self.cmax,
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
        """Create a PNG map overlay of the flood map in EPSG:3857.

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
            logger.error(
                "Dataset is not initialized. Call make() or read() before map_overlay()."
            )
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
                discrete_colors=self.discrete_colors,
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

            if self.discrete_colors:
                self.legend = {}
                self.legend["title"] = "Flood Depth (m)"
                self.legend["contour"] = []

                if isinstance(self.color_values, str):
                    color_values = []
                    color_values.append(
                        {"color": "lightgreen", "lower_value": 0.1, "upper_value": 0.3}
                    )
                    color_values.append(
                        {"color": "yellow", "lower_value": 0.3, "upper_value": 1.0}
                    )
                    color_values.append(
                        {"color": "#FFA500", "lower_value": 1.0, "upper_value": 2.0}
                    )
                    color_values.append({"color": "red", "lower_value": 2.0})
                else:
                    color_values = self.color_values

                for cv in color_values:
                    legend_item = {}
                    legend_item["color"] = cv["color"]
                    if "upper_value" in cv and "lower_value" in cv:
                        legend_item["lower_value"] = cv["lower_value"]
                        legend_item["upper_value"] = cv["upper_value"]
                        legend_item["text"] = (
                            f"{cv['lower_value']}--{cv['upper_value']} m"
                        )
                    elif "upper_value" in cv:
                        legend_item["upper_value"] = cv["upper_value"]
                        legend_item["text"] = f"{cv['lower_value']}- m"
                    else:
                        legend_item["lower_value"] = cv["lower_value"]
                        legend_item["text"] = f"{cv['lower_value']}+ m"
                    self.legend["contour"].append(legend_item)

            else:
                self.legend = {}
                self.legend["title"] = "Flood Depth (m)"
                self.legend["cmin"] = self.cmin
                self.legend["cmax"] = self.cmax
                self.legend["cmap"] = self.cmap

            return True

        except Exception as e:
            logger.exception(e)
            return False

    def plot(
        self,
        pngfile: str,
        zoom: int | None = None,
        title: str = "Flood Depth (m)",
        color_values: list[dict] | None = None,
        cmap: str = "Blues",
        vmin: float = 0.0,
        vmax: float = 5.0,
        lon_lim: list[float] | None = None,
        lat_lim: list[float] | None = None,
        width: float = 10.0,
        background: str = "EsriWorldImagery",
    ) -> None:
        """Plot the flood map with a basemap and save to PNG.

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
            flood depth color scheme is used.
        cmap : str
            Matplotlib colormap for continuous coloring.
        vmin : float
            Minimum value for color mapping.
        vmax : float
            Maximum value for color mapping.
        lon_lim : list[float] | None
            Longitude limits ``[lon_min, lon_max]``.
        lat_lim : list[float] | None
            Latitude limits ``[lat_min, lat_max]``.
        width : float
            Figure width in inches.
        background : str
            Basemap provider: ``"osm"`` or ``"EsriWorldImagery"``.
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

        if color_values is None:
            discrete_colors = False
        else:
            discrete_colors = True
            if isinstance(color_values, str):
                color_values = []
                color_values.append(
                    {"color": "lightgreen", "lower_value": 0.1, "upper_value": 0.3}
                )
                color_values.append(
                    {"color": "yellow", "lower_value": 0.3, "upper_value": 1.0}
                )
                color_values.append(
                    {"color": "#FFA500", "lower_value": 1.0, "upper_value": 2.0}
                )
                color_values.append({"color": "red", "lower_value": 2.0})

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

            cmap = ListedColormap(colors)
            bounds = list(range(1, len(colors) + 2))

            norm = BoundaryNorm(bounds, cmap.N)

            classified.plot(ax=ax, cmap=cmap, norm=norm, add_colorbar=False)

            legend_elements = []
            for i, color_value in enumerate(color_values):
                legend_elements.append(
                    Patch(facecolor=color_value["color"], label=labels[i])
                )
            plt.legend(handles=legend_elements, title="Flood Depth", loc="lower right")

        else:
            da_3857.plot(
                ax=ax,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                add_colorbar=True,
                cbar_kwargs={"label": "Flood Depth (m)"},
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


def get_appropriate_overview_level(
    src: rasterio.io.DatasetReader, max_pixel_size: float
) -> int:
    """Determine the appropriate rasterio overview level for a target resolution.

    Parameters
    ----------
    src : rasterio.io.DatasetReader
        An open rasterio dataset.
    max_pixel_size : float
        Maximum desired pixel size in metres.

    Returns
    -------
    int
        The overview level index (0 = native resolution).
    """
    original_resolution = src.res
    if src.crs.is_geographic:
        original_resolution = (
            original_resolution[0] * 111000,
            original_resolution[1] * 111000,
        )

    overview_levels = src.overviews(1)

    if not overview_levels:
        return 0

    resolutions = [
        (original_resolution[0] * factor, original_resolution[1] * factor)
        for factor in overview_levels
    ]

    selected_overview = 0
    for i, (x_res, y_res) in enumerate(resolutions):
        if x_res <= max_pixel_size and y_res <= max_pixel_size:
            selected_overview = i
        else:
            break

    return selected_overview


def get_rgb_data_array(
    da: xr.DataArray,
    cmap: str,
    cmin: float | None = None,
    cmax: float | None = None,
    color_values: list[dict] | None = None,
    discrete_colors: bool = False,
) -> xr.DataArray:
    """Convert an xarray DataArray to an RGBA DataArray using a colormap.

    Supports both continuous colormaps and discrete color value ranges.

    Parameters
    ----------
    da : xr.DataArray
        Input 2-D data array.
    cmap : str
        Matplotlib colormap name for continuous coloring.
    cmin : float | None
        Minimum value for normalization. Defaults to data minimum.
    cmax : float | None
        Maximum value for normalization. Defaults to data maximum.
    color_values : list[dict] | None
        Discrete color definitions with ``lower_value``, ``upper_value``,
        and ``rgb`` keys.
    discrete_colors : bool
        If True and ``color_values`` is provided, use discrete coloring
        via named color strings.

    Returns
    -------
    xr.DataArray
        RGBA DataArray with shape ``(4, height, width)`` and dtype uint8.
    """
    ny, nx = da.shape
    if color_values is not None:
        zz = da.to_numpy()

    if discrete_colors:
        if isinstance(color_values, str):
            color_values = []
            color_values.append(
                {"color": "lightgreen", "lower_value": 0.1, "upper_value": 0.3}
            )
            color_values.append(
                {"color": "yellow", "lower_value": 0.3, "upper_value": 1.0}
            )
            color_values.append(
                {"color": "#FFA500", "lower_value": 1.0, "upper_value": 2.0}
            )
            color_values.append({"color": "red", "lower_value": 2.0})

        rgba = np.zeros((ny, nx, 4), "uint8")
        for color_value in color_values:
            lower = color_value.get("lower_value", -np.inf)
            upper = color_value.get("upper_value", np.inf)
            inr = np.logical_and(zz >= lower, zz < upper)
            valid = np.logical_and(inr, ~np.isnan(zz))

            if "rgb" in color_value:
                rgba[valid, 0] = color_value["rgb"][0]
                rgba[valid, 1] = color_value["rgb"][1]
                rgba[valid, 2] = color_value["rgb"][2]
            elif "color" in color_value:
                color_rgba = cm.colors.to_rgba(color_value["color"])
                rgba[valid, 0] = int(color_rgba[0] * 255)
                rgba[valid, 1] = int(color_rgba[1] * 255)
                rgba[valid, 2] = int(color_rgba[2] * 255)
            rgba[valid, 3] = 255

    else:
        if cmap is None:
            raise ValueError("Either color_values or cmap must be provided")

        if cmin is None:
            cmin = da.min()
        if cmax is None:
            cmax = da.max()

        if cmin == cmax:
            cmin = cmax - 1.0
            cmax = cmax + 1.0

        normed = (da - cmin) / (cmax - cmin)

        cmap_obj = plt.get_cmap(cmap)

        rgba = cmap_obj(normed)

        rgba = (rgba[:, :, :] * 255).astype("uint8")

    rgb_da = xr.DataArray(
        np.moveaxis(rgba, -1, 0),
        dims=("band", "y", "x"),
        coords={"band": [0, 1, 2, 3], "y": da.y, "x": da.x},
        attrs=da.attrs,
    )

    rgb_da.rio.write_crs(da.rio.crs, inplace=True)

    return rgb_da


def reproject_bbox(
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
    crs_src: str,
    crs_dst: str,
    buffer: float = 0.0,
) -> tuple[float, float, float, float]:
    """Reproject a bounding box between coordinate reference systems.

    Parameters
    ----------
    xmin : float
        Minimum x (or longitude).
    ymin : float
        Minimum y (or latitude).
    xmax : float
        Maximum x (or longitude).
    ymax : float
        Maximum y (or latitude).
    crs_src : str
        Source CRS string (e.g. ``"EPSG:4326"``).
    crs_dst : str
        Destination CRS string.
    buffer : float
        Fractional buffer to expand the bounding box before reprojection.

    Returns
    -------
    tuple[float, float, float, float]
        Reprojected bounding box ``(xmin, ymin, xmax, ymax)``.
    """
    transformer = Transformer.from_crs(crs_src, crs_dst, always_xy=True)

    dx = (xmax - xmin) * buffer
    dy = (ymax - ymin) * buffer
    xmin -= dx
    xmax += dx
    ymin -= dy
    ymax += dy

    x0, y0 = transformer.transform(xmin, ymin)
    x1, y1 = transformer.transform(xmax, ymin)
    x2, y2 = transformer.transform(xmax, ymax)
    x3, y3 = transformer.transform(xmin, ymax)

    xs = [x0, x1, x2, x3]
    ys = [y0, y1, y2, y3]

    return min(xs), min(ys), max(xs), max(ys)


def make_flood_map_tiles(
    valg: np.ndarray,
    index_path: str,
    png_path: str,
    topo_path: str,
    option: str = "deterministic",
    zoom_range: list[int] | None = None,
    color_values: list[dict] | None = None,
    caxis: list[float] | None = None,
    zbmax: float = -999.0,
    merge: bool = True,
    depth: float | None = None,
    quiet: bool = False,
) -> None:
    """Generate flood map PNG tiles from water level data and index/topo tiles.

    Parameters
    ----------
    valg : np.ndarray
        Water level values (1-D array indexed by cell index, or a list
        of CDF interpolators for probabilistic mode).
    index_path : str
        Directory containing index tile PNG files.
    png_path : str
        Output directory for the generated flood map tiles.
    topo_path : str
        Directory containing topobathy tile PNG files.
    option : str
        Tile generation mode: ``"deterministic"`` or ``"probabilistic"``.
    zoom_range : list[int] | None
        Two-element list ``[min_zoom, max_zoom]``. Auto-detected if None.
    color_values : list[dict] | None
        Discrete color definitions with ``lower_value``, ``upper_value``,
        and ``rgb`` keys.
    caxis : list[float] | None
        Color axis range ``[vmin, vmax]``. Auto-detected if None.
    zbmax : float
        Maximum bed level; flood in areas below this is suppressed.
    merge : bool
        Whether to merge new tiles with existing ones.
    depth : float | None
        Water depth offset for probabilistic mode.
    quiet : bool
        Whether to suppress progress output.
    """
    if isinstance(valg, list):
        pass
    else:
        valg = valg.transpose().flatten()

    if not caxis:
        caxis = []
        caxis.append(np.nanmin(valg))
        caxis.append(np.nanmax(valg))

    # First do highest zoom level, then derefine from there
    if not zoom_range:
        levs = fo.list_folders(os.path.join(index_path, "*"), basename=True)
        zoom_range = [999, -999]
        for lev in levs:
            zoom_range[0] = min(zoom_range[0], int(lev))
            zoom_range[1] = max(zoom_range[1], int(lev))

    izoom = zoom_range[1]

    if not quiet:
        logger.info(f"Processing zoom level {izoom}")

    index_zoom_path = os.path.join(index_path, str(izoom))

    png_zoom_path = os.path.join(png_path, str(izoom))
    fo.mkdir(png_zoom_path)

    for ifolder in fo.list_folders(os.path.join(index_zoom_path, "*")):
        path_okay = False
        ifolder = os.path.basename(ifolder)
        index_zoom_path_i = os.path.join(index_zoom_path, ifolder)
        png_zoom_path_i = os.path.join(png_zoom_path, ifolder)

        for jfile in fo.list_files(os.path.join(index_zoom_path_i, "*.png")):
            jfile = os.path.basename(jfile)
            j = int(jfile[:-4])

            index_file = os.path.join(index_zoom_path_i, jfile)
            png_file = os.path.join(png_zoom_path_i, f"{j}.png")

            ind = png2int(index_file, -1)
            ind = ind.flatten()

            if option == "probabilistic":
                bathy_file = os.path.join(topo_path, str(izoom), ifolder, f"{j}.png")
                if not os.path.exists(bathy_file):
                    continue
                zb = png2elevation(bathy_file).flatten()
                zs = zb + depth

                valt = valg[ind](zs)
                valt[ind < 0] = np.nan

            else:
                bathy_file = os.path.join(topo_path, str(izoom), ifolder, f"{j}.png")
                if not os.path.exists(bathy_file):
                    continue
                zb = png2elevation(bathy_file).flatten()

                noval = np.where(ind < 0)
                ind[ind < 0] = 0
                valt = valg[ind]

                valt = valt - zb
                valt[valt < 0.10] = np.nan
                valt[zb < zbmax] = np.nan
                valt[noval] = np.nan

            if color_values:
                rgb = np.zeros((256 * 256, 4), "uint8")

                for color_value in color_values:
                    inr = np.logical_and(
                        valt >= color_value["lower_value"],
                        valt < color_value["upper_value"],
                    )
                    rgb[inr, 0] = color_value["rgb"][0]
                    rgb[inr, 1] = color_value["rgb"][1]
                    rgb[inr, 2] = color_value["rgb"][2]
                    rgb[inr, 3] = 255

                rgb = rgb.reshape([256, 256, 4])
                if not np.any(rgb > 0):
                    continue
                im = Image.fromarray(rgb)

            else:
                valt = valt.reshape([256, 256])
                valt = (valt - caxis[0]) / (caxis[1] - caxis[0])
                valt[valt < 0.0] = 0.0
                valt[valt > 1.0] = 1.0
                im = Image.fromarray(cm.jet(valt, bytes=True))

            if not path_okay:
                if not os.path.exists(png_zoom_path_i):
                    fo.mkdir(png_zoom_path_i)
                    path_okay = True

            if os.path.exists(png_file):
                if merge:
                    im0 = Image.open(png_file)
                    rgb = np.array(im)
                    rgb0 = np.array(im0)
                    isum = np.sum(rgb, axis=2)
                    rgb[isum == 0, :] = rgb0[isum == 0, :]
                    im = Image.fromarray(rgb)

            im.save(png_file)

    # Now make tiles for lower level by merging

    for izoom in range(zoom_range[1] - 1, zoom_range[0] - 1, -1):
        if not quiet:
            logger.info(f"Processing zoom level {izoom}")

        index_zoom_path = os.path.join(index_path, str(izoom))

        if not os.path.exists(index_zoom_path):
            continue

        png_zoom_path = os.path.join(png_path, str(izoom))
        png_zoom_path_p1 = os.path.join(png_path, str(izoom + 1))
        fo.mkdir(png_zoom_path)

        for ifolder in fo.list_folders(os.path.join(index_zoom_path, "*")):
            path_okay = False
            ifolder = os.path.basename(ifolder)
            i = int(ifolder)
            index_zoom_path_i = os.path.join(index_zoom_path, ifolder)
            png_zoom_path_i = os.path.join(png_zoom_path, ifolder)

            for jfile in fo.list_files(os.path.join(index_zoom_path_i, "*.png")):
                jfile = os.path.basename(jfile)
                j = int(jfile[:-4])

                png_file = os.path.join(png_zoom_path_i, f"{j}.png")

                rgb = np.zeros((256, 256, 4), "uint8")

                i0 = i * 2
                i1 = i * 2 + 1
                j0 = j * 2 + 1
                j1 = j * 2

                tile_name_00 = os.path.join(png_zoom_path_p1, str(i0), f"{j0}.png")
                tile_name_10 = os.path.join(png_zoom_path_p1, str(i0), f"{j1}.png")
                tile_name_01 = os.path.join(png_zoom_path_p1, str(i1), f"{j0}.png")
                tile_name_11 = os.path.join(png_zoom_path_p1, str(i1), f"{j1}.png")

                okay = False

                # Lower-left
                if os.path.exists(tile_name_00):
                    okay = True
                    rgb0 = np.array(Image.open(tile_name_00))
                    rgb[128:256, 0:128, :] = rgb0[0:255:2, 0:255:2, :]
                # Upper-left
                if os.path.exists(tile_name_10):
                    okay = True
                    rgb0 = np.array(Image.open(tile_name_10))
                    rgb[0:128, 0:128, :] = rgb0[0:255:2, 0:255:2, :]
                # Lower-right
                if os.path.exists(tile_name_01):
                    okay = True
                    rgb0 = np.array(Image.open(tile_name_01))
                    rgb[128:256, 128:256, :] = rgb0[0:255:2, 0:255:2, :]
                # Upper-right
                if os.path.exists(tile_name_11):
                    okay = True
                    rgb0 = np.array(Image.open(tile_name_11))
                    rgb[0:128, 128:256, :] = rgb0[0:255:2, 0:255:2, :]

                if okay:
                    im = Image.fromarray(rgb)

                    if not path_okay:
                        if not os.path.exists(png_zoom_path_i):
                            fo.mkdir(png_zoom_path_i)
                            path_okay = True

                    if os.path.exists(png_file):
                        if merge:
                            im0 = Image.open(png_file)
                            rgb = np.array(im)
                            rgb0 = np.array(im0)
                            isum = np.sum(rgb, axis=2)
                            rgb[isum == 0, :] = rgb0[isum == 0, :]
                            im = Image.fromarray(rgb)

                    im.save(png_file)


def make_flood_map_overlay_v2(
    valg: np.ndarray,
    index_path: str,
    topo_path: str,
    zmax_minus_zmin: np.ndarray | None = None,
    mean_depth: np.ndarray | None = None,
    npixels: list[int] = [1200, 800],
    hmin: float = 0.10,
    dzdx_mild: float = 0.01,
    lon_range: list[float] | None = None,
    lat_range: list[float] | None = None,
    option: str = "deterministic",
    color_values: list[dict] | None = None,
    caxis: list[float] | None = None,
    zbmax: float = -999.0,
    merge: bool = True,
    depth: float | None = None,
    quiet: bool = False,
    file_name: str | None = None,
) -> tuple[list[float], list[float], list[float]] | tuple[None, None]:
    """Generate a single flood map overlay PNG from tiles at an auto-selected zoom.

    Parameters
    ----------
    valg : np.ndarray
        Water level values (1-D array or DataArray).
    index_path : str
        Directory containing index tile PNG files.
    topo_path : str
        Directory containing topobathy tile PNG files.
    zmax_minus_zmin : np.ndarray | None
        Per-cell elevation range for slope filtering.
    mean_depth : np.ndarray | None
        Per-cell mean water depth (volume / area).
    npixels : list[int]
        Target output size ``[width, height]`` in pixels.
    hmin : float
        Minimum water depth threshold.
    dzdx_mild : float
        Slope threshold below which mean_depth overrides pixel depth.
    lon_range : list[float] | None
        Longitude range ``[lon_min, lon_max]``.
    lat_range : list[float] | None
        Latitude range ``[lat_min, lat_max]``.
    option : str
        Mode: ``"deterministic"`` or ``"probabilistic"``.
    color_values : list[dict] | None
        Discrete color definitions.
    caxis : list[float] | None
        Color axis range. Auto-detected if None.
    zbmax : float
        Maximum bed level for flood masking.
    merge : bool
        Whether to merge with existing tiles.
    depth : float | None
        Water depth offset for probabilistic mode.
    quiet : bool
        Whether to suppress progress output.
    file_name : str | None
        Output PNG file path.

    Returns
    -------
    tuple[list[float], list[float], list[float]] | tuple[None, None]
        ``([lon_min, lon_max], [lat_min, lat_max], caxis)`` on success,
        or ``(None, None)`` on failure.
    """
    try:
        if isinstance(valg, list):
            logger.info("valg is a list!")
        elif isinstance(valg, xr.DataArray):
            valg = valg.to_numpy()
            if mean_depth is not None:
                mean_depth = mean_depth.to_numpy()
            if zmax_minus_zmin is not None:
                zmax_minus_zmin = zmax_minus_zmin.to_numpy()
        else:
            valg = valg.transpose().flatten()
            if mean_depth is not None:
                mean_depth = mean_depth.transpose().flatten()
            if zmax_minus_zmin is not None:
                zmax_minus_zmin = zmax_minus_zmin.transpose().flatten()

        if mean_depth is not None and zmax_minus_zmin is not None:
            mean_depth[(zmax_minus_zmin > dzdx_mild)] = np.nan

        max_zoom = 0
        levs = fo.list_folders(os.path.join(index_path, "*"), basename=True)
        for lev in levs:
            max_zoom = max(max_zoom, int(lev))

        for izoom in range(max_zoom + 1):
            ix0, it0 = deg2num(lat_range[1], lon_range[0], izoom)
            ix1, it1 = deg2num(lat_range[0], lon_range[1], izoom)
            if (ix1 - ix0 + 1) * 256 > npixels[0] and (it1 - it0 + 1) * 256 > npixels[
                1
            ]:
                break

        index_zoom_path = os.path.join(index_path, str(izoom))

        nx = (ix1 - ix0 + 1) * 256
        ny = (it1 - it0 + 1) * 256
        zz = np.empty((ny, nx))
        zz[:] = np.nan

        if not quiet:
            logger.info(f"Processing zoom level {izoom}")

        index_zoom_path = os.path.join(index_path, str(izoom))

        for i in range(ix0, ix1 + 1):
            ifolder = str(i)
            index_zoom_path_i = os.path.join(index_zoom_path, ifolder)

            for j in range(it0, it1 + 1):
                index_file = os.path.join(index_zoom_path_i, f"{j}.png")

                if not os.path.exists(index_file):
                    continue

                ind = png2int(index_file, -1)

                if option == "probabilistic":
                    bathy_file = os.path.join(
                        topo_path, str(izoom), ifolder, f"{j}.png"
                    )

                    if not os.path.exists(bathy_file):
                        continue

                    zb = np.fromfile(bathy_file, dtype="f4")
                    zs = zb + depth

                    valt = valg[ind](zs)
                    valt[ind < 0] = np.nan

                else:
                    bathy_file = os.path.join(
                        topo_path, str(izoom), ifolder, f"{j}.png"
                    )
                    if not os.path.exists(bathy_file):
                        continue

                    zb = png2elevation(bathy_file)

                    valt = valg[ind]

                    valt = valt - zb

                    if mean_depth is not None:
                        mean_depth_p = mean_depth[ind]
                        valt[~np.isnan(mean_depth_p)] = mean_depth_p[
                            ~np.isnan(mean_depth_p)
                        ]

                    valt[valt < hmin] = np.nan
                    valt[zb < zbmax] = np.nan

                ii0 = (i - ix0) * 256
                ii1 = ii0 + 256
                jj0 = (j - it0) * 256
                jj1 = jj0 + 256
                zz[jj0:jj1, ii0:ii1] = valt

        if color_values:
            zz = zz.flatten()
            rgb = np.zeros((ny * nx, 4), "uint8")
            for color_value in color_values:
                inr = np.logical_and(
                    zz >= color_value["lower_value"], zz < color_value["upper_value"]
                )
                rgb[inr, 0] = color_value["rgb"][0]
                rgb[inr, 1] = color_value["rgb"][1]
                rgb[inr, 2] = color_value["rgb"][2]
                rgb[inr, 3] = 255
            im = Image.fromarray(rgb.reshape([ny, nx, 4]))

        else:
            if not caxis:
                caxis = []
                caxis.append(np.nanmin(valg))
                caxis.append(np.nanmax(valg))

            zz = (zz - caxis[0]) / (caxis[1] - caxis[0])
            zz[zz < 0.0] = 0.0
            zz[zz > 1.0] = 1.0
            im = Image.fromarray(cm.jet(zz, bytes=True))

        if file_name:
            im.save(file_name)

        lat1, lon0 = num2deg(ix0, it0, izoom)
        lat0, lon1 = num2deg(ix1 + 1, it1 + 1, izoom)

        return [lon0, lon1], [lat0, lat1], caxis

    except Exception as e:
        logger.exception(e)
        return None, None
