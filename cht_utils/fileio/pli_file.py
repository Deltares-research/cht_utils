"""Read and write Deltares PLI/POL files (polyline and polygon data).

PLI files store polylines; POL files store polygons. Both use the Tekal
block format internally.
"""

from typing import Optional

import geopandas as gpd
import pandas as pd
import shapely.geometry

import cht_utils.fileio.tekal as tek


class Polyline:
    """Simple polyline wrapper with x/y coordinate arrays.

    Parameters
    ----------
    x : array-like
        X-coordinates.
    y : array-like
        Y-coordinates.
    """

    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y


def read_pli_file(file_name: str) -> list:
    """Read all polylines from a PLI file.

    Parameters
    ----------
    file_name : str
        Path to the PLI file.

    Returns
    -------
    list of Polyline
        List of polyline objects.
    """
    polylines = []
    D = tek.tekal(file_name)
    D.info()
    for j in range(len(D.blocks)):
        m = D.read(j)
        polylines.append(Polyline(m[0, :, 0], m[1, :, 0]))
    return polylines


def pli2gdf(
    file_name: str,
    crs: Optional[int] = None,
    name_string: str = "name",
) -> gpd.GeoDataFrame:
    """Convert a PLI file to a GeoDataFrame of LineStrings.

    Parameters
    ----------
    file_name : str
        Path to the PLI file.
    crs : int or None
        Coordinate reference system (EPSG code).
    name_string : str
        Column name for the polyline name attribute.

    Returns
    -------
    gpd.GeoDataFrame
    """
    D = tek.tekal(file_name)
    D.info()
    gdf_list = []
    for j in range(len(D.blocks)):
        name = D.blocks[j].name
        if isinstance(name, bytes):
            name = name.decode("utf-8")
        m = D.read(j)
        x, y = m[0, :, 0], m[1, :, 0]
        line = shapely.geometry.LineString(list(zip(x, y)))
        gdf_list.append({name_string: name, "geometry": line})
    return gpd.GeoDataFrame(gdf_list, crs=crs)


def pli2geojson(
    file_name_in: str,
    file_name_out: str,
    crs: Optional[int] = None,
) -> None:
    """Convert a PLI file to GeoJSON.

    Parameters
    ----------
    file_name_in : str
        Input PLI file path.
    file_name_out : str
        Output GeoJSON file path.
    crs : int or None
        Coordinate reference system (EPSG code).
    """
    gdf = pli2gdf(file_name_in, crs=crs)
    gdf.to_file(file_name_out, driver="GeoJSON")


def pol2gdf(
    file_name: str,
    crs: Optional[int] = None,
    header: bool = True,
) -> gpd.GeoDataFrame:
    """Convert a POL file to a GeoDataFrame of Polygons.

    Parameters
    ----------
    file_name : str
        Path to the POL file.
    crs : int or None
        Coordinate reference system (EPSG code).
    header : bool
        If ``False``, read as plain whitespace-separated x/y columns.

    Returns
    -------
    gpd.GeoDataFrame
    """
    gdf_list = []
    if not header:
        df = pd.read_csv(
            file_name,
            index_col=False,
            header=None,
            delim_whitespace=True,
            names=["x", "y"],
        )
        poly = shapely.geometry.Polygon(list(zip(df.x.values, df.y.values)))
        gdf_list.append({"geometry": poly})
    else:
        D = tek.tekal(file_name)
        D.info()
        for j in range(len(D.blocks)):
            m = D.read(j)
            x, y = m[0, :, 0], m[1, :, 0]
            poly = shapely.geometry.Polygon(list(zip(x, y)))
            gdf_list.append({"geometry": poly})
    return gpd.GeoDataFrame(gdf_list, crs=crs)


def pol2geojson(
    file_name_in: str,
    file_name_out: str,
    crs: Optional[int] = None,
) -> None:
    """Convert a POL file to GeoJSON.

    Parameters
    ----------
    file_name_in : str
        Input POL file path.
    file_name_out : str
        Output GeoJSON file path.
    crs : int or None
        Coordinate reference system (EPSG code).
    """
    gdf = pol2gdf(file_name_in, crs=crs)
    gdf.to_file(file_name_out, driver="GeoJSON")


def gdf2pli(
    gdf: gpd.GeoDataFrame,
    file_name: str,
    header: bool = True,
    name_string: str = "name",
    add_point_name: bool = False,
) -> None:
    """Write a GeoDataFrame of LineStrings to a PLI file.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame with LineString geometries.
    file_name : str
        Output PLI file path.
    header : bool
        Write Tekal block headers.
    name_string : str
        Column containing the polyline name.
    add_point_name : bool
        Append vertex names (``NAME_0001``, etc.) to each coordinate line.
    """
    fmt = "{:12.6f}" if gdf.crs.is_geographic else "{:12.1f}"

    with open(file_name, "w") as fid:
        for index, row in gdf.iterrows():
            coords = list(row["geometry"].coords)
            nrp = len(coords)
            if header:
                hdr = row.get(name_string, f"BL{index + 1:04d}")
                fid.write(f"{hdr}\n{nrp} 2\n")
            for ip, (x, y) in enumerate(coords):
                line = fmt.format(x) + fmt.format(y)
                if add_point_name:
                    line += f" {hdr}_{ip + 1:04d}"
                fid.write(line + "\n")


def gdf2pol(
    gdf: gpd.GeoDataFrame,
    file_name: str,
    header: bool = True,
) -> None:
    """Write a GeoDataFrame of Polygons to a POL file.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame with Polygon geometries.
    file_name : str
        Output POL file path.
    header : bool
        Write Tekal block headers.
    """
    fmt = "{:12.6f}" if gdf.crs.is_geographic else "{:12.1f}"

    with open(file_name, "w") as fid:
        for index, row in gdf.iterrows():
            coords = list(row["geometry"].exterior.coords)
            nrp = len(coords)
            if header:
                fid.write(f"BL{index + 1:04d}\n{nrp} 2\n")
            for x, y in coords:
                fid.write(fmt.format(x) + fmt.format(y) + "\n")
