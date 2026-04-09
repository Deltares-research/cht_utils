"""Generate probability maps and merge ensemble NetCDF output.

Supports both regular grids and quadtree (unstructured mesh) datasets.
"""

from typing import List, Optional

import numpy as np
import xarray as xr

from cht_utils.fileops import delete_file


def prob_floodmaps(
    file_list: List[str],
    variables: List[str],
    prcs: List[float],
    delete: bool = False,
    output_file_name: Optional[str] = None,
) -> None:
    """Generate probability flood maps from ensemble NetCDF files.

    Parameters
    ----------
    file_list : List[str]
        Paths to ensemble member NetCDF files.
    variables : List[str]
        Variable names to process.
    prcs : List[float]
        Percentiles to compute (0-100 scale).
    delete : bool
        Delete input files after reading.
    output_file_name : str or None
        Output NetCDF file path.
    """
    ds_concat = []
    for file in file_list:
        dsin = xr.open_dataset(file)
        ds_concat.append(dsin[variables])
        if delete:
            delete_file(file)

    combined_ds = xr.concat(ds_concat, dim="ensemble")

    out_qq = {}
    for v in variables:
        combined_ds[v] = xr.where(np.isnan(combined_ds[v]), 0, combined_ds[v])
        out_qq[v] = np.percentile(combined_ds[v], prcs, axis=0)
        out_qq[v] = np.where(out_qq[v] == 0, np.nan, out_qq[v])

    dsin = xr.open_dataset(file_list[0])
    for i, p in enumerate(prcs):
        for v in variables:
            dsin[f"{v}_{p}"] = xr.DataArray(out_qq[v][i], dims=dsin[v].dims)
            dsin[f"{v}_{p}"].attrs["long_name"] = f"{v}_{p}"

    try:
        delete_file(output_file_name)
    except Exception:
        pass
    dsin.to_netcdf(path=output_file_name)
    dsin.close()


def merge_nc_his(
    file_list: List[str],
    variables: List[str],
    prcs: Optional[List[float]] = None,
    delete: bool = False,
    output_file_name: Optional[str] = None,
) -> None:
    """Merge ensemble history files and compute quantiles.

    Parameters
    ----------
    file_list : List[str]
        Paths to ensemble member NetCDF files.
    variables : List[str]
        Variable names to merge.
    prcs : List[float] or None
        Quantile probabilities (0-1 scale). Defaults to ``[0.05, 0.5, 0.95]``.
    delete : bool
        Delete input files after reading.
    output_file_name : str or None
        Output NetCDF file path.
    """
    if prcs is None:
        prcs = [0.05, 0.5, 0.95]
    if len(file_list) == 0:
        print("his-file list is empty")
        return

    nens = len(file_list)
    ds = xr.open_dataset(file_list[0])
    if "runtime" in ds.dims:
        ds = ds.drop_dims("runtime")
    if "structures" in ds.dims:
        ds = ds.drop_dims("structures")

    dimensions = list(ds.dims.keys())
    dimensions = ["time"] + [d for d in dimensions if d != "time"]
    new_dimensions = dimensions + ["ensemble"]
    ens = range(nens)

    for v in variables:
        ds[v] = xr.DataArray(
            data=np.zeros(
                tuple(ds.dims[d] if d in ds.dims else nens for d in new_dimensions)
            ),
            dims=new_dimensions,
            coords={d: ds[d] if d in ds.dims else ens for d in new_dimensions},
        )

    for iens, file in enumerate(file_list):
        dsin = xr.open_dataset(file)
        if "runtime" in dsin.dims:
            dsin = dsin.drop_dims("runtime")
        if "structures" in dsin.dims:
            dsin = dsin.drop_dims("structures")
        for v in variables:
            ds[v][:, :, iens] = dsin[v].transpose("time", ...)

    for v in variables:
        arr = ds[v].fillna(-999.0).quantile(prcs, dim="ensemble")
        for ip, p in enumerate(prcs):
            ds[f"{v}_{round(p * 100)}"] = arr[ip, :, :]

    try:
        delete_file(output_file_name)
    except Exception:
        pass
    ds.to_netcdf(path=output_file_name)
    ds.close()


def merge_nc_map(
    file_list: List[str],
    variables: List[str],
    prcs: Optional[List[float]] = None,
    delete: bool = False,
    output_file_name: Optional[str] = None,
) -> None:
    """Merge ensemble map files and compute quantiles.

    Supports both regular grids and quadtree (unstructured) meshes.
    Uses sorting instead of ``xr.quantile`` for performance.

    Parameters
    ----------
    file_list : List[str]
        Paths to ensemble member NetCDF files.
    variables : List[str]
        Variable names to merge.
    prcs : List[float] or None
        Quantile probabilities (0-1 scale). Defaults to ``[0.9]``.
    delete : bool
        Delete input files after reading.
    output_file_name : str or None
        Output NetCDF file path.
    """
    if prcs is None:
        prcs = [0.9]

    nens = len(file_list)
    ds = xr.open_dataset(file_list[0])

    if "mesh2d_face_nodes" in ds:
        _merge_quadtree(ds, file_list, variables, prcs, nens)
    else:
        _merge_regular(ds, file_list, variables, prcs, nens)

    # Remove original variables
    to_drop = ["zs", "zsmax", "cumprcp", "cuminf", "qinf", "hm0max"]
    for v in to_drop:
        if v in ds:
            ds = ds.drop(v)

    try:
        delete_file(output_file_name)
    except Exception:
        pass
    ds.to_netcdf(path=output_file_name)
    ds.close()


def _merge_quadtree(
    ds: xr.Dataset,
    file_list: List[str],
    variables: List[str],
    prcs: List[float],
    nens: int,
) -> None:
    """Merge quadtree (unstructured) map data."""
    npoints = ds["nmesh2d_face"].size
    for v in variables:
        attrs = ds[v].attrs
        timedim = ds[v].dims[0]
        ntimes = ds.sizes[timedim]

        arr = xr.DataArray(
            data=np.zeros((ntimes, npoints, nens)),
            dims=(timedim, "nmesh2d_face", "ensemble"),
            coords={timedim: ds[timedim], "ensemble": range(nens)},
        )
        ds[v] = arr

        for iens, file in enumerate(file_list):
            dsin = xr.open_dataset(file)
            ds[v][:, :, iens] = dsin[v]
            dsin.close()

        sorted_arr = np.sort(ds[v].values, axis=2)
        for ip, p in enumerate(prcs):
            indx = min(max(int(np.ceil(p * nens)) - 1, 0), nens - 1)
            da = xr.DataArray(
                data=sorted_arr[:, :, indx],
                dims=(timedim, "nmesh2d_face"),
                coords={timedim: ds[timedim]},
            )
            da.attrs = attrs
            ds[f"{v}_{round(p * 100)}"] = da


def _merge_regular(
    ds: xr.Dataset,
    file_list: List[str],
    variables: List[str],
    prcs: List[float],
    nens: int,
) -> None:
    """Merge regular grid map data."""
    for v in variables:
        attrs = ds[v].attrs
        timedim = ds[v].dims[0]
        ntimes = ds.sizes[timedim]

        arr = xr.DataArray(
            data=np.zeros((ntimes, ds.x.shape[0], ds.x.shape[1], nens)),
            dims=(timedim, "n", "m", "ensemble"),
            coords={
                timedim: ds[timedim],
                "x": ds.x,
                "y": ds.y,
                "ensemble": range(nens),
            },
        )
        ds[v] = arr

        for iens, file in enumerate(file_list):
            dsin = xr.open_dataset(file)
            ds[v][:, :, :, iens] = dsin[v]
            dsin.close()

        sorted_arr = np.sort(ds[v].values, axis=3)
        for ip, p in enumerate(prcs):
            indx = min(max(int(np.ceil(p * nens)) - 1, 0), nens - 1)
            da = xr.DataArray(
                data=sorted_arr[:, :, :, indx],
                dims=(timedim, "n", "m"),
                coords={timedim: ds[timedim], "x": ds.x, "y": ds.y},
            )
            da.attrs = attrs
            ds[f"{v}_{round(p * 100)}"] = da
