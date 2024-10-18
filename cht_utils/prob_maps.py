# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 13:57:00 2023

@author: roelvink
"""

import xarray as xr
import numpy as np
import cht_utils.fileops as fo
import numpy as np

def prob_floodmaps(file_list, variables, prcs, delete=False, output_file_name=None):
    
    out_qq = {}
    ds_concat = []

    for file in file_list:
        dsin = xr.open_dataset(file)
        ds_concat.append(dsin[variables])
        if delete:
            fo.delete_file(file)

    combined_ds = xr.concat(ds_concat, dim='ensemble')

    for v in variables:
        combined_ds[v] = xr.where(np.isnan(combined_ds[v]), 0, combined_ds[v])
        out_qq[v] = np.percentile(combined_ds[v], prcs, axis=0)
        out_qq[v] = np.where(out_qq[v] == 0, np.nan, out_qq[v])

    dsin = xr.open_dataset(file_list[0])

    for i, v in enumerate(prcs):
        for vv in variables:
            dsin[vv + "_" + str(v)] = xr.DataArray(out_qq[vv][i], dims=dsin[vv].dims)
            dsin[vv + "_" + str(v)].attrs['long_name'] = vv + "_" + str(v)

    try:
        fo.delete_file(output_file_name)
    except:
        pass

    dsin.to_netcdf(path=output_file_name)
    dsin.close()

def merge_nc_his(file_list, variables, prcs=None, delete=False, output_file_name=None):

    if prcs is None:
        prcs = [0.05, 0.5, 0.95]
    
    if len(file_list)==0:
        print('his-file list is empty')
        return
    nens = len(file_list)

    # Read first file    
    ds = xr.open_dataset(file_list[0])
    if 'runtime' in ds.dims:
        ds = ds.drop_dims('runtime') #hurrywave
    if "structures" in ds.dims:
        ds = ds.drop_dims('structures') #sfincs

    # Read dimensions from first file
    dimensions = list(ds.dims.keys())

    # Ensure 'time' comes first in the list of dimensions
    dimensions = ['time'] + [dim for dim in dimensions if dim != 'time']

    # Create ensemble dimension
    new_dimensions = list(dimensions + ['ensemble'])
    ens = range(nens)

    # Make new data array filled with zeros with dimensions time, stations and ens
    for v in variables:
        ds[v] = xr.DataArray(
                data=np.zeros(tuple(ds.dims[dim] if dim in ds.dims else nens for dim in new_dimensions)),
                dims=new_dimensions,
                coords={dim: ds[dim] if dim in ds.dims else ens for dim in new_dimensions})

    for iens, file in enumerate(file_list):
        dsin = xr.open_dataset(file)

        if 'runtime' in dsin.dims:
            dsin = dsin.drop_dims('runtime') #hurrywave
        if "structures" in dsin.dims:
            dsin = dsin.drop_dims('structures') #sfincs

        for v in variables:
            ds[v][:,:,iens] = dsin[v].transpose('time', ...)

    for v in variables:
        ds[v].fillna(-999.0) 
        arr = ds[v].fillna(-999.0).quantile(prcs, dim="ensemble")
        for ip, p in enumerate(prcs):
            ds[v + "_" + str(round(p*100))] = arr[ip,:,:] 
        
    try:
        fo.delete_file(output_file_name)
    except:
        pass
    ds.to_netcdf(path= output_file_name)
    ds.close()

def merge_nc_map(file_list, variables, prcs=None, delete=False, output_file_name=None):

    if prcs is None:
        prcs = [0.9]

    nens = len(file_list)

    # Read first file (this is the xarray dataset that will be used to store the data merged data)

    ds = xr.open_dataset(file_list[0])

    if "mesh2d_face_nodes" in ds:

        # Quadtree

        npoints = ds["nmesh2d_face"].size

        for v in variables:

            # Get the attributes of the variable
            attrs = ds[v].attrs

            # Create a new dataarray with dimensions time, x, y, ens
            arr = xr.DataArray(data=np.zeros((len(ds.timemax), npoints, nens)),
                               dims=('timemax', 'nmesh2d_face', 'ensemble'),
                               coords={'timemax': ds.timemax, 'ensemble': range(nens)})
            ds[v] = arr

            for iens, file in enumerate(file_list):
                dsin = xr.open_dataset(file)
                ds[v][:,:,iens] = dsin[v]
                dsin.close()
            
            # arr = ds[v].quantile(prcs, dim="ensemble", skipna=True)
            # for ip, p in enumerate(prcs):
            #     ds[v + "_" + str(round(p*100))] = arr[ip,:,:]

            # Quantile method takes insanely long ...
            # Try simple sorting and indexing instead (round up to be conservative)
            arr = np.sort(ds[v], axis=2)
            for ip, p in enumerate(prcs):
                indx = min(max(int(np.ceil(p * nens)) - 1, 0), nens - 1)
                da = xr.DataArray(data=arr[:,:,indx],
                                  dims=('timemax', 'nmesh2d_face'),
                                  coords={'timemax': ds.timemax})
                da.attrs = attrs
                ds[v + "_" + str(round(p*100))] = da

    else:

        # Regular grid

        # Loop over variables to merge  
        for v in variables:

            # Get the attributes of the variable
            attrs = ds[v].attrs

            # Create a dataarray with dimensions time, x, y, ens
            arr = xr.DataArray(data=np.zeros((len(ds.timemax), np.shape(ds.x)[0], np.shape(ds.x)[1], nens)),
                               dims=('timemax', 'n', 'm', 'ensemble'),
                               coords={'timemax': ds.timemax, 'x': ds.x, 'y': ds.y, 'ensemble': range(nens)})
            ds[v] = arr

            for iens, file in enumerate(file_list):
                dsin = xr.open_dataset(file)
                ds[v][:,:,:,iens] = dsin[v]
                dsin.close()

            # arr = ds[v].quantile(prcs, dim="ensemble", skipna=True)
            # for ip, p in enumerate(prcs):
            #     ds[v + "_" + str(round(p*100))] = arr[ip,:,:,:]

            # Quantile method takes insanely long ...
            # Try simple sorting and indexing instead (round up to be conservative)
            arr = np.sort(ds[v], axis=3)
            for ip, p in enumerate(prcs):
                indx = min(max(int(np.ceil(p * nens)) - 1, 0), nens - 1)
                da = xr.DataArray(data=arr[:,:,indx],
                                  dims=('timemax', 'n', 'm'),
                                  coords={'timemax': ds.timemax, 'x': ds.x, 'y': ds.y})
                da.attrs = attrs
                ds[v + "_" + str(round(p*100))] = da

    # Remove the original variables
    to_drop = ["zs", "zsmax", "cumprcp", "cuminf", "qinf", "hm0max"]
    for v in to_drop:
        if v in ds:
            ds = ds.drop(v)
        
    try:
        fo.delete_file(output_file_name)
    except:
        pass

    ds.to_netcdf(path=output_file_name)
    ds.close()


def merge_nc_map_1by1(file_list, variables, prcs=None, delete=False, output_file_name=None, nbins=10):

    if prcs is None:
        prcs = [0.9]

    nens = len(file_list)

    # Read first file    
    ds = xr.open_dataset(file_list[0])

    # Read time and stations from first file
    # time = ds.timemax
    # x = ds.x
    # y = ds.y
    # ens = range(nens)

    d0 = np.zeros((len(ds.timemax), np.shape(ds.x)[0], np.shape(ds.x)[1]))
    i0 = np.zeros((len(ds.timemax), np.shape(ds.x)[0], np.shape(ds.x)[1]), nbins).astype(int)
    
    # Loop through variables to get min and max at all time steps
    for v in variables:
        minstr = v + "_min"
        maxstr = v + "_max"
        ds[minstr] = xr.DataArray(data=d0-999.0,
                                  dims=('timemax', 'n', 'm'),
                                  coords={'timemax': ds.timemax, 'x': ds.x, 'y': ds.y})
        ds[maxstr] = xr.DataArray(data=d0+999.0,
                                  dims=('timemax', 'n', 'm'),
                                  coords={'timemax': ds.timemax, 'x': ds.x, 'y': ds.y})

    # Now loop through all ensemble members to get min and max
    for iens, file in enumerate(file_list):
        dsin = xr.open_dataset(file)
        for v in variables:
            minstr = v + "_min"
            maxstr = v + "_max"
            ds[minstr] = np.min(dsin[v], ds[minstr])
            ds[maxstr] = np.min(dsin[v], ds[maxstr])

    # Loop through variables to get bins        
    for v in variables:
        minstr = v + "_min"
        maxstr = v + "_max"
        ds[v + "_dbin"] = (ds[maxstr] - ds[minstr]) / (nbins - 1)
        ds[v + "_ibin"] = xr.DataArray(data=i0,
                                       dims=('timemax', 'n', 'm', 'bin'),
                                       coords={'timemax': ds.timemax, 'x': ds.x, 'y': ds.y, 'ibin': range(nbins)})
     
    # Now again loop through all ensemble members to determine number of occurences in bins
    for iens, file in enumerate(file_list):
        dsin = xr.open_dataset(file)
        for v in variables:
            minstr = v + "_min"
            maxstr = v + "_max"
            ibin   = (dsin[v] - ds[v + "_min"]) / ds[v + "_dbin"]
            ds[v + "_ibin"][ibin] += 1

    # Make new data array filled with zeros with dimensions time, stations and ens
    for v in variables:
        d = np.zeros((len(ds.timemax), np.shape(ds.x)[0], np.shape(ds.x)[1], nens))
        # Create a dataarray with dimensions time, x, y, ens
        arr = xr.DataArray(data=d,
                           dims=('timemax', 'n', 'm', 'ensemble'),
                           coords={'timemax': ds.timemax, 'x': ds.x, 'y': ds.y, 'ensemble': range(nens)})
        ds[v] = arr

    for iens, file in enumerate(file_list):
        dsin = xr.open_dataset(file)
        for v in variables:
            ds[v][:,:,:,iens] = dsin[v]

    for v in variables:
        arr = ds[v].fillna(-999.0).quantile(prcs, dim="ensemble", skipna=False)
        for ip, p in enumerate(prcs):
            ds[v + "_" + str(round(p*100))] = arr[ip,:,:,:]

    # Remove the original variables
    to_drop = ["zs", "zsmax", "cumprcp", "cuminf", "qinf"]
    for v in to_drop:
        if v in ds:
            ds = ds.drop(v)
        
    try:
        fo.delete_file(output_file_name)
    except:
        pass
    ds.to_netcdf(path= output_file_name)
    ds.close()
