# -*- coding: utf-8 -*-
"""
FileIO scripts for TEKAL pol and pli files

Created on Mon Sep  6 12:08:35 2021

@author: ormondt
"""

import geopandas as gpd
import pandas as pd
import shapely

import .tekal as tek

class Polyline:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
class PliFile:
    
    def __init__(self, file_name):
        
        self.file_name=file_name
        self.read()

    def read(self):
        
        D = tek.tekal(self.file_name)
        D.info()
        m=D.read(0)
        x = m[0,:,0]
        y = m[1,:,0]        
        self.x = x
        self.y = y

def read_pli_file(file_name):
    
    polylines = []

    D = tek.tekal(file_name)
    D.info()
    for j in range(len(D.blocks)):
        m=D.read(j)
        x = m[0,:,0]
        y = m[1,:,0]
        polylines.append(Polyline(x, y))
        
    return polylines    

def pli2gdf(file_name, crs=None):
    D = tek.tekal(file_name)
    D.info()
    gdf_list = []
    for j in range(len(D.blocks)):
        m = D.read(j)
        x = m[0,:,0]
        y = m[1,:,0]
        line = shapely.geometry.LineString(list(zip(x,y)))
        d = {"geometry": line}
        gdf_list.append(d)
    if crs is not None:     
        gdf = gpd.GeoDataFrame(gdf_list, crs=crs)
    else:
        gdf = gpd.GeoDataFrame(gdf_list)
    return gdf

def pli2geojson(file_name_in, file_name_out, crs=None):
    gdf = pli2gdf(file_name_in, crs=None)
    gdf.to_file(file_name_out, driver='GeoJSON')


def pol2gdf(file_name, crs=None, header=True):
    gdf_list = []
    if not header:
        df = pd.read_csv(file_name, index_col=False, header=None,
              delim_whitespace=True, names=['x', 'y'])
        line = shapely.geometry.Polygon(list(zip(df.x.values, df.y.values)))
        d = {"geometry": line}
        gdf_list.append(d)
    else:
        D = tek.tekal(file_name)
        D.info()
        for j in range(len(D.blocks)):
            m = D.read(j)
            x = m[0,:,0]
            y = m[1,:,0]
            line = shapely.geometry.Polygon(list(zip(x,y)))
            d = {"geometry": line}
            gdf_list.append(d)

    if crs is not None:     
        gdf = gpd.GeoDataFrame(gdf_list, crs=crs)
    else:
        gdf = gpd.GeoDataFrame(gdf_list)
    return gdf

def pol2geojson(file_name_in, file_name_out, crs=None):
    gdf = pol2gdf(file_name_in, crs=None)
    gdf.to_file(file_name_out, driver='GeoJSON')

def gdf2pli(gdf, file_name, header=True):
    if gdf.crs.is_geographic:
        fid = open(file_name, "w")
        for index, row in gdf.iterrows():
            nrp = len(row["geometry"].coords)
            if header:
                fid.write("BL" + str(index + 1).zfill(4) + "\n")
                fid.write(str(nrp) + " " + "2\n")
            for ip in range(nrp):
                x = row["geometry"].coords[ip][0]
                y = row["geometry"].coords[ip][1]
                string = f'{x:12.6f}{y:12.6f}\n'
                fid.write(string)
        fid.close()
    else:
        fid = open(file_name, "w")
        for index, row in gdf.iterrows():
            nrp = len(row["geometry"].coords)
            if header:
                fid.write("BL" + str(index).zfill(4) + "\n")
                fid.write(str(nrp) + " " + "2\n")
            for ip in range(nrp):
                x = row["geometry"].coords[ip][0]
                y = row["geometry"].coords[ip][1]
                string = f'{x:12.1f}{y:12.1f}\n'
                fid.write(string)
        fid.close()

def gdf2pol(gdf, file_name, header=True):
    if gdf.crs.is_geographic:
        fid = open(file_name, "w")
        for index, row in gdf.iterrows():
            nrp = len(row["geometry"].exterior.coords)
            if header:
                fid.write("BL" + str(index + 1).zfill(4) + "\n")
                fid.write(str(nrp) + " " + "2\n")
            for ip in range(nrp):
                x = row["geometry"].exterior.coords[ip][0]
                y = row["geometry"].exterior.coords[ip][1]
                string = f'{x:12.6f}{y:12.6f}\n'
                fid.write(string)
        fid.close()
    else:
        fid = open(file_name, "w")
        for index, row in gdf.iterrows():
            nrp = len(row["geometry"].exterior.coords)
            if header:
                fid.write("BL" + str(index).zfill(4) + "\n")
                fid.write(str(nrp) + " " + "2\n")
            for ip in range(nrp):
                x = row["geometry"].exterior.coords[ip][0]
                y = row["geometry"].exterior.coords[ip][1]
                string = f'{x:12.1f}{y:12.1f}\n'
                fid.write(string)
        fid.close()
