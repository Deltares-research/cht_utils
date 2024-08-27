# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 11:29:26 2022

@author: mvano
"""
import os
import argparse
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union, transform
import geopandas as gpd
import fiona
import pyproj

def get_watersheds_in_bbox(gdb_path, output_path, bbox, layer, name):

    gdb_file = os.path.join(gdb_path, "WBD_National_GDB.gdb")
    
    if not output_path:
        output_path = os.getcwd()
    layers = fiona.listlayers(os.path.join(gdb_path, "WBD_National_GDB.gdb"))

    print('layer=' + layer)
    gdf = gpd.read_file(gdb_file, layer=layer)
    
    if layer=="WBDHU2":
        hucstr = "huc2"
    elif layer=="WBDHU4":
        hucstr = "huc4"
    elif layer=="WBDHU6":
        hucstr = "huc6"
    elif layer=="WBDHU8":
        hucstr = "huc8"
    elif layer=="WBDHU10":
        hucstr = "huc10"
    elif layer=="WBDHU12":
        hucstr = "huc12"
    elif layer=="WBDHU14":
        hucstr = "huc14"
    elif layer=="WBDHU16":
        hucstr = "huc16"
    
    hucs = gdf[hucstr]
    watersheds = gdf["geometry"]
    
    fid = open(os.path.join(output_path, name + "_watersheds_" + layer + ".pol"), "w")
    
    for j, watershed in enumerate(watersheds):
        if watershed.bounds[0]<bbox[2] and watershed.bounds[2]>bbox[0] and watershed.bounds[1]<bbox[3] and watershed.bounds[3]>bbox[1]:
            gs = watershed.geoms
            for polygon in gs:
                crds = polygon.exterior.coords
                huc = hucs[j]
                fid.write(huc + "\n")
                fid.write(str(len(crds)) + ' 2\n')
                for crd in crds:
                    lon = crd[0]
                    lat = crd[1]
                    fid.write(f'{lon:11.6f}' + ' ' + f'{lat:11.6f}' + '\n')
    
    fid.close()

def create_include_polygon(gdb_path, huc_file, output_path, name, layer="WBDHU8", dbuf=0.0, crs=None):

    if not output_path:
        output_path = os.getcwd()

    if crs:
        try:
            # Try epsg code
            crs = int(crs)            
        except:
            pass
        wgs84 = pyproj.CRS('EPSG:4326')
        crs1  = pyproj.CRS(crs)
        if crs1.is_geographic:
            # Convert metres to degrees
            dbuf = dbuf / 100000.0
        # Construct transformer    
        transformer = pyproj.Transformer.from_crs(wgs84,
                                                  crs1,
                                                  always_xy=True).transform

    # Read list of huc8 IDs    
    huc_list = []
    fid = open(huc_file, 'r')
    lines = fid.readlines()
    for line in lines:
        huc_list.append(line.strip())

    # Read gdb file and make list of polygons
    gdb_file = os.path.join(gdb_path, "WBD_National_GDB.gdb")
    source = fiona.open(gdb_file, driver='OpenFileGDB', layer=layer)
    polys = []
    wnames = []

    if layer=="WBDHU2":
        hucstr = "huc2"
    elif layer=="WBDHU4":
        hucstr = "huc4"
    elif layer=="WBDHU6":
        hucstr = "huc6"
    elif layer=="WBDHU8":
        hucstr = "huc8"
    elif layer=="WBDHU10":
        hucstr = "huc10"
    elif layer=="WBDHU12":
        hucstr = "huc12"
    elif layer=="WBDHU14":
        hucstr = "huc14"
    elif layer=="WBDHU16":
        hucstr = "huc16"
    
    for f in source:
        if f["properties"][hucstr] in huc_list:
            print(f["properties"][hucstr])
            geom = f["geometry"]
            wnames.append(f["properties"]["name"])
            for crds in geom["coordinates"]:
                p = Polygon(crds[0])
                polys.append(p)
    source.close()    

    # Merge polygons
    merged = unary_union(polys)            

    crsgeo = True
    if crs:    
        merged = transform(transformer, merged)        
        if crs1.is_projected:
            crsgeo = False

    # Apply buffer
    if dbuf>0.0:
        merged = merged.buffer(dbuf, resolution=16)
        merged = merged.simplify(dbuf)

    filename_pol = os.path.join(output_path, name + ".pol")
    filename_txt = os.path.join(output_path, name + "_watershed_names.txt")

    # Write text file with watershed names
    fit = open(filename_txt, "w")
    for wname in wnames:
        fit.write(wname + '\n')
    fit.close()

    # Write polygon file

    fid = open(filename_pol, "w")

    if isinstance(merged, Polygon):

        crds = merged.exterior.coords
        fid.write("BL01\n")
        fid.write(str(len(crds)) + ' 2\n')
        for crd in crds:
            lon = crd[0]
            lat = crd[1]
            if crsgeo:
                fid.write(f'{lon:11.6f}' + ' ' + f'{lat:11.6f}' + '\n')
            else:
                fid.write(f'{lon:11.1f}' + ' ' + f'{lat:11.1f}' + '\n')

    else:
        
        for geom in merged.geoms:
            crds = geom.exterior.coords
            fid.write("BL01\n")
            fid.write(str(len(crds)) + ' 2\n')
            for crd in crds:
                lon = crd[0]
                lat = crd[1]
                if crsgeo:
                    fid.write(f'{lon:11.6f}' + ' ' + f'{lat:11.6f}' + '\n')
                else:
                    fid.write(f'{lon:11.1f}' + ' ' + f'{lat:11.1f}' + '\n')
            
    fid.close()
        
    return merged

def main():

        # e.g. python watersheds get c:\\work\\projects\\delftdashboard\\wbd -o . -bbox -79.0 35.0 -75.0 41.0 -l WBDHU8 -name chesapeake 
        # e.g. python watersheds merge -wpath c:\\work\\projects\\delftdashboard\\wbd id_list.txt --outputpath . --buffer 200.0 --crs "WGS 84" 

	# create parser object
    parser = argparse.ArgumentParser(description = "Bathymetry/roughness dataset tiler")

	# defining arguments for parser object
    parser.add_argument("command")
	
    parser.add_argument("wpath", type = str,
                        help = "GDB path with WBD database.")

    parser.add_argument("-i","--idfile", type = str,
                        help = "File with list of IDs to be merged.")
	
    parser.add_argument("-o", "--outputpath", type = str,
                         metavar = "path", default = None,
                         help = "Output path.")
	
    parser.add_argument("-b", "--bbox", type = float, nargs = 4,
	                    metavar = ("lonmin", "latmin", "lonmax", "latmax"),
	                    default = [-999.0, -999.0, 999.0, 999.0],
                        help = "Bounding box coordinates(e.g. -79.0 35.0 -75.0 41.0.")
	
    parser.add_argument("-l", "--layer", type = str,
                        metavar = "layer", default = "WBDHU8",
                        help = "WBD layer (e.g. WBDHU8).")

    parser.add_argument("-n", "--name", type = str,
                        metavar = "name", default = "wbd",
                        help = "Output name.")

    parser.add_argument("-c", "--crs", type = str,
                        metavar = "crs", default = "WGS 84",
                        help = "Name of CRS (e.g. WGS 84).")
    
    parser.add_argument("-buffer", "--buffer", type = float,
                        metavar = "buffer", default = 0.0,
                        help = "Buffer around watershed polygon in metres.")

	# parse the arguments from standard input
    args = parser.parse_args()
	
	# calling functions depending on type of argument
    if args.command == "get":
        get_watersheds_in_bbox(args.wpath, args.outputpath, args.bbox, args.layer, args.name)
    elif args.command == "merge":
        create_include_polygon(args.wpath, args.idfile, args.outputpath, args.name, dbuf=args.buffer, crs=args.crs, layer=args.layer)
    else:
        print("Error! Unknown command " + args.command)

if __name__ == "__main__":
	# calling the main function
	main()
        