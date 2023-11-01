"""
This module will extract raster values at points of interest.

One could do this in ArcGIS Pro, but for multiband rasters, the process is 
cumbersome.  This script takes paths to a shapefile with points of interest
and a raster as inputs.  The script currently just prints the results, however,
I suggest you use pandas to log results and save to CSV.  Not knowing how this
script will be used, I left it to simply print.

To use the script, you may either call it from the command line or import
tif_extracter as a submodule.

Command line example:

python tif_extracter.py "C:\path\to\pts.shp" "C:\path\to\raster.tif"

"""

import sys
from osgeo import gdal,ogr


def extract_points(point_path, raster_path):
    # Load raster data
    raster = gdal.Open(raster_path) 
    band_count = raster.RasterCount 
    transform = raster.GetGeoTransform()

    # Load shapefile data
    ds = ogr.Open(point_path)
    lyr = ds.GetLayer()

    # Iterate through all points
    for feat in lyr:
        # Get raster indices of point
        geom = feat.GetGeometryRef()
        mx,my = geom.GetX(), geom.GetY()
        px = int((mx - transform[0]) / transform[1])
        py = int((my - transform[3]) / transform[5])
        
        # Get values at raster indices
        for i in range(band_count):
            tmp_band = raster.GetRasterBand(int(i + 1))
            val = tmp_band.ReadAsArray(px,py,1,1)[0][0]

            # Currently set up just to print output.  I suggest using pandas 
            # to log info and save to csv.
            print(f'{tmp_band.GetDescription()}: {val}')

if __name__ == '__main__':
    points = sys.argv[1]
    raster = sys.argv[2]
    extract_points(points, raster)
