import xarray as xr
import rioxarray as riox
import os
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator as rgi
import utm
import rasterio
import contextily as cx
import numpy as np

class AEM3DSheet:
    """ Class of AEM3D Lake Data Sheet Manipulation Methods """
    
    def open ( file, twist = True) : 
        ds = xr.open_dataset(file)
        if (twist) :     
            ds = AEM3DSheet.coordtwist(ds)
        else :
            ds = ds
        ds = AEM3DSheet.set_crs(ds)
        return ds
    
    def set_crs(ds):
        loc_utm = utm.from_latlon(ds.attrs['latitude0'], ds.attrs['longitude0'])
        aem3d_crs = rasterio.crs.CRS.from_proj4(f'+proj=tmerc +lon_0=-75 +k_0=0.9996 +x_0={500000-loc_utm[0]} +y_0={0-loc_utm[1]}')
        ds.rio.write_crs(aem3d_crs, inplace=True)
        return ds

    def coordtwist (ds):
        newds = ds.transpose('T', 'X', 'Y',...)
        newds = newds.rename({'X': 'y', 'Y': 'x'})
        newds = newds.assign_coords(y=-newds.coords['y'])
        newds = newds.rio.set_spatial_dims(x_dim='x', y_dim='y')
        return newds
    
    def cringeterpolate(meshda, cellsize=50, interpmethod="nearest") :
        
        # create the interpolator object
        interper = rgi(points = [meshda['y'],meshda['x']], values = meshda.values, method=interpmethod, bounds_error=False,fill_value = None)

        # generate uniform grid to build
        xreg = list(range(int(cellsize/2), int(meshda['x'].values[-1]), cellsize))
        yreg = list(range(int(-cellsize/2), int(meshda['y'].values[-1]),-cellsize))

        # initialize size of dataarray to hold interpolated data
        interpda = xr.DataArray(coords=(yreg,xreg))
        #  the brute force double loop builds the array that the mapping function currently works with
        for yi in range(0,len(interpda['dim_0'])-1) :
            for xi in range(0,len(interpda['dim_1'])-1) :
                interpda[yi,xi] = interper((interpda['dim_0'][yi],interpda['dim_1'][xi]))
    
        interpda = interpda.rename({'dim_0': 'y', 'dim_1': 'x'})
        interpda = interpda.rio.set_spatial_dims(x_dim='x', y_dim='y')

        interpda = interpda.rio.write_crs(meshda.rio.crs)
    
        return interpda
    
    def interpolate(meshda, cellsize=50, interpmethod="nearest") :
        
        # create the interpolator object
        interper = rgi(points = [meshda['y'],meshda['x']], values = meshda.values, method=interpmethod, bounds_error=False,fill_value = None)

        # generate uniform grid to build
        xreg = list(range(int(cellsize/2), int(meshda['x'].values[-1]), cellsize))
        yreg = list(range(int(-cellsize/2), int(meshda['y'].values[-1]),-cellsize))

        # points list for vectorizing call to interpolater
        plist = []
        for yi in range(0,len(yreg)) :
            for xi in range(0,len(xreg)) :
                plist.append([yreg[yi],xreg[xi]])
                
        # vectorized interpolater invocation
        interplist=interper(plist)   
        # reshape to 2D array
        interprast=interplist.reshape(len(yreg),len(xreg))
        
        interpda = xr.DataArray(
                    data=interprast,  # Your NumPy array of values
                    coords={
                        "x": xreg,  # Dimension name and corresponding coordinates
                        "y": yreg  # Dimension name and corresponding coordinates
                        },                   
                    dims=("y", "x") # Specify the dimension names in the correct order
        )         

        interpda = interpda.rio.set_spatial_dims(x_dim='x', y_dim='y')
        interpda = interpda.rio.write_crs(meshda.rio.crs)
    
        return interpda
    
    def reproject(da, crs):
        targetcrs = rasterio.crs.CRS.from_string(crs)
        sliceproj = da.rio.reproject(targetcrs, resampling=rasterio.enums.Resampling.nearest)
        return sliceproj
    
    def plot (da, crs=None, base=True, title=None, **kwargs) :
        """ Plot a sheet variable for a time slice 
            Default to using the dataset CRS
            Reprojection to passed CRS supported
            Default to using time slice 0"""
            
        if (crs != None) :
            # reproject original to target crs
            targetcrs = rasterio.crs.CRS.from_string(crs)
            daslice = AEM3DSheet.reproject(da, crs)
        else :
            targetcrs = da.rio.crs
            daslice = da
         
        fig, axs = plt.subplots()
        # check if there's a _FillValue attribute and mask it if there is
        if ('_FillValue' in daslice.attrs) :
            print('Masking by _FillValue')
            daslice.where(daslice<daslice.attrs['_FillValue']).plot(ax=axs, **kwargs)
        else :
            daslice.plot(ax=axs, **kwargs)
        if (base) :
            cx.add_basemap(axs, crs=targetcrs)
        
        plt.title(title)
        return(axs)