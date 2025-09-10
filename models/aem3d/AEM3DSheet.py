import xarray as xr
import rioxarray as riox
import os
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator as rgi
import utm
import rasterio
import contextily as cx
import numpy as np
from .AEM3DDatetimeCoder import AEM3DDatetimeCoder

class AEM3DSheet:
	""" Class of AEM3D Lake Data Sheet Manipulation Methods """
	
	def open ( file, twist = True):
		ds = xr.open_dataset(file, decode_times={'T': AEM3DDatetimeCoder(msec_resolution=False)})
		if (twist):
			ds = AEM3DSheet.coordtwist(ds)
		else:
			ds = ds
		ds = AEM3DSheet.set_crs(ds)
		return ds
	
	def set_crs(ds):
		loc_utm = utm.from_latlon(ds.attrs['latitude0'], ds.attrs['longitude0'])
		aem3d_crs = rasterio.crs.CRS.from_proj4(f'+proj=tmerc +lon_0=-75 +k_0=0.9996 +x_0={500000-loc_utm[0]} +y_0={0-loc_utm[1]}')
		ds.rio.write_crs(aem3d_crs, inplace=True)
		return ds

	def coordtwist (ds):
		newds = ds.transpose('T', 'X', 'Y', ...)
		newds = newds.rename({'X': 'y', 'Y': 'x'})
		newds = newds.assign_coords(y=-newds.coords['y'])
		newds = newds.rio.set_spatial_dims(x_dim='x', y_dim='y')
		return newds
	
	def cringeterpolate(meshda, cellsize=50, interpmethod="nearest"):
		
		# create the interpolator object
		interper = rgi(points = [meshda['y'],meshda['x']], values = meshda.values, method=interpmethod, bounds_error=False,fill_value = None)

		# generate uniform grid to build
		xreg = list(range(int(cellsize/2), int(meshda['x'].values[-1]), cellsize))
		yreg = list(range(int(-cellsize/2), int(meshda['y'].values[-1]), -cellsize))

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
	
	def interpolate(meshda, cellsize=50, interpmethod="nearest"):
		
		# create the interpolator object
		interper = rgi(points = [meshda['y'],meshda['x']], values = meshda.values, method=interpmethod, bounds_error=False,fill_value = None)

		# generate uniform grid to build
		xreg = list(range(int(cellsize/2), int(meshda['x'].values[-1]), cellsize))
		yreg = list(range(int(-cellsize/2), int(meshda['y'].values[-1]),-cellsize))

		# points list for vectorizing call to interpolater
		plist = []
		for yi in range(0,len(yreg)):
			for xi in range(0,len(xreg)):
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

	def interpolate_stack(da_stack: xr.DataArray, cellsize=50, interpmethod="nearest") -> xr.DataArray:
		'''
		Interpolates a stack of AEM3D rectilinear spatial grids to uniform grids. Performs exactly the same as AEM3DSheet.interpolate(),
		except this function efficiently interpolates every grid in a stack as opposed to a single grid, and is xarray compatible in that
		it perserves meta data like dataset attributes.

		Args:
		-- da_stack: the stack of AEM3D TCHLA time slices; each 'slice' is a 2D grid of the lake representing TCHLA data at one time stamp
		-- cellsize: size in meters of each regular grid cell
		-- interpmethod: method for interpolation forwarded to scipy.interpolate.RegularGridInterpolator()

		Returns:
		An xarray DataArray with interpolated regular spatial grids
		'''

		# extract data for interpolation from the original data array
		values = da_stack.values
		x = da_stack['x'].values
		y = da_stack['y'].values
	
		# Determine output grid
		xreg = np.arange(int(cellsize/2), int(x[-1]), cellsize)
		yreg = np.arange(int(-cellsize/2), int(y[-1]), -cellsize)
		xx, yy = np.meshgrid(xreg, yreg)
		points = np.column_stack([yy.ravel(), xx.ravel()])

		# Loop over each 2D slice in the stack
		grid_list = []
		for i in range(values.shape[0]):
			interpolator = rgi(
				points=[y, x],
				values=values[i],
				method=interpmethod,
				bounds_error=False,
				fill_value=np.nan
			)
			interp_vals = interpolator(points)
			interp_grid = interp_vals.reshape(len(yreg), len(xreg))
			grid_list.append(interp_grid)

		# Stack along axis 0
		interp_data = np.stack(grid_list, axis=0)

		# create the new interpolated data array
		interp_da = xr.DataArray(interp_data, coords={'x': xreg, 'y': yreg, 'date':da_stack['date'].values}, dims=('date', 'y', 'x'), name=da_stack.name)

		### add metadata back to interpolated data array from original data arrray
		# add date attributes back
		interp_da['date'].attrs = da_stack['date'].attrs

		# recompute min and max values for x and y, use original values for the other attributes
		x_attrs = {'MinValue' : da_stack['x'].values.min(), 'MaxValue' : da_stack['x'].values.max()}
		x_attrs = x_attrs | {k:v for k, v in da_stack['x'].attrs.items() if k not in ['MinValue', 'MaxValue']}
		y_attrs = {'MinValue' : da_stack['y'].values.min(), 'MaxValue' : da_stack['y'].values.max()}
		y_attrs = y_attrs | {k:v for k, v in da_stack['y'].attrs.items() if k not in ['MinValue', 'MaxValue']}

		# now add x and y attributes back
		interp_da['x'].attrs = x_attrs
		interp_da['y'].attrs = y_attrs
		
		# add whole dataset attributes (except min and max, since they will be different)
		interp_da.attrs = {k:v for k, v in da_stack.attrs.items() if k not in ['MinValue', 'MaxValue']}

		# add the original CRS back to the interpolated datasest
		interp_da = interp_da.rio.set_spatial_dims(x_dim='x', y_dim='y')
		interp_da = interp_da.rio.write_crs(da_stack.rio.crs)

		return interp_da

	
	def reproject(da, crs):
		targetcrs = rasterio.crs.CRS.from_string(crs)
		sliceproj = da.rio.reproject(targetcrs, resampling=rasterio.enums.Resampling.nearest)
		return sliceproj
	
	def plot (da, crs=None, base=True, title=None, **kwargs):
		""" Plot a sheet variable for a time slice 
			Default to using the dataset CRS
			Reprojection to passed CRS supported
			Default to using time slice 0"""
			
		if (crs != None):
			# reproject original to target crs
			targetcrs = rasterio.crs.CRS.from_string(crs)
			daslice = AEM3DSheet.reproject(da, crs)
		else:
			targetcrs = da.rio.crs
			daslice = da
		 
		fig, axs = plt.subplots()
		# check if there's a _FillValue attribute and mask it if there is
		if ('_FillValue' in daslice.attrs):
			print('Masking by _FillValue')
			daslice.where(daslice<daslice.attrs['_FillValue']).plot(ax=axs, **kwargs)
		else:
			daslice.plot(ax=axs, **kwargs)
		if (base):
			cx.add_basemap(axs, crs=targetcrs)
		
		plt.title(title)
		return(axs)