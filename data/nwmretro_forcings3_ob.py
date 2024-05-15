import datetime as dt
from more_itertools import chunked
import numpy as np
import pandas as pd
import pyproj
from s3fs import S3FileSystem, S3Map
import xarray as xr
from lib import get_hour_diff, parse_to_datetime
from aiohttp import ServerDisconnectedError


# stick with print statements
# Question: In terms of best practicee, should we use print() or logger.info() statements for log messages

# Data source: https://ciroh-nwm-zarr-retrospective-data-copy.s3.amazonaws.com/index.html


def get_data(start_date,
			 end_date,
			 locations,
			 variables):
	'''
	A function to download and process NWM Retro Forcings data to return nested dictionary of pandas series for each variable, for each location.

		Args:
		-- start_date (str, date, or datetime) [req]: the start date for which to grab data
		-- end_date (str, date, or datetime) [req]: the end date for which to grab data
		-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get data for.
		-- variables (dict) [req]: a dictionary of variables to download. Keys should be user-defined var names, value should be dataset-specific var names
		
		Returns:
		NWM retrospective forcings timeseries for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
		are variables names and values are the respective data in a Pandas Series object.
	'''
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	print(f"BEGIN NWM RETRO FORCING DATA GRAB FOR DATES [{start_date.strftime('%Y-%m-%d %H:%M:%S')} TO {end_date.strftime('%Y-%m-%d %H:%M:%S')})")
	print(f'VARIABLES TO EXTRACT: {list(variables.values())}')

	# NOTE: there are 24 timesteps for each day, 00-23
 	# define the amazon web bucket you want to use
	bucket = 's3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/forcing/'
	fs = S3FileSystem(anon=True)
	
	dates = [start_date+(dt.timedelta(hours=1)*i) for i in range(0, get_hour_diff(start_date, end_date))]
	
    # Convert lat, lon pairs to x,y in NWM forcings projection
    # define the lat lon crs
	wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
	# define the nwm crs - Lambert Conformal Conic Projection
	wrf_proj = pyproj.Proj(proj='lcc',
						lat_1=30.,
						lat_2=60., 
						lat_0=40., lon_0=-97., # Center point
						a=6370000, b=6370000)

	# transform lat, lon into X, Y
	to_wrf_x_y = pyproj.Transformer.from_crs(wgs_proj.crs, wrf_proj.crs)
	
	# Create an empty list to store the results for x and y
	y_indices = []
	x_indices = []

	return_data = {}
	for loc, coords in locations.items():
		return_data[loc] = {}

		lat, lon = coords
		X, Y = to_wrf_x_y.transform(lon, lat)

		# Append the result to the list
		x_indices.append(X)
		y_indices.append(Y)

		print(f'FOR LOCATION: {loc}')
		print(f"\tLatitude, Longitude {(lat, lon)} is {(X, Y)} in NWM V3.0 forcing X,Y projection")

	zarr_variable_fn = {'RAINRATE': 'precip',
						'T2D': 't2d'}
	
	for var, nwm_name in variables.items():
		ds = xr.open_dataset(S3Map(f"s3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/forcing/{zarr_variable_fn[nwm_name]}.zarr", s3=fs), engine='zarr')
		# print(f"Reducing time slice to {start_date.strftime('%Y-%m-%dT%H:%M:%S')} - {end_date.strftime('%Y-%m-%dT%H:%M:%S')}")
		ds = ds.sel(time=slice(start_date.strftime('%Y-%m-%dT%H:%M:%S'), end_date.strftime('%Y-%m-%dT%H:%M:%S')))
		ds = ds.sel({'x': x_indices, 'y': y_indices}, method = "nearest")
		# print("ds after select")
		# print(ds)

		locations_xy = {loc: (x,y) for loc,x,y in zip(locations.keys(), ds.coords['x'].values,  ds.coords['y'].values)}
		# print(locations_xy)

		# convert to dataframe
		while True:
			try:
				print(f'ATTEMPTING DATA ACQUISITION FOR {var}...')
				df = ds.drop_vars(['crs']).to_dataframe()
				break
			except ServerDisconnectedError as sde:
				print("Server disconnected, retrying in 5 seconds...")
			except Exception as e:
				raise(e)
		ds.close()
		print(f"DATA ACQUISITION COMPLETE... ACQUIRED {df.shape} DATAFRAME")

		# set proper index and drop dataset dimension columns
		# indexed_df = df.reset_index().set_index(['time','x','y'])
		# group by x and y coordinates
		location_groups = df.groupby(["x", "y"])
		# extract the locations of interest and drop x,y coordinates after extraction
		for name, xy in locations_xy.items():
			print(f"PROCESSING location {name} at {xy}...\n")
			# Get the group for each location and rename variables
			this_group = location_groups.get_group(xy).droplevel(['x','y']).rename({nwm_name: var}, axis='columns')[var]
			# Drop duplicated time stamps (from duplicated x,y coordinates) and drop end_date to make it exclusive in the range [start, end)
			return_data[name][var] = this_group[~this_group.index.duplicated(keep='first')].drop(end_date.strftime('%Y-%m-%dT%H:%M:%S'))
	
	return return_data