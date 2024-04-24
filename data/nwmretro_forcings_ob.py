import datetime as dt
from kerchunk.combine import MultiZarrToZarr
from more_itertools import chunked
import numpy as np
import pandas as pd
import pyproj
from s3fs import S3FileSystem
import xarray as xr
from lib import get_hour_diff, parse_to_datetime

# stick with print statements
# Question: In terms of best practicee, should we use print() or logger.info() statements for log messages

# Data source: https://ciroh-nwm-zarr-retrospective-data-copy.s3.amazonaws.com/index.html


def get_data(start_date,
			 end_date,
			 locations,
			 variables,
			 return_type='dict'):
	'''
	A function to download and process NWM Retro Forcings data to return nested dictionary of pandas series for each variable, for each location.

		Args:
		-- start_date (str, date, or datetime) [req]: the start date for which to grab data
		-- end_date (str, date, or datetime) [req]: the end date for which to grab data
		-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get data for.
		-- variables (dict) [req]: a dictionary of variables to download.
		-- return_type (string) [opt]: string indicating which format to return data in. Default is "dict", which will return data in a nested dict format:
										{locationID1:{
											var1_name:pd.Series,
											var2_name:pd.Series,
											...},
										locationID2:{...},
										...
										}
										Alternative return type is "dataframe", which smashes all data into a single dataframe muliIndex'd by station ID, then timestamp
		
		Returns:
		
	'''
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	print(f'BEGIN NWM RETRO FORCING DATA GRAB FOR DATE {start_date.strftime("%m/%d/%Y:%H")} TO {end_date.strftime("%m/%d/%Y:%H")}')
	print(f'VARIABLES TO EXTRACT: {list(variables.values())}')

	# NOTE: there are 24 timesteps for each day, 00-23
 	# define the amazon web bucket you want to use
	bucket = 's3://ciroh-nwm-zarr-retrospective-data-copy/noaa-nwm-retrospective-2-1-zarr-pds/forcing/'
	s3 = S3FileSystem(anon=True)
	
	dates = [start_date+(dt.timedelta(hours=1)*i) for i in range(0, get_hour_diff(start_date, end_date))]
	bucket_parts = [f'/{d}' for d in bucket.split("/")]
	files = [f'https:/{bucket_parts[2]}.s3.amazonaws.com{"".join(bucket_parts[3:])}{d.year}/{d:%Y%m%d%H}.LDASIN_DOMAIN1.json' for d in dates]
	
	# Build a multiZarr only for getting the structure of the data
	jsonlist = files[0:5]
	mzz = MultiZarrToZarr(jsonlist,
		remote_protocol='s3',
		remote_options={'anon':True},
		concat_dims=['valid_time'])
	d = mzz.translate()
	backend_args = {"consolidated": False, "storage_options": {"fo": d}, "consolidated": False}
	ds = xr.open_dataset("reference://", engine="zarr", backend_kwargs=backend_args)
	# Remove extra Time dimension of size 1
	ds = ds.squeeze(dim='Time')

	# Add spatial metadata to the dataset
	# Load the metadata dataset using xarray and add spatial metadata to it.
	ds_meta = xr.open_dataset('http://thredds.hydroshare.org/thredds/dodsC/hydroshare/resources/2a8a3566e1c84b8eb3871f30841a3855/data/contents/WRF_Hydro_NWM_geospatial_data_template_land_GIS.nc', engine='netcdf4')
	x = ds_meta.x.values
	y = ds_meta.y.values
	ds = ds.rename({'valid_time': 'time', 'south_north':'y', 'west_east':'x'})
	X, Y = np.meshgrid(x, y)
	# define the input crs - Lambert Conformal Conic Projection
	wrf_proj = pyproj.Proj(proj='lcc',
						lat_1=30.,
						lat_2=60., 
						lat_0=40.0000076293945, lon_0=-97., # Center point
						a=6370000, b=6370000)

	# define the output crs
	wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
	# transform X, Y into Lat, Lon
	transformer = pyproj.Transformer.from_crs(wrf_proj.crs, wgs_proj.crs)
	lon, lat = transformer.transform(X, Y)

	ds = ds.assign_coords(lon = (['y', 'x'], lon))
	ds = ds.assign_coords(lat = (['y', 'x'], lat))
	ds = ds.assign_coords(x = x)
	ds = ds.assign_coords(y = y)
	ds.x.attrs['axis'] = 'X'
	ds.x.attrs['standard_name'] = 'projection_x_coordinate'
	ds.x.attrs['long_name'] = 'x-coordinate in projected coordinate system'
	ds.x.attrs['resolution'] = 1000.  # cell size
	ds.y.attrs['axis'] = 'Y' 
	ds.y.attrs['standard_name'] = 'projection_y_coordinate'
	ds.y.attrs['long_name'] = 'y-coordinate in projected coordinate system'
	ds.y.attrs['resolution'] = 1000.  # cell size
	ds.lon.attrs['units'] = 'degrees_east'
	ds.lon.attrs['standard_name'] = 'longitude' 
	ds.lon.attrs['long_name'] = 'longitude'
	ds.lat.attrs['units'] = 'degrees_north'
	ds.lat.attrs['standard_name'] = 'latitude'
	ds.lat.attrs['long_name'] = 'latitude'
	
	# add crs to netcdf file
	ds.rio.write_crs(ds_meta.crs.attrs['spatial_ref'],
				  	 inplace=True).rio.set_spatial_dims(x_dim="x",
                                       					y_dim="y",
                                       					inplace=True,
                                       					).rio.write_coordinate_system(inplace=True)
	
	# Create an empty list to store the results for x and y
	y_indices = []
	x_indices = []

	locations_xy = {}
	for loc, coords in locations.items():
		lat, lon = coords
		# Lets compute the difference between the desired lat longs against every lat long
		lon_diff = abs(ds['lon'] - lon)
		lat_diff  = abs(ds['lat'] - lat)
		
		# Now filter the values with a minumum difference
		min_lon_diff = lon_diff.where(lon_diff == lon_diff.min(), drop=True)
		min_lat_diff  = lat_diff.where(lat_diff == lat_diff.min(), drop=True)

		# print(min_lon_diff)
		# At this point, we have the index points which we can use to filter the data points at that lat long. However, we need to further refine the 
		# index points. Lets exactly get those points

		y_index = np.where(ds['y'].values == min_lat_diff['y'].values.item())[0][0]
		x_index = np.where(ds['x'].values == min_lon_diff['x'].values.item())[0][0]
			
		# Append the result to the list
		y_indices.append(y_index)
		x_indices.append(x_index)
		
		# define locations dict with x, y coords
		locations_xy[loc] = (x_index, y_index)

		print(f'FOR LOCATION: {loc}')
		print(f"\tLatitude, Longitude is {(lat, lon)} and nearest x, y coordinates are {(x_index, y_index)}")

	mzz = MultiZarrToZarr(files,
	remote_protocol='s3',
	remote_options={'anon':True},
	concat_dims=['valid_time'])

	d = mzz.translate()
	backend_args = {'consolidated': False, 'storage_options': {'fo': d}}
	ds = xr.open_dataset('reference://', engine='zarr', backend_kwargs=backend_args)

	# create list of variable names to drop based of variables dict
	vars_to_drop = [v for v in list(ds.data_vars.keys()) if v not in list(variables.values())]

	# Remove extra Time dimension of size 1 and remove unwated variables
	ds = ds.squeeze(dim='Time').drop_vars(vars_to_drop)
	# create x and y coordinates from dimensions
	ds = ds.assign_coords(x = ds.west_east)
	ds = ds.assign_coords(y = ds.south_north)
	# rename valid_time coord to time
	ds = ds.rename({'valid_time': 'time'})
	# select the locations we need
	ds = ds.sel({'south_north':y_indices, 'west_east':x_indices})
	# convert to dataframe
	while True:
		try:
			print('ATTEMPTING DATA ACQUISITION...')
			df = ds.to_dataframe()
			break
		except ServerDisconnectedError as sde:
			print("Server disconnected, retrying in 5 seconds...")
		except Exception as e:
			raise(e)
	
	print("DATA ACQUISITION COMPLETE")

	# set proper index and drop dataset dimension columns
	indexed_df = df.reset_index().set_index(['time','x','y']).drop(['south_north','west_east'], axis=1)
	# group by x and y coordinates
	location_groups = indexed_df.groupby(["x", "y"])
	# extract the locations of interest and drop x,y coordinates after extraction
	location_dataframes = {name : location_groups.get_group(xy).droplevel(['x','y']) for name, xy in locations_xy.items()}
	
	# ensure return_type is a valid value
	if return_type not in ['dict', 'dataframe']:
		raise ValueError(f"'{return_type}' is not a valid return_type. Please use 'dict' or 'dataframe'")
	elif return_type == 'dict':
		# created nested dictionary of pd.Series for each variable for each location
		nwm_retro_forcings = {location:{name:data.dropna() for name, data in loc_df.T.iterrows()} for location, loc_df in location_dataframes.items()}
	elif return_type == 'dataframe':
		raise Exception("'dataframe' option not implemented yet. Please use return_type = 'dict'")
	
	return nwm_retro_forcings