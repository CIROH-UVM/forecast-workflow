import datetime as dt
from .utils import add_units, get_hour_diff, parse_to_datetime
import pandas as pd
import zarr
import xarray as xr

'''
Data Aquisition Module for NOAA Analysis of Record for Calibration (AORC) Dataset
Amazon Bucket (cuurently implmented): https://registry.opendata.aws/noaa-nws-aorc/
NOAA NWS AORC archive: (not implemented): https://hydrology.nws.noaa.gov/aorc-historic/AORC_NERFC_4km/NERFC_precip_partition/
'''

def get_data(start_date,
			 end_date,
			 locations,
			 variables):
	'''
	A function to download and process AORC data to return nested dictionary of pandas series for each variable, for each location.

		Args:
		-- start_date (str, date, or datetime) [req]: the start date for which to grab data.
		-- end_date (str, date, or datetime) [req]: the end date for which to grab data.
		-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get data for.
		-- variables (dict) [req]: a dictionary of variables to download. Keys should be user-defined names, values should be dataset-specific var names.
	
		Returns:
		AORC timeseries data for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
		are variables names and values are the respective data in a Pandas Series object.
	'''
	# complete lsit of variables available from the AORC bucket
	# variables = {'Total Precipitation':'APCP_surface',
	# 			 'Air Temperature':'TMP_2maboveground',
	# 			 'Specific Humidity':'SPFH_2maboveground',
	# 			 'Downward Long-Wave Radiation Flux':'DLWRF_surface',
	# 			 'Downward Short-Wave Radiation Flux':'DSWRF_surface',
	# 			 'Pressure':'PRES_surface',
	# 			 'U-Component of Wind':'UGRD_10maboveground',
	# 			 'V-Component of Wind':'VGRD_10maboveground'}

	# standardize datetime inputs
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	# remove time zone info dor Datetime index, as AORC zarr files don't specify timezone
	dates = pd.DatetimeIndex([start_date.replace(tzinfo=None)+(dt.timedelta(hours=1)*i) for i in range(0, get_hour_diff(start_date, end_date))])
	years = list(range(start_date.year, end_date.year+1))

	# define the AORC bucket
	bucket = 'noaa-nws-aorc-v1-1-1km'

	# define the fileset for each year requested
	
	# Original Version by Noah Becakge... But ZARR no longer likes S3Map
	# s3_out = s3fs.S3FileSystem(anon=True)
	# fileset = [s3fs.S3Map(
	# 			root=f"s3://{bucket}/{dataset_year}.zarr", s3=s3_out, check=False
	# 		) for dataset_year in years]
	
	# Suggested fix at https://github.com/zarr-developers/zarr-python/issues/2706
	#   But, doesn't allow for anonymous access???
	# fs = fsspec.filesystem('s3', asynchronous=True)
	# fileset = [zarr.storage.FsspecStore(fs, 
	# 	path=f"{bucket}/{dataset_year}.zarr") for dataset_year in years]
	
	# Current suggested solution at https://zarr.readthedocs.io/en/v3.0.0/user-guide/storage.html#remote-store
	#   Also, see https://zarr.readthedocs.io/en/latest/user-guide/v3_migration.html
	fileset = [zarr.storage.FsspecStore.from_url(
		url=f's3://{bucket}/{dataset_year}.zarr',
		# read_only=True,
		storage_options={'anon': True}) for dataset_year in years]
	
	ds_multi_year = xr.open_mfdataset(fileset, engine='zarr')

	# creating lists of lats and lons to grab
	lats = [coord[0] for coord in locations.values()]
	lons= [coord[1] for coord in locations.values()]

	# create list of variable names to drop based of variables dict
	vars_to_drop = [v for v in list(ds_multi_year.data_vars.keys()) if v not in list(variables.values())]

	# filter the dataset for locations and dates
	ds = ds_multi_year.sel(latitude=lats, longitude=lons, time=dates, method='nearest').drop_vars(vars_to_drop)

	# making a dictionary of units for each variable
	var_units = {user_var:ds[ds_var].units for user_var, ds_var in variables.items()}

	# get the dataset-approximated lat and lon values, useful for later grouping
	approx_lats = ds['latitude'].values
	approx_lons = ds['longitude'].values

	# group dataset by latitude
	grouped_ds = ds.groupby('latitude')
	# select and concat only the datasets that correspond to the lat/lon pairs that we need
	filtered_ds = xr.concat([grouped_ds[lat].sel(longitude=lon) for lat, lon in zip(approx_lats, approx_lons)], dim='latitude')

	# convert to a flat dataframe and rename variables to user-defined names
	filtered_df = filtered_ds.to_dataframe().reset_index().set_index(['latitude','longitude']).rename(columns = {ds_name : user_name for user_name, ds_name in variables.items()})

	# now group dataframe by (lat, lon) pairs
	grouped_df = filtered_df.groupby(["latitude", "longitude"])

	# inverted locations dataframe with the approximate coords (from dataset) as keys and user-defined location names as values
	approx_coords = {approx_coords:name for name, approx_coords in zip(locations.keys(), list(zip(approx_lats, approx_lons)))}
		
	# now get dataframes for each location
	# .drop(['latitude','longitude'], axis=1) add this line to the below group manipulations if you want to get rid of lat/lon columns
	location_dataframes = {approx_coords[name]: group.reset_index().set_index('time') for name, group in grouped_df}

	# created nested dictionary of pd.Series for each variable for each location
	aorc_data = {location:{name:data.dropna() for name, data in loc_df.drop(['latitude','longitude'], axis=1).T.iterrows()} for location, loc_df in location_dataframes.items()}

	# add units to nested series dictionary
	add_units(aorc_data, var_units)

	# return the AORC data
	return aorc_data