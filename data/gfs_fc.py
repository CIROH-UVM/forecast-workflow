# import cfgrib
import datetime as dt
from glob import glob
from lib import (download_data,
				 multithreaded_download,
				 multithreaded_loading,
				 parse_to_datetime,
				 get_hour_diff,
				 generate_hours_list)
import more_itertools
import numpy as np
import os
import tempfile as tf
import time
import pandas as pd
import subprocess as sp
import xarray as xr

############# Functions for Processing GFS grib data ############################

def process_gfs_data(location_dict,
					 gfs_data_dir,
					 num_threads):
	"""
	Will load NWM data and return a dictionary of timestamped streamflow data (pd.DataFrame) for each staion name in location_dict

	Args:
	-- location_dict (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to load meterological data for.
	-- gfs_data_dir (str) [req]: the directroy in which all of the GFS grib2 files are stored.
	-- num_threads (int) [req]: number of threds to use for reading grib2 files

	Returns:
	a dictionary of the GFS data in the format {station_ID:pd.Dataframe}
	"""

	station_dict = {}
	grib_file_list = sorted(glob(f'{gfs_data_dir}/gfs.t00z.pgrb2.0p25.f[0-9][0-9][0-9]'))
	# seems like the optimal number of threads to use is 2-4 - any more actually slows down the function
	datasets_dict = multithreaded_loading(cfgrib.open_datasets, grib_file_list, num_threads)
	for grib_file, grib_datasets in datasets_dict.items():
		# grib_datasets = cfgrib.open_datasets(grib_file)
		datasets_indices_to_drop = []
		coords_to_drop = ["step", "atmosphere", "heightAboveGround", "surface", "time", "t"]
		grib_datasets = [ds.drop_vars(coords_to_drop, errors='ignore') for ds in grib_datasets]
		for i, ds in enumerate(grib_datasets):
			for var_name, var in ds.variables.items():
				if var_name == 'tcc' and var.attrs['GRIB_stepType'] == 'avg' and var.attrs['GRIB_typeOfLevel'] == 'atmosphere':
					datasets_indices_to_drop.append(i)
				if var_name == 'prate' and var.attrs['GRIB_stepType'] == 'avg' and var.attrs['GRIB_typeOfLevel'] == 'surface':
					datasets_indices_to_drop.append(i)
		grib_datasets = [grib_datasets[i] for i, _ in enumerate(grib_datasets) if i not in datasets_indices_to_drop]
		merged_ds = xr.merge(grib_datasets)
		# create a downward short-wave radiation flux for the .f000 files (since they don't have one) and set to 0
		if grib_file.endswith('.f000'):
			merged_ds['dswrf'] = 0
		locations_ds = isolate_loc_rows(ds=remap_longs(merged_ds), loc_dict=location_dict)
		locations_df = locations_ds.to_dataframe().reset_index()
		location_groups = locations_df.groupby(["latitude", "longitude"])
		location_dataframes = {name: group for name, group in location_groups}
		append_timestamp(sta_dict=station_dict, loc_dict=location_dict, loc_dfs=location_dataframes)
	return station_dict

def calibrate_columns(df,
					  grib_to_expected_names = {'valid_time':'time', 't2m': 'T2', 'tcc':'TCDC', 'dswrf':'SWDOWN', 'u10':'U10', 'v10':'V10', 'r2':'RH2', 'prate':'RAIN', 'cpofp':'CPOFP'}):
	"""
	Renames and reorders the variable names for GFS dataframes to the expected naming/order convention

	Args:
	-- df (dataframe) [req]: dataframe to modify
	-- grib_to_expected_names (dict) [opt]: dictionary mapping default grib variable names to desired names

	Returns
	The dataframe to append with calibrated columns and index
	"""
	# renaming the columns
	df.rename(grib_to_expected_names, axis=1, inplace=True)
	# Make time the index
	# 20231211 - Add timezone suffix set to UTC
	df.set_index('time', inplace=True)
	df = df.tz_localize('UTC')
	return df

def append_timestamp(sta_dict, loc_dict, loc_dfs):
	"""
	An in-place function that appends each timestamp row (f000, f001, etc) to the station dataframe dict

	Args:
	-- sta_dict (dict) [req]: the dictionary that will contain a dataframe for each station ID
	-- loc_dict (dict) [req]: the dictionary where keys are the station ID and value is the corresponding lat/long tuple
	-- loc_dfs (dict) [req]: the dictionary where keys are the lat/long tuple and values are the corresponding data for the given timestamp (f00, f001, ect.) for that location
	"""
	for stationID, loc in loc_dict.items():
		df_to_append = loc_dfs[loc].drop_duplicates().drop(columns = ['latitude','longitude'])
		# reorder & rename columns to expected convention
		df_to_append = calibrate_columns(df_to_append)
		if stationID in sta_dict:
			sta_dict[stationID] = pd.concat([sta_dict[stationID], df_to_append]).sort_index()
		else:
			sta_dict[stationID] = df_to_append

def dict_to_csv(loc_dict={}, location_dataframes={}):
	for station in loc_dict:
		location = loc_dict[station]
		filename = f"{station}_{location[0]}_{location[1]}.csv"
		location_dataframes[station].to_csv(filename)

def execute(cmd):
	popen = sp.Popen(cmd, stdout=sp.PIPE, universal_newlines=True)
	for stdout_line in iter(popen.stdout.readline, ""):
		yield stdout_line
	popen.stdout.close
	return_code = popen.wait()
	if return_code:
		raise sp.CalledProcessError(return_code, cmd)
		
# pulls just the desired lat/long pairs out of the datasets and isolat
def isolate_loc_rows(ds, loc_dict):
	# initialize empty list to store ds for each station in loc_dict
	station_ds_list = []
	for coords in loc_dict.values():
		# grab the ds for the lat/long pair
		station_ds = ds.sel({"latitude": coords[0], "longitude": coords[1]})
		# add the pulled ds to the station ds list
		station_ds_list.append(station_ds)
	# concatenate all of the station datasets pulled
	concat_stations_ds = xr.concat(station_ds_list, dim="latitude")
	return concat_stations_ds

### In-place function that transforms the longitude indices from 0-360 t0 -180-180
def remap_longs(ds):
	map_function = lambda lon: (lon - 360) if (lon > 180) else lon
	vector_fcn = np.vectorize(map_function)
	longitudes = ds.coords["longitude"].values
	remapped_longitudes = vector_fcn(longitudes)
	remapped_ds = ds.assign_coords(longitude=remapped_longitudes)
	return remapped_ds

################## Function for downloading GFS data ##############################

def download_gfs(dates=dt.datetime.today().strftime("%Y%m%d"),
				 hours=generate_hours_list(168, 'gfs'),
				 grib_data_dir="/data/forecastData/gfs"):
	"""
	Downloads the grib files for the specified dates and hours

	Args:
	-- dates (list of strs) [opt]: list of dates to download gribs for. Default is just the current date
	-- hours (list of strs) [opt]: list of forecast hours to download gribs for. Default is 7 day forecast, or 168 hours
	-- log (logger) [required]: logger track info and download progress
	-- grib_data_dir (str) [opt]: directory in which to store gfs grib files
	
	"""
	print(f'TASK INITIATED: Download {int(hours[-1])}-hour GFS forecasts for the following dates: {dates}')
	for d in dates:
		print(f'DOWNLOADING GFS DATA FOR DATE {d}')
		date_dir = os.path.join(grib_data_dir, f'gfs.{d}/00/atmos')
		if not os.path.exists(date_dir):
			os.makedirs(date_dir)
		for h in hours:
			grib_destination = os.path.join(grib_data_dir, date_dir, f'gfs.t00z.pgrb2.0p25.f{h}')
			grib_url = f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.{d}%2F00%2Fatmos&file=gfs.t00z.pgrb2.0p25.f{h}&var_CPOFP=on&var_DSWRF=on&var_PRATE=on&var_RH=on&var_TCDC=on&var_TMP=on&var_UGRD=on&var_VGRD=on&lev_2_m_above_ground=on&lev_10_m_above_ground=on&lev_surface=on&lev_entire_atmosphere=on&subregion=&toplat=47.5&leftlon=280&rightlon=293.25&bottomlat=40.25"
			download_data(url=grib_url, filepath=grib_destination)
	print('TASK COMPLETE: GFS DOWNLOAD')

def download_gfs_threaded(date,
				 		  hours,
				 		  gfs_data_dir,
						  num_threads):
	"""
	Downloads the grib files for the specified date and hours

	Args:
	-- date (str) [req]: date to download gribs for.
	-- hours (list of strs) [req]: list of forecast hours to download gribs for. Default is 7 day forecast, or 168 hours
	-- gfs_data_dir (str) [req]: directory in which to store gfs grib files
	-- num_threads (int) [req]: number of threads to use.
	"""
	download_list = []
	print(f'TASK INITIATED: Download {int(hours[-1])}-hour GFS forecasts for the following date: {date}')
	# making the directroy in the download function so that it can run independently without error
	if not os.path.exists(gfs_data_dir):
		os.makedirs(gfs_data_dir)
	for h in hours:
		grib_fpath = os.path.join(gfs_data_dir, f'gfs.t00z.pgrb2.0p25.f{h}')
		# if the grib file isn't downloaded already, then download it
		if not os.path.exists(grib_fpath):
			grib_url = f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.{date}%2F00%2Fatmos&file=gfs.t00z.pgrb2.0p25.f{h}&var_CPOFP=on&var_DSWRF=on&var_PRATE=on&var_RH=on&var_TCDC=on&var_TMP=on&var_UGRD=on&var_VGRD=on&lev_2_m_above_ground=on&lev_10_m_above_ground=on&lev_surface=on&lev_entire_atmosphere=on&subregion=&toplat=47.5&leftlon=280&rightlon=293.25&bottomlat=40.25"
			# Ending False is to now use a google bucket -- that's not an option for GFS
			download_list.append((grib_url, grib_fpath, False))
		else:
			print(f'Skipping download; {os.path.basename(grib_fpath)} found at: {gfs_data_dir}')
	# if download_list isn't empty:
	if download_list:
		# The true max hits/min allowed before getting blocked by NOMADS is 120/min, but they've requested users stick to 50/min
		max_hits_per_min = 50
		# break our list into chunks of 50 elements
		download_chunks = list(more_itertools.chunked(download_list, max_hits_per_min))
		for i, dwnld_chunk in enumerate(download_chunks):
			# download one chunked list at a time
			multithreaded_download(dwnld_chunk, num_threads)
			# wait a minute after every download request, except the last one b/c there's no more downloads after the last chunk
			if i < len(download_chunks)-1:
				print('waiting a minute to avoid being flagged by NOMADS')
				time.sleep(60)
	print('TASK COMPLETE: GFS DOWNLOAD')


########################################### GFS get_data() ###########################################

def get_data(forecast_datetime,
			 end_datetime,
			 locations,
			 data_dir=tf.gettempdir(),
			 dnwld_threads=int(os.cpu_count()/2),
			 load_threads=2,
			 return_type='dict'):
	"""
	Download specified GFS forecast data and return nested dictionary of pandas series fore each variable, for each location.

	Args:
	-- forecast_datetime (str, date, or datetime) [req]: the start date and time (00, 06, 12, 18) of the forecast to download. Times are assumed to be UTC time.
	-- end_datetime (str, date, or datetime) [req]: the end date and time for the forecast. GFS forecasts 16-days out for a given start date.
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to download forecast data for.
	-- data_dir (str) [opt]: directory to store donwloaded data. Defaults to OS's default temp directory.
	-- dwnld_threads (int) [opt]: number of threads to use for downloads. Default is half of OS's available threads.
	-- load_threads (int) [opt]: number of threads to use for reading data. Default is 2 for GFS, since file reads are already pretty fast.
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
	GFS forecast data for the given locations in the format specified by return_type
	"""
	forecast_datetime = parse_to_datetime(forecast_datetime)
	end_datetime = parse_to_datetime(end_datetime)

	# grabbing the cycle (12am, 6am, 12pm, 6pm) from the datetime
	forecast_cycle = f"{forecast_datetime.hour:02d}"

	# we need the end_datetime to match the time of forecast_Datetime in order to calculate the number of forecast hours to download accurately
	end_datetime = dt.datetime.combine(end_datetime.date(), forecast_datetime.time()).replace(tzinfo=dt.timezone.utc)
	# calculate the number of hours of forecast data to grab. I.e. for a 5 day forecast, hours would be 120
	forecast_hours = generate_hours_list(get_hour_diff(forecast_datetime, end_datetime), 'gfs')
	forecast_date = forecast_datetime.strftime("%Y%m%d")

	# define the directory for storing GFS data
	# directory should be made in download function, so that it can be used independently if necessary
	gfs_date_dir = os.path.join(data_dir, f'gfs/gfs.{forecast_date}/{forecast_cycle}/atmos')

	# downloading the GFS data with multithreading; 
	download_gfs_threaded(date=forecast_date,
						  hours=forecast_hours,
						  gfs_data_dir=gfs_date_dir,
						  num_threads=dnwld_threads)
	
	# Now, load and process the data
	loc_data = process_gfs_data(location_dict=locations,
								gfs_data_dir=gfs_date_dir,
								num_threads=load_threads)
	# ensure return_type is a valid value
	if return_type not in ['dict', 'dataframe']:
		raise ValueError(f"'{return_type}' is not a valid return_type. Please use 'dict' or 'dataframe'")
	elif return_type == 'dict':
		# created nested dictionary of pd.Series for each variable for each location
		gfs_data = {location:{name:data.dropna().astype('float') for name, data in loc_df.T.iterrows()} for location, loc_df in loc_data.items()}
	elif return_type == 'dataframe':
		raise Exception("'dataframe' option not implemented yet. Please use return_type = 'dict'")
	
	# return the GFS data
	return gfs_data