import datetime as dt
from data.gfs_fc import isolate_loc_rows, remap_longs
from glob import glob
from lib import (multithreaded_download,
				 multithreaded_loading,
				 parse_to_datetime,
				 get_hour_diff,
				 generate_hours_list)
import more_itertools
import os
import tempfile as tf
import pandas as pd
import re
import time
import xarray as xr

def append_timestamp(sta_dict, loc_dict, loc_dfs, timestamp):
	"""
	An in-place function that appends each timestamp row (f000, f001, etc) to the station dataframe dict

	Args:
	-- sta_dict (dict) [req]: the dictionary that will contain a dataframe for each station ID
	-- loc_dict (dict) [req]: the dictionary where keys are the station ID and value is the corresponding lat/long tuple
	-- loc_dfs (dict) [req]: the dictionary where keys are the lat/long tuple and values are the corresponding data for the given timestamp (f00, f001, ect.) for that location
	"""
	for stationID, loc in loc_dict.items():
		df_to_append = loc_dfs[loc].drop_duplicates().drop(columns = ['latitude','longitude'])
		df_to_append['time'] = timestamp
		df_to_append.set_index('time', inplace=True)
		if stationID in sta_dict:
			sta_dict[stationID] = pd.concat([sta_dict[stationID], df_to_append]).sort_index()
		else:
			sta_dict[stationID] = df_to_append

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
	print(f'TASK INITIATED: Download {int(hours[-1])}-hour GFS forecasts for the following date: {date.strftime("%Y%m%d")}')
	# making the directroy in the download function so that it can run independently without error
	if not os.path.exists(gfs_data_dir):
		os.makedirs(gfs_data_dir)
	# grabbing the cycle (12am, 6am, 12pm, 6pm) from the datetime
	fc_cycle = f"{date.hour:02d}"
	for h in hours:
		fpath = os.path.join(gfs_data_dir, f'gfs.0p25.{date.strftime("%Y%m%d")}{fc_cycle}.f{h}.grib2.nc')
		# if the grib file isn't downloaded already, then download it
		if not os.path.exists(fpath):
			ts_date = date.replace(tzinfo=None) + dt.timedelta(days=int(int(h)/24), hours=int(h)%24)
			if h == '000':
				url = f'https://thredds.rda.ucar.edu/thredds/ncss/grid/files/g/ds084.1/{date.year}/{date.strftime("%Y%m%d")}/gfs.0p25.{date.strftime("%Y%m%d")}{fc_cycle}.f{h}.grib2?var=u-component_of_wind_height_above_ground&var=v-component_of_wind_height_above_ground&var=Temperature_height_above_ground&var=Relative_humidity_height_above_ground&var=Per_cent_frozen_precipitation_surface&var=Precipitation_rate_surface&var=Total_cloud_cover_entire_atmosphere&north=47.5&west=280&east=293.25&south=40.25&horizStride=1&time_start={ts_date.isoformat()}Z&time_end={ts_date.isoformat()}Z&&&accept=netcdf4-classic'
			elif int(h) % 6 == 3:
				url = f'https://thredds.rda.ucar.edu/thredds/ncss/grid/files/g/ds084.1/{date.year}/{date.strftime("%Y%m%d")}/gfs.0p25.{date.strftime("%Y%m%d")}{fc_cycle}.f{h}.grib2?var=Downward_Short-Wave_Radiation_Flux_surface_3_Hour_Average&var=u-component_of_wind_height_above_ground&var=v-component_of_wind_height_above_ground&var=Temperature_height_above_ground&var=Relative_humidity_height_above_ground&var=Per_cent_frozen_precipitation_surface&var=Precipitation_rate_surface&var=Total_cloud_cover_entire_atmosphere&north=47.5&west=280&east=293.25&south=40.25&horizStride=1&time_start={ts_date.isoformat()}Z&time_end={ts_date.isoformat()}Z&&&accept=netcdf4-classic'
			elif int(h) % 6 == 0:
				url = f'https://thredds.rda.ucar.edu/thredds/ncss/grid/files/g/ds084.1/{date.year}/{date.strftime("%Y%m%d")}/gfs.0p25.{date.strftime("%Y%m%d")}{fc_cycle}.f{h}.grib2?var=Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average&var=u-component_of_wind_height_above_ground&var=v-component_of_wind_height_above_ground&var=Temperature_height_above_ground&var=Relative_humidity_height_above_ground&var=Per_cent_frozen_precipitation_surface&var=Precipitation_rate_surface&var=Total_cloud_cover_entire_atmosphere&north=47.5&west=280&east=293.25&south=40.25&horizStride=1&time_start={ts_date.isoformat()}Z&time_end={ts_date.isoformat()}Z&&&accept=netcdf4-classic'
			# Ending False is to now use a google bucket -- that's not an option for GFS
			download_list.append((url, fpath, False))
		else:
			print(f'Skipping download; {os.path.basename(fpath)} found at: {gfs_data_dir}')
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

def process_gfs_data(date,
					 location_dict,
					 gfs_data_dir,
					 num_threads):
	"""
	Will load NWM data and return a dictionary of timestamped streamflow data (pd.DataFrame) for each staion name in location_dict

	Args:
	-- date (datetime) [req]: datetime object representing the forecast date.
	-- location_dict (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to load meterological data for.
	-- gfs_data_dir (str) [req]: the directroy in which all of the GFS grib2 files are stored.
	-- num_threads (int) [req]: number of threds to use for reading grib2 files

	Returns:
	a dictionary of the GFS data in the format {station_ID:pd.Dataframe}
	"""

	station_dict = {}
	file_list = sorted(glob(f'{gfs_data_dir}/gfs.0p25.{date.strftime("%Y%m%d")}{date.hour:02d}.f[0-9][0-9][0-9].grib2.nc'))
	dataset_dict = multithreaded_loading(xr.open_dataset, file_list, num_threads)
	for fname, ds in dataset_dict.items():
		# making a timestamp for each file based on the filename. Doing this b/c the NetCDFs have no reliable (aka conisistently named) timestamp column.
		match = re.search(r'f(\d{3})', fname)
		hours_out = int(match.group(1))
		ts = (date + dt.timedelta(days=int(hours_out/24), hours=hours_out%24)).replace(tzinfo=dt.timezone.utc)
		# grabbing just the lat/long coords we need based on locations
		loc_ds = isolate_loc_rows(remap_longs(ds), location_dict)
		# selecting the hight levelss we need, and then dropping a bunch of junk
		loc_ds = loc_ds.sel({'height_above_ground2':10, 'height_above_ground3':2}).drop(['reftime','LatLon_721X1440-0p13S-180p00E','height_above_ground2','height_above_ground3','height_above_ground4'])
		# convert to a flat datafranme
		loc_df = loc_ds.to_dataframe().reset_index()
		# grab only the columns we need, except swdown b/c the4 name of that column is inconsistent
		filtered_df = loc_df.filter(items=['latitude','longitude','Relative_humidity_height_above_ground','Total_cloud_cover_entire_atmosphere','Precipitation_rate_surface','Temperature_height_above_ground','u-component_of_wind_height_above_ground','v-component_of_wind_height_above_ground','Per_cent_frozen_precipitation_surface'], axis=1)
		# rename columns appropriately
		filtered_df.rename({'Temperature_height_above_ground': 'T2', 'Total_cloud_cover_entire_atmosphere':'TCDC', 'u-component_of_wind_height_above_ground':'U10', 'v-component_of_wind_height_above_ground':'V10', 'Relative_humidity_height_above_ground':'RH2', 'Precipitation_rate_surface':'RAIN', 'Per_cent_frozen_precipitation_surface':'CPOFP'}, axis=1, inplace=True)
		# for .f000, set swdon to 0
		if fname.endswith('.f000.grib2.nc'):
			filtered_df['SWDOWN'] = 0
		# for all other timestamps, grab the data and add it to the df
		else:
			swdown = loc_df.filter(regex='Downward_Short-Wave_Radiation_Flux_surface*', axis=1)
			filtered_df['SWDOWN'] = swdown
		# group by lat/long
		location_groups = filtered_df.groupby(["latitude", "longitude"])
		location_dataframes = {name: group for name, group in location_groups}
		append_timestamp(sta_dict=station_dict, loc_dict=location_dict, loc_dfs=location_dataframes, timestamp=ts)
	return station_dict

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
	forecast_hours = generate_hours_list(get_hour_diff(forecast_datetime, end_datetime), 'gfs', archive=True)
	forecast_date = forecast_datetime.strftime("%Y%m%d")

	# define the directory for storing GFS data
	# directory should be made in download function, so that it can be used independently if necessary
	gfs_date_dir = os.path.join(data_dir, f'gfs/{forecast_datetime.strftime("%Y")}/gfs.{forecast_date}/{forecast_cycle}/atmos')

	# downloading the GFS data with multithreading; 
	download_gfs_threaded(date=forecast_datetime,
						  hours=forecast_hours,
						  gfs_data_dir=gfs_date_dir,
						  num_threads=dnwld_threads)
	
	# Now, load and process the data
	loc_data = process_gfs_data(date=forecast_datetime,
							    location_dict=locations,
								gfs_data_dir=gfs_date_dir,
								num_threads=load_threads)
	# ensure return_type is a valid value
	if return_type not in ['dict', 'dataframe']:
		raise ValueError(f"'{return_type}' is not a valid return_type. Please use 'dict' or 'dataframe'")
	elif return_type == 'dict':
		# created nested dictionary of pd.Series for each variable for each location
		gfs_data = {location:{name:data for name, data in loc_df.T.iterrows()} for location, loc_df in loc_data.items()}
	elif return_type == 'dataframe':
		raise Exception("'dataframe' option not implemented yet. Please use return_type = 'dict'")
	
	# return the GFS data
	return gfs_data