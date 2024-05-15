import datetime as dt
from data.gfs_fc import isolate_loc_rows, remap_longs
from glob import glob
from .utils import (multithreaded_download,
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
		#print(f'shape before: {loc_dfs[loc].shape}')
		df_to_append = loc_dfs[loc].drop_duplicates().drop(columns = ['latitude','longitude'])
		#print(f'shape after: {df_to_append.shape}')
		df_to_append['time'] = timestamp
		df_to_append.set_index('time', inplace=True)
		if stationID in sta_dict:
			sta_dict[stationID] = pd.concat([sta_dict[stationID], df_to_append]).sort_index()
		else:
			sta_dict[stationID] = df_to_append

def download_gfs_threaded(date,
				 		  hours,
				 		  gfs_data_dir,
						  variables,
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
			
			# Create url string for 'normal' variables
			# Using += to build up URL: https://waymoot.org/home/python_string/
			reg_variable_string = ''
			reg_vars_dict = {x:y for x,y in variables.items() if not '_$3OR6' in y}
			for short_name, gfs_name in reg_vars_dict.items():
				reg_variable_string += '&var=' + gfs_name
			# Take off leading '&'
			if reg_variable_string != '':
				reg_variable_string = reg_variable_string[1:]
			# The variables that are 3/6 Hour averages need special handling
			ave_variable_string = ''
			ave_vars_dict = {x:y for x,y in variables.items() if '_$3OR6' in y}
			for short_name, gfs_name in ave_vars_dict.items():
				if h == '000':
					ave_variable_suffix_string = ''
				elif int(h) % 6 == 3:
					ave_variable_suffix_string = '_3_Hour_Average'
				elif int(h) % 6 == 0:
					ave_variable_suffix_string = '_6_Hour_Average'
				else:
					ave_variable_suffix_string = ''
				if ave_variable_suffix_string != '':
					ave_variable_string += '&var=' + gfs_name.split('_$3OR6')[0] + ave_variable_suffix_string
			url = f'https://thredds.rda.ucar.edu/thredds/ncss/grid/files/g/ds084.1/{date.year}/{date.strftime("%Y%m%d")}/gfs.0p25.{date.strftime("%Y%m%d")}{fc_cycle}.f{h}.grib2?{reg_variable_string}{ave_variable_string}&north=47.5&west=280&east=293.25&south=40.25&horizStride=1&time_start={ts_date.isoformat()}Z&time_end={ts_date.isoformat()}Z&&&accept=netcdf4-classic'

			# Ending False is to not use a google bucket -- that's not an option for GFS
			download_list.append((url, fpath, False))
		else:
			print(f'Skipping download; {os.path.basename(fpath)} found at: {gfs_data_dir}')
	# if download_list isn't empty:
	if download_list:
		# The true max hits/min allowed before getting blocked by NOMADS is 120/min, but they've requested users stick to 50/min
		max_hits_per_min = 65
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
					 variables_dict,
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
	# Checking the filelist for bad data. Most gfs netcdfs are around 200kb, so checking for any less thean 1 kb should suffice
	bad_files = []
	for file in file_list:
		# if a bad file is larger than 1000 bytes, it will get through and multithreaded_loading will fail)
		if os.path.getsize(file) < 1000:
			print(f"BAD FILE: {file}")
			print(f"FILE SIZE: {os.path.getsize(file)} BYTES")
			print("REMOVING FILE")
			# remove bad file from list; this is usually fine as AEM3D can interpolate over a missing timestep here and there
			bad_files.append(file)
	file_list = [f for f in file_list if f not in bad_files]
	# print(f"FILE LIST AFTER REMOVALS: {file_list}")
	dataset_dict = multithreaded_loading(xr.open_dataset, file_list, num_threads)
	for fname, ds in dataset_dict.items():
		# making a timestamp for each file based on the filename. Doing this b/c the NetCDFs have no reliable (aka conisistently named) timestamp column.
		match = re.search(r'f(\d{3})', fname)
		hours_out = int(match.group(1))
		ts = (date + dt.timedelta(days=int(hours_out/24), hours=hours_out%24)).replace(tzinfo=dt.timezone.utc)
		# grabbing just the lat/long coords we need based on locations
		loc_ds = isolate_loc_rows(remap_longs(ds), location_dict)
		# selecting the hight levelss we need, and then dropping a bunch of junk
		#print(loc_ds.sel({'height_above_ground2':10, 'height_above_ground3':2}).drop(['LatLon_721X1440-0p13S-180p00E','height_above_ground2','height_above_ground3']).to_dataframe())
		#loc_ds.to_dataframe().filter(regex='Total_cloud_cover_entire_atmosphere').to_csv(f"{fname}_TCDC.csv")
		#loc_ds.to_dataframe().filter(regex='Downward_Short-Wave_Radiation_Flux_surface').to_csv(f"{fname}_SWDOWN.csv")
		if 'bounds_dim' in loc_ds.dims.keys():
			#print('Selecting time_bounds dimension 0')
			loc_ds = loc_ds.sel({'bounds_dim':0})
		#print(f"Selecting on {loc_ds['Temperature_height_above_ground'].dims}, {loc_ds['u-component_of_wind_height_above_ground'].dims}")
		## TODO: do better on finding dims... should be [1], but something above is messing with them
		loc_ds = loc_ds.sel(
			{loc_ds['Temperature_height_above_ground'].dims[2]         :  2,
			 loc_ds['u-component_of_wind_height_above_ground'].dims[2] : 10,}).drop(
			['reftime','LatLon_721X1440-0p13S-180p00E',
	  		 loc_ds['Temperature_height_above_ground'].dims[2],
			 loc_ds['u-component_of_wind_height_above_ground'].dims[2],
			 loc_ds['Relative_humidity_height_above_ground'].dims[2]])
		# convert to a flat datafranme
		loc_df = loc_ds.to_dataframe().reset_index()
		# grab only the columns we need, except SWDOWN and TCDC b/c the4 name of that column is inconsistent
		filtered_df = loc_df.filter(items=['latitude','longitude','Relative_humidity_height_above_ground','Precipitation_rate_surface','Temperature_height_above_ground','u-component_of_wind_height_above_ground','v-component_of_wind_height_above_ground','Per_cent_frozen_precipitation_surface'], axis=1)
		# rename columns appropriately
		filtered_df.rename({'Temperature_height_above_ground': 'T2', 'u-component_of_wind_height_above_ground':'U10', 'v-component_of_wind_height_above_ground':'V10', 'Relative_humidity_height_above_ground':'RH2', 'Precipitation_rate_surface':'RAIN', 'Per_cent_frozen_precipitation_surface':'CPOFP'}, axis=1, inplace=True)
		# Now, grab SWDOWN and TCDC with their inconsistent names
		ave_vars_dict = {x:y for x,y in variables_dict.items() if '_$3OR6' in y}
		for short_name, gfs_name in ave_vars_dict.items():
			value = loc_df.filter(regex=gfs_name.split('_$3OR6')[0])
			if value.shape[1] == 0:
				filtered_df[short_name] = pd.NA
			else:
				filtered_df[short_name] = value
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
			 useTCDCInstant = False):
	"""
	Download specified GFS forecast data and return nested dictionary of pandas series fore each variable, for each location.

	Args:
	-- forecast_datetime (str, date, or datetime) [req]: the start date and time (00, 06, 12, 18) of the forecast to download. Times are assumed to be UTC time.
	-- end_datetime (str, date, or datetime) [req]: the end date and time for the forecast. GFS forecasts 16-days out for a given start date.
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to download forecast data for.
	-- data_dir (str) [opt]: directory to store donwloaded data. Defaults to OS's default temp directory.
	-- dwnld_threads (int) [opt]: number of threads to use for downloads. Default is half of OS's available threads.
	-- load_threads (int) [opt]: number of threads to use for reading data. Default is 2 for GFS, since file reads are already pretty fast.
	-- useTCDCInstant (bool) [opt]: wether to use instantaneous var for cloud cover or rolling average value.

	Returns:
	GFS forecast data for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	"""
	variables_dict = {
		'T2'    :'Temperature_height_above_ground',
		'TCDC'  :'Total_cloud_cover_entire_atmosphere_$3OR6_Hour_Average',
		'U10'   :'u-component_of_wind_height_above_ground',
		'V10'   :'v-component_of_wind_height_above_ground',
		'RH2'   :'Relative_humidity_height_above_ground',
		'RAIN'  :'Precipitation_rate_surface',
		'CPOFP' :'Per_cent_frozen_precipitation_surface',
		'SWDOWN':'Downward_Short-Wave_Radiation_Flux_surface_$3OR6_Hour_Average'
	}
		
	forecast_datetime = parse_to_datetime(forecast_datetime)
	end_datetime = parse_to_datetime(end_datetime)

	if useTCDCInstant:
		variables_dict['TCDC'] = 'Total_cloud_cover_entire_atmosphere'

	# grabbing the cycle (12am, 6am, 12pm, 6pm) from the datetime
	forecast_cycle = f"{forecast_datetime.hour:02d}"

	# we need the end_datetime to match the time of forecast_Datetime in order to calculate the number of forecast hours to download accurately
	# end_datetime = dt.datetime.combine(end_datetime.date(), forecast_datetime.time()).replace(tzinfo=dt.timezone.utc)
	# calculate the number of hours of forecast data to grab. I.e. for a 5 day forecast, hours would be 120
	print(f"forecast_datetime: {forecast_datetime} | end_datetime: {end_datetime} | hour_diff: {get_hour_diff(forecast_datetime, end_datetime)}")
	forecast_hours = generate_hours_list(get_hour_diff(forecast_datetime, end_datetime), 'gfs', archive=True)
	forecast_date = forecast_datetime.strftime("%Y%m%d")

	# define the directory for storing GFS data
	# directory should be made in download function, so that it can be used independently if necessary
	gfs_date_dir = os.path.join(data_dir, f'gfs/{forecast_datetime.strftime("%Y")}/gfs.{forecast_date}/{forecast_cycle}/atmos')

	# downloading the GFS data with multithreading; 
	download_gfs_threaded(date=forecast_datetime,
						  hours=forecast_hours,
						  gfs_data_dir=gfs_date_dir,
						  variables=variables_dict,
						  num_threads=dnwld_threads)
	
	# Now, load and process the data
	loc_data = process_gfs_data(date=forecast_datetime,
							    location_dict=locations,
								variables_dict=variables_dict,
								gfs_data_dir=gfs_date_dir,
								num_threads=load_threads)
	
	gfs_data = {location:{name:data.dropna().astype('float') for name, data in loc_df.T.iterrows()} for location, loc_df in loc_data.items()}
	
	# return the GFS data
	return gfs_data