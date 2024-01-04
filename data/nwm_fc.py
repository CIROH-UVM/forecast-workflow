import datetime as dt
from glob import glob
from lib import (download_data,
				 multithreaded_download,
				 multithreaded_loading,
				 parse_to_datetime,
				 get_hour_diff,
				 generate_hours_list)
import numpy as np
import os
import pandas as pd
import tempfile as tf
import xarray as xr

# Adapted from python notebooks at https://www.hydroshare.org/resource/5949aec47b484e689573beeb004a2917/

def get_forecast_File(Url, download_dir='.'):
	"""	
	A Function to download the given forecast url from NOAA NWM.

	Args:
	-- Url (str) [req]: The URL of the forecast file to be downloaded, should be string.
	-- download_dir (str) [opt]: The directory where the file should be saved. 

	Returns:
	A String Path of the Downloaded File. 
	
	"""
	FileName = os.path.basename(Url)
	# Lets construct the complete file path
	FilePath = os.path.join(download_dir, FileName)
	# Lets create download_path if it doesn't exist yet
	if not os.path.exists(download_dir):
		os.makedirs(download_dir)    
	# Lets make sure the file is not already downloaded - If yes, we will not download it again.
	if os.path.exists(FilePath):
		print(f'Skipping download: "{FilePath}" already exists')
		return None
	   
	print(f'Downloading {FileName}')

	# Lets do the download
	# r = requests.get(Url, allow_redirects=True)
	# # Time to save the file to the download_dir. 
	# open(FilePath, 'wb').write(r.content)
	
	# sh.curl('-o', FilePath,'-C','-', Url)

	download_data(Url, FilePath)
	
	return FilePath

# Next we will define a function which will call these two function and download the data. 
def download_nwm_threaded(date,
						  hours,
						  nwm_data_dir,
						  num_threads,
						  fname_template,
						  use_google_bucket=False):
	"""
	Downloads the NWM medium range mem1 CONUS forecast NetCDF files for the specified date and hours

	!!!  NOMADS has a 120 hits / minute limit that this code will quickly overwhelm.  We need to build in a rate limiter.
	TODO: Ideas for rate limiting code?  Maybe break the file list into 120 file sections and submit those to the
	      parallel download code once every 1 minute or so?

	Args:
	-- date (str) [req]: date to download gribs for. Default is just the current date
	-- hours (list of strs) [req]: list of forecast hours to download gribs for. Default is 7 day forecast, or 168 hours
	-- nwm_data_dir (str) [req]: directory in which to store nwm netcdf files
	-- num_threads (int) [req]: number of threads to use.
	-- fname_template (str) [req]: the filename template for the NWM forecast; includes info about forecat cycle, type, and member
	-- use_google_bucket (bool) [opt]: flag determining wheteher or not to use google buckets for nwm download as opposed to NOMADs site
	"""
	# Lets create an empty list to store the complete path of downloaded files
	download_list = []
	print(f'TASK INITIATED: Download {int(hours[-1])}-hour NWM hydrology forecasts for the following date: {date}')
	if not os.path.exists(nwm_data_dir):
		os.makedirs(nwm_data_dir)
	for h in hours:
		netcdf_fname = f'{fname_template}.f{h}.conus.nc'
		netcdf_fpath = os.path.join(nwm_data_dir, netcdf_fname)
		# if the netcdf file isn't downloaded already, then download it
		if not os.path.exists(netcdf_fpath):
			# Getting the URL
			if use_google_bucket:
				base_name = 'gs://national-water-model/nwm.'
			else:
				base_name = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.'
			netcdf_url = f'{base_name}{nwm_data_dir.split("nwm.")[-1]}/{netcdf_fname}'
			download_list.append((netcdf_url, netcdf_fpath, use_google_bucket))
		else:
			print(f'Skipping download; {os.path.basename(netcdf_fpath)} found at: {netcdf_fpath}')
	if download_list:
		multithreaded_download(download_list, num_threads)
	print('TASK COMPLETE: NWM DOWNLOAD')

def process_nwm_data(location_dict,
					 nwm_data_dir,
					 num_threads,
					 fname_template):
	"""
	Will load NWM data and return a dictionary of timestamped streamflow data (pd.DataFrame) for each staion name in location_dict

	Args:
	-- location_dict (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to load streamflow data for.
	-- nwm_data_dir (str) [req]: the directroy in which all of the NWM netCDF files are stored.
	-- num_threads (int) [req]: number of threds to use for reading netCDF files
	-- fname_template (str) [req]: the generic filenmame for the specific netCDF's being loaded, up to the timestamp component of the file name, '.f###'.

	Returns:
	a dictionary of the streamflow data in the format {reach_name:pd.Dataframe}
	"""
	
	# Get the filenames from download_base_path - in chronological order
	download_files = sorted(glob(f'{nwm_data_dir}/{fname_template}.f[0-9][0-9][0-9].conus.nc'))
	
	reach_dict = {reach_name:{'time':[], 'streamflow':[]} for reach_name, reach_id in location_dict.items()}
   
	# load the NWM data with multithreading
	# dataset_list will be in the same order as file list submitted, download_files
	dataset_dict = multithreaded_loading(xr.open_dataset, download_files, num_threads)

	# Time to read the data files
	for data in dataset_dict.values():
		# Time to extract the relevent streamflow and timestamp data from the dataset
		for reach_name, reach_id in location_dict.items():
			stream_value = np.float64(data.sel(feature_id=int(reach_id)).streamflow.values)
			timestamp = data.coords['time'].values[0]
			# Time to append it in results dict
			reach_dict[reach_name]['time'].append(timestamp)
			reach_dict[reach_name]['streamflow'].append(stream_value)

	# now get the data in the proper seies format
	results_dict = {}
	for reach_name, reach_dict in reach_dict.items():
		reach_df  = pd.DataFrame(data=reach_dict['streamflow'], index=reach_dict['time'], columns=['streamflow'])
		reach_df = reach_df.tz_localize('UTC')
		reach_df.index.name = 'time'
		results_dict[reach_name] = reach_df
	
	return results_dict

def get_data(forecast_datetime,
			 end_datetime,
			 locations,
			 data_dir=tf.gettempdir(),
			 dwnld_threads=int(os.cpu_count()/2),
			 load_threads=int(os.cpu_count()/2),
			 google_buckets=False,
			 return_type='dict'):
	"""

	A function to download and process NWM hydrology forecast data to return nested dictionary of pandas series fore each variable, for each location.
	
	Args:
	-- forecast_datetime (str, date, or datetime) [req]: the start date and time (00, 06, 12, 18) of the forecast to download. Times are assumed to be UTC time.
	-- end_datetime (str, date, or datetime) [req]: the end date and time for the forecast. GFS forecasts 16-days out for a given start date.
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to download forecast data for.
	-- data_dir (str) [opt]: directory to store donwloaded data. Defaults to OS's default temp directory.
	-- dwnld_threads (int) [opt]: number of threads to use for downloads. Default is half of OS's available threads.
	-- load_threads (int) [opt]: number of threads to use for reading data. Default is 2 for GFS, since file reads are already pretty fast.
	-- forecast_cycle (str) [req]: The starting time for the forecasts. valid values are 00, 06, 12, 18
	-- forecast_type (str) [req]: The type of forecast (default is 'medium_range').
	-- forecast_member (str) [req]: The member of the forecast model.
	-- google_buckets (bool) [opt]: Flag determining wheteher or not to use google buckets for nwm download as opposed to NOMADs site.
	-- return_type (string) [opt]: string indicating which format to return data in. Default is "dict", which will return data in a nested dict format:
									{locationID1:{
										var1_name:pd.Series,
										var2_name:pd.Series,
										...},
									locationID2:{...},
									...
									}
									Alternative return type is "dataframe", which smashes all data into a single dataframe muliIndex'd by station ID, then timestamp
	save_csv              : Flag indicating whether or not dataframes should be saved as CSV files.
	forecast_member        : The member of the forecast model.
	
	Returns:
	NWM data in the format specified by return_type
	"""
	forecast_datetime = parse_to_datetime(forecast_datetime)
	end_datetime = parse_to_datetime(end_datetime)

	# the type of forecast to grab should be a parameter in the future, though right now the function can't accomodate all fc types produced by NWM
	# To get any kind of data at https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.20231218/ would probably require a new, comprehensive function
	forecast_cycle = f"{forecast_datetime.hour:02d}"
	forecast_type = 'medium_range'
	forecast_member = '1'
	netcdf_template = f'nwm.t{forecast_cycle}z.{forecast_type}.channel_rt_{forecast_member}'


	# we need the end_datetime to match the time of forecast_datetime in order to calculate the number of forecast hours to download accurately
	end_datetime = dt.datetime.combine(end_datetime.date(), forecast_datetime.time())
	# calculate the number of hours of forecast data to grab. I.e. for a 5 day forecast, hours would be 120
	forecast_hours = generate_hours_list(get_hour_diff(forecast_datetime, end_datetime), 'nwm')
	forecast_date = forecast_datetime.strftime("%Y%m%d")
	
	# define the directory for storing NWM data
	# directory should be made in download function, so that it can be used independently if necessary
	nwm_date_dir = os.path.join(data_dir, f'nwm/nwm.{forecast_date}/{forecast_type}_mem{forecast_member}')

	# make the CSV directory - this may be outdated/unneeded at this point
	nwm_csv_dir = f'nwm/nwm.{forecast_date}'
	if not os.path.exists(nwm_csv_dir):
		os.makedirs(nwm_csv_dir)

	# downloading the NWM data with multithreading; 
	download_nwm_threaded(date=forecast_date,
						  hours=forecast_hours,
						  nwm_data_dir=nwm_date_dir,
						  num_threads=dwnld_threads,
						  fname_template=netcdf_template,
						  use_google_bucket=google_buckets)

	# Now, load and process the data
	reach_data = process_nwm_data(location_dict=locations,
							   	  nwm_data_dir=nwm_date_dir,
								  num_threads=load_threads,
								  fname_template=netcdf_template)
	
	# ensure return_type is a valid value
	if return_type not in ['dict', 'dataframe']:
		raise ValueError(f"'{return_type}' is not a valid return_type. Please use 'dict' or 'dataframe'")
	elif return_type == 'dict':
		# created nested dictionary of pd.Series for each variable for each location
		nwm_data = {reach:{name:data for name, data in reach_df.T.iterrows()} for reach, reach_df in reach_data.items()}
	elif return_type == 'dataframe':
		raise Exception("'dataframe' option not implemented yet. Please use return_type = 'dict'")
	
	# return the NWM data
	return nwm_data