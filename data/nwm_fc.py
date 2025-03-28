import datetime as dt
from glob import glob
import data.utils as utils
import numpy as np
import os
import pandas as pd
import sh
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

	utils.download_data(Url, FilePath)
	
	return FilePath

def prepForDownloads(reference_date, member, download_dir):
	'''
	Prepares the directory structure, file name template, and reference date for downloading NWM forecast data, independent of the download method (GCS or NOMADs).
	
	Args:
	-- reference_date (str, date, or datetime) [req]: the launch date and time of the forecast to download. Time should specify model cycle (ex.'2015010112')
	-- member (str) [req]: The member type of NWM forecast to get (medium_range_mem1, long_range_mem3, short_range, etc).
	-- download_dir (str) [req]: directory to store donwloaded data.

	Returns:
	reference_date (datetime): The parsed reference datetime object
	netcdf_template (str): The file name template for the forecast files
	nwm_date_dir (str): The precise directory where the forecast files will be stored. Copies the NWM dir structure.
	'''
	# parse the reference date
	reference_date = utils.parse_to_datetime(reference_date)

	# break up the member string for file name and path construction
	mem_num = member.split('range')[-1]
	mem_type = member.split('_')[0]

	# define the forecast file name template that we want to look for - note that we are interested in channel_rt files as these contain streamflow data
	netcdf_template = f'nwm.t{reference_date.strftime("%H")}z.{mem_type}_range.channel_rt'

	# define the directory for storing NWM data. This directroy structure mimics the NWM GCS bucket structure
	nwm_date_dir = os.path.join(download_dir, f'nwm/{reference_date.strftime("%Y")}/nwm.{reference_date.strftime("%Y%m%d")}/{mem_type}_range{mem_num}')

	# print(f'TASK INITIATED: Download {hours} hours of {member} NWM hydrologic forecast for the following date and model cycle: {reference_date.strftime("%Y%m%d.t%Hz")}')
	# if not os.path.exists(nwm_date_dir):
	# 	os.makedirs(nwm_date_dir)

	return reference_date, netcdf_template, nwm_date_dir

def download_nwm(reference_date, member, hours='all', gcs=True, download_dir=tf.gettempdir(), num_threads=int(os.cpu_count()/2)):
	'''
	Downloads NWM forecast data from either the Google Cloud Storage (GCS) or NOMADS server. Designed to download one forecast product at a time
	NWM GCS Bucket: https://console.cloud.google.com/storage/browser/national-water-model
	NWM NOMADS server: https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/

	Args:
	-- reference_date (str, date, or datetime) [req]: the launch date and time of the forecast to download. Time should specify model cycle (ex.'2015010112')
	-- member (str) [req]: The member type of NWM forecast to get (medium_range_mem1, long_range_mem3, short_range, etc).
	-- hours (str or int) [opt]: the number of hours of forecast data to download. Default is 'all', which downloads all available files in the bucket.
	-- gcs (bool) [opt]: Flag determining wether or not to use google buckets for nwm download as opposed to NOMADs site. Default is True
	-- download_dir (str) [opt]: directory to store donwloaded data. Defaults to OS's default temp directory.
	-- num_threads (int) [opt]: number of threads to use for downloads. Default is half of OS's available threads.
	'''
	# prepare variables for NWM downloads, which includes 1) parsing the refernce date to a datetime object
	# 2) breaking up the member string for file name and path construction, and enables us to then 3) define the directory for storing NWM data
	reference_date, netcdf_template, nwm_date_dir = prepForDownloads(reference_date, member, download_dir)

	print(f'TASK INITIATED: Download {hours} hours of {member} NWM hydrologic forecast for the following date and model cycle: {reference_date.strftime("%Y%m%d.t%Hz")}')
	if not os.path.exists(nwm_date_dir):
		os.makedirs(nwm_date_dir)

	if gcs:
		forecast_page = os.path.join('gs://national-water-model/', f'nwm.{reference_date.strftime("%Y%m%d")}', member, f'{netcdf_template}*')
	else: forecast_page = os.path.join('https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/', f'nwm.{reference_date.strftime("%Y%m%d")}', f'{member}/')
	
	print(f'Acquiring data from: {forecast_page}')

	if gcs:
		# run ls command using gsutil to get the list of files in the bucket that match our file name template
		output = sh.gsutil(['ls', '-l', forecast_page])
		# Split output into lines
		lines = output.split("\n")
		# Extract All GCS file paths (last column in each line)
		all_server_paths = [line.split()[-1] for line in lines if line.strip() and not line.startswith("TOTAL:")]
	else:
		# scrape the forecast page for netcdf files based on the constructed filename template
		all_server_paths = utils.scrapeUrlsFromHttp(forecast_page, pattern=netcdf_template)


	# Filter the bucket paths to only include files that correspond with the number of hours of data requested
	if hours == 'all':
		# if all hours were requested, then just use all the files in the bucket
		bucket_paths = all_server_paths
	# otherwise filter by looking at the timestep number (fxxx) in the file name
	else: bucket_paths = [file for file in all_server_paths if int(file.split('.')[-3].split('f')[-1]) <= hours]

	# create local file paths where the downloaded data will be stored, mimicking the NWM GCS bucket structure
	file_paths = [os.path.join(nwm_date_dir, os.path.basename(path)) for path in bucket_paths]

	# zip the bucket paths and file paths together
	all_files = list(zip(bucket_paths, file_paths))

	download_list = []
	for bucket_file, local_file in all_files:
		# if the netcdf file isn't downloaded already, then download it
		if not os.path.exists(local_file):
			download_list.append((bucket_file, local_file, gcs))
		else:
			print(f'Skipping download; {os.path.basename(local_file)} found at: {local_file}')
	if download_list:
		utils.multithreaded_download(download_list, num_threads)
		pass
	print('TASK COMPLETE: NWM DOWNLOAD')
	# return a list of just the filepaths that were downloaded
	return [download_list[i][1] for i, _ in enumerate(download_list)]

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
		utils.multithreaded_download(download_list, num_threads)
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
	dataset_dict = utils.multithreaded_loading(xr.open_dataset, download_files, num_threads)

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

def get_data(start_date,
			 end_date,
			 locations,
			 forecast_type,
			 data_dir=tf.gettempdir(),
			 dwnld_threads=int(os.cpu_count()/2),
			 load_threads=int(os.cpu_count()/2),
			 google_buckets=False,
			 archive=False):
	"""
	A function to download and process NWM hydrology forecast data to return nested dictionary of pandas series fore each variable, for each location.
	
	Args:
	-- start_date (str, date, or datetime) [req]: the start date and time (00, 06, 12, 18) of the forecast to download. Times are assumed to be UTC time.
	-- end_date (str, date, or datetime) [req]: the end date and time for the forecast. GFS forecasts 16-days out for a given start date.
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to download forecast data for.
	-- forecast_type (str) [req]: The type of forecast (medium_range_mem1, long_range_mem3, short_range, etc).
	-- forecast_cycle (str) [req]: The starting time for the forecasts. valid values are 00, 06, 12, 18.
	-- data_dir (str) [opt]: directory to store donwloaded data. Defaults to OS's default temp directory.
	-- dwnld_threads (int) [opt]: number of threads to use for downloads. Default is half of OS's available threads.
	-- load_threads (int) [opt]: number of threads to use for reading data. Default is 2 for GFS, since file reads are already pretty fast.
	-- google_buckets (bool) [opt]: Flag determining wether or not to use google buckets for nwm download as opposed to NOMADs site.
	-- archive (bool) [opt]: Flag determining wether or not data you are grabbing is older than the last two days (relevant for NWM only)
	
	Returns:
	NWM timeseries for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	"""
	start_date = utils.parse_to_datetime(start_date)
	end_date = utils.parse_to_datetime(end_date)

	# the type of forecast to grab should be a parameter in the future, though right now the function can't accomodate all fc types produced by NWM
	# To get any kind of data at https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.20231218/ would probably require a new, comprehensive function
	forecast_cycle = f"{start_date.hour:02d}"
	forecast_member = forecast_type.split('range')[-1]
	forecast_type = forecast_type.split('_')[0]
	netcdf_template = f"nwm.t{forecast_cycle}z.{forecast_type}_range.channel_rt{forecast_member.replace('mem','')}"

	# determining what source we are downlaoding data from in order to create proper hours list
	if google_buckets:
		source = 'buckets'
		# Set archive = true if forecast older than 11/1/2020.  Might need to go even later... need to test
		if start_date < utils.parse_to_datetime('20201101'):
			archive = True
	else: source = 'nwm'

	# calculate the number of hours of forecast data to grab. I.e. for a 5 day forecast, hours would be 120
	forecast_hours = utils.generate_hours_list(utils.get_hour_diff(start_date, end_date), source, forecast_type, archive)
	forecast_date = start_date.strftime("%Y%m%d")

	# define the directory for storing NWM data
	# directory should be made in download function, so that it can be used independently if necessary
	nwm_date_dir = os.path.join(data_dir, f'nwm/{start_date.strftime("%Y")}/nwm.{forecast_date}/{forecast_type}_range{forecast_member}')

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
	
	nwm_data = {reach:{name:data for name, data in reach_df.T.iterrows()} for reach, reach_df in reach_data.items()}
	
	utils.add_units(nwm_data, {'streamflow':'m3 s-1'})

	# return the NWM data
	return nwm_data