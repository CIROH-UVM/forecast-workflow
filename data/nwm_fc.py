from glob import glob
import data.utils as utils
import numpy as np
import os
import re
import sh
import tempfile as tf
import xarray as xr

# Adapted from python notebooks at https://www.hydroshare.org/resource/5949aec47b484e689573beeb004a2917/

'''
This module provides functional tools to download and process NWM forecast product data from either NOMADs (which holds the last 2 days' NWM runs) or Google Cloud Services (which has NWM data back to 2018)

NOMADS: https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/
NWM GCS Bucket: https://console.cloud.google.com/storage/browser/national-water-model;tab=objects?invt=AbtP-w&prefix=&forceOnObjectsSortingFiltering=false

NOTE: This module provides support for aqcuiring any NWM forecast member (long_range_mem1, medium_range_mem2, short_range, etc) at any time cycle (00, 06, 12, 18). However, from these forecast products, only the `channel_rt` files are downloaded and processed.
This is because the 'channel_rt' files contain streamflow data, which is the most common type of data requested from the NWM in our use cases. Other data files, i.e. 'land' and/or 'reservoir' are not currently accessible with this version of the module. However, this module
could be modified in the future if desired to also aqcuire these files. Please reach out to the repository owners to request this feature.

Variables available in 'channel_rt' files:
	- "streamflow"
	- "velocity"
	- "nudge"
'''

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
	mem_type = member.split('_')[0]
	# define the forecast file name template that we want to look for - note that we are interested in channel_rt files as these contain streamflow data
	netcdf_template = f'nwm.t{reference_date.strftime("%H")}z.{mem_type}_range.channel_rt'
	# if the member is short_range, which has no members, there is no member number
	if member != 'short_range':
		# but if there is a member number (medium and long range), then include JUST an underscore plus the integer for the fname template (i.e '_1')
		netcdf_template = netcdf_template+'_'+member.split('range')[-1][-1]
	
	# define the directory for storing NWM data. This directroy structure mimics the NWM GCS bucket structure
	dir_structure = f'nwm/{reference_date.strftime("%Y")}/nwm.{reference_date.strftime("%Y%m%d")}/{mem_type}_range'
	if member != 'short_range':
		dir_structure = dir_structure+member.split('range')[-1]
	nwm_date_dir = os.path.join(download_dir, dir_structure)
	
	return reference_date, netcdf_template, nwm_date_dir

def download_nwm(reference_date, member, hours='all', gcs=True, download_dir=tf.gettempdir(), num_threads=int(os.cpu_count()/2)):
	'''
	Downloads NWM forecast data from either the Google Cloud Storage (GCS) or NOMADS server. Designed to download one forecast product at a time
	NWM GCS Bucket: https://console.cloud.google.com/storage/browser/national-water-model
	NWM NOMADS server: https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/

	Args:
	-- reference_date (str, date, or datetime) [req]: the launch date and time of the forecast to download. Time should specify model cycle (ex.'2015010112')
	-- member (str) [req]: The member type of NWM forecast to get (medium_range_mem1, long_range_mem3, short_range, etc).
	-- hours (str, int, or list) [opt]: the number of hours of forecast data to download. Default is 'all', which downloads all available files in the bucket. A list of integers
		indicating what hours to get is also acceptable (i.e. [12,14,16,...,20])
	-- gcs (bool) [opt]: Flag determining wether or not to use google buckets for nwm download as opposed to NOMADs site. Default is True
	-- download_dir (str) [opt]: directory to store downloaded data. Defaults to OS's default temp directory.
	-- num_threads (int) [opt]: number of threads to use for downloads. Default is half of OS's available threads.
	'''
	# prepare variables for NWM downloads, which includes 1) parsing the refernce date to a datetime object
	# 2) breaking up the member string for file name and path construction, and enables us to then 3) define the directory for storing NWM data
	reference_date, netcdf_template, nwm_date_dir = prepForDownloads(reference_date, member, download_dir)

	print(f'TASK INITIATED: Download {member} NWM hydrologic forecast for the following date and model cycle: {reference_date.strftime("%Y%m%d.t%Hz")}')
	if not os.path.exists(nwm_date_dir):
		os.makedirs(nwm_date_dir)

	# make te forecat page to scan for avaialable data, depending on the source (GCS vs NOMADS)
	if gcs:
		forecast_page = os.path.join('gs://national-water-model/', f'nwm.{reference_date.strftime("%Y%m%d")}', member)
	else: forecast_page = os.path.join('https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/', f'nwm.{reference_date.strftime("%Y%m%d")}', f'{member}/')
	
	print(f'Scanning: {forecast_page} for available data...')

	# Use gsutil to list all of the files available for download
	if gcs:
		# run ls command using gsutil to get the list of files in the bucket that match our file name template
		output = sh.gsutil(['ls', '-l', f'{forecast_page}/{netcdf_template}*'])
		# Split output into lines
		lines = output.split("\n")
		# Extract All GCS file paths (last column in each line)
		all_server_paths = [line.split()[-1] for line in lines if line.strip() and not line.startswith("TOTAL:")]
	# or if the source is NOMADs, scrap the links listed on the forecast product page using BeautifulSoup
	else:
		# scrape the forecast page for netcdf files based on the constructed filename template
		all_server_paths = utils.scrapeUrlsFromHttp(forecast_page, pattern=netcdf_template)


	# Filter the bucket paths to only include files that correspond with the number of hours of data requested
	if hours == 'all':
		# if all hours were requested, then just use all the files in the bucket
		bucket_paths = all_server_paths
	# if a list of hours to get was passed, then only filter those files that are in the hours list
	elif isinstance(hours, list):
		bucket_paths = [file for file in all_server_paths if int(file.split('.')[-3].split('f')[-1]) in hours]
	# if hours is passed as an int, get all the files up to and including hours
	elif isinstance(hours, int):
		bucket_paths = [file for file in all_server_paths if int(file.split('.')[-3].split('f')[-1]) <= hours]
	else: raise TypeError(f"'hours' argument must be an int, list of ints, or 'all':{hours}")

	# print(bucket_paths)
	# report what forecast timesteps are being selected
	first_ts = re.search(r'\.f(\d{3})\.', bucket_paths[0]).group(1)
	last_ts = re.search(r'\.f(\d{3})\.', bucket_paths[-1]).group(1)
	print(f"Seeking the following forecast timesteps: {first_ts} to {last_ts}")

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
		# execute multithreaded downloading
		utils.multithreaded_download(download_list, num_threads)
	print('TASK COMPLETE: NWM DOWNLOAD')
	# return a list of just the filepaths that were downloaded
	return [download_list[i][1] for i, _ in enumerate(download_list)]

def process_nwm(start_ts,
				end_ts,
				locations,
				variables,
				nwm_data_dir,
				fname_template,
				format='dictionary',
				end_date_exclusive=True,
				num_threads = int(os.cpu_count()/2)):
	'''
	Loads and process NWM data froma a single forecast product and returns the data in a format specified by the user. Data can be selected by start and end time, location, and variable.

	Args:
	-- start_ts (str, date, or datetime) [req]: The start date (and time) for which to slice the forecast data
	-- end_ts (str, date, or datetime) [req]: The end date (and time) for which to slice the forecast data
	-- locations (dict or list) [req]: a dictionary or list of locations to pull out of the forecast files. When a dict is passed in the format {"user-name":"gauge_ID"},
		then the function will use those user-defined names when returning data. Otherwise, default location IDs found in the dataset are used.
	-- variables (dict or list) [req]: a dictionary or list of variables to pull out of the 'channel_rt' forecast files. When a dictionary is passed in the format ("user-name":"variable-name"),
		then the function will use those user-defined names when returning data. Otherwise, default variable names found in the dataset are used. 
	-- nwm_data_dir (str) [req]: the directory in which the NWM forecast files to process are located
	-- fname_template (str) [req]: the generic filenmame for the specific netCDF's being loaded, up to the timestamp component of the file name, '.f###'. For instance,
		for a medium_range_mem1 forecast launched at the t06z model cycle, the valid fname_template would be "nwm.t06z.medium_range.channel_rt_1"
	-- format (str) [opt]: The format in which to return the data. Default is 'dictionary', which returns a nested dictionary of pandas series. Other valid option is 'xarray', which returns an xarray dataset.
	-- end_date_exclusive (bool) [opt]: Whether to exclude the ending timestamp from the time series. Defaults to True.
	-- num_threads (int) [req]: number of threds to use for reading netCDF files. Defauls to half of avaialbe CPU's to speed up processing

	Returns:
	a dictionary of the streamflow data in the format {reach_name:pd.Dataframe}
	'''
	# reconstructing the forecast member name based on the file name template
	# using the file name tempate rather than the nwm_data_dir becasue unlike the directory names, the file name template must have the relevant information about the forecast product
	fname_template_parts = fname_template.split('.')
	member = fname_template_parts[2]+'_mem'+fname_template_parts[-1][-1]
	cycle = fname_template_parts[1]

	print(f'TASK INITIATED: Process {member} {cycle} NWM hydrologic forecast data located at: {nwm_data_dir}')
	# print(fname_template)
	# ensure start and end ts are datetime objects
	start_ts = utils.parse_to_datetime(start_ts)
	end_ts = utils.parse_to_datetime(end_ts)
	
	# Get the filenames from download_base_path - in chronological order
	download_files = sorted(glob(os.path.join(nwm_data_dir, f'{fname_template}.f[0-9][0-9][0-9].conus.nc')))
	# print(download_files)
	
	# load the NWM data with multithreading
	# dataset_list will be in the same order as file list submitted, download_files
	print("Loading data...")
	dataset_dict = utils.multithreaded_loading(xr.open_dataset, download_files, num_threads)
	print("Loading complete.")

	# type checking for locations and variables. Must be a list
	if not isinstance(locations, (list, dict)):
		raise TypeError(f"Expected 'locations' to be a list or dictionary, but got type: {type(locations)}")
	if not isinstance(variables, (list, dict)):
		raise TypeError(f"Expected 'variables' to be a list or dictionary, but got type: {type(variables)}")

	# if locations (or variables) is a list, make it a dict with identical keys and values
	# We want to do this so users can pass either of these arguments as a list or dict
	locations = ({reach_id:reach_id for reach_id in locations} if isinstance(locations, list) else locations)
	variables = ({var:var for var in variables} if isinstance(variables, list) else variables)
	
	# create list of dataset-specic locations (i.e. the feature_ids) and variables (i.e. the variable names as defined in the datset)
	reaches_to_extract = [int(reach) for reach in locations.values()]
	vars_to_extract = [var for var in variables.values()]

	print(f"Extracting the following Gauge IDs (feature_id): {reaches_to_extract}")
	print(f"Extracing the following variables: {vars_to_extract}")

	print(f"Slicing data to get timestamps from {start_ts.strftime('%m-%d-%Y %H:%M:%S')} to {end_ts.strftime('%m-%d-%Y %H:%M:%S')}")
	timeslices = []

	# print(dataset_dict)
	# iterate through every netcdf file (i.e. forecast timeslice)
	for ds in dataset_dict.values():
		# we expect each file to be a singular timeslice of a forecast and so therfore should have exactly one timestamp
		if len(ds.time.values) != 1:
			raise ValueError(f"Multiple time steps found in dataset. Expected only one time step. Found {ds.time.values}")
		# if the timestamp comes before the start date, then skip it
		if ds.time.values[0] < np.datetime64(start_ts.replace(tzinfo=None)):
			continue
		# if end_date_exclusive is true, then exclude the end date as well
		if end_date_exclusive:
			if ds.time.values[0] >= np.datetime64(end_ts.replace(tzinfo=None)):
				continue
		elif ds.time.values[0] > np.datetime64(end_ts.replace(tzinfo=None)):
			continue
		# select reaches of interest
		ds = ds.sel(feature_id=reaches_to_extract)
		# list of variable names we want to drop
		drop_vars = list(set(ds.data_vars) - set(vars_to_extract))
		# drop variables that we don't want
		ds = ds.drop_vars(drop_vars)
		# print(f"Time step addded to concat list: {ds.time.values[0]}")
		timeslices.append(ds)

	# concatenate all of the timeslices into a single dataset
	ds = xr.concat(timeslices, dim='time')
	# localize timestamps to UTC
	# ds["time"] = ds["time"].to_index().tz_localize("UTC")

	# if 'locations' (or 'variables') was originally passed as a dictionary, then the keys shouldn't match the values (i.e. the user is giving us a mapping)
	# if every key DOES match the corresponding value, then the argument was probably passed as a list
	# we only create a 'locations' coord and rename variables IF a mapping (i.e. a dict) was originally passsed
	if any([k != v for k,v in locations.items()]):
		# if the user provided a dictionary of custom location names to use, then we will create a "locations" coordinate
		ds = ds.assign_coords(location=("feature_id", list(locations.keys())))
		# and then we swap with "feature_id" so that the "locations" are the primary dimension after "time"
		ds = ds.swap_dims({"feature_id":"location"})
		# assign an attribute to explain what the "location" coordinate is
		ds.coords['location'].attrs.update({"description":"user-defined location names", "comment":"location names correspond to feature_id's, respectively"})
	if any([k != v for k,v in variables.items()]):
		# if the user provided a dictionary of custom variable names to use, then ;et's rename the vars in the dataset
		ds = ds.rename({v:k for k,v in variables.items()})

	# returned the filtered xarray dataset if that is the format requested
	if format == 'xarray':
		print('TASK COMPLETE: PROCESS NWM')
		return ds
	elif format == 'dictionary':
		# create the classic nested dictionary of pandas series if that is the format requested
		nwm_data = {}
		for reach_name, reach_id in locations.items():
			# again, this if statement serves as our check as to whether or not 'locations' was originally passed as a dict;
			# if it was, then a location coordinate was created and we need to access that one, not feature_id
			if any([k != v for k,v in locations.items()]):
				loc_ds = ds.sel(location=reach_name)
			else: loc_ds = ds.sel(feature_id=int(reach_id))
			nwm_data[reach_name] = {}
			# we can just use the keys of 'variables' at this point because 'ds' variable names should match 'var_names'
			# whether the data variables or not
			for var_name in variables.keys():
				# Extract series and localize to UTC
				nwm_data[reach_name][var_name] = loc_ds[var_name].to_pandas().tz_localize("UTC")

		# Add units safely; i.e. if "units" attribute does not exist, just use the variable name as the series' name
		var_units = {var_name: getattr(ds.data_vars[var_name], "units", var_name) for var_name in variables.keys()}
		utils.add_units(nwm_data, var_units=var_units)
		print('TASK COMPLETE: PROCESS NWM')
		return nwm_data

def get_data(start_date,
			 end_date,
			 member,
			 locations,
			 variables,
			 reference_date=None,
			 data_dir=tf.gettempdir(),
			 format='dictionary',
			 gcs=True,
			 end_date_exclusive=True,
			 dwnld_threads=int(os.cpu_count()/2),
			 load_threads=int(os.cpu_count()/2)):
	'''
	A function to download and process NWM hydrology forecast data to return nested dictionary of pandas series fore each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: The start date for which to retrieve data.
	-- end_date (str, date, or datetime) [req]: The end date for which to retrieve data.
	-- member (str) [req]: The member type of NWM forecast to get (medium_range_mem1, long_range_mem3, short_range, etc).
	-- locations (dict) [req]: A dictionary mapping location names (or station IDs) to latitude-longitude tuples.
	-- variables (dict) [req]: A dictionary mapping user-defined variable names to dataset-specific variable names.
	-- reference_date (str, date, or datetime) [opt]: The forecast reference time, i.e., the date and time at which 
		the forecast was initialized. Defaults to start_date if None.
	-- data_dir (str) [opt]: directory to store donwloaded data. Defaults to OS's default temp directory.
	-- format (str) [opt]: The format in which to return the data. Default is 'dictionary', which returns a nested dictionary of pandas series. Other valid option is 'xarray', which returns an xarray dataset.
	-- gcs (bool) [opt]: Flag determining wether or not to use google buckets for nwm download as opposed to NOMADs site. Default is True
	-- end_date_exclusive (bool) [opt]: Whether to exclude the end date from the time series. Defaults to True.
	-- dwnld_threads (int) [opt]: number of threads to use for downloads. Default is half of OS's available threads.
	-- load_threads (int) [opt]: number of threads to use for reading data. Default is 2 for GFS, since file reads are already pretty fast.

	Returns:
	NWM timeseries for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	'''
	# Validate and process forecast dates
	start_date, end_date, reference_date = utils.validate_forecast_times(start_date, end_date, reference_date)

	reference_date, netcdf_template, nwm_date_dir = prepForDownloads(reference_date, member, data_dir)
	# print(netcdf_template)

	# Calculate which forecast timesteps to download based on reference_date, start_date, and end_date
	# this allows for a precise selection of data from any given forecast
	ref_to_end_hours = int((end_date - reference_date).total_seconds() / 3600)
	ref_to_start_hours = int((start_date - reference_date).total_seconds() / 3600)
	num_hours = list(range(ref_to_start_hours, ref_to_end_hours+1))

	# now download the data
	downloaded_list = download_nwm(reference_date, member, hours=num_hours, gcs=gcs, download_dir=data_dir, num_threads=dwnld_threads)

	# process the data
	nwm_data = process_nwm(start_date, end_date, locations, variables, nwm_date_dir, netcdf_template, format=format, end_date_exclusive=end_date_exclusive, num_threads=load_threads)

	return nwm_data