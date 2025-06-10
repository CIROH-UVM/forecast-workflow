from glob import glob
import datetime as dt
import data.utils as utils
import numpy as np
import os
import re
import sh
import tempfile as tf
import xarray as xr

'''
This module provides functional tools to download and process NWM forecast forcings data from the NWM Google Cloud Storage (GCS) bucket.

NWM GCS Bucket: https://console.cloud.google.com/storage/browser/national-water-model

NOTE: This module provides support for aqcuiring medium and short range NWM forecast forcings data only. Forcings for other forecast projects
are not supported at this time but may be added by request

For avaiable forcings variables, see the NWM_FORCING_VARS dictionary below.
'''

# all variables included in the NWM forecast forcings files
# keys are the 'long_name' attributes of the variables as described in the netCDF files
# values are the actual data variable names used in the NWM forcing netCDF files
NWM_FORCING_VARS = {'10-m U-component of wind':'U2D',
					'10-m V-component of wind':'V2D',
					'Surface downward long-wave radiation flux':'LWDOWN',
					'Surface Precipitation Rate':'RAINRATE',
					'2-m Air Temperature':'T2D',
					'2-m Specific Humidity':'Q2D',
					'Surface Pressure':'PSFC',
					'Surface downward short-wave radiation flux':'SWDOWN'}

def prepForDownloads(
	reference_date: str | dt.date | dt.datetime,
	member: str,
	download_dir: str
) -> tuple[dt.datetime, str, str]:
	'''
	Prepares the directory structure, file name template, and reference date for downloading NWM forcing data from the NWM Google Bucket:
	https://console.cloud.google.com/storage/browser/national-water-model

	Args:
	-- reference_date: The launch date and time of the forecast that you want to download forcings for. Time should specify model cycle (e.g., '2015010112').
	-- member: The NWM forecast member for which you to get forcings for (Currently accepts 'forcing_medium_range', 'forcing_short_range').
	-- download_dir: Directory to store downloaded data.

	Returns:
	tuple:
		- reference_date: The parsed reference datetime object.
		- netcdf_template: The file name template for the forecast forcings files.
		- nwm_date_dir: The precise directory where the forecast forcings will be stored. Mimics the NWM Google Bucket directory structure.
	'''
	reference_date = utils.parse_to_datetime(reference_date)

	# define the forecast file name template that we want to look for - note that we are interested in channel_rt files as these contain streamflow data
	netcdf_template = f'nwm.t{reference_date.strftime("%H")}z.{member}.forcing'

	# define the directory for storing NWM data. This directory structure mimics the NWM GCS bucket structure
	dir_structure = f'nwm/{reference_date.strftime("%Y")}/nwm.{reference_date.strftime("%Y%m%d")}/forcing_{member}'

	nwm_date_dir = os.path.join(download_dir, dir_structure)

	return reference_date, netcdf_template, nwm_date_dir

def download_nwm_forcings(
	reference_date: str | dt.date | dt.datetime,
	member: str,
	hours: str | int | list[int] = 'all',
	download_dir: str = tf.gettempdir(),
	num_threads: int = int(os.cpu_count() / 2)
) -> list[str]:
	'''
	Downloads NWM forecast forcings data from the Google Cloud Storage (GCS). Designed to download one forecast product at a time
	NWM GCS Bucket: https://console.cloud.google.com/storage/browser/national-water-model

	Args:
	-- reference_date: The launch date and time of the forecast to get forcings for. Time should specify model cycle (e.g., '2015010112').
	-- member: The NWM forecast member for which you to get forcings for (Currently accepts 'medium_range', 'short_range').
	-- hours: The number of hours of forcing data to download. Default is 'all', which downloads all available files in the bucket. A list of integers indicating what hours to get is also acceptable (e.g., [12, 14, 16, ..., 20]).
	-- download_dir: Directory to store downloaded data. Defaults to OS's default temp directory.
	-- num_threads: Number of threads to use for downloads. Default is half of OS's available threads.

	Returns:
	list: A list of file paths for the downloaded files.
	'''
	# prepare variables for NWM forcings downloads, which includes 1) parsing the refernce date to a datetime object
	# 2) breaking up the member string for file name and path construction, and enables us to then 3) define the directory for storing NWM data
	reference_date, netcdf_template, nwm_date_dir = prepForDownloads(reference_date, member, download_dir)

	print(f'TASK INITIATED: Download {member} NWM forcings for the following date and model cycle: {reference_date.strftime("%Y%m%d.t%Hz")}')
	if not os.path.exists(nwm_date_dir):
		os.makedirs(nwm_date_dir)

	# make the forecat page to scan for avaialable data
	forecast_page = os.path.join('gs://national-water-model/', f'nwm.{reference_date.strftime("%Y%m%d")}', f'forcing_{member}')

	print(f'Scanning: {forecast_page} for available data...')

	# Use gsutil to list all of the files available for download
	# run ls command using gsutil to get the list of files in the bucket that match our file name template
	output = sh.gsutil(['ls', '-l', f'{forecast_page}/{netcdf_template}*'])
	# Split output into lines
	lines = output.split("\n")
	# Extract All GCS file paths (last column in each line)
	all_server_paths = [line.split()[-1] for line in lines if line.strip() and not line.startswith("TOTAL:")]

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
	print(f"Seeking the following forcings timesteps: {first_ts} to {last_ts}")

	# create local file paths where the downloaded data will be stored, mimicking the NWM GCS bucket structure
	file_paths = [os.path.join(nwm_date_dir, os.path.basename(path)) for path in bucket_paths]
	# zip the bucket paths and file paths together
	all_files = list(zip(bucket_paths, file_paths))

	download_list = []
	for bucket_file, local_file in all_files:
		# if the netcdf file isn't downloaded already, then download it
		if not os.path.exists(local_file):
			download_list.append((bucket_file, local_file, True))
		else:
			print(f'Skipping download; {os.path.basename(local_file)} found at: {local_file}')
	if download_list:
		# execute multithreaded downloading
		utils.multithreaded_download(download_list, num_threads)

	print('TASK COMPLETE: NWM FORCINGS DOWNLOAD')
	# return a list of just the filepaths that were downloaded
	return [download_list[i][1] for i, _ in enumerate(download_list)]