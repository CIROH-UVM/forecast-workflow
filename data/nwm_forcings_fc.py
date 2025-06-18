import datetime as dt
import data.utils as utils
import numpy as np
import os
import pyproj
import re
import rioxarray
import sh
import tempfile as tf
import xarray as xr
from glob import glob

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

def subset_by_bbox(ds: xr.Dataset, bbox: dict) -> xr.Dataset:
	'''
	Subsets the NWM forcings xarray dataset by a bounding box.

	Args:
	-- ds: The xarray dataset to subset. 
	-- bbox: A dictionary containing the bounding box with keys 'min_lat', 'max_lat', 'min_lon', 'max_lon'. CRS of bounding box is assumed to be WGS84 (EPSG:4326).

	Returns:
	xr.Dataset: The subsetted dataset.
	'''
	# check that the bounding box has the required keys
	required_keys = ['min_lat', 'max_lat', 'min_lon', 'max_lon']
	missing_keys = [key for key in required_keys if key not in bbox]
	if missing_keys:
		raise ValueError(f"Bounding box is missing the following key(s): {missing_keys}. Please check your input.")
	
	ref_time = ds.coords['reference_time']
	# need to specify to rioxarray that the bounding box is in WGS84 coordinates
	wgs84 = rioxarray.crs.CRS.from_string("EPSG:4326")
	# the crs arg tells rioxarray that the bounding box coordinates are in WGS84
	# rioxarray infers the CRS of the dataset and reprojects the bounding box to the dataset's CRS
	clipped_ds = ds.rio.clip_box(minx=bbox['min_lon'], miny=bbox['min_lat'], maxx=bbox['max_lon'], maxy=bbox['max_lat'], crs=wgs84)
	
	# asssign the reference time to the clipped dataset
	# the reference time is dropped by the clip_box function because rioxarray drops non-index, non-spatial, non-dimension coordinates
	# The reason is that under the hood, rioxarray/rasterio performs operations like masking/clipping that do not know or care about additional coordinates that are not tied to raster grid layout.
	clipped_ds = clipped_ds.assign_coords(reference_time=ref_time)

	return clipped_ds

def subset_by_points(ds: xr.Dataset, points: list[tuple[float, float]]) -> xr.Dataset:
	'''
	Filters the NWM forcings dataset to only include the geospatial coordinates requested

	Args:
	-- ds: The dataset to subset.
	-- points: List of (latitude, longitude) pairs in WGS84 (EPSG:4326) coordinates.

	Returns:
	xr.Dataset: Dataset containing data at the nearest grid cell to each requested point.
		Note: Exact coordinate matches are not guaranteed; the nearest available grid points are used.
	'''
	# pull the CRS WKT string from the datset
	nwm_forcings_crs_wkt = ds['crs'].attrs['crs_wkt']
	# Create a Transformer to reproject from the NWM forcings CRS to WGS 1984
	wgs_to_nwm = pyproj.Transformer.from_crs('EPSG:4326', nwm_forcings_crs_wkt, always_xy=True)
	# Transform points - note that the resultant lcc_corners tuples will be in the order (x, y) (derived from (lon, lat))
	proj_points = [wgs_to_nwm.transform(lon, lat) for (lat, lon) in points]
	x, y = [pair[0] for pair in proj_points], [pair[1] for pair in proj_points]
	# now select the x and y coordinates from the dataset
	return ds.sel(x=x, y=y, method='nearest')

def validate_locations(locations: dict | None) -> str | None:
	'''
	Validates the locations argument to ensure it is either None or a dictionary with the correct keys.

	Args:
	-- locations: A nested dictionary containing a single 1st-level key, either 'bbox' or 'points' for subsetting the data (see examples below). Or None if no subsetting is required.
		{'points' : [(45.21, -75.89), (42.34, -72.68)]}
		{'bbox' : {'min_lat' : 42.34, 'max_lat' : 45.21, 'min_lon' : -75.89, 'max_lon' : -72.68}}

	Raises:
	TypeError: If locations is not None or a dictionary.
	TypeError: If 'points' is not a list of tuples.
	ValueError: If both 'bbox' and 'points' are specified in the dictionary.

	Returns:
	str | None: Returns 'bbox' if a bounding box is specified, 'points' if a list of points is specified, or None if no subsetting is required.
	'''
	# type checking for locations. locations must be None or a dict, with a single key, either 'bbox' for bounding box or 'points' for a list of points.
	if locations is None:
		# if locations is None, then we will not subset the data by location
		return None

	# if locations is not None, then it must be a dictionary with either 'bbox' or 'points' as keys
	elif not isinstance(locations,  (dict)):
		raise TypeError(f"Expected 'locations' to be a dictionary or None, but got type: {type(locations)}")

	#  we don't want to allow both 'bbox' and 'points' to be specified at the same time, so we check for that
	elif 'bbox' in locations and 'points' in locations:
		raise ValueError("Only one of 'bbox' or 'points' can be specified in 'locations'. Please choose one to subset forcing data with.")

	# if a bounding box is specified, check that it has the required keys and subset
	elif 'bbox' in locations:
		# subset the data by the bounding box
		# ds = subset_by_bbox(ds, locations['bbox'])
		return 'bbox'
	# if a list of points is specified, check that it is a list of tuples and subset
	elif 'points' in locations:
		if not isinstance(locations['points'], list) or not all(isinstance(point, tuple) and len(point) == 2 for point in locations['points']):
			raise TypeError("'points' argument must be a list of tuples with (lat, lon) pairs.")
		# subset the data by the points
		# ds = subset_by_points(ds, locations['points'])
		return 'points'
	
def parse_variables(variables: dict | list) -> dict:
	'''
	Parses and validate the variables argument.

	Args:
	-- variables: A dictionary or list of variables to extract from the dataset. If a list is passed, it will be converted to a dictionary with identical keys and values.

	Raises:
	ValueError: If the variables argument is a string that is not 'all'.
	TypeError: If the variables argument is not a dictionary, list, or 'all'.

	Returns:
	dict: A dictionary where the values are always the dataset-specifc variable names. Keys are either the same as values (if a list or 'all' was passed) or user-defined names (if a dictionary was passed).
	'''
	if isinstance(variables, str):
		if variables.lower() == 'all':
			# if the user wants all variables, then use the NWM_FORCING_VARS dictionary to get all of the variables
			variables = list(NWM_FORCING_VARS.values())
		else:
			raise ValueError(f"Invalid string value for 'variables': {variables}.")
	elif not isinstance(variables, (dict, list)):
		raise TypeError(f"Invalid type for 'variables': {type(variables)}. Expected dict, list, or 'all'.")

	# if variables is a list, convert it to a dictionary with identical keys and values
	variables = ({var:var for var in variables} if isinstance(variables, list) else variables)

	return variables

def process_nwm_forcings(
	start_ts: str | dt.date | dt.datetime,
	end_ts: str | dt.date | dt.datetime,
	nwm_data_dir: str,
	fname_template: str,
	locations: dict | None = None,
	variables: dict | list | str = 'all',
	end_date_exclusive: bool = True,
	num_threads: int = int(os.cpu_count() / 2)
) -> dict | xr.Dataset:
	'''
	Loads and processes NWM forcing data. Slices the dataset by location, time, and variables of interest.

	Args:
	-- start_ts: The start date (and time) for which to slice the forecast data.
	-- end_ts: The end date (and time) for which to slice the forecast data.
	-- nwm_data_dir: The directory in which the NWM forcing files to process are located.
	-- fname_template: The generic filename for the specific netCDF files being loaded, up to the timestamp component of the file name, '.f###'. For instance, for forcing_medium_range, the valid fname_template would be "nwm.t00z.medium_range.forcing".
	-- locations: A dictionary containing either bounding box information or a list of points to extract from the gridded forcings dataset. Default value of None does not spatially subset the data. See validate_locations() for more details.
		Note that the bounding box must be in latitude, longitude (WGS 1984), but the dataset is NOT reprojected and will maintain its orginal CRS.
	-- variables: A dictionary or list of variables to pull out of the forcing files. When a dictionary is passed in the format {"user-name":"variable-name"}, the function will use those user-defined names when returning data. Otherwise, the variable names found in the dataset are used. Default value 'all' keeps all variables
	-- end_date_exclusive: Whether to exclude the ending timestamp from the time series. Defaults to True.
	-- num_threads: Number of threads to use for reading netCDF files. Defaults to half of available CPUs to speed up processing.

	Returns:
	xr.Dataset or dict: If locations is None or a bounding box, returns an xarray.Dataset. If locations specifies points, returns a nested dictionary where 1st-level keys are points, 2nd-level keys are variables, and 2nd-level values are pandas.Series
	'''
	# reconstructing the forcing member name based on the file name template
	# using the file name tempate rather than the nwm_data_dir becasue unlike the directory names, the file name template must have the relevant information about the forcings product
	fname_template_parts = fname_template.split('.')
	member = fname_template_parts[-2]
	cycle = fname_template_parts[1]

	print(f'TASK INITIATED: Process {member} {cycle} NWM forecast forcing data located at: {nwm_data_dir}')

	##### PARSE AND VALIDATE INPUTS #####
	# print(fname_template)
	# ensure start and end ts are datetime objects
	start_ts = utils.parse_to_datetime(start_ts)
	end_ts = utils.parse_to_datetime(end_ts)

	# parse variables arg
	variables = parse_variables(variables)
	# validate locations arg and determine which method to use for subsetting
	loc_method = validate_locations(locations)

	# Get all of the filenames from download_base_path - in chronological order
	download_files = sorted(glob(os.path.join(nwm_data_dir, f'{fname_template}.f[0-9][0-9][0-9].conus.nc')))
	# print(download_files)

	##### LOADING NETCDF FILES #####
	# load the NWM data with multithreading
	# dataset_list will be in the same order as file list submitted, download_files
	print("Loading data...")
	dataset_dict = utils.multithreaded_loading(xr.open_dataset, download_files, num_threads)
	print("Loading complete.")

	print(f"Slicing data to get timestamps from {start_ts.strftime('%m-%d-%Y %H:%M:%S')} to {end_ts.strftime('%m-%d-%Y %H:%M:%S')}")
	timeslices = []

	vars_to_extract = [var for var in variables.values()]
	print(f"Extracing the following variables: {vars_to_extract}")

	# boolean flag to help us print a lcoation subsetting message ONLY once in the following loop
	loc_message = True

	##### TIME SLICING, VARIABLE SELECTION, & LOCATION SUBSETTING #####
	# iterate through every netcdf file (i.e. forecast timeslice)
	for i, ds in enumerate(dataset_dict.values()):
		# print(i)
		
		### Time Slicing ###
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

		### Variable Selection ###
		# writing the CRS first so that rioxarray can use it later for clipping
		# The CRS is initally stored as a data varaiable and consequently will get dropped when we drop the other variables unless we write it to the dataset's coordinates
		ds = ds.rio.write_crs(ds['crs'].attrs['esri_pe_string'])
		# list of variable names we want to drop
		drop_vars = list(set(ds.data_vars) - set(vars_to_extract))
		# drop variables that we don't want
		ds = ds.drop_vars(drop_vars)

		### Location Subsetting ###
		# subsetting by location really before concatenating all of the timeslices really speeds up the concatenation
		# One might argue that logically it makes more since to do the location subsetting at the end, but making the dataset lighter before concatenation I'd say is more consequential
		match loc_method:
			case 'bbox':
				# print a log message only during the first loop iteration... no need to repeat it for every timeslice
				if loc_message:
					print(f"Subsetting forcings grid by bounding box: {locations['bbox']}")
					loc_message = False
				# if a bounding box is specified, then subset the data by the bounding box
				ds = subset_by_bbox(ds, locations['bbox'])
			case 'points':
				# if a list of points is specified, then subset the data by the points
				if loc_message:
					print(f"Extracting the following points from the forcings grid: {locations['points']}")
					loc_message = False
				ds = subset_by_points(ds, locations['points'])
			case None:
				# if no locations are specified, then do not subset the data
				pass
		
		# print(f"Time step addded to concat list: {ds.time.values[0]}")
		timeslices.append(ds)
	
	print("Concatenating timeslices...")
	# concatenate all of the timeslices into a single dataset
	ds = xr.concat(timeslices, dim='time')
	print("Concatenation complete.")

	### Adding UTC Time Zone ###
	# saving attributes for later restoration
	time_attrs = ds['time'].attrs
	reftime_attrs = ds['reference_time'].attrs
	# creating new coords that are the time and reference_time coordinates to UTC
	ds = ds.assign_coords(timeutc=ds['time'].to_index().tz_localize("UTC"), reference_timeutc=ds['reference_time'].to_index().tz_localize("UTC"))
	# this gets rid of the new temporary coordinates that were created just created and replaces the original time and reference_time values with the new UTC values
	ds = ds.set_index(time='timeutc', reference_time='reference_timeutc')
	# now we just restore the original attributes
	ds.coords['time'].attrs = time_attrs
	ds.coords['reference_time'].attrs = reftime_attrs

	### Rename Dataset Variables if Applicable ###
	if any([k != v for k,v in variables.items()]):
		# if the user provided a dictionary of custom variable names to use, then ;et's rename the vars in the dataset
		ds = ds.rename({v:k for k,v in variables.items()})

	# reassigning so we can return one variable no matter the location method
	nwm_forcings_data = ds

	# point data extraction and nested dictionary creation
	if loc_method == 'points':
		# create a transformer to convert coords in the forcings CRS to WGS 1984
		nwm_to_wgs = pyproj.Transformer.from_crs(ds['crs'].attrs['crs_wkt'], 'EPSG:4326', always_xy=True)
		# make coordinate pairs of the x,y coordinates in the dataset (order is seeminly preserved this way)
		nearest_xy = [(x, y) for x, y in zip(ds['x'].values, ds['y'].values)]

		# get a dictionary of units
		var_units = {var : ds.data_vars[var].units for var in ds.data_vars}
		nwm_forcings_data = {}
		# iterate through the requested coords and the nearest xy grid coords
		for requested_coords, xy_coords in zip(locations['points'], nearest_xy):
			# print(requested_coords, xy_coords)
			nwm_forcings_data[requested_coords] = {}
			# print(f"x = {xy_coords[0]}, y = {xy_coords[1]}")
			# select the data for the given point location
			point_ds = ds.sel(x=xy_coords[0], y=xy_coords[1])
			for var in point_ds.data_vars:
				# convert each variable into a series
				nwm_forcings_data[requested_coords][var] = point_ds[var].to_pandas()

		# add unit info to series' names
		utils.add_units(nwm_forcings_data, var_units)

		# now convert the dataset's nearest x,y coordinates into lat,lon. These are the lat,lon coords in the dataset that are nearest to the requested lat,lon points
		nearest_xy_to_wgs = [nwm_to_wgs.transform(x, y) for (x, y) in nearest_xy]
		# create a dictionary where the keys are the user-requested lat,lon coords, and values are the nearest la,lon values in the dataset
		req_to_act_coords = {requested_coord : (lat, lon) for requested_coord, (lon, lat) in zip(locations['points'], nearest_xy_to_wgs)}
		# this list comp will add the actual grid coords to the dictionary for each requested point 
		[nwm_forcings_data[requested].update(grid_coords=nearest) for requested, nearest in req_to_act_coords.items()]

	print('TASK COMPLETE: PROCESS NWM FORCINGS DATA')
	return nwm_forcings_data