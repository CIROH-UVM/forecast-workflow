import data.utils as utils
import os
import datetime as dt
import tempfile as tf
import xarray as xr
import pandas as pd

# IMPORTANT: to get additional variables from CFS, you must add said variable's info to this dict
# Keys: CFS variable name as displayed on NetCDF Subset Service for given variable file
# fname: CFS variable file prefix as found on NCEI THREDDS
# abrv: abbreviation for variable name found in files downloaded from fileServer service
CFS_VAR_NAMES = {
	'Temperature_height_above_ground':
		{'fname':'tmp2m',
   		 'abrv':'t2m'},
	'Total_cloud_cover_entire_atmosphere_single_layer':
		{'fname':'tcdcclm',
   		 'abrv':'tcc'},
	'u-component_of_wind_height_above_ground':
		{'fname':'uwnd10m',
   		 'abrv':'u10'},
	'v-component_of_wind_height_above_ground':
		{'fname':'vwnd10m',
   		 'abrv':'v10'},
	'Specific_humidity_height_above_ground':
		{'fname':'q2m',
		 'abrv':'sh2'},
	'Precipitation_rate_surface':
		{'fname':'prate',
   		 'abrv':'prate'},
	'Snow_Phase_Change_Heat_Flux_surface':
		{'fname':'snohf',
   		 'abrv':'snohf'},
	'Downward_Short-Wave_Radiation_Flux_surface':
		{'fname':'dswsfc',
   		 'abrv':'dswrf'},
	'Pressure_surface':
		{'fname':'pressfc',
   		 'abrv':'sp'}
}

def construct_download_paths(reference_date, variables, data_dir, ncss=False):
	'''
	Constructs CFS download URLs and file paths based on reference date and download method.
	
	Args:
	-- reference_date (datetime) [req]: The forecast reference time.
	-- variables (dict) [req]: A dictionary mapping user-defined variable names to dataset-specific variable names.
	-- data_dir (str) [req]: Directory to store downloaded data.
	-- ncss (bool) [opt]: Whether to use the NetCDF Subset Service (NCSS). Defaults to False.
	
	Returns:
	-- tuple (list, list): file_urls, file_paths
	'''
	# define method of acquiring data from THREDDS server; by default, use fileServer method, which simply downloads the whole forecast file (no server-side subsetting)
	method = 'ncss' if ncss else 'fileServer'
	# base THREDDS url for the full 9-month operational CFSv2 forecast product
	base_url = f'https://www.ncei.noaa.gov/thredds/{method}/model-cfs_v2_for_ts'
	# define date path for forecast url/download directory structure
	date_path = f'{reference_date.year}/{reference_date.strftime("%Y%m")}/{reference_date.strftime("%Y%m%d")}/{reference_date.strftime("%Y%m%d%H")}'
	# create date-specific directory tree for storing files
	date_dir = os.path.join(data_dir, 'cfs', date_path)
	# make the date directory for CFS data storeage if they do not exist already
	os.makedirs(date_dir, exist_ok=True)
	print(f'\tData will be stored at: {date_dir}')
	
	# create list of CFS filenames to download/store
	grib_files = []
	for user_name, long_name in variables.items():
		var_prefix = CFS_VAR_NAMES[long_name]['fname']
		# adjust prefix for wind (u & v components stored in  the same file)
		if 'wnd10m' in var_prefix: var_prefix = 'wnd10m'
		file_name = f'{var_prefix}.01.{reference_date.strftime("%Y%m%d%H")}.daily.grb2'
		grib_files.append(file_name)
	# there will be duplicate grib files for wind since u and v components are separate vars in the same file
	# so remove the duplicate:
	grib_files = list(set(grib_files))
	
	file_urls = [f"{base_url}/{date_path}/{f}" for f in grib_files]
	
	# add an ncss tag to beginning of file name so users know that said file has been pre-sliced
	if ncss: grib_files = ["ncss_"+f for f in grib_files]
	file_paths = [os.path.join(date_dir, f) for f in grib_files]
	
	return file_urls, file_paths

def download_full_cfs(reference_date, variables, data_dir, num_threads=int(os.cpu_count()/2)):
	'''
	Downloads full CFS forecast files using the fileServer method.
	
	Args:
	-- reference_date (datetime) [req]: The forecast reference time.
	-- variables (dict) [req]: A dictionary of requested variables.
	-- data_dir (str) [req]: Directory to store downloaded files.
	-- num_threads (int) [opt]: Number of threads for downloading. Defaults to half CPU cores.
	
	Returns:
	-- list: Paths to downloaded files.
	'''
	print('TASK INITIATED: Download CFS 9-month operational forecast from NCEI THREDDS fileServer')
	print(f'\tFORECAST REFERENCE TIME: {reference_date}')
	print(f'\tFORECAST CYCLE: {reference_date.hour:02d}')
	# construct full grib file URLS and full paths to store data at
	file_urls, file_paths = construct_download_paths(reference_date, variables, data_dir, ncss=False)
	# execute downloads
	execute_downloads(file_urls, file_paths, num_threads)
	print('TASK COMPLETE: CFS DOWNLOAD')
	return file_paths

def download_ncss(reference_date, variables, locations, start_date, end_date, data_dir, num_threads=int(os.cpu_count()/2)):
	'''
	Downloads subsetted CFS forecast data using the NCSS method.
	
	Args:
	-- reference_date (datetime) [req]: The forecast reference time.
	-- variables (dict) [req]: A dictionary of requested variables.
	-- locations (dict) [req]: Mapping of location names to lat/lon tuples.
	-- start_date (datetime) [req]: Start date for data retrieval.
	-- end_date (datetime) [req]: End date for data retrieval.
	-- data_dir (str) [req]: Directory to store downloaded files.
	-- num_threads (int) [opt]: Number of threads for downloading. Defaults to half CPU cores.
	
	Returns:
	-- list: Paths to downloaded files.
	'''
	print('TASK INITIATED: Download CFS 9-month operational forecast from NCEI THREDDS NetCDF Subset Service')
	print(f'\tFORECAST REFERENCE TIME: {reference_date}')
	print(f'\tFORECAST CYCLE: {reference_date.hour:02d}')
	print(f'\tSTART DATE: {start_date}')
	print(f'\tEND DATE: {end_date}')

	# remapping longitude using map_function
	remapped_locs = remap_to_cfs_coords(locations)
	# now get the bounding box for the netcdf subset slice
	min_lat, max_lat, min_lon, max_lon = utils.get_bounding_box(remapped_locs)
	# construct base URLs for NCSS downloads, and full file paths
	file_urls, file_paths = construct_download_paths(reference_date, variables, data_dir, ncss=True)
	
	file_urls = [
		f"{url}?var=all&north={max_lat}&west={min_lon}&east={max_lon}&south={min_lat}" 
		f"&horizStride=1&time_start={start_date.strftime('%Y-%m-%dT%H')}%3A00%3A00Z&time_end={end_date.strftime('%Y-%m-%dT%H')}%3A00%3A00Z"
		"&timeStride=1&vertStride=1"
		for url in file_urls
	]
	execute_downloads(file_urls, file_paths, num_threads)
	print('TASK COMPLETE: CFS DOWNLOAD')
	return file_paths

def execute_downloads(file_urls, file_paths, num_threads):
	'''
	Executes the download process using multi-threading.
	
	Args:
	-- file_urls (list) [req]: List of URLs to download.
	-- file_paths (list) [req]: List of local file paths to save downloaded files.
	-- num_threads (int) [req]: Number of threads to use for parallel downloading.
	'''
	# making a download list for multithreaded downloading
	download_list = []
	for i, (url, fpath) in enumerate(zip(file_urls, file_paths)):
		if not os.path.exists(fpath):
			download_list.append((url, fpath, False))
		else: print(f'Skipping download; {os.path.basename(fpath)} found at: {fpath}')
	# print(f"Download list:")
	# print(download_list)
	if download_list:
		utils.multithreaded_download(download_list, num_threads)

def process_full_cfs(ds, variables, locations):
	'''
	A function to process a full CFS forecast dataset downloaded from the THREDDS fileServer and extract time-series data for specified locations and variables.  
	Filters the dataset to include only requested variables, maps user-provided locations to the nearest available coordinates in the dataset,  
	and returns a nested dictionary of Pandas Series objects for each variable at each location.

		Args:
		-- ds (xarray.Dataset) [req]: The dataset containing complete CFS forecast data.
		-- variables (dict) [req]: A dictionary mapping user-defined variable names to dataset-specific variable names (keys of CFS_VAR_NAMES).
		-- locations (dict) [req]: A dictionary of locations to extract data for. Keys are user-defined location names, values are (lat, lon) tuples.

		Returns:
		A nested dictionary where:
		- 1st-level keys are user-defined location names.
		- 2nd-level keys are user-defined variable names.
		- Values are Pandas Series objects containing time-series data for the corresponding variable at the location.
		
		Notes:
		- Maps user-provided location coordinates to the nearest available dataset coordinates.
		- Drops unnecessary coordinates for efficiency.
		- Converts extracted data into a Pandas DataFrame and localizes time indices to UTC.
		- Groups data by unique latitude-longitude pairs for optimized extraction.
		- Attaches unit metadata to the extracted data.
	'''
	print("TASK INITIATED: Process CFS data")
	# making a dictionary of units for each variable... to use later
	var_units = {user_var:ds[CFS_VAR_NAMES[long_name]['abrv']].units for user_var, long_name in variables.items()}

	# get lists of remapped lat/lons
	remapped_locs = remap_to_cfs_coords(locations)
	lats, lons = zip(*remapped_locs.values())

	# filter by latitude and longitude values
	ds = ds.sel(latitude=list(lats), longitude=list(lons), method='nearest')

	# get the dataset-approximated lat and lon values, report lat/lon approximation
	approx_lats, approx_lons = ds['latitude'].values, ds['longitude'].values
	approx_locs = {loc_name:(ap_lat,ap_lon) for (loc_name, ap_lat, ap_lon) in zip(locations.keys(), approx_lats, approx_lons)}
	print('NOTE: location coordinates adjusted to nearest values available in dataset:')
	for (loc_name, user_coords), approx_coords in zip(locations.items(), approx_locs.values()):
		print(f'\tFor {loc_name}: {user_coords} -> {approx_coords}')

	# create new coordinate, based off of the reference time and steps (timedeltas, already an index)
	# this new coordinate is the unlocalized datetime index
	ds = ds.assign_coords(time=('step',(pd.Timestamp(ds.time.values) + pd.to_timedelta(list(ds.step.values)))))
	# drop all coordinates except for the new time coordinate we just created
	coords_to_drop = [coord for coord in ds.coords._names if coord not in ['latitude', 'longitude', 'time', 'step']]
	ds = ds.drop(coords_to_drop)
	# set the time index as the 'step' dimension on the data
	ds = ds.swap_dims({'step':'time'})
	# squeeze the data variabes from 3d to 2d ((n, 1, 1) -> (n,) where n is the number of timesteps), for faster loading into pandas
	ds_squeezed = ds.map(lambda var: var.squeeze())
	# drop duplactes lat/lons after selecting, which occurs when two different locations have the same lat/lon values
	ds_squeezed = ds_squeezed.drop_duplicates(dim='longitude').drop_duplicates(dim='latitude')
	# convert to dataframe - will take 2-3 minutes
	df = ds_squeezed.to_dataframe()
	# reset the index and make 'time' the only index (latitutde and longitude are now columns)
	df = df.reset_index()
	df = df.set_index('time')
	# localize datetime index to UTC time
	df.index = df.index.tz_localize(dt.timezone.utc)

	# we want the set of coordinate pairs to pull out of the dataframe, for efficiency (don't need to operate on the same data multiple times)
	# we want a set because there might be the identical coord pairs for multiple locations
	unique_coords = set(approx_locs.values())
	# make a dictionary of unique groups of coordinates pair dataframes
	loc_groups = {group_name: group for group_name, group in df.groupby(['latitude', 'longitude']) if group_name in unique_coords}

	# now build the CFS nested dictionary
	cfs_data = {}
	for loc_name, coords in approx_locs.items():
		cfs_data[loc_name] = {}
		for user_name, long_name in variables.items():
			cfs_data[loc_name][user_name] = loc_groups[coords][CFS_VAR_NAMES[long_name]['abrv']]
	utils.add_units(cfs_data, var_units)
	return cfs_data

def process_ncss(ds, variables, locations):
	'''
	A function to process NetCDF Subset Service (NCSS) CFS dataset and extract time-series data for specified locations and variables.  
	Filters the dataset to include only requested variables, maps user-provided locations to the nearest available coordinates in the dataset,  
	and returns a nested dictionary of Pandas Series objects for each variable at each location.

		Args:
		-- ds (xarray.Dataset) [req]: The dataset containing CFS forecast data.
		-- variables (dict) [req]: A dictionary mapping user-defined variable names to dataset-specific variable names (keys of CFS_VAR_NAMES).
		-- locations (dict) [req]: A dictionary of locations to extract data for. Keys are user-defined location names, values are (lat, lon) tuples.

		Returns:
		A nested dictionary where:
		- 1st-level keys are user-defined location names.
		- 2nd-level keys are user-defined variable names.
		- Values are Pandas Series objects containing time-series data for the corresponding variable at the location.
		
		Notes:
		- Drops unrequested variables from the dataset.
		- Adjusts user-provided location coordinates to the nearest available values in the dataset.
		- Converts extracted data into Pandas Series, ensuring time indices are localized to UTC.
		- Attaches unit metadata to the extracted data.
	'''
	print("TASK INITIATED: Process CFS data")
	# create a list of vars to drop
	vars_to_drop = [v for v in list(ds.data_vars.keys()) if v not in list(variables.values())]
	# making a dictionary of units for each variable
	var_units = {user_var:ds[long_name].units for user_var, long_name in variables.items()}

	# get lists of remapped lat/lons
	remapped_locs = remap_to_cfs_coords(locations)
	lats, lons = zip(*remapped_locs.values())
	# filter by lat/lons, drop unrequested vars
	ds = ds.sel(lat=list(lats), lon=list(lons), method='nearest').drop_vars(vars_to_drop)

	# get the dataset-approximated lat and lon values, report lat/lon approximation
	approx_lats, approx_lons = ds['lat'].values, ds['lon'].values
	approx_locs = {loc_name:(ap_lat,ap_lon) for (loc_name, ap_lat, ap_lon) in zip(locations.keys(), approx_lats, approx_lons)}
	print('NOTE: location coordinates adjusted to nearest values available in dataset:')
	for (loc_name, user_coords), approx_coords in zip(locations.items(), approx_locs.values()):
		print(f'\tFor {loc_name}: {user_coords} -> {approx_coords}')

	# group dataset by latitude
	grouped_ds = ds.groupby('lat')
	# get grouped datasets
	location_datasets = [grouped_ds[ap_lat].sel(lon=ap_lon) for ap_lat, ap_lon in zip(approx_lats, approx_lons)]
	cfs_data = {}
	for _, (loc_ds, (loc_name, coords)) in enumerate(zip(location_datasets, remapped_locs.items())):
		# convert to dataframe, set time as index, drop uneeded coords, drop duplicates and rename columns
		loc_df = loc_ds.to_dataframe().reset_index().set_index('time').drop(['height_above_ground'], axis=1).drop_duplicates().rename(columns = {ds_name : user_name for user_name, ds_name in variables.items()})
		loc_series = {}
		for user_name, ds_name in variables.items():
			var_series = loc_df[user_name].dropna()
			var_series = var_series[~var_series.index.duplicated()]
			# localize datetime index to UTC time
			var_series.index = var_series.index.tz_localize(dt.timezone.utc)
			loc_series[user_name] = var_series
		cfs_data[loc_name] = loc_series

	utils.add_units(cfs_data, var_units)
	return cfs_data

def remap_to_cfs_coords(locations):
	'''
	Given a standard location dictionary, will  return a tuple containing two paired lists of -90 - 90 latitudes and 0-360 longitudes
	
	Returns:
	--lats: list of latitudes
	--lons: corresponding list of longitudes
	'''
	# function to convert longitude on the -180 to 180 scale to 0 to 360 (Higher longitude makes east side of boundary box, lower value is west boundary)
	map_function = lambda lon: 360 + lon if (lon < 0) else lon
	# remapping longitude using map_function
	remapped_locs = {key:(lat, map_function(lon)) for key, (lat, lon) in locations.items()}
	return remapped_locs

def get_data(start_date,
			 end_date,
			 locations,
			 variables,
			 reference_date=None,
			 ncss=False,
			 end_date_exclusive=True,
			 data_dir=tf.gettempdir(),
			 num_threads=int(os.cpu_count()/2)):
	'''
	A function to download and process CFS Operational 9-month forecast data, returning a nested dictionary of Pandas Series  
	for each requested variable at each location. Supports both full dataset downloads and NetCDF Subset Service (NCSS)  
	for optimized data retrieval.

	CFS operational 9-month forecast data is available from 04/01/2011 to present.  
	See available data here: https://www.ncei.noaa.gov/thredds/catalog/model-cfs_v2_for_ts/catalog.html  
	CFS products overview: https://www.ncei.noaa.gov/products/weather-climate-models/climate-forecast-system  

		Args:
		-- start_date (str, date, or datetime) [req]: The start date for which to retrieve data.
		-- end_date (str, date, or datetime) [req]: The end date for which to retrieve data.
		-- locations (dict) [req]: A dictionary mapping location names (or station IDs) to latitude-longitude tuples.
		-- variables (dict) [req]: A dictionary mapping user-defined variable names to dataset-specific variable names.
		-- reference_date (str, date, or datetime) [opt]: The forecast reference time, i.e., the date and time at which  
		the forecast was initialized. Defaults to start_date if None.
		-- ncss (bool) [opt]: Whether to use the NetCDF Subset Service (NCSS) to retrieve data.  
		Defaults to False. NCSS improves efficiency by pre-filtering data but may not be suitable for shared network drives.
		-- end_date_exclusive (bool) [opt]: Whether to exclude the end date from the time series. Defaults to True.
		-- data_dir (str) [opt]: The directory to store downloaded data. Defaults to the OS's temporary directory.
		-- num_threads (int) [opt]: The number of threads to use for downloading data. Defaults to half of the available CPU cores.

		Returns:
		A nested dictionary where:
		- 1st-level keys are user-defined location names.
		- 2nd-level keys are user-defined variable names.
		- Values are Pandas Series objects containing time-series data for the corresponding variable at the location.

		Notes:
		- Determines the appropriate data retrieval method (fileServer vs. NCSS) based on the `ncss` flag.
		- Parses forecast reference time and validates it against CFS forecast cycles (0, 6, 12, 18 UTC).
		- Remaps user-specified locations to the nearest available dataset grid points.
		- Supports multi-threaded downloading for faster data retrieval.
		- Filters time-series data based on the specified start and end dates.
		- Attaches unit metadata to the extracted data.
	'''
	# Validate and process forecast dates
	start_date, end_date, reference_date = utils.validate_forecast_cycle(start_date, end_date, reference_date)

	# Download data
	if ncss:
		file_paths = download_ncss(reference_date, variables, locations, start_date, end_date, data_dir, num_threads)
	else:
		file_paths = download_full_cfs(reference_date, variables, data_dir, num_threads)
	
	# Load dataset
	ds = xr.open_mfdataset(file_paths, engine=('netcdf4' if ncss else 'cfgrib'))
	
	# Process dataset
	cfs_data = process_ncss(ds, variables, locations) if ncss else process_full_cfs(ds, variables, locations)
	
	print("Slicing CFS timeseries:")
	print(f'\tSTART DATE: {start_date}')
	print(f'\tEND DATE: {end_date}')
	
	# Filter time series according to start and end date
	for loc_name, var_dict in cfs_data.items():
		for var_name, var_series in var_dict.items():
			var_dict[var_name] = var_series.loc[start_date:end_date]
			if end_date_exclusive:
				var_dict[var_name].drop(end_date, inplace=True)
	
	return cfs_data