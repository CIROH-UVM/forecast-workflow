from data.utils import add_units, multithreaded_download, parse_to_datetime
import os
import tempfile as tf
import xarray as xr


# IMPORTANT: to get additional variables from CFS, you must add said variable's info to this dict
# Keys: CFS variable name as displayed on NetCDF Subset Service for given variable file
# Values: CFS variable file prefix as found on NCEI THREDDS
cfs_variables_dict = {
	'Temperature_height_above_ground': 'tmp2m',
	'Total_cloud_cover_entire_atmosphere_single_layer': 'tcdcclm',
	'u-component_of_wind_height_above_ground': 'uwnd10m',
	'v-component_of_wind_height_above_ground': 'vwnd10m',
	'Specific_humidity_height_above_ground': 'q2m',
	'Precipitation_rate_surface': 'prate',
	'Snow_Phase_Change_Heat_Flux_surface': 'snohf',
	'Downward_Short-Wave_Radiation_Flux_surface': 'dswsfc',
	'Pressure_surface':'pressfc'
}

def get_data(start_date,
			 end_date,
			 locations,
			 variables,
			 reference_date = None,
			 data_dir = tf.gettempdir(),
			 num_threads = int(os.cpu_count()/2)):
	'''
	A function to download and process CFS Operational 9-month forecast data to return nested dictionary of pandas series for each variable, for each location.
	CFS operational 9-month forecast is available from 04/01/2011 to present. See available data here: https://www.ncei.noaa.gov/thredds/catalog/model-cfs_v2_for_ts/catalog.html
	CFS products overview: https://www.ncei.noaa.gov/products/weather-climate-models/climate-forecast-system

		Args:
		-- start_date (str, date, or datetime) [req]: the start date for which to grab data.
		-- end_date (str, date, or datetime) [req]: the end date for which to grab data.
		-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get data for.
		-- variables (dict) [req]: a dictionary of variables to download. Keys should be user-defined names, values should be dataset-specific var names (keys of cfs_variables_dict).
		-- reference_date (str, date, or datetime) [opt]: The forecast reference time, i.e. the date and time at which the forecast is launched. Set to start_date if None.
		-- data_dir (str) [opt]: directory to store downloaded data. Defaults to OS's default temp directory.
		
		Returns:
		CFS timeseries data for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
		are variables names and values are the respective data in a Pandas Series object.
	'''
	# base THREDDS url for the forecast product 
	ncss_base_url = 'https://www.ncei.noaa.gov/thredds/ncss/model-cfs_v2_for_ts'

	##### Parsing arguments #####
	# add variable file prefixes to variables dictionary
	variables = {user_name:(cfs_variables_dict[long_name], long_name) for user_name, long_name in variables.items()}

	print('TASK INITIATED: Download CFS 9-month operational forecast')
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)
	'''
	if no reference time is passed, reference time is inferrred from start_date
		if start date does not have a valid forecast cycle, an error is raised
	if a reference time is passed, then it is considered the forecast reference time and cycle and start date is used to slice forecast timeseries
	if an invalid reference time is passed, an error is raised
	'''
	# if no reference time is passed, set it to start date
	if reference_date is None:
		reference_date = start_date
		print(f'Forecast reference time and cycle inferred from start_date: {reference_date}')
	# if reference time was passed, parse it
	else: reference_date = parse_to_datetime(reference_date)
	# if the reference time does not have a valid forecast cycle, raise error
	if reference_date.hour not in [0, 6, 12, 18]:
		raise ValueError(f'Invalid forecast reference time: {reference_date}. Reference time hour indicates forecast cycle and must be 0, 6, 12, or 18')
	# if start date is before reference time, raise an error. With reference time passed, start date is used to slice forecast timeseries
	if start_date < reference_date:
		raise ValueError(f'Forecast start date: {start_date} comes before forecast reference time: {reference_date}')
	print(f'\tFORECAST REFERENCE TIME: {reference_date}')
	print(f'\tFORECAST CYCLE: {reference_date.hour:02d}')
	print(f'\tSTART DATE: {start_date}')
	print(f'\tEND DATE: {end_date}')

	# function to convert longitude on the -180 to 180 scale to 0 to 360 (Higher longitude makes east side of boundary box, lower value is west boundary)
	map_function = lambda lon: 360 + lon if (lon < 0) else lon
	# remapping longitude using map_function
	remapped_locs = {key:(lat, map_function(lon)) for key, (lat, lon) in locations.items()}
	lats, lons = zip(*remapped_locs.values())
	# getting mins and maxs for boundary box
	# if there's only one lat/lon value requested pad either side of lat/lon to get boundaries
	if len(set(lons)) == 1:
		# just grab the first value, they should all be the same
		max_lon, min_lon = lons[0]+0.25, lons[0]-0.25
		# padding each boundary with a quarter degree to ensure the requested lat/lons are in the boundary box (cannot assume boundary is inclusive)
	else: max_lon, min_lon = max(lons)+0.25, min(lons)-0.25
	if len(set(lats)) == 1:
		# just grab the first value, they should all be the same
		max_lat, min_lat = lats[0]+0.25, lats[0]-0.25
	else: max_lat, min_lat = max(lats)+0.25, min(lats)-0.25

	##### Downloading data #####
	# define date path for forecast url/download directory structure
	date_path = f'{reference_date.year}/{reference_date.strftime("%Y%m")}/{reference_date.strftime("%Y%m%d")}/{reference_date.strftime("%Y%m%d%H")}'
	# create date-specific directory tree for storing files
	date_dir = os.path.join(data_dir, 'cfs', date_path)
	# date-specific url for THREDDS server
	ncss_date_url = f'{ncss_base_url}/{date_path}'
	print(f'\tData will be stored at: {date_dir}')
	# make the date directory for CFS data storeage if they do not exist already
	if not os.path.exists(date_dir):
		os.makedirs(date_dir)

	# create list of CFS filenames to download/store
	grib_files = []
	for var_prefix, var_longname in variables.values():
		# adjust prefix for wind (u & v components stored in  the same file)
		if 'wnd10m' in var_prefix: var_prefix = 'wnd10m'
		grib_files.append(f'{var_prefix}.01.{reference_date.strftime("%Y%m%d%H")}.daily.grb2')
	# there will be duplicate grib files for wind since u and v components are separate vars in the dame file
	# so remove the duplicate:
	grib_files = list(set(grib_files))

	# construct full NCSS urls for downloading
	ncss_file_urls = [
	f'{ncss_date_url}/{grib_file}?var=all'
	f'&north={max_lat}&west={min_lon}&east={max_lon}&south={min_lat}'
	f'&horizStride=1'
	f'&time_start={start_date.strftime("%Y-%m-%dT%H")}%3A00%3A00Z'
	f'&time_end={end_date.strftime("%Y-%m-%dT%H")}%3A00%3A00Z'
	f'&timeStride=1&vertStride=1'
	for grib_file in grib_files]

	file_paths = [f'{date_dir}/{file}' for file in grib_files]

	# making a download list for mulitthreaded downloading
	download_list = []
	for i, (url, fpath) in enumerate(zip(ncss_file_urls, file_paths)):
		if not os.path.exists(fpath):
			download_list.append((url, fpath, False))
		else: print(f'Skipping download; {os.path.basename(fpath)} found at: {date_dir}')

	# download if the list is not empty
	if download_list:
		multithreaded_download(download_list, num_threads)
	print('TASK COMPLETE: CFS DOWNLOAD')

	##### Dataset Processing #####
	ds = xr.open_mfdataset(file_paths, engine='netcdf4')

	# create a list of vars to drop
	vars_to_drop = [v for v in list(ds.data_vars.keys()) if v not in list(zip(*variables.values()))[1]]
	# making a dictionary of units for each variable
	var_units = {user_var:ds[ds_var].units for user_var, (var_prefix, ds_var) in variables.items()}
	# filter by lat/lons, drop unrequested vars
	ds = ds.sel(lat=list(lats), lon=list(lons), method='nearest').drop_vars(vars_to_drop)

	# get the dataset-approximated lat and lon values, report lat/lon approximation
	approx_lats, approx_lons = ds['lat'].values, ds['lon'].values
	approx_locs = {loc_name:(ap_lat,ap_lon) for (loc_name, ap_lat, ap_lon) in zip(locations.keys(), approx_lats, approx_lons)}
	print('NOTE: location coordinates adjusted to nearest values available in dataset:')
	for (loc_name, user_coords), approx_coords in zip(locations.items(), approx_locs.values()):
		print(f'For {loc_name}: {user_coords} -> {approx_coords}')

	# group dataset by latitude
	grouped_ds = ds.groupby('lat')
	# get grouped datasets
	location_datasets = [grouped_ds[ap_lat].sel(lon=ap_lon) for ap_lat, ap_lon in zip(approx_lats, approx_lons)]
	cfs_data = {}
	for _, (loc_ds, (loc_name, coords)) in enumerate(zip(location_datasets, remapped_locs.items())):
		# convert to dataframe, set time as index, drop uneeded coords, drop duplicates and rename columns
		loc_df = loc_ds.to_dataframe().reset_index().set_index('time').drop(['height_above_ground'], axis=1).drop_duplicates().rename(columns = {ds_name : user_name for user_name, (var_prefix, ds_name) in variables.items()})
		loc_series = {}
		for user_name in variables.keys():
			var_series = loc_df[user_name].dropna()
			var_series = var_series[~var_series.index.duplicated()]
			# localize datetime index to UTC time
			var_series.index = var_series.index.tz_localize(dt.timezone.utc)
			loc_series[user_name] = var_series
		cfs_data[loc_name] = loc_series

	add_units(cfs_data, var_units)
	return cfs_data