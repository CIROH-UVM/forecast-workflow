import datetime as dt
import concurrent.futures
import pandas as pd
import sh
'''
A module containing helper and utility functions for the data acquisiton modules
'''

def add_units(series_data, var_units):
	'''
	In-place function to add units to final series data returned by data modules. Units are stored 
	as the name attribute for each series in series_data.

	Args:
	-- series_data (dict) [req]: a nested dictionary of Pandas Series (returned by data modules).
	-- var_units (dict) [req]: a dictionary where keys are user-defined variable names and values are units

	Returns:
	None; function directly modifies series_data.
	'''
	for series_dict in series_data.values():
		for user_var_name, series in series_dict.items():
			series_dict[user_var_name] = series.rename(f'{user_var_name} ({var_units[user_var_name]})')

def download_data(url, filepath, use_google_bucket=False):
	"""
	Downloads, via curl, a file from a url to a specified filepath. Prints download progress.

	Args:
	-- url (str) [required]: full url to the file that is to be downloaded
	-- filepath (str) [required]: filepath to the download directory. should include desired name for file to be downloaded (ex /data/files/somedata.csv)

	"""
	# need to add curl output to list the old-fashioned way (no list comp) in case curl download is interupted/fails
	progress = []
	# '-#' turns on progress bar, '-o' writes download to filepath
	# _iter = '_err' allows us to iterate over the stderr produced by the call (for some reason progress bar is directed to stderr, not stdout)
	# _err_bufsize = 0 returns stderr one character at a time rather than as one very long string
	if use_google_bucket:
		print(f"Downloading file from URL: {url} using gsutil")
		call = sh.gsutil(['-m', 'cp', url, filepath], _iter = True, _out_bufsize = 100)
	else:
		print(f"Downloading file from URL: {url} using curl")
		call = sh.curl(['-#','-o', filepath, url], _iter = 'err', _err_bufsize = 0)
	try:
		for line in call:
			progress.append(line)
	except Exception as e:
		print("Download failed:")
		# if downloads fail, stop the run
		raise e
	else:
		# join strings, split by carriage return
		progress = ''.join(progress).split('\r')
		# log the last progress bar update (should be 100% ir download was finished, otherwise will show the point at which download ceased)
		print(progress[-1])
		print(f"Download Complete: {filepath}\n")

def fahr_to_celsius(fahren):
	'''
	Utility function to convert temperature data from Farenheit to Celsius

	Args:
	-- fahren (int, float, array-like) [req]: interger/float or array-like object to convert to Celcius from Fahrenheit
	
	Returns:
	New object in Celcius
	'''
	celsius = (fahren - 32) * 5 / 9
	return celsius
		
def generate_date_strings(start_date, end_date):
	"""
	Creates a list of date strings in the format "YYYMMDD" for the range provided by start_date and end_date

	Args:
	-- start_date (obj) [req]: First date in the date range. Should be a string in the format "YYYYMMDD", a dt.date, or dt.datetime object
	-- end_date (obj) [req]: Last date in the date range. Should be a string in the format "YYYYMMDD", a dt.date, or dt.datetime object

	Returns:
	A list of date strings in the format "YYYMMDD" 
	"""
	date_strings = []

	# if not a datetime object, convert start_date to one
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	for i in range(((end_date-start_date).days)+1):
		delta = dt.timedelta(days=i)
		date_strings.append((start_date+delta).strftime('%Y%m%d'))
	
	return date_strings

def generate_hours_list(num_hours, source, forecast_type='medium', archive=False):
	"""
	Creates a list of forecast hours to be downloaded

	Args:
	-- num_hours (int) [req]: how many hours of forecast data you want.
	-- saource (str) [req]: string indicating the forecast source. Valid values are 'gfs', 'nwm', and 'archive'.
	-- forecast_type (str) [opt]: The type of NWM forecast to grab. Valid values are 'short', 'medium', or 'long'.

	Returns:
	A list of hours in the format needed for the specifed forecast source. Or None if no valid source value is passed.
	"""
	match source:
		case 'gfs':
			if archive:
				return [f"{hour:03}" for hour in range(0, num_hours, 3)]
			if num_hours <= 120:
				return [f"{hour:03}" for hour in range(0, num_hours)]
			else:
				return [f"{hour:03}" for hour in range(0, 120)] + [f"{hour:03}" for hour in range(120, num_hours, 3)]
		case 'nwm':
			if forecast_type=='short' or forecast_type=='medium':
				return [f"{hour:03}" for hour in range(1, num_hours)]
			if forecast_type=='long':
				return [f"{hour:03}" for hour in range(6, num_hours, 6)]
		case 'buckets':
			if forecast_type=='short':
				return [f"{hour:03}" for hour in range(1, num_hours)]
			if forecast_type=='medium':
				if archive:
					return [f"{hour:03}" for hour in range(3, num_hours, 3)]
				else: return [f"{hour:03}" for hour in range(1, num_hours)]
			if forecast_type=='long':
				return [f"{hour:03}" for hour in range(6, num_hours, 6)]

def get_hour_diff(start_date, end_date):
	"""
	Utility function to get the number of hours 
	"""
	time_difference = end_date - start_date
	hours = time_difference.total_seconds() / 3600
	return int(hours)

def multithreaded_download(download_list, n_threads):
	"""
	Implements download_data() with mulithreading to speed up downloads

	Args:
	-- download_list (list) [required]: a list containing tuples of (url_to_download, download_destination_path)
	-- num_threads (int) [opt]: number of threads to use.
	"""
	urls, paths, use_google_buckets = zip(*download_list)
	print(f'Using {n_threads} threads for downloads')
	with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as executor:
		executor.map(download_data, urls, paths, use_google_buckets)

def multithreaded_loading(load_func, file_list, n_threads):
	"""
	The idea here is to pass the function that reads datasets (xr.open_dataset, cfgrib) with multithreading, 
	and return a dictionary with the structure {fname : dataset}

	Args:
	-- load_func (function) [required]: the function to use to open the datasets
	-- file_list (list) [required]: list of file names to load
	-- n_threads (int) [req]: number of threads to use for mulithreading

	Returns:
	a dictionary where every key is a filename in file_list, and each corresponding value is the datastruct loaded for that file
	"""
	with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as executor:
		datasets = executor.map(load_func, file_list)
	# returns datasets in the order in which they were submited; order of file_list will be preserved
	# dataset_list = [d for d in datasets]
	dataset_list = []
	# try to capture errors from failed file loads, but continue on since a missing file or too won't hurt AEM3D much
	for d in datasets:
		try:
			dataset_list.append(d)
		except Exception as e:
			print(e)
	return dict(zip(file_list, dataset_list))

def parse_to_datetime(date):
	"""
	Utility function to convert date or string to datetime obj, in UTC time

	Args:
	-- date (str, date, or datetime) [req]: the date to be parsed

	Returns:
	a datetime object representing the passed date
	"""
	# if already a datetime obj, just return the same obj; be sure to convert to UTC
	if isinstance(date, dt.datetime):
		return date.replace(tzinfo=dt.timezone.utc)
	# if a date obj, convert to datetime
	elif isinstance(date, dt.date):
		# return as a datetime, time set to midnight by default
		return dt.datetime.combine(date, dt.datetime.min.time()).replace(tzinfo=dt.timezone.utc)
	elif isinstance(date, str):
		if date == "":
			return None
		try:
			date = dt.datetime.strptime(date, "%Y%m%d").replace(tzinfo=dt.timezone.utc)
			return date
		except ValueError as ve:
			raise ValueError(f'{ve}. Please enter date strings in the format "YYYYMMDD"')
		except Exception as e:
			raise e

def smash_to_dataframe(series_data):
	'''
	A quick helper function that smashes the nested series dictionary returned by a get_data() into a single df.
	If timeseries do not align for each series in series_data, NaN's will be introduced during cocatenation to fill missing values.

	Args:
	-- series_data (dict) [req]: a nested dictionary of Pandas Series to be smashed into a dataframe
	
	Returns:
	A single timeseries dataframe for all of the variables and locations in series_data. Indexed by time (location name is a column).
	'''
	df_dict = {}
	for locname, series_dict in series_data.items():
		loc_df = pd.concat([pd.DataFrame(series) for series in series_dict.values()], axis=1)
		loc_df.insert(0, 'location', locname)
		df_dict[locname] = loc_df
	return pd.concat(df_dict.values())

def strip_non_numeric_chars(raw_series):
	'''
	Processes a Pandas Series to strip any non-numeric character, except a single decimal point and/or negative sign, from a string.
	
	Args:
	-- raw_series (Pandas Series) [req]: the series to parse

	Returns:
	A copy of the passed series with all non-numeric characters removed
	'''
	corrected_series = raw_series.copy()
	for i, value in enumerate(raw_series.astype(str)):
		corrected_value = ''
		# Check if the string is numeric after accounting for one decimal point
		if not value.replace('.', '', 1).replace('-','',1).isdigit() and value != 'nan':
			# Remove non-numeric characters and rejoin digits
			corrected_value = ''.join(filter(str.isdigit, value))
			if not corrected_value:
				corrected_value = 'nan'
			# print(f'Value at index {i} corrected from {value} to {corrected_value}')
			corrected_series[i] = corrected_value
	return corrected_series.astype(float)
