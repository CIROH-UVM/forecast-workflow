import datetime as dt
import concurrent.futures
import pandas as pd
import sh
'''
A module containing helper and utility functions for the data acquisiton modules
'''

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

def smash_dataframes(series_data):
	'''
	A quick helper function that smashes the nested series dictionary returned by a get_data() into a single df

		Args:
		-- series_data (dict) [req]: a nested dictionary of Pandas Series to be smashed into a dataframe
		Returns:
		A single dataframe containing all of the data for each. Indexed by time, latitude and longitude are columns.
	'''
	df_dict = {locname : pd.concat([s for s in series_dict.values()], axis=1) for locname, series_dict in series_data.items()}
	df = pd.concat(df_dict.values())
	df.index = pd.MultiIndex.from_product([df.index.unique(), df_dict.keys()], names = ['time','location'])
	return df
