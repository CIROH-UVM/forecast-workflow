import datetime as dt
import concurrent.futures
import pandas as pd
import sh
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
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

def combine_timeseries(primary, secondary, interval):
    """
    Combine two time series into one, ensuring all original rows in the primary
    time series are preserved and adding rows from the secondary time series if 
    their datetime index does not exist within the primary time series.

    Args:
    -- primary (pd.Series): The primary time series with a datetime index.
    -- secondary (pd.Series): The secondary time series with a datetime index.
    -- interval (str): The rounding interval for the datetime index (e.g., 'T' for minutes).

    Returns:
    pd.Series: A combined time series with source timesteps, preserving the primary 
               time series data and adding data from the secondary time series where 
               necessary.
    """
    # Convert the primary and secondary series to dataframes
    primary_df = primary.to_frame()
    secondary_df = secondary.to_frame()

    # Add a column for the original timestamps for each dataframe
    primary_df['timestamps'] = primary_df.index
    secondary_df['timestamps'] = secondary_df.index

    # Round the index to the specified interval to ensure alignment
    primary_df.index = primary_df.index.round(interval)
    secondary_df.index = secondary_df.index.round(interval)

    # Remove duplicate indices, which may be introduced by timestep rounding
    primary_df = primary_df[~primary_df.index.duplicated(keep='first')]
    secondary_df = secondary_df[~secondary_df.index.duplicated(keep='first')]

    # Combine the two dataframes, using the primary dataframe as the base
    combined_df = primary_df.combine_first(secondary_df)

    # Restore the original timestamps as the index
    combined_df.index = combined_df['timestamps'].rename('time')

    # Drop the temporary 'timestamps' column
    combined_df = combined_df.drop(columns=['timestamps'])

    # Convert the final dataframe back to a series
    return combined_df.squeeze()

def complete_dt_index(tsindx, start_t=None, end_t=None, gap_dur=None, **kwargs):
	'''
    Create a complete `pandas.DatetimeIndex` with no temporal gaps based on an uncertain or incomplete time series.
	Reports any gaps found in the passed index.

    This function generates a new DatetimeIndex that fills in temporal gaps from an existing DatetimeIndex. If the user
    does not provide a specific start time, end time, or gap duration, the function infers these from the input index. 
    The result is a continuous time series from `start_t` to `end_t` with a constant interval, `gap_dur`.

	Args:
	-- tsindx (pd.DatetimeIndex) [req]: The original datetime index that may contain temporal gaps. The function will infer the appropriate interval between timestamps if `gap_dur` is not provided.
    -- start_t (pd.timestamp, dt.datetime, None) [opt]: Expected start time of the timeseries. If None, uses the first timestamp in `tsindx`
    -- end_t (pd.timestamp, dt.datetime, None) [opt]: Expected end time of timeseries. If None, uses the last timestamp in `tsindx`. 
    -- gap_dur (pd.Timedelta or None) [opt]: The longest acceptable time interval between timestamps in the series. If None, eses the most common time interval found in tsindx
    -- **kwargs: Additional keyword arguments passed to `pd.date_range()`

	Returns:
    -- complete_index (pd.DatetimeIndex): 
            A complete datetime index from `start_t` to `end_t` with a constant interval of `gap_dur`. The new index 
            will have no temporal gaps, and the interval between consecutive timestamps will be uniform.
			Returns None if no gaps are found in the index.
	'''
	if gap_dur is None:
		# timedelta representing the smallest (temporally shortest) mode interval between timestamps in the index. 
		gap_dur = tsindx.to_series().diff().dropna().mode().min()
		print(f"Interval inferred from index: {gap_dur}")

	# If no expected start or end timestamps are passed, infer from index
	if start_t is None:
		start_t = tsindx[0]
	if end_t is None:
		end_t = tsindx[-1]
	
	# flag determining if any gaps at all were found
	gaps_found = False
	for i, ts in enumerate(tsindx):
		# For beginning of series index, check for a gap between the expected start ts and actual first ts
		if i == 0:
			b = start_t
		# for all other cases, check for gap with previous timestamp
		else:
			b = tsindx[i-1]
		delta = ts - b
		if delta > gap_dur:
			gaps_found = True
			print(f"Temporal gap in DatetimeIndex found: {delta}")
			print(f"\t{b} jumps to {ts}")

	# Finally, check if the last ts is the expected end ts or if there's an end_gap
	delta = end_t - ts
	if delta > gap_dur:
		gaps_found = True
		print(f"Temporal gap in DatetimeIndex found: {delta}")
		print(f"\t{ts} jumps to {end_t}")

	if gaps_found:
		# convert the timedelta into a DateOffset (which can be passed to 'freq' param in pd.date_range())
		toffset = pd.DateOffset(seconds=gap_dur.total_seconds())
		complete_index = pd.date_range(start_t, end_t, freq=toffset, **kwargs)
		return complete_index
	else: 
		print("No gaps were found in the index.")
		return None

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

def dfprint(df):
	'''
	Helper function to pretty-print an entire pandas dataframe in a Jupyter notebook
	Taken from: https://stackoverflow.com/questions/19124601/pretty-print-an-entire-pandas-series-dataframe

	Args:
	-- df (pd.DataFrame) [req]: the dataframe to print
	'''
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
		print(df)

		
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

def get_dt_index_gaps(tsindx, start_t=None, end_t=None, max_acceptable_td=None):
	'''
	Analyzes a pandas DatetimeIndex for temporal gaps. 'Gaps' are definied as a time delta between any two adjacent indices that is greater than max_acceptable_td.

	Args:
	-- tsindx (pd.DatetimeIndex) [req]: The original datetime index that may contain temporal gaps.
    -- start_t (pd.timestamp, dt.datetime, None) [opt]: Expected start time of the timeseries. If None, uses the first timestamp in `tsindx`
    -- end_t (pd.timestamp, dt.datetime, None) [opt]: Expected end time of timeseries. If None, uses the last timestamp in `tsindx`. 
    -- max_acceptable_td (pd.Timedelta or None) [opt]: The longest acceptable time interval between timestamps in the series. If None, uses the most common time interval (mode) found in tsindx

	Returns:
    gap_list (list): a list of tuples where each tuple is the pair of indices that form the gap
	'''
	# create empty gap list - will have list of tuples that mark the timestamps identifying the gaps
	gap_list = []

	if max_acceptable_td is None:
		# timedelta representing the smallest (temporally shortest) mode interval between timestamps in the index. 
		max_acceptable_td = tsindx.to_series().diff().dropna().mode().min()
		print(f"Max allowable time difference between two indices inferred from index: {max_acceptable_td}")

	# If no expected start or end timestamps are passed, infer from index
	if start_t is None:
		start_t = tsindx[0]
	if end_t is None:
		end_t = tsindx[-1]
	
	# flag determining if any gaps at all were found
	for i, ts in enumerate(tsindx):
		# For beginning of series index, check for a gap between the expected start ts and actual first ts
		if i == 0:
			b = start_t
		# for all other cases, check for gap with previous timestamp
		else:
			b = tsindx[i-1]
		delta = ts - b
		if delta > max_acceptable_td:
			print(f"Temporal gap in DatetimeIndex found: {delta}")
			print(f"\t{b} jumps to {ts}")
			gap_list.append((b, ts))

	# Finally, check if the last ts is the expected end ts or if there's an end_gap
	delta = end_t - ts
	if delta > max_acceptable_td:
		print(f"Temporal gap in DatetimeIndex found: {delta}")
		print(f"\t{ts} jumps to {end_t}")
		gap_list.append((ts, end_t))

	return gap_list

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
	Utility function to convert a date and time, represented by a variety of object types, to a datetime obj, in UTC time

	Args:
	-- date (str, date, dt.datetime, pd.Timestamp) [req]: the date to be parsed

	Returns:
	a datetime object representing the passed date
	"""
	# if already a datetime obj, just return the same obj; be sure to convert to UTC
	# pd.Timestamps are also considered an isntance of dt.datetimes, so will be parsed here
	if isinstance(date, dt.datetime):
		return date.replace(tzinfo=dt.timezone.utc)
	# if a date obj, convert to datetime
	elif isinstance(date, dt.date):
		# return as a datetime, time set to midnight by default
		return dt.datetime.combine(date, dt.datetime.min.time()).replace(tzinfo=dt.timezone.utc)
	elif isinstance(date, str):
		if date == "":
			return None
		# define acceptable string formats
		formats = ["%Y%m%d", "%Y%m%d%H"]
		for date_format in formats:
			try:
				date = dt.datetime.strptime(date, date_format).replace(tzinfo=dt.timezone.utc)
				return date
			except ValueError:
				continue
		raise ValueError(f'Invalid date string: {date}. Please enter date strings in the format "YYYYMMDD" or "YYYYMMDDHH".')

def plot_ts(series_list, scale='auto', **kwargs):
	"""
	Makes a single plot with one or many Pandas.Series. Assumes series are timeseries with Pandas.DatetimeIndex

	Args:
	-- series_list (list) [req]: a list of Pandas.Series to plot
	-- scale (str) [opt]: determines the x-axis tick mark scale. Options are 'hours', 'days', or 'auto', which will try to determine the best scale for the x-axis based on the series' indices'.
	-- kwargs [opt]: various kwargs to pass to matpllotlib plotting functions. 
			- interval (int): interval for x-axis tick marks, whether the scale be days, hours, etc
			- labels (list): list of labels to use for each series, respectively, in plt.plot()
			- colors (list): list of colros to use for each series, respectively, in plt.plot()
	
	Returns:
	the figure object created. 
	""" 
	fig, ax = plt.subplots(figsize=(10, 6))
	
	# Extract label from kwargs if it exists
	interval = kwargs.pop('interval', None)
	labels = kwargs.pop('labels', None)
	colors = kwargs.pop('colors', None)

	for i, series in enumerate(series_list):
		kwargs_single = dict()
		if labels: kwargs_single['label'] = labels[i]
		if colors: kwargs_single['color'] = colors[i]
		ax.plot(series.index, series.values, **kwargs_single)
	
	if labels:
		ax.legend()
	
	ax.set_title(f'{series.name} from {series.index[0].strftime("%m-%d-%Y")} to {series.index[-1].strftime("%m-%d-%Y")}')
	ax.set_ylabel(series.name)
	ax.set_xlabel('Date')

	# Format x-axis: show only month and day, and increase tick frequency
	if scale == 'days':
		ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))  # Month-Day format
		ax.xaxis.set_major_locator(mdates.DayLocator(interval=interval))
	elif scale == 'hours':
		ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))  # Month-Day format
		ax.xaxis.set_major_locator(mdates.HourLocator(interval=interval))
	elif scale == 'auto':
		# Auto format axis by default
		autolocator = mdates.AutoDateLocator()
		autoformatter = mdates.AutoDateFormatter(autolocator)
		ax.xaxis.set_major_formatter(autoformatter)
		ax.xaxis.set_major_locator(autolocator)

	# Rotate the x-axis labels
	ax.tick_params(axis='x', rotation=45)

	ax.grid(True)
	fig.tight_layout()

	# Remove plt.show() and add plt.close() to prevent immediate display
	plt.close(fig)
	return fig

def plot_nested_dict(data):
	"""
	Makes a figure with n subplots, where n is the number of vaiables in each location dict. Must have the same variables for each location
	
	Parameters:
		data (dict): A nested dictionary where the outer keys are categories 
					 and the inner keys are variables with their corresponding values.
	"""
	# Get the list of variables from the first category
	locations = list(data.keys())
	variables = list(data[locations[0]].keys())

	n = len(variables)  # Number of variables
	fig, axes = plt.subplots(n, 1, figsize=(10, n*4))  # Create subplots

	# define a color map
	cmap = plt.get_cmap("tab10") # 'Set1', 'Set2'
	colors = [cmap(i) for i in range(len(locations))]
	for i, var in enumerate(variables):
		for j, loc in enumerate(locations):
			series = data[loc][var]
			axes[i].plot(series.index, series.values, label=loc, color=colors[j])
			axes[i].set_title(var)
			axes[i].set_xlabel('Datetime')
			axes[i].set_ylabel(series.name)
			axes[i].legend()
			axes[i].grid()

	plt.tight_layout()  # Adjust layout to prevent overlap
	return fig

def report_gaps(series):
	"""
	This function identifies gaps (NaN values) within a pandas Series with a DatetimeIndex.
	Prints the start, end, and duration of each gap.
	
	Args:
	-- series (pd.Series) [req]: A pandas Series with a DatetimeIndex.
	
	Returns:
	None
	"""
	# Ensure the index is a DatetimeIndex
	if not isinstance(series.index, pd.DatetimeIndex):
		raise ValueError("The input series must have a DatetimeIndex.")
	
	# Find indices where values are NaN
	is_nan = series.isna()

	# Identify where the gaps (blocks of consecutive NaNs) are
	gap_start = None

	# were any gaps found at all? by default no
	gap_found = False
	for i in range(len(is_nan)):
		if is_nan.iloc[i] and gap_start is None:
			# Start of a new gap
			gap_found = True
			gap_start = series.index[i]
		elif not is_nan.iloc[i] and gap_start is not None:
			# End of the current gap
			gap_end = series.index[i - 1]
			duration = gap_end - gap_start
			# if the timedelta between the start and end of the gap is zero, then it reprewsents a single missing data point
			if duration.to_pytimedelta().total_seconds() == 0:
				print(f"Missing data at {gap_start}: {series.loc[gap_start]}")
			else:
				print(f"Data Gap Found: {duration}")
				print(f"\tFirst missing timestamp: {gap_start}")
				print(f"\tLast missing timestamp:  {gap_end}")
			gap_start = None
	# If the series ends with a gap
	if gap_start is not None:
		gap_end = series.index[-1]
		duration = gap_end - gap_start
		print(f"Data Gap Found: {duration}")
		print(f"\tFirst missing timestamp: {gap_start}")
		print(f"\tLast missing timestamp:  {gap_end}")
	if not gap_found:
		print("No gaps found in the timeseries")

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
