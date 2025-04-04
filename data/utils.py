from bs4 import BeautifulSoup
import datetime as dt
import concurrent.futures
import pandas as pd
import requests
import sh
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
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

def calculate_rh(psfc, q2, t2):
	'''
	Calculate a relative humidty series. Input series must have identical indices.
	Source: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity

	Args:
	-- psfc (pd.Series) [req]: surface pressure series, in pascals (Pa).
	-- q2 (pd.Series) [req]: specifc humidity series at 2m (dimensionless)
	-- t2 (pd.Series) [req]: temperature at 2m series (k)
	'''
	# reference temperature (k), equal to 0 degrees celsius (freezing point of water)
	T0 = 273.15
	# Calculate rh using element-wise operations
	rh = psfc * q2 * 0.263 * np.exp((17.67 * (t2 - T0)) / (t2 - 29.65)) ** -1
	# cap rhum at 100%
	cap_rhum_at_100 = lambda x: 100 if x > 100 else x
	rh = rh.apply(cap_rhum_at_100).rename('RH2 (%)')
	return rh

def calculate_rh_LeBauer(psfc, q2, t2):
	'''
	Calculate a relative humidty series using equations authored by David LeBauer. Input series must have identical indices.
	Source: https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html#eq:1

	Args:
	-- psfc (pd.Series) [req]: surface pressure series, in pascals (Pa).
	-- q2 (pd.Series) [req]: specifc humidity series at 2m (dimensionless)
	-- t2 (pd.Series) [req]: temperature at 2m series (k)
	'''
	##### Original Code #####
	#   Relative Humidity : WRF, Calculated from Q2
	#
	#
	# rhum = dict()
	# for zone in climate[cs['RHUM']]['T2'].keys():
	#     # Using an equation ported from metuils.r that was authored by David LeBauer
	#     # es = 6.112 * np.exp((17.67 * t2) / (t2 + 243.5))
	#     # constants synched with https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html

	#     q2 = climate[cs['RHUM']]['Q2'][zone]
	#     psfc = climate[cs['RHUM']]['PSFC'][zone]
	#     t2 = climate[cs['RHUM']]['T2'][zone]
	#     es = 6.1078 * np.exp((17.2694 * t2) / (t2 + 237.30))    
	# Saturation vapor pressure over water
	#     e = q2 * psfc * 0.01 / (0.378 * q2 + 0.622)        # factor of 0.01 is converting Pascal to mbar
	#     rhum_temp = e/es
	#     rhum_temp[rhum_temp>1] = 1
	#     rhum_temp[rhum_temp<0] = 0
	#     rhum[zone] = rhum_temp

	es = 6.1078 * np.exp((17.2694 * (t2 - 273.16)) / (t2 - 35.86))
	# Saturation vapor pressure over water
	# factor of 0.01 is converting Pascal to mbar
	e = (q2 * psfc * 0.01) / (0.378 * q2 + 0.622)
	rhum = (e/es) * 100
	rhum[rhum>100] = 100
	rhum[rhum<0] = 0
	rhum.name = 'RH2 (%)'
	return rhum

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
    -- gap_dur (pd.Timedelta or None) [opt]: The longest acceptable time interval between timestamps in the series. If None, uses the most common (mode) time interval found in tsindx
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

def get_bounding_box(locations):
	'''
	Computes the geospatial bounding box for a given dictionary of location coordinates. IOW calculates min and max for lat/lon values given 
	the coordinates passed, plus adds a quarter degree of padding to ensure location is included in the box.
	
	Args:
	-- locations (dict) [req]: A dictionary mapping location names to latitude-longitude tuples.
	
	Returns:
	-- tuple (float, float, float, float): min_lat, max_lat, min_lon, max_lon
	'''
	lats, lons = zip(*locations.values())
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
	
	return min_lat, max_lat, min_lon, max_lon
		
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
	-- source (str) [req]: string indicating the forecast source. Valid values are 'gfs', 'nwm', and 'archive'.
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
	-- scale (str) [opt]: determines the x-axis tick mark scale. Options are:
	 	- 'years'
		- 'months'
		- 'weeks'
		- 'days'
		- 'hours'
		- 'auto' : will try to determine the best scale for the x-axis based on the series' indices'
	-- kwargs [opt]: various kwargs to pass to matpllotlib plotting functions. 
		- interval (int): interval for x-axis tick marks, whether the scale be days, hours, etc
		- labels (list): list of labels to use for each series, respectively, in plt.plot()
		- colors (list): list of colors to use for each series, respectively, in plt.plot()
		- title (str): a custom title for the plot
	
	Returns:
	the figure object created. 
	""" 
	fig, ax = plt.subplots(figsize=(10, 6))
	
	# Extract label from kwargs if it exists
	interval = kwargs.pop('interval', 1)
	labels = kwargs.pop('labels', None)
	colors = kwargs.pop('colors', None)
	title = kwargs.pop('title', None)

	for i, series in enumerate(series_list):
		kwargs_single = dict()
		if labels: kwargs_single['label'] = labels[i]
		if colors: kwargs_single['color'] = colors[i]
		ax.plot(series.index, series.values, **kwargs_single)
	
	if labels:
		ax.legend()
	
	if title:
		ax.set_title(title)
	else: ax.set_title(f'{series_list[0].name} from {series_list[0].index[0].strftime("%m-%d-%Y")} to {series_list[0].index[-1].strftime("%m-%d-%Y")}')
	
	ax.set_ylabel(series.name)
	ax.set_xlabel('Date')

	# Format x-axis: show only month and day, and increase tick frequency
	match scale:
		case 'years':
			ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
			ax.xaxis.set_major_locator(mdates.YearLocator(base=interval))  # Set ticks at year intervals
		case 'months':
			ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
			ax.xaxis.set_major_locator(mdates.MonthLocator(interval=interval))  # Set ticks at month intervals
		case 'weeks':
			ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
			ax.xaxis.set_major_locator(mdates.WeekdayLocator(interval=interval))  # Set ticks at week intervals
		case 'days':
			ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
			ax.xaxis.set_major_locator(mdates.DayLocator(interval=interval))
		case 'hours':
			ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))
			ax.xaxis.set_major_locator(mdates.HourLocator(interval=interval))
		case 'auto':
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

def plot_nested_dict(data, **kwargs):
	"""
	Makes a figure with n subplots, where n is the number of vaiables in each location dict. Must have the same variables for each location

	Parameters:
		-- data (dict): A nested dictionary where the outer keys are categories 
						and the inner keys are variables with their corresponding values.
		-- kwargs [opt]: various kwargs to pass to matpllotlib plotting functions. 
			- interval (int): interval for x-axis tick marks, whether the scale be days, hours, etc (forwarded to plot_ts() in single-variable cases only)
			- labels (list): list of labels to use for each series, respectively, in plt.plot()
			- colors (list): list of colors to use for each series, respectively, in plt.plot()
	"""
	# Get the list of variables from the first category
	locations = list(data.keys())
	variables = list(data[locations[0]].keys())

	n = len(variables)  # Number of variables

	# if there's more than 1 variable to plot, we need to make subplots
	if n > 1:
		fig, axes = plt.subplots(n, 1, figsize=(10, n*4))  # Create subplots
		# check to see if colors (list) was passed. If not, set to None
		colors = kwargs.pop('colors', None)
		if colors is None:
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
	
	# otherwise, we just need one plot
	else:
		# returns a list of series to be plotted
		# NOTE: indexing the first position is allowable because it is assumed that there is one variable per location
		# that is why we are making a single timeseries plot
		series_list = [list(vardict.values())[0] for vardict in data.values()]
		# if 'labels' were passed in plot_nested_dict, use those. otherwise, use the default location names from the nested dict
		if 'labels' in kwargs:
			labels = kwargs.pop('labels', None)
		else: labels = locations
		# pass any other kwargs to plot_ts
		return plot_ts(series_list, labels=labels, **kwargs)

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
			gap_start_ind = series.index.get_loc(gap_start)
			# duration = series.index[i] - series.index[gap_start_ind - 1]
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

def scrapeUrlsFromHttp(url, pattern=None):
	'''
	Fetch and list file URLs from an HTTP directory listing page.
	
	Args:
	-- url (str): The URL of the directory listing page.
	--pattern (str, optional): A substring to filter file URLs. Defaults to None.
		
	Returns:
	A list of matching file URLs.
	'''
	try:
		# Get page content
		response = requests.get(url)
		response.raise_for_status()  # Raise an error for bad status codes
		
		# Parse HTML content
		soup = BeautifulSoup(response.text, "html.parser")
		
		# Find all links in the page
		file_links = [url + a["href"] for a in soup.find_all("a", href=True) if not a["href"].startswith("?")]

	except requests.RequestException as e:
		print(f"Error fetching URL: {e}")
		return []

	if pattern:
		file_links = [link for link in file_links if pattern in link]
	return file_links

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

def validate_forecast_times(start_date, end_date, reference_date=None):
	'''
	Validates and processes forecast dates to ensure compliance with model forecast cycles.
	
	Args:
	-- start_date (str, date, or datetime) [req]: The start date for which to slice data.
	-- end_date (str, date, or datetime) [req]: The end date for which to slice data.
	-- reference_date (str, date, or datetime) [opt]: The forecast reference time, i.e. the date and time cycle at which the particular forecast was launched.
		Defaults to start_date if None.

	Returns:
	-- tuple (datetime, datetime, datetime): Processed start_date, end_date, and reference_date.
	'''
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)
	
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
	
	return start_date, end_date, reference_date