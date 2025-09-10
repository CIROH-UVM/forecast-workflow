from data import utils
import datetime as dt
import numpy as np
import pandas as pd
import requests
import warnings

'''
var_units represents a dictionary of LCD variable names and their corresponding units.
Only the variables listed here have been tested and are sure to work in this module.

LCD metadata and attributes, inlcude standard and metric units, can be found at:
https://www.ncei.noaa.gov/access/services/support/v3/datasets.json
'''
global standard_var_units
standard_var_units = {'HourlySkyConditions':'%',
					  'HourlyPrecipitation':'inches',
					  'HourlyDryBulbTemperature':'\N{DEGREE SIGN}F',
					  'HourlyRelativeHumidity':'%',
					  'HourlyWindSpeed':'mph',
					  'HourlyWindDirection':'\N{DEGREE SIGN}E of N'}

global metric_var_units
metric_var_units = {'HourlySkyConditions':'%',
					'HourlyPrecipitation':'mm',
					'HourlyDryBulbTemperature':'\N{DEGREE SIGN}C',
					'HourlyRelativeHumidity':'%',
					'HourlyWindSpeed':'m/s',
					'HourlyWindDirection':'\N{DEGREE SIGN}E of N'}

'''

Var-specific processing functions:

'''
def splitsky ( instring ) :
	thestring = str(instring)
	tokens = thestring.split(':')   # looking for at least one cover code marker :
	if len(tokens) > 1 :
		token = tokens[-2].split()[-1]  # keep the cover code just before last : marker
	else:
		token = ' '     # any records without a : cover code separator end up as a blank and discarded later
	# Remove X:10 (Obscured Sky) Observations
	if token == 'X':
		token = ' '
	return token

def sky2prop (theskycode) :
  skypropmap = {'CLR': 0.000, 'FEW': 0.250, 'SCT': 0.5000, 'BKN': 0.875, 'OVC': 1.000, 'VV': 1.000, ' ': 1.000}
  theprop = skypropmap[str(theskycode)]
  return theprop

def process_clouds(cloud_series):
	cloud_series = cloud_series.apply(splitsky)
	# Remove those that don't convert to skycode... junk entries
	cloud_series = cloud_series[cloud_series != ' ']
	cloud_series = cloud_series.apply(sky2prop)
	return cloud_series

def leavenotrace (precip) :
	if str(precip) == 'T' :
		return '0.00'
	# Some suspect values (marked with s) contain junk data with 2 decimal points
	# (i.e. STN 72617014742 DATE 2020-01-16T09:54:00)
	elif len(str(precip).split('.')) > 2 :
		return 'NaN'
	else :
		#return precip
		return str(precip).replace('s', '')

def process_rain(precip_df, user_name):
	# First replace 'T's for trace precip with 0.0
	#  leavenotrace also removes 's' notations on some precip values
	#  Also, convert to float
	precip_df[user_name] = precip_df['HourlyPrecipitation'].apply(leavenotrace).astype('float')
	# Then, dump rows with NaN for RAIN
	precip_df = precip_df[~precip_df['RAIN'].isna()]
	# add units - gathered from LCD documentation
	precip_df = precip_df.assign(Units='inches')
	return precip_df

# def process_air_temp(temp_df, user_name):
# 	# create a new column with user-defined var name, and strip special indicator characters
# 	temp_df[user_name] = leavenosuspect(temp_df['HourlyDryBulbTemperature'])
# 	# uncomment the below line to convert temperature to celsius
# 	# temp_df[user_name] = fahr_to_celsius(temp_df[user_name])
# 	# drop NA's, based on corrected column
# 	temp_df = temp_df[~temp_df[user_name].isna()]
# 	# adding a units column
# 	temp_df = temp_df.assign(Units='\N{DEGREE SIGN}F')
# 	return temp_df

# def process_relhum(relhum_df, user_name):
# 	# get rid of duplicate timestamps, keep first instance of each duplicated index
# 	relhum_df = relhum_df[~relhum_df.index.duplicated(keep='first')]
# 	# remove special characters from rel hum column
# 	relhum_df[user_name] = leavenosuspect(relhum_df['HourlyRelativeHumidity'])
# 	# now drop any NA's that remain in non-duplicated timestamps
# 	relhum_df = relhum_df[~relhum_df[user_name].isna()]
# 	relhum_df = relhum_df.assign(Units='%')
# 	return relhum_df

'''
Some variables may need unique special methods to parse the data (such as splitsky() for clouds for example). Others do not.
Regardless, the general procedure for cleaning the raw data returned by the API is written partially in  clean_raw_df() and then get_data(), and is in order as follows:

1. Remove duplicate values in the 'DATE' column, comparing 'REPORT_TPYE' to decide between which duplicates to keep
2. Set the index of the df to the 'DATE' column as UTC datetime timestamps, rename to 'time'
3. Remove special indicator characters ['*', 's'] from the data
4. Cast series values as np.float64 type
5. Remove Nan's from the data
6. Add a unit the series' name
'''
def clean_raw_df(raw_df):
	'''
	1. Remove duplicate values in the 'DATE' column
	2. Set the index of the df to the 'DATE' column as UTC datetime timestamps, rename to 'time'
	'''
	# all of the duplicated Dates
	all_dup_rows = raw_df.loc[raw_df['DATE'].duplicated(keep=False)]
	# Keep FM-16 reports, as we will prefer those over FM-15 reports when they exist for duplicated timestamps, since they are specially updated reports
	# we will also keep SOD reports over SOM reports when they exist for duplicated timestamps
	'''
	FM-15 = METAR Aviation routine weather report
	FM-16 = SPECI Aviation selected special weather report
	SOD = Summary of day report from U.S. ASOS or AWOS station
	SOM = Summary of month report from U.S. ASOS or AWOS station

	Source: https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
	'''
	drop_indices = []
	for ts, dups in all_dup_rows.groupby(['DATE']):
		# print(f"for duplicate Timestamp: {ts}")
		rtypes = dups['REPORT_TYPE'].values
		# print(f"\treport types: {rtypes}")
		# If there is a FM-16 report in the duplicate pair, keep it; drop the other report
		if 'FM-16' in rtypes:
			dups = dups[dups['REPORT_TYPE'] != 'FM-16']
		# else if there's no FM-16, keep the FM-15 report
		elif 'FM-15' in rtypes:
			dups = dups[dups['REPORT_TYPE'] != 'FM-15']
		# next highest priority is FM-12
		elif 'FM-12' in rtypes:
			dups = dups[dups['REPORT_TYPE'] != 'FM-12']
		# next highest priority is SOD
		elif 'SOD  ' in rtypes:
			dups = dups[dups['REPORT_TYPE'] != 'SOD  ']
		# theoretically, duplicate groups should be in pairs, so there should be just one left aftering dropping the rows we intend to keep
		if len(dups) > 1:
			raise IndexError(f"Duplicates still detected after filtering: {dups.values}")
		# print(f"DUPS AFTER FILTERING: {dups['REPORT_TYPE']}")
		# the indices that remain after filtering correspond to the rows we want to drop
		drop_indices.append(dups.index[0])
	# print(drop_indices)
	# drop the row for the report that we don't want to keep
	no_dup_df = raw_df.drop(index=drop_indices)
	# check to ensure that no duplicates remain
	if any(no_dup_df['DATE'].duplicated(keep=False)):
		raise IndexError(f"Duplicate timestamps detected: {no_dup_df[no_dup_df['DATE'].duplicated(keep=False)]}")
	# set the index to the DATE column as datetime timestamps, rename to 'time'
	no_dup_df.set_index(pd.to_datetime(no_dup_df['DATE'], utc=True), inplace=True)
	no_dup_df.index.rename('time', inplace=True)
	return no_dup_df

def scrubSpecialChars(raw_series, inplace=False):
	'''
	Processes raw LCD data Series to strip away 's' and '*' characters. These are special indicator characters appended to weather data values
	from the Burlington Airport weather station (and LCD datasets broadly). For more info on special indicator characters, see https://www.ncei.noaa.gov/pub/data/cdo/documentation/LCD_documentation.pdf
	
	Args:
	-- raw_series (Pandas Series) [req]: the raw Series to scrub

	Returns:
	A copy of the passed series with 's' and '*' characters removed
	'''
	# print('Removing special indicator chars from raw data...')
	corrected_series = raw_series.copy()
	for i, value in enumerate(raw_series.astype(str)):
		corrected_value = ''
		# correct suspect values
		# 's' character means suspect value; keep the value, ditch the char
		if 's' in value:
			corrected_value = value.replace('s','')
			corrected_series.iloc[i] = corrected_value
		# correct missing values
		# '*' shows up by itself and has no value associated with it
		if '*' in value:
			corrected_value = value.replace('*', 'nan')
			corrected_series.iloc[i] = corrected_value
		# removing variable wind (or could replace with different value)
		if value == 'VRB':
			corrected_value = value.replace('VRB', 'nan')
			corrected_series.iloc[i] = corrected_value
		# if no correction is needed, use original value
		if corrected_value:
			# optional log message below
			# print(f'Value at index {i} corrected from {value} to {corrected_value}')
			pass
		else: corrected_series.iloc[i] = np.float64(value)
	if inplace:
		raw_series[:] = corrected_series.astype(float)
		return None
	else: return corrected_series.astype(float)


# def create_final_df(df, colToKeep, index):
# 	# print(f'Creating final df for {colToKeep}')
# 	# print(df[colToKeep])
# 	# print(df['Units'])
# 	# 20231211 - set index as datetime with timezone suffix set to UTC; data is in UTC, tz_localize() tells pandas that
# 	return pd.DataFrame(data={colToKeep: df[colToKeep].to_numpy(), 'Units':df['Units'].to_numpy()}, index=pd.DatetimeIndex(data=pd.to_datetime(df[index]), name='time')).tz_localize('UTC')

def lcdRequest(startDate, endDate, var_list, station_id, units='standard'):
	# join all requested variables by a comma for the API call
	varstring = (',').join(var_list)
	# put this in loop since this fails frequently
	returnValue = None
	# note that 'T00:00:00Z' is added to startDate (and similar appendage) to endDate in order to grab data for UTC time
	while(returnValue is None):
		requeststring = 'https://www.ncei.noaa.gov/access/services/data/v1'+\
								'?dataset=local-climatological-data'+\
								'&stations='+\
									station_id+\
								'&startDate='+\
									str(startDate.date())+'T00:00:00Z'+\
								'&endDate='+\
									str(endDate.date()-dt.timedelta(days=1))+'T23:59:00Z'+\
								'&dataTypes='+\
									varstring+\
								'&format=json'+\
								'&reportTypes=FM-15'+\
								'&units='+\
									units
		print(requeststring)
		result = requests.get(requeststring)
		# Old way to test for valid response
		# if len(result.text) > 10:
		# 	resultReceived = True
		
		# New way to test for valid response
		try:
			returnValue = result.json()
		except:
			print('NOAA Local Climatological Data Request Failed... Will retry')
			print(result.text)
	# logger.info('result.text')
	# logger.info(result.text)

	return pd.DataFrame(returnValue)

def get_data(start_date,
			  end_date,
			  locations,
			  variables,
			  units='standard'):
	"""
	A function to download and process NOAA Local Climatological Data data to return nested dictionary of pandas series for each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab LCD data.
	-- end_date (str, date, or datetime) [req]: the end date for which to grab LCD data.
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get data for.
	-- variables (dict) [req]: a dictionary of variables to download, where keys are user-defined variable names and values are LCD-specific variable names.
								Currently only tested for variables listed in global var_units. 
	-- units (str) [opt]: specifies unit convention for the data request. Options are 'standard' for standard US units, or 'metric' for metric units.
		
	Returns:
	NOAA Local Climatological Data timeseries for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	"""
	start_date = utils.parse_to_datetime(start_date)
	end_date = utils.parse_to_datetime(end_date)
	lcd_data = {loc:{} for loc in locations.keys()}
	for station, id in locations.items():
		print(f'Requesting data for {station}: {id}')
		raw_df = lcdRequest(start_date, end_date, var_list=list(variables.values()), station_id=id, units=units)
		if raw_df.empty:
			warnings.warn(f"Above API request returns empty dataframe")
		else:
			clean_df = clean_raw_df(raw_df)
			var_dict = {}
			for user_name, lcd_name in variables.items():
				print(f"Processing {lcd_name} data...")
				try:
					var_series = clean_df[lcd_name]
				except Exception as e:
					# print(e)
					warnings.warn(f"The following variable was not available for the above API request: {lcd_name}")
					continue
				match lcd_name:
					case 'HourlySkyConditions':
						var_series = process_clouds(var_series)
					case 'HourlyPrecipitation':
						var_series = var_series.apply(leavenotrace)
				# scrub special indicator characters from data as described in LCD documentation
				# casts remaining values to float64
				var_series = scrubSpecialChars(var_series)
				# wind processing requires that the series already be cast as floats
				match lcd_name:
					case 'HourlyWindDirection':
						# drop rows where wind direction = 0 (little to no wind)
						var_series = var_series[var_series != 0]
					case 'HourlyWindSpeed':
						# Remove any outliers in the data (removing wind speed greater than 300 mph seems reasonable)
						var_series = var_series[var_series<300]
				# now drop nan's from series
				# var_series.dropna(inplace=True)
				# add unit information to the series' name
				match units:
					case 'standard':
						var_units = standard_var_units
					case 'metric':
						var_units = metric_var_units
				var_series.rename(f"{user_name} ({var_units[lcd_name]})", inplace=True)
				# assign series to variable key
				var_dict[user_name] = var_series
			lcd_data[station] = var_dict
	return lcd_data