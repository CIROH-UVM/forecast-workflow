import datetime as dt
import requests
import pandas as pd
from .utils import parse_to_datetime, fahr_to_celsius

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

def leavenosuspect(raw_temp):
	'''
	Processes air temp Series to strip away 's' and '*' characters from the data. These are special indicator characters appended to some air temp values
	from the Burlington Airport weather station. For more info on special indicator characters, see https://www.ncei.noaa.gov/pub/data/cdo/documentation/LCD_documentation.pdf
	
	Args:
	-- raw_temp (Pandas Series) [req]: the raw temp Series to parse

	Returns:
	A copy of the passed series with 's' and '*' characters removed
	'''
	print('Removing special indicator chars from raw data...')
	corrected_temp = raw_temp.copy()
	for i, value in enumerate(raw_temp.astype(str)):
		corrected_value = ''
		# correct suspect values
		if 's' in value:
			corrected_value = value.replace('s','')
			corrected_temp[i] = corrected_value
		# correct missing values
		if '*' in value:
			corrected_value = value.replace('*', 'nan')
			corrected_temp[i] = corrected_value
		# if no correction is needed, use original value
		if corrected_value:
			print(f'Value at index {i} corrected from {value} to {corrected_value}')
		else: corrected_temp[i] = value
	return corrected_temp.astype(float)

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

def create_final_df(df, colToKeep, index):
	# print(f'Creating final df for {colToKeep}')
	# print(df[colToKeep])
	# print(df['Units'])
	# 20231211 - set index as datetime with timezone suffix set to UTC; data is in UTC, tz_localize() tells pandas that
	return pd.DataFrame(data={colToKeep: df[colToKeep].to_numpy(), 'Units':df['Units'].to_numpy()}, index=pd.DatetimeIndex(data=pd.to_datetime(df[index]), name='time')).tz_localize('UTC')

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

def process_clouds(cloud_df, user_name):
	cloud_df['skycode'] = cloud_df['HourlySkyConditions'].apply(splitsky)
	# Remove those that don't convert to skycode... junk entries
	cloud_df = cloud_df[cloud_df['skycode'] != ' ']
	cloud_df[user_name] = cloud_df['skycode'].apply(sky2prop).astype('float')
	# add units
	cloud_df = cloud_df.assign(Units='%')
	return cloud_df

def process_air_temp(temp_df, user_name):
	# create a new column with user-defined var name, and strip special indicator characters
	temp_df[user_name] = leavenosuspect(temp_df['HourlyDryBulbTemperature'])
	# uncomment the below line to convert temperature to celsius
	# temp_df[user_name] = fahr_to_celsius(temp_df[user_name])
	# drop NA's, based on corrected column
	temp_df = temp_df[~temp_df[user_name].isna()]
	# adding a units column
	temp_df = temp_df.assign(Units='\N{DEGREE SIGN}F')
	return temp_df

def process_relhum(relhum_df, user_name):
	# get rid of duplicate timestamps, keep first instance of each duplicated index
	relhum_df = relhum_df[~relhum_df.index.duplicated(keep='first')]
	# remove special characters from rel hum column
	relhum_df[user_name] = leavenosuspect(relhum_df['HourlyRelativeHumidity'])
	# now drop any NA's that remain in non-duplicated timestamps
	relhum_df = relhum_df[~relhum_df[user_name].isna()]
	relhum_df = relhum_df.assign(Units='%')
	return relhum_df

def retrieve_data(startDate, endDate, variable, station_id):
	# put this in loop since this fails frequently
	returnValue = None
	# note that 'T00:00:00Z' is added to startDate (and similar appendage) to endDate in order to grab data for UTC time
	while(returnValue is None):
		requeststring = 'https://www.ncei.noaa.gov/access/services/data/v1/'+\
								'?dataset=local-climatological-data'+\
								'&stations='+\
									station_id+\
								'&startDate='+\
									str(startDate)+'T00:00:00Z'+\
								'&endDate='+\
									str(endDate-dt.timedelta(days=1))+'T23:59:00Z'+\
								'&dataTypes='+\
									variable+\
								'&format=json' 
		print(requeststring)
		result = requests.get(requeststring)
		# Old way to test for valid response
		# if len(result.text) > 10:
		# 	resultReceived = True
		
		# New way to test for valid response
		try:
			returnValue = result.json()
		except:
			print("NOAA Local Climatological Data Request Failed... Will retry")
			print(result.text)
	# logger.info('result.text')
	# logger.info(result.text)

	return pd.DataFrame(returnValue)
	

def get_data(start_date,
			 end_date,
			 locations,
			 variables):
	"""
	A function to download and process NOAA Local Climatological Data data to return nested dictionary of pandas series for each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab LCD data.
	-- end_date (str, date, or datetime) [req]: the end date for which to grab LCD data.
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get USGS data for.
	-- variables (dict) [req]: a dictionary of variables to download, where keys are user-defined variable names and values are LCD-specific variable names.
								Currently only works for LCD variables HourlyPrecipitation, HourlySkyCondtions, and HourlyDryBulbTemperature. 
		
	Returns:
	NOAA Local Climatological Data timeseries for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	"""

	# 20231211 - Do not adjust passed dates to a previous day, that's a caller concern if that data buffer is needed
	start_date = parse_to_datetime(start_date).date()
	end_date = parse_to_datetime(end_date).date()
	lcd_data = {loc:{} for loc in locations.keys()}

	for station, id in locations.items():
		# define an empty return dict for each station
		returnDict = {}
		# define dicionaries to hold raw and processed lcd data
		raw_data = {}
		processed_data = {}
		# retrieve all of the data in the variables dictionary
		for user_name, lcd_name in variables.items():
			raw_data[user_name] = retrieve_data(start_date, end_date, lcd_name, id)
			# keep first instance of any duplicates in the raw data
			print('Removing duplicate timestamps from raw data, keeping first instance of each duplication...')
			raw_data[user_name] = raw_data[user_name].loc[~raw_data[user_name]['DATE'].duplicated(keep='first')]
			try:
				print(f"{lcd_name} for station ID {id}")
				print(raw_data[user_name])
			except Exception as e:
				print(f'{type(e)}:{e}')
			# add more cases as more variables come in and need specific processing
			match lcd_name:
				case 'HourlySkyConditions':
					processed_data[user_name] = process_clouds(raw_data[user_name], user_name)
				case 'HourlyPrecipitation':
					processed_data[user_name] = process_rain(raw_data[user_name], user_name)
				case 'HourlyDryBulbTemperature':
					processed_data[user_name] = process_air_temp(raw_data[user_name], user_name)
				case 'HourlyRelativeHumidity':
					processed_data[user_name] = process_relhum(raw_data[user_name], user_name)
			
			returnDict[user_name] = create_final_df(processed_data[user_name], user_name, 'DATE')
		# return processed_data
		# created nested dictionary of pd.Series for each variable for each location
		lcd_data[station] = {var:df[var].rename(f'{var} ({df["Units"].iloc[0]})') for var, df in returnDict.items()}

	return lcd_data