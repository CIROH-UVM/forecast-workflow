import datetime as dt
import requests
import pandas as pd
from lib import parse_to_datetime, logger

def splitsky ( instring ) :
	thestring = str(instring)
	tokens = thestring.split(':')   # looking for at least one cover code marker :
	if len(tokens) > 1 :
		token = tokens[-2].split()[-1]  # keep the cover code just before last : marker
	else:
		token = " "     # any records without a : cover code separator end up as a blank
	return token

def sky2prop (theskycode) :
  skypropmap = {'CLR': 0.000, 'FEW': 0.250, 'SCT': 0.5000, 'BKN': 0.875, 'OVC': 1.000, 'VV': 1.000, ' ': 1.000}
  theprop = skypropmap[str(theskycode)]
  return theprop

def leavenotrace (precip) :
	if str(precip) == 'T' :
		return '0.00'
	else :
		#return precip
		return str(precip).replace('s', '')

def create_final_df(df, colToKeep, index):
	# print(f'Creating final df for {colToKeep}')
	# print(df[colToKeep])
	
	# 20231211 - set index as datetime with timezone suffix set to UTC; data is in UTC, tz_localize() tells pandas that
	return pd.DataFrame(data={colToKeep: df[colToKeep].to_numpy()}, index=pd.DatetimeIndex(data=pd.to_datetime(df[index]), name='time')).tz_localize('UTC')

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
			 locations={"BTV":"72617014742"},
			 return_type='dict'):
	"""
	A function to download and process NOAA Local Climatological Data data to return nested dictionary of pandas series for each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab LCD data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab LCD data
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get USGS data for.
	-- return_type (string) [opt]: string indicating which format to return data in. Default is "dict", which will return data in a nested dict format:
									{locationID1:{
										var1_name:pd.Series,
										var2_name:pd.Series,
										...},
									locationID2:{...},
									...
									}
									Alternative return type is "dataframe", which smashes all data into a single dataframe muliIndex'd by station ID, then timestamp
	
	Returns:
	NOAA Local Climatological Data (total cloud cover and precipitation currently) for the dat range and locations provided
	"""
	# end_date = date.today()
	#d = datetime.timedelta(days = 90)
	#start_date = end_date - d
	# 20231211 - Do not adjust passed dates to a previous day, that's a caller concern if that data buffer is needed
	start_date = parse_to_datetime(start_date).date()
	end_date = parse_to_datetime(end_date).date()
	lcd_data = {loc:{} for loc in locations.keys()}
	
	# requeststring = 'https://www.ncei.noaa.gov/access/services/data/v1/'+\
	#                         '?dataset=local-climatological-data'+\
	#                         '&stations=72617014742'+\
	#                         '&startDate='+\
	#                          str(start_date)+\
	#                         '&endDate='+\
	#                          str(end_date)+\
	#                         '&dataTypes=HourlyPrecipitation,HourlySkyConditions'+\
	#                         '&format=json' 
	# print(requeststring)
	# result = requests.get(requeststring)

	# df = pd.DataFrame(result.json())
	
	for station, id in locations.items():
		cloud_df = retrieve_data(start_date, end_date, 'HourlySkyConditions', id)
		precip_df = retrieve_data(start_date, end_date, 'HourlyPrecipitation', id)
		
		try:
			logger.info('cloud_df in btv_met')
			logger.info(cloud_df)
			logger.info('precip_df in btv_met')		
			logger.info(precip_df)
		except Exception as e:
			print(f'{type(e)}:{e}')
		
		returnDict = {}

		cloud_df['skycode'] = cloud_df['HourlySkyConditions'].apply(splitsky)
		cloud_df['TCDC'] = cloud_df['skycode'].apply(sky2prop).astype('float')
		# Remove those that don't convert to skycode... junk entries
		cloud_df = cloud_df[cloud_df['skycode'] != ' ']
		
		# First replace 'T's for trace precip with 0.0
		#  leavenotrace also removes 's' notations on some precip values
		#  Also, convert to float
		precip_df['RAIN'] = precip_df['HourlyPrecipitation'].apply(leavenotrace).astype('float')
		# Then, dump rows with NaN for RAIN
		precip_df = precip_df[~precip_df['RAIN'].isna()]

		returnDict['TCDC'] = create_final_df(cloud_df, 'TCDC', 'DATE')
		returnDict['RAIN'] = create_final_df(precip_df, 'RAIN', 'DATE')

		# ensure return_type is a valid value
		if return_type not in ['dict', 'dataframe']:
			raise ValueError(f"'{return_type}' is not a valid return_type. Please use 'dict' or 'dataframe'")
		elif return_type == 'dict':
			# created nested dictionary of pd.Series for each variable for each location
			lcd_data[station] = {var:df[var] for var, df in returnDict.items()}
		elif return_type == 'dataframe':
			raise Exception("'dataframe' option not implemented yet. Please use return_type = 'dict'")
	return lcd_data
