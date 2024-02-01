from lib import parse_to_datetime
import requests
import pandas as pd

# USGS parameter codes
# https://help.waterdata.usgs.gov/codes-and-parameters/parameters
# https://help.waterdata.usgs.gov/parameter_cd?group_cd=PHY
# Streamflow, mean. daily in cubic ft / sec: '00060',
# Streamflow, instantaneous cubic ft / sec: '00061',
# Gage Height, feet: '00065'

def USGSstreamflow_function(id, parameter, start, end):
	"""
	Form url for USGS station request and return formatted dataframe for the given station

	Args:
	-- id (str) [req]: station ID to get data for
	-- paramter (str) [req]: parameter code of data to get
	-- start (datetime) [req]: start datetime
	-- end (datetime) [req]: end datetime

	Returns:
	A dataframe of USGS streamflow data indexed by timestamp
	"""
	# put in while loop to ensure the request doesn't fail
	returnValue = None
	# for more info on how to format URL requests, see:
	# https://waterservices.usgs.gov/docs/instantaneous-values/instantaneous-values-details/#url-format
	while(returnValue is None):
		gage = requests.get('https://waterservices.usgs.gov/nwis/iv/'
							 '?format=json'
							f'&sites={id}'
#                   		 f'&period={period}'
							f'&startDT={start.strftime("%Y-%m-%d")}T00:00Z'
							f'&endDT={end.strftime("%Y-%m-%d")}T23:59Z'                     
							f'&parameterCd={parameter}'
							)
		# print(gage.text)
		try:
			returnValue = gage.json()
			values = gage.json()['value']['timeSeries'][0]['values'][0]['value']
		except:
			print("USGS Observational Hydrology Data Request Failed... Will retry")
			# print(gage.text)
	df = pd.DataFrame(values)
	# notice timezone is localized to UTC
	# print(df['dateTime'])
	# df = df.set_index('dateTime')
	# df = df.set_index(pd.to_datetime(df['dateTime']))
	# df = df.drop(['dateTime','qualifiers'],axis =1)
	# df.columns = ['streamflow']
	# df.to_csv(id+"_flow.csv", sep=',')
	# 'US/Eastern' is the other option, but what about fall daylight savings "fall back"
	#return pd.DataFrame(data={'streamflow': df['value'].values}, index=pd.to_datetime(df['dateTime'], utc=True).dt.tz_convert('Etc/GMT+4').dt.tz_localize(None))
	# 20231211 - set index as datetime with timezone suffix set to UTC
	# utc=True converts the dateTime column to UTC, since dateTime is already UTC-localized
	station_df =  pd.DataFrame(data={'streamflow': df['value'].astype(float).values}, index=pd.to_datetime(df['dateTime'], utc=True))
	station_df.index.name = 'time'
	# print(station_df)
	return station_df

def get_data(start_date,
			 end_date,
			 locations,
			 return_type='dict'):
	"""
	A function to download and process USGS observational hydrology data to return nested dictionary of pandas series fore each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab USGS data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab USGS data
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
	USGS observed streamflow data for the given stations in the format specified by return_type
	"""
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	# 04294000 (MS), 04292810 (J-S), 04292750 (Mill)
	parameter = '00060'
	# Get 90 Days Prior
	period = 'P90D'
	returnVal = {}

	# 20231211 - do not adjust passed dates to a previous day. that is a caller concern if that additional data buffer is needed.
	for station, id in locations.items():
		returnVal[station] = USGSstreamflow_function(id,
												parameter,
												start_date,
												end_date)
	
		# ensure return_type is a valid value
		if return_type not in ['dict', 'dataframe']:
			raise ValueError(f"'{return_type}' is not a valid return_type. Please use 'dict' or 'dataframe'")
		elif return_type == 'dict':
			# created nested dictionary of pd.Series for each variable for each location
			usgs_data = {station:{name:data for name, data in station_df.T.iterrows()} for station, station_df in returnVal.items()}
		elif return_type == 'dataframe':
			raise Exception("'dataframe' option not implemented yet. Please use return_type = 'dict'")

	return usgs_data

