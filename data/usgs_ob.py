import datetime as dt
from .utils import parse_to_datetime
import requests
import pandas as pd

# USGS parameter codes
# https://help.waterdata.usgs.gov/codes-and-parameters/parameters
# https://help.waterdata.usgs.gov/parameter_cd?group_cd=PHY
# Streamflow, mean. daily in cubic ft / sec: '00060',
# Streamflow, instantaneous cubic ft / sec: '00061',
# Gage Height, feet: '00065'
# Lake Elevation above NGVD, ft: '62614''

def USGSgetvars_function(id, variables, start, end, serv='iv'):
	"""
	Form url for USGS station request and return formatted dataframe for the given station

	Args:
	-- id (str) [req]: station ID to get data for
	-- variables (dictionary) : ['column name' : 'usgs var code']
	-- paramter (str) [req]: parameter code of data to get
	-- start (datetime) [req]: start datetime
	-- end (datetime) [req]: end datetime
	-- serv (str) [opt]: what USGS service to get data from. Default is instanteous values service. For more options, see https://waterservices.usgs.gov/docs/

	Returns:
	A dataframe of USGS streamflow data indexed by timestamp
	"""
	# put in while loop to ensure the request doesn't fail
	returnValue = None

	# daily values service does not accept timezones, but instantaneous values service does
	if serv == 'dv':
		start_tz = ''
		end_tz = ''
	else:
		start_tz = 'T00:00Z'
		end_tz = 'T23:59Z'
	parameter = variables[list(variables)[0]]	# extract first variable code from passed dictionary
	# for more info on how to format URL requests, see:
	# https://waterservices.usgs.gov/docs/instantaneous-values/instantaneous-values-details/#url-format
	while(returnValue is None):
		gage = requests.get(f'https://waterservices.usgs.gov/nwis/{serv}/'
							 '?format=json'
							f'&sites={id}'
								# f'&period={period}'
							f'&startDT={start.strftime("%Y-%m-%d")}{start_tz}'
							f'&endDT={(end-dt.timedelta(days=1)).strftime("%Y-%m-%d")}{end_tz}'                     
							f'&parameterCd={parameter}'
							)
		# print(gage.text)
		try:
			returnValue = gage.json()
			values = gage.json()['value']['timeSeries'][0]['values'][0]['value']
		except:
			print("USGS Observational Hydrology Data Request Failed... Will retry")
			print(gage.text)
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
	station_df =  pd.DataFrame(data={list(variables)[0]: df['value'].astype(float).values}, index=pd.to_datetime(df['dateTime'], utc=True))
	station_df.index.name = 'time'
	# print(station_df)
	return station_df

def get_data(start_date,
			 end_date,
			 locations,
			 variables={'streamflow':'00060'},
			 service='iv'):
	"""
	A function to download and process USGS observational hydrology data to return nested dictionary of pandas series fore each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab USGS data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab USGS data
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get USGS data for.
	-- variables (dict) [req]: a dictionary of variables to download, where keys are user-defined variable names and values are dataset-specific variable names.
	-- service (str) [opt]: what USGS service to get data from. Default is instanteous values service. For more options, see https://waterservices.usgs.gov/docs/
	
	Returns:
	USGS observed streamflow data for the given stations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	"""
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	# 04294000 (MS), 04292810 (J-S), 04292750 (Mill)

	# Get 90 Days Prior
	period = 'P90D'
	returnVal = {}

	# 20231211 - do not adjust passed dates to a previous day. that is a caller concern if that additional data buffer is needed.
	for station, id in locations.items():
		returnVal[station] = USGSgetvars_function(id,
												variables,
												start_date.date(),
												end_date.date(),
												service)
	
	usgs_data = {station:{name:data for name, data in station_df.T.iterrows()} for station, station_df in returnVal.items()}
	
	return usgs_data