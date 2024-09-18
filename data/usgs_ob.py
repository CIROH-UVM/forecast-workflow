import datetime as dt
from .utils import parse_to_datetime, add_units
import requests
import pandas as pd

# USGS parameter codes
# https://help.waterdata.usgs.gov/codes-and-parameters/parameters
# https://help.waterdata.usgs.gov/parameter_cd?group_cd=PHY
# Streamflow, mean. daily in cubic ft / sec: '00060',
# Streamflow, instantaneous cubic ft / sec: '00061',
# Gage Height, feet: '00065'
# Lake Elevation above NGVD, ft: '62614''

def USGSgetvars_function(id, variables, start, end, service='iv'):
	"""
	Form url for USGS station request and return formatted dataframe for the given station

	Args:
	-- id (str) [req]: station ID to get data for
	-- variables (dictionary) : {'column name' : 'usgs var code'}
	-- paramter (str) [req]: parameter code of data to get
	-- start (datetime) [req]: start datetime
	-- end (datetime) [req]: end datetime
	-- service (str) [opt]: what USGS service to get data from. Default is instanteous values service. For more options, see https://waterservices.usgs.gov/docs/

	Returns:
	A dataframe of USGS streamflow data indexed by timestamp
	"""
	# daily values service does not accept timezones, but instantaneous values service does
	if service == 'dv':
		start_tz = ''
		end_tz = ''
		utc = False
	elif service == 'iv':
		start_tz = 'T00:00Z'
		end_tz = 'T23:59Z'
		utc = True
	else: raise ValueError(f'Invalid service requested: "{service}"')

	loc_data = {}

	for var, parameter in variables.items():
		returnValue = None
		# for more info on how to format URL requests, see:
		# https://waterservices.usgs.gov/docs/instantaneous-values/instantaneous-values-details/#url-format
		while(returnValue is None):
			url = ''.join(f'https://waterservices.usgs.gov/nwis/{service}/'
						'?format=json'
						f'&sites={id}'
						# f'&period={period}'
						f'&startDT={start.strftime("%Y-%m-%d")}{start_tz}'
						f'&endDT={(end-dt.timedelta(days=1)).strftime("%Y-%m-%d")}{end_tz}'                     
						f'&parameterCd={parameter}'
						)
			gage = requests.get(url)
			print(f"\tAqcuiring USGS data from: {url}")
			try:
				returnValue = gage.json()
				values = gage.json()['value']['timeSeries'][0]['values'][0]['value']
			except:
				if returnValue is None:
					print("USGS Observational Hydrology Data Request Failed... Will retry")
				else: raise ValueError(f"Bad request... ensure the data you are requesting is available from station: {id}")
			# get the units for the var from the API
			unit = gage.json()['value']['timeSeries'][0]['variable']['unit']['unitCode']
		var_name = f'{var} ({unit})'

		df = pd.DataFrame(values)
		# localize timestamps to utc time zone IFF they instantaneous data was collected
		station_df =  pd.DataFrame(data={var: df['value'].astype(float).values}, index=pd.to_datetime(df['dateTime'], utc=utc))
		station_df.index.name = 'time'

		loc_data[var] = station_df[var].rename(var_name)

	return loc_data

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
	usgs_data = {}

	# 20231211 - do not adjust passed dates to a previous day. that is a caller concern if that additional data buffer is needed.
	for station, id in locations.items():
		print(f'Station: {station} ({id})')
		usgs_data[station] = USGSgetvars_function(id,
												variables,
												start_date.date(),
												end_date.date(),
												service)
	
	return usgs_data