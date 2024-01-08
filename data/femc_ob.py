from lib import parse_to_datetime
import os
import pandas as pd
import datetime as dt

def get_data(start_date,
			 end_date,
			 locations = {'CR':'ColReefQAQC'},
			 return_type='dict'):
	"""
	A function to download and process observational meterological data from UVM FEMC (Forest Ecosysytem Monitoring Cooperative - https://www.uvm.edu/femc/) to return nested dictionary of pandas series for each variable, for each location.
	
	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab FEMC data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab FEMC data
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get FEMC data for.
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
	FEMC obsrvational meterological data for the specifed data range and locations, in the format specified by return_type
	"""
	# ensure start and end dates are datetime objects
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)
	femc_data = {loc:{} for loc in locations.keys()}

	for loc, loc_value in femc_data.items():
		with open(os.path.join("/data/forecastData/colchesterReefFEMC", "Z0080_CR_QAQC.csv")) as file:
			cached_data = pd.read_csv(file, delimiter=",", header=0, index_col=0, parse_dates=True) # Make DateTime as index
			cached_df = pd.DataFrame(cached_data)

			# Also get the URL for the most recent Colchester Reef data in the CSV. Concat the two frames.
			url="https://uvm.edu/femc/MetData/ColReefQAQC/CR_QAQC_latest.csv"
			loaded_df = pd.DataFrame(pd.read_csv(url, delimiter=",", header=0, index_col=0, parse_dates=True))

			#print(loaded_df)

			#df = pd.concat([cached_df, loaded_df], axis=0, ignore_index=True)   # ignore dupe time indexes (overlap)
			df = pd.concat([cached_df, loaded_df], axis=0)
			# Drop duplicate time indices from overlap between 2 .csvs -- And sort
			df = df[~df.index.duplicated(keep='first')].sort_index()

			# Before: Keep the last 90 days (24h*4quarters*90days=8640) plus buffer
			# df = df.tail(8800)
			# Now: Keep the days from start_date to present

			# 12231211 - remove shift of start date to the previous day. That's an AEM3D specific concern
			df = df[df.index > dt.datetime.combine(start_date, dt.datetime.min.time())]

			# Drop rows after midnight today to prevent duplicate entries with forecast
			df = df[df.index < dt.datetime.combine(end_date, dt.datetime.min.time())]

			# generate index as datetime with timezone suffix in UTC time
			df.set_index(df.index.tz_localize('EST').tz_convert('UTC'), inplace = True)
			
			# Subset to the columns of interest
			df = df[['38m_AIRTEMP','PYRANOM','38m_RELHUMID', 'NRG_38m_MEAN_RESULTANT_WINDSPEED','NRG_38m_MEAN_WIND_DIRECTION']]
			df.columns = ['T2', 'SWDOWN', 'RH2', 'WSPEED', 'WDIR']
			df.index.name = 'time'

		# ensure return_type is a valid value
		if return_type not in ['dict', 'dataframe']:
			raise ValueError(f"'{return_type}' is not a valid return_type. Please use 'dict' or 'dataframe'")
		elif return_type == 'dict':
			# created nested dictionary of pd.Series for each variable for each location
			femc_data[loc] = {var:df[var] for var in df.columns}
		elif return_type == 'dataframe':
			raise Exception("'dataframe' option not implemented yet. Please use return_type = 'dict'")
			
	return femc_data
