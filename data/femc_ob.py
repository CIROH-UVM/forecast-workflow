from .utils import parse_to_datetime
import numpy as np
import os
import pandas as pd
import datetime as dt
from zoneinfo import ZoneInfo

def get_data(start_date,
			 end_date,
			 locations = {'CR':'ColReefQAQC'},
			 data_dir = '/data/forecastData/'):
	"""
	A function to download and process observational meterological data from UVM FEMC (Forest Ecosysytem Monitoring Cooperative - https://www.uvm.edu/femc/) to return nested dictionary of pandas series fore each variable, for each location.
	
	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab FEMC data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab FEMC data
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get FEMC data for.
	-- data_dir (str) [opt]: Directory in which the cached Colchester Reef CSV can be found
	
	Returns:
	FEMC obsrvational meterological data for the specifed data range and locations, in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	"""
	# ensure start and end dates are datetime objects
	start_date = parse_to_datetime(start_date)
	# subtracting a day so that 
	end_date = parse_to_datetime(end_date)

	femc_data = {loc:{} for loc in locations.keys()}

	for loc, loc_value in femc_data.items():
		with open(os.path.join(data_dir, 'colchesterReefFEMC/Z0080_CR_QAQC.csv')) as file:
			cached_data = pd.read_csv(file, delimiter=",", header=0, index_col=0, parse_dates=True) # Make DateTime as index
			cached_df = pd.DataFrame(cached_data)

			# Also get the URL for the most recent Colchester Reef data in the CSV. Concat the two frames.
			url="https://uvm.edu/femc/MetData/ColReefQAQC/CR_V2_QAQC_latest.csv"
			print(url)
			loaded_df = pd.DataFrame(pd.read_csv(url, delimiter=",", header=0, index_col=0, parse_dates=True))

			#df = pd.concat([cached_df, loaded_df], axis=0, ignore_index=True)   # ignore dupe time indexes (overlap)
			df = pd.concat([cached_df, loaded_df], axis=0)
			# Drop duplicate time indices from overlap between 2 .csvs -- And sort
			df = df[~df.index.duplicated(keep='first')].sort_index()

			# establish the default timezone as UTC
			df.set_index(df.index.tz_localize('EST').tz_convert('UTC'), inplace = True)

			# data returned should be start_Date includsive, end_date exclusive
			df = df[df.index >= start_date]
			df = df[df.index < end_date]
			
			# Subset to the columns of interest
			df = df[['38m_AIRTEMP','PYRANOM','38m_RELHUMID', 'NRG_38m_MEAN_RESULTANT_WINDSPEED','NRG_38m_MEAN_WIND_DIRECTION']]
			df.columns = ['T2', 'SWDOWN', 'RH2', 'WSPEED', 'WDIR']
			df.index.name = 'time'

		femc_data[loc] = {var:df[var] for var in df.columns}

			
	return femc_data