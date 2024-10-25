from .utils import parse_to_datetime, add_units
import dask.dataframe as dd
import numpy as np
import os
import pandas as pd
import datetime as dt
import pytz

"""
Data Acquisition module for the UVM FEMC (Forest Ecosysytem Monitoring Cooperative - https://www.uvm.edu/femc/)

First timestamp in V2 dataset: 2024-07-26 20:10:00
Last timestamp in v1 BEFORE first timestamp in V2: 2024-07-26 20:00:00
"""

##### Global Variables #####
# define the timestamp in the v2 data - manually code, but should not change
# unless FEMC folks add data to the beginning of the v2 dataset
FIRST_V2_TS = '2024-07-26 20:10:00'

V1_URL = 'https://water.w3.uvm.edu/data/colchester-reef-v1.csv'
V2_URL = 'https://water.w3.uvm.edu/data/colchester-reef-v2.csv'
LATEST_URL = "https://uvm.edu/femc/MetData/ColReefQAQC/CR_V2_QAQC_latest.csv"

# to get units, go to this link: https://www.uvm.edu/femc/data/archive/project/colchester-reef-meteorological-monitoring/dataset/colchester-reef-met-site/metadata#notes
# and click on 'View XML File'. Then Ctrl + f 'units' to manually find units for each var
FEMC_UNITS = {'T2':'°C',
			'SWDOWN':'W/m²',
			'RH2':'%',
			'WSPEED':'m/s',
			'WDIR':'degrees'
}
############################

def load_femc_data(start_date, end_date):
	"""
	Loads observational meteorological data from UVM FEMC (Forest Ecosystem Monitoring Cooperative) 
	in Colchester Reef, handling v1, v2 snd latest datasets, which are split by time.

	Args:
	-- start_date (str, date, or datetime) [req]: The start date for which to load the FEMC data.
	-- end_date (str, date, or datetime) [req]: The end date for which to load the FEMC data.
	
	Returns:
	A pandas DataFrame containing the meteorological data sliced between the start and end dates. 
	If both v1 and v2 datasets are required (i.e., the date range crosses between v1 and v2), the datasets 
	are concatenated. If only v1 or v2 is required, only that dataset is returned. 
	
	The resulting DataFrame is indexed by timestamps and sorted chronologically.

	Notes:
	-- The function uses lazy loading (via Dask) for v1 data due to its large size. The v2 dataset is loaded using pandas.
	-- Assumes that data past 2024-07-26 20:10:00 is contained in the v2 dataset.
	-- Ensures the date range is valid and that the end date is not in the future.
	"""
	# parse start and end dates again so that function can be used independently of get_data()
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	# set to eastern time zone
	start_date_est = start_date.astimezone(pytz.timezone('EST'))
	end_date_est = end_date.astimezone(pytz.timezone('EST'))

	# format as timestamps for easy indexing
	start_date_est_ts = start_date_est.strftime("%Y-%m-%d %H:%M:%S")
	end_date_est_ts = end_date_est.strftime("%Y-%m-%d %H:%M:%S")

	# add eastern time zone so we can compare with start/end_date_est
	firstV2dt = dt.datetime.strptime(FIRST_V2_TS, "%Y-%m-%d %H:%M:%S").replace(tzinfo=pytz.timezone('EST'))

	# assert assumptions for that need to be true for the rest of the function's logic to work
	assert start_date_est < end_date_est, 'start date must come before end date'
	today = dt.datetime.today().astimezone(pytz.timezone('EST'))
	assert end_date_est <= today, 'end date cannot be in the future'

	# start with some None objects
	v1_sliced = None
	v2_sliced = None
	v1v2_df = None

	# need to load v1 data if start date is before beginning of v2
	if start_date_est < firstV2dt:
		v1_to_v2_colnames = {'38m_AIRTEMP':'35m_AIRTEMP',
					 '38m_RELHUMID':'35m_RELHUMID',
					 'NRG_38m_MEAN_RESULTANT_WINDSPEED':'35m_MEAN_WINDSPEED',
					 'NRG_38m_MEAN_WIND_DIRECTION':'35m_MEAN_WIND_DIRECTION'}
		# if the start date is before when the v2 dataset begins, we need the v1 dataset
		# use a dask dataframe to load lazily, since the v1 dataset is huge
		v1 = dd.read_csv(V1_URL, delimiter=",", header=0, parse_dates=True)
		v1 = v1.set_index(dd.to_datetime(v1['RECORDTIME'])).drop('RECORDTIME', axis=1)
		v1 = v1.rename(columns=v1_to_v2_colnames)
		v1_sliced_lazy = v1.loc[start_date_est_ts:end_date_est_ts]
		v1_sliced = v1_sliced_lazy.compute()
		# if the end date is before the beginning of the v2 dataset, the we just need v1
		if end_date_est < firstV2dt: 
			print("returning just v1 data")
			return v1_sliced

	# load v2 csv if start date is before the last 20 days
	if start_date_est < (today - dt.timedelta(days=20)):
		# loading in the v2 dataset with pandas since its not too big
		v2 =  pd.read_csv(V2_URL, delimiter=",", header=0, index_col=0, parse_dates=True)
		# slice v2 as we will need it in all future cases
		v2_sliced = v2.loc[start_date_est_ts:end_date_est_ts]
		### need to check for v2 cache gap here
		# concat v2 with v1 if v1 exists
		if v1_sliced is not None:
			v1v2_df = pd.concat([v2_sliced, v1_sliced[:FIRST_V2_TS]]).sort_index()
		# if the end date is before the last 20 days, we have all the data we need already
		if end_date_est <= (today - dt.timedelta(days=20)):
			if v1v2_df is not None:
				print("returning v1 and v2 concated data")
				# return v1_sliced[:FIRST_V2_TS], v2_sliced
				return v1v2_df
			else: 			
				print('returning just v2 data')
				return v2_sliced

	# Load the last 21 days from url only, IF the end date is within the last 20 days (should be all other cases)
	print(f"Loading lastest Colchester Reef Data from: {LATEST_URL}")
	latest = pd.DataFrame(pd.read_csv(LATEST_URL, delimiter=",", header=0, index_col=0, parse_dates=True))
	latest_sliced = latest.loc[start_date_est_ts:end_date_est_ts]
	# if v1 and v2 were combined, concat with the latest
	if v1v2_df is not None:
		print("concating all three datasets")
		# combine first will use indices from the first df, and any from that second that aren't in the first 
		df = v1v2_df.combine_first(latest_sliced)
		return df
		# return v1v2_df, latest
	# if v1 is None but v2 is sliced, concat with latest
	elif v2_sliced is not None:
		print('returning v2 and latest concated data')
		df = v2_sliced.combine_first(latest_sliced)
		return df
	# otherwise of v1 or v2 were not used, just return v1
	else: 
		print('returning just latest data')
		return latest_sliced

def get_data(start_date,
			 end_date,
			 locations = {'CR':'ColReefQAQC'},
			 variables = {'T2':'T2',
				 		  'SWDOWN':'SWDOWN',
						  'RH2':'RH2',
						  'WSPEED':'WSPEED',
						  'WDIR':'WDIR'}):
	"""
	A function to download and process observational meterological data from UVM FEMC (Forest Ecosysytem Monitoring Cooperative - https://www.uvm.edu/femc/) to return nested dictionary of pandas series fore each variable, for each location.
	
	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab FEMC data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab FEMC data
	-- locations (dict) [opt]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get FEMC data for.
	-- variables (dict) [opt]: a dictionary of variables to get; keys can be whatever you want to call the variables, but the values must be the variable abbreviations as seen in default dictionary
	
	Returns:
	FEMC obsrvational meterological data for the specifed data range and locations, in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
	"""
	# ensure start and end dates are datetime objects
	start_date = parse_to_datetime(start_date)

	end_date = parse_to_datetime(end_date)

	femc_data = {loc:{} for loc in locations.keys()}

	for loc, _ in femc_data.items():

		print(f'Aquiring FEMC data from {start_date.strftime("%m-%d-%Y")} to {end_date.strftime("%m-%d-%Y")}')
		print(f"Getting the following variables: {variables.keys()}")

		# load data
		df = load_femc_data(start_date, end_date)

		# set index to UTC
		df.set_index(df.index.tz_localize('EST').tz_convert('UTC'), inplace = True)

		# data returned should be start_Date includsive, end_date exclusive
		df = df[df.index >= start_date]
		df = df[df.index < end_date]

		# Subset to the columns of interest
		# these column names represent the v2 naming conventions
		df = df[['35m_AIRTEMP','PYRANOM','35m_RELHUMID', '35m_MEAN_WINDSPEED','35m_MEAN_WIND_DIRECTION']]
		df.columns = ['T2', 'SWDOWN', 'RH2', 'WSPEED', 'WDIR']
		df.index.name = 'time'

		# map units to user-variable names
		units = {var_name:FEMC_UNITS[var] for var_name, var in variables.items()}

		femc_data[loc] = {var_name:df[var] for var_name, var in variables.items()}
		
		# add units to series
		add_units(femc_data, units)
		
	return femc_data