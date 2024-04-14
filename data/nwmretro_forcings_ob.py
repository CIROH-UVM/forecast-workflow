import datetime as dt
from kerchunk.combine import MultiZarrToZarr
from more_itertools import chunked
import numpy as np
import pandas as pd
import pyproj
from s3fs import S3FileSystem
import xarray as xr
# from lib import logger

# stick with print statements
# Question: In terms of best practicee, should we use print() or logger.info() statements for log messages

# Data source: https://ciroh-nwm-zarr-retrospective-data-copy.s3.amazonaws.com/index.html


def get_data(start_date,
			 end_date,
			 locations,
			 variables,
			 return_type='dict'):
	'''
	A function to download and process NWM Retro Forcings data to return nested dictionary of pandas series for each variable, for each location.

		Args:
		-- start_date (str, date, or datetime) [req]: the start date for which to grab data
		-- end_date (str, date, or datetime) [req]: the end date for which to grab data
		-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get data for.
		-- variables (dict) [req]: a dictionary of variables to download.
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
		
	'''
	pass
	return

