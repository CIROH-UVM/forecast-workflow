import datetime as dt
from .utils import parse_to_datetime
import requests
import pandas as pd

# Candada parameter codes 
#	Débit (m³/s) :	streamflow
# 	030424 : Pike 
# 	030425 : Rock



def get_data(start_date,
			 end_date,
			 locations={'Pike':'030424','Rock':'030425'},
			 variables={'streamflow':'Débit (m³/s)'}):
	"""
	A function to download and process Canadian observational hydrology data to return nested dictionary of pandas series fore each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab Canadian Instantaneous data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab Canadian Instantaneous data
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get Canadian Instantaneous data for.
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
	Canadian Instantaneous observed streamflow data for the given stations in the format specified by return_type
	"""
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)
	yearlist = list(range(start_date.year, end_date.year+1)) # the list of all years data is requested for
	print('ca_flow.get_data years to process:')
	print(yearlist)
	print('ca_flow.get_data() locations')
	print(locations.items())

	caflow_data = {}

	# 20231211 - do not adjust passed dates to a previous day. that is a caller concern if that additional data buffer is needed.
	for station, id in locations.items():

		df = pd.DataFrame()

		for year in yearlist :
			locurl = 'https://www.cehq.gouv.qc.ca/depot/historique_donnees_instantanees/'+id+'_Q_'+str(year)+'.txt'
			df = pd.concat([df,pd.read_csv(locurl, delimiter=r'\s{2,}', index_col= 'Date', header=15, encoding='ISO-8859-1', parse_dates=True, engine='python')])


		df = df[pd.to_numeric(df['Débit (m³/s)'], errors='coerce').notnull()]	# remove all rows with non-float flow values

		# convert timestamp to UTC with zone suffix
		df.set_index(df.index.tz_localize('EST').tz_convert('UTC'), inplace = True)
	
		# trim dates to the requested range
		df = df[df.index >= start_date]
		df = df[df.index < end_date]

		# Subset to the columns of interest
		df = df[['Débit (m³/s)']]
		df.columns = ['streamflow']
		df.index.name = 'time'
		#print(df)

		# created nested dictionary of pd.Series for each variable for each location
		caflow_data[station] = {var:df[var] for var in df.columns}

	return caflow_data

