import datetime as dt
from .utils import parse_to_datetime, add_units
import requests
import pandas as pd
from urllib.error import HTTPError

# Candada parameter codes 
#	Débit (m³/s) :	streamflow
# 	030424 : Pike 
# 	030425 : Rock

def get_daily(id):
	locurl = 'https://www.cehq.gouv.qc.ca/depot/historique_donnees/fichier/'+id+'_Q.txt'
	try:
		df = pd.read_csv(locurl, delimiter=r'\s{2,}', index_col= 'Date', header=19, encoding='ISO-8859-1', parse_dates=True, engine='python')
		print(f'Getting data from {locurl}')
	except Exception as e:
		print(f"{e}. {locurl}")
	return df

def get_instantaneous(id, yearlist):
	df = pd.DataFrame()
	for year in yearlist :
		locurl = 'https://www.cehq.gouv.qc.ca/depot/historique_donnees_instantanees/'+id+'_Q_'+str(year)+'.txt'
		try:
			new_year = pd.read_csv(locurl, delimiter=r'\s{2,}', index_col= 'Date', header=16, encoding='ISO-8859-1', parse_dates=True, engine='python')
			print(f'Getting data from {locurl}')
		except Exception as e:
			print(f"{e}. {locurl}")
			continue
		df = pd.concat([df,new_year])
	return df

def get_data(start_date,
			 end_date,
			 locations={'Pike':'030424','Rock':'030425'},
			 variables={'streamflow':'Débit (m³/s)'},
			 service='iv'):
	"""
	A function to download and process Canadian observational hydrology data to return nested dictionary of pandas series fore each variable, for each location.

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for which to grab Canadian Instantaneous data
	-- end_date (str, date, or datetime) [req]: the end date for which to grab Canadian Instantaneous data
	-- locations (dict) [req]: a dictionary (stationID/name:IDValue/latlong tuple) of locations to get Canadian Instantaneous data for.
	-- variables (dict) [req]: a dictionary of variables to download, where keys are user-defined variable names and values are dataset-specific variable names.
	-- service (str) [opt]: what frequency of data to get. Default is 'iv', or instantaneous data (15-min frequency). Other option is 'dv', which returns daily data. For more info, see https://www.cehq.gouv.qc.ca/hydrometrie/historique_donnees/fiche_station.asp?NoStation=030425
	
	Returns:
	Canadian Instantaneous or Daily observed streamflow data for the given stations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object.
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
		
		print(f"ACQUIRING DATA FOR STATION: {station} ({id})")
		match service:
			case 'iv':
				df = get_instantaneous(id, yearlist)
				# convert timestamp to UTC with zone suffix
				df.set_index(df.index.tz_localize('EST').tz_convert('UTC'), inplace = True)
				# trim dates to the requested range
				df = df[df.index >= start_date]
				df = df[df.index < end_date]
			case 'dv':
				df = get_daily(id)
				# trim dates to the requested range
				# Time zone info not relevant for daily data
				df = df[df.index >= start_date.replace(tzinfo=None)]
				df = df[df.index < end_date.replace(tzinfo=None)]
		
		df = df[pd.to_numeric(df['Débit (m³/s)'], errors='coerce').notnull()]	# remove all rows with non-float flow values

		# Subset to the columns of interest
		df = df[['Débit (m³/s)']]
		df.columns = [list(variables.keys())[0]]
		df.index.name = 'time'
		#print(df)

		# created nested dictionary of pd.Series for each variable for each location
		caflow_data[station] = {var:df[var].astype('float64') for var in df.columns}

		# add unit info
		# assumes variables dict will have one entry, for streamflow, discharge - whatever a ueser wants to call it
		add_units(caflow_data, {list(variables.keys())[0]:'m³/s'})

	return caflow_data

