import requests
import json
import pandas as pd
import numpy as np
import datetime as dt
from datetime import date
from lib import *

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
	return pd.DataFrame(data={colToKeep: df[colToKeep].to_numpy()}, index=pd.DatetimeIndex(data=pd.to_datetime(df[index]), name='time'))

def retrieve_data(startDate, endDate, variable):
	# put this in loop since this fails frequently
	resultReceived = False
	while(not resultReceived):
		requeststring = 'https://www.ncei.noaa.gov/access/services/data/v1/'+\
								'?dataset=local-climatological-data'+\
								'&stations=72617014742'+\
								'&startDate='+\
									str(startDate)+\
								'&endDate='+\
									str(endDate)+\
								'&dataTypes='+\
									variable+\
								'&format=json' 
		print(requeststring)
		result = requests.get(requeststring)
		if len(result.text) > 10:
			resultReceived = True
	
	# logger.info('result.text')
	# logger.info(result.text)

	return pd.DataFrame(result.json())
	

def get_data (ForecastStartDate, SpinupStartDate) :

		# endday = date.today()
		#d = datetime.timedelta(days = 90)
		#startday = endday - d
		endday = ForecastStartDate - dt.timedelta(days=1)
		startday = SpinupStartDate - dt.timedelta(days=1)

		# requeststring = 'https://www.ncei.noaa.gov/access/services/data/v1/'+\
		#                         '?dataset=local-climatological-data'+\
		#                         '&stations=72617014742'+\
		#                         '&startDate='+\
		#                          str(startday)+\
		#                         '&endDate='+\
		#                          str(endday)+\
		#                         '&dataTypes=HourlyPrecipitation,HourlySkyConditions'+\
		#                         '&format=json' 
		# print(requeststring)
		# result = requests.get(requeststring)

		# df = pd.DataFrame(result.json())
		
		cloud_df = retrieve_data(startday, endday, 'HourlySkyConditions')
		precip_df = retrieve_data(startday, endday, 'HourlyPrecipitation')
		
		# logger.info('cloud_df in btv_met')
		# logger.info(cloud_df)
		# logger.info('precip_df in btv_met')		
		# logger.info(precip_df)
		
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

		return returnDict
