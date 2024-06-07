#  Creates the AEM3D Lake Model Input Contents
#
#  Time Series
#       Flow from RHESSys or SWAT
#       Weather from WRF
#         Net Longwave Radiation
#         Precipitation
#         Temperature
#		  Wind
#         Humidity
#       Lake Level (flow and temp related)
#       Salinity
#       Tracers
#       Call waterquality.py to wq files
#
#  Control File

import copy
from lib import *
from data import (femc_ob,
				  nwm_fc, 
				  usgs_ob,
				  gfs_fc_thredds,
				  lcd_ob,
				  caflow_ob
)

from data.utils import combine_timeseries

from .waterquality import *
from .AEM3D import *

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import datetime as dt

AEM3D_DEL_T = 300

# initialize empty dict to store information for each subplot to be made
SUBPLOT_PACKAGES = {}

# set up 2x5 figure to add subplots to
FIG, AXES = plt.subplots(nrows=2, ncols=5, figsize=(20, 10), sharex=True, layout='constrained')
FIG.supxlabel('Datetime')
# fig.suptitle()
# Rotate x-axis tick labels by 45 degrees
for i in range(0,5):
	for label in AXES[1,i].get_xticklabels():
		label.set_rotation(45)
# FIG.delaxes(AXES[0,4])

def get_climate_zone_keys(dict):
	return [zone for zone in dict.keys() if int(int(zone) / 100) == 4]

def add_plot(labelled_data, ylabel, title, fc_start, row, col, axis):
	ax = axis[row,col]
	# Plot data for location 1
	for lab, data in labelled_data.items():
		time_sliced_data = data[data.index <= pd.Timestamp(fc_start + dt.timedelta(days=7))]
		time_sliced_data = time_sliced_data[time_sliced_data.index > pd.Timestamp(fc_start - dt.timedelta(days=14))]
		ax.plot(time_sliced_data.index, time_sliced_data, label=lab)

	# Add labels and title
	# ax.set_xlabel('Datetime')
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	ax.grid(True)
	ax.legend()

	# Adding vert line at point where data goes from USGS to NWM
	ax.axvline(x=pd.Timestamp(fc_start), color='r', linestyle='--')

def make_figure(plot_package):
	for kwargs in plot_package.values():
		add_plot(**kwargs)

def print_df(df):
	logger.info('\n'
				f'{df}\n'
				f'column_dtypes:\n{df.dtypes}\n'
				f'index_dtype: {df.index.dtype}\n'
				)

def colsToHeader(columns):
	return_string = ''
	for col in columns:
		return_string += ' ' + str(col)
	return return_string

def remove_nas(series):
	# logger.info('Size with nas')
	# logger.info(len(series))
	# new_series = series[~series.isna()]
	# logger.info('Size without nas')
	# logger.info(len(new_series))    
	# return new_series
	return series[~series.isna()]

def ordinalnudgerow(rowtonudge, columntonudge, nudgeframe):
	# apply a proportional nudge value that is specific to the day of year in the passed row
	#   rowtonudge - row from a dataframe with TIME column and a value column
	#   columntonudge - the specific column name to apply the data adjustement to
	#   nudgeframe - dataframe with days of year and proportions (NudgeObs) to apply to the data that needs nudging

	DOY = int(rowtonudge['TIME'].split('.')[0][-3:])    # pull the day of year from ordinal date
	nudged = rowtonudge[columntonudge] * nudgeframe.loc[nudgeframe['Day of Year'] == DOY]['NudgeObs']
	return nudged.reset_index(drop=True)

def writeFile(filename, bayid, zone, varName, dataSeries):
	with open(filename, mode='w', newline='') as output_file:

		# output the header text
		output_file.write('!-----------------------------------------------------!\n')
		output_file.write('! Written by AEM3D_prep_IAM                           !\n')
		output_file.write(f'! Bay ID: {bayid}                                 !\n')
		#output_file.write('! Bay Source: ' + bs_name + '                         !\n')
		output_file.write('!-----------------------------------------------------!\n')
		output_file.write('1 data sets\n')
		output_file.write('0 seconds between data\n')
		output_file.write(f'0            {zone}\n')
		output_file.write(f'TIME         {varName}\n')

		dataSeries.to_csv(path_or_buf = output_file, float_format='%.3f', sep=' ', index=True, header=False)

'''
def writeLongwaveRadiationDownward(climate, THEBAY):
	###########################################################################################
	#
	#   Longwave Radiation Downward : GLW
	#

	for zone in climate['AEMLW'].keys():
		filename = f'LWRADIN_{zone}.dat'
		logger.info('Generating Bay Longwave Radiation Downward File: '+filename)
		writeFile(
			os.path.join(THEBAY.infile_dir, filename),
			THEBAY.bayid,
			zone,
			"LW_RAD_IN",
			index_to_ordinal_date(climate['AEMLW'][zone]))
		THEBAY.addfile(fname=filename)
'''

'''
def writeCloudCover(climate, THEBAY):
	###########################################################################################
	#
	#   Total Cloud Cover Entire Atmosphere : TCDC
	#

	for zone in climate['AEMLW'].keys():
		filename = f'CLOUDS_{zone}.dat'
		logger.info('Generating Bay Cloud Cover File: '+filename)
		writeFile(
			os.path.join(THEBAY.infile_dir, filename),
			THEBAY.bayid,
			zone,
			"CLOUDS",
			index_to_ordinal_date(climate['AEMLW'][zone]))
		THEBAY.addfile(fname=filename)
'''

#########################################################################
#
#   Import Hydrology Flow - Generate Flow Files
#
##########################################################################


def getflowfiles(whichbay, settings):
	'''
	getflowfiles : Get hydrology model flow for Bay Inflow
		Most information contained in passed Bay Object
	'''

	THEBAY = whichbay
	logger.info('Loading Hydrology Flow Data')

	observedHydro = None
	if settings['hydrology_dataset_observed'] == 'USGS_IV':
		observedHydro = usgs_ob.get_data(start_date = settings['spinup_date'] - dt.timedelta(days=1),
								 		 end_date = settings['forecast_start'],
										 locations = {"MS":'04294000',
			   									 	  "J-S":'04292810',
			   										  "Mill":'04292750'})
		# Convert USGS streamflow from cubic ft / s to cubic m / s
		for location in observedHydro.keys():
			observedHydro[location]['streamflow'] = observedHydro[location]['streamflow'] * 0.0283168

	else:
		raise ValueError(f"'{settings['hydrology_dataset_observed']}' is not a valid observational hydrology dataset")
	

	# pull observed Rock and Pike flow data from candadian service site
	observedca = caflow_ob.get_data(start_date = settings['spinup_date'] - dt.timedelta(days=1),
								 		 end_date = settings['forecast_start'],
										 locations = {"Rock":'030424',
			   										  "Pike":'030425'})

	forecastHydro = None
	if settings['hydrology_dataset_forecast'] == 'USGS_IV':
		forecastHydro = usgs_ob.get_data(start_date = settings['forecast_start'],
								 		 end_date = settings['forecast_end'] + dt.timedelta(days=1),
										 locations = {"MS":'04294000',
			   									 	  "J-S":'04292810',
			   										  "Mill":'04292750'})
		# Convert USGS streamflow from cubic ft / s to cubic m / s
		for location in forecastHydro.keys():
			forecastHydro[location]['streamflow'] = forecastHydro[location]['streamflow'] * 0.0283168

		forecastca = caflow_ob.get_data(start_date = settings['forecast_start'],
								 		end_date = settings['forecast_end'] + dt.timedelta(days=1),
										locations = {"Rock":'030424',
			   										 "Pike":'030425'})

	
	elif(settings['hydrology_dataset_forecast'] == 'NOAA_NWM_PROD'):
		forecastHydro = nwm_fc.get_data(forecast_datetime = settings['forecast_start'],
						end_datetime = settings['forecast_end'] + dt.timedelta(days=1),
						locations = {"MS":"166176984",
								     "J-S":"4587092",
								     "Mill":"4587100"},
						forecast_type = settings['nwm_forecast_member'],
						data_dir = settings['data_dir'],
						load_threads = 1,
						google_buckets = True)

		#############  speculative coding ahead
		#	approximate Rock and Pike forecast from scaled Missisquoi forecast
		#forecastca['Rock'] = forecastHydro['MS']['streamflow'] * THEBAY.sourcemap['21']['prop']
		#forecastca['Pike'] = forecastHydro['MS']['streamflow'] * THEBAY.sourcemap['22']['prop']
		####  need to tweak bay object sourcemap so that later the Rock and Pike are adjusted from concat'ed flows

	else:
		raise ValueError(f"'{settings['hydrology_dataset_forecast']}' is not a valid hydrology forecast dataset")
			
	
	# Need to adjust for column names
	flowdf = pd.concat([observedHydro['MS']['streamflow'], forecastHydro['MS']['streamflow']]).rename_axis('time').astype('float').to_frame()
	mlflow = pd.concat([observedHydro['Mill']['streamflow'], forecastHydro['Mill']['streamflow']]).rename_axis('time').astype('float').to_frame()
	jsflow = pd.concat([observedHydro['J-S']['streamflow'], forecastHydro['J-S']['streamflow']]).rename_axis('time').astype('float').to_frame()
	# TODO add the concats for the pike and rock 


	# these df's may not be the same length, there may be missing data from USGS gauges
	logger.info(flowdf)
	logger.info(mlflow)
	logger.info(jsflow)

	# sometimes there are negative sensor readings, let's get rid of those if they're there
	flowdf =  flowdf[flowdf['streamflow'] >= 0].reindex(flowdf.index, method='nearest')
	mlflow =  mlflow[mlflow['streamflow'] >= 0].reindex(mlflow.index, method='nearest')
	jsflow =  jsflow[jsflow['streamflow'] >= 0].reindex(jsflow.index, method='nearest')

	flow_data = {'Missisquoi':flowdf['streamflow'],
				 'Mill':mlflow['streamflow'],
				 'Jewett-Stevens':jsflow['streamflow']}

	global SUBPLOT_PACKAGES
	global AXES
	SUBPLOT_PACKAGES['streamflow'] = {'labelled_data':flow_data,
								   	  'ylabel':'Streamflow (m/s^3)',
									  'title':f"{settings['hydrology_dataset_observed']} vs. {settings['hydrology_dataset_forecast']}",
									  'fc_start':settings['forecast_start'],
									  'row':0,
									  'col':0,
									  'axis':AXES}


	flowdf.columns = ['msflow']
	## Have to merge the other two because USGS sometimes doesn't return all dates for all gauges
	jsflow.columns = ['jsflow']
	flowdf = pd.merge(flowdf, jsflow, on='time', how='inner')
	mlflow.columns = ['mlflow']
	flowdf = pd.merge(flowdf, mlflow, on='time', how='inner')

	flowdf['ordinaldate'] = flowdf.index.to_series().apply(datetimeToOrdinal)

	logger.info(flowdf)

	# Store flow series in bay object for later
	THEBAY.flowdf = flowdf[['ordinaldate', 'msflow', 'mlflow', 'jsflow']].copy()

	logger.info('Daily Flow Data Scaled')
	# Scale Additional Inflows from Predicted Inflow
	#       sourcelist has list of water source IDs
	#       sourcemap defines the source name (wshed) and proportion of hydromodel output flow
	#		the adjust value is defined from InlandSeaModel_Notes Document describing how model calibration was done


	# TODO - rework at least the Rock and Pike handling to work off the series built above rather than a proportion of Missisquoi
	#		better yet, use series for all the sources for consistent coding.   series should all be local as built just above
	for baysource in THEBAY.sourcelist :
		logger.info('Generating Bay Source File for Id: '+baysource)
		bs_prop = THEBAY.sourcemap[baysource]['prop'] * THEBAY.sourcemap[baysource]['adjust']  # proportion of input file for this source
		wshed = THEBAY.sourcemap[baysource]['wshed']            # get column name of watershed flow source for this stream
		flowdf[baysource] = flowdf[wshed] * bs_prop             # scale source from hydromodel flow (some are split)
		bs_name = THEBAY.sourcemap[baysource]['name']
		filename = bs_name + '_Flow.dat'
		logger.info('Bay Source File to Generate: '+filename)
		# Write Inflow File
		# open the file in output directory
		pathedfile = os.path.join(THEBAY.infile_dir, filename)
		with open(pathedfile, mode='w', newline='') as output_file:
			THEBAY.addfile(fname=filename)    # remember generated file names
			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ THEBAY.bayid + '                         !\n')
			output_file.write('! Bay Source: ' + bs_name + '                         !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('1 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write('0    ' + baysource + '\n')
			output_file.write('TIME      INFLOW\n')
			# output the ordinal date and flow value time dataframe columns
			flowdf.to_csv(path_or_buf = output_file, columns= ['ordinaldate', baysource], float_format='%.3f',
			sep=' ', index=False, header=False)

##
#       End of Flow Data Import
#
##

#########################################################################
#
#   All Things Meteorology Related Handled Herein
#
##########################################################################

def adjustCRTemp(air_data):
	# S.E.T. - 20240206 - Adjust Colchester Reef Observed temp by month for the Missisquoi Bay Zone (401)
	# the adjustment, by month, to apply to CR data for use in MissBay
	tempadjust = {
    	1 : -2.6,
    	2 : -2.4,
    	3 : -0.7,
    	4 : 0.8,
    	5 : 1.3,
    	6 : 0.8,
    	7 : -0.4,
    	8 : -0.8,
    	9 : -1.0,
    	10 : -1.2,
    	11 : -1.6,
    	12 : -2.4 }

	adjustedCRtemp = pd.Series()
	for row in range(air_data.shape[0]):
		adjustedCRtemp[row] = air_data.iloc[row] + tempadjust[air_data.index[row].month]
	adjustedCRtemp = adjustedCRtemp.set_axis(air_data.index)
	return adjustedCRtemp


# Class for adjusting shortwave radiation

class ShortwaveNudger:

	nudge_df = None

	@classmethod
	def initialize(cls,nudgefile):
		if type(cls.nudge_df) != pd.core.frame.DataFrame:
			cls.nudge_df = pd.read_csv(nudgefile,index_col="Day of Year")

	@classmethod
	def nudgeDF(cls, swdownobj):

		# the day of year for each data record to nudge
		#dayofyear = dt.strftime(swdownDF.index, '%-j')
		dayofyear = swdownobj.index.strftime('%j').astype(int).to_series()

		#print('Days of Year ',dayofyear)
		#print(cls.nudge_df)

		multiplier = cls.nudge_df.iloc[dayofyear]['Ratio (MB/CR)'] # the multiplier for each data record
		#print('Multipliers ',multiplier)

		#swdownobj_nudged = swdownobj.reset_index(drop=True) * multiplier.reset_index(drop=True)
		swdownobj_nudged = swdownobj.multiply(multiplier.to_numpy()) 
 
		#swdownready = pd.Series(swdownobj_nudged.array,swdownobj.index)
		return swdownobj_nudged


def adjustFEMCLCD(whichbay, dataset):
	# BTV rain adjustment
	dataset['403']['RAIN'] = remove_nas(dataset['403']['RAIN']) * 0.6096
	# define a function to set relative humidity values to 100 if greater than 100
	# seems to be a problem in observations prior to 6/6/2019 in colchesterReefFEMC/Z0080_CR_QAQC.csv
	cap_rhum_at_100 = lambda x: 100 if x > 100 else x
	for zone in get_climate_zone_keys(dataset):
		# air temp and swr adjustments
		if zone == '401':
			dataset[zone]['T2'] = remove_nas(adjustCRTemp(dataset[zone]['T2']))

			ShortwaveNudger.initialize(os.path.join(whichbay.template_dir, 'SolarRadiationFactor_MB.csv'))
			dataset[zone]['SWDOWN'] = ShortwaveNudger.nudgeDF(dataset[zone]['SWDOWN'])
		else :
			dataset[zone]['T2'] = remove_nas(dataset[zone]['T2'])

		# wind speed adjustments
		if zone == '403':
			dataset[zone]['WSPEED'] = remove_nas(dataset[zone]['WSPEED']) * 0.75
		else :
			dataset[zone]['WSPEED'] = remove_nas(dataset[zone]['WSPEED']) * 0.65
		
		# Removing NAs for wind direction, relative humidity, and short-wave radiation
		dataset[zone]['WDIR'] = remove_nas(dataset[zone]['WDIR'])
		dataset[zone]['RH2'] = remove_nas(dataset[zone]['RH2'].apply(cap_rhum_at_100))
		# logger.info("REL_HUM ABOIVE 100:")
		# logger.info(dataset[zone]['RH2'][dataset[zone]['RH2'] > 100])
		# logger.info("The abover should be an empty series")
		dataset[zone]['SWDOWN'] = remove_nas(dataset[zone]['SWDOWN'])
	return dataset

def lakeht_est(whichbay, dataset):
	# Estimate predicted Lake Level based on rolling averages of Missisquoi flow and Airtemp

	THEBAY = whichbay
	########################################### Formerly...
	# # TODO: Which to use?  streamflow_cms or streamflow_adj_cms
	# streamflow_unadj_cms = flowdf['msflow'] / 0.98
	# flowdf['flowmean_07'] = streamflow_unadj_cms.rolling(window=7,min_periods=1).mean() # moving average over 7 days
	# flowdf['flowmean_30'] = streamflow_unadj_cms.rolling(window=30,min_periods=1).mean() # moving average over 30 days
	# flowdf['flowmean_60'] = streamflow_unadj_cms.rolling(window=60,min_periods=1).mean() # moving average over 60 days

	######################################### Now...
	lakeLevel_df = pd.DataFrame({'ordinaldate': THEBAY.flowdf['ordinaldate'],
								 'msflow': THEBAY.flowdf['msflow'],
								 'flowmean_07': THEBAY.flowdf['msflow'].rolling(window=7,min_periods=1).mean(),
								 'flowmean_30': THEBAY.flowdf['msflow'].rolling(window=30,min_periods=1).mean(),
								 'flowmean_60': THEBAY.flowdf['msflow'].rolling(window=60,min_periods=1).mean()},
								 index=pd.DatetimeIndex(THEBAY.flowdf.index, name='time')
								 )
	logger.info(f'lakelevel_df pre merge')
	logger.info(print_df(lakeLevel_df))
	
	###################
	# Need to get parity with timestamps between flow and temp
	## Formerly...
	# # drop any null leap days generated by resample and convert to F
	# wrfdailyraw = air_temp['403'].resample('D').mean()
	# wrfdailyF = 32 + 1.8 * wrfdailyraw.dropna(axis=0, inplace=False, how='any')
	# TemperatureF = wrfdailyF.to_numpy()

	## Now...
	### Need to merge these df's together on time stamp
	### Moved from inner merge to outer merge 3/7/2024 to account for mis-matched time stamps between time series and
	###   data gaps in either time series
	# lakeLevel_df = pd.merge(lakeLevel_df, air_temp['403'], on='time', how='inner')
	#lakeLevel_df = pd.merge(lakeLevel_df, air_temp['403'], on='time', how='outer').interpolate(method="time").dropna()
	logger.info(f'LakeLevel Dataset')
	logger.info(f'{dataset}')
	lakeLevel_df = pd.merge(lakeLevel_df, dataset['403']['T2'].to_frame(), on='time', how='outer').interpolate(method="time").dropna()

	logger.info(f'lakelevel_df after merge')
	logger.info(print_df(lakeLevel_df))

	print('Calculating Lake Levels')
	lakeLevel_df['LakeLevel'] = 94.05887 + 0.007910834 * lakeLevel_df['T2'] + \
		7.034478e-05 * lakeLevel_df['msflow'] + \
		0.003396492 * lakeLevel_df['flowmean_07'] + \
		0.01173037 * lakeLevel_df['flowmean_30'] + \
		0.0258206 * lakeLevel_df['flowmean_60']

	print('Calculating Lake Level Bias')
	# Next, apply the bias correction from the bias correction quadratic regression on the residuals against the observed lake level:
	bias_correction = -358.51020205 + \
		7.16150850 * lakeLevel_df['LakeLevel'] + \
		-0.03570562 * np.power(lakeLevel_df['LakeLevel'],2)

	lakeLevel_df['LakeLevel_corrected'] = lakeLevel_df['LakeLevel'] + bias_correction


	# AEM3D lake input defined as meters above 93ft - do the math
	lakeLevel_df['LakeLevel_delta'] = (lakeLevel_df['LakeLevel_corrected'] - 93) * 0.3048

	return lakeLevel_df['LakeLevel_delta']

def adjustGFS(whichbay, dataset):
	# slice GFS data at forecast start date, so that AEM3D uses observed dataset up to forecast start date
	for zone in get_climate_zone_keys(dataset):
		# GFS rain adjustment
		dataset[zone]['RAIN'] = dataset[zone]['RAIN'] * 86.4
		# GFS temperature adjustment
		dataset[zone]['T2'] = dataset[zone]['T2']-273.15
		# GFS TCDC adjustment
		# Divide GFS TCDC by 100 to get true percentage (see new adjustments function)
		dataset[zone]['TCDC'] = dataset[zone]['TCDC']/100.0
		# GFS windspeed adjustments
		dataset[zone]['WSPEED'] = np.sqrt(np.square(dataset[zone]['U10']) + np.square(dataset[zone]['V10']))
		# GFS wind direction adjustments
		dataset[zone]['WDIR'] = 180 + np.arctan2(dataset[zone]['U10'], dataset[zone]['V10']) * 180 / np.pi

	#dataset['300']['LAKEHT'] = lakeht_est(whichbay, dataset)
	dataset['300'] = {'LAKEHT' : lakeht_est(whichbay, dataset)}

	return dataset

def genclimatefiles(whichbay, settings):

	global SCENARIO

	THEBAY = whichbay   # passed object defining bay characteristics
	year=THEBAY.year   # starting year pulled from IAMBAY class object (lib.py)

	#
	#   Read WRF climate data
	#
	logger.info('Processing Meterological Data')

	##### adjusting spinup date for LCD dataset due to weird blackout dates between 12/26/2021 - 12/31/2021. Data is there if you request an earlier date, just directly request these dates
	'''
	# blackout begins on this date
	blackout_start = dt.datetime(2021,12,26, tzinfo = dt.timezone.utc)
	# blackout is over by this date
	blackout_end = dt.datetime(2022,1,1, tzinfo = dt.timezone.utc)
	blackout_dates = [blackout_start + dt.timedelta(days=d) for d in range((blackout_end - blackout_start).days)]
	
	adjusted_spinup = settings['spinup_date'] - dt.timedelta(days=1)
	if adjusted_spinup in blackout_dates:
		adjusted_spinup = blackout_start - dt.timedelta(days=1)
	'''
	# NEW METHOD because blackout seems to be for 2020 as well
	# Move a day back so AEM3D is happy
	adjusted_spinup = settings['spinup_date'] - dt.timedelta(days=1)
	# Keep moving a day back until out of blackout region
	while (adjusted_spinup.year in [2018, 2019, 2020, 2021] and adjusted_spinup.month == 12 and adjusted_spinup.day > 25):
		adjusted_spinup = adjusted_spinup - dt.timedelta(days=1)
	print(f"ADJUSTED SPINUP: {adjusted_spinup}")
	#####

	observedClimate = {}
	##### GRABBING OBSERVATIONAL CLIMATE DATA #####
	if settings['weather_dataset_observed'] == 'NOAA_LCD+FEMC_CR':
		observedClimateBTV = lcd_ob.get_data(start_date = adjusted_spinup,
										end_date = settings['forecast_start'],
										locations = {"401":"72617014742"},
										variables = {'TCDC':'HourlySkyConditions',
					   								 'RAIN':'HourlyPrecipitation'})
		
		observedClimateCR = femc_ob.get_data(start_date = adjusted_spinup,
										end_date = settings['forecast_start'],
										locations = {'401':'ColReefQAQC'},
										data_dir = settings['data_dir'])
		
		###### FEMC+LCD DATA ADJUSTMENTS HERE ######
		# make additional location dictionaries here rather than call for the same location multiple times in get_data()'s
		for key in ['402', '403']:
			observedClimateBTV[key] = copy.deepcopy(observedClimateBTV['401'])
			observedClimateCR[key] = copy.deepcopy(observedClimateCR['401'])

		# For MB (zone 401), we actually want cloud cover from Franklin Airport
		# Franklin airport ID: 00152
		# common code prefix for vermont stations: 726170
		observedClimateFSO = lcd_ob.get_data(start_date = adjusted_spinup,
									  	  end_date = settings['forecast_start'],
									  	  locations = {"401":"72049400152"},
										  variables = {'TCDC':'HourlySkyConditions'})

		logger.info("Observed TCDC BTV:")
		logger.info(observedClimateBTV['401']['TCDC'].info())
		# now overwrite cloud cover for 401 - combine franklin cloud cover with that from BTV to fill in data gaps
		observedClimateBTV['401']['TCDC'] = combine_timeseries(primary = observedClimateFSO['401']['TCDC'],
															   secondary = observedClimateBTV['401']['TCDC'],
															   interval = '20min')
		logger.info("Observed TCDC FSO:")
		logger.info(observedClimateBTV['401']['TCDC'].info())
		# cool new way to combine dictionaries (python >= 3.9)for zone, ds in climateObsCR.items():
		# Both FEMC and LCD data grabbers now need mathing location keys to work since they are being combine
		for zone in observedClimateCR.keys():
			observedClimate[zone] = observedClimateCR[zone] | observedClimateBTV[zone]
		
		# fetch observed Richelieu height for use as Lake Height
		# Richelieu data represents observed lake height
		getvars = {'LAKEHT':'62614'}
		observedlake = usgs_ob.get_data(start_date = settings['spinup_date'] - dt.timedelta(days=1),
								 		 end_date = settings['forecast_start'],
										 locations = {"RL":'04295000'},
                                          variables = getvars)
		# adjust height reference to 93 ft and convert to meters
		observedlake['RL']['LAKEHT'] = (observedlake['RL']['LAKEHT']-93) * 0.3048
		# store observed lake height in bay object for later concat with predicted height
		observedClimate['300'] = observedlake['RL']


		observedClimate = adjustFEMCLCD(THEBAY, observedClimate)

	else:
		raise ValueError(f"'{settings['weather_dataset_observed']}' is not a valid observational weather dataset")

	forecastClimate = {}
	##### GRABBING FORECASTED CLIMATE DATA #####
	if settings['weather_dataset_forecast'] == 'NOAA_LCD+FEMC_CR':
		# Adjust end_date for FEMC data gap 2021-06-14 21:15:00+00:00 -> 2021-06-18 10:45:00+00:00
		adjusted_end_date = settings['forecast_end'] + dt.timedelta(days=1)
		while adjusted_end_date.year == 2021 and adjusted_end_date.month == 6 and adjusted_end_date.day in [16, 17, 18]:
			adjusted_end_date = adjusted_end_date + dt.timedelta(days=1)
		print(f"ADJUSTED END_DATE: {adjusted_end_date}")

		forecastClimateBTV = lcd_ob.get_data(start_date = settings['forecast_start'],
								end_date = adjusted_end_date,
								locations = {"401":"72617014742"},
								variables = {'TCDC':'HourlySkyConditions',
					   						 'RAIN':'HourlyPrecipitation'})
		
		forecastClimateCR = femc_ob.get_data(start_date = settings['forecast_start'],
										end_date = adjusted_end_date,
										locations = {'401':'ColReefQAQC'},
										data_dir = settings['data_dir'])
		
		
		# copy data from 401 for 402, 403. Make sure it's a deep copy, not memeory reference
		###### FEMC+LCD DATA ADJUSTMENTS HERE ######
		# make additional location dictionaries here rather than call for the same location multiple times in get_data()'s
		for key in ['402', '403']:
			forecastClimateBTV[key] = copy.deepcopy(forecastClimateBTV['401'])
			forecastClimateCR[key] = copy.deepcopy(forecastClimateCR['401'])

		# For MB (zone 401), we actually want cloud cover from Franklin Airport
		# Franklin airport ID: 00152
		# common code prefix for vermont stations: 726170
		# try Franklin data grab - 7-day data grab might return empty json, so in that case, use BTV data
		
		# boolean switch - if true, we are using franklin aiport cloud data for the forecast period
		fso_cloud_forecast = True
		try:
			# trying Franklin grab, and if that doesn't work...
			logger.info("Trying to get Forecast Period cloud data from Franklin airport (72049400152)...")
			forecastClimateFSO = lcd_ob.get_data(start_date = settings['forecast_start'],
									  	  end_date = adjusted_end_date,
									  	  locations = {"401":"72049400152"},
										  variables = {'TCDC':'HourlySkyConditions'})
		except Exception as e:
			# use Burlington data
			logger.warning("Franklin airport cloud data grab for Forecast Period Failed.")
			logger.warning(e)
			logger.warning("Getting Forecast Period cloud data from Burlington instead...")
			fso_cloud_forecast = False
			
		logger.info("forecast TCDC BTV:")
		logger.info(forecastClimateBTV['401']['TCDC'].info())

		# only need to combine franklin and burlington cloud data if franklin data grab worked
		# otherwise, no need to combine, will just use BTV cloud data
		if fso_cloud_forecast:
			# now overwrite cloud cover for 401 - combine franklin cloud cover with that from BTV to fill in data gaps
			forecastClimateBTV['401']['TCDC'] = combine_timeseries(primary = forecastClimateFSO['401']['TCDC'],
																secondary = forecastClimateBTV['401']['TCDC'],
																interval = '20min')
			logger.info("forecast TCDC FSO:")
			logger.info(forecastClimateBTV['401']['TCDC'].info())

		# cool new way to combine dictionaries (python >= 3.9)for zone, ds in climateObsCR.items():
		# Both FEMC and LCD data grabbers now need mathing location keys to work since they are being combine
		for zone in forecastClimateCR.keys():
			forecastClimate[zone] = forecastClimateCR[zone] | forecastClimateBTV[zone]

		# Richelieu data represents forecast lake height
		getvars = {'LAKEHT':'62614'}
		forecastlake = usgs_ob.get_data(start_date = settings['forecast_start'],
								 		end_date = adjusted_end_date,
										locations = {"RL":'04295000'},
                                        variables = getvars)
		# adjust height reference to 93 ft and convert to meters
		forecastlake['RL']['LAKEHT'] = (forecastlake['RL']['LAKEHT']-93) * 0.034478e-05
		# store observed lake height in bay object for later concat with predicted height
		forecastClimate['300'] = forecastlake['RL']

		forecastClimate = adjustFEMCLCD(THEBAY, forecastClimate)

	elif settings['weather_dataset_forecast'] == 'NOAA_GFS':
		############## Use this bit to load forecast climate from .csvs previously created above
		if settings['csv']:
			logger.info("Loading GFS From CSVs")
			forecastClimate = {}
			for zone in ['401', '402', '403']:
				forecastClimate[zone] = pd.read_csv(
							os.path.join(settings['data_dir'], f'gfs{zone}.csv'),
							index_col='time',
							parse_dates=True)
		else:
		############## Use this bit to load forecast climate from original GRIB files and create .csvs for quick loading later
			logger.info(f"Begin GFS get_data()")
			# determine the timedelta adjustment for GFS based of NWM forecast memeber
			member = settings['nwm_forecast_member'][-1]
			# if there is no member num (such as 'short_range/') default to 1 so that there is no timedelta adjust
			if not member.isdigit():
				member = '1'
			# NOTE forecast_end adjustment: between forecast_start and end, we only want 7.5 days of data.
				# forecast_start will be adjusted based on the NWM member, but end will stay the same
				# this means that we will be grabbing more data ofr higher members
				# ex. member 7 uses a GFS adjustment of forecast_start-1.5 days, so will grab 9 days total (1.5+7.5)
			forecastClimate = gfs_fc_thredds.get_data(forecast_datetime = settings['forecast_start'] - dt.timedelta(hours=(6*(int(member)-1))),
													end_datetime = settings['forecast_end'] + dt.timedelta(hours=12),
													locations = {'401': (45.00, -73.25),
																'402': (44.75, -73.25),
																'403': (44.75, -73.25)},
													data_dir = settings['data_dir'],
													load_threads = 1)
			
			###### GFS DATA ADJUSTMENTS HERE ######
			# slice GFS data at forecast start date, so that AEM3D uses observed dataset up to forecast start date
			for zone in forecastClimate.keys():
				for var in forecastClimate[zone].keys():
					forecastClimate[zone][var] = forecastClimate[zone][var].loc[pd.to_datetime(settings['forecast_start'].replace(tzinfo=dt.timezone.utc)):]
			forecastClimate = adjustGFS(THEBAY, forecastClimate)

	logger.info('Climate Observed Data')
	logger.info(observedClimate)
	logger.info('Climate Forecast Data (Zone 401)')
	logger.info(forecastClimate)

	air_temp = {"401": pd.concat([observedClimate['401']['T2'],forecastClimate['401']['T2']]),
				"402": pd.concat([observedClimate['402']['T2'],forecastClimate['402']['T2']]),
				"403": pd.concat([observedClimate['403']['T2'],forecastClimate['403']['T2']])}

	global SUBPLOT_PACKAGES
	global AXES
	SUBPLOT_PACKAGES['air temp'] = {'labelled_data':air_temp,
								   	'ylabel':'Air Temperature at 2m (C)',
									'title':f"{settings['weather_dataset_observed']} vs. {settings['weather_dataset_forecast']}",
									'fc_start':settings['forecast_start'],
									'row':0,
									'col':1,
									'axis':AXES}
	
	logger.info(print_df(air_temp['401']))

	'''
	# OLD METHOD FOR WATER TEMP - BASED ON ZONE 403 FOR ALL BAYSOURCES
	# Use air temp at zone 403 (ILS)
	wtr_temp = air_temp['403'].rolling(window=96,min_periods=1).mean() # moving average over 4 days

	wtr_temp.loc[wtr_temp<0] = 0  # no subfreezing water
	wtr_temp = wtr_temp + 0.75    # 0.75 correction based on WQS Docs 2021.05.27
	'''

	wtr_temp_dict = {
		'201': '401',
		'202': '401',
		'203': '401',
		'204': '401',
		'21':  '401',
		'22':  '401',
		'17':  '402',
		'19':  '402'
	}

	wtr_temp_zones = {'401' : air_temp['401'].resample("15min").interpolate("time").rolling(window=96,min_periods=1).mean(),
				   	  '402' : air_temp['402'].resample("15min").interpolate("time").rolling(window=96,min_periods=1).mean()}
	
	# General water temp nudges
	for zone, data in wtr_temp_zones.items():
		wtr_temp_zones[zone].loc[wtr_temp_zones[zone] < 0.0] = 0.0  # no subfreezing water
		wtr_temp_zones[zone] = wtr_temp_zones[zone] + 0.75          # 0.75 correction based on WQS Docs 2021.05.27

	# Create water temp dictionary for bay object for later use in wq calcs
	# THEBAY.tempdf = wrfdf[['ordinaldate', 'wtr_temp']].copy()
	THEBAY.wtr_temp_dict = {}
	for inflow, climZone in wtr_temp_dict.items():
		THEBAY.wtr_temp_dict[inflow] = wtr_temp_zones[climZone]

	#
	#       Write the temperature file for each source of the bay
	#
	for baysource in THEBAY.sourcelist :

		bs_name = THEBAY.sourcemap[baysource]['name']
		filename = bs_name + '_Temp.dat'
		logger.info('Generating Bay Source Temperature File: '+filename)

		# Write Temp File
		# open the file in output directory
		pathedfile = os.path.join(THEBAY.infile_dir, filename)
		with open(pathedfile, mode='w', newline='') as output_file:

			THEBAY.addfile(fname=filename)        # remember generated bay files

			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ THEBAY.bayid + '                            !\n')
			output_file.write('! Bay Source: ' + bs_name + '                         !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('1 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write('0    '+baysource+'\n')
			output_file.write('TIME      WTR_TEMP\n')

			# output the ordinal date and temp dataframe columns
			index_to_ordinal_date(wtr_temp_zones[wtr_temp_dict[baysource]]).to_csv(path_or_buf = output_file, float_format='%.3f',
			sep=' ', index=True, header=False)

	#
	#   end water temp file

	######################################################################################
	#
	#   Rain (WRF RAIN) and Snow
	#
	bay_rain = dict()
	bay_snow = dict()

	### TODO: Remove hardcoded single RAIN zone
	for zone in ['403']:
		
		# ['RAIN'] for BTV to get to a series
		# AEM3D wants meters/day
		# BTV is in inches/hr so thats * 24 to get hr -> day and * 0.0254 to convert inches to meters
		#   for a total of 0.6096
		#   LCD Ref: https://www.ncei.noaa.gov/pub/data/cdo/documentation/LCD_documentation.pdf
		# GFS is in kg/m^2/s. kg/m^2 is mm, so really mm/s. To convert, * 86400 to get sec -> day,
		#   and / 1000 to get mm to m, so, in all, * 86.4
		#   GFS Ref: https://www.nco.ncep.noaa.gov/pmb/products/gfs/gfs.t00z.pgrb2.0p25.f003.shtml
		bay_rain[zone] = pd.concat([observedClimate[zone]['RAIN'], forecastClimate[zone]['RAIN']])
		# give series a name... beacuse pandas wants one for the merge below
		bay_rain[zone] = bay_rain[zone].rename('RAIN')

		print(air_temp[zone])
		TEMP = air_temp[zone].reindex(bay_rain[zone].index, method='nearest')

		logger.info('TEMP after reindex')
		logger.info(print_df(TEMP))

		# Original Try: Use https://www.ncdc.noaa.gov/sites/default/files/attachments/Estimating_the_Water_Equivalent_of_Snow.pdf
		#   and fit a quadratic regression through it
		#   SNOW = e^(2.3129520 + -0.0962303*TEMP(C) + -0.0009388*TEMP(C)*TEMP(C)) * PRECIP
		#   This was too much though, so looking at the calibration input, divide in half
		#   to get closer to calibration input
		#snowcoeff = np.exp(2.3129520 - 0.0962303*TEMP - 0.0009388*TEMP*TEMP)/2.0
		
		# Take Two: Using GHCN data from BTV (Burlington Airport)
		# See Pat's rainToSnow.R
		# Calculate ratio using SNOW/PRCP, only take 1991 - present and days where SNOW and PRCP > 0
		# Tried using TMAX, but TAVG ((TMAX + TMIN)/2) got a curve closer to NASA table and makes more sense
		snowcoeff = np.exp(2.1413626 - 0.1921400*TEMP - 0.0079924*TEMP*TEMP)

		logger.info('bay_rain')
		logger.info(print_df(bay_rain[zone]))
		
		# Also a series, name = T2
		logger.info('snowcoeff')
		logger.info(print_df(snowcoeff))

		#####################

		# Merge bay_rain and snowcoeff by time stamps so we know we have times for both
		snowcalc_df = pd.merge(bay_rain[zone], snowcoeff, on='time', how='inner')
		logger.info('snowcalc_df')
		logger.info(print_df(snowcalc_df))
		bay_snow[zone] = snowcalc_df['RAIN'] * snowcalc_df['T2']
		
		# remove duplicate indices from bay_snow, bay_rain and TEMP
		# duplicated timestamps was creating some concat and .loc issues down the line
		# removing these duplicated timestamps seems to work nicely
		bay_rain[zone] = bay_rain[zone][~bay_rain[zone].index.duplicated()]
		bay_snow[zone] = bay_snow[zone][~bay_snow[zone].index.duplicated()]
		TEMP = TEMP[~TEMP.index.duplicated()]

		# If too warm, no snow - Use -2.0 C to get mean snowfall for 2017-2020 close to calibration data
		#  NOAA table above uses 34 F (1.1 C)
		bay_snow[zone].loc[TEMP > -2.0] = 0.0
		# In the WQS AEM3D calibration SNOW / RAIN data, RAIN looks like snow equivalent
		#   and SNOW is depth of snow, so we should have values for both when it snows
		# bay_rain[zone].loc[bay_snow[zone] > 0.0] = 0.0

		logger.info('bay_snow before write')
		logger.info(print_df(bay_snow[zone]))        
		logger.info('bay_rain before write')
		logger.info(print_df(bay_rain[zone])) 

		logger.info(f'Snow Stats for zone {zone}')
		logger.info(bay_snow[zone].describe(percentiles=[]))
		logger.info(f'Rain Stats for zone {zone}')
		logger.info(bay_rain[zone].describe(percentiles=[]))		

	SUBPLOT_PACKAGES['bay rain'] = {'labelled_data':bay_rain,
									'ylabel':'Rainfall (m/day)',
									'title':'Rainfall for Inland Sea',
									'fc_start':settings['forecast_start'],
									'row':0,
									'col':2,
									'axis':AXES}
	
	SUBPLOT_PACKAGES['bay snow'] = {'labelled_data':bay_snow,
									'ylabel':'Snow Fall (m/day)',
									'title':'Snow Fall for Inland Sea',
									'fc_start':settings['forecast_start'],
									'row':0,
									'col':3,
									'axis':AXES}

	#
	# Write Precip Files for Bay
	#
	for zone in bay_rain.keys():
		#################################
		# Hardcode zone to 0 for now... only one precip station
		filename = f'PRECIP_0.dat'
		logger.info('Generating Bay Precipitation File: '+filename)

		with open(os.path.join(THEBAY.infile_dir, filename), mode='w', newline='') as output_file:

			THEBAY.addfile(fname=filename)        # remember generated bay files

			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ THEBAY.bayid + '                         !\n')
			#output_file.write('! Bay Source: ' + bs_name + '                        !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('2 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write('0    0    0\n')
			output_file.write('TIME      RAIN     SNOW\n')

			# output the ordinal date and flow value time dataframe columns
			pd.concat([
					index_to_ordinal_date(bay_rain[zone]),
					index_to_ordinal_date(bay_snow[zone])],
					axis=1).to_csv(
					path_or_buf = output_file,
					float_format='%.3f',
					sep=' ',
					index=True, header=False)


	###########################################################################################
	#
	#   Longwave Radiation Measure : CLOUDS, LW_RAD_IN, or LW_RAD_NET
	#

	# Use Cloud Cover
	cloud_plot_data = {}
	for zone in get_climate_zone_keys(forecastClimate):

		filename = f'CLOUDS_{zone}.dat'
		logger.info('Generating Bay Cloud Cover File: '+filename)

		full_cloud_series = pd.concat([observedClimate[zone]['TCDC'],forecastClimate[zone]['TCDC']])
		cloud_plot_data[zone] = full_cloud_series

		cloud_series = index_to_ordinal_date(full_cloud_series)
		logger.info(f'TCDC for zone {zone}') 
		logger.info(print_df(cloud_series))

		writeFile(
			os.path.join(THEBAY.infile_dir, filename),
			THEBAY.bayid,
			zone,
			"CLOUDS",
			cloud_series
		)
		THEBAY.addfile(fname=filename)

	SUBPLOT_PACKAGES['cloud cover'] = {'labelled_data':cloud_plot_data,
										'ylabel':'Total Cloud Cover (%)',
										'title':'Cloud Cover',
										'fc_start':settings['forecast_start'],
										'row':1,
										'col':0,
										'axis':AXES}

	###########################################################################################
	#
	#   Wind Speed and Direction
	#           U10 and V10 variables converted to Earth compass direction and speed
	#

	windspd = dict()
	winddir = dict()
	for zone in get_climate_zone_keys(forecastClimate):

		# # a bit of vector math to combine the East(U) and North(V) wind components
		# windspd[zone] = np.sqrt(
		# 	np.square(forecastClimate[zone]['U10']) +
		# 	np.square(forecastClimate[zone]['V10'])
		# )
		windspd[zone] = pd.concat([observedClimate[zone]['WSPEED'], forecastClimate[zone]['WSPEED']])

		# #  a bit of trig to map the wind vector components into a direction
		# #  𝜙 =180+(180/𝜋)*atan2(𝑢,𝑣)
		# winddir[zone] = 180 + np.arctan2(
		# 		forecastClimate[zone]['U10'],
		# 		forecastClimate[zone]['V10']
		# 	) * 180 / np.pi
		winddir[zone] = pd.concat([observedClimate[zone]['WDIR'], forecastClimate[zone]['WDIR']])

		logger.info(f'WINDSP for zone {zone}')
		logger.info(print_df(windspd[zone]))        
		logger.info(f'WINDDIR for zone {zone}')
		logger.info(print_df(winddir[zone]))        

	SUBPLOT_PACKAGES['wind speed'] = {'labelled_data':windspd,
										'ylabel':'Wind Speed at 10m (m/s)',
										'title':'Wind Speed',
										'fc_start':settings['forecast_start'],
										'row':1,
										'col':1,
										'axis':AXES}
	
	SUBPLOT_PACKAGES['wind direction'] = {'labelled_data':winddir,
										'ylabel':'Wind Direction at 10m (degrees clockwise from North)',
										'title':'Wind Direction',
										'fc_start':settings['forecast_start'],
										'row':1,
										'col':2,
										'axis':AXES}

	# Write Wind Speed and Direction File
	#
	for zone in windspd.keys():
		filename = f'WS_WD_{zone}.dat'
		logger.info('Generating Wind Speed and Direction File: '+filename)
		# logger.info(index_to_ordinal_date(windspd[zone]))
		# logger.info(index_to_ordinal_date(winddir[zone]))
		with open(os.path.join(THEBAY.infile_dir, filename), mode='w', newline='') as output_file:

			THEBAY.addfile(fname=filename)        # remember generated bay files

			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ THEBAY.bayid + '                            !\n')
			#output_file.write('! Bay Source: ' + bs_name + '                         !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('2 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write(f'0            {zone}          {zone}\n')
			output_file.write('TIME         WIND_SPEED	  WIND_DIR\n')

			# output the ordinal date and flow value time dataframe columns
			'''
			pd.concat([
				index_to_ordinal_date(windspd[zone]),
				index_to_ordinal_date(winddir[zone])],
				axis=1).drop_na().sort().to_csv(
				path_or_buf = output_file,
				float_format='%.3f',
				sep=' ',
				index=True, header=False)
			'''
			# now interpolating instead of dropping rows with missing data afger calibration
			# converting index to ordinal date last, as we need a datetime index in order to interpolate
			index_to_ordinal_date(pd.concat([
				windspd[zone],
				winddir[zone]],
				axis=1).sort_index().interpolate(method='time')).to_csv(
				path_or_buf = output_file,
				float_format='%.3f',
				sep=' ',
				index=True, header=False)
	#
	#   end wind file



	###########################################################################################
	#
	#   Relative Humidity : WRF, Calculated from Q2
	#

	rhum = {}
	rhum['403'] = pd.concat([observedClimate['403']['RH2'] * .01, forecastClimate['403']['RH2'] * .01])   
	
	logger.info(f'RH2')
	logger.info(print_df(rhum['403']))

	SUBPLOT_PACKAGES['relative humidity'] = {'labelled_data':rhum,
										'ylabel':'Relative Humidity at 2m (%)',
										'title':'Relative Humidity',
										'fc_start':settings['forecast_start'],
										'row':1,
										'col':3,
										'axis':AXES}

	#   Write Relative Humidity File
	for zone in rhum.keys():
		filename = f'RELHUM_{zone}.dat'
		logger.info('Generating Rel Humidity File: '+filename)
		writeFile(
			os.path.join(THEBAY.infile_dir, filename),
			THEBAY.bayid,
			zone,
			"REL_HUM",
			index_to_ordinal_date(rhum[zone]))
		THEBAY.addfile(fname=filename)


	#   Write Air Temp File
	for zone in air_temp.keys():
		filename = f'AIRTEMP_{zone}.dat'
		logger.info('Generating Air Temp File: '+filename)
		writeFile(
			os.path.join(THEBAY.infile_dir, filename),
			THEBAY.bayid,
			zone,
			"AIR_TEMP",
			index_to_ordinal_date(air_temp[zone]))
		THEBAY.addfile(fname=filename)


	# Write ShortwaveRad File
	
	# SET 20240326 - Depricated by move of nudge to adjustFEMCLCD
	# Load the CSV that nudges SWR values based on Day Of Year
	#swradjust = pd.read_csv(os.path.join(THEBAY.template_dir, 'SolarRadiationFactor_MB.csv'))
	## Use the passed forecast date to only nudge the observed values before the forecast
	#swradjust['NudgeObs'] = swradjust['Ratio (MB/CR)']  # copy entire original nudge column
	#dayofyear = int(settings['forecast_start'].strftime('%j'))
	#swradjust.loc[swradjust['Day of Year'] >= dayofyear,'NudgeObs'] = 1.0   # don't adjust forecast data

	swdown_plot_data = {}
	for zone in get_climate_zone_keys(forecastClimate):

		filename = f'SOLAR_{zone}.dat'
		logger.info('Generating Short Wave Radiation File: '+filename)
		full_swdown_series = pd.concat([observedClimate[zone]['SWDOWN'], forecastClimate[zone]['SWDOWN']])
		swdown_plot_data[zone] = full_swdown_series
		swdown_series = index_to_ordinal_date(full_swdown_series)
	
		# SET 20240326 - Depricated by move of nudge to adjustFEMCLCD
		# build a dataframe to nudge the data series
		#swdown_df = pd.concat([swdown_series.index.to_series(),swdown_series], axis = 1)
		#swdown_df.columns = ['TIME', 'SWDOWN']
		#swdown_df['SWNUDGED'] = swdown_df.apply(ordinalnudgerow, args = ('SWDOWN', swradjust), axis = 1)

		logger.info(f'SWDOWN for zone {zone}')
		logger.info(print_df(swdown_series))
	
		# SET 20240326 - Depricated by move of nudge to adjustFEMCLCD
		#writeFile(
		#	os.path.join(THEBAY.infile_dir, filename),
		#	THEBAY.bayid,
		#	zone,
		#	"SOLAR_RAD",
		#	swdown_df['SWNUDGED']
		#)

		writeFile(
			os.path.join(THEBAY.infile_dir, filename),
			THEBAY.bayid,
			zone,
			"SOLAR_RAD",
			swdown_series
		)

		THEBAY.addfile(fname=filename)

	SUBPLOT_PACKAGES['short wave rad'] = {'labelled_data':swdown_plot_data,
										'ylabel':'Short-wave Downward Radition (W/m^2)',
										'title':'Short-Wave Down Rad (Un-nudged)',
										'fc_start':settings['forecast_start'],
										'row':1,
										'col':4,
										'axis':AXES}

	###
	#
	#   Generate Lake Level Regression from Temp and Flow Derived Terms
	#       Calculate lake level based on previous regression model results.
	#		Ported from original IAM/tools/EFDC_file_prep/EFDC_file_prep.R
	#       The regression model equation for lake level is:
	#       LakeLevel = 94.05887 +
	#           0.007910834(temperature) +
	#           7.034478e-05 (discharge) +
	#           0.003396492(discharge7days) +
	#           0.01173037(discharge30days) +
	#           0.0258206(discharge60days)
	#
	###
	logger.info('Generating Lake Levels')

	# Calculate misalignment between first predicted and last observered, and adjust entire prediction that amount
	htoffset = observedClimate['300']['LAKEHT'].iloc[-1] - forecastClimate['300']['LAKEHT'].iloc[0]
	logger.info(f"Last observed = {observedClimate['300']['LAKEHT'].index[-1]}, {observedClimate['300']['LAKEHT'].iloc[-1]}")
	logger.info(f"First estimated = {forecastClimate['300']['LAKEHT'].index[0]}, {forecastClimate['300']['LAKEHT'].iloc[0]}")
	logger.info('Lake Height Prediction Offset = ' + str(htoffset))
	forecastClimate['300']['LAKEHT'] = forecastClimate['300']['LAKEHT'] + htoffset 
	lake_height = pd.concat([observedClimate['300']['LAKEHT'],forecastClimate['300']['LAKEHT']])
	lake_ht_series = index_to_ordinal_date(lake_height)

	#
	#   write out lake level file
	#
	logger.info('Lake Level (m) above 93ft')
	logger.info(lake_ht_series.describe(percentiles=[]))

	filename = 'Lake_Level.dat'
	logger.info('Writing Lake Level File '+filename)

	writeFile(
			os.path.join(THEBAY.infile_dir, filename),
			THEBAY.bayid,
			'300',
			"HEIGHT",
			lake_ht_series
	)
	THEBAY.addfile(fname=filename)

	logger.info(f'lake_height:{lake_height}')
	# adding Lake Level plot package
	SUBPLOT_PACKAGES['lake level'] = {'labelled_data':{'300':lake_height},
									'ylabel':'Lake Level (m) above 93ft',
									'title':'Lake Level',
									'fc_start':settings['forecast_start'],
									'row':0,
									'col':4,
									'axis':AXES}
	
	# pathedfile = os.path.join(THEBAY.infile_dir, filename)
	# with open(pathedfile, mode='w', newline='') as output_file:

	# 	THEBAY.addfile(fname=filename)        # remember generated bay files

	# 	# output the header text
	# 	output_file.write(
	# 		'!-----------------------------------------------------!\n')
	# 	output_file.write(
	# 		'! Written by AEM3D_prep_IAM                           !\n')
	# 	output_file.write('! Bay ID: ' + THEBAY.bayid +
	# 					  '                            !\n')
	# 	output_file.write(
	# 		'! values in (m) above 93 ft                           !\n')
	# 	output_file.write(
	# 		'!-----------------------------------------------------!\n')
	# 	output_file.write('1 data sets\n')
	# 	output_file.write('0 seconds between data\n')
	# 	output_file.write('0    300\n')
	# 	output_file.write('TIME	  HEIGHT\n')

	# 	# output the ordinal date and flow value time dataframe columns for AEM3D and as .csv
	# 	# get series output csv write syntax from other series file writer
	# 	lake_height
	# 	#			  sep=' ', index=False, header=False)
	# 	# lakeLevel_df.to_csv(path_or_buf='lakeheight.csv', float_format='%.3f', sep=' ', index=False, header=True)


	##
	#
	#   End of Lake Level From Temp and Flow Generation
	#
	##


def gensalinefile(theBay):
	#
	#   Salinity data time series
	#
	generate_file_from_template('salinity.template.txt',
								'Inflows_Salinity.dat',
								theBay,
								{'source_id_list': '  '.join(theBay.sourcelist),
								 'firstdate': theBay.FirstDate,
								 'lastdate': theBay.LastDate
								})


def genboundaryfile(theBay):
	#
	#   P boundary condition data time series
	#
	generate_file_from_template('bc_p.template.txt',
								'OpenBC_P.dat',
								theBay,
								{'firstdate': theBay.FirstDate,
								 'lastdate': theBay.LastDate
								})
		

def gentracerfiles(theBay):
	#
	# Write Tracer Files - fixed series, use templates and set $year
	#

	for templateFile in [f for f in os.listdir(theBay.template_dir) if 'Tracer' in f]:
		generate_file_from_template(templateFile,
									os.path.splitext(templateFile)[0]+'.dat',
									theBay,
									{'firstdate': theBay.FirstDate,
									 'lastdate': theBay.LastDate
									})
	
	generate_file_from_template('tracer_release.template',
								'tracer_release.dat',
								theBay,
								{'firstdate': theBay.FirstDate,
								 'lastdate': theBay.LastDate
								},
								'update_file')
	#
	# End of Tracer File Generation
	#


def gencntlfile(theBay, settings):

	# Calculate 1 hour in iterations
	hourIter = int(86400 / AEM3D_DEL_T / 24)

	output_iters = settings['output_write_iteration_hours']*hourIter
	simstart = settings['spinup_date']	# initial sim start time, pending restart override
	restartswitch = 0		# begin with assumption there will be no restart file

	# aem3d control file parser wants a value for restart file key, even if restart inactive
	#		alternate method would be to programmtically comment the line out
	restart_read_base = "norestartinput"	# the base (pathless) name for incomming restart file

	##
	# Check and setup for an incomming restart file
	if settings['restart_read_file'] != "" :
		if not os.path.exists(settings['restart_read_file']):
				logger.critical(
				'File {} does not exist. Exiting.'.format(
				settings['restart_read_file']
 			)
		)
		#exit_with_code_based_on_host()

		# copy restart file to runtime infiles directory
		os.system('cp '+ settings['restart_read_file'] + " " + theBay.infile_dir)

		restart_read_base = os.path.basename(settings['restart_read_file'])  # strip path off

		with open(settings['restart_read_file'], 'rb') as restartfile:
			# extract restart file's timestamp for use in setting interation count

			# unpack original julian date from binary array
			#  number of full days since Jan 1, 4713BC with fraction of day since noon as decimal
			#  3rd record of the unf file after two character (Len=256) vars
			#  Julian date is Fortran double

			# {4bytecount}{256 chars}{copy of 4bytecount}
			# {4bytecount}{256 chars}{copy of 4bytecount}
			# {4bytecount}{Fortran double value for date}

			restartfile.seek((256+8)*2+4)   # skip first two records and a count
			lastsimtime = np.fromfile(restartfile, dtype='float64',count=1)[0]
			epoch = pd.to_datetime(0, unit='s').to_julian_date()  #first datetime (1970)
			simdatetime = pd.to_datetime(lastsimtime - epoch, unit='D', utc=True)

			simstart = simdatetime
			print('Restart time as Datetime ', simdatetime)

		output_iters = output_iters + 1  	# when using restart and extra iteration needed to line up first save
		restartswitch = 1 		# turn on restart in control file
	

	# Calculate iterations: Time between forecast end and start of sim
	iterations = int((settings['forecast_end'] - simstart).total_seconds() / AEM3D_DEL_T)
	logger.info(f'Configuring AEM3D to run {iterations} iterations')


	##
	#	Check and setup for possible restart file generation
	#		If no write_start specified, no saves to schedule
	if settings['output_write_start_datetime'] != None :
		# Convert write_start datetime to a sim iteration count
		output_start_iter = 1 + int((settings['output_write_start_datetime'] - simstart).total_seconds() / AEM3D_DEL_T)
	else :
		output_start_iter = iterations + 1 	# set output start after max sim iterations

	
	with open(os.path.join(theBay.template_dir, 'aem3dcntl.template.txt'), 'r') as file:
		template = Template(file.read())

	# control file is written to runtime directory
	pathedfile = os.path.join(theBay.run_dir, 'run_aem3d.dat')
	with open(pathedfile, 'w') as output_file:
		output_file.write(template.substitute(**{
			'start_date': datetimeToOrdinal(simstart),
			'del_t': AEM3D_DEL_T,
			# number of 300s steps in a 364 days (year minus 1 day, because 1st day is a start, not a step)
			# 27936 for 97 days (97*24*12)
			'iter_max': iterations,
			'hour': hourIter,
			'eighthours': (hourIter * 8),
			'daysthirty': (hourIter * 24 * 30),
			'output_start' : output_start_iter,
			'output_iters': output_iters,
			'userestart' : restartswitch,
			'restart_to_read' : restart_read_base
			}))

		# update the control file with all generated input files
		output_file.write(
			'! --------------------------------------------------------------- !\n')
		output_file.write(
			'! Input file names                                                 !\n')
		output_file.write('\''+os.path.split(theBay.infile_dir)[1] +
						'\'                                     infile_dir\n')
		output_file.write(
			' sparsedata_aem3d.unf                              3D_data_file\n')
		output_file.write(
			' usedata_aem3d.unf                             preprocessor_file\n')
		# output_file.write(
		#     ' datablock.xml                                     datablock_file\n')

		for f, t in theBay.bayfiles:        # for each filename, filetype
			output_file.write(' '+f+'              '+t+'\n')

		output_file.write(
					'!  End                                                            !\n')
# end of control file generation


def gendatablockfile(theBay, settings):

	# Calculate 1 hour in iterations
	hourIter = int(86400 / AEM3D_DEL_T / 24)
	# Calculate iteration for forecast start: Time between forecast date and spinup start
	forecastStartIter = int((settings['forecast_start'] - (settings['spinup_date'] + dt.timedelta(days=1))).total_seconds() / AEM3D_DEL_T) + 1

	logger.info(f'Configuring datablock.xml file')
	generate_file_from_template('datablock.xml.template',
								'datablock.xml',
								theBay,
								{'hour': hourIter,
								 'sixhours': (hourIter * 6),
								 'day': (hourIter * 24),
								 'spinupEnd': (forecastStartIter - 1),
								 'forecastStart': forecastStartIter
								},
								'datablock_file')
# end of datablock.xml generation
 
def AEM3D_prep_IAM(theBay, settings):

	logger.info(f'Processing Bay: {theBay.bayid} for year {theBay.year}')

	# get flow files from hydrology model data, including observed lake level
	getflowfiles(theBay, settings)

	# generate climate files including lake levels
	genclimatefiles(theBay, settings)

	# generate salinity file
	gensalinefile(theBay)

	# generate boundary condition file
	genboundaryfile(theBay)
	
	# generate tracer files
	gentracerfiles(theBay)

	# generate the water quality files (waterquality.py script)
	genwqfiles(theBay)

	# generate the datablock.xml file
	gendatablockfile(theBay, settings)

	# generate control file
	gencntlfile(theBay, settings)

	global FIG
	global SUBPLOT_PACKAGES
	make_figure(SUBPLOT_PACKAGES)
	FIG.savefig('weatherVars.png')

	return 0
