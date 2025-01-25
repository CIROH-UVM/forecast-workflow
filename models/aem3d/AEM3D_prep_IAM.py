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
				  nwmv3_retro_fc,
				  usgs_ob,
				  gfs_fc_thredds,
				  cfs_fc,
				  lcd_ob,
				  caflow_ob,
				  utils
)

from .waterquality import *
from .AEM3D import *
from .Aem3dForcings import *

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import datetime as dt
import pytz
import warnings

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

def add_plot(labelled_data, ylabel, title, fc_start, fc_end, row, col, axis):
	ax = axis[row,col]
	# Plot data for location 1
	for lab, data in labelled_data.items():
		data.index = pd.to_datetime(data.index, utc=True)
		# print("data.index BEFORE")
		# print(data.index)
		
		# if there is no time component to the index...
		if all(d.time() == dt.time(0) for d in data.index) and not data.index.empty:
			# add a time component
			data.index = data.index.map(lambda x: x.replace(hour=12)).tz_localize('UTC')
		# print("data.index AFTER")
		# print(data.index)
		time_sliced_data = data[data.index <= (fc_end + dt.timedelta(days=1))]
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

def scaleRockandPikeQ(missisquoi_q):
	#	approximate Rock and Pike forecast from scaled Missisquoi forecast
	# TODO: Devise a better plan to reference the scaling factors than hardcoding them here
	forecastca = {'RK':{},'PK':{}}
	forecastca['RK']['streamflow'] = missisquoi_q * 0.038  # scale Missisquoi to Rock
	forecastca['PK']['streamflow'] = missisquoi_q * 0.228  # scale Missisquoi to Pike
	return forecastca

def backfillCaFlowsSpinup(ca_data, settings):
	'''
	Fill in missing dates in the Canadian instantaneous data spinup series with daily values.
	Only fills in dates that are completely missing from the instantaneous timerseries with daily means from Canadian daily values (dv) service.
	Note that the nested ditionary ca_data is modified inplace.

	Args:
	-- ca_data (nested dict) [req]: instantaneous canadian streamflow timeseries in standard get_data() format
	-- settings (dict) [req]: dictionary of run settings

	Returns:
	None; modifies ca_data inplace 
	'''
	ca_gauges = {"PK":'030424',"RK":'030425'}
	for loc, iv_dict in ca_data.items():
		print(f"Backfilling {loc} IV streamflow with daily values...")
		# getting daily data, for spinup period, for specific gauge
		spinup = settings['spinup_date'] - dt.timedelta(days=1)
		forecast_start = settings['forecast_start']
		dv_data = caflow_ob.get_data(start_date = spinup,
									end_date = forecast_start,
									locations = {loc:ca_gauges[loc]},
									service = 'dv')
		# note that daily values dates are reported in ET time
		dv = dv_data[loc]['streamflow']
		# note that the index of iv will be in UTC time
		iv = iv_dict['streamflow']
		# convert iv timestamps from UTC to ET before getting dates
		# then get rid of time info (normalize), then drop subsequent date duplicates
		iv_daily_et_index = iv.index.map(lambda t: t.tz_convert(pytz.timezone('America/New_York'))).normalize().drop_duplicates()
		# Now get an index of just the iv dates so we can compare to the dv.index, which is also just dates
		iv_dates = iv_daily_et_index.map(lambda t: t.date())
		# missing days will be difference in the two indices, assuming dv is missing no dates
		missing_days = dv.index.difference(iv_dates)
		# subset the daily data to be only the days that are missing form instantaneous data
		dv_days = dv.loc[missing_days]
		
		# some logging messages to check my math... realized the index comparison was wrong before
		print(f"\nSpinup period is {(forecast_start - spinup).days} days long")
		print(f"There are {len(iv_dates)} dates with at least one instantaneous value in ca_data")
		print(f"Therefore there should be {(forecast_start - spinup).days - len(iv_dates)} daily value dates injected")
		print(f"In fact there were {len(dv_days)} dates injected")
		if len(missing_days) + len(iv_dates) != (forecast_start - spinup).days:
			warnings.warn("Faulty index comparison; number of missing dates plus exisiting dates does not add up to total timespan. Daily values could potentially have missing dates.", UserWarning)
		print()

		# option logger messages for debugging purposes 
		# print("Daily values to be injected into instantaneous timeseries:")
		# print(dv_days)
		# No inherent time zone info to daily values
		# so, reset tz to ET, set hour to Noon, convert tz to UTC
		dv_days.index = dv_days.index.map(lambda t: t.replace(hour=12, tzinfo=pytz.timezone('America/New_York')).tz_convert(dt.UTC))
		# now combine iv and dv data
		combined = pd.concat([iv, dv_days]).sort_index()
		# set the streamflow dict of ca_data to be the new combined series
		iv_dict['streamflow'] = combined

def getRichelieuLakeHt(start, end, service='dv'):
	# this is the date and time (in UTC) where "Current Observations" begin for the Richelieu gage
	# for Richelieu gage lake height calls prior to this date, we will use the "Daily Statistics" service ("dv")
	# https://waterdata.usgs.gov/nwis/inventory?agency_code=USGS&site_no=04295000
	rl_current_obs_begin = dt.datetime(2007, 10, 1, 5, tzinfo=dt.timezone.utc)
	rl = {"RL":'04295000'}
	lakeht = {'LAKEHT':'62614'}

	# This logic is only relevant if iv service Richelieu lake height data is requested; dv service request should be fine
	# if the end of the timeseries comes before the beginning of iv obs, use dv service
	if end < rl_current_obs_begin or service == 'dv':
		dv_data_chunks = []
		# use the "dv" service
		print(f"Requested date range ends before the beginning of the iv observations for Lake Height at the Richelieu gage ({rl_current_obs_begin})")
		print("Or you've explicitly requested daily means... either way getting daily means:")
		# current_start = start
		# counter = 1
		# while current_start < end:
		# 	print(f"ITER #: {counter}")
		# 	current_end = min(current_start + dt.timedelta(days=365), end)  # End of the 1-year chunk
			
		# 	# Fetch data for the current 1-year period
		# 	print("starting data get")
		# 	dv_data_chunk = usgs_ob.get_data(
		# 		start_date=current_start,
		# 		end_date=current_end,
		# 		locations=rl,
		# 		variables=lakeht,
		# 		service='dv')["RL"]["LAKEHT"]
		
		# 	# Add the chunk to the list
		# 	dv_data_chunks.append(dv_data_chunk)
    
		# 	# Update the start date for the next iteration
		# 	current_start = current_end
		# 	counter += 1
		# print("Out of the loop")
		# dv_data = {'RL':{"LAKEHT":pd.concat(dv_data_chunks)}}
		# print("data concatenated")
		dv_data = usgs_ob.get_data(start_date = start,
								end_date = end,
								locations = rl,
								variables = lakeht,
								service='dv')
		return dv_data
	# don't logically need the below end >= rl_current_obs_begin condition, since if the above if statement is false, then this one must be true
	elif end >= rl_current_obs_begin and start >= rl_current_obs_begin:
		# use the "iv" service
		iv_data = usgs_ob.get_data(start_date = start,
								end_date = end,
								locations = rl,
								variables = lakeht,
								service=service)
		return iv_data
	# this elif could just be an else, since if start >= rl_current_obs_begin is false, then this elif must always be true
	elif start < rl_current_obs_begin:
		print(f"Requested date range starts before the beginning of the iv observations for Lake Height at the Richelieu gage ({rl_current_obs_begin}) and ends after it")
		print("Getting both daily means and iv observations, then concating:")
		dv_data = usgs_ob.get_data(start_date = start,
								end_date = rl_current_obs_begin,
								locations = rl,
								variables = lakeht,
								service='dv')
		iv_data = usgs_ob.get_data(start_date = rl_current_obs_begin,
								end_date = end,
								locations = rl,
								variables = lakeht,
								service=service)
		# Make dv_series timezone-aware in order to combine it with iv_series;, hour also set to noon
		dv_data["RL"]["LAKEHT"].index = dv_data["RL"]["LAKEHT"].index.map(lambda x: x.replace(hour=12)).tz_localize('UTC')

		# Combine the series
		combined_data = {"RL": {"LAKEHT":pd.concat([dv_data["RL"]["LAKEHT"], dv_data["RL"]["LAKEHT"]])}}
		return combined_data

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
	if settings['hydrology_dataset_spinup'] == 'USGS_IV':
		observedHydro = usgs_ob.get_data(start_date = settings['spinup_date'] - dt.timedelta(days=1),
								 		 end_date = settings['forecast_start'],
										 locations = {"MS":'04294000',
			   									 	  "JS":'04292810',
			   										  "ML":'04292750'})
		# Convert USGS streamflow from cubic ft / s to cubic m / s
		for location in observedHydro.keys():
			observedHydro[location]['streamflow'] = observedHydro[location]['streamflow'].rename('streamflow (m³/s)') * 0.0283168

		# pull observed Rock and Pike flow data from candadian service site
		observedca = caflow_ob.get_data(start_date = settings['spinup_date'] - dt.timedelta(days=1),
								 		end_date = settings['forecast_start'],
										locations = {"PK":'030424',
			   										 "RK":'030425'})
		# backfill missing IV data with DV data - Right now we are just doing this for th spinup periop only
		backfillCaFlowsSpinup(observedca, settings)

		# add canadian reaches to observedHydro
		observedHydro.update(observedca)

	# this is a redundant value check... settings values are validated in get_args()
	# else: raise ValueError(f"'{settings['hydrology_dataset_spinup']}' is not a valid observational hydrology dataset")

	# the complete implementation of the HydroForcings class is forthcoming...
	# in the full implmentation, the below if statement logic will be unnecessary - you'll just need to call get_hydro_data()
	# Data adjustment methods, like backfillCaFlowsSpinup(), still need to be adapted for HydroForcings object
	if observedHydro is None:
		hydroSpinup = HydroForcings(start_date=settings['spinup_date'] - dt.timedelta(days=1),
								end_date=settings['forecast_start'],
								source=settings['hydrology_dataset_spinup'],
								period='spinup',
								dir=settings['data_dir'])
		hydroSpinup.get_hydro_data()

		observedHydro = hydroSpinup.access_data()

	forecastHydro = None

	# define reaches for NWM
	nwm_reaches = {"MS":"166176984",
				   "JS":"4587092",
				   "ML":"4587100"}
	
	if settings['hydrology_dataset_forecast'] == 'USGS_IV':
		forecastHydro = usgs_ob.get_data(start_date = settings['forecast_start'],
								 		 end_date = settings['forecast_end'] + dt.timedelta(days=1),
										 locations = {"MS":'04294000',
			   									 	  "JS":'04292810',
			   										  "ML":'04292750'})
		# Convert USGS streamflow from cubic ft / s to cubic m / s
		for location in forecastHydro.keys():
			forecastHydro[location]['streamflow'] = forecastHydro[location]['streamflow'].rename('streamflow (m³/s)') * 0.0283168

		forecastca = caflow_ob.get_data(start_date = settings['forecast_start'],
								 		end_date = settings['forecast_end'] + dt.timedelta(days=1),
										locations = {"PK":'030424',
			   										 "RK":'030425'})
		# update forecastHydro with Canadian data
		forecastHydro.update(forecastca)

	
	elif settings['hydrology_dataset_forecast'] == 'NOAA_NWM_PROD':
		forecastHydro = nwm_fc.get_data(forecast_datetime = settings['forecast_start'],
						end_datetime = settings['forecast_end'] + dt.timedelta(days=1),
						locations = nwm_reaches,
						forecast_type = settings['nwm_forecast_member'],
						data_dir = settings['data_dir'],
						load_threads = 1,
						google_buckets = True)
		# calculate streamflow for Rock and Pike based off Missisquoi streamflow
		forecastca = scaleRockandPikeQ(forecastHydro['MS']['streamflow'])
		# update forecastHydro with Canadian data
		forecastHydro.update(forecastca)
	
	elif settings['hydrology_dataset_forecast'] == 'NOAA_NWMV3_RETRO':
		forecastHydro = nwmv3_retro_fc.get_data(start_date = settings['forecast_start'],
										  		end_date = settings['forecast_end'] + dt.timedelta(days=1),
												locations = nwm_reaches)
		# calculate streamflow for Rock and Pike based off Missisquoi streamflow
		forecastca = scaleRockandPikeQ(forecastHydro['MS']['streamflow'])
		# update forecastHydro with Canadian data
		forecastHydro.update(forecastca)

	# Partial HydroForcings object implementation:
	if forecastHydro is None:
		hydroForecast = HydroForcings(start_date=settings['forecast_start'],
									  end_date=settings['forecast_end'] + dt.timedelta(days=1),
									  source=settings['hydrology_dataset_forecast'],
									  period='forecast',
									  dir=settings['data_dir'])
		hydroForecast.get_hydro_data()

		forecastHydro = hydroForecast.access_data()



			
	# Build Dictionary of Series with to adjusted column names
	iv_flows = {'MS':{},'ML':{},'JS':{},'PK':{},'RK':{}}	# initialize empty dictionary
	iv_flows['MS']['streamflow'] = pd.concat([observedHydro['MS']['streamflow'], forecastHydro['MS']['streamflow']]).rename_axis('time').astype('float')
	iv_flows['ML']['streamflow'] = pd.concat([observedHydro['ML']['streamflow'], forecastHydro['ML']['streamflow']]).rename_axis('time').astype('float')
	iv_flows['JS']['streamflow'] = pd.concat([observedHydro['JS']['streamflow'], forecastHydro['JS']['streamflow']]).rename_axis('time').astype('float')
	iv_flows['PK']['streamflow'] = pd.concat([observedHydro['PK']['streamflow'], forecastHydro['PK']['streamflow'] ]).rename_axis('time').astype('float')
	iv_flows['RK']['streamflow'] = pd.concat([observedHydro['RK']['streamflow'], forecastHydro['RK']['streamflow'] ]).rename_axis('time').astype('float')
	
	print('Structure of instantaneous flows dictionary')
	# print(iv_flows)

	logger.info(iv_flows)

	
	flow_data = {'Missisquoi':iv_flows['MS']['streamflow'],
				 'Mill':iv_flows['ML']['streamflow'],
				 'Jewett-Stevens':iv_flows['JS']['streamflow'],
				 'Pike':iv_flows['PK']['streamflow'],
				 'Rock':iv_flows['RK']['streamflow']}

	global SUBPLOT_PACKAGES
	global AXES
	SUBPLOT_PACKAGES['streamflow'] = {'labelled_data':flow_data,
								   	  'ylabel':'Streamflow (m/s^3)',
									  'title':f"{settings['hydrology_dataset_spinup']} vs. {settings['hydrology_dataset_forecast']}",
									  'fc_start':settings['forecast_start'],
									  'fc_end':settings['forecast_end'],
									  'row':0,
									  'col':0,
									  'axis':AXES}
	
	# save the daily values flow data in bay for later use in wqcalcs
	# THEBAY.flowdict = dv_flows			

	logger.info('Instantaneous Flow Data Scaled')
	# Scale Additional Inflows from Predicted Inflow
	#       sourcelist has list of water source IDs
	#       sourcemap defines the source name (wshed) and proportion and adjust of hydromodel output flow
	#		the adjust value is defined from InlandSeaModel_Notes Document describing how model calibration was done

	for baysource in THEBAY.sourcelist :
		logger.info('Generating Bay Source File for Id: '+baysource)
		wshed = THEBAY.sourcemap[baysource]['wshed']            # get column name of watershed flow source for this stream
		theflow = pd.Series()
		theflow['streamflow'] = iv_flows[wshed]['streamflow'] * THEBAY.sourcemap[baysource]['adjust']  # adjust gauge to loc
		theflow['streamflow'] =  theflow['streamflow'][theflow['streamflow'] >= 0]	# scrub out any negative values
		theflow['streamflow'].index = theflow['streamflow'].index.to_series().apply(datetimeToOrdinal) # ordinal dates index

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
			# output the ordinal date and adjusted flow value series
			theflow['streamflow'].to_csv(path_or_buf = output_file, float_format='%.3f',
				sep=' ', index=True, header=False)


##
#       End of Flow Data Import
#
##

#########################################################################
#
#   Import Hydrology Daily Flow - Generate Daily Values Streamflow Timeseries
#
##########################################################################

def getdailyflows(whichbay, settings):
	'''
		getdailyflows: Get daily flows for water quality calculations.
			dv_flows is stored in whichbay.flowdict
	'''

	logger.info("Loading Daily Streamflow")

	THEBAY = whichbay

	# define empty daily flows object - this should be final timeseries to be assigned to whichbay.flowdict
	dv_flows = {}

	# Assign US and canadian gauges
	us_gauges = {"MS":'04294000',
				 "JS":'04292810',
				 "ML":'04292750'}
	ca_gauges = {"PK":'030424',
				 "RK":'030425'}
	all_gauges = us_gauges | ca_gauges

	logger.info("Getting daily streamflow for spinup period...")
	logger.info(f"Dataset: {settings['hydrology_dataset_spinup']}")

	# getting spinup period daily streamflows
	spinupDailyFlows = None

	# pad the spinup date to ensure we get data on the specified spinup date
	adjusted_spinup_date = settings['spinup_date'] - dt.timedelta(days=1)

	if settings['hydrology_dataset_spinup'] == 'USGS_IV':
		# get spinup daily values from USGS
		usgs_spinup_dv = usgs_ob.get_data(start_date = adjusted_spinup_date,
										  end_date = settings['forecast_start'],
										  locations = us_gauges,
										  service = 'dv')
		
		# Convert USGS streamflow from cubic ft / s to cubic m / s
		for location in usgs_spinup_dv.keys():
			usgs_spinup_dv[location]['streamflow'] = usgs_spinup_dv[location]['streamflow'].rename('streamflow (m³/s)') * 0.0283168

		# get spinup daily values from canada
		caflow_spinup_dv = caflow_ob.get_data(start_date = adjusted_spinup_date,
											  end_date = settings['forecast_start'],
											  locations = ca_gauges,
											  service = 'dv')
		
		# combine spinup daily values into one nested dict
		spinupDailyFlows = usgs_spinup_dv | caflow_spinup_dv

		# set timestamps to noon ET, then convert to UTC
		for loc, q_dict in spinupDailyFlows.items():
			q = q_dict['streamflow']
			q.index = q.index.map(lambda t: t.replace(hour=12, tzinfo=None).tz_localize('America/New_York').tz_convert('UTC'))
			q_dict['streamflow'] = q

	if spinupDailyFlows is None:
		hydroSpinupDv = HydroForcings(start_date=adjusted_spinup_date,
									  end_date=settings['forecast_start'],
									  source=settings['hydrology_dataset_spinup'],
									  period='spinup',
									  dir=settings['data_dir'],
									  service='dv')
		hydroSpinupDv.get_hydro_data()

		spinupDailyFlows = hydroSpinupDv.access_data()

	logger.info("Spinup Period Daily Streamflow:")
	logger.info(spinupDailyFlows)

	logger.info("Getting daily streamflow for forecast period...")
	logger.info(f"Dataset: {settings['hydrology_dataset_forecast']}")


	# getting forecast period daily streamflows
	forecastDailyFlows = None

	# pad the forecast end date to ensure we get data on the specified forecast end date
	adjusted_forecast_enddate = settings['forecast_end'] + dt.timedelta(days=1)

	if settings['hydrology_dataset_forecast'] == 'USGS_IV':
		# get forecast daily values from USGS
		usgs_forecast_dv = usgs_ob.get_data(start_date = settings['forecast_start'],
											end_date = adjusted_forecast_enddate,
											locations = us_gauges,
											service = 'dv')
		
		# Convert USGS streamflow from cubic ft / s to cubic m / s
		for location in usgs_forecast_dv.keys():
			usgs_forecast_dv[location]['streamflow'] = usgs_forecast_dv[location]['streamflow'].rename('streamflow (m³/s)') * 0.0283168
		
		# get forecast daily values from canada
		caflow_forecast_dv = caflow_ob.get_data(start_date = settings['forecast_start'],
												end_date = adjusted_forecast_enddate,
												locations = ca_gauges,
												service = 'dv')
		
		# combine forecast daily values into one nested dict
		forecastDailyFlows = usgs_forecast_dv | caflow_forecast_dv

		# set timestamps to noon ET, then convert to UTC
		for loc, q_dict in forecastDailyFlows.items():
			q = q_dict['streamflow']
			q.index = q.index.map(lambda t: t.replace(hour=12, tzinfo=None).tz_localize('America/New_York').tz_convert('UTC'))
			q_dict['streamflow'] = q
		
	elif(settings['hydrology_dataset_forecast'] == 'NOAA_NWM_PROD'):
		# get forecast streamflow from NWM
		nwm_forecast = nwm_fc.get_data(forecast_datetime = settings['forecast_start'],
									   end_datetime = adjusted_forecast_enddate,
									   locations = {"MS":"166176984",
													"JS":"4587092",
													"ML":"4587100"},
									   forecast_type = settings['nwm_forecast_member'],
									   data_dir = settings['data_dir'],
									   load_threads = 1,
									   google_buckets = True)
		
		#	approximate Rock and Pike forecast from scaled Missisquoi forecast
		# TODO: Devise a better plan to reference the scaling factors than hardcoding them here
		nwm_forecast_ca = {'RK':{},'PK':{}}
		nwm_forecast_ca['RK']['streamflow'] = nwm_forecast['MS']['streamflow'] * 0.038  # scale Missisquoi to Rock
		nwm_forecast_ca['PK']['streamflow'] = nwm_forecast['MS']['streamflow'] * 0.228  # scale Missisquoi to Pike

		# combine forecast daily values into one nested dict
		forecastDailyFlows = nwm_forecast | nwm_forecast_ca

		# creating daily averages for NWM streamflow
		# NOTE: Do not convert timestamps to ET before daily averaging. Just leave it in UTC
		for loc, q_dict in forecastDailyFlows.items():
			q = q_dict['streamflow']
			# converting NWM timestamps from UTC to Eastern time zone
			# q.index = q.index.tz_convert('America/New_York')
			# daily streamflow averages
			q_daily_ave = q.resample('D').mean().rename('streamflow (m³/s)')
			# set hour to noon
			q_daily_ave.index = q_daily_ave.index.map(lambda t: t.replace(hour=12))
			# set hour to noon (eastern time) and then convert to UTC time once again
			# q_daily_ave.index = q_daily_ave.index.map(lambda t: t.replace(hour=12).tz_convert(dt.UTC))
			# Now drop the row for the date prior to forecast_start, which was created when converting UTC timestamps to ET
			# adding timezone info to forecast_start in order to compare with timestamps
			# q_dict['streamflow'] = q_daily_ave.loc[settings['forecast_start'].replace(tzinfo=dt.UTC)]
			q_dict['streamflow'] = q_daily_ave

	if forecastDailyFlows is None:
		hydroForecastDv = HydroForcings(start_date=settings['forecast_start'],
										end_date=adjusted_forecast_enddate,
									  	source=settings['hydrology_dataset_forecast'],
									  	period='forecast',
									  	dir=settings['data_dir'],
									  	service='dv')
		hydroForecastDv.get_hydro_data()

		forecastDailyFlows = hydroForecastDv.access_data()

	logger.info("Forecast Period Daily Streamflow:")
	logger.info(forecastDailyFlows)

	# now combine spinup and forecast period flows
	# this will combine spinup and forecast period series and fix timestamps
	for loc in all_gauges.keys():
		# Build Dictionary of Series with adjusted column names
		dv_combined = pd.concat([spinupDailyFlows['MS']['streamflow'], forecastDailyFlows['MS']['streamflow']]).rename_axis('time').astype('float')
		# remove any duplicate timestamps
		dv_combined = dv_combined[~dv_combined.index.duplicated(keep='first')]
		dv_dict = {'streamflow':dv_combined}
		dv_flows[loc] = dv_dict
	
	logger.info("Final Daily Streamflow Timeseries")
	logger.info(dv_flows)

	# set flowdict to be the daily flows dictionary
	THEBAY.flowdict = dv_flows

#########################################################################
#
#   All Things Meteorology Related Handled Herein
#
##########################################################################

# def adjustCRTemp(air_data):
# 	# S.E.T. - 20240206 - Adjust Colchester Reef Observed temp by month for the Missisquoi Bay Zone (401)
# 	# the adjustment, by month, to apply to CR data for use in MissBay
# 	tempadjust = {
# 		1 : -2.6,
# 		2 : -2.4,
# 		3 : -0.7,
# 		4 : 0.8,
# 		5 : 1.3,
# 		6 : 0.8,
# 		7 : -0.4,
# 		8 : -0.8,
# 		9 : -1.0,
# 		10 : -1.2,
# 		11 : -1.6,
# 		12 : -2.4 }
# 	print("line 715")
# 	adjustedCRtemp = pd.Series()
# 	print("air_data shape:")
# 	print(air_data.shape)
# 	for row in range(air_data.shape[0]):
# 		print(row)
# 		adjustedCRtemp[row] = air_data.iloc[row] + tempadjust[air_data.index[row].month]
# 	print("out of loop")
# 	adjustedCRtemp = adjustedCRtemp.set_axis(air_data.index)
# 	print("line 720")
# 	return adjustedCRtemp

def adjustCRTemp(air_data):
    # S.E.T. - 20240206 - Adjust Colchester Reef Observed temp by month for the Missisquoi Bay Zone (401)
    # the adjustment, by month, to apply to CR data for use in MissBay
    # P.J.C. - 20250124 - Updated to use pandas vectorized math so performance is much improved
   
    tempadjust_dict = {
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
 
    tempadjust = pd.Series(tempadjust_dict)
    return air_data + tempadjust[air_data.index.month].set_axis(air_data.index)


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

		# print(cls.nudge_df)

		multiplier = cls.nudge_df.loc[dayofyear]['Ratio (MB/CR)'] # the multiplier for each data record
		#print('Multipliers ',multiplier)

		#swdownobj_nudged = swdownobj.reset_index(drop=True) * multiplier.reset_index(drop=True)
		swdownobj_nudged = swdownobj.multiply(multiplier.to_numpy()) 
 
		#swdownready = pd.Series(swdownobj_nudged.array,swdownobj.index)
		return swdownobj_nudged

def femcRhumGapfill(femc_rh, start, end, allowed_gap_size=dt.timedelta(hours=2), figname='FEMCRelHumGapFilled.png'):
	'''
	Identifies gaps in FEMC relative humidity series and fills them with LCD relative humidity data.
	No adjustments or smoothing methods are applied to the LCD data. A plot of the gap-filled timeseries is
	created and saved whenever function is called.

	Args:
	-- femc_rh (pd.Series) [req]: FEMC relative humidity series
	-- start (dt.Datetime) [req]: the expected start date of the series
	-- end (dt.Datetime) [req]: the expected end date of the series
	-- allowed_gap_size (dt.Timedelta) [opt]: the largest acceptable duration of time between any two adjacent indices
	-- fignmae (str) [opt]: name for the figure to be saved

	Returns:
	Saves a figure, called 'FEMCRelHumGapFilled.png' that visualies the gap-filled FEMC series.
	-- rh_processed (pd.Series): the gap-filled FEMC relative humidity series
	'''
	# define some parameters for LCD data call
	rhum = {'RH2':'HourlyRelativeHumidity'}
	btv = {"BTV":"72617014742"}

	# first off, drop Na's
	rh_processed = femc_rh.dropna()

	# it's possible entire series in empty after removing NAN's (gap spans entire series)
	# in that case... just get LCD data and return, simple as that
	# Even if femc_rh has just one non-Na value, the rest of the gap processing code should work as expected
	if rh_processed.empty:
		print(f"FEMC Relhum data is completely missing from {start.strftime('%m-%d-%Y %H%M%S')} to {end.strftime('%m-%d-%Y %H%M%S')}")
		print("Getting BTV airport relhum data instead...\n")
		rh_processed = lcd_ob.get_data(start, end, btv, rhum)['BTV']['RH2']
		# now make a plot of the data
		fig = utils.plot_ts([rh_processed], scale='days', labels=['LCD'], colors=['orange'])
		fig.savefig(figname)
		return rh_processed

	# now, get a list of gaps in the series (assuming the series is not empty after dropping Na's)
	gaps = utils.get_dt_index_gaps(rh_processed.index, start_t=start, end_t=end, max_acceptable_td=allowed_gap_size)

	# check if that gap list is empty; if not, process and fill the gaps
	if gaps:
		# get lcd rhum data for the time period
		lcd_rh = lcd_ob.get_data(start, end, btv, rhum)['BTV']['RH2']
		# empty list to hold all of the chunks that will be stitched into the FEMC timeseries
		gap_plugs = []
		for left, right in gaps:
			# get the lcd data to plug the gap with
			lcd_plug_data = lcd_rh.loc[left:right]
			# pad the lcd data with the gap-bordering values from the FEMC data
			plug = rh_processed.loc[left:right].combine_first(lcd_plug_data)
			# combine all the lcd chunks into a list of series
			gap_plugs.append(plug)

		# combine femc data with all of the lcd chunks
		rh_processed = rh_processed.combine_first(pd.concat(gap_plugs))
		# remove duplicate indices that can occur if a FEMC data point is the right end of one gap and left end of another
		rh_processed = rh_processed.loc[~rh_processed.index.duplicated(keep='first')]

		# come up with some paramters for the plot function
		# each chunk in gap_plugs is plotted as an individual series. We give them the same color and only one label
		# so that they appear to be the same continuous series in the plot.
		color_list = ['blue'] + ['orange' for _ in range(len(gap_plugs))]
		full_labels = ['FEMC', 'LCD'] + [None for _ in range(len(gap_plugs) - 1)]

		# want to only plot the portion of the FEMC timeseries that involves gap-filling
		gap_zone_start = gaps[0][0] - dt.timedelta(days=7)
		gap_zone_end = gaps[-1][-1] + dt.timedelta(days=7)
		full_data = [rh_processed.loc[gap_zone_start:gap_zone_end]] + gap_plugs

		# now make a plot
		fig = utils.plot_ts(full_data, scale='auto', labels=full_labels, colors=color_list)
		fig.savefig(figname)

		# quick catch to see if any Nan's were introduced with all the pandas shuffling above... they shouldn't have been
		if rh_processed.hasnans: warnings.warn("NaN's detected in FEMC RelHum series after Gap Filling")
		# check for any gaps remaining after filling
		any_gaps_left = utils.get_dt_index_gaps(rh_processed.index, start_t=start, end_t=end, max_acceptable_td=allowed_gap_size)
		if any_gaps_left:
			warnings.warn('Gaps still present in FEMC relhum after femcRhumGapfill():')
			print(any_gaps_left)
		else: print("FEMC Relhum gaps successfully filled.")
	else: print("No Gaps detected in FEMC Relhum data")
	return rh_processed

def adjustFEMCLCD(whichbay, dataset, start_dt, end_dt, figure_name):
	# define a function to set relative humidity values to 100 if greater than 100
	# seems to be a problem in observations prior to 6/6/2019 in colchesterReefFEMC/Z0080_CR_QAQC.csv
	cap_rhum_at_100 = lambda x: 100 if x > 100 else x

	# BTV rain adjustment
	dataset['403']['RAIN'] = remove_nas(dataset['403']['RAIN']) * 0.6096

	for zone in get_climate_zone_keys(dataset):
		# air temp and swr adjustments
		if zone == '401':
			dataset[zone]['T2'] = remove_nas(adjustCRTemp(dataset[zone]['T2']).rename('T2', inplace=True))

			ShortwaveNudger.initialize(os.path.join(whichbay.template_dir, 'SolarRadiationFactor_MB.csv'))
			dataset[zone]['SWDOWN'] = ShortwaveNudger.nudgeDF(dataset[zone]['SWDOWN'])
		else:
			dataset[zone]['T2'] = remove_nas(dataset[zone]['T2'].rename('T2', inplace=True))

		# wind speed adjustments
		if zone == '403':
			dataset[zone]['WSPEED'] = remove_nas(dataset[zone]['WSPEED']) * 0.75
		else:
			dataset[zone]['WSPEED'] = remove_nas(dataset[zone]['WSPEED']) * 0.65
		
		# Removing NAs for wind direction, relative humidity, and short-wave radiation
		dataset[zone]['WDIR'] = remove_nas(dataset[zone]['WDIR'])
		# identify and process gaps in FEMC Relhum data
		# gap size is set to 2 hours (default); it's a parameter we could tweak if needed
		dataset[zone]['RH2'] = femcRhumGapfill(dataset[zone]['RH2'], start_dt, end_dt, allowed_gap_size=dt.timedelta(hours=2), figname=figure_name)
		# remove Nan's and apply relhum cap function; Intentionally keeping this seperate from femcRhumGapfill() function
		dataset[zone]['RH2'] = remove_nas(dataset[zone]['RH2']).apply(cap_rhum_at_100)
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
	#lakeLevel_df = pd.DataFrame({'ordinaldate': THEBAY.flowdf['ordinaldate'],
	#							 'msflow': THEBAY.flowdf['msflow'],
	#							 'flowmean_07': THEBAY.flowdf['msflow'].rolling(window=7,min_periods=1).mean(),
	#							 'flowmean_30': THEBAY.flowdf['msflow'].rolling(window=30,min_periods=1).mean(),
	#							 'flowmean_60': THEBAY.flowdf['msflow'].rolling(window=60,min_periods=1).mean()},
	#							 index=pd.DatetimeIndex(THEBAY.flowdf.index, name='time')
	#							 )
	lakeLevel_df = pd.DataFrame({'ordinaldate': THEBAY.flowdict['MS']['streamflow'].index.to_series().apply(datetimeToOrdinal),
								 'msflow': THEBAY.flowdict['MS']['streamflow'],
								 'flowmean_07': THEBAY.flowdict['MS']['streamflow'].rolling(window=7,min_periods=1).mean(),
								 'flowmean_30': THEBAY.flowdict['MS']['streamflow'].rolling(window=30,min_periods=1).mean(),
								 'flowmean_60': THEBAY.flowdict['MS']['streamflow'].rolling(window=60,min_periods=1).mean()},
								 index=pd.DatetimeIndex(THEBAY.flowdict['MS']['streamflow'].index, name='time')
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
	lakeLevel_df['LakeLevel'] = 94.05887 + 0.007910834 * lakeLevel_df[dataset['403']['T2'].name] + \
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

def adjustNOAAFSProducts(whichbay, dataset, settings):
	# slice GFS data at forecast start date, so that AEM3D uses observed dataset up to forecast start date
	for zone in get_climate_zone_keys(dataset):
		# calculate relative humidity for CFS only, BEFORE adjusting other variables (TEMP namely)
		if settings['weather_dataset_forecast'] == 'NOAA_CFS':
			dataset[zone]['RH2'] = cfs_fc.calculate_rh(psfc = dataset[zone]['PRSFC'], q2 = dataset[zone]['SH2'], t2 = dataset[zone]['T2'])
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
	if settings['weather_dataset_spinup'] == 'NOAA_LCD+FEMC_CR':
		observedClimateBTV = lcd_ob.get_data(start_date = adjusted_spinup,
										end_date = settings['forecast_start'],
										locations = {"401":"72617014742"},
										variables = {'TCDC':'HourlySkyConditions',
					   								 'RAIN':'HourlyPrecipitation'})
		
		observedClimateCR = femc_ob.get_data(start_date = adjusted_spinup,
										end_date = settings['forecast_start'],
										locations = {'401':'ColReefQAQC'})
		
		###### FEMC+LCD DATA ADJUSTMENTS HERE ######
		# make additional location dictionaries here rather than call for the same location multiple times in get_data()'s
		for key in ['402', '403']:
			observedClimateBTV[key] = copy.deepcopy(observedClimateBTV['401'])
			observedClimateCR[key] = copy.deepcopy(observedClimateCR['401'])


		# For MB (zone 401), we actually want cloud cover from Franklin Airport
		# Franklin airport ID: 00152
		# common code prefix for vermont stations: 726170
		# try Franklin data grab - 7-day data grab might return empty json, so in that case, use BTV data
		
		# boolean switch - if true, we are using franklin aiport cloud data for the forecast period
		fso_cloud_forecast = True
		try:
			# trying Franklin grab, and if that doesn't work...
			logger.info("Trying to get Spinup Period cloud data from Franklin airport (72049400152)...")
			observedClimateFSO = lcd_ob.get_data(start_date = adjusted_spinup,
											end_date = settings['forecast_start'],
											locations = {"401":"72049400152"},
											variables = {'TCDC':'HourlySkyConditions'})
		except Exception as e:
			# use Burlington data
			logger.warning("Franklin airport cloud data grab for Spinup Period Failed.")
			logger.warning(e)
			logger.warning("Getting Spinup Period cloud data from Burlington instead...")
			fso_cloud_forecast = False

		logger.info("Observed TCDC BTV:")
		logger.info(observedClimateBTV['401']['TCDC'].info())

		# only need to combine franklin and burlington cloud data if franklin data grab worked
		# otherwise, no need to combine, will just use BTV cloud data
		if fso_cloud_forecast:
			# now overwrite cloud cover for 401 - combine franklin cloud cover with that from BTV to fill in data gaps
			observedClimateBTV['401']['TCDC'] = utils.combine_timeseries(primary = observedClimateFSO['401']['TCDC'],
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
		observedlake = getRichelieuLakeHt(start=adjusted_spinup,
										  end=settings['forecast_start'])
		# observedlake = usgs_ob.get_data(start_date = adjusted_spinup,
		# 						 		 end_date = settings['forecast_start'],
		# 								 locations = {"RL":'04295000'},
		# 								 variables = {'LAKEHT':'62614'})
		# adjust height reference to 93 ft and convert to meters
		observedlake['RL']['LAKEHT'] = (observedlake['RL']['LAKEHT']-93) * 0.3048
		# store observed lake height in bay object for later concat with predicted height
		observedClimate['300'] = observedlake['RL']

		# have to pass the expected start and end dates for the observed climate Series for femcRhumGapfill()
		observedClimate = adjustFEMCLCD(THEBAY, observedClimate, start_dt=adjusted_spinup, end_dt=settings['forecast_start'], figure_name='ObservedFEMCRelHumGapFilled.png')

	else:
		raise ValueError(f"'{settings['weather_dataset_spinup']}' is not a valid observational weather dataset")

	forecastClimate = {}

	# define lat/lon coords for bay zones (needed for GFS and CFS)
	zone_coords = {'401': (45.00, -73.25),
				   '402': (44.75, -73.25),
				   '403': (44.75, -73.25)}
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
										locations = {'401':'ColReefQAQC'})
		
		
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
			forecastClimateBTV['401']['TCDC'] = utils.combine_timeseries(primary = forecastClimateFSO['401']['TCDC'],
																secondary = forecastClimateBTV['401']['TCDC'],
																interval = '20min')
			logger.info("forecast TCDC FSO:")
			logger.info(forecastClimateBTV['401']['TCDC'].info())

		# cool new way to combine dictionaries (python >= 3.9)for zone, ds in climateObsCR.items():
		# Both FEMC and LCD data grabbers now need mathing location keys to work since they are being combine
		for zone in forecastClimateCR.keys():
			forecastClimate[zone] = forecastClimateCR[zone] | forecastClimateBTV[zone]

		# Richelieu data represents forecast lake height
		forecastlake = getRichelieuLakeHt(start=settings['forecast_start'],
										  end=adjusted_end_date)
		# forecastlake = usgs_ob.get_data(start_date = settings['forecast_start'],
		# 						 		end_date = adjusted_end_date,
		# 								locations = {"RL":'04295000'},
		# 								variables = {'LAKEHT':'62614'})
		# adjust height reference to 93 ft and convert to meters
		print("Out of getRichelieuLakeHt()...")
		forecastlake['RL']['LAKEHT'] = (forecastlake['RL']['LAKEHT']-93) * 0.3048
		# store observed lake height in bay object for later concat with predicted height
		forecastClimate['300'] = forecastlake['RL']

		# passing the expected start and end dates for the forecast climate series
		forecastClimate = adjustFEMCLCD(THEBAY, forecastClimate, start_dt=settings['forecast_start'], end_dt=adjusted_end_date, figure_name='ForecastFEMCRelHumGapFilled.png')

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
													locations = zone_coords,
													data_dir = settings['data_dir'],
													load_threads = 1)
			
			###### GFS DATA ADJUSTMENTS HERE ######
			# slice GFS data at forecast start date, so that AEM3D uses observed dataset up to forecast start date
			for zone in forecastClimate.keys():
				for var in forecastClimate[zone].keys():
					forecastClimate[zone][var] = forecastClimate[zone][var].loc[pd.to_datetime(settings['forecast_start'].replace(tzinfo=dt.timezone.utc)):]
			forecastClimate = adjustNOAAFSProducts(THEBAY, forecastClimate, settings)

	# Handling CFS as the weather dataset for the forecast
	elif settings['weather_dataset_forecast'] == 'NOAA_CFS':
		logger.info(f"Begin CFS get_data()")
		forecastClimate = cfs_fc.get_data(start_date = settings['forecast_start'],
										  end_date = settings['forecast_end'],
										  locations = zone_coords,
										  variables = {'T2':'Temperature_height_above_ground',
													   'TCDC':'Total_cloud_cover_entire_atmosphere_single_layer',
													   'U10':'u-component_of_wind_height_above_ground',
													   'V10':'v-component_of_wind_height_above_ground',
													   'SH2':'Specific_humidity_height_above_ground',
													   'RAIN':'Precipitation_rate_surface',
													   'SWDOWN':'Downward_Short-Wave_Radiation_Flux_surface',
													   'PRSFC':'Pressure_surface'},
										   data_dir = settings['data_dir'])
		####### CFS Data Adjustments ######
		forecastClimate = adjustNOAAFSProducts(THEBAY, forecastClimate, settings)


	logger.info('Climate Observed Data')
	logger.info(observedClimate)
	logger.info('Climate Forecast Data (Zone 401)')
	logger.info(forecastClimate)

	air_temp = {"401": pd.concat([observedClimate['401']['T2'],forecastClimate['401']['T2']]).rename('T2'),
				"402": pd.concat([observedClimate['402']['T2'],forecastClimate['402']['T2']]).rename('T2'),
				"403": pd.concat([observedClimate['403']['T2'],forecastClimate['403']['T2']]).rename('T2')}

	global SUBPLOT_PACKAGES
	global AXES
	SUBPLOT_PACKAGES['air temp'] = {'labelled_data':air_temp,
								   	'ylabel':'Air Temperature at 2m (C)',
									'title':f"{settings['weather_dataset_spinup']} vs. {settings['weather_dataset_forecast']}",
									'fc_start':settings['forecast_start'],
									'fc_end':settings['forecast_end'],
									'row':0,
									'col':1,
									'axis':AXES}
	
	logger.info("Air Temp for Zone 401")
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

		print(f"Air Temp for zone {zone}")
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
									'fc_end':settings['forecast_end'],
									'row':0,
									'col':2,
									'axis':AXES}
	
	SUBPLOT_PACKAGES['bay snow'] = {'labelled_data':bay_snow,
									'ylabel':'Snow Fall (m/day)',
									'title':'Snow Fall for Inland Sea',
									'fc_start':settings['forecast_start'],
									'fc_end':settings['forecast_end'],
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
										'fc_end':settings['forecast_end'],
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
										'fc_end':settings['forecast_end'],
										'row':1,
										'col':1,
										'axis':AXES}
	
	SUBPLOT_PACKAGES['wind direction'] = {'labelled_data':winddir,
										'ylabel':'Wind Direction at 10m (degrees clockwise from North)',
										'title':'Wind Direction',
										'fc_start':settings['forecast_start'],
										'fc_end':settings['forecast_end'],
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
										'fc_end':settings['forecast_end'],
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
										'fc_end':settings['forecast_end'],
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
									'fc_end':settings['forecast_end'],
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
	### NOTE: adding 12 to iterations so that the last forecast timestamp is 6/8 00:00 instead of 6/7 23:00 (given forecast start is 6/1 00:00)\
	### each AEM3D iteration is 5 minutes, so 12 iterations is an hour, which is how long we want to extend the forecast to get to a full 7 days
	iterations = int((settings['forecast_end'] - simstart).total_seconds() / AEM3D_DEL_T) + 12
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
	# forecastStartIter = int((settings['forecast_start'] - (settings['spinup_date'] + dt.timedelta(days=1))).total_seconds() / AEM3D_DEL_T) + 1
	forecastStartIter = int((settings['forecast_start'] - settings['spinup_date']).total_seconds() / AEM3D_DEL_T) + 1

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

	# get daily flows for later water quality calculations
	getdailyflows(theBay, settings)

	# generate climate files including lake levels
	genclimatefiles(theBay, settings)

	# generate salinity file
	gensalinefile(theBay)

	# generate boundary condition file
	genboundaryfile(theBay)
	
	# generate tracer files
	gentracerfiles(theBay)

	# generate the water quality files (waterquality.py script)
	genwqfiles(theBay, settings)

	# generate the datablock.xml file
	gendatablockfile(theBay, settings)

	# generate control file
	gencntlfile(theBay, settings)

	global FIG
	global SUBPLOT_PACKAGES
	# make_figure(SUBPLOT_PACKAGES)
	# FIG.savefig('weatherVars.png')

	return 0
