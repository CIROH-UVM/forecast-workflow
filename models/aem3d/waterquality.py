#!/usr/bin/env python3

#  Creates the AEM3D Lake Model Water Quality Related Files
#
#  Time Series
#     Nutrient Series
#     Chlorophyl Series
#
#  Control File

from lib import *
from models.islesRF.cqrf import add_features
import joblib
import glob
import os
from string import Template
import numpy as np
import pandas as pd
from .AEM3D import *
from .Aem3dForcings import *

'''
List of WQ control Files present in wq_files from AEM3D_prep_file
	carbon.dat
	cyano.dat
	fdiat.dat
	light.dat
	nitrogen.dat
	organic_matter.dat
	oxygen.dat
	phosphorus.dat
	ssol1.dat
	ssol2.dat
	ts_bottom_51.dat
	ts_surface_51.dat

Sample of files to create in infiles (one of each for each river source)
	Pike_River_WQ_C.dat     Constant for all sources
	Pike_River_WQ_DO.dat    Based on Bay Temp
	Pike_River_WQ_N.dat
	Pike_River_WQ_P.dat
	Pike_River_WQ_Phyto.dat
	Pike_River_WQ_Si.dat     Constant for all sources
	Pike_River_WQ_TSS.dat
Specifications per Dr. Clelia Luisa Marti : June 2020

DO Dissolved Oxygen (mg O /L)
Set dissolved oxygen equal to the saturated value for the water temperature
DO = (14.652 - 4.1022e-1 * WTR_TEMP + 7.991e-3 * (WTR_TEMP2) –
7.77774e-5 * (WTR_TEMP3))
(Rich 1973 from Henderson-Sellers 1984, Engineering Limnology)


Q Discharge (m^3/s)

TP Total Phosphorous (mg P /L)
PO4 Filterable Reactive Phosphorous (mg P /L)
DOPL Dissolved Organic Phosphorous (Labile) (mg P /L)
POPL Particulate Organic Phosphorous (Labile) (mg P /L)

Equation to estimate TP
TP = 0.011001 * z2 + 0.073104 * z + 0.091528
where z = (Q - 159.78)/132.64

PO4 = 0.1050 * TP
DOPL = 0.23275 * TP
POPL = 0.30875 * TP


TN Total Nitrogen (mg N /L)
NH4 Ammonium (mg N /L)
NO3 Nitrate (mg N /L)
DONL Dissolved Organic Nitrogen (Labile) (mg N /L)
PONL Particulate Organic Nitrogen (Labile) (mg N /L)

Equation to estimate TN
TN = 0.00407 * z2  + 0.12853 * z + 0.75675
where z = (Q-166.3734)/138.1476
NH4 = 0.1 * TN
NO3 = 0.4 * TN
DONL = 0.2375 * TN
PONL = 0.2375 * TN

DOCL Dissolved Organic Carbon (Labile) (mg C /L)
POCL Particulate Organic Carbon (Labile) (mg C /L)
DOCL = 4.750 (Set to a constant value)
POCL  = 4.750 (Set to a constant value)

SiO2 Silica (mg Si /L)
SiO2 = 2.5 (Set to a constant value)

SSOL1 Suspended Solids Group 1 (mg /L)
Equation to estimate SSOL1
SSOL1 = 17.7558 * z2 + 59.5663 * z + 48.0127
where z = (Q - 170.0903)/ 142.1835

'''

def gencarbonfile(theBay):
	#
	#   Carbon data time series
	#
	generate_file_from_template('wq_c.template.txt',
								'WQ_C.dat',
								theBay,
								{'source_id_list': '  '.join(theBay.sourcelist),
								 'firstdate': theBay.FirstDate,
								 'lastdate': theBay.LastDate
								})    


def gensilicafile(theBay):
	#
	#   Silica data time series
	#     constant value for all sources
	#     set start and stop dates in template that has all sources
	#
	generate_file_from_template('wq_si.template.txt',
								'WQ_Si.dat',
								theBay,
								{'source_id_list': '  '.join(theBay.sourcelist),
								 'firstdate': theBay.FirstDate,
								 'lastdate': theBay.LastDate
								})        


def genwqfiles (theBay, settings):

	THEBAY = theBay

	# print('That is some Quality Water, right there.')

	# print('Copy Bay Flow')
	#flowdf = theBay.flowdf.copy()       # get flow dataframe from bay object
	
	flows = theBay.flowdict    # get dictionary of series of flow data
   
	#
	#       Write the water temperature file for each source of the bay
	#
	for baysource in THEBAY.sourcelist :

		# Dissolved Oxygen based on Water Temp
		DO  = 14.652 \
			- 4.1022e-1  * THEBAY.wtr_temp_dict[baysource] \
			+ 7.991e-3   * np.power(THEBAY.wtr_temp_dict[baysource], 2) \
			- 7.77774e-5 * np.power(THEBAY.wtr_temp_dict[baysource], 3)
		
		bs_name = THEBAY.sourcemap[baysource]['name']
		filename = "WQ_" + bs_name + '_DO.dat'
		logger.info('Generating Bay Source DO File: '+filename)

		# open the file in output directory
		pathedfile = os.path.join(THEBAY.infile_dir, filename)
		with open(pathedfile, mode='w', newline='') as output_file:

			THEBAY.addfile(fname=filename)        # remember generated bay files

			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ THEBAY.bayid + '                         !\n')
			output_file.write('! Bay Source: ' + bs_name + '                         !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('1 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write('0    '+baysource+'\n')
			output_file.write('TIME      DO\n')

			# output the ordinal date and temp dataframe columns
			index_to_ordinal_date(DO).to_csv(path_or_buf = output_file, float_format='%.3f',
			sep=' ', index=True, header=False)

	##
	#
	#   Generate Phyto related series, FDIAT and CYANO
	#
	##

	#
	# Write Phyto Files - fixed series, use templates and set $year
	#
	for templateFile in [f for f in os.listdir(theBay.template_dir) if 'WQ_Phyto' in f]:
		generate_file_from_template(templateFile,
									os.path.splitext(templateFile)[0]+'.dat',
									theBay,
									{'year': theBay.FirstDate[0:4]})
	#
	# End of Phyto File Generation
	#

	#
	#   Generate Source Specific P, N, TSS files
	#

	for baysource in theBay.sourcelist :

		bs_name = theBay.sourcemap[baysource]['name']

		#
		# Bring in P reduction csv
		#         Adapted from RCA file prep bcseries.r

		# TODO: Check for p_reductions to be present, if not, provide sample where p_redux = 1.0
		
		print('NOT NOT NOT Reading P Reduction CSV')
		p_redux = 1.0

		'''
		#p_reductions = read.csv('p_reductions.csv', row.names='year', check.names = FALSE)
		p_reductions = pd.read_csv('p_reductions.csv', delimiter=',', index_col='year')
		#print('P Reduction CSV read')
		#print(p_reductions.columns)
		#print(p_reductions)

		# Get the one line from the p_reduction table for the current year
		theyears_p_reductions = p_reductions.loc[THEBAY.year]

		logger.info(' Selecting for reduxP = ', SCENARIO.reduxP)

		## Convert percent reduction for current scenario to a "percent of" and store (reduxP is scenario string)
		if (SCENARIO.reduxP == 'no_redux') :
				p_redux = (100 - theyears_p_reductions[['0_redux']]) / 100
		else:
				p_redux = (100 - theyears_p_reductions[[SCENARIO.reduxP]]) / 100
		'''
		#p_redux = (100 - theyears_p_reductions[['40_redux']]) / 100       # Fix This hardcoded redux reference
		# logger.info(f'Scaling Flow by P_Redux: {p_redux}')

		phosdf = pd.DataFrame()
		#phosdf['ordinaldate'] = flowdf['ordinaldate']
		phosdf['ordinaldate'] = flows[THEBAY.sourcemap[baysource]['wshed']]['streamflow'].index.to_series().apply(datetimeToOrdinal) # ordinal dates index

		# Get Q for use in equations
		# Don't use 'adjust' factor to calculate nutrient fluxes
		#   per 2024 06 email discussion with Peter Isles
		# ! Use Q (streamflow) estimates / observations
		#   at USGS gauge / nutrient monitoring locations
		Q = flows[THEBAY.sourcemap[baysource]['wshed']]['streamflow'] 
		# Clip low flows to a minimum for log calculations
		# TODO -  consider using a Bay resource variable to clamp min flows (as suggested below)
		#logQnoz = Q.apply(lambda x: np.log10(THEBAY.sourcemap[baysource]['minCQflow']) if x < THEBAY.sourcemap[baysource]['minCQflow'] else np.log10(x))
		if bs_name.startswith('MissisquoiRiver'):
			logQ = Q.apply(lambda x: np.log10(10.0) if x < 10.0 else np.log10(x))
		else:
			logQ = Q.apply(lambda x: np.log10(0.1) if x < 0.1 else np.log10(x))

		# Changed default 20240607 by pjc

		# if phosdf.isna().any().any():
		# 	print(phosdf[phosdf.isna().any(axis=1)])
		# 	raise Exception(f"NA's detected in phos df for bs_name: {bs_name}")

		#TODO: Implement BREE2021Seg
		#TODO: Move all this junk to THEBAY, choose cqVersion at THEBAY creation
		# 8/19/2024 - cqVersion is now a attribute of THEBAY object and is determined by SETTINGS
		cqVersion = THEBAY.cqVersion

		# I think it's appropriate to have this be hard-coded, as it shouldn't change... unless someone runs our workflow outside of the VACC
		rf_models_dir = '/netfiles/ciroh/models/cq_randomforest_peterisles/current/'
		# rf_models_dir = '/netfiles/ciroh/nbeckage/rfcq/'

		print(f"Q DATAFRAME:")
		print(Q)
		# make wQ a settings, have peter's stuff be an option
		if cqVersion == 'Clelia':
			# Clelia TP Concentration - Discharge Relationship
			#   Same for ALL ILS inputs!!
			zTP = (flows['MS']['streamflow'] - 159.78) / 132.64
			phosdf['TP'] =  ( 0.011001 * np.power(zTP,2) + 0.073104 * zTP + 0.091528 ) * p_redux
		elif cqVersion == 'islesRF':
			rf_tp_dir = os.path.join(rf_models_dir, "TP")
			Q_features = add_features(pd.DataFrame(Q))
			if bs_name.startswith('MissisquoiRiver'):
				rf_tp = joblib.load(os.path.join(rf_tp_dir, 'missisquoi_TP.joblib'))
				phosdf['TP'] = rf_tp.predict(Q_features)
			elif bs_name.startswith('RockRiver'):
				rf_tp = joblib.load(os.path.join(rf_tp_dir, 'rock_TP.joblib'))
				phosdf['TP'] = rf_tp.predict(Q_features)
			elif bs_name.startswith('PikeRiver'):
				rf_tp = joblib.load(os.path.join(rf_tp_dir, 'pike_TP.joblib'))
				phosdf['TP'] = rf_tp.predict(Q_features)
			elif bs_name.startswith('MillRiver'):
				rf_tp = joblib.load(os.path.join(rf_tp_dir, 'mill_TP.joblib'))
				phosdf['TP'] = rf_tp.predict(Q_features)
			elif bs_name.startswith('JewettStevens'):
				rf_tp = joblib.load(os.path.join(rf_tp_dir, 'jewett_TP.joblib'))
				phosdf['TP'] = rf_tp.predict(Q_features)
	
		elif cqVersion == '202406Calibration':
			if bs_name.startswith('MissisquoiRiver'):
				zTP = (Q - 159.78) / 132.64
				phosdf['TP'] = 0.091528 + 0.073104*zTP + 0.011001*(zTP*zTP) * p_redux
			elif bs_name.startswith('RockRiver'):
				phosdf['TP'] = 0.065 + 0.0115*Q - 0.000168*(Q*Q) * p_redux
			elif bs_name.startswith('PikeRiver'):
				phosdf['TP'] = 0.016 + 0.0016*Q - 0.0000015*(Q*Q) * p_redux
			# 2024 Calibration used VT DEC data for Mill and JS
			#   For now, these are the BREE2021Quad Equations
			#   But, we should switch to Takis' new 2024 equations ASAP
			elif bs_name.startswith('MillRiver'):
				phosdf['TP'] = np.power(10, (1.7935 + 0.4052 * logQ + 0.1221 * logQ * logQ)) * p_redux / 1000
			elif bs_name.startswith('JewettStevens'):
				phosdf['TP'] = np.power(10, (2.2845 + 0.5185 * logQ + 0.1995 * logQ * logQ)) * p_redux / 1000
			else:
				raise Exception(f'CQ Equation for baysource={bs_name} not found for cqVersion={cqVersion}') 
			
		elif cqVersion == 'BREE2021Quad':
			# Takis' new quadratic CQ equations by stream derived from historical data
			# TODO Fix these to use Q from above instead of raw flowdf[]
			if (
			bs_name.startswith('MissisquoiRiver') or 
			bs_name.startswith('RockRiver') or
			bs_name.startswith('PikeRiver')
			):
				phosdf['TP'] = np.power(10, (1.6884 - 0.7758 * logQ + 0.3952 * logQ * logQ)) * p_redux / 1000
			elif bs_name.startswith('JewettStevens'):
				phosdf['TP'] = np.power(10, (2.2845 + 0.5185 * logQ + 0.1995 * logQ * logQ)) * p_redux / 1000
			elif bs_name.startswith('MillRiver'):
				phosdf['TP'] = np.power(10, (1.7935 + 0.4052 * logQ + 0.1221 * logQ * logQ)) * p_redux / 1000
			else:
				raise Exception(f'CQ Equation for baysource={bs_name} not found for cqVersion={cqVersion}')
		elif cqVersion.startswith("READ_CSV"):
			### NOTE: NutrientForcings class NOT implmented -- that will require restructuring a lot of genwqfiles()
			tpForcings = NutrForcings(start_date=ordinalToDatetime(THEBAY.FirstDate),
								 		  end_date=ordinalToDatetime(THEBAY.LastDate),
										  source=cqVersion.split(":")[0],
										  dir=cqVersion.split(":")[-1])
			tp_csv_dir = os.path.join(settings["data_dir"], cqVersion.split(":")[-1], "TP")
			if bs_name.startswith("Missisquoi"):
				fname = "Conc_TP_Missisquoi"
				phosdf['TP'] = pd.read_csv(os.path.join(tp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Rock"):
				fname = "Conc_TP_Rock"
				phosdf['TP'] = pd.read_csv(os.path.join(tp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Pike"):
				fname = "Conc_TP_Pike"
				phosdf['TP'] = pd.read_csv(os.path.join(tp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Mill"):
				fname = "Conc_TP_Mill"
				phosdf['TP'] = pd.read_csv(os.path.join(tp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Jewett"):
				fname = "Conc_TP_Jewett"
				phosdf['TP'] = pd.read_csv(os.path.join(tp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
		else: raise Exception(f'cqVersion {cqVersion} not defined')

		# Check for NAs -- will throw an error and log data if NA's are returned
		#   Gonna keep this code in case we need to debug similar issues again
		if phosdf['TP'].isna().any():
			logger.info("Q indices that became NA")
			logger.info(Q[phosdf.isna().any(axis=1)])
			logger.info("logQ indices that became NA's:")
			logger.info(logQ[phosdf.isna().any(axis=1)])
			logger.info("phosdf with NA's:")
			logger.info(phosdf[phosdf.isna().any(axis=1)])
			e = Exception(f"NA's detected in phosdf['TP'] bs_name: {bs_name}")
			logger.exception(e)
			raise e        

		# for other phosphorus concentrations, islesRF should use the same formulas as 202406Calibration
		if cqVersion == '202406Calibration' or cqVersion == 'islesRF':
			if bs_name.startswith('MissisquoiRiver'):
				phosdf['PO4'] = 0.01401 * np.power(phosdf['TP'],0.319)
				phosdf['DOPL'] = 0.03268 * np.power(phosdf['TP'],0.319)
				phosdf['POPL'] = phosdf['TP'] - phosdf['PO4'] - phosdf['DOPL']
			elif(
			bs_name.startswith('RockRiver') or
			bs_name.startswith('PikeRiver') or
			bs_name.startswith('JewettStevens') or
			bs_name.startswith('MillRiver')
			):
				phosdf['PO4'] = 0.16 * phosdf['TP']
				phosdf['DOPL'] = 0.36 * phosdf['TP']
				phosdf['POPL'] = 0.48 * phosdf['TP']
			else:
				raise Exception(f'baysource={bs_name} not found when calculating phosphorus species for 202406Calibration')
		
		elif cqVersion.startswith("READ_CSV"):
			dp_csv_dir = os.path.join(settings["data_dir"], cqVersion.split(":")[-1], "DP")
			if bs_name.startswith("Missisquoi"):
				fname = "Conc_DP_Missisquoi"
				phosdf['DP'] = pd.read_csv(os.path.join(dp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Rock"):
				fname = "Conc_DP_Rock"
				phosdf['DP'] = pd.read_csv(os.path.join(dp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Pike"):
				fname = "Conc_DP_Pike"
				phosdf['DP'] = pd.read_csv(os.path.join(dp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Mill"):
				fname = "Conc_DP_Mill"
				phosdf['DP'] = pd.read_csv(os.path.join(dp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Jewett"):
				fname = "Conc_DP_Jewett"
				phosdf['DP'] = pd.read_csv(os.path.join(dp_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			
			# These equations come from email communication Peter Isles', 1/22/2025 (RE:Status of longer-term AEM3D Hindcasts)
			phosdf['PO4'] = phosdf['DP'] * 0.3
			phosdf['DOPL'] = phosdf['DP'] * 0.7
			print("phosdf['TP']")
			print(phosdf['TP'])
			print("phosdf['DP']")
			print(phosdf['DP'])
			# Replace values with 0 where the below condition is False (i.e. where DP > TP)
			phosdf['POPL'] = (phosdf['TP'] - phosdf['DP']).where(phosdf['DP'] <= phosdf['TP'], 0)

		# Same for cqVersion = Clelia or BREE2021
		#   Updated 2021.05.27 per WQS Docs
		else:
			if bs_name.startswith('MissisquoiRiver'):
				phosdf['PO4'] = 0.01401 * np.power(phosdf['TP'],0.319)
				phosdf['DOPL'] = 0.03268 * np.power(phosdf['TP'],0.319)
				phosdf['POPL'] = phosdf['TP'] - phosdf['PO4'] - phosdf['DOPL']
			elif(
			bs_name.startswith('RockRiver') or
			bs_name.startswith('PikeRiver') or
			bs_name.startswith('JewettStevens') or
			bs_name.startswith('MillRiver')
			):
				phosdf['PO4'] = 0.1050 * phosdf['TP']
				phosdf['DOPL'] = 0.23275 * phosdf['TP']
				phosdf['POPL'] = 0.30875 * phosdf['TP']
			else:
				raise Exception(f'baysource={bs_name} not found when calculating phosphorus species for {cqVersion}')

		print('Check of PO4 dataframe')
		print(phosdf['PO4'])

		# Nitrogen Series
		nitdf = pd.DataFrame()
		#nitdf['ordinaldate'] = flowdf['ordinaldate']
		nitdf['ordinaldate'] = flows[THEBAY.sourcemap[baysource]['wshed']]['streamflow'].index.to_series().apply(datetimeToOrdinal) # ordinal dates index

		if cqVersion == '202406Calibration':
			if bs_name.startswith('PikeRiver'):
				nitdf['TN'] = 0.95 + 0.004*Q + 0.7*np.power(2.72, (-Q/11))
			elif bs_name.startswith('RockRiver'):
				nitdf['TN'] = 2.0 + 0.021*Q + 1.45*np.power(2.72, (-(Q+0.9)/1.6))
			elif (
			bs_name.startswith('MissisquoiRiver') or
			bs_name.startswith('JewettStevens') or
			bs_name.startswith('MillRiver')
			):
				zTN = (Q - 166.3734) / 138.1476
				nitdf['TN'] = 0.00407 * (zTN*zTN)  + 0.12853 * zTN + 0.75675            
		
		elif cqVersion == 'islesRF':
			rf_tn_dir = os.path.join(rf_models_dir, "TN")
			if bs_name.startswith('MissisquoiRiver'):
				rf_tn = joblib.load(os.path.join(rf_tn_dir, 'missisquoi_TN.joblib'))
				nitdf['TN'] = rf_tn.predict(Q_features)
			elif bs_name.startswith('RockRiver'):
				rf_tn = joblib.load(os.path.join(rf_tn_dir, 'rock_TN.joblib'))
				nitdf['TN'] = rf_tn.predict(Q_features)
			elif bs_name.startswith('PikeRiver'):
				rf_tn = joblib.load(os.path.join(rf_tn_dir, 'pike_TN.joblib'))
				nitdf['TN'] = rf_tn.predict(Q_features)
			elif bs_name.startswith('MillRiver'):
				rf_tn = joblib.load(os.path.join(rf_tn_dir, 'mill_TN.joblib'))
				nitdf['TN'] = rf_tn.predict(Q_features)
			elif bs_name.startswith('JewettStevens'):
				rf_tn = joblib.load(os.path.join(rf_tn_dir, 'jewett_TN.joblib'))
				nitdf['TN'] = rf_tn.predict(Q_features)

		elif cqVersion.startswith("READ_CSV"):
			tn_csv_dir = os.path.join(settings["data_dir"], cqVersion.split(":")[-1], "TN")
			if bs_name.startswith("Missisquoi"):
				fname = "Conc_TN_Missisquoi"
				nitdf['TN'] = pd.read_csv(os.path.join(tn_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Rock"):
				fname = "Conc_TN_Rock"
				nitdf['TN'] = pd.read_csv(os.path.join(tn_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Pike"):
				fname = "Conc_TN_Pike"
				nitdf['TN'] = pd.read_csv(os.path.join(tn_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Mill"):
				fname = "Conc_TN_Mill"
				nitdf['TN'] = pd.read_csv(os.path.join(tn_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
			if bs_name.startswith("Jewett"):
				fname = "Conc_TN_Jewett"
				nitdf['TN'] = pd.read_csv(os.path.join(tn_csv_dir, f"{fname}.csv"), parse_dates=True, index_col='Date')[fname]
		else:
			zTN = (Q - 166.3734) / 138.1476
			nitdf['TN'] = 0.00407 * (zTN*zTN)  + 0.12853 * zTN + 0.75675
			
		# islesRf and READ_CSV (for Peter's CSVs) should use the same equations as 202406Calibration for all other nitrogren species
		if cqVersion == '202406Calibration' or cqVersion == 'islesRF' or cqVersion.startswith("READ_CSV"):
			# using the same equations for nitrogen species for every trib
			nitdf['NH4'] = 0.1 * nitdf['TN']
			nitdf['NO3'] = 0.3 * nitdf['TN']
			nitdf['DONL'] = 0.1 * nitdf['TN']
			nitdf['PONL'] = 0.5 * nitdf['TN']   
		else:
			nitdf['NH4'] = 0.1 * nitdf['TN']            # Updated 2021.05.27 per WQS Docs
			nitdf['NO3'] = 0.3 * nitdf['TN']
			nitdf['DONL'] = 0.1 * nitdf['TN']
			nitdf['PONL'] = 0.475 * nitdf['TN']

		# Suspended Solids Series
		ssdf = pd.DataFrame()
		#ssdf['ordinaldate'] = flowdf['ordinaldate']
		ssdf['ordinaldate'] = flows[THEBAY.sourcemap[baysource]['wshed']]['streamflow'].index.to_series().apply(datetimeToOrdinal) # ordinal dates index

		zSS = (Q - 170.0903) / 142.1835
		ssdf['SSOL1'] = 17.7558 * np.power(zSS,2)  + 59.5663 * zSS + 48.0127

		#
		# Write Phosophorus File
		#
		filename = bs_name + '_WQ_P.dat'
		logger.info('Generating Bay Source File: '+filename)

		# open the file in output directory
		pathedfile = os.path.join(theBay.infile_dir, filename)
		with open(pathedfile, mode='w', newline='') as output_file:

			theBay.addfile(fname=filename)    # remember generated file names

			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ theBay.bayid + '                         !\n')
			output_file.write('! Bay Source: ' + bs_name + '                         !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('3 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write('0    ' + baysource + '  ' + baysource + '  ' + baysource + '\n')
			output_file.write('TIME      PO4	DOPL	POPL\n')

			# output the ordinal date and P values time dataframe columns
			phosdf.to_csv(path_or_buf = output_file, columns= ['ordinaldate', 'PO4', 'DOPL', 'POPL'], float_format='%.3f',
			sep=' ', index=False, header=False)

		#
		# Write Nitrogen File
		#
		filename = bs_name + '_WQ_N.dat'
		logger.info('Generating Bay Source File: '+filename)

		# open the file in output directory
		pathedfile = os.path.join(theBay.infile_dir, filename)
		with open(pathedfile, mode='w', newline='') as output_file:

			theBay.addfile(fname=filename)    # remember generated file names

			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ theBay.bayid + '                         !\n')
			output_file.write('! Bay Source: ' + bs_name + '                         !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('4 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write('0    ' + baysource + '  ' + baysource + '  ' + baysource + '  ' + baysource +'\n')
			output_file.write('TIME     NH4	NO3	DONL	PONL\n')

			# output the ordinal date and N values time dataframe columns
			nitdf.to_csv(path_or_buf = output_file, columns= ['ordinaldate', 'NH4', 'NO3', 'DONL', 'PONL'], float_format='%.3f',
			sep=' ', index=False, header=False)

		#
		# Write Total Suspended Solids File
		#
		filename = bs_name + '_WQ_TSS.dat'
		logger.info('Generating Bay Source File: '+filename)

		# open the file in output directory
		pathedfile = os.path.join(theBay.infile_dir, filename)
		with open(pathedfile, mode='w', newline='') as output_file:

			theBay.addfile(fname=filename)    # remember generated file names

			# output the header text
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('! Written by AEM3D_prep_IAM                           !\n')
			output_file.write('! Bay ID: '+ theBay.bayid + '                         !\n')
			output_file.write('! Bay Source: ' + bs_name + '                         !\n')
			output_file.write('!-----------------------------------------------------!\n')
			output_file.write('1 data sets\n')
			output_file.write('0 seconds between data\n')
			output_file.write('0    ' + baysource +'\n')
			output_file.write('TIME     SSOL1 \n')

			# output the ordinal date and TSS value time dataframe columns
			ssdf.to_csv(path_or_buf = output_file, columns= ['ordinaldate', 'SSOL1'], float_format='%.3f',
			sep=' ', index=False, header=False)

	gencarbonfile(theBay)   # one file for constant carbon data series
	gensilicafile(theBay)   # one file for constant silica data series
