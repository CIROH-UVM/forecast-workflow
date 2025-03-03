from collections import deque
import datetime as dt
import logging
from models.aem3d.AEM3D import ordinalToDatetime
import os
import pandas as pd
import warnings
import xarray as xr
import json
import argparse

'''
This purpose of this script is to Loop through some FEE runs (specififed by the paramters below) and report any failed runs.
If a failed run is detected, the error is diagnosed and reported. This script creates two files by default; (1) a comprehensive log file that
reports error messages detected in failed runs, and (2) a concise csv summary file of failed runs and their causes.
The script is designed to be run in a debugging environment, and is not intended for production use.

Settings in the debug_params.json file:
- root_dir (str): the root directory of the FEE runs
- cq_paradigms (list of str): a list of CQ paradigms to loop through
- scenarios (str or list of str): 'all' or a list of scenarios to loop through
- dir_years (str or list of str): 'all' or a list of years to loop through
- launch_start (str): the start date of the FEE runs in MMDD format
- launch_end (str): the end date of the FEE runs in MMDD format
- file_path (str): the path to the debug_params.json file
- last_lines_to_read (int): the number of lines to read from the end of the .out file
- log_name (str): the name of the log file to create

TO RUN THIS DEBUG SCRIPT:
1. Create a custom debug parameters json file with the settings above
2. Run the script with the --conf argument pointing to the debug parameters json file
3. The script will output a log file and a csv file with the error messages detected in failed runs
'''

def aem3d_fatal(fatal_msg, err_report):
	'''
    Function that parses error messages from AEM3D model runs, identifying fatal errors and extracting relevant details.

    Args:
    -- fatal_msg (list) [req]: A list of strings representing the error message lines from the AEM3D model output.

    Returns:
    -- err_report (dict): A dictionary containing details of the error with the following keys:
        - 'error' (str): Always set to 'aem3d_fatal'.
        - 'err_type' (str or None): The type of fatal error identified (e.g., "Insufficient data", "run_aem3d.dat missing"). Extend this function to parse more AEM3D error types as they come up.
        - 'file' (str or None): The file involved in the error (if applicable).
        - 'last_date_ord' (str or None): The last date in ordinal form before the error occurred (if applicable).
        - 'last_date_dt' (str or None): The last date in human-readable datetime format before the error occurred (if applicable).
    '''
	err_report['error'] = 'aem3d_fatal'
	for i, line in enumerate(fatal_msg):
		print(f'\t\t\t\t {line}')
		# add more if clauses to deal with other AEM3D fatal error types as they come up
		if "Insufficient data" in line:
			err_report['err_type'] = "Insufficient data"
			err_report['file'] = line.split('file ')[-1]
			ord_date = fatal_msg[i+1].split('read: ')[-1].split(' ')[0]
			err_report['last_date_ord'] = ord_date
			err_report['last_date_dt'] = ordinalToDatetime(ord_date).strftime("%m/%d/%Y")
			# print(err_report)
		elif "run_aem3d.dat missing" in line:
			err_report['err_type'] = "run_aem3d.dat missing"
			err_report['file'] = "run_aem3d.dat"
	return err_report

def get_timestamps(ds):
	'''
	Function that pulls the timestamps out of an "aem3d-run/outfiles/nc/sheet_surface_forecast.nc" file.
	May work for other AEM3D netcdf output files if they are structured similariliy to the surface sheets.

	Args:
	-- ds (xarray.Dataset) [req]: the dataset to get timestamps for

	Returns:
	-- timestamps (list): a list of datetime objects representing the timestamps of the netcdf file
	'''
	timestamps = []
	for i in range(len(ds.T)):
		ts = dt.datetime(year=int(ds['Year'].values[i]),
						month=int(ds['Month'].values[i]),
						day=int(ds['Day'].values[i]),
						hour=int(ds['Hour'].values[i]),
						minute=int(ds['Minute'].values[i]),
						second=int(ds['Second'].values[i]))
		timestamps.append(ts)
	return timestamps

def prep_IAM_err(prep_err, err_report):
	err_report['error'] = 'prep_IAM'
	prnt = False
	for i, line in enumerate(prep_err):
		if "During handling of the above exception" in line: prnt = True
		if prnt: print(f'\t\t\t\t {line}')
		if "AEM3D_prep_worker - None" in line: prnt = False
	return err_report

def parse_failed_run(filelines, report):
	'''
    Function that parses the lines of the .out file to detect what kind of error caused the run to fail.
	Passes on the .out file lines to more specific error parser functions once a braod kind of error is determined

    Args:
    -- filelines (list) [req]: A list of strings representing the lines of the .out file from a workflow run.

    Returns:
    -- report (dict or None): A dictionary containing the parsed fatal error report if a fatal error is found, or None
    	if no fatal error is detected. The structure of the report is the same as the output of the aem3d_fatal() function.
    '''
	# for line in filelines:
	# 	print(line)
	for i, line in enumerate(filelines):
		# check if AEM3D threw a fatal error
		if "***** FATAL ERROR ****" in line:
			fatal_error = [line.rstrip(' \n') for line in filelines[i:i+9] if line != '\n']
			report = aem3d_fatal(fatal_error, report)
			break
		if "AEM3D_prep_IAM.py failed." in line:
			prep_error = [line.rstrip(' \n') for line in filelines[i:] if line != '\n']
			report = prep_IAM_err(prep_error, report)
			break
	return report

def report_df(failure_dict):
	'''
    Function that converts a nested dictionary of failure data into a formatted pandas DataFrame. The DataFrame is structured
    by extracting key components and organizing fields for easier analysis and reporting.

    Args:
    -- failure_dict (dict) [req]: A nested dictionary containing failure data, where keys are hierarchical strings
       (e.g., 'cq_paradigm.scenario.year.date.field') and values are the corresponding data.

    Returns:
    -- df (pandas.DataFrame): A pandas DataFrame indexed by 'cq_paradigm', 'scenario', 'year', and 'date' with
       the fields 'error', 'err_type', 'file', 'last_date_ord', and 'last_date_dt' as columns. None if failure_dict contains no data
    '''
	# Use pandas json_normalize to flatten the dictionary
	df = pd.json_normalize(failure_dict, sep='.')
	# Transpose the DataFrame so that the keys form the columns
	df = df.T.reset_index()
	# if the df is not empty, continue
	if not df.empty:
		# Split the flattened keys into individual components
		df[['cq_paradigm', 'scenario', 'year', 'date', 'field']] = df['index'].str.split('.', expand=True)
		# Drop the old 'index' column
		df = df.drop(columns=['index'])
		# Rearrange the DataFrame for better readability
		df = df[['cq_paradigm', 'scenario', 'year', 'date', 'field', 0]]
		# Rename the last column to 'value'
		df.columns = ['cq_paradigm', 'scenario', 'year', 'date', 'Field', 'Value']
		# Pivot the DataFrame to have one row per 'Year' and 'Date', and fields as columns
		df = df.pivot_table(index=['cq_paradigm', 'scenario','year', 'date'], columns='Field', values='Value', aggfunc='first').reset_index()
		df.columns.name = None
		# set the new multi-index
		df = df.set_index(['cq_paradigm', 'scenario', 'year', 'date'])
		order = ['error', 'err_type', 'file', 'last_date_ord', 'last_date_dt']
		df = df[order]
	else: df = None

	return df

def load_debug_parameters():
	parser = argparse.ArgumentParser(description='command-line arguments for executing the feeDebug.py script.')
	parser.add_argument('--conf', type=str, help='The path to the debug parameters JSON file.')
	args = parser.parse_args()
	with open(args.conf, 'r') as file:
		return json.load(file)

def setup_logging(log_name):
	logging.basicConfig(
		level=logging.INFO,           # Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
		format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
	)
	logger = logging.getLogger()
	logger.addHandler(logging.FileHandler(f'{log_name}.log', 'w'))
	return logger

def main():
	# define scenarios and dir_years as global variables, because we are modifying them in main()
	# NOTE: if you're just accessing the value of variables defined outsidd the local scope (function)
	# you don't need to define them as globals; only if you plan on modifying them do you need to declar them as global
	# global scenarios
	# global dir_years

	# determine the error message to use based on fee version
	# NOTE: prior to v3, there was no universal failed run message, but the most common failure was an AEM3D one
	# common_error_msg = 'Exception: AEM3D failure.'
	common_error_msg = 'ERROR: RUN FAILED'

	##### DEBUG PARAMETERS #####
	debug_params = load_debug_parameters()
	root_dir = debug_params["root_dir"]
	cq_paradigms = debug_params["cq_paradigms"]
	dir_years = debug_params["dir_years"]
	scenarios = debug_params["scenarios"]
	launch_start = debug_params["launch_start"]
	launch_end = debug_params["launch_end"]
	file_path = debug_params["file_path"]
	last_lines_to_read = debug_params["last_lines_to_read"]
	log_name = debug_params["log_name"]
	
	logger = setup_logging(log_name)
	print = logger.info

	# Print the debug parameters
	print(f"##### DEBUG PARAMETERS #####")
	print(f"ROOT DIR: {root_dir}")
	print(f"CQ PARADIGMs: {cq_paradigms}")
	print(f"SCENARIOS: {scenarios}")
	print(f"YEARS: {dir_years}")
	print(f"FIRST DATE: {launch_start}")
	print(f"LAST DATE: {launch_end}\n")

	fail_dict = {}

	for cq in cq_paradigms:
		cq_path = os.path.join(root_dir, cq)
		print(f"CURRENT CQ PARADIGM: {cq}")
		fail_dict[cq] = {}
		if scenarios == 'all':
			scenarios = os.listdir(cq_path)
		for scen in scenarios:
			print(f"\t SCENARIO: {scen}")
			scen_path = os.path.join(cq_path, scen)
			if dir_years == 'all':
				with os.scandir(scen_path) as entries:
					dir_years_list = [entry.name for entry in entries if entry.is_dir()]
				dir_years_list = [str(year) for year in sorted([int(yr) for yr in dir_years_list], reverse=True)]
			else: dir_years_list = dir_years
			fail_dict[cq][scen] = {}
			for dr_yr in dir_years_list:
				print(f"\t\t YEAR: {dr_yr}")
				fail_dict[cq][scen][dr_yr] = {}
				launch_start_date = dt.datetime.strptime(launch_start, "%m%d").replace(year=int(dr_yr))
				launch_end_date = dt.datetime.strptime(launch_end, "%m%d").replace(year=int(dr_yr))
				# print(launch_start_date)
				date_range = [launch_start_date + dt.timedelta(days=x) for x in range((launch_end_date - launch_start_date).days + 1)]
				for date in date_range:
					# fail_dict[cq][scen][dr_yr][date.strftime("%Y%m%d")] = None
					run_dir = os.path.join(root_dir, cq, scen, dr_yr, date.strftime("%Y%m%d"))

					# skip to try the next run dir if the current one does not exist
					if not os.path.exists(run_dir): continue

					# list the directory and filter to get just the outfile
					outfile = [f for f in os.listdir(run_dir) if f.endswith(".out")]
					# the list above should have a length of exactly 1 (there should be one outfile)
					# raise a warning if for some reason ther is more than 1 outfile
					if len(outfile) != 1:
						warnings.warn(f"Multiple outfiles detected: {outfile}")
					outfile = outfile[0]

					# create full file path
					full_out_path = os.path.join(run_dir, outfile)
					# print(full_out_path)
					# open file and read the lines at the end of the file
					with open(full_out_path, 'r') as file:
						lastlines = list(deque(file, maxlen=last_lines_to_read))

					if lastlines[-1].rstrip('\n') == common_error_msg:
						print(f"\t\t\t DATE: {date.strftime('%m/%d/%Y')}")
						print(f"\t\t\t PATH: {full_out_path}")
						error_report = {'error':'run_failure', 'err_type':'unkown', 'file':outfile, 'last_date_ord':-1, 'last_date_dt':-1}
						# print("\nPASSING TO PARSE_FAILED_RUNS\n")
						error_report = parse_failed_run(lastlines, error_report)
						# print(error_report)
						fail_dict[cq][scen][dr_yr][date.strftime('%Y%m%d')] = error_report
					# # optional message
					# else: print("\t\t\t\t No run failure detected")

	# generate the concise report table
	report = report_df(fail_dict)
	
	# don't try to write df if it is None (no data to report)
	if report is not None:
		report.to_csv(f'{log_name}.csv')

	print("DEBUG SCRIPT COMPLETE")

if __name__ == '__main__':
	main()