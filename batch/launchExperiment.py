import json
import os
import subprocess
import datetime as dt
import time
from string import Template
import sys
import models.aem3d.get_args as get_args
import random

'''
Script to Launch Model Runs for a given hindcast scenario
'''
# directory in which this .py script is executed from - should be a scenario dir according to FEE protocol
orig = os.getcwd()

# read in the experiment-spoecific config file from the command line
experiment_conf_file = sys.argv[1]

# read in the experiment-specific config file
with open(experiment_conf_file) as experiment_config:
	experiment_launch_params = json.load(experiment_config)

# parse start and end year args from command line if passed
try:
	start_year = int(sys.argv[2])
except IndexError:
	start_year = int(experiment_launch_params['start_year'])
try:
	end_year = int(sys.argv[3])
except IndexError:
	end_year = int(experiment_launch_params['end_year'])
# parse test arg from command line if passed
# should be an integer specifying the number of runs to randomly launch for a given year
try:
	test = int(sys.argv[4])
except IndexError:
	test = False

# load in the default run config file
config = get_args.load_defaults()

# load in the default job script template
with open('/users/n/b/nbeckage/ciroh/forecast-workflow/batch/submitTEMPLATE.sh') as f:
	src = Template(f.read())

### 6/11/24 - I think this is outdated - orig will be the scenario dir (same dir as experiemnt_config)
# root_dir = '/netfiles/ciroh/7dayHABsHindcast/'
# scenario_dir = os.path.join(root_dir, f"{experiment_launch_params['scenario']}/")

if start_year > end_year:
	years = list(reversed([str(y) for y in range(end_year, start_year+1)]))
else: years = [str(y) for y in range(start_year, end_year+1)]

# print(years)
# sys.exit()

for year in years:
	start_dt = dt.datetime.strptime(year+experiment_launch_params['start_date'], '%Y%m%d%H')
	end_dt = dt.datetime.strptime(year+experiment_launch_params['end_date'], '%Y%m%d')
	delta = end_dt - start_dt
	scenario_dir_year = os.path.join(orig, f'{year}/')
	dates = [start_dt + dt.timedelta(days=d) for d in range(delta.days+1)]
	# if testing, takes a random sample of the dates
	if test:
		dates = random.sample(dates, test)
	spinup_month_and_day = dt.datetime.strptime(experiment_launch_params['spinup_date'], "%m%d")
	# get the length of the forecast in days (7 days, 30, etc)
	forecast_days = int(experiment_launch_params['forecast_days'])
	for date in dates:
		scenario_dir_date = os.path.join(scenario_dir_year,date.strftime('%Y%m%d.t%Hz'))
		# if the run exists, skip it
		if os.path.exists(scenario_dir_date):
			continue
		os.makedirs(scenario_dir_date)
		os.chdir(scenario_dir_date)
		
		config['spinup_date'] = dt.datetime(date.year, spinup_month_and_day.month, spinup_month_and_day.day).strftime('%Y%m%d')
		config['forecast_start'] = date.strftime('%Y%m%d%H')
		config['forecast_end'] = (date + dt.timedelta(days=forecast_days)).strftime('%Y%m%d')
		config['data_dir'] = experiment_launch_params['data_dir']
		config['weather_dataset_spinup'] = experiment_launch_params['weather_dataset_spinup']
		config['weather_dataset_forecast'] = experiment_launch_params['weather_dataset_forecast']
		config['hydrology_dataset_spinup'] = experiment_launch_params['hydrology_dataset_spinup']
		config['hydrology_dataset_forecast'] = experiment_launch_params['hydrology_dataset_forecast']
		config['nwm_forecast_member'] = experiment_launch_params['nwm_forecast_member']
		config['aem3d_command_path'] = experiment_launch_params['aem3d_command_path']
		config['cqVersion'] = experiment_launch_params['cqVersion']

		# wrtie the run-specific config file
		with open('configuration.json', 'w') as config_file:
			json.dump(config, config_file, indent=2)
			config_file.write('\n')

		# define the job params for the run
		job_params = {'job_name':f'{experiment_launch_params["job_name_prefix"]}{date.strftime("%Y%m%d.t%Hz")}',
			  		  'run_dir':scenario_dir_date}
		job_script = src.substitute(job_params)

		# write the run-specifc submit script
		with open('submit.sh', 'w') as submit_script:
			submit_script.write(job_script)

		print(f"Launching run in: {scenario_dir_date}")
		sys.stdout.flush()
		# Uncomment to submit
		subprocess.run(['sbatch ' 'submit.sh'], shell=True)
		# wait so as to not overwhelm the VACC with GFS / NWM file reads for AEM3D_prep_IAM
		time.sleep(45)
os.chdir(orig)

print("Experiment Launch Complete.")