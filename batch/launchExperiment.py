import json
import os
import subprocess
import datetime as dt
import time
from string import Template
import sys

'''
Script to Launch Model Runs for a given hindcast scenario
'''
# directory in ehcih this .py script is located
orig = os.getcwd()

# read in the experiment-spoecific config file from the command line
experiment_conf_file = sys.argv[1]

# read in the experiment-specific config file
with open(experiment_conf_file) as experiment_config:
	experiment_launch_params = json.load(experiment_config)

# load in the default run config file
with open("/gpfs1/home/n/b/nbeckage/ciroh/forecast-workflow/default_settings.json") as default_json_file:
	config = json.load(default_json_file)

# load in the default job script template
with open('/netfiles/ciroh/7dayHABsHindcast/submitTEMPLATE.sh') as f:
	src = Template(f.read())

root_dir = '/netfiles/ciroh/7dayHABsHindcast/'
scenario_dir = os.path.join(root_dir, f"{experiment_launch_params['scenario']}/")

years = list(reversed([str(y) for y in range(int(experiment_launch_params['end_year']), int(experiment_launch_params['start_year'])+1)]))

for year in years:
	start_dt = dt.datetime.strptime(year+experiment_launch_params['start_date'], '%Y%m%d')
	end_dt = dt.datetime.strptime(year+experiment_launch_params['end_date'], '%Y%m%d')
	delta = end_dt - start_dt
	scenario_dir_year = os.path.join(scenario_dir, f'{year}/')
	dates = [start_dt + dt.timedelta(days=d) for d in range(delta.days+1)]
	for date in dates:
		scenario_dir_date = os.path.join(scenario_dir_year,date.strftime('%Y%m%d'))
		# if the run exists, skip it
		if os.path.exists(scenario_dir_date):
			continue
		os.makedirs(scenario_dir_date)
		os.chdir(scenario_dir_date)

		config['spinup_date'] = dt.datetime(date.year, 1, 1).strftime('%Y%m%d')
		config['forecast_start'] = date.strftime('%Y%m%d')
		config['forecast_end'] = (date + dt.timedelta(days=7)).strftime('%Y%m%d')
		config['root_dir'] = root_dir
		config['weather_dataset_observed'] = experiment_launch_params['weather_dataset_observed']
		config['weather_dataset_forecast'] = experiment_launch_params['weather_dataset_forecast']
		config['hydrology_dataset_observed'] = experiment_launch_params['hydrology_dataset_observed']
		config['hydrology_dataset_forecast'] = experiment_launch_params['hydrology_dataset_forecast']
		config['nwm_forecast_member'] = experiment_launch_params['nwm_forecast_member']

		# wrtie the run-specific config file
		with open('configuration.json', 'w') as config_file:
			json.dump(config, config_file, indent=2)
			config_file.write('\n')

		# define the job params for the run
		job_params = {'job_name':f'{experiment_launch_params["job_name_prefix"]}{date.strftime("%Y%m%d")}',
			  		  'run_dir':scenario_dir_date}
		job_script = src.substitute(job_params)

		# write the run-specifc submit script
		with open('submit.sh', 'w') as submit_script:
			submit_script.write(job_script)

		print(scenario_dir_date)
		# Uncomment to submit
		subprocess.run(['sbatch ' 'submit.sh'], shell=True)
		# wait so as to not overwhelm the VACC with GFS / NWM file reads for AEM3D_prep_IAM
		time.sleep(30)
os.chdir(orig)
