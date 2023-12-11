import argparse
from datetime import date, datetime, timedelta
import json
import os

"""
Settings parser for AEM3D_prep_IAM.py. To get settings, simply import get_args(), define the path to your default
settings json file, and call get_args():

	from get_args import get_args

	fpath = 'path/to/default_settings.json'

	settings = get_args(fpath)

Now when you call AEM3D_prep_IAM.py, you can include args in the command-line call:

python -m models.aem3d.AEM3D_prep_IAM --conf config.json

To add new arguments:
	1. add arg and default value to default_settings.json
	2. 

**NOTE: blending variable and dataset values are currently not being validated in check_values(). Establish a list of 
valid dataset/blending var strings in order to do this.
"""

def check_keys(settings_dict):
	SETTINGS_KEYS = get_settings_keys()
	for key in settings_dict:
		if key not in SETTINGS_KEYS:
			raise KeyError(f'invalid settings key: {key}')
		
def check_values(settings_dict):
	# forecast start date must be datetime object
	if not isinstance(settings_dict['forecast_start'], datetime):
		settings_dict['forecast_start']  = make_datetime(settings_dict['forecast_start'])
	# forecast end date must be datetime object
	if not isinstance(settings_dict['forecast_end'], datetime):
		settings_dict['forecast_end'] = make_datetime(settings_dict['forecast_end'])
	settings_dict['spinup_date'] = make_datetime(settings_dict['spinup_date'])
	# # check blending variable to make sure it's a valid var name
	# if settings_dict['blending_variable'] not in valid_vars:
	# 	raise ValueError(f'{settings_dict["blending_variable"]} is an invalid variable name.')
	# blending ratio should be between 0 and 1
	if not 0 <= settings_dict['blending_ratio'] <= 1.0:
		raise ValueError(f"'{settings_dict['blending_ratio']}' is not a valid blending ratio; must be between 0 and 1.")
	# # datasets should be validated
	# if settings_dict['weather_dataset_observed'] not in valid_wdo:
	# 	raise ValueError(f"'{settings_dict['weather_dataset_observed']}' is not a valid observational weather dataset")
	# if settings_dict['weather_dataset_forecast'] not in valid_wdf:
	# 	raise ValueError(f"'{settings_dict['weather_dataset_forecast']}' is not a valid forecast weather dataset")
	# if settings_dict['hydrology_dataset_observed'] not in valid_hdo:
	# 	raise ValueError(f"'{settings_dict['hydrology_dataset_observed']}' is not a valid observational hydrology dataset")
	# if settings_dict['hydrology_dataset_forecast'] not in valid_hdf:
	# 	raise ValueError(f"'{settings_dict['hydrology_dataset_forecast']}' is not a valid forecast hydrology dataset")
	# don't need to check csv flag, as default is false, and if the flag is passed, it will change to true. Trying to pass a
	#  string or number for --csv will throw an error as an unrecognized bool
	# check to see if root dir exists
	# don't need to check --old_dirs flag, as default is false, and if the flag is passed, it will change to true. Trying to pass a
	#  string or number for --old_dirs will throw an error as an unrecognized bool
	
# setting up command-line argument parser
def get_cmdln_args():
	parser = argparse.ArgumentParser(description="command-line arguments for running the AEM3D-based HABs Forecast",
									epilog="more details and documentation to come soon")
	# adding an optional argument for a config file
	parser.add_argument('--conf', type=str, help="name of the configuration file containing model run settings")
	parser.add_argument('--fc_start', type=str, help='forecast start date')
	parser.add_argument('--fc_end', type=str, help='forecast end date')
	parser.add_argument('--spinup', type=str, help='model spinup date')
	parser.add_argument('--bl_var', type=str, help='name of the input variable to be blended')
	parser.add_argument('--bl_ratio', type=float, help='blending ratio')
	parser.add_argument('--wdo', type=str, help='observed weather dataset to use for model run')
	parser.add_argument('--wdf', type=str, help='forecasted weather dataset to use for model run')
	parser.add_argument('--hdo', type=str, help='observed hydrological dataset to use for model run')
	parser.add_argument('--hdf', type=str, help='forecasted hydrological dataset to use for model run')
	parser.add_argument('--csv', action='store_true', help="flag determining whether or not to use GFS/NWM CSV's. Default is False.")
	parser.add_argument('--root', type=str, help='root dir containing forecastScripts, forecastRuns, forecastData')
	parser.add_argument('--old_dirs', action='store_false', help="flag determining whether or not to use post-update dir structure. Default is True. IOW, pass this flag to use old dir structure")

	args = parser.parse_args()

	return args

def get_settings_keys(default_fpath):
	return list(load_json(default_fpath).keys())

def load_config(config_fpath):
	config_settings = load_json(config_fpath)
	check_keys(config_settings)
	return config_settings

# reads default settings json file and returns dictionary of default settings
def load_defaults(default_fpath):
	today = date.today()
	defaults = load_json(default_fpath)
	# manually set the current date and 7 days from today - can't do this programmatically in json
	defaults['forecast_start'] = datetime(today.year, today.month, today.day)
	defaults['forecast_end'] = datetime(today.year, today.month, today.day) + timedelta(days=7)
	return defaults

def load_json(fpath):
	with open(fpath) as file:
		data = json.load(file)
	return data

def make_datetime(date_str):
	try:
		datetime_obj = datetime.strptime(date_str, '%m/%d/%Y')
		return datetime_obj
	except ValueError:
		raise ValueError(f"'{date_str}' is not a valid date. Make sure dates are in the format: 'DD/MM/YYYY'")
	except Exception as e:
		raise e

def process_args(args, default_fpath):
	SETTINGS_KEYS = get_settings_keys(default_fpath)
	# convert args from a Namespace to a dict
	args = vars(args)
	# we don't need the configuration file setting anymore
	args.pop('conf')
	# create a mapping from json-format setting names to command-line setting names
	zipped_key_map = zip(SETTINGS_KEYS, args)
	key_map = list(zipped_key_map)
	# update the keys of the command-line arg dict to have the same key names as the json settings dict
	renamed_args = {setting[0]: args[setting[1]] for setting in key_map}
	# we only want to update settings with command-line args if they were passed
	passed_args = {key:renamed_args[key] for key in renamed_args if renamed_args[key] is not None}
	return passed_args


def get_args(default_fpath):
	"""
	Main method to read and parse settings for forecast-workflow.
	Hierachy of settings is as follows:
		- Command-line args override configuration file args
		- Configuration file args override default args
		- Default args are loaded by, of course, default

	Args:
	-- default_fpath (str) [required]: path of the default settings configuration file. Said file is stored in ../forecast-workflow/ currently.
	
	Returns:
	a dictionary of settings
	"""
	# establish default settings first
	settings = load_defaults(default_fpath)
	# load command-line arguments
	args = get_cmdln_args()
	# load settings from configuration file if config file is given
	if args.conf is not None:
		# loading config file settings
		custom_settings = load_config(args.conf)
		# updating settings dict with config file settings
		settings.update(custom_settings)
	# processing cmd line args to match cmd line arg names with json arg names
	cmd_args = process_args(args, default_fpath)
	# update settings with cmd line arguments
	settings.update(cmd_args)
	# check settings to ensure all values are valid
	check_values(settings)
	return settings

def main():
	# absolute filepath for the default_settings.json file
	default_settings_fpath = "/data/users/n/b/nbeckage/forecast-workflow/default_settings.json"
	print(f'SETTINGS:{get_args(default_settings_fpath)}')
	# settings = load_defaults(default_settings_fpath)
	# args = get_cmdln_args()
	# cmd_args = process_args(args, default_settings_fpath)
	# print(cmd_args)

if __name__ == '__main__':
	main()