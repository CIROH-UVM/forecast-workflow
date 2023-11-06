# from models.aem3d.default_settings import defaults
import argparse
from datetime import datetime, timedelta
import json
import os

##### Hardcoded params ######
# absolute filepath for the default_settings.json file
default_settings_fpath = "/data/users/n/b/nbeckage/forecast-workflow/default_settings.json"

def check_keys(settings_dict):
	SETTINGS_KEYS = get_settings_keys()
	for key in settings_dict:
		if key not in SETTINGS_KEYS:
			raise KeyError(f'invalid settings key: {key}')
		
def check_settings(sd):
	if not isinstance(sd['forecast_start'], datetime):
		sd['forecast_start'] = datetime.strptime(sd['forecast_start'], '%m/%d/%Y')
	if not isinstance(sd['forecast_end'], datetime):
		sd['forecast_end'] = datetime.strptime(sd['forecast_end'], '%m/%d/%Y')
	sd['spinup_date'] = datetime.strptime(sd['spinup_date'], '%m/%d/%Y')


# setting up command-line argument parser
def get_args():
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

	args = parser.parse_args()

	return args

def get_settings_keys():
	return list(load_json(default_settings_fpath).keys())

def load_json(fpath):
	with open(fpath) as file:
		data = json.load(file)
	return data

def load_config(config_fpath):
	config_settings = load_json(config_fpath)
	check_keys(config_settings)
	return config_settings

# reads default settings json file and returns dictionary of default settings
def load_defaults():
	today = datetime.today()
	defaults = load_json(default_settings_fpath)
	# manually set the current date and 7 days from today - can't do this programmatically in json
	defaults['forecast_start'] = today
	defaults['forecast_end'] = today + timedelta(days=7)
	return defaults

def process_args(args):
	SETTINGS_KEYS = get_settings_keys()
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

def main():
	config_fpath = "/data/users/n/b/nbeckage/forecast-workflow/custom_settings.json"
	# establish default settings first
	settings = load_defaults()
	# load command-line arguments
	args = get_args()
	# load settings from configuration file if given
	if args.conf is not None:
		custom_settings = load_config(args.conf)
		settings.update(custom_settings)
	cmd_args = process_args(args)
	settings.update(cmd_args)
	check_settings(settings)
	print(settings)


if __name__ == '__main__':
	main()

"""

def get_config():
	config_fpath = '/data/users/n/b/nbeckage/forecast-workflow/test_settings.cfg'
	config_file = configparser.ConfigParser()
	# read the configuration file
	config_file.read(config_fpath)
	print(config_file.sections())
	# the settings section should be the first and only section
	settings_section = config_file.sections()[0]
	print(settings_section)
	# pull the settings section and make it a dictionary
	config_settings = dict(config_file[settings_section])
	print(config_settings)

"""