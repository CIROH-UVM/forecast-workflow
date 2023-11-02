from default_config import defaults
import argparse
import configparser
from datetime import date, timedelta
import os

SETTINGS_KEYS = ['forecast_start',
				'forecast_end',
				'spinup_start',
				'blending_variable',
				'blending_ratio',
				'weather_dataset_observed',
				'weather_dataset_forecast',
				'hydrology_dataset_observed',
				'hydrology_dataset_forecast']

# setting up command-line argument parser
def get_args():
	parser = argparse.ArgumentParser(description="command-line arguments for running the AEM3D-based HABs Forecast",
									epilog="more details and documentation to come soon")
	# adding an optional argument for a config file
	parser.add_argument('--config', type=str, help="name of the configuration file containing model run settings")
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

def main():
	args = get_args()
	if args.config is None:
		print(f'config args:{args.config}')
	else: print('no config args passed')

if __name__ == '__main__':
	main()
