from datetime import date, timedelta

# Default configuration for aem3d

today = date.today()

defaults = {
	'forecast_start' : today,
	'forecast_end' : today + timedelta(days=7),
	'spinup_date' : date(year=today.year, month=1, day=1),
	'blending_variable' : '',
	'blending_ratio' : 1.0,
	'weather_dataset_observed' : 'NOAA_LCD+FEMC_CR',
	'weather_dataset_forecast' : 'NOAA_GFS',
	'hydrology_dataset_observed' : 'USGS_IV',
	'hydrology_dataset_forecast' : 'NOAA_NWM_PROD'
}

# print(defaults)