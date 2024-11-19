from data.utils import parse_to_datetime
import datetime as dt
import fsspec
import xarray as xr

def get_data(start_date, end_date, locations, variables = {'streamflow':'streamflow'}):
	'''
	Get channel output (CHRTOUT files) from NWMv3.0 retrospective forecasts hosted on Amazon Web Services (https://registry.opendata.aws/nwm-archive/).
	NWM retrospective forecasts only produce one long-term retrospective run per NWM verison. Source: https://onlinelibrary.wiley.com/doi/10.1111/1752-1688.13184

	Args:
	-- start_date (str, date, or datetime) [req]: the start date for the forecat data grab
	-- end_date (str, date, or datetime) [req]: the end date for the forecast data grab
	-- locations (dict) [req]: a dictionary (reachName:reachID) of NWM reaches to get data for.
	-- variables (dict) [req]: a dictionary of variables to download. Keys should be user-defined var names, value should be dataset-specific var names (qBtmVertRunoff, qbucket, etc).
		Available variables in CHRTOUT files:
		 - qBtmVertRunoff: Runoff from bottom of soil to bucket (m3)
		 - qBucket: Flux from gw bucket (m3 s-1)
		 - qSfcLatRunoff: Runoff from terrain routing (m3 s-1)
		 - q_lateral: Runoff into channel reach (m3 s-1)
		 - streamflow: River Flow (m3 s-1)
		 - velocity: River Velocity (m s-1)
		
	Returns:
	NWMv3.0 retrospective forecast CHRTOUT data timeseries for the given locations in a nested dict format where 1st-level keys are user-provided location names and 2nd-level keys
	are variables names and values are the respective data in a Pandas Series object. 
	'''
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	print(f"BEGIN NWMv3.0 RETRO FORECAST CHRTOUT GET FROM {start_date.strftime('%Y-%m-%d %H:%M:%S')} TO {end_date.strftime('%Y-%m-%d %H:%M:%S')}")

	# NOTE: there are 24 timesteps for each day, 00-23
	# define the amazon web bucket you want to use
	bucket = 's3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/chrtout.zarr/'

	# open entire chrtout dataset
	ds = xr.open_zarr(fsspec.get_mapper(bucket, anon=True), consolidated=True)

	# can filter by time first, since we are getting the same time slice for each location and var
	ds = ds.sel(time=slice(start_date.strftime('%Y-%m-%dT%H:%M:%S'), end_date.strftime('%Y-%m-%dT%H:%M:%S')))

	# get list of reach ids for 
	reach_ids = [int(id) for id in locations.values()]

	### dataset coordinates:
	# -- gage_id - NHD Gage Event ID from SOURCE_FEA field in Gages feature class
	# -- feature_id - NHDPlusv2 ComIDs within CONUS, arbitrary Reach IDs outside of CONUS# filter by reah ID

	# NOTE: dataset also has latitude, longitude, elevation, and gage_id coordinates. To select by these coords,
	#		a list of coord values (ids, etc) is needed. Can select by mutiple coords in one .sel() statement.
	# although if you did add more coords, you would need a new dict like 'locations' for events, etc
	# and then of course implement in downstream code
	# IF you did want this module to get data by other coords (events, (lat,lon) pairs, etc), one way would be
	# to add new dicts as mentioned above, then return 3-layed dict (add 'gages', 'events', 'locations', etc as new top layer)
	ds = ds.sel(feature_id = reach_ids)

	# extract units from datset
	units = {varname : ds[varname].units for varname in variables.values()}

	# get a list of dataset var names
	vars_to_get = list(variables.values())

	print(f"Getting the following variables: {vars_to_get}")
	print(f"For the following reaches: {locations}")
	# now convert filtered dataset to flat dataframe
	df = ds[vars_to_get].to_dataframe().reset_index().loc[:, ['time','feature_id']+vars_to_get].set_index('time')

	# create dict for series data
	nwmretro_q = {locname : {} for locname in locations.keys()}

	# group df by reach_id
	reach_groups = df.groupby('feature_id')

	inverted_locations = {id : locname for locname, id in locations.items()}

	# first, iterate through reach groups (locations)
	for reach_id, reach_df in reach_groups:
		# get the user-name for the specific reach
		reach_name = inverted_locations[str(reach_id)]
		# then, iterate through the variables list for each reach
		for varname, var in variables.items():
			# get variable series from df
			reach_series = reach_df[var]
			# localize series indiex to UTC time
			reach_series.index = reach_series.index.tz_localize(dt.timezone.utc)
			# rename series to include units and drop last row to make series end_date exclusive
			reach_series = reach_series.rename(f"{var} ({units[var]})").iloc[:-1]
			nwmretro_q[reach_name].update({varname:reach_series})


	print("NWMv3.0 RETRO FORECAST CHRTOUT GET COMPLETE")
	return nwmretro_q