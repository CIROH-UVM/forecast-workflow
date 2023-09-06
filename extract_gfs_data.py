from gfs_download_fcns import *

# locations to pull forcast data for
loc_list = [(45.0, -73.25), (44.75, -73.25)]

# initialize empty dictionary which will hold dataframes for each lat/long pair in loc_list
location_dataframes = {}

##### Varaibles Needed:
# ENTIRE ATMOSPHERE, instant + avg
# - Total Cloud Cover [%]
# SURFACE, avg
# - Downward Short-Wave Radiation Flux [W/m^2]
# 2 M ABOVE GROUND, instant
# - Temperature [K]
# - Relative Humidity [%]
# 10 M ABOVE GROUND, instant
# - U-Component of Wind [m/s]
# - V-Component of Wind [m/s]

##### Variables sucessfully downloaded:
# SURFACE, instant
# - "Temperature, K"
# - "Precipitation rate, kg m**-2 s**-1"e
# - "Snow depth, m"
# - "Categorical snow, (Code table 4.222)"


### RUN PARAMS
# home directory in which to run the program
run_dir = '/netfiles/ciroh/nbeckage/gfs_data/stInstant_lvSurface/'
if not os.path.exists(run_dir): os.makedirs(run_dir)
os.chdir(run_dir)

# dates to get
dates = generate_date_strings('20230828',1)

# forecast daysa ahead to get
hours = generate_hours_list(0)

# define stepType
stepType = 'instant'

# cfgrib.dataset.DatasetBuildError: multiple values for unique key, try re-open the file with one of:
#     filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'}
#     filter_by_keys={'stepType': 'avg', 'typeOfLevel': 'surface'}
#     filter_by_keys={'stepType': 'accum', 'typeOfLevel': 'surface'}

# define typeOfLevel
# typeOfLevel = 'surface'
typeOfLevel = 'heightAboveGround'

# cfgrib.dataset.DatasetBuildError: multiple values for unique key, try re-open the file with one of:
#     filter_by_keys={'typeOfLevel': 'meanSea'}
#     filter_by_keys={'typeOfLevel': 'hybrid'}
#     filter_by_keys={'typeOfLevel': 'atmosphere'}
#     filter_by_keys={'typeOfLevel': 'surface'}
#     filter_by_keys={'typeOfLevel': 'planetaryBoundaryLayer'}
#     filter_by_keys={'typeOfLevel': 'isobaricInPa'}
#     filter_by_keys={'typeOfLevel': 'isobaricInhPa'}
#     filter_by_keys={'typeOfLevel': 'heightAboveGround'}
#     filter_by_keys={'typeOfLevel': 'depthBelowLandLayer'}
#     filter_by_keys={'typeOfLevel': 'heightAboveSea'}
#     filter_by_keys={'typeOfLevel': 'atmosphereSingleLayer'}
#     filter_by_keys={'typeOfLevel': 'lowCloudLayer'}
#     filter_by_keys={'typeOfLevel': 'middleCloudLayer'}
#     filter_by_keys={'typeOfLevel': 'highCloudLayer'}
#     filter_by_keys={'typeOfLevel': 'cloudCeiling'}
#     filter_by_keys={'typeOfLevel': 'heightAboveGroundLayer'}
#     filter_by_keys={'typeOfLevel': 'tropopause'}
#     filter_by_keys={'typeOfLevel': 'maxWind'}
#     filter_by_keys={'typeOfLevel': 'isothermZero'}
#     filter_by_keys={'typeOfLevel': 'highestTroposphericFreezing'}
#     filter_by_keys={'typeOfLevel': 'pressureFromGroundLayer'}
#     filter_by_keys={'typeOfLevel': 'sigmaLayer'}
#     filter_by_keys={'typeOfLevel': 'sigma'}
#     filter_by_keys={'typeOfLevel': 'potentialVorticity'}

### Past this line back into aggregate_df()
# , backend_kwargs={'filter_by_keys': {'stepType': stepType,'typeOfLevel': typeOfLevel}}


# 1. pull the data at the specified data/time ranges
# pull_gribs(dates, hours)

# 2. open each file, read it as a dataframe, process it, and it it to location_dataframes
aggregate_df(dates, hours, loc_list, location_dataframes, stepType, typeOfLevel)

# write out csv's
# dict_to_csv(loc_list, location_dataframes)