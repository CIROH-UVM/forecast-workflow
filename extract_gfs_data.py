from gfs_download_fcns import *

# locations to pull forcast data for
loc_dict = {'401':(45.0, -73.25),
            '402':(44.75, -73.25),
            '403':(44.75, -73.25)}

##### Varaibles Needed:
### listed in gfs grib2 'longnames' format
# ENTIRE ATMOSPHERE, instant
# - 'Total Cloud Cover, %'
# SURFACE, avg
# - 'Downward short-wave radiation flux, W m**-2'
# 2 M ABOVE GROUND, instant
# - '2 metre temperature, K'
# - '2 metre relative humidity, %'
# 10 M ABOVE GROUND, instant
# - '10 metre U wind component, m s**-1'
# - '10 metre V wind component, m s**-1'
# SURFACE, instant
# - 'Precipitation rate, kg m**-2 s**-1'
# - 'Percent frozen precipitation, %'

### to do
# - Convert step to date and time
# - grab forecasts t00 instead of t12
# - time 000 for short wave
# - convert execute to sh 


### RUN PARAMS
# home directory in which to run the program and store data
run_dir = '/netfiles/ciroh/nbeckage/gfs_data/'
if not os.path.exists(run_dir): os.makedirs(run_dir)
os.chdir(run_dir)

# determine dates for which to get forecast data
# date_str = datetime.today().strftime("%Y%m%d")
date_str = '20230906'
dates = generate_date_strings(date_str,2)

# determine how many days worth of forecasts to get
hours = generate_hours_list(7)[0:3]

##### Future reference, pasted stepType and typeOfLevel params
### - I couldn't find any way to see all of these params besides in these error messages
# cfgrib.dataset.DatasetBuildError: multiple values for unique key, try re-open the file with one of:
#     filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'}
#     filter_by_keys={'stepType': 'avg', 'typeOfLevel': 'surface'}
#     filter_by_keys={'stepType': 'accum', 'typeOfLevel': 'surface'}
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


##### RUNNING THE SCRIPT
# uncomment each function as needed

# 1. pull the data at the specified data/time ranges
# pull_gribs(dates, hours)

# 2. open each grib2 file, read the dataframes, merge them, and then return a dict of dfs for each station
master_dict = aggregate_df_dict(dates, hours, loc_dict)

# 3. write out csv's for each dataframe in master_dict
dict_to_csv(loc_dict, master_dict)