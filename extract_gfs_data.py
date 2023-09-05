from gfs_download_fcns import *

# locations to pull forcast data for
loc_list = [(45.0, -73.25), (44.75, -73.25)]

# initialize empty dictionary which will hold dataframes for each lat/long pair in loc_list
location_dataframes = {}

# home directory in which to run the program
home_dir = '/netfiles/ciroh/nbeckage/gfs_data'
os.chdir(home_dir)

dates = generate_date_strings('20230827',10)
hours = generate_hours_list(7)

# 1. pull the data at the specified data/time ranges
pull_gribs(dates, hours)

# 2. open each file, read it as a dataframe, process it, and it it to location_dataframes
aggregate_df(dates, hours, loc_list, location_dataframes)

# write out csv's
dict_to_csv(loc_list, location_dataframes)