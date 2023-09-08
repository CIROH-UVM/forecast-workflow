from datetime import datetime, timedelta
import os
import pandas as pd
import subprocess as sp
import xarray as xr

### global vars that wouldn't change based on user; may change depending on what forecast data is being pulled
# root url where the past 10 day forecasts subfolders are located:
gfs_root = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/'
fc_time = '/12/atmos/'
fc_file = 'gfs.t12z.pgrb2.0p25.f'
# where the raw grib2 files will be stored
fc_data_dir = 'raw_fc_data/'
# where any csv files created for each location/station will be stored
location_data_dir = "loc_data/"


def aggregate_df_dict(dates = [], hours = [], loc_dict = {}):
    # initialize a dictionary to store a dataframe for each station in
    data_dict = {}
    
    if not os.path.exists(fc_data_dir): os.makedirs(fc_data_dir)
    os.chdir(fc_data_dir)
    
    for d in dates:
        date_dir = 'gfs.'+d+fc_time
        os.chdir(date_dir)
        for h in hours:
            hour_file = fc_file+h
            df_list = get_subgrib_df_list(hour_file)
            vars_df = merge_subgrib_dfs(df_list)
            extract_locs(vars_df, loc_dict, data_dict)
        os.chdir('../../../')
    os.chdir('../')
    return data_dict

def dict_to_csv(loc_dict = {}, location_dataframes = {}):
    # make directory for location data
    if not os.path.exists(location_data_dir): os.makedirs(location_data_dir)
    os.chdir(location_data_dir)
    for station in loc_dict:
        location = loc_dict[station]
        filename = f"{station}_{location[0]}_{location[1]}.csv"
        location_dataframes[station].to_csv(filename)
    return

def execute(cmd):
    popen = sp.Popen(cmd, stdout=sp.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close
    return_code = popen.wait()
    if return_code:
        raise sp.CalledProcessError(return_code, cmd)
    return

### Given the transformed dataframe and a list of lat/long tuples to extract, returns new df containing just the rows for each lat/long pair
# -- df : the grib2 df post-long-transform
# -- loc_dict : dict of of station names and corresponding lat/long tuples; will pull said coords from the datafram
# given the transformed dataframe, adds the new rows to the respective composite dataframes in the df dictionary
def extract_locs(df, loc_dict = {}, location_dataframes = {}):
    for station in loc_dict:
        coords = loc_dict[station]
        extracted_df = pd.DataFrame(df.loc[coords]).T.set_index('step')
        # if the dictionary does not already have a key for each location, then initialize that key; should only be true when extracting the 1st df
        if len(location_dataframes) != len(loc_dict):
            location_dataframes[station] = extracted_df
        else:
            location_dataframes[station] = pd.concat([location_dataframes[station], extracted_df])
    return

def extract_subgrib(fname = '', args = {}):
    ds = xr.open_dataset(fname, engine="cfgrib", backend_kwargs={'filter_by_keys': args})
    longnames = ['time','step',args['typeOfLevel'], 'valid_time']
    for v in ds:
        longnames.append("{}, {}".format(ds[v].attrs["long_name"], ds[v].attrs["units"]))
    df = ds.to_dataframe()
    df.columns=longnames
    remap_longs(df)
    return df

### Creates a list of forecast dates to download
# -- start_date : first date of forecast data to download
# -- num_dates : the number of dates ahead or behind the start date you want to download
# -- cast : str switch that determines whether to grab n dates ahead ("fore") or behind ("hind")
def generate_date_strings(start_date, num_dates, cast="fore"):
    date_strings = []
    current_date = datetime.strptime(start_date, "%Y%m%d")

    for _ in range(num_dates):
        date_strings.append(current_date.strftime("%Y%m%d"))
        if cast == "hind":
            current_date -= timedelta(days=1)
        else:
            current_date += timedelta(days=1)

    return date_strings

def get_subgrib_df_list(fname = ''):
    args_list = [{'typeOfLevel':'atmosphere', 'stepType':'instant'},
             {'typeOfLevel':'heightAboveGround', 'topLevel':10},
             {'typeOfLevel':'heightAboveGround', 'topLevel':2},
             {'typeOfLevel':'surface', 'stepType':'avg'},
             {'typeOfLevel':'surface', 'stepType':'instant'}]
    # temporary list to store un-extracted dataframes
    df_list = []
    for args in args_list:
        subgrib_df = extract_subgrib(fname, args)
        df_list.append(subgrib_df)
    return df_list


### Creates a list of forecast hours to be downloaded
# -- num_days : how many days out of forecast data you want to download (i.e. 5, 7, 10, etc)
# - may want to change this fcn so that it can generate a more precise hour list (i.e. specify how many hours out you want fc data)
def generate_hours_list(num_days):
    hours_list = []
    hour = 0
    hours_list.append(f"{hour:03}")
    for day in range(1, num_days + 1):
        if day <= 5:
            for h in range(24):
                hour += 1
                hours_list.append(f"{hour:03}")
        else:
            for h in range(8):
                hour += 3
                hours_list.append(f"{hour:03}")
    return hours_list

def merge_subgrib_dfs(df_list = []):
    # defining a list of the data that we want to extract for the forecast model
    # - to grab more data, simply add them to the list
    # - note that you may have to add an additional args dict to args_list (See get_subgrib_df_list()) in order to get the data you want
    cols_to_keep = ['2 metre temperature, K',
                'Total Cloud Cover, %',
                'Downward short-wave radiation flux, W m**-2',
                '10 metre U wind component, m s**-1',
                '10 metre V wind component, m s**-1',
                '2 metre relative humidity, %',
                'Precipitation rate, kg m**-2 s**-1',
                'Percent frozen precipitation, %']
    # variable names required for model input file
    var_names = ['time', 'step', 'T2', 'TCDC', 'SWDOWN' , 'U10', 'V10', 'RH2', 'RAIN', 'CPOFP']
    # drop the surface average precip rate; it has the same col name as the surface instant precip col, which is problematic b/c we just need the latter
    # ! easily breakable line ! 
    df_list[3] = df_list[3].drop('Precipitation rate, kg m**-2 s**-1', axis=1)
    # combine all dfs in list
    all_cols_df = pd.concat(df_list, axis=1)
    # keep time and step indices at beginning of df
    ts_indices = all_cols_df.iloc[:,[0,1]]
    # get the vars we need, in order
    vars_df = all_cols_df[cols_to_keep]
    merged_df = pd.concat([ts_indices, vars_df], axis=1)
    merged_df.columns = var_names
    return merged_df

### Downloads gfs data into directories that mirrors the GFS directory structure
# -- dates : list of forecast dates to download
# -- hours : list of forecast hours to download
def pull_gribs(dates = generate_date_strings('20230828',1), hours = generate_hours_list(0)):
    # make the subdirectory for the raw data
    if not os.path.exists(fc_data_dir): os.makedirs(fc_data_dir)
    os.chdir(fc_data_dir)
    
    for d in dates:
        date_dir = 'gfs.'+d+fc_time
        if not os.path.exists(date_dir): os.makedirs(date_dir)
        os.chdir(date_dir)
        date_url = gfs_root+date_dir+fc_file
        for h in hours:
            hour_url = date_url + h
            if not os.path.exists(fc_file+h):
                print("Downloading file from URL:",hour_url)
                for path in execute(['curl', '--connect-timeout','120','-O', hour_url]):
                    print(path, end="")
                print("Download Complete:",date_dir+fc_file+h,"\n")
                # pull_call = sp.run(['curl', '-O', hour_url], capture_output = True, check = True)
        os.chdir('../../../')
    os.chdir('../')
    return


### In-place function that transforms the longitude indices from 0-360 t0 -180-180
def remap_longs(df):
    map_function = lambda lon: (lon - 360) if (lon > 180) else lon
    longitudes = df.index.levels[1]
    remapped_longitudes = longitudes.map(map_function)
    df.index = df.index.set_levels(remapped_longitudes, level="longitude")
    return
