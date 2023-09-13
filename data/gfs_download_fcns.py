from datetime import datetime, timedelta
import os
import pandas as pd
import subprocess as sp
import xarray as xr

### global vars that wouldn't change based on user; may change depending on what forecast data is being pulled
# root url where the past 10 day forecasts subfolders are located:
gfs_root = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/'
fc_time = '/00/atmos/'
fc_file = 'gfs.t00z.pgrb2.0p25.f'
# where the raw grib2 files will be stored
fc_data_dir = '/data/forecastData/gfs/raw_fc_data/'
# where any csv files created for each location/station will be stored
location_data_dir = "loc_data/"

# Define alias for aggreate_df_dict
def get_data(dates = [], hours = [], loc_dict = {}):
    return aggregate_df_dict(dates, hours, loc_dict)

def aggregate_df_dict(dates = [], hours = [], loc_dict = {}):
    # initialize a dictionary to store a dataframe for each station in
    data_dict = {}
    
    # Don't think we need this anymore now that downloads are seperate
    # if not os.path.exists(fc_data_dir): os.makedirs(fc_data_dir)
    origDir = os.getcwd()
    os.chdir(fc_data_dir)
    
    for d in dates:
        date_dir = 'gfs.'+d+fc_time
        os.chdir(date_dir)
        for h in hours:
            hour_file = fc_file+h
            print("aggregating subgribs for file",hour_file)
            df_list = get_subgrib_df_list(hour_file)
            # df_list = get_subgrib_df_list(os.path.join(fc_data_dir, date_dir, hour_file))
            print("succesfully aggregated subgribs")
            vars_df = merge_subgrib_dfs(df_list, hour_file)
            extract_locs(vars_df, loc_dict, data_dict)
            print("successfully merged and extracted loc data from subgribs\n")
        os.chdir(fc_data_dir)
    os.chdir(origDir)
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
        extracted_df = pd.DataFrame(df.loc[coords]).T.set_index('time')
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
# ---- should be a string in the format 'YYYYMMDD' or a datetime object
# -- num_dates : the number of dates ahead or behind the start date you want to download
# -- cast : str switch that determines whether to grab n dates ahead ("fore") or behind ("hind")
def generate_date_strings(start_date, num_dates = 1, cast="fore"):
    date_strings = []
    # if not a datetime object, convert start_date to one
    if not isinstance(start_date, datetime):
        # strptime(str, format) converts a string to a datetime object
        current_date = datetime.strptime(start_date, "%Y%m%d")
    else: current_date = start_date
    for _ in range(num_dates):
        date_strings.append(current_date.strftime("%Y%m%d"))
        if cast == "hind":
            current_date -= timedelta(days=1)
        else:
            current_date += timedelta(days=1)

    return date_strings

def get_subgrib_df_list(fname = ''):
    # drop surface avg args for the first (f000) file
    if fname[-3:] == "000":
        args_list = [{'typeOfLevel':'atmosphere', 'stepType':'instant'},
         {'typeOfLevel':'heightAboveGround', 'topLevel':10},
         {'typeOfLevel':'heightAboveGround', 'topLevel':2},
         {'typeOfLevel':'surface', 'stepType':'instant'}]
    else:
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
# -- num_hours: how many hours of forecast data you want: note that date goes up to 16 days out
# -- archive: boolean flag indicating if data is going to be pulled from archives; if true returns only step = 3 list
def generate_hours_list(num_hours = 168, archive = False):
    if archive:
        return [f"{hour:03}" for hour in range(0, num_hours+1, 3)]
    if not archive:
        if num_hours <= 120:
            return [f"{hour:03}" for hour in range(0,num_hours+1)]
        else: return [f"{hour:03}" for hour in range(0,120)] + [f"{hour:03}" for hour in range(120, num_hours+1, 3)]

def merge_subgrib_dfs(df_list = [], fname = ''):
    cols_to_keep = ['2 metre temperature, K',
                'Total Cloud Cover, %',
                'Downward short-wave radiation flux, W m**-2',
                '10 metre U wind component, m s**-1',
                '10 metre V wind component, m s**-1',
                '2 metre relative humidity, %',
                'Precipitation rate, kg m**-2 s**-1',
                'Percent frozen precipitation, %']
    var_names = ['time', 'T2', 'TCDC', 'SWDOWN' , 'U10', 'V10', 'RH2', 'RAIN', 'CPOFP']
    
    # if not the first file, drop surface avg precip rate
    if not fname[-3:] == "000":
        # drop the surface average precip rate; it has the same col name as the surface instant precip col, which is problematic b/c we just need the latter
        df_list[3] = df_list[3].drop('Precipitation rate, kg m**-2 s**-1', axis=1)
    # combine all dfs in list
    all_cols_df = pd.concat(df_list, axis=1)
    # keep valid_time indices at beginning of df
    ts_indices = all_cols_df.iloc[:,[3]]
    # create frozen precip if it does not exist
    if not 'Downward short-wave radiation flux, W m**-2' in all_cols_df:
        all_cols_df['Downward short-wave radiation flux, W m**-2'] = 0
    # get the vars we need, in order
    vars_df = all_cols_df[cols_to_keep]
    merged_df = pd.concat([ts_indices, vars_df], axis=1)
    merged_df.columns = var_names
    
    # print(all_cols_df)
    # print(ts_indices)
    # print(vars_df)
        
    #merged_df['forecastTime'] = pd.to_datetime(merged_df['time']) + pd.to_timedelta(merged_df['step'])

    return merged_df

### Downloads gfs data into directories that mirrors the GFS directory structure
# -- dates : list of forecast dates to download
# -- hours : list of forecast hours to download
def pull_gribs(dates = generate_date_strings('20230828',1), hours = generate_hours_list(0)):
    
    # make the subdirectory for the raw data
    if not os.path.exists(fc_data_dir): os.makedirs(fc_data_dir)
    origDir = os.getcwd()
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
        os.chdir(fc_data_dir)
    os.chdir(origDir)
    return


### In-place function that transforms the longitude indices from 0-360 t0 -180-180
def remap_longs(df):
    map_function = lambda lon: (lon - 360) if (lon > 180) else lon
    longitudes = df.index.levels[1]
    remapped_longitudes = longitudes.map(map_function)
    df.index = df.index.set_levels(remapped_longitudes, level="longitude")
    return
