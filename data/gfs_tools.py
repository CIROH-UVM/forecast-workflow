import cfgrib
from datetime import datetime, timedelta
import glob
from lib import download_data
import numpy as np
import os
import pandas as pd
import subprocess as sp
import xarray as xr

### global vars that wouldn't change based on user; may change depending on what forecast data is being pulled
gfs_file_template = "gfs.t00z.pgrb2.0p25.f"
# where the raw grib2 files will be stored
gfs_dir = "/data/forecastData/gfs/"

############# Functions for Processing GFS grib data ############################

# Define alias for aggreate_df_dict
def get_data(gfs_dir = f'/data/forecastData/gfs/gfs.{datetime.today().strftime("%Y%m%d")}/00/atmos/',
             location_dict = {"401": (45.0, -73.25),
                              "402": (44.75, -73.25),
                              "403": (44.75, -73.25)
                             }
):
    return aggregate_station_df_dict(gfs_dir=gfs_dir, location_dict=location_dict)

def aggregate_station_df_dict(gfs_dir = f'/data/forecastData/gfs/gfs.{datetime.today().strftime("%Y%m%d")}/00/atmos/',
							location_dict = {"401": (45.0, -73.25),"402": (44.75, -73.25),"403": (44.75, -73.25)}
							):
	station_dict = {}
	grib_file_list = glob.glob(f'{gfs_dir}gfs.t00z.pgrb2.0p25.f[0-9][0-9][0-9]')
	for grib_file in grib_file_list:
		grib_datasets = cfgrib.open_datasets(grib_file)
		datasets_indices_to_drop = []
		coords_to_drop = ["step", "atmosphere", "heightAboveGround", "surface", "time", "t"]
		grib_datasets = [ds.drop_vars(coords_to_drop, errors='ignore') for ds in grib_datasets]
		for i, ds in enumerate(grib_datasets):
			for var_name, var in ds.variables.items():
				if var_name == 'tcc' and var.attrs['GRIB_stepType'] == 'avg' and var.attrs['GRIB_typeOfLevel'] == 'atmosphere':
					datasets_indices_to_drop.append(i)
				if var_name == 'prate' and var.attrs['GRIB_stepType'] == 'avg' and var.attrs['GRIB_typeOfLevel'] == 'surface':
					datasets_indices_to_drop.append(i)
		grib_datasets = [grib_datasets[i] for i, _ in enumerate(grib_datasets) if i not in datasets_indices_to_drop]
		merged_ds = xr.merge(grib_datasets)
		# create a downward short-wave radiation flux for the .f000 files (since they don't have one) and set to 0
		if grib_file.endswith('.f000'):
			merged_ds['dswrf'] = 0
		locations_ds = isolate_loc_rows(ds=remap_longs(merged_ds), loc_dict=location_dict)
		locations_df = locations_ds.to_dataframe().reset_index()
		location_groups = locations_df.groupby(["latitude", "longitude"])
		location_dataframes = {name: group for name, group in location_groups}
		append_timestamp(sta_dict=station_dict, loc_dict=location_dict, loc_dfs=location_dataframes)
	return station_dict

# This loop will append each timestamp row (f000, f001, etc) to the station dataframe dict
def append_timestamp(sta_dict, loc_dict, loc_dfs):
	for stationID, loc in loc_dict.items():
		df_to_append = loc_dfs[loc].drop_duplicates().drop(columns = ['latitude','longitude'])
		if stationID in sta_dict:
			sta_dict[stationID] = pd.concat([sta_dict[stationID], df_to_append])
		else:
			sta_dict[stationID] = df_to_append

def dict_to_csv(loc_dict={}, location_dataframes={}):
    for station in loc_dict:
        location = loc_dict[station]
        filename = f"{station}_{location[0]}_{location[1]}.csv"
        location_dataframes[station].to_csv(filename)

def execute(cmd):
    popen = sp.Popen(cmd, stdout=sp.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close
    return_code = popen.wait()
    if return_code:
        raise sp.CalledProcessError(return_code, cmd)

### Creates a list of forecast dates to download
# -- start_date : first date of forecast data to download
# ---- should be a string in the format 'YYYYMMDD' or a datetime object
# -- num_dates : the number of dates ahead or behind the start date you want to download
# -- cast : str switch that determines whether to grab n dates ahead ("fore") or behind ("hind")
def generate_date_strings(start_date, num_dates=1, cast="fore"):
    date_strings = []
    # if not a datetime object, convert start_date to one
    if not isinstance(start_date, datetime):
        # strptime(str, format) converts a string to a datetime object
        current_date = datetime.strptime(start_date, "%Y%m%d")
    else:
        current_date = start_date
    for _ in range(num_dates):
        date_strings.append(current_date.strftime("%Y%m%d"))
        if cast == "hind":
            current_date -= timedelta(days=1)
        else:
            current_date += timedelta(days=1)

    return date_strings

### Creates a list of forecast hours to be downloaded
# -- num_hours: how many hours of forecast data you want: note that date goes up to 16 days out
# -- archive: boolean flag indicating if data is going to be pulled from archives; if true returns only step = 3 list
def generate_hours_list(num_hours=168, archive=False):
    if archive:
        return [f"{hour:03}" for hour in range(0, num_hours + 1, 3)]
    if not archive:
        if num_hours <= 120:
            return [f"{hour:03}" for hour in range(0, num_hours + 1)]
        else:
            return [f"{hour:03}" for hour in range(0, 120)] + [
                f"{hour:03}" for hour in range(120, num_hours + 1, 3)
            ]
        
# pulls just the desired lat/long pairs out of the datasets and isolat
def isolate_loc_rows(ds, loc_dict):
    # initialize empty list to store ds for each station in loc_dict
    station_ds_list = []
    for coords in loc_dict.values():
        # grab the ds for the lat/long pair
        station_ds = ds.sel({"latitude": coords[0], "longitude": coords[1]})
        # add the pulled ds to the station ds list
        station_ds_list.append(station_ds)
    # concatenate all of the station datasets pulled
    concat_stations_ds = xr.concat(station_ds_list, dim="latitude")
    return concat_stations_ds

### In-place function that transforms the longitude indices from 0-360 t0 -180-180
def remap_longs(ds):
    map_function = lambda lon: (lon - 360) if (lon > 180) else lon
    vector_fcn = np.vectorize(map_function)
    longitudes = ds.coords["longitude"].values
    remapped_longitudes = vector_fcn(longitudes)
    remapped_ds = ds.assign_coords(longitude=remapped_longitudes)
    return remapped_ds

################## Function for downloading GFS data ##############################

### will download the grib files for the specified dates and hours
# -- dates (list of strs) [opt]: list of dates to download gribs for. Default is just the current date
# -- hours (list of strs) [opt]: list of forecast hours to download gribs for. Default is 7 day forecast, or 168 hours
# -- log (logger) [required]: logger track info and download progress
# -- grib_data_dir (str) [opt]: directory in which to store gfs grib files
def download_gfs(log, dates=generate_date_strings(start_date=datetime.today().strftime("%Y%m%d"), num_dates=1), hours=generate_hours_list(168), grib_data_dir="/data/forecastData/gfs"):
	log.info(f'TASK INITIATED: Download {int(hours[-1])}-hour GFS forecasts for the following dates: {dates}')
	for d in dates:
		log.info(f'DOWNLOADING GFS DATA FOR DATE {d}')
		date_dir = os.path.join(grib_data_dir, f'gfs.{d}/00/atmos')
		if not os.path.exists(date_dir):
			os.makedirs(date_dir)
		for h in hours:
			grib_destination = os.path.join(grib_data_dir, date_dir, f'gfs.t00z.pgrb2.0p25.f{h}')
			grib_url = f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.{d}%2F00%2Fatmos&file=gfs.t00z.pgrb2.0p25.f{h}&var_CPOFP=on&var_DSWRF=on&var_PRATE=on&var_RH=on&var_TCDC=on&var_TMP=on&var_UGRD=on&var_VGRD=on&lev_2_m_above_ground=on&lev_10_m_above_ground=on&lev_surface=on&lev_entire_atmosphere=on&subregion=&toplat=47.5&leftlon=280&rightlon=293.25&bottomlat=40.25"
			download_data(url=grib_url, filepath=grib_destination, log=log)
	log.info('TASK COMPLETE: GFS DOWNLOAD')

############################# OLD FUNCTIONS - NOT TO BE USED ##########################

# def aggregate_df_dict(
#     files = [],
#     dates=[datetime.today().strftime("%Y%m%d")],
#     loc_dict={"401": (45.0, -73.25), "402": (44.75, -73.25), "403": (44.75, -73.25)},
# ):
#     arg_vars = {
#         ("atmosphere", "instant"): ["tcc"],
#         ("heightAboveGround", 10): ["u10", "v10"],
#         ("heightAboveGround", 2): ["t2m", "r2"],
#         ("surface", "avg"): ["dswrf"],
#         ("surface", "instant"): ["cpofp", "prate"],
#     }
#     drop_coords = ["step", "atmosphere", "heightAboveGround", "surface"]
#     var_names = ["T2", "TCDC", "SWDOWN", "U10", "V10", "RH2", "RAIN", "CPOFP"]
#     # initialize a dictionary to store a dataframe for each station in
#     data_dict = {}

#     # Don't think we need this anymore now that downloads are seperate
#     # if not os.path.exists(fc_data_dir): os.makedirs(fc_data_dir)
#     origDir = os.getcwd()
#     # move into raw_fc_Data/
#     os.chdir(fc_data_dir)

#     # for every date in the list of date strings dates
#     for d in dates:
#         # create and navigate to the dir containing the grib2 files for the current date d
#         date_dir = "gfs." + d + fc_time
#         os.chdir(date_dir)
#         print(
#             "aggregating all gribs for date {}:".format(
#                 datetime.strptime(d, "%Y%m%d").strftime("%Y-%m-%d")
#             )
#         )

#         # this chunk should go in aggregate_df_dict, as it creates a dataset for each arg combo, 5 total
#         # initalize empty dataset list
#         total_ds_list = []
#         # for each arg combo in args_list, open all grib2 files and add the aggregate dataset to the list
#         for av in arg_vars:
#             args = {"typeOfLevel": av[0]}
#             # if the type of level is heigh above ground, then the 2nd arg is topLevel, not stepType
#             if args["typeOfLevel"] == "heightAboveGround":
#                 args["topLevel"] = av[1]
#             else:
#                 args["stepType"] = av[1]
#             print("\topening gribs with args dict: {}".format(args))
#             preprocess_with_args = functools.partial(
#                 preprocesses, date=d, loc_dict=loc_dict
#             )
#             files=glob.glob("gfs.t00z.pgrb2.0p25.f[0-9][0-9][0-9]")
#             ds = xr.open_mfdataset(
#                 files,
#                 engine="cfgrib",
#                 compat="override",
#                 coords="minimal",
#                 parallel=True,
#                 preprocess=preprocess_with_args,
#                 combine="nested",
#                 concat_dim="valid_time",
#                 backend_kwargs={"indexpath": "", "filter_by_keys": args},
#             )
#             total_ds_list.append(ds)

#         print("successfully opened all gribs")
#         # enumerate through arg_vars again to only keep the vars we want from each arg combo
#         total_filtered_ds_list = []

#         # we want to enumerate so that we have an index i to access the datasets in total_ds_list
#         # keep_vars will be the list of vars to keep; the values of the arg_vars dict
#         for i, keep_vars in enumerate(arg_vars.values()):
#             # if a var is not in the keep list, add it to the drop list
#             drop_vars = [
#                 var for var in total_ds_list[i].data_vars if var not in keep_vars
#             ]
#             # create new dataset that has dropped the vars in the drop list
#             ds_filtered = total_ds_list[i].drop_vars(drop_vars)
#             # add new filtered dataseets to a new list
#             total_filtered_ds_list.append(ds_filtered)

#         # now concat all of the filtered datasets into one
#         final_ds = xr.concat(
#             total_filtered_ds_list, dim="latitude", coords="minimal", compat="override"
#         ).drop_vars(drop_coords)
#         # convert to a dataframe and drop the time column - we only need valid_time, which contains the time series
#         final_df = final_ds.to_dataframe().drop("time", axis=1)
#         # Rename valid_time (already an index, along with latitude) to time
#         #   to keep consistent with naming in old extraction algorithm
#         final_df.index = final_df.index.set_names('time', level=1)

#         print("successfully aggregated all datasets")

#         # reset index to drop all duplicate rows, then reset and sort index
#         final_df = (
#             final_df.reset_index()
#             .drop_duplicates()
#             .set_index(["time", "latitude", "longitude"])
#             .sort_index()
#         )
#         # group by unique valid time, lat/long combinations, collapse null values
#         final_df = final_df.groupby(
#             ["time", "latitude", "longitude"], as_index=True
#         ).first()
#         # rename columns according to specified variable names
#         final_df.columns = [
#             "TCDC",
#             "U10",
#             "V10",
#             "T2",
#             "RH2",
#             "SWDOWN",
#             "CPOFP",
#             "RAIN",
#         ]
#         # re order columns according to specified var order
#         final_df = final_df.reindex(columns=var_names)

#         # now groupby unique lat/long pairs for location data extraction
#         loc_groups = final_df.groupby(["latitude", "longitude"])
#         # create a dict where tje key is the lat/long pair and value is corresponding df
#         loc_dfs = {name: group for name, group in loc_groups}

#         # for every station in loc_dict...
#         for station in loc_dict:
#             coords = loc_dict[station]
#             # pull the dataframe for the station, drop lat/long indices
#             df = (
#                 loc_dfs[coords]
#                 .reset_index(["latitude", "longitude"])
#                 .drop(["latitude", "longitude"], axis=1)
#             )
#             # if there is already a dataframe for the station (i.e. from a previous date)
#             if station in data_dict:
#                 # add the new date df to the last date df
#                 data_dict[station] = pd.concat([data_dict[station], df])
#             # if there is no dict entry for the station, enter the df
#             else:
#                 data_dict[station] = df

#         print("successfully merged and extracted loc data to dictionary\n")

#         # TODO: Can we remove this os.chdir on the next line?
#         os.chdir(fc_data_dir)
#     os.chdir(origDir)
#     return data_dict

# # pulls just the desired lat/long pairs out of the datasets and isolat
# def concat_locs(ds, loc_dict):
#     # initialize empty list to store ds for each station in loc_dict
#     station_ds_list = []
#     for coords in loc_dict.values():
#         # grab the ds for the lat/long pair
#         station_ds = ds.sel({"latitude": coords[0], "longitude": coords[1]})
#         # add the pulled ds to the station ds list
#         station_ds_list.append(station_ds)
#     # concatenate all of the station datasets pulled
#     concat_stations_ds = xr.concat(station_ds_list, dim="latitude")
#     return concat_stations_ds

# ### calls the curl command to the terminal to download a file from the web
# # -- url (str) [required]: complete url of the file to be downloaded
# # -- file_path (str) [required]: complete* path to the destination file. *should inlcude the desired name of the download file at the end of the path
# # -- silent (bool) [optional]: switch to turn off progress report. Default is false
# def curl(url, file_path, silent = False):
# 	if not silent:
# 		print(f'Downloading file from: {url}')
# 		print(f'\tto destination: {file_path}')
# 		# -C - enables curl to pickup download where it left off in case of connection disruption
# 	sh.curl('-o',file_path, '-C', '-', url)
# 	if not silent: print('Download complete\n')

# ### Given the transformed dataframe and a list of lat/long tuples to extract, returns new df containing just the rows for each lat/long pair
# # -- df : the grib2 df post-long-transform
# # -- loc_dict : dict of of station names and corresponding lat/long tuples; will pull said coords from the datafram
# # given the transformed dataframe, adds the new rows to the respective composite dataframes in the df dictionary
# def extract_locs(df, loc_dict={}, location_dataframes={}):
#     for station in loc_dict:
#         coords = loc_dict[station]
#         extracted_df = pd.DataFrame(df.loc[coords]).T.set_index("time")
#         # if the dictionary does not already have a key for each location, then initialize that key; should only be true when extracting the 1st df
#         if len(location_dataframes) != len(loc_dict):
#             location_dataframes[station] = extracted_df
#         else:
#             location_dataframes[station] = pd.concat(
#                 [location_dataframes[station], extracted_df]
#             )

# def extract_subgrib(fname="", args={}):
#     ds = xr.open_dataset(
#         fname, engine="cfgrib", backend_kwargs={"filter_by_keys": args}
#     )
#     longnames = ["time", "step", args["typeOfLevel"], "valid_time"]
#     for v in ds:
#         longnames.append(
#             "{}, {}".format(ds[v].attrs["long_name"], ds[v].attrs["units"])
#         )
#     df = ds.to_dataframe()
#     df.columns = longnames
#     remap_longs(df)
#     return df

# def get_subgrib_df_list(fname=""):
#     # drop surface avg args for the first (f000) file
#     if fname[-3:] == "000":
#         args_list = [
#             {"typeOfLevel": "atmosphere", "stepType": "instant"},
#             {"typeOfLevel": "heightAboveGround", "topLevel": 10},
#             {"typeOfLevel": "heightAboveGround", "topLevel": 2},
#             {"typeOfLevel": "surface", "stepType": "instant"},
#         ]
#     else:
#         args_list = [
#             {"typeOfLevel": "atmosphere", "stepType": "instant"},
#             {"typeOfLevel": "heightAboveGround", "topLevel": 10},
#             {"typeOfLevel": "heightAboveGround", "topLevel": 2},
#             {"typeOfLevel": "surface", "stepType": "avg"},
#             {"typeOfLevel": "surface", "stepType": "instant"},
#         ]
#     # temporary list to store un-extracted dataframes
#     df_list = []
#     for args in args_list:
#         subgrib_df = extract_subgrib(fname, args)
#         df_list.append(subgrib_df)
#     return df_list

# def merge_subgrib_dfs(df_list=[], fname=""):
#     cols_to_keep = [
#         "2 metre temperature, K",
#         "Total Cloud Cover, %",
#         "Downward short-wave radiation flux, W m**-2",
#         "10 metre U wind component, m s**-1",
#         "10 metre V wind component, m s**-1",
#         "2 metre relative humidity, %",
#         "Precipitation rate, kg m**-2 s**-1",
#         "Percent frozen precipitation, %",
#     ]
#     var_names = ["time", "T2", "TCDC", "SWDOWN", "U10", "V10", "RH2", "RAIN", "CPOFP"]

#     # if not the first file, drop surface avg precip rate
#     if not fname[-3:] == "000":
#         # drop the surface average precip rate; it has the same col name as the surface instant precip col, which is problematic b/c we just need the latter
#         df_list[3] = df_list[3].drop("Precipitation rate, kg m**-2 s**-1", axis=1)
#     # combine all dfs in list
#     all_cols_df = pd.concat(df_list, axis=1)
#     # keep valid_time indices at beginning of df
#     ts_indices = all_cols_df.iloc[:, [3]]
#     # create frozen precip if it does not exist
#     if not "Downward short-wave radiation flux, W m**-2" in all_cols_df:
#         all_cols_df["Downward short-wave radiation flux, W m**-2"] = 0
#     # get the vars we need, in order
#     vars_df = all_cols_df[cols_to_keep]
#     merged_df = pd.concat([ts_indices, vars_df], axis=1)
#     merged_df.columns = var_names

#     # print(all_cols_df)
#     # print(ts_indices)
#     # print(vars_df)

#     # merged_df['forecastTime'] = pd.to_datetime(merged_df['time']) + pd.to_timedelta(merged_df['step'])

#     return merged_df

# # mainly isolates the desired lat/long pairs and returns datasets containing only those coords rather than whole globe
# # also creates a ds with valid time set to the current date at midnight and downward short wave UV to 0; this deals with the issue of f000 grib2 datasets not having
# # any data for the arg combo {'typeOfLevel':'surface', 'stepType':'avg'}
# def preprocesses(ds, date, loc_dict):
#     # if the ds is empty...
#     if len(ds) == 0:
#         # create valid_time coord and set to the date, at midnight
#         dt = np.datetime64(
#             datetime.strptime(date, "%Y%m%d").replace(
#                 hour=0, minute=0, second=0, microsecond=0
#             )
#         )
#         # Suppress the warning message
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             ds = ds.assign_coords(valid_time=dt)
#         # add a column for downward short-wave radiation and set to 0
#         ds["dswrf"] = 0
#         return ds
#     else:
#         ds = concat_locs(remap_longs(ds), loc_dict)
#         return ds

# ### Downloads gfs data into directories that mirrors the GFS directory structure
# # -- dates : list of forecast dates to download
# # -- hours : list of forecast hours to download
# def pull_gribs(
#     dates=generate_date_strings("20230828", 1), hours=generate_hours_list(0)
# ):
#     # make the subdirectory for the gfs data
#     if not os.path.exists(fc_data_dir):
#         os.makedirs(fc_data_dir)
#     origDir = os.getcwd()
#     os.chdir(fc_data_dir)

#     for d in dates:
#         date_dir = "gfs." + d + fc_time
#         if not os.path.exists(date_dir):
#             os.makedirs(date_dir)
#         os.chdir(date_dir)
#         date_url = gfs_root + date_dir + fc_file
#         for h in hours:
#             hour_url = date_url + h
#             if not os.path.exists(fc_file + h):
#                 print("Downloading file from URL:", hour_url)
#                 for path in execute(
#                     ["curl", "--silent", "--connect-timeout", "120", "-O", hour_url]
#                 ):
#                     print(path, end="")
#                 print("Download Complete:", date_dir + fc_file + h, "\n")
#                 # pull_call = sp.run(['curl', '-O', hour_url], capture_output = True, check = True)
#         os.chdir(fc_data_dir)
#     os.chdir(origDir)