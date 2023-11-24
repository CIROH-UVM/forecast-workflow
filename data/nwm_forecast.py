from datetime import datetime, timedelta
from lib import download_data
import os
import pandas as pd
import xarray as xr

# Adapted from python notebooks at https://www.hydroshare.org/resource/5949aec47b484e689573beeb004a2917/

# StartDate = '20230907'
StartDate = datetime.now().strftime('%Y%m%d')
StartTimestep = '00'
forecast_files_path = '/data/forecastData/nwm'

## Example Uses
#### To Download
# download_files = download_forecast_files(ForecastStartDate=StartDate, ForecastStartTimestep=StartTimestep, download_dir=forecast_files_path)

#### To generate dict with data
# reach_Series_Data = get_data(ForecastStartDate=StartDate, ForecastStartTimestep=StartTimestep, download_files=download_files)


def GetForecastFileName(ForecastStartDate = '20230907', ForecastStartTimestep='00', 
                        ForecastType = 'medium_range', ForecastMember='1', TimeStep = '001'):
  
    """
    
    A Function to construct a URL for Downloading a Specific Forecast File from NOAA NWM
    
    Args:
    ForecastStartDate     : The date for which the forecast is needed. Default is '20230907'.
    ForecastStartTimestep : The hour at which the forecast starts. Default is '00' (midnight).
    ForecastType          : Specifies the type of forecast, Default is medium_range.
    ForecastMember        : Represents which member of the forecast model we want. Default is '1' ---> 'medium_range_mem1'. 
    TimeStep              : Represents the specific forecast file within the range. 

    Return:
    Complete URL for the File. 
    
    """
    BaseName = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.'
    return BaseName + ForecastStartDate + '/medium_range_mem' + ForecastMember + '/nwm.t' + ForecastStartTimestep + 'z.medium_range.channel_rt_' + ForecastMember + '.f' + TimeStep + '.conus.nc'

  
def GetForecastFile(Log, Url, download_dir='.'):
    """
    
    A Function to download the given forecast url from NOAA NWM.

    Args:
    Log          : logging.Logger object - to log messages for data download
    Url          : The URL of the forecast file to be downloaded, should be string.
    download_dir: The directory where the file should be saved. 

    Returns:
    A String Path of the Downloaded File. 
    
    """
    FileName = os.path.basename(Url)
    # Lets construct the complete file path
    FilePath = os.path.join(download_dir, FileName)
    # Lets create download_path if it doesn't exist yet
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)    
    # Lets make sure the file is not already downloaded - If yes, we will not download it again.
    if os.path.exists(FilePath):
        Log.info(f'Skipping download: "{FilePath}" already exists')
        return None
       
    print(f'Downloading {FileName}')

    # Lets do the download
    # r = requests.get(Url, allow_redirects=True)
    # # Time to save the file to the download_dir. 
    # open(FilePath, 'wb').write(r.content)
    
    # sh.curl('-o', FilePath,'-C','-', Url)

    download_data(Url, FilePath, Log)
    
    return FilePath

# Next we will define a function which will call these two function and download the data. 
def download_forecast_files(Log, ForecastStartDate, ForecastStartTimestep='00', 
                            ForecastType='medium_range', ForecastMember='1', data_dir='forecastData/nwm/'):
    """
    
    A Function to download the forecast data for a given start date and start time step. It will first call the GetForecastFileName()
    to get the URL of the Filename to be downloaded. That URL will further be passed to GetForecastFile() to download the file. 

    Args:
    Log                  : logging.Logger object - to log messages for data download
    ForecastStartDate    : The starting date for which forecasts are to be downloaded.
    ForecastStartTimestep: The starting time for the forecasts.
    ForecastType         : The type of forecast (default is 'medium_range').
    ForecastMember       : The member of the forecast model.
    data_dir             : Path in which to build the NWM data download subdirectroy structure.

    Returns:
    List of paths of downloaded files.
    
    """
    # First define the timestamps
    time_stamps = ['%03d' % (i+1) for i in range(240)]  # For 10 days time Stamps

    Log.info(f'TASK INITIATED: Download {int(time_stamps[-1])}-hour NWM hydrology forecasts for the following date: {ForecastStartDate[4:6]}/{ForecastStartDate[6:8]}/{ForecastStartDate[0:4]}')
    # Lets create an empty list to store the complete path of downloaded files
    download_files = []
    for time_stamp in time_stamps:
        # Getting the URL
        url = GetForecastFileName(ForecastStartDate=ForecastStartDate, TimeStep=time_stamp)
        # create the directroy structure in which to download the data - mirrors the NWM URL structure
        download_dir = os.path.join(data_dir, f'nwm.{ForecastStartDate}/{ForecastType}_mem{ForecastMember}')
        # Lets download the url file - The function will return the filepath which we will append to download_files
        file_path = GetForecastFile(Log=Log, Url=url, download_dir=download_dir)
        if file_path is not None:
            download_files.append(file_path)
    Log.info('TASK COMPLETE: NWM DOWNLOAD')
    # Lets return the list
    return download_files

# Lastly, lets read all the downloaded files and process them into a nice Dictionary, where the Key will be the Reach Name and 
# values will be Pandas Series
def get_data(ForecastStartDate, ForecastStartTimestep, data_dir='forecastData/nwm/', save_csv=False, ForecastType='medium_range', ForecastMember='1', directory_flag = True):
    """

    A Function to Process the Downloaded data. It will read each individual file, extract the Stream Value and Add it to 
    the Dictionary Values against Reach Names.
    
    Args:
    ForecastStartDate     : The date for which the forecast is needed - This is needed to extract the Datetime value.
    ForecastStartTimestep : The starting time for the forecasts - This is needed to extract the Datetime value. 
    data_dir              : Path in which to build the NWM data download subdirectroy structure.
    save_csv              : Flag indicating whether or not dataframes should be saved as CSV files.
    ForecastType          : The type of forecast (default is 'medium_range').
    ForecastMember        : The member of the forecast model.
    
    """
    # Lets define an empty dictionary to store the Results.
    results  = {}

    if directory_flag:
        nwm_download_path = f'nwm.{ForecastStartDate}/{ForecastType}_mem{ForecastMember}'
        nwm_csv_save_path = f'nwm.{ForecastStartDate}'
    else:
        nwm_download_path = f'{ForecastStartDate}/'
        nwm_csv_save_path = f'{ForecastStartDate}'

    # Append ForecastStartDate to download_base_path
    download_base_path = os.path.join(data_dir, nwm_download_path)

    # Get the filenames from download_base_path
    download_files = sorted([f for f in os.listdir(download_base_path) if f.endswith('.nc')])

    # Lets initiate an another dictionary where key is reach name from reaches and value is an empty list. 
    # We will add the StreamFlow values later to this list.
    reaches = {
    
                166176984: "MS",
                4587092: "J-S",
                4587100: "Mill"
              }
    
    data_dict = {reach_name:[] for reach_id, reach_name in reaches.items()}
   

    # Lets initiate an empty list to store the timestamps info
    timestamps = []

    # Time to read the data files
    for file in download_files:
        data = xr.open_dataset(os.path.join(download_base_path, file))

        # Time to extract the relevent reach id data from data
        for reach_id, reach_name in reaches.items():
            stream_value  = float(data.sel(feature_id=reach_id).streamflow.values)
            # Time to append it in data_dictionary
            data_dict[reach_name].append(stream_value)

        # Now lets extract the Timestamp data from file.
        time_stamp = os.path.basename(file).split('f')[-1].split('.')[0]
        hour_value = int(time_stamp)
        timestamps.append(datetime.strptime(ForecastStartDate + ForecastStartTimestep, '%Y%m%d%H') + timedelta(hours=hour_value))
        print(f'NWMFile: {file} | TimeStamp: {datetime.strptime(ForecastStartDate + ForecastStartTimestep, "%Y%m%d%H") + timedelta(hours=hour_value)}')

    # At step we have all the data need to get the final format - which {'Reach_Name':pd.Series(streamflow, index=timestamp)}
    for reach_name, series_data in data_dict.items():
        results[reach_name] = pd.DataFrame(data={'streamflow': series_data}, index=timestamps)
        if save_csv:
            filename = os.path.join(data_dir, f'{nwm_csv_save_path}/nwm_{reach_name}.csv')
            results[reach_name].to_csv(filename)

    # Simply return the results
    return results