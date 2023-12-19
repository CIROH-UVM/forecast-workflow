import requests
import json
import pandas as pd
import datetime as dt

# USGS parameter codes
# https://help.waterdata.usgs.gov/codes-and-parameters/parameters
# https://help.waterdata.usgs.gov/parameter_cd?group_cd=PHY
# Streamflow, mean. daily in cubic ft / sec: '00060',
# Streamflow, instantaneous cubic ft / sec: '00061',
# Gage Height, feet: '00065'

def USGSstreamflow_function(station_id, parameter, startDate, endDate):
    gage = requests.get('https://waterservices.usgs.gov/nwis/iv/'
                      '?format=json'
                     f'&sites={station_id}'
#                     f'&period={period}'
                     f'&startDT={startDate.strftime("%Y-%m-%d")}'
                     f'&endDT={endDate.strftime("%Y-%m-%d")}'                     
                     f'&parameterCd={parameter}'
                     )
    #print(gage.text)
    values = gage.json()['value']['timeSeries'][0]['values'][0]['value']
    df = pd.DataFrame(values)
    # print(df)
    # df = df.set_index('dateTime')
    # df = df.set_index(pd.to_datetime(df['dateTime']))
    # df = df.drop(['dateTime','qualifiers'],axis =1)
    # df.columns = ['streamflow']
    # df.to_csv(station_id+"_flow.csv", sep=',')
    # 'US/Eastern' is the other option, but what about fall daylight savings "fall back"
    #return pd.DataFrame(data={'streamflow': df['value'].values}, index=pd.to_datetime(df['dateTime'], utc=True).dt.tz_convert('Etc/GMT+4').dt.tz_localize(None))
    # 20231211 - set index as datetime with timezone suffix set to UTC
    return pd.DataFrame(data={'streamflow': df['value'].values}, index=pd.to_datetime(df['dateTime'], utc=True))

def get_data(ForecastStartDate, SpinupStartDate, station_ids = ['04294000', '04292810', '04292750']):

    # 04294000 (MS), 04292810 (J-S), 04292750 (Mill)
    parameter = '00060'
    # Get 90 Days Prior
    period = 'P90D'
    returnVal = {}

    # 20231211 - do not adjust passed dates to a previous day. that is a caller concern if that additional data buffer is needed.
    for station_id in station_ids:
        returnVal[station_id] = USGSstreamflow_function(station_id,
                                                        parameter,
                                                        SpinupStartDate,
                                                        ForecastStartDate
                                                        )
    
    return returnVal
