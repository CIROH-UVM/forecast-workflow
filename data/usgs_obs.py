import requests
import json
import pandas as pd

# USGS parameter codes
# Streamflow: '00060',
# Gage Height: '00065'

def USGSstreamflow_function(station_id, parameter, period):
    gage = requests.get('https://waterservices.usgs.gov/nwis/iv/'
                      '?format=json'
                     f'&sites={station_id}'
                     f'&period={period}'
                     f'&parameterCd={parameter}'
                     )
    values = gage.json()['value']['timeSeries'][0]['values'][0]['value']
    df = pd.DataFrame(values)
    # print(df)
    # df = df.set_index('dateTime')
    # df = df.set_index(pd.to_datetime(df['dateTime']))
    # df = df.drop(['dateTime','qualifiers'],axis =1)
    # df.columns = ['streamflow']
    # df.to_csv(station_id+"_flow.csv", sep=',')
    return pd.DataFrame(data={'streamflow': df['value'].values}, index=pd.to_datetime(df['dateTime']).dt.tz_localize(None))

def get_data(station_ids = ['04294000', '04292810', '04292750']):

    # 04294000 (MS), 04292810 (J-S), 04292750 (Mill)
    parameter = '00060'
    # Get 90 Days Prior
    period = 'P90D'
    returnVal = {}

    for station_id in station_ids:
        returnVal[station_id] = USGSstreamflow_function(station_id, parameter, period)
    
    return returnVal
