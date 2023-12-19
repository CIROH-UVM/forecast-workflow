import os
import csv
import pandas as pd
import datetime as dt

def get_data(ForecastStartDate, SpinupStartDate):
    with open(os.path.join("/data/forecastData/colchesterReefFEMC", "Z0080_CR_QAQC.csv")) as file:
        hcr_data = pd.read_csv(file, delimiter=",", header=0, index_col=0, parse_dates=True) # Make DateTime as index
        hcr_df = pd.DataFrame(hcr_data)

        # Also get the URL for the most recent Colchester Reef data in the CSV. Concat the two frames.
        url="https://uvm.edu/femc/MetData/ColReefQAQC/CR_QAQC_latest.csv"
        scr_df = pd.DataFrame(pd.read_csv(url, delimiter=",", header=0, index_col=0, parse_dates=True))

        #print(scr_df)

        #cr_df = pd.concat([hcr_df, scr_df], axis=0, ignore_index=True)   # ignore dupe time indexes (overlap)
        cr_df = pd.concat([hcr_df, scr_df], axis=0)
        # Drop duplicate time indices from overlap between 2 .csvs -- And sort
        cr_df = cr_df[~cr_df.index.duplicated(keep='first')].sort_index()

        # Before: Keep the last 90 days (24h*4quarters*90days=8640) plus buffer
        # cr_df = cr_df.tail(8800)
        # Now: Keep the days from SpinupStartDate to present

        # 12231211 - remove shift of start date to the previous day. That's an AEM3D specific concern
        cr_df = cr_df[cr_df.index > dt.datetime.combine(SpinupStartDate, dt.datetime.min.time())]

        # Drop rows after midnight today to prevent duplicate entries with forecast
        cr_df = cr_df[cr_df.index < dt.datetime.combine(ForecastStartDate, dt.datetime.min.time())]

        # generate index as datetime with timezone suffix in UTC time
        cr_df.set_index(cr_df.index.tz_localize('EST').tz_convert('UTC'), inplace = True)
        
        # Subset to the columns of interest
        cr_df = cr_df[['38m_AIRTEMP','PYRANOM','38m_RELHUMID', 'NRG_38m_MEAN_RESULTANT_WINDSPEED','NRG_38m_MEAN_WIND_DIRECTION']]
        cr_df.columns = ['T2', 'SWDOWN', 'RH2', 'WSPEED', 'WDIR']
    return cr_df
