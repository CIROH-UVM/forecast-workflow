import os
import csv
import pandas as pd

def get_data():
    with open(os.path.join("/data/forecastData/colchesterReefFEMC", "Z0080_CR_QAQC.csv")) as file:
        cr_data = pd.read_csv(file, delimiter=",", header=0, index_col=0) # Make DateTime as index
        cr_df = pd.DataFrame(cr_data)
        # Keep the last 90 days (24h*4quarters*90days=8640) plus buffer
        cr_df = cr_df.tail(8800)
        # Subset to the columns of interest
        cr_df = cr_df[['38m_AIRTEMP','PYRANOM','38m_RELHUMID', 'NRG_38m_MEAN_RESULTANT_WINDSPEED','NRG_38m_MEAN_WIND_DIRECTION']]
        cr_df.columns = ['T2C', 'SWDOWN', 'RH2', 'WSPEED', 'WDIR']
    return cr_df
