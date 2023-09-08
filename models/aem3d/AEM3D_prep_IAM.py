#  Creates the AEM3D Lake Model Input Contents
#
#  Time Series
#       Flow from RHESSys or SWAT
#       Weather from WRF
#         Net Longwave Radiation
#         Precipitation
#         Temperature
#		  Wind
#         Humidity
#       Lake Level (flow and temp related)
#       Salinity
#       Tracers
#       Call waterquality.py to wq files
#
#  Control File

from ...lib import *
from ...NWM_forecast import *
from ...USGS_obs import *
from ...gfs_download_fcns import *
from ...colchesterReefData import *
from ...btvData import *
from .waterquality import *

from sh import cp, tar, mkdir, mv, Rscript
from string import Template
import pandas as pd
import numpy as np
import glob
import os
import datetime


def colsToHeader(columns):
    return_string = ''
    for col in columns:
        return_string += ' ' + str(col)
    return return_string

def datetimeToOrdinal(date):

    dayofyear = date.strftime('%j')

    # Left pad dayofyear to length 3 by zeros
    yearday = str(date.year) + dayofyear.zfill(3)

    totseconds = date.hour * 3600 + \
                 date.minute * 60 + \
                 date.second
    fracsec = totseconds / datetime.timedelta(days=1).total_seconds()  #Fraction of the day's seconds

    ordinaldate = yearday + str(fracsec)[1:6].ljust(5,'0')  # add the percentage seconds since noon
    return ordinaldate

## Deprecating this in favor of more generic datetimeToOrdinal() and then using .apply()

# def pandasDatetimeToOrdinal(pd_dt_index):

#     # This doesn't work for dayofyear because Dec 31 in leap year is 366
#     #dayofyear = pd.Series(pd_dt_index.strftime('%j'))

#     year = pd.unique(pd_dt_index.year)
#     if not len(year) == 1:
#         raise Exception('DateTimeIndex must have only one year for conversion to Ordinal')
#     year = year[0]

#     # So... make own dayofyear sequence from hourly dt sequence, adjusted for indexing from 1
#     dayofyear = 1 + np.floor(pd.Series(range(len(pd_dt_index)))/24).astype(int)
#     # Left pad dayofyear to length 3 by zeros
#     yearday = str(year) + dayofyear.apply(str).apply(lambda x: x.zfill(3))

#     totseconds = pd.Series(pd_dt_index.hour * 3600 + \
#                 pd_dt_index.minute * 60 + \
#                 pd_dt_index.second)
#     fracsec = totseconds / pd.Timedelta(days=1).total_seconds()  #Fraction of the day's seconds

#     ordinaldate = yearday + fracsec.astype(str).str.slice(start=1,stop=6)  # add the percentage seconds since noon
#     return ordinaldate.str.pad(width=12, side='right', fillchar='0')  # pretty right edge


def seriesIndexToOrdinalDate(series):
    # Now, using datetimeToOrdinal()
    ordinaldate = series.index.to_series().apply(datetimeToOrdinal)

    #ordinaldate = pandasDatetimeToOrdinal(series.index)
    #ordinaldate = pd.Series(wrfdf['ordinaldate'].array, index = wrfdf['wrftime'])
    return pd.Series(series.array, index = ordinaldate)

def writeFile(filename, bayid, zone, varName, dataSeries):
    with open(filename, mode='w', newline='') as output_file:

        # output the header text
        output_file.write('!-----------------------------------------------------!\n')
        output_file.write('! Written by AEM3D_prep_IAM                           !\n')
        output_file.write(f'! Bay ID: {bayid}                                 !\n')
        #output_file.write('! Bay Source: ' + bs_name + '                         !\n')
        output_file.write('!-----------------------------------------------------!\n')
        output_file.write('1 data sets\n')
        output_file.write('0 seconds between data\n')
        output_file.write(f'0            {zone}\n')
        output_file.write(f'TIME         {varName}\n')

        dataSeries.to_csv(path_or_buf = output_file, float_format='%.3f', sep=' ', index=True, header=False)


# Doesn't work... need to now calculate in climate_lib
# def writeLongwaveRadiationNet(climate, THEBAY):
#     ###########################################################################################
#     #
#     #   Net Longwave Radiation : GLW - LWUP
#     #

#     netrad_df = dict()
#     for zone in climate['GLW'].keys():
#         # Net Longwave = Downward Longwave at Ground (GLW) - Longwave Rad Up (LWUP)
#         netrad_df[zone] = climate['GLW'][zone] - climate['LWUP'][zone]

#     for zone in netrad_df.keys():
#         filename = f'LWRADNET_{zone}.dat'
#         logger.info('Generating Bay Longwave Radiation Net File: '+filename)
#         writeFile(
#             os.path.join(THEBAY.infile_dir, filename),
#             THEBAY.bayid,
#             zone,
#             "LW_RAD_NET",
#             seriesIndexToOrdinalDate(netrad_df[zone]))
#         THEBAY.addfile(fname=filename)


def writeLongwaveRadiationDownward(climate, THEBAY):
    ###########################################################################################
    #
    #   Longwave Radiation Downward : GLW
    #

    for zone in climate['AEMLW'].keys():
        filename = f'LWRADIN_{zone}.dat'
        logger.info('Generating Bay Longwave Radiation Downward File: '+filename)
        writeFile(
            os.path.join(THEBAY.infile_dir, filename),
            THEBAY.bayid,
            zone,
            "LW_RAD_IN",
            seriesIndexToOrdinalDate(climate['AEMLW'][zone]))
        THEBAY.addfile(fname=filename)


def writeCloudCover(climate, THEBAY):
    ###########################################################################################
    #
    #   Longwave Radiation Downward : GLW
    #

    for zone in climate['AEMLW'].keys():
        filename = f'CLOUDS_{zone}.dat'
        logger.info('Generating Bay Cloud Cover File: '+filename)
        writeFile(
            os.path.join(THEBAY.infile_dir, filename),
            THEBAY.bayid,
            zone,
            "CLOUDS",
            seriesIndexToOrdinalDate(climate['AEMLW'][zone]))
        THEBAY.addfile(fname=filename)


#########################################################################
#
#   Import Hydrology Flow - Generate Flow Files
#
##########################################################################


def getflowfiles(whichbay):
    '''
    getflowfiles : Get hydrology model flow file(s) for Bay Inflow
        Most information contained in passed Bay Object
    '''

    THEBAY = whichbay

    #
    #   InLandSea has 3 *basin.daily files that provide input flow
    #
    logger.info('Processing Hydrology Flow Data')
    #print (flowFiles)

    #
    #   2 years spinup and 10 years of output in the _daily file
    #   Trim extracted data to just the year of interest
    #   Don't use time stamps, because rhessys generates unwanted leap days
    #   Generate 365 records for the year to match the leapless climate data
    #
    #records2skip = (2 + THEBAY.year - year_to_decade(THEBAY.year))*365      # skip spinup and prior years

    ######### TODO: Instead of from file below, get from data gathering functions
    
    # dict by id: 04294000 (MS), 04292810 (J-S), 04292750 (Mill)
    #observedUSGS = get_USGS_data(station_ids = ['04294000', '04292810', '04292750'])
    #forecastNWM = xxxx
    
    # Need to adjust for column names
    #flowdf = pd.concat([observedUSGS['04294000'], forecastNWM['MS'])
    #mlflow = pd.concat([observedUSGS['04292750'], forecastNWM['Mill'])
    #jsflow = pd.concat([observedUSGS['04292810'], forecastNWM['J-S'])

    ##############  Remove filebased df initialization
    '''
    #  Some empty dataframes that should get filled in by found .daily files
    flowdf = pd.DataFrame()
    mlflow = pd.DataFrame()
    jsflow = pd.DataFrame()

    for fname in flowFiles:
        if fname[0:2] == 'MS' :
            # read Missisquoi daily
            #flowdf = pd.read_csv(fname, delimiter=' ', skiprows=[i for i in range(1,records2skip)], nrows=365)
            allflow = pd.read_csv(fname, delimiter=' ')
            flowdf = allflow.loc[allflow['year'] == THEBAY.year].copy().reset_index(drop=True)
            logger.info('Missisquoi Daily File Read: '+fname)

        elif fname[0:2] == 'ML' :
            # read Mill Daily
            #mlflow = pd.read_csv(fname, delimiter=' ', skiprows=[i for i in range(1,records2skip)], nrows=365)
            allflow = pd.read_csv(fname, delimiter=' ')
            mlflow = allflow.loc[allflow['year'] == THEBAY.year].copy().reset_index(drop=True)
            logger.info('Mill Daily File Read: '+fname)

        elif fname[0:2] == 'JS' :
            # read Jewitt/Stevens daily
            #jsflow = pd.read_csv(fname, delimiter=' ', skiprows=[i for i in range(1,records2skip)], nrows=365)
            allflow = pd.read_csv(fname, delimiter=' ')
            jsflow = allflow.loc[allflow['year'] == THEBAY.year].copy().reset_index(drop=True)
            logger.info('Jewitt/Stevens Daily File Read: '+fname)

        else :
            # don't recognize the flow data - throw message
            logger.info('Unrecognized River Source Flow File '+fname)

    if flowdf.empty:
        logger.info('Missisquoi Flow not read')
        sys.exit(1)
    '''

    ############### Should also be done for us
    '''
    # Combine Columns into a datetime object (for inspection, only. forcing dates to use in output
    flowdf['date_given'] = pd.to_datetime(flowdf[['year','month','day']])   # build time-stamp from separate columns
    logger.info('Hydrology Start Date: ' + flowdf['date_given'][0].strftime('%Y-%m-%d'))
    # Force the Ordinal Date as year requested with 365 days
    #   day of year is array index, adjusted for indexing from 1, left padded by zeros
    flowdf['dayofyear'] = flowdf.index + 1
    flowdf['yearday'] = str(THEBAY.year) + flowdf['dayofyear'].apply(str).apply(lambda x: x.zfill(3))
    flowdf['ordinaldate'] = flowdf['yearday'] + '.5000'  # add the fractional day (noon is at .5000 for ordinal)
    logger.info('First Ordinal Date: ' + flowdf['ordinaldate'][0])
    THEBAY.FirstDate = flowdf['ordinaldate'][0]      # Update Bay with start and stop dates
    THEBAY.LastDate = flowdf['ordinaldate'].iloc[-1]
    '''

    ############## TODO: Hopefully, converted... but, need to rename columns
    #flowdf.columns = ['msflow']
    #flowdf['jsflow'] = jsflow.iloc[:,0]
    #flowdf['mlflow'] = mlflow.iloc[:,0]
    #flowdf['ordinaldate'] = flowdf.index.to_series().apply(datetimeToOrdinal)

    '''
    # Convert Streamflow units (per EFDC_file_prep precedent)
    #       RHESSYS .daily output is mm/m^2 of basin
    #       generating cubic meters per second
    #       scale up slightly to reflect portion of watershed downstream of swanton gauge
    # Convert streamflow value to cubic meters per second
    #flowdf['msflow'] = flowdf['streamflow'] /1000 * 2199524000 / 3600 / 24
    # Scale to include ungaged proportion of the watershed not captured by the gage station at swanton
    #flowdf['msflow'] = flowdf['msflow'] / 0.98
    flowdf['msflow'] = flowdf['streamflow_adj_cms'] # Unit conversion and adjustment now done in basin.daily

    # get flow from Mill and JewittStevens also
    if jsflow.empty:
        logger.info('Jewett/Stevens Flow not read. Scaling from Missisquoi Flow.')
        flowdf['jsflow'] = flowdf['msflow'] * 0.038   # J/S about Flow of Rock, 0.038 of Missisquoi
    else:
        flowdf['jsflow'] = jsflow['streamflow_adj_cms']

    if mlflow.empty:
        logger.info('Mill Flow not read. Scaling from Missisquoi Flow.')
        flowdf['mlflow'] = flowdf['msflow'] * 0.012   # Mill about 1/3 of Rock
    else:
        flowdf['mlflow'] = mlflow['streamflow_adj_cms']
    '''

    # Store flow series in bay object for later
    THEBAY.flowdf = flowdf[['ordinaldate', 'msflow', 'mlflow', 'jsflow']].copy()

    logger.info('Daily Flow Data Scaled')
    # Scale Additional Inflows from Predicted Inflow
    #       sourcelist has list of source IDs
    #       sourcemap defines the source name and proportion of hydromodel output flow
    for baysource in THEBAY.sourcelist :
        logger.info('Generating Bay Source File for Id: '+baysource)
        bs_prop = THEBAY.sourcemap[baysource]['prop']       # proportion of input file for this source
        wshed = THEBAY.sourcemap[baysource]['wshed']        # get column name of watershed flow source for this stream
        flowdf[baysource]= flowdf[wshed] * bs_prop          # scale source from hydromodel flow (some are split)
        bs_name = THEBAY.sourcemap[baysource]['name']
        filename = bs_name + '_Flow.dat'
        logger.info('Bay Source File to Generate: '+filename)
        # Write Inflow File
        # open the file in output directory
        pathedfile = os.path.join(THEBAY.infile_dir, filename)
        with open(pathedfile, mode='w', newline='') as output_file:
            THEBAY.addfile(fname=filename)    # remember generated file names
            # output the header text
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('! Written by AEM3D_prep_IAM                           !\n')
            output_file.write('! Bay ID: '+ THEBAY.bayid + '                         !\n')
            output_file.write('! Bay Source: ' + bs_name + '                         !\n')
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('1 data sets\n')
            output_file.write('0 seconds between data\n')
            output_file.write('0    ' + baysource + '\n')
            output_file.write('TIME      INFLOW\n')
            # output the ordinal date and flow value time dataframe columns
            flowdf.to_csv(path_or_buf = output_file, columns= ['ordinaldate', baysource], float_format='%.3f',
            sep=' ', index=False, header=False)

##
#       End of Flow Data Import
#
##

#########################################################################
#
#   All Things Meteorology Related Handled Herein
#
##########################################################################


def genclimatefiles(whichbay):

    global SCENARIO

    THEBAY = whichbay   # passed object defining bay characteristics
    year=THEBAY.year   # starting year pulled from IAMBAY class object (lib.py)

    #
    #   Read WRF climate data
    #
    logger.info('Processing Meterological Data')

    ######################## Have to get new climate
    # TODO: rename columns if necessary
    #climateForecast = aggregate_df(dates, hours, loc_list, location_dataframes, stepType, typeOfLevel)
    #climateObsBTV = getBTV()
    #climateObsCR = getCol()

    '''
    climate = [
        cl.get_climate_data(
        SCENARIO,
        config['vars'],
        cl.HOURLY, [year],
        config['zones'])
        for config in THEBAY.climateZones
    ]
    # Unlist climate
    climate = {k : v for d in climate for k, v in d.items()}
    '''

    #
    #   Generate the Ordinal Date (yeardayofyear.percentseconds)
    #

    #wrfdf = pd.DataFrame(climate[0]['T2']['403'].index, columns = ['wrftime'])
    #wrftime = climate[0]['T2']['403'].index

    #   Cannot use DateTime DayOfYear function, as it generates 366 days in a Leap Year
    #wrfdf['yearday'] = wrfdf['wrftime'].dt.strftime('%Y%j')   # Concatenates year and day of year in one string

    #   day of year is array index/24, adjusted for indexing from 1, left padded by zeros
    # wrfdf['dayofyear'] = 1 + np.floor(wrfdf.index/24)   # groups of 24 hourly records with same D of Y
    # wrfdf['dayofyear'] = wrfdf['dayofyear'].astype(int)
    # wrfdf['yearday'] = str(THEBAY.year) + wrfdf['dayofyear'].apply(str).apply(lambda x: x.zfill(3))

    # wrfdf['totseconds']= wrfdf['wrftime'].dt.hour * 3600 + \
    #                      wrfdf['wrftime'].dt.minute * 60 + \
    #                      wrfdf['wrftime'].dt.second

    # wrfdf['fracsec'] = wrfdf['totseconds']/pd.Timedelta(days=1).total_seconds()  #Fraction of the day's seconds

    # wrfdf['ordinaldate'] = wrfdf['yearday'] + wrfdf['fracsec'].astype(str).str.slice(start=1,stop=6)  # add the percentage seconds since noon
    # wrfdf['ordinaldate'] = wrfdf['ordinaldate'].str.pad(width=12, side='right', fillchar='0')  # pretty right edge

    #wrfdf['ordinaldate'] = pandasDatetimeToOrdinal(climate[0]['T2']['403'].index)
    #ordinaldate = pd.Series(wrfdf['ordinaldate'].array, index = wrfdf['wrftime'])

    # Remove this... only used in two other spots... use function instead
    #ordinaldate = pandasDatetimeToOrdinal(climate['T2']['403'].index)
    #print(ordinaldate)

    ################################################################################
    #
    #   Moving average over 4 days of Air Temp (WRF T2)
    #
    #print('Whole dataset TEMP Shape')
    #print(wrf_data.variables['T2'].shape)
    # TODO New air_temp -- adjust window below if not hourly
    # air_temp = {"401": pd.concat([climateObsCR['T2'],climateForecast['401']['T2']),
    #             "402": pd.concat([climateObsCR['T2'],climateForecast['401']['T2']),
    #             "403": pd.concat([climateObsCR['T2'],climateForecast['401']['T2'])
    # }

    # air_temp = climate['T2']    # temp at 2m
    #wrfdf['wtr_temp'] = wrfdf['air_temp'].rolling(window=4,min_periods=1).mean() # moving average over 4 days
    # Use air temp at zone 403 (ILS)
    wtr_temp = air_temp['403'].rolling(window=96,min_periods=1).mean() # moving average over 4 days

    wtr_temp.loc[wtr_temp<0] = 0  # no subfreezing water
    wtr_temp = wtr_temp + 0.75    # 0.75 correction based on WQS Docs 2021.05.27

    # Store temp series in bay object for later use in wq calcs
    # THEBAY.tempdf = wrfdf[['ordinaldate', 'wtr_temp']].copy()
    THEBAY.tempdf = pd.DataFrame({'ordinaldate' : wtr_temp.index.to_series().apply(datetimeToOrdinal), 'wtr_temp' : wtr_temp.array})


    #
    #       Write the temperature file for each source of the bay
    #
    for baysource in THEBAY.sourcelist :

        bs_name = THEBAY.sourcemap[baysource]['name']
        filename = bs_name + '_Temp.dat'
        logger.info('Generating Bay Source Temperature File: '+filename)

        # Write Temp File
        # open the file in output directory
        pathedfile = os.path.join(THEBAY.infile_dir, filename)
        with open(pathedfile, mode='w', newline='') as output_file:

            THEBAY.addfile(fname=filename)        # remember generated bay files

            # output the header text
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('! Written by AEM3D_prep_IAM                           !\n')
            output_file.write('! Bay ID: '+ THEBAY.bayid + '                            !\n')
            output_file.write('! Bay Source: ' + bs_name + '                         !\n')
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('1 data sets\n')
            output_file.write('0 seconds between data\n')
            output_file.write('0    '+baysource+'\n')
            output_file.write('TIME      WTR_TEMP\n')

            # output the ordinal date and temp dataframe columns
            seriesIndexToOrdinalDate(wtr_temp).to_csv(path_or_buf = output_file, float_format='%.3f',
            sep=' ', index=True, header=False)

    #
    #   end water temp file

    ######################################################################################
    #
    #   Rain (WRF RAIN) and Snow
    #
    bay_rain = dict()
    bay_snow = dict()
    for zone in climate['RAIN'].keys():
        bay_rain[zone] = climate['RAIN'][zone] * 24 / 1000.0     # want daily cummulative rate in meters
        TEMP = pd.to_numeric(air_temp['403'])

        # Original Try: Use https://www.ncdc.noaa.gov/sites/default/files/attachments/Estimating_the_Water_Equivalent_of_Snow.pdf
        #   and fit a quadratic regression through it
        #   SNOW = e^(2.3129520 + -0.0962303*TEMP(C) + -0.0009388*TEMP(C)*TEMP(C)) * PRECIP
        #   This was too much though, so looking at the calibration input, divide in half
        #   to get closer to calibration input
        #snowcoeff = np.exp(2.3129520 - 0.0962303*TEMP - 0.0009388*TEMP*TEMP)/2.0
        
        # Take Two: Using GHCN data from BTV (Burlington Airport)
        # See Pat's rainToSnow.R
        # Calculate ratio using SNOW/PRCP, only take 1991 - present and days where SNOW and PRCP > 0
        # Tried using TMAX, but TAVG ((TMAX + TMIN)/2) got a curve closer to NASA table and makes more sense
        snowcoeff = np.exp(2.1413626 - 0.1921400*TEMP - 0.0079924*TEMP*TEMP)

        # Calculate snow from snowcoeff
        bay_snow[zone] = bay_rain[zone] * snowcoeff
        
        # If too warm, no snow - Use -2.0 C to get mean snowfall for 2017-2020 close to calibration data
        #  NOAA table above uses 34 F (1.1 C)
        bay_snow[zone].loc[air_temp['403'] > -2.0] = 0.0
        # In the WQS AEM3D calibration SNOW / RAIN data, RAIN looks like snow equivalent
        #   and SNOW is depth of snow, so we should have values for both when it snows
        # bay_rain[zone].loc[bay_snow[zone] > 0.0] = 0.0

        # print('Snow Stats for zone ', zone)
        # print(bay_snow[zone].describe(percentiles=[]))
        # print('Rain Stats for zone ', zone)
        # print(bay_rain[zone].describe(percentiles=[]))

    #
    # Write Precip Files for Bay
    #
    for zone in bay_rain.keys():
        filename = f'PRECIP_{zone}.dat'
        logger.info('Generating Bay Precipitation File: '+filename)

        with open(os.path.join(THEBAY.infile_dir, filename), mode='w', newline='') as output_file:

            THEBAY.addfile(fname=filename)        # remember generated bay files

            # output the header text
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('! Written by AEM3D_prep_IAM                           !\n')
            output_file.write('! Bay ID: '+ THEBAY.bayid + '                         !\n')
            #output_file.write('! Bay Source: ' + bs_name + '                        !\n')
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('2 data sets\n')
            output_file.write('0 seconds between data\n')
            output_file.write('0    0    0\n')
            output_file.write('TIME      RAIN     SNOW\n')

            # output the ordinal date and flow value time dataframe columns
            pd.concat([
                    seriesIndexToOrdinalDate(bay_rain[zone]),
                    seriesIndexToOrdinalDate(bay_snow[zone])],
                    axis=1).to_csv(
                    path_or_buf = output_file,
                    float_format='%.3f',
                    sep=' ',
                    index=True, header=False)

    #
    #   end precip file


    ###########################################################################################
    #
    #   Longwave Radiation Measure : CLOUDS, LW_RAD_IN, or LW_RAD_NET
    #

    if SCENARIO.gcm.startswith("bree."):
        #writeLongwaveRadiationNet(climate, THEBAY)
        writeLongwaveRadiationDownward(climate, THEBAY)
    elif SCENARIO.gcm == "era5":
        #writeCloudCover(climate, THEBAY)
        # Found LWDOWN in ERA5 -- strd
        writeLongwaveRadiationDownward(climate, THEBAY)



    ###########################################################################################
    #
    #   Wind Speed and Direction
    #           U10 and V10 variables converted to Earth compass direction and speed
    #

    # pjc 6/2 -- You are going to have to figure out a way to do this for multiple columns
    #            of U10 and V10.  I recommend setting this to a new DataFrame
    #            (without the iloc[] at the end) and then
    #            each column of that new dataframe is the windspeed series for each location
    #            !! Similar issue for winddir and rhum
    windspd = dict()
    winddir = dict()
    for zone in climate['V10'].keys():
        # a bit of vector math to combine the East(U) and North(V) wind components
        windspd[zone] = np.sqrt(
            np.square(climate['U10'][zone]) +
            np.square(climate['V10'][zone])
        )

        #  a bit of trig to map the wind vector components into a direction
        #  𝜙 =180+(180/𝜋)*atan2(𝑢,𝑣)
        winddir[zone] = 180 + np.arctan2(
                climate['U10'][zone],
                climate['V10'][zone]
            ) * 180 / np.pi

    # Write Wind Speed and Direction File
    #
    for zone in windspd.keys():
        filename = f'WS_WD_{zone}.dat'
        logger.info('Generating Wind Speed and Direction File: '+filename)

        with open(os.path.join(THEBAY.infile_dir, filename), mode='w', newline='') as output_file:

            THEBAY.addfile(fname=filename)        # remember generated bay files

            # output the header text
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('! Written by AEM3D_prep_IAM                           !\n')
            output_file.write('! Bay ID: '+ THEBAY.bayid + '                            !\n')
            #output_file.write('! Bay Source: ' + bs_name + '                         !\n')
            output_file.write('!-----------------------------------------------------!\n')
            output_file.write('2 data sets\n')
            output_file.write('0 seconds between data\n')
            output_file.write(f'0            {zone}          {zone}\n')
            output_file.write('TIME         WIND_SPEED	  WIND_DIR\n')

            # output the ordinal date and flow value time dataframe columns
            pd.concat([
                seriesIndexToOrdinalDate(windspd[zone]),
                seriesIndexToOrdinalDate(winddir[zone])],
                axis=1).to_csv(
                path_or_buf = output_file,
                float_format='%.3f',
                sep=' ',
                index=True, header=False)
    #
    #   end wind file



    ###########################################################################################
    #
    #   Relative Humidity : WRF, Calculated from Q2
    #
    #
    # rhum = dict()
    # for zone in climate[cs['RHUM']]['T2'].keys():
    #     # Using an equation ported from metuils.r that was authored by David LeBauer
    #     # es = 6.112 * np.exp((17.67 * t2) / (t2 + 243.5))
    #     # constants synched with https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html

    #     q2 = climate[cs['RHUM']]['Q2'][zone]
    #     psfc = climate[cs['RHUM']]['PSFC'][zone]
    #     t2 = climate[cs['RHUM']]['T2'][zone]

    #     es = 6.1078 * np.exp((17.2694 * t2) / (t2 + 237.30))    # Saturation vapor pressure over water
    #     e = q2 * psfc * 0.01 / (0.378 * q2 + 0.622)        # factor of 0.01 is converting Pascal to mbar
    #     rhum_temp = e/es
    #     rhum_temp[rhum_temp>1] = 1
    #     rhum_temp[rhum_temp<0] = 0
    #     rhum[zone] = rhum_temp
    rhum = climate['RH2']

    #   Write Relative Humidity File
    #
    for zone in rhum.keys():
        filename = f'RELHUM_{zone}.dat'
        logger.info('Generating Rel Humidity File: '+filename)
        writeFile(
            os.path.join(THEBAY.infile_dir, filename),
            THEBAY.bayid,
            zone,
            "REL_HUM",
            seriesIndexToOrdinalDate(rhum[zone]))
        THEBAY.addfile(fname=filename)


    #   Write Air Temp File
    #
    for zone in air_temp.keys():
        filename = f'AIRTEMP_{zone}.dat'
        logger.info('Generating Air Temp File: '+filename)
        writeFile(
            os.path.join(THEBAY.infile_dir, filename),
            THEBAY.bayid,
            zone,
            "AIR_TEMP",
            seriesIndexToOrdinalDate(air_temp[zone]))
        THEBAY.addfile(fname=filename)


    # Write ShortwaveRad File
    #
    swdown = climate['SWDOWN'] # * 0.875         # Scale Solar based on matching observed for 2018

    for zone in swdown.keys():
        filename = f'SOLAR_{zone}.dat'
        logger.info('Generating Short Wave Radiation File: '+filename)
        writeFile(
            os.path.join(THEBAY.infile_dir, filename),
            THEBAY.bayid,
            zone,
            "SOLAR_RAD",
            seriesIndexToOrdinalDate(swdown[zone]))
        THEBAY.addfile(fname=filename)


    ###
    #
    #   Generate Lake Level Regression from Temp and Flow Derived Terms
    #       Calculate lake level based on previous regression model results.
    #		Ported from original IAM/tools/EFDC_file_prep/EFDC_file_prep.R
    #       The regression model equation for lake level is:
    #       LakeLevel = 94.05887 +
    #           0.007910834(temperature) +
    #           7.034478e-05 (discharge) +
    #           0.003396492(discharge7days) +
    #           0.01173037(discharge30days) +
    #           0.0258206(discharge60days)
    #
    ###
    logger.info('Generating Lake Levels')
    flowdf = THEBAY.flowdf.copy()       # get flow dataframe from bay object

    # TODO: Which to use?  streamflow_cms or streamflow_adj_cms
    streamflow_unadj_cms = flowdf['msflow'] / 0.98
    flowdf['flowmean_07'] = streamflow_unadj_cms.rolling(window=7,min_periods=1).mean() # moving average over 7 days
    flowdf['flowmean_30'] = streamflow_unadj_cms.rolling(window=30,min_periods=1).mean() # moving average over 30 days
    flowdf['flowmean_60'] = streamflow_unadj_cms.rolling(window=60,min_periods=1).mean() # moving average over 60 days

    # print('Missisqoui Flow Entering Bay')
    # print(flowdf['msflow'].describe(percentiles=[]))
    # print(flowdf['flowmean_07'].describe(percentiles=[]))
    # print(flowdf['flowmean_30'].describe(percentiles=[]))
    # print(flowdf['flowmean_60'].describe(percentiles=[]))

    #  need a temperature series in Fahrenheit
    #		reconcile wrf hourly series with flow daily series
    #print('Resampling Hourly Temp')

    #wrftemp = wrfdf[['wrftime', 'air_temp']]
    #wrfdailyraw = wrftemp.set_index('wrftime').resample('D').mean()
    
    # drop any null leap days generated by resample and convert to F
    wrfdailyraw = air_temp['403'].resample('D').mean()
    wrfdailyF = 32 + 1.8 * wrfdailyraw.dropna(axis=0, inplace=False, how='any')
    TemperatureF = wrfdailyF.to_numpy()

    # print('Air Temp Raw Stats (C)')
    # print(air_temp['403'].describe(percentiles=[]))
    # print(wrfdailyraw.describe(percentiles=[]))
    # print('WrfDailyF shape', wrfdailyF.shape)
    # print(wrfdailyF)
    # print(TemperatureF.describe(percentiles=[]))

    # print('Calculating Lake Levels')
    flowdf['LakeLevel'] = 94.05887 + 0.007910834 * TemperatureF + \
        7.034478e-05 * flowdf['msflow'] + \
        0.003396492 * flowdf['flowmean_07'] + \
        0.01173037 * flowdf['flowmean_30'] + \
        0.0258206 * flowdf['flowmean_60']

    # print(flowdf['LakeLevel'].describe(percentiles=[]))

    # print('Calculating Lake Level Bias')
    # Next, apply the bias correction from the bias correction quadratic regression on the residuals against the observed lake level:
    bias_correction = -358.51020205 + \
        7.16150850 * flowdf['LakeLevel'] + \
        -0.03570562 * np.power(flowdf['LakeLevel'],2)

    flowdf['LakeLevel_corrected'] = flowdf['LakeLevel'] + bias_correction


    #	AEM3D lake input defined as meters above 93ft - do the math
    flowdf['LakeLevel_delta'] = (flowdf['LakeLevel_corrected'] - 93) * 0.3048

    #
    #   write out lake level file
    #
    print('Lake Level (m) above 93ft')
    print(flowdf['LakeLevel_delta'].describe(percentiles=[]))



    filename = 'Lake_Level.dat'
    logger.info('Writing Lake Level File '+filename)

    pathedfile = os.path.join(THEBAY.infile_dir, filename)
    with open(pathedfile, mode='w', newline='') as output_file:

        THEBAY.addfile(fname=filename)        # remember generated bay files

        # output the header text
        output_file.write(
            '!-----------------------------------------------------!\n')
        output_file.write(
            '! Written by AEM3D_prep_IAM                           !\n')
        output_file.write('! Bay ID: ' + THEBAY.bayid +
                          '                            !\n')
        output_file.write(
            '! values in (m) above 93 ft                           !\n')
        output_file.write(
            '!-----------------------------------------------------!\n')
        output_file.write('1 data sets\n')
        output_file.write('0 seconds between data\n')
        output_file.write('0    300\n')
        output_file.write('TIME	  HEIGHT\n')

        # output the ordinal date and flow value time dataframe columns
        flowdf.to_csv(path_or_buf=output_file, columns=['ordinaldate', 'LakeLevel_delta'], float_format='%.3f',
                      sep=' ', index=False, header=False)

        flowdf.to_csv(path_or_buf='lakeheight.csv', float_format='%.3f', sep=' ', index=False, header=True)

    ##
    #
    #   End of Lake Level From Temp and Flow Generation
    #
    #


def gensalinefile(whichbay):
    #
    #   Salinity data time series
    #
    global SCENARIO

    THEBAY = whichbay

    with open('TEMPLATES/salinity.template.txt', 'r') as file:
        template = Template(file.read())

    pathedfile = os.path.join(THEBAY.infile_dir, 'Inflows_Salinity.dat')
    with open(pathedfile, 'w') as output_file:

        # remember generated bay files
        THEBAY.addfile(fname='Inflows_Salinity.dat')

        output_file.write(template.substitute(**{
            'source_id_list': '  '.join(THEBAY.sourcelist),
            'firstdate': THEBAY.FirstDate,
            'lastdate': THEBAY.LastDate
        }))


def genboundaryfile(whichbay):
    #
    #   P boundary condition data time series
    #
    global SCENARIO

    THEBAY = whichbay

    with open('TEMPLATES/bc_p.template.txt', 'r') as file:
        template = Template(file.read())

    pathedfile = os.path.join(THEBAY.infile_dir, 'OpenBC_P.dat')
    with open(pathedfile, 'w') as output_file:

        # remember generated bay files
        THEBAY.addfile(fname='OpenBC_P.dat',ftype='boundary_condition_file')

        output_file.write(template.substitute(**{
            'firstdate': THEBAY.FirstDate,
            'lastdate': THEBAY.LastDate
        }))
        

def gentracerfiles(whichbay):
  #
  # Write Tracer Files - fixed series, use templates and set $year
  #
  THEBAY = whichbay

  tracerfiles = glob.glob('TEMPLATES/*racer*.template')
  # Loop through list of Phyto Templates, replacing $year with specific year
  for thisfile in tracerfiles:        # for each filename

    with open(thisfile,'r') as file:
      tracerStr = file.read()
    template = Template(tracerStr)

    #	pull name out of the path/filename.template and add dat extention
    thisfilename = os.path.splitext(os.path.split(thisfile)[1])[0]+'.dat'

    pathedfile = os.path.join(THEBAY.infile_dir, thisfilename)
    logger.info('Writing Tracer File, '+pathedfile)

    with open(pathedfile, 'w') as output_file:
        output_file.write(template.substitute(**{
            'nextyear': str(THEBAY.year+1),
            'year': str(THEBAY.year)
            }))

    # Change name of phyto templates to phyto dat
    if (thisfilename == 'tracer_release.dat') :
        THEBAY.addfile(fname=thisfilename,ftype='update_file')
    else :
        THEBAY.addfile(fname=thisfilename)
  #
  # End of Tracer File Generation
  #


def gencntlfile(whichbay):
    global SCENARIO

    THEBAY = whichbay

    with open('TEMPLATES/aem3dcntl.template.txt', 'r') as file:
        template = Template(file.read())

    #	control file is written to runtime directory
    pathedfile = os.path.join(THEBAY.run_dir, 'run_aem3d.dat')
    with open(pathedfile, 'w') as output_file:
        output_file.write(template.substitute(**{
            'start_date': THEBAY.FirstDate,
            # number of 300s steps in a 364 days (year minus 1 day, because 1st day is a start, not a step)
            # 27936 for 97 days (97*24*12)
            'iter_max': 27936
            }))

        # update the control file with all generated input files
        output_file.write(
            '! --------------------------------------------------------------- !\n')
        output_file.write(
            '! Input file names                                                 !\n')
        output_file.write('\''+os.path.split(THEBAY.infile_dir)[1] +
                        '\'                                     infile_dir\n')
        output_file.write(
            ' sparsedata_aem3d.unf                              3D_data_file\n')
        output_file.write(
            ' usedata_aem3d.unf                             preprocessor_file\n')
        output_file.write(
            ' datablock.xml                                     datablock_file\n')

        for f, t in THEBAY.bayfiles:        # for each filename, filetype
            output_file.write(' '+f+'              '+t+'\n')

        output_file.write(
                    '!  End                                                            !\n')


# end of control file generation

def AEM3D_prep_IAM(whichbay):

    THEBAY = whichbay
    #logger.info(f'Processing Bay: {THEBAY.bayid} for year {THEBAY.year}')

    # get flow files from hydrology model data
    getflowfiles(THEBAY)

    # generate climate files including lake levels
    genclimatefiles(THEBAY)

    # generate salinity file
    gensalinefile(THEBAY)

    # generate boundary condition file
    genboundaryfile(THEBAY)
    
    # generate tracer files
    gentracerfiles(THEBAY)

    # generate the water quality files (waterquality.py script)
    genwqfiles(whichbay=THEBAY)

    # generate control file
    gencntlfile(THEBAY)

    return 0
