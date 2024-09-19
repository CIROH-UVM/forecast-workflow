import datetime as dt
import pandas as pd

def datetimeToOrdinal(date):

	dayofyear = date.strftime('%j')

	# Left pad dayofyear to length 3 by zeros
	yearday = str(date.year) + dayofyear.zfill(3)

	totseconds = date.hour * 3600 + \
				 date.minute * 60 + \
				 date.second
	fracsec = totseconds / dt.timedelta(days=1).total_seconds()  #Fraction of the day's seconds

	ordinaldate = yearday + str(fracsec)[1:6].ljust(5,'0')  # add the percentage seconds since noon
	return ordinaldate

def ordinalToDatetime(ordinaldate):
    # Split the ordinal date into year, day of year, and fraction of the day
    year = int(ordinaldate[:4])  # First 4 characters are the year
    day_of_year = int(ordinaldate[4:7])  # Next 3 characters are the day of the year
    fraction_of_day = float('0' + ordinaldate[7:])  # Remaining part is the fraction of the day

    # Convert year and day of year into a date
    base_date = dt.datetime(year, 1, 1) + dt.timedelta(days=day_of_year - 1)

    # Calculate total seconds in the day and multiply by the fraction of the day
    total_seconds_in_day = 86400  # 24 * 3600 seconds in a day
    seconds_since_midnight = fraction_of_day * total_seconds_in_day

    # Add the seconds to the base_date
    final_datetime = base_date + dt.timedelta(seconds=seconds_since_midnight)

    return final_datetime

def index_to_ordinal_date(pandas_data):
	return pandas_data.rename(index=datetimeToOrdinal)
