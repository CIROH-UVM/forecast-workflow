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

def index_to_ordinal_date(pandas_data):
	return pandas_data.rename(index=datetimeToOrdinal)
