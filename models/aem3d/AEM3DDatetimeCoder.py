import xarray
import datetime as dt
import numpy as np

class AEM3DDatetimeCoder(xarray.coders.CFDatetimeCoder):
	def add_timedelta_to_origin(self, delta, msec_resolution = False):
		### Tried using np's datetime64 functions, but timedelta64 requires an int as a first argument
		#return np.datetime64('-4713-01-01') + np.timedelta64(int(delta), 'D')
		# Days from Julian Start to Jan 1, 1 AD: 4714*365 + 816 Leap Days = 1721426
		return_value = dt.datetime(1,1,1,12,0) + dt.timedelta(days = delta-1721426)
		if not msec_resolution:
			return_value = (return_value + dt.timedelta(seconds = .5)).replace(microsecond=0)
		return return_value

	def __init__(self, use_cftime = None, time_unit = "ns", msec_resolution = False):
		self.msec_resolution = msec_resolution
		super().__init__(use_cftime, time_unit)

	def decode(self, variable, name = None):
		data = np.vectorize(self.add_timedelta_to_origin)(variable.data, self.msec_resolution)
		return xarray.Variable(dims='T', data=data, attrs={'long_name': 'Datetime'}, encoding={'dtype': 'np.datetime64'})
		
		### Other options I tried... but need to return an xarray Variable
		#return pd.to_datetime("01-01-1970") + pd.to_timedelta(variable, unit='seconds')
		#return super().decode(variable, name)	