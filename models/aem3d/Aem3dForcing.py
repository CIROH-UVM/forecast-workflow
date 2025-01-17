
from data import (usgs_ob,
				  caflow_ob,
				  utils)
import datetime as dt
from get_args import get_args
import os
import pandas as pd

# I don't intend for the parent class, Aem3dForcings, to really be that functional;
# instead, it just has common attributes and methods that the child classes can inherit
class Aem3dForcings:
	def __init__(self,
				data=None,
				start_date=None,
				end_date=None,
				locations=None,
				variables=None,
				source=None,
				period=None,
				dir=None):
		
		# most basic instance of the class should reflect the basic parameters of a data module call
		self.data = data

		# ensure that the dates are datetime objects
		self.start_date = utils.parse_to_datetime(start_date)
		self.end_date = utils.parse_to_datetime(end_date)

		# set the rest of the attributes
		self.locations = locations
		self.variables = variables
		self.source = source
		self.period = period
		self.dir = dir

	def __str__(self):
		outstr = "".join([f"{self.__class__}\n",
						  f"Data: {self.data}\n",
						  f"Start Date: {self.start_date}\n",
						  f"End Date: {self.end_date}\n",
						  f"Locations: {self.locations}\n",
						  f"Variables: {self.variables}\n",
						  f"Source: {self.source}\n",
						  f"Period: {self.period}\n",
						  f"Directory: {self.dir}\n"])
		return outstr

	### Copy method
	def copy(self):
		return Aem3dForcings(self.data,
							self.start_date,
							self.end_date,
							self.locations,
							self.variables,
							self.source,
							self.period,
							self.dir)

	### Getter methods
	# calling this method access data because get_data is reserved for the data modules
	def access_data(self):
		return self.data
	
	### Getters for dates
	def get_start_date(self):
		return self.start_date
	
	def get_end_date(self):
		return self.end_date

	def get_dates(self):
		return self.start_date, self.end_date
	
	# rest of the getters
	def get_locations(self):
		return self.locations
	
	def get_variables(self):
		return self.variables	
	
	def get_source(self):
		return self.source

	def get_period(self):
		return self.period
	
	def get_dir(self):
		return self.dir
	
	### Setter methods
	def set_data(self, data):
		self.data = data

	def set_dates(self, start_date=None, end_date=None):
		if not isinstance(start_date, dt.datetime) and start_date is not None:
			raise TypeError(f"Start date must be a datetime object: {start_date}")
		if not isinstance(end_date, dt.datetime) and end_date is not None:
			raise TypeError(f"End date must be a datetime object: {end_date}")
		if start_date is not None and isinstance(start_date, dt.datetime):
			self.start_date = start_date
		if end_date is not None and isinstance(end_date, dt.datetime):
			self.end_date = end_date

	def set_locations(self, locations):
		self.locations = locations

	def set_variables(self, variables):
		self.variables = variables

	def set_source(self, source):
		self.source = source

	def set_period(self, period):
		if period not in ["spinup", "forecast"]:
			raise ValueError("Period must be 'spinup' or 'forecast'")
		self.period = period

	def set_dir(self, dir):
		self.dir = dir


class HydroForcings(Aem3dForcings):
	##### Hardcoding paramters relevant to AEM3D workflow #####
	global us_gauges
	us_gauges = {"MS":'04294000',
				 "JS":'04292810',
				 "ML":'04292750'}
	
	global ca_gauges
	ca_gauges = {"PK":'030424',
				 "RK":'030425'}
	
	global reaches
	reaches = {"US":us_gauges,
			   "CA":ca_gauges}

	streamflow = {'streamflow':'00060'}

	global reachnames
	reachnames = {"MS":"Missisquoi",
				  "ML":"Mill",
				  "JS":"Jewett",
				  "PK":"Pike",
				  "RK":"Rock"}
	
	### Initial Hydroings
	def __init__(self,
				data=None,
				start_date=None,
				end_date=None,
				locations=reaches,
				variables=streamflow,
				source=None,
				period=None,
				dir=None,
				service="iv"):
		# Call the parent class constructor
		super().__init__(data, start_date, end_date, locations, variables, source, period, dir)
		self.service = service
	
	def copy(self):
		copy_instance = super().copy()
		copy_instance.service = self.service
		return copy_instance
	
	def __str__(self):
		return super().__str__() + f"Service: {self.service}\n"

	# get method for service, which is unique to Hydroings
	def get_service(self):
		return self.service

	def set_service(self, service):
		if service not in ["iv", "dv"]:
			raise ValueError("Service must be 'iv' (instantaneous values) or 'dv' (daily values)")
		self.service = service
	
	# the function to call to get streamflow, wether from USGS, NWM, or CSV's
	def get_hydro_data(self):
		print("ACQUIRING HYDROLOGY DATA")
		print(f"\tPERIOD: {self.period}")
		print(f"\tSTART DATE: {self.start_date}")
		print(f"\tEND DATE: {self.end_date}")
		print(f"\tLOCATIONS: {self.locations}")
		print(f"\tVARIABLES: {self.variables}")
		print(f"\tSOURCE: {self.source.split(':')[0]}")
		if self.start_date is None or self.end_date is None:
			raise ValueError("Start and end date must be provided to get hydrology data")
		if self.locations is None:
			raise ValueError("Locations must be provided to get hydrology data")

		if self.source == "USGS_IV":
			usgs_data = usgs_ob.get_data(start_date=self.start_date,
										end_date=self.end_date,
										locations=self.locations["US"],
										variables=self.variables,
										service=self.service)
			# convert usgs streamflow to m³/s
			print(f"Converting USGS streamflow to m³/s...")
			for loc in usgs_data.keys():
				usgs_data[loc]['streamflow'] = usgs_data[loc]['streamflow'].rename('streamflow (m³/s)') * 0.0283168
			ca_data = caflow_ob.get_data(start_date=self.start_date,
										end_date=self.end_date,
										locations=self.locations["CA"])
			
			# now combine the dictionaries and set as data
			q = usgs_data | ca_data
			self.set_data(q)

			# This function is not ready to be implemented yet
			# if self.period == "spinup":
			# 	self.backfillCaFlowsSpinup()
		elif self.source.startswith("READ_CSV"):
			csv_dir = self.source.split(":")[-1]
			data_dir = self.dir
			self.set_dir(os.path.join(data_dir, csv_dir))
			print(f"\tHYDRO CSV DIR: {self.dir}")
			self.parse_hydro_csvs()

	def parse_hydro_csvs(self):
		'''
		Loads streamflow data from CSV files and sets resulting nested dict of Pandas Series as the data attribute of the HydroForcing object.
		CSV files should be stored in the directory specified by the "data_dir" setting in default_Settings.json. Naming and format
		for the streamflow CSVs follow the conventions establishded by Peter Isles' csv files. This and related functions could be made more general
		in the future, in part by changing the naming conventiond in the code below
		'''
		csv_q = {}
		for country, loc_dict in self.locations.items():
			# print(country)
			for loc in loc_dict.keys():
				fname = f"Q_{reachnames[loc]}"
				# print(fname)
				q = pd.read_csv(os.path.join(self.dir, f"{fname}.csv"), index_col='Date', parse_dates=True)
				new_index = q.index.map(lambda x: x.replace(hour=14))
				q.index = new_index.tz_localize('UTC')
				csv_q[loc] = {"streamflow":q[fname].rename("Streamflow (m³/s)")}
		self.set_data(csv_q)

if __name__ == "__main__":

	SETTINGS = get_args()
	# Test the class
	hydro = HydroForcings(period='spinup',
						  dir=SETTINGS["data_dir"],
						  source=SETTINGS['hydrology_dataset_spinup'])
	
	start = dt.datetime(2024, 1, 1)
	end = dt.datetime(2024, 9, 1)
	hydro.set_dates(start_date=start, end_date=end)
	hydro.get_hydro_data()

	print(hydro.access_data())