import os
import sh
import sys
import logging, logging.config
import concurrent.futures
from contextlib import contextmanager
import datetime as dt
from string import Template
import inspect


class IAMBAY:
	'''
	IAMBAY Class contains the most common information used to describe a lake bay.

		bayid - Which Bay
		sourcelist - list of flow source IDs
		hydromodel - which hydrology model feeds the bay
	'''

	def __init__(self, bayid='ILS', year=2000, **kwargs):
		'''
		Init.
		'''
		self.bayid = bayid
		self.year = year
		self.FirstDate = str(year)+'001.000'  # default first date of data
		self.LastDate =  str(year)+'365.000'  # default last date of data
		self.wrfgridxy = (1,1)      # a default grid selection
		self.infile_dir = 'AEM3D-inputs/infiles'            # directory for boundary condition inputs
		self.template_dir = 'AEM3D-inputs/TEMPLATES'            # directory for boundary condition inputs
		self.run_dir = 'AEM3D-inputs'     # directory to contain everything the lake model needs to run
		self.bayfiles = []                           # initial (empty) list of boundary condition files for bay modeling
		self.flowdf =  None       # eventually to contain bay input flow time series dataframe
		self.tempdf = None        # will contain wtr_temp time series dataframe



		# list of water source IDs for multiple bays
		sourcelist_MAP = {
			'MB': {'201', '202', '203', '204', '21', '22'},
			'STA': {'17', '19'},
			'ILS': {'201', '202', '203', '204', '21', '22', '17', '19' }
		}

		# which WRF grid cell corresponds to lake data
		wrfgrid_MAP = {
			'MB': (32,45),
			'STA': (31,38),
			'ILS': (31,38)
		}
		
		climateZones_MAP = {
			'MB': [{}],
			'STA': [{}],
			'ILS': [{'vars': ['T2', 'AEMLW', 'SWDOWN', 'U10', 'V10'],
					 'zones': [{
						'name': '401',
						'desc': 'MB',
						'lat': 45.035531,
						'lon': -73.127665},
						{
						'name': '402',
						'desc': 'SAB',
						'lat': 44.790926,
						'lon': -73.155695},
						{
						'name': '403',
						'desc': 'IS',
						'lat': 44.778645,
						'lon': -73.17273}]},
					{'vars': ['RH2', 'RAIN'],
					 'zones': [{
						 'name': '0',
						 'desc': 'Unified',
						 'lat': 44.778645,
						 'lon': -73.17273}]}
					 ]}
					 
		#
		#  bay source list structure
		#      map of bay names
		#           bayname: map of source IDs to watershed and the flow proportion of watershed for that sourceID
		#               Watershed designations:
		#                   MS - Missisquoi (rhessys or swat)
		#                   ML - Mill (swat)
		#                   JS - Jewett Stevens (swat)
		#					prop: proportion of the associated watershed value to use for the source's flow
		#					adjust: adjust to flow as defined when bay sim model was calibrated (1.0 for no adjust)

		sourcelist_STRUCT = {
			'MB': {
				'201': {'name': 'MissisquoiRiverDC', 'wshed': 'msflow', 'prop': 0.57, 'adjust': 1.0},      # proportions of Miss River Flow
				'202': {'name': 'MissisquoiRiverNE', 'wshed': 'msflow', 'prop': 0.15, 'adjust': 1.0},
				'203': {'name': 'MissisquoiRiverCE', 'wshed': 'msflow', 'prop': 0.14, 'adjust': 1.0},
				'204': {'name': 'MissisquoiRiverNW', 'wshed': 'msflow', 'prop': 0.14, 'adjust': 1.0},
				'21':  {'name': 'RockRiver',         'wshed': 'msflow', 'prop': 0.038, 'adjust': 2.17},     # .03 of Total Miss Bay Inflow
				'22':  {'name': 'PikeRiver',         'wshed': 'msflow', 'prop': 0.228, 'adjust': 1.23}      # .18 of Total Miss Bay Inflow
			},
			'STA': {
				'17': {'name': 'MillRiver', 'wshed': 'mlflow', 'prop': 1.0, 'adjust': 1.04},                # Mill River flow
				'19': {'name': 'JewettStevens', 'wshed': 'jsflow', 'prop': 1.0, 'adjust': 4.28}             # JewettStevens flow
			},
			'ILS': {
				'201': {'name': 'MissisquoiRiverDC', 'wshed': 'msflow', 'prop': 0.57, 'adjust': 1.0},      # proportions of Miss River Flow
				'202': {'name': 'MissisquoiRiverNE', 'wshed': 'msflow', 'prop': 0.15, 'adjust': 1.0},
				'203': {'name': 'MissisquoiRiverCE', 'wshed': 'msflow', 'prop': 0.14, 'adjust': 1.0},
				'204': {'name': 'MissisquoiRiverNW', 'wshed': 'msflow', 'prop': 0.14, 'adjust': 1.0},
				'21':  {'name': 'RockRiver',         'wshed': 'msflow', 'prop': 0.038, 'adjust': 2.17},    # .03 of Total Miss Bay Inflow
				'22':  {'name': 'PikeRiver',         'wshed': 'msflow', 'prop': 0.228, 'adjust': 1.23},    # .18 of Total Miss Bay Inflow
				'17':  {'name': 'MillRiver',         'wshed': 'mlflow', 'prop': 1.0, 'adjust': 1.04},      #  about 1/3 volume of Rock
				'19':  {'name': 'JewettStevens',     'wshed': 'jsflow', 'prop': 1.0, 'adjust': 4.28}       #  about the volume of Rock
			}
		}



		# hydrology model specifications by bay ID
		hydromodel_MAP = {
			'MB': 'rhessys', # rhessys model
			'STA': 'swat',  # swat model
			'ILS': 'swat' # swat and rhessys models needed
		}


		if(bayid == 'MB' or bayid == 'STA' or bayid == 'ILS'):
			self.sourcelist = sourcelist_MAP.get(bayid)     # map of flow source IDs for specified bay
			self.sourcemap = sourcelist_STRUCT.get(bayid)   # structure giving name and proportion for sources
			self.hydromodel = hydromodel_MAP.get(bayid)     # the hydrology model associated with bay sources
			self.wrfgridxy = wrfgrid_MAP.get(bayid)         # which wrf climate grid coordinates
			self.climateZones = climateZones_MAP.get(bayid) # climate zones for weather inputs
		else:
			raise Exception('Bay ID is neither "MB", "STA", nor "ILS" ')


		self.validate()

	def validate(self):
		if True == False:
			raise Exception('The Bay is Busted')

	#
	#   method to keep track of files used to model the bay
	#       fname : unpathed name of file
	#       ftype : descriptor of file type (default = "boundary_condition_file")
	def addfile(self, fname, ftype = 'boundary_condition_file'):

		self.bayfiles.append((fname, ftype))
	##
	#       End of IAMBAY Class
	##


def generate_file_from_template(templateFile, outFile, theBay, sub_dict, outFileType='boundary_condition_file'):
	
	logger.info(f'Writing out {outFile} from template file {templateFile}')

	with open(os.path.join(theBay.template_dir, templateFile), 'r') as file:
		template = Template(file.read())

	with open(os.path.join(theBay.infile_dir, outFile), 'w') as output_file:
		output_file.write(template.substitute(**sub_dict))
	
	# remember generated bay files
	theBay.addfile(fname=outFile, ftype=outFileType)


@contextmanager
def cd(path):
	'''
	A context manager for changing directories.

	Will change the current working directory to path, and then once the
	context exits, change back to the original working directory, optionally
	re-raising an exception that was generated along the way.
	'''
	origin = os.getcwd()
	os.chdir(path)

	try:
		yield
	except Exception as e:
		os.chdir(origin)
		raise e

	os.chdir(origin)

class IAMLogger(logging.Logger):

	logging_initialized = False

	@classmethod
	def setup_logging(cls, name=None, default_level=logging.INFO):
		log_name = 'workers.lib'
		calling_package = get_calling_package()
		if(calling_package is not None):
			log_name = calling_package
		logging_config = {
			"version": 1,
			"disable_existing_loggers": False,
			"formatters": {
				"simple": {
					"format": "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
				}
			},
			"handlers": {
				"console": {
					"class": "logging.StreamHandler",
					"level": "DEBUG",
					"formatter": "simple",
					"stream": "ext://sys.stdout"
				},
				"info_file_handler": {
					"class": "logging.handlers.RotatingFileHandler",
					"level": "INFO",
					"formatter": "simple",
					"filename": f"{log_name}.log",
					"maxBytes": 10485760,
					"backupCount": 20,
					"encoding": "utf8"
				}
			},
			"root": {
				"level": "INFO",
				"handlers": ["console", "info_file_handler"]
			}
		}
		logging.config.dictConfig(logging_config)

		logger = logging.getLogger(log_name)
		logger.info('IAMlogging initialized')
		cls.logging_initialized = True

	@classmethod
	def getLogger(cls, name):
		if(not cls.logging_initialized):
			cls.setup_logging(name)

		return logging.getLogger(name)


class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, level):
       self.logger = logger
       self.level = level
       self.linebuf = ''

    def write(self, buf):
       for line in buf.rstrip().splitlines():
          self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass

def check_frame(frame):
	if(frame.f_globals['__package__'] is None):
		return False
	if(frame.f_globals['__package__'].startswith('importlib') or
	   frame.f_globals['__package__'] == 'workers'):
		return True
	return False

# New get_calling_package()
def get_calling_package():
	"""
	Returns the package name of the original caller in the call stack
	"""
	stack = inspect.stack()
	# iterating through the stack in reverse order to check the bottom first and work back up
	for frame_info in reversed(stack):
		# I don't know if this is a bullet-proof check, but it seems to work in my tests
		if frame_info.code_context is not None:
			module = os.path.splitext(os.path.basename(frame_info.filename))[0]
			return module
		
def report_stack():
	"""
	Helper function to provide execution frame information for the current call stack.
	"""
	stack = inspect.stack()
	print('Frame order: first is top of stack, last is bottom')
	for i, frame_info in enumerate(stack):
		print(f'For frame number {i+1}')
		print(f'\tframe: {frame_info.frame}')
		print(f'\tframe filename: {frame_info.filename}')
		print(f'\tframe function: {frame_info.function}')
		print(f'\tframe code context: {frame_info.code_context}')

def parse_to_datetime(date):
	"""
	Utility function to convert date or string to datetime obj, in UTC time

	Args:
	-- date (str, date, or datetime) [req]: the date to be parsed

	Returns:
	a datetime object representing the passed date
	"""
	# if already a datetime obj, just return the same obj; be sure to convert to UTC
	if isinstance(date, dt.datetime):
		return date.replace(tzinfo=dt.timezone.utc)
	# if a date obj, convert to datetime
	elif isinstance(date, dt.date):
		# return as a datetime, time set to midnight by default
		return dt.datetime.combine(date, dt.datetime.min.time()).replace(tzinfo=dt.timezone.utc)
	elif isinstance(date, str):
		if date == "":
			return None
		try:
			date = dt.datetime.strptime(date, "%Y%m%d").replace(tzinfo=dt.timezone.utc)
			return date
		except ValueError as ve:
			raise ValueError(f'{ve}. Please enter date strings in the format "YYYYMMDD"')
		except Exception as e:
			raise e
		
def get_hour_diff(start_date, end_date):
	"""
	Utility function to get the number of hours 
	"""
	time_difference = end_date - start_date
	hours = time_difference.total_seconds() / 3600
	return int(hours)

def generate_hours_list(num_hours, source, forecast_type='medium', archive=False):
	"""
	Creates a list of forecast hours to be downloaded

	Args:
	-- num_hours (int) [req]: how many hours of forecast data you want.
	-- saource (str) [req]: string indicating the forecast source. Valid values are 'gfs', 'nwm', and 'archive'.
	-- forecast_type (str) [opt]: The type of NWM forecast to grab. Valid values are 'short', 'medium', or 'long'.

	Returns:
	A list of hours in the format needed for the specifed forecast source. Or None if no valid source value is passed.
	"""
	match source:
		case 'gfs':
			if archive:
				return [f"{hour:03}" for hour in range(0, num_hours, 3)]
			if num_hours <= 120:
				return [f"{hour:03}" for hour in range(0, num_hours)]
			else:
				return [f"{hour:03}" for hour in range(0, 120)] + [f"{hour:03}" for hour in range(120, num_hours, 3)]
		case 'nwm':
			if forecast_type=='short' or forecast_type=='medium':
				return [f"{hour:03}" for hour in range(1, num_hours)]
			if forecast_type=='long':
				return [f"{hour:03}" for hour in range(6, num_hours, 6)]
		case 'buckets':
			if forecast_type=='short':
				return [f"{hour:03}" for hour in range(1, num_hours)]
			if forecast_type=='medium':
				if archive:
					return [f"{hour:03}" for hour in range(3, num_hours, 3)]
				else: return [f"{hour:03}" for hour in range(1, num_hours)]
			if forecast_type=='long':
				return [f"{hour:03}" for hour in range(6, num_hours, 6)]


		
def generate_date_strings(start_date, end_date):
	"""
	Creates a list of date strings in the format "YYYMMDD" for the range provided by start_date and end_date

	Args:
	-- start_date (obj) [req]: First date in the date range. Should be a string in the format "YYYYMMDD", a dt.date, or dt.datetime object
	-- end_date (obj) [req]: Last date in the date range. Should be a string in the format "YYYYMMDD", a dt.date, or dt.datetime object

	Returns:
	A list of date strings in the format "YYYMMDD" 
	"""
	date_strings = []

	# if not a datetime object, convert start_date to one
	start_date = parse_to_datetime(start_date)
	end_date = parse_to_datetime(end_date)

	for i in range(((end_date-start_date).days)+1):
		delta = dt.timedelta(days=i)
		date_strings.append((start_date+delta).strftime('%Y%m%d'))
	
	return date_strings

def download_data(url, filepath, use_google_bucket=False):
	"""
	Downloads, via curl, a file from a url to a specified filepath. Prints download progress.

	Args:
	-- url (str) [required]: full url to the file that is to be downloaded
	-- filepath (str) [required]: filepath to the download directory. should include desired name for file to be downloaded (ex /data/files/somedata.csv)

	"""
	# need to add curl output to list the old-fashioned way (no list comp) in case curl download is interupted/fails
	progress = []
	# '-#' turns on progress bar, '-o' writes download to filepath
	# _iter = '_err' allows us to iterate over the stderr produced by the call (for some reason progress bar is directed to stderr, not stdout)
	# _err_bufsize = 0 returns stderr one character at a time rather than as one very long string
	if use_google_bucket:
		print(f"Downloading file from URL: {url} using gsutil")
		call = sh.gsutil(['-m', 'cp', url, filepath], _iter = True, _out_bufsize = 100)
	else:
		print(f"Downloading file from URL: {url} using curl")
		call = sh.curl(['-#','-o', filepath, url], _iter = 'err', _err_bufsize = 0)
	try:
		for line in call:
			progress.append(line)
	except Exception as e:
		print("Download failed:")
		# if downloads fail, stop the run
		raise e
	else:
		# join strings, split by carriage return
		progress = ''.join(progress).split('\r')
		# log the last progress bar update (should be 100% ir download was finished, otherwise will show the point at which download ceased)
		print(progress[-1])
		print(f"Download Complete: {filepath}\n")

def multithreaded_download(download_list, n_threads):
	"""
	Implements download_data() with mulithreading to speed up downloads

	Args:
	-- download_list (list) [required]: a list containing tuples of (url_to_download, download_destination_path)
	-- num_threads (int) [opt]: number of threads to use.
	"""
	urls, paths, use_google_buckets = zip(*download_list)
	print(f'Using {n_threads} threads for downloads')
	with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as executor:
		executor.map(download_data, urls, paths, use_google_buckets)

def multithreaded_loading(load_func, file_list, n_threads):
	"""
	The idea here is to pass the function that reads datasets (xr.open_dataset, cfgrib) with multithreading, 
	and return a dictionary with the structure {fname : dataset}

	Args:
	-- load_func (function) [required]: the function to use to open the datasets
	-- file_list (list) [required]: list of file names to load
	-- n_threads (int) [req]: number of threads to use for mulithreading

	Returns:
	a dictionary where every key is a filename in file_list, and each corresponding value is the datastruct loaded for that file
	"""
	with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as executor:
		datasets = executor.map(load_func, file_list)
	# returns datasets in the order in which they were submited; order of file_list will be preserved
	# dataset_list = [d for d in datasets]
	dataset_list = []
	# try to capture errors from failed file loads, but continue on since a missing file or too won't hurt AEM3D much
	for d in datasets:
		try:
			dataset_list.append(d)
		except Exception as e:
			print(e)
	return dict(zip(file_list, dataset_list))

IAMLogger.setup_logging()
# initializing the logger; Logger name will be the calling package, IOW the py file that is originally called from the command line.
# i.e. if "python -m models.aem3d.AEM3D_prep_worker" is called, logger will be named "AEM3D_prep_worker"
logger = logging.getLogger(get_calling_package())
# implementing StreamToLogger will forward standard output (such as from print statements) to the logger
sys.stdout = StreamToLogger(logger, logging.INFO)