from datetime import datetime, timedelta
from glob import glob
import os
import sh

### Moves all forecast data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the folder in which the forecast data are located
# -- destination (str) [required]: absolute path to the folder in which the forecast data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source folder. current date included in count.
def archive_forecasts(source = "/data/forecastRuns/", destination = "/netfiles/ciroh/forecastArchive/", past_n = 10):
	print("Forecast archive process initiated")
	today = datetime.today()
	past5 = ['7dayforecast-'+(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]
	# navigate to the forecast source dir
	os.chdir(source)
	all_dates = glob('7dayforecast-*')
	to_archive = [date for date in all_dates if date not in past5]
	print('the following folders in {} will be moved to {}:\n{}'.format(source, destination, to_archive))
	for folder in to_archive:
		print('archiving {}'.format(folder))
		sh.mv(folder, destination)
		print('successfully archived {} in {}'.format(folder, destination))
	print('Forecast archive process complete')
	return

### Moves all gfs data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the folder in which the gfs data are located
# -- destination (str) [required]: absolute path to the folder in which the gfs data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source folder. current date included in count.
def archive_gfs(source, destination, past_n = 5):
	print("GFS archive process initiated")
	today = datetime.today()
	past5 = ['gfs.'+(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]
	# navigate to the gfs source dir
	os.chdir(source)
	# create glob to catch all gfs subfolders
	all_dates = glob('gfs.*')
	to_archive = [date for date in all_dates if date not in past5]
	print('the following folders in {} will be moved to {}:\n{}'.format(source, destination, to_archive))
	for folder in to_archive:
		print('archiving {}'.format(folder))
		sh.mv(folder, destination)
		print('successfully archived {} in {}'.format(folder, destination))
	print('GFS archive process complete')
	return

### Moves all nwm data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the folder in which the nwm data are located
# -- destination (str) [required]: absolute path to the folder in which the nwm data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source folder. current date included in count.
def archive_nwm(source, destination, past_n = 5):
	print("NWM archive process initiated")
	today = datetime.today()
	past5 = [(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]
	# navigate to the nwm source dir
	os.chdir(source)
	all_dates = glob(today.strftime('%Y')+'*')
	to_archive = [date for date in all_dates if date not in past5]
	print('the following folders in {} will be moved to {}:\n{}'.format(source, destination, to_archive))
	for folder in to_archive:
		print('archiving {}'.format(folder))
		sh.mv(folder, destination)
		print('successfully archived {} in {}'.format(folder, destination))
	print('NWM archive process complete')
	return