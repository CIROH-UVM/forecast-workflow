from datetime import datetime, timedelta
from glob import glob
import os
import sh

def remove_directories(dirs):
	for dir in dirs:
		print('Deleting {}'.format(dir))
		sh.rm('-rf', dir)
		print('Successfully deleted {}'.format(dir))


def archive_directories(source_dirs, destination_dir):
	for source_dir in source_dirs:
		print('Archiving {}'.format(source_dir))
		sh.rsync('-a', source_dir, destination_dir)
		print('Successfully archived {} in {}'.format(source_dir, destination_dir))


### Moves all forecast data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the folder in which the forecast data are located
# -- destination (str) [required]: absolute path to the folder in which the forecast data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source folder. current date included in count.
def archive_forecasts(source = "/data/forecastRuns/", destination = "/netfiles/ciroh/forecastArchive/", past_n = 10):
	print("Forecast archive process initiated")
	today = datetime.today()
	past5 = ['7dayforecast-'+(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]
	
	#### Old file list method
	# # navigate to the forecast source dir
	# os.chdir(source)
	# all_dates = glob('7dayforecast-*')
	# to_archive = [date for date in all_dates if date not in past5]
	#### New file list method
	to_archive = [dir for dir in os.listdir(source) if dir not in past5]
	
	#print('The following folders in {} will be archived to {}:\n{}'.format(source, destination, to_archive))
	print(f'The following folders in {source} will be archived to {destination}:\n{to_archive}')
	archive_directories([os.path.join(source, directory) for directory in to_archive], destination)
	# for folder in to_archive:
	# 	print('Archiving {}'.format(os.path.join(source, folder)))
	# 	#sh.rsync('-a', os.path.join(source, folder), destination)
	# 	print('Successfully archived {} in {}'.format(os.path.join(source, folder), destination))
	print('Forecast archive process complete')

	print('The following directories will be deleted:\n{}'.format(to_archive))
	remove_directories([os.path.join(source, directory) for directory in to_archive])
	print('Forecast directory deleting complete')


### Moves all gfs data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the folder in which the gfs data are located
# -- destination (str) [required]: absolute path to the folder in which the gfs data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source folder. current date included in count.
def archive_gfs(source, destination, past_n = 5):
	print("GFS archive process initiated")
	today = datetime.today()
	past5 = ['gfs.'+(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]
	
	#### Old file list method
	# # navigate to the gfs source dir
	# os.chdir(source)
	# # create glob to catch all gfs subfolders
	# all_dates = glob('gfs.*')
	# to_archive = [date for date in all_dates if date not in past5]
	#### New file list method
	to_archive = [dir for dir in os.listdir(source) if dir not in past5]

	print('The following folders in {} will be archived to {}:\n{}'.format(source, destination, to_archive))
	archive_directories([os.path.join(source, directory) for directory in to_archive], destination)
	# for folder in to_archive:
	# 	print('Archiving {}'.format(os.path.join(source, folder)))
	# 	#sh.rsync('-a', os.path.join(source, folder), destination)
	# 	print('Successfully archived {} in {}'.format(os.path.join(source, folder), destination))
	print('GFS archive process complete')

	print('The following directories will be deleted:\n{}'.format(to_archive))
	remove_directories([os.path.join(source, directory) for directory in to_archive])
	print('GFS directory deleting complete')


### Moves all nwm data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the folder in which the nwm data are located
# -- destination (str) [required]: absolute path to the folder in which the nwm data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source folder. current date included in count.
def archive_nwm(source, destination, past_n = 5):
	print("NWM archive process initiated")
	today = datetime.today()
	past5 = [(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]
	
	#### Old file list method
	# # navigate to the nwm source dir
	# os.chdir(source)
	# all_dates = glob(today.strftime('%Y')+'*')
	# to_archive = [date for date in all_dates if date not in past5]
	#### New file list method
	to_archive = [dir for dir in os.listdir(source) if dir not in past5]

	print('The following folders in {} will be archived to {}:\n{}'.format(source, destination, to_archive))
	archive_directories([os.path.join(source, directory) for directory in to_archive], destination)
	# for folder in to_archive:
	# 	print('Archiving {}'.format(os.path.join(source, folder)))
	# 	#sh.rsync('-a', os.path.join(source, folder), destination)
	# 	print('successfully archived {} in {}'.format(os.path.join(source, folder), destination))
	print('NWM archive process complete')

	print('The following directories will be deleted:\n{}'.format(to_archive))
	remove_directories([os.path.join(source, directory) for directory in to_archive])
	print('NWM directory deleting complete')
