from datetime import datetime, timedelta
import os
import sh

def remove_directories(dirs):
	for dir in dirs:
		print(f'Deleting {dir}')
		sh.rm('-rf', dir)
		print(f'Successfully deleted {dir}')


def archive_directories(source_dirs, destination_dir):
	for source_dir in source_dirs:
		print(f'Archiving {source_dir}')
		sh.rsync('-a', source_dir, destination_dir)
		print(f'Successfully archived {source_dir} in {destination_dir}')

### Moves all forecast data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the directory in which the forecast data are located
# -- destination (str) [required]: absolute path to the directory in which the forecast data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source directory. current date included in count.
def archive_forecasts(source = "/data/forecastRuns/", destination = "/netfiles/ciroh/forecastArchive/", past_n = 10):
	print("Forecast archive process initiated")
	today = datetime.today()
	past5 = ['7dayforecast-'+(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]
	
	#### New file list method
	to_archive = [dir for dir in os.listdir(source) if dir not in past5]
	
	print(f'The following directories in {source} will be archived to {destination}:\n{to_archive}')
	archive_directories([os.path.join(source, directory) for directory in to_archive], destination)
	print('Forecast archive process complete')

	print(f'The following directories will be deleted:\n{to_archive}')
	remove_directories([os.path.join(source, directory) for directory in to_archive])
	print('Forecast directory deleting complete')


### Moves all gfs data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the directory in which the gfs data are located
# -- destination (str) [required]: absolute path to the directory in which the gfs data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source directory. current date included in count.
def archive_gfs(source, destination, past_n = 5):
	print("GFS archive process initiated")
	today = datetime.today()
	past5 = ['gfs.'+(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]

	#### New file list method
	to_archive = [dir for dir in os.listdir(source) if dir not in past5]

	print(f'The following directories in {source} will be archived to {destination}:\n{to_archive}')
	archive_directories([os.path.join(source, directory) for directory in to_archive], destination)
	print('GFS archive process complete')

	print(f'The following directories will be deleted:\n{to_archive}')
	remove_directories([os.path.join(source, directory) for directory in to_archive])
	print('GFS directory deleting complete')


### Moves all nwm data except those from the past n days to a specified archive directory
# -- source (str) [required]: absolute path to the directory in which the nwm data are located
# -- destination (str) [required]: absolute path to the directory in which the nwm data will be archived
# -- past_n (int) [optional]: the number of past days to keep in source directory. current date included in count.
def archive_nwm(source, destination, past_n = 5):
	print("NWM archive process initiated")
	today = datetime.today()
	past5 = [(today - timedelta(days=i)).strftime('%Y%m%d') for i in range(past_n)]

	#### New file list method
	to_archive = [dir for dir in os.listdir(source) if dir not in past5]

	print(f'The following directories in {source} will be archived to {destination}:\n{to_archive}')
	archive_directories([os.path.join(source, directory) for directory in to_archive], destination)
	print('NWM archive process complete')

	print(f'The following directories will be deleted:\n{to_archive}')
	remove_directories([os.path.join(source, directory) for directory in to_archive])
	print('NWM directory deleting complete')
