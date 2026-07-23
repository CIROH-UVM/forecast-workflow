from datetime import datetime, timedelta
import os
import sh

def remove_directories(dirs):
	'''
	Recursively deletes each directory in the given list (via `rm -rf`).

	Args:
	-- dirs (list of str) [req]: absolute paths of the directories to delete.
	'''
	for dir in dirs:
		print(f'Deleting {dir}')
		sh.rm('-rf', dir)
		print(f'Successfully deleted {dir}')


def archive_directories(source_dirs, destination_dir):
	'''
	Copies each source directory into a destination directory (via `rsync -a`, preserving attributes).

	Args:
	-- source_dirs (list of str) [req]: absolute paths of the directories to archive.
	-- destination_dir (str) [req]: absolute path of the directory to archive them into.
	'''
	for source_dir in source_dirs:
		print(f'Archiving {source_dir}')
		sh.rsync('-a', source_dir, destination_dir)
		print(f'Successfully archived {source_dir} in {destination_dir}')

def archive_forecasts(source, destination, past_n = 10):
	'''
	Moves all forecast data except those from the past n days to a specified archive directory.

	Args:
	-- source (str) [req]: absolute path to the directory in which the forecast data are located.
	-- destination (str) [req]: absolute path to the directory in which the forecast data will be archived.
	-- past_n (int) [opt]: the number of past days to keep in the source directory. Current date included in count. Defaults to 10.
	'''
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


def archive_gfs(source, destination, past_n = 5):
	'''
	Moves all GFS data except those from the past n days to a specified archive directory.

	Args:
	-- source (str) [req]: absolute path to the directory in which the GFS data are located.
	-- destination (str) [req]: absolute path to the directory in which the GFS data will be archived.
	-- past_n (int) [opt]: the number of past days to keep in the source directory. Current date included in count. Defaults to 5.
	'''
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


def archive_nwm(source, destination, past_n = 5):
	'''
	Moves all NWM data except those from the past n days to a specified archive directory.

	Args:
	-- source (str) [req]: absolute path to the directory in which the NWM data are located.
	-- destination (str) [req]: absolute path to the directory in which the NWM data will be archived.
	-- past_n (int) [opt]: the number of past days to keep in the source directory. Current date included in count. Defaults to 5.
	'''
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
