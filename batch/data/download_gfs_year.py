import datetime as dt
from data.gfs_fc_thredds import download_gfs_threaded
from lib import generate_hours_list
import os

# define year for download
year = 2023
# start date for data grab - note that forecast cycle info is grabbed form the 'date' object in the active for loop
start_dt = dt.datetime(year,6,1)
end_dt = dt.datetime(year,10,31)
# number of forecast days to download - this is the amount of forecast data we want for all members
# from forecast_start to forecast_end
# member 7, for example, will grab 9 total days of data (1.5 for member timedelta adjustment, 7.5 from forecast_start to end)
fc_days = 7.5

# define yearly gfs dir
year_dir = f"/netfiles/ciroh/7dayHABsHindcast/hindcastData/gfs/{year}/"
delta = end_dt - start_dt

# download gfs data for members 1-7
for member in range(1, 8):
	# define the timedelta for the forecast cycle adjust
	cycle_delta = dt.timedelta(hours=6*(member-1))
	n_hours = int(fc_days*24 + 6*(member-1))
	# create hours list with cycle adjustment
	hours = generate_hours_list(num_hours=n_hours, source='gfs', archive=True)
	dates = [(start_dt - cycle_delta) + dt.timedelta(days=d) for d in range(delta.days+1)]
	print(f"Downloading GFS data from THREDDS for year: {year}")
	print(f"\t START DATE: {dates[0].strftime('%Y%m%d')}")
	print(f"\t END DATE: {dates[-1].strftime('%Y%m%d')}")
	print(f"\t FORECAST CYCLE: {dates[0].hour:02d}")
	print(f"\t FORECAST HOURS: {n_hours}")
	for date in dates:
		download_dir = os.path.join(year_dir, f"gfs.{date.strftime('%Y%m%d')}/{date.hour:02d}/atmos")
		print(download_dir)
		download_gfs_threaded(date=date, hours=hours, gfs_data_dir=download_dir, num_threads=8)
