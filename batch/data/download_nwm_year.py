import datetime as dt
from data.nwm_fc import download_nwm_threaded
from lib import generate_hours_list
import os

# define year for download
year = 2023
# start date for data grab - note that forecast cycle info is grabbed form the 'date' object in the active for loop
start_dt = dt.datetime(2023,6,1)
end_dt = dt.datetime(2023,10,31)
# forecast type
fc_type = "medium_range_mem"
# days of forecast data to download (for generate_hours_list())
fc_days = 7.5
# define yearly gfs dir
year_dir = f"/netfiles/ciroh/7dayHABsHindcast/hindcastData/nwm/{year}/"
delta = end_dt - start_dt

# download nwm data for members 1-7
for member in range(1, 8):
	# hours list will be the same for NWM, regardless of member
	hours = generate_hours_list(num_hours=int(24*fc_days), source='nwm', forecast_type='medium')
	dates = [(start_dt + dt.timedelta(days=d)) for d in range(delta.days+1)]
	print(f"Downloading NWM {fc_type}{member} for year: {year}")
	print(f"\t START DATE: {dates[0].strftime('%Y%m%d')}")
	print(f"\t END DATE: {dates[-1].strftime('%Y%m%d')}")
	print(f"\t FORECAST CYCLE: {dates[0].hour:02d}")
	print(f"\t FORECAST HOURS: {24*fc_days}")
	for date in dates:
		forecast_cycle = f"{date.hour:02d}"
		forecast_type = fc_type.split('_')[0]
		netcdf_template = f"nwm.t{forecast_cycle}z.{forecast_type}_range.channel_rt_{member}"
		download_dir = os.path.join(year_dir, f"nwm.{date.strftime('%Y%m%d')}/{fc_type}{member}")
		print(download_dir)
		download_nwm_threaded(date=date, hours=hours, nwm_data_dir=download_dir, num_threads=8, fname_template=netcdf_template, use_google_bucket=True)
