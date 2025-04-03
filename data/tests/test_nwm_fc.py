import os
import pytest
import datetime as dt
from data.nwm_fc import prepForDownloads, download_nwm, process_nwm, get_data

TEST_DATA_DIR = "/gpfs2/scratch/nbeckage/test_data/"

# @pytest.mark.parametrize("reference_date, xreference_date, member, netcdf_template", [
# 	("20240101", dt.datetime(2024, 1, 1, tzinfo=dt.timezone.utc), "long_range_mem3", "nwm.t00z.long_range.channel_rt_3"),
# 	("2024010118", dt.datetime(2024, 1, 1, 18, tzinfo=dt.timezone.utc), "medium_range_mem5", "nwm.t18z.medium_range.channel_rt_5"),
# 	(dt.date(2024, 1, 1), dt.datetime(2024, 1, 1, tzinfo=dt.timezone.utc), "short_range", "nwm.t00z.short_range.channel_rt"),
# 	(dt.datetime(2024, 1, 1, 12), dt.datetime(2024, 1, 1, 12, tzinfo=dt.timezone.utc), "medium_range_mem1", "nwm.t12z.medium_range.channel_rt_1")
# ])
# def test_prepForDownloads(reference_date, xreference_date, member, netcdf_template):
# 	ref, template, _ = prepForDownloads(reference_date, member, TEST_DATA_DIR)
# 	assert ref == xreference_date
# 	assert template == netcdf_template

def test_download_nwm():
	download_nwm("2024010112", "long_range_mem1", hours=list(range(20,40)), gcs=True, download_dir=TEST_DATA_DIR)

def test_process_nwm():
	pass

def test_get_data():
	pass