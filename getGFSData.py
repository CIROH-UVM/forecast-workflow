# 1. Create loop to download all of the relevant .grib2 files (or whatever other file extension turns out to be inmportant)
# 2. Create a pandas dataframe that represents a specific lat/long (location), and each row is forecast data x hours out from the inital time
import os
import xarray as xr

ds = xr.open_dataset("../GFS_Data/08_18_2023/gfs.t00z.pgrb2.0p25.f/gfs.t00z.pgrb2.0p25.f000", engine="cfgrib", backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})
print(ds)

longnames = ["time","step","surface","valid_time"]
for v in ds:
    print("{}, {}, {}".format(v, ds[v].attrs["long_name"], ds[v].attrs["units"]))
    longnames.append("{}, {}".format(ds[v].attrs["long_name"], ds[v].attrs["units"]))

df = ds.to_dataframe()
df.columns=longnames

print(df)
print(df.info())
print(df.shape)