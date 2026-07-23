---
icon: lucide/cloud
---

# NOAA GFS Forecast (THREDDS)

Downloads and processes NOAA's GFS 0.25° forecast for a fixed Lake Champlain bounding box, using server-side NetCDF Subset Service (NCSS) requests against [UCAR RDA's THREDDS server](https://thredds.rda.ucar.edu/thredds/catalog/files/g/ds084.1/catalog.html). This is the recommended module for GFS acquisition — it supersedes [`gfs_fc`](gfs_fc.md), which is deprecated.

## `get_data()`

Downloads and processes GFS forecast data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `forecast_datetime` — the forecast launch date and cycle (00/06/12/18 UTC).
- `end_datetime` — the end date/time of the forecast window.
- `locations` — a dictionary mapping station name to a `(lat, lon)` tuple.
- `data_dir` — directory to store downloaded data. Defaults to the OS temp directory.
- `dnwld_threads` / `load_threads` — thread counts for downloading and reading data, respectively.
- `useTCDCInstant` — whether to use the instantaneous total-cloud-cover variable instead of the default 3/6-hour rolling average. Defaults to `False`.

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names (`T2`, `TCDC`, `U10`, `V10`, `RH2`, `RAIN`, `CPOFP`, `SWDOWN`), with values as Pandas Series.

## `download_gfs_threaded()`

Builds NCSS request URLs for each forecast hour (handling the special-cased 3/6-hour-average variables like cloud cover and shortwave radiation) and downloads them, throttled to stay within rate limits.

**Parameters:** `date`, `hours`, `gfs_data_dir`, `variables`, `num_threads`.

## `process_gfs_data()`

Loads the downloaded NetCDF files (skipping any truncated/bad files), isolates the requested locations (reusing [`isolate_loc_rows()`/`remap_longs()`](gfs_fc.md) from `gfs_fc`), renames variables to the standard convention, and assembles a per-station DataFrame indexed by time.

**Parameters:** `date`, `location_dict`, `variables_dict`, `gfs_data_dir`, `num_threads`.

**Returns:** a dict of `{station_ID: pd.DataFrame}`.

??? info "Internal helpers"
    - `append_timestamp()` — appends a single forecast hour's data into the running per-station DataFrame.
