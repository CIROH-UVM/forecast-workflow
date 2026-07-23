---
icon: lucide/cloud
---

# NOAA GFS Forecast (NOMADS)

Downloads and processes NOAA's GFS 0.25° forecast for a fixed Lake Champlain bounding box, using the [NOMADS GRIB filter service](https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl).

!!! warning "Deprecated"
    `get_data()` raises a `DeprecationWarning`. [`gfs_fc_thredds`](gfs_fc_thredds.md) is the recommended module for acquiring GFS forecast data going forward.

## `get_data()`

Downloads and processes GFS forecast data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `forecast_datetime` — the forecast launch date and cycle (00/06/12/18 UTC).
- `end_datetime` — the end date/time of the forecast window (GFS forecasts out to 16 days).
- `locations` — a dictionary mapping station name to a `(lat, lon)` tuple.
- `data_dir` — directory to store downloaded data. Defaults to the OS temp directory.
- `dnwld_threads` / `load_threads` — thread counts for downloading and reading grib files, respectively.
- `return_type` — `'dict'` (default) for a nested Series dict, or `'dataframe'` (not yet implemented).

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series.

## `download_gfs_threaded()`

Downloads the GRIB2 files for a forecast date across the requested hours, throttled to stay within NOMADS's requested rate limit (50 requests/minute, in chunks with a 60-second pause between them).

**Parameters:** `date`, `hours`, `gfs_data_dir`, `num_threads`.

## `process_gfs_data()`

Loads the downloaded GRIB2 files, isolates the requested locations, renames variables to the workflow's standard convention (`T2`, `TCDC`, `SWDOWN`, `U10`, `V10`, `RH2`, `RAIN`, `CPOFP`), and assembles a per-station DataFrame indexed by time.

**Parameters:** `location_dict`, `gfs_data_dir`, `num_threads`.

**Returns:** a dict of `{station_ID: pd.DataFrame}`.

??? info "Internal helpers"
    - `download_gfs()` — an older, non-throttled, non-threaded download function; superseded by `download_gfs_threaded()` and not currently called by `get_data()`.
    - `calibrate_columns()` / `append_timestamp()` — rename/reorder variable columns to the standard convention and append each timestep's data into the running per-station DataFrame.
    - `dict_to_csv()` — writes per-station DataFrames out to CSV (not currently called elsewhere).
    - `execute()` — generic subprocess runner that yields stdout line by line.
    - `isolate_loc_rows()` — selects the grid cells nearest each requested `(lat, lon)` out of a GFS dataset.
    - `remap_longs()` — remaps GFS's native 0–360 longitude convention to -180–180.
