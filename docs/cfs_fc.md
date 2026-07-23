---
icon: lucide/cloud-sun
---

# NOAA CFS (Climate Forecast System) 9-Month Forecast

Downloads and processes NOAA's [CFSv2 operational 9-month forecast](https://www.ncei.noaa.gov/products/weather-climate-models/climate-forecast-system) (available from 2011-04-01 to present) from the [NCEI THREDDS server](https://www.ncei.noaa.gov/thredds/catalog/model-cfs_v2_for_ts/catalog.html). Two retrieval strategies are supported: downloading full forecast files (`fileServer`), or requesting a server-side spatial/temporal slice via the NetCDF Subset Service (`ncss`) — useful for reducing bandwidth, though not always suitable for shared network drives.

## `get_data()`

Downloads and processes CFS forecast data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` / `end_date` — the date range to retrieve data for.
- `locations` — a dictionary mapping location name to a `(lat, lon)` tuple.
- `variables` — a dictionary mapping user-defined variable names to CFS long names (keys of `CFS_VAR_NAMES`, e.g. `'Temperature_height_above_ground'`, `'Precipitation_rate_surface'`).
- `reference_date` — the forecast reference time (initialization time). Defaults to `start_date` if `None`. Must align to a CFS forecast cycle (0, 6, 12, or 18 UTC).
- `ncss` — whether to use the NetCDF Subset Service instead of downloading full files. Defaults to `False`.
- `end_date_exclusive` — whether to exclude `end_date` from the returned series. Defaults to `True`.
- `data_dir` — directory to store downloaded data. Defaults to the OS temp directory.
- `num_threads` — number of threads for downloading. Defaults to half the available CPU cores.

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series with units attached.

## `download_full_cfs()` / `download_ncss()`

Download the full forecast files, or a server-subsetted slice, respectively, for a given reference date and variable set.

**Parameters:** `reference_date`, `variables`, `data_dir`, `num_threads` (both); `download_ncss()` additionally takes `locations`, `start_date`, `end_date` to build the spatial/temporal subset request.

**Returns:** a list of local file paths for the downloaded files.

## `process_full_cfs()` / `process_ncss()`

Extract per-location, per-variable time series from a downloaded full-file dataset or NCSS dataset, respectively, mapping requested `(lat, lon)` locations to the nearest available grid coordinates.

**Parameters:** `ds` (the opened xarray Dataset), `variables`, `locations`.

**Returns:** a nested dict of Pandas Series (see `get_data()`), with units attached.

??? info "Internal helpers"
    - `construct_download_paths()` — builds THREDDS download URLs and local file paths for a given reference date and variable set.
    - `execute_downloads()` — runs the actual multi-threaded download given a list of URLs and destination paths.
    - `remap_to_cfs_coords()` — converts `(lat, lon)` tuples from -180/180 to the 0/360 longitude convention CFS uses.
