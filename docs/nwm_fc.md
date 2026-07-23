---
icon: lucide/droplet
---

# National Water Model Acquisition Module

Downloads and processes National Water Model (NWM) streamflow forecast data from either [Google Cloud Storage](https://console.cloud.google.com/storage/browser/national-water-model) (GCS, holds NWM data back to 2018) or [NOMADS](https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/) (holds only the last 2 days' runs). Supports any forecast member (`medium_range_mem1`, `long_range_mem3`, `short_range`, `analysis_assim`, etc.) at any model cycle, but only downloads `channel_rt` files, which contain `streamflow`, `velocity`, and `nudge` variables — `land` and `reservoir` files are not currently supported.

## `get_data()`

Downloads and processes NWM hydrology forecast data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` / `end_date` — the date range to retrieve data for.
- `member` — the NWM forecast member to get (e.g. `medium_range_mem1`, `long_range_mem3`, `short_range`, `analysis_assim`).
- `locations` — a dictionary mapping location name to reach/gauge ID (`feature_id`).
- `variables` — a dictionary mapping user-defined variable names to `channel_rt` variable names (`streamflow`, `velocity`, `nudge`).
- `reference_date` — the forecast reference time (launch date/cycle). Defaults to `start_date` if `None`. For `analysis_assim`, this is required and represents the corresponding short-range forecast launch time.
- `data_dir` — directory to store downloaded data. Defaults to the OS temp directory.
- `format` — `'dictionary'` (default) for a nested Series dict, or `'xarray'` for an `xr.Dataset`.
- `gcs` — whether to use the GCS bucket (`True`, default) instead of NOMADS.
- `end_date_exclusive` — whether to exclude `end_date` from the returned series. Defaults to `True`.
- `dwnld_threads` / `load_threads` — thread counts for downloading and reading data, respectively.

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series — or an `xr.Dataset` if `format='xarray'`.

## `download_nwm()`

Downloads the `channel_rt` NetCDF files for a single forecast product (member + reference date) from GCS or NOMADS.

**Parameters:** `reference_date`, `member`, `hours` (`'all'`, an int, or a list of ints), `gcs`, `download_dir`, `num_threads`.

**Returns:** a list of downloaded file paths.

## `process_nwm()`

Loads and slices NWM forecast NetCDF files by time, location, and variable, returning the result in the requested format.

**Parameters:** `start_ts`, `end_ts`, `locations`, `variables`, `nwm_data_dir`, `fname_template`, `format`, `end_date_exclusive`, `num_threads`.

**Returns:** a nested dict of Pandas Series, or an `xr.Dataset`, matching `get_data()`'s return convention.

??? info "Internal helpers"
    - `prepForDownloads()` — parses the reference date and builds the file name template and local directory structure for a given member.