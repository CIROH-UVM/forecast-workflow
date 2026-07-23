---
icon: lucide/cloud-rain
---

# National Water Model Forecast Forcings

Downloads and processes the meteorological forcing inputs (precipitation rate, air temperature, wind, humidity, pressure, radiation — see `NWM_FORCING_VARS`) used to drive NWM `medium_range` and `short_range` forecasts, from the [NWM Google Cloud Storage bucket](https://console.cloud.google.com/storage/browser/national-water-model). Unlike most acquisition modules, this one supports two output shapes: spatial subsetting by bounding box (returns a gridded `xr.Dataset`) or by a set of named points (returns the usual nested dict of Pandas Series).

## `get_data()`

Downloads and processes NWM forecast forcings data.

**Parameters:**

- `start_date` / `end_date` — the date range to retrieve data for.
- `member` — the NWM forecast member to get forcings for (`'medium_range'`, `'short_range'`, `'analysis_assim'`, `'analysis_assim_extend'`).
- `locations` — `None` (no spatial subsetting), `{'bbox': {'min_lat':..., 'max_lat':..., 'min_lon':..., 'max_lon':...}}`, or `{'points': {name: (lat, lon), ...}}`. See `validate_locations()`.
- `variables` — a dict/list of variables to extract, or `'all'` (default) to keep every forcing variable.
- `reference_date` — the forecast reference time. Defaults to `start_date` if `None`.
- `data_dir` — directory to store downloaded data. Defaults to the OS temp directory.
- `end_date_exclusive` — whether to exclude `end_date` from the returned series. Defaults to `True`.
- `dwnld_threads` — number of threads for downloading. Defaults to half the available CPU cores.

**Returns:**

- If `locations` is `None` or a bounding box: an `xr.Dataset` covering the requested area/variables. If `locations` specifies points: a nested dict where 1st-level keys are point names and 2nd-level keys are variable names, with values as Pandas Series.

## `download_nwm_forcings()`

Downloads the forcing NetCDF files for a given forecast member, reference date, and set of forecast hours from the GCS bucket.

**Parameters:** `reference_date`, `member`, `hours`, `download_dir`, `num_threads`.

**Returns:** a list of downloaded file paths.

## `process_nwm_forcings()`

Loads the downloaded forcing files, applies variable/location subsetting during load (via `preprocess_forcings_datasets()` as an `xarray.open_mfdataset` preprocessor for efficiency), localizes timestamps to UTC, and — for point subsetting — converts the result into a nested Series dictionary.

**Parameters:** `start_ts`, `end_ts`, `nwm_date_dir`, `member`, `reference_date`, `locations`, `variables`, `end_date_exclusive`.

**Returns:** an `xr.Dataset` or nested dict, matching `get_data()`'s return convention.

??? info "Internal helpers"
    - `prepForDownloads()` — builds the file name template and local directory structure for a given member/reference date.
    - `subset_by_bbox()` / `subset_by_points()` — perform the actual bounding-box clip (via `rioxarray`) or nearest-point selection (via a CRS transform into the forcing grid's projection).
    - `validate_locations()` — validates the `locations` argument and determines which subsetting method (`'bbox'`, `'points'`, or `None`) applies.
    - `parse_variables()` — normalizes the `variables` argument (list, dict, or `'all'`) into a `{user_name: dataset_name}` dict.
    - `preprocess_forcings_datasets()` — per-file preprocessing hook applying variable and location subsetting before concatenation.
    - `estimate_chunk_sizes_auto()` — estimates dask chunk sizes for loading based on available system memory.
