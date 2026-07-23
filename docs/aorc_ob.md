---
icon: lucide/cloud-rain
---

# NOAA AORC Meteorological Reanalysis

Downloads gridded, hourly meteorological reanalysis data from NOAA's Analysis of Record for Calibration (AORC) v1.1 dataset, reading directly from the [public Zarr store on AWS S3](https://registry.opendata.aws/noaa-nws-aorc/) (anonymous access, no download step required). For each requested location, the nearest available grid cell is selected.

## `get_data()`

Downloads and processes AORC reanalysis data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` — the start date for which to grab data.
- `end_date` — the end date for which to grab data.
- `locations` — a dictionary mapping location name to a `(lat, lon)` coordinate tuple.
- `variables` — a dictionary mapping user-defined variable names to AORC dataset variable names, e.g. `APCP_surface` (total precipitation), `TMP_2maboveground` (air temperature), `SPFH_2maboveground` (specific humidity), `DLWRF_surface`/`DSWRF_surface` (long/short-wave radiation), `PRES_surface` (pressure), `UGRD_10maboveground`/`VGRD_10maboveground` (wind components).

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series.
