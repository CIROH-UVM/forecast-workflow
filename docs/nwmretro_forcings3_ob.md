---
icon: lucide/cloud-rain
---

# NWM v3.0 Retrospective Forcings

Downloads the meteorological forcing inputs (precipitation rate, 2-m air temperature) used to drive version 3.0 of the National Water Model's retrospective run, reading directly from the [public v3.0 retrospective Zarr store on AWS](https://noaa-nwm-retrospective-3-0-pds.s3.amazonaws.com/index.html#CONUS/zarr/forcing/). Each requested lat/lon location is reprojected into the forcing grid's Lambert Conformal Conic coordinates and matched to the nearest grid cell.

## `get_data()`

Downloads and processes NWM v3.0 retrospective forcings data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` — the start date for which to grab data.
- `end_date` — the end date for which to grab data.
- `locations` — a dictionary mapping location name to a `(lat, lon)` coordinate tuple.
- `variables` — a dictionary mapping user-defined variable names to dataset variable names. Supported values are `'RAINRATE'` (precipitation rate) and `'T2D'` (2-m air temperature).

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series (end-date exclusive).

!!! note "Related module"
    For the older NWM v2.1 retrospective forcings, see [`nwmretro_forcings_ob`](nwmretro_forcings_ob.md).
