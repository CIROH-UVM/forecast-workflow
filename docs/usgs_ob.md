---
icon: lucide/droplet
---

# USGS Streamflow Observations

Downloads observed streamflow, gage height, and related hydrology data from the [USGS NWIS web services](https://waterservices.usgs.gov/docs/), supporting both instantaneous and daily-values services.

## `get_data()`

Downloads and processes USGS observational hydrology data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` — the start date for which to grab data.
- `end_date` — the end date for which to grab data.
- `locations` — a dictionary mapping station name to USGS site ID.
- `variables` — a dictionary mapping user-defined variable names to [USGS parameter codes](https://help.waterdata.usgs.gov/parameter_cd?group_cd=PHY) (e.g. `'00060'` for mean daily streamflow, `'00061'` for instantaneous streamflow, `'00065'` for gage height). Defaults to `{'streamflow':'00060'}`.
- `service` — which USGS service to use: `'iv'` (instantaneous values, default) or `'dv'` (daily values). See the [USGS docs](https://waterservices.usgs.gov/docs/) for other options.

**Returns:**

- A nested dict where 1st-level keys are user-provided station names and 2nd-level keys are variable names, with values as Pandas Series named `"{var} ({unit})"`.

## `USGSgetvars_function()`

Builds the USGS NWIS request URL for a single station and set of variables, and returns the parsed response. Retries on failure.

**Parameters:**

- `id` — USGS site ID.
- `variables` — a dictionary mapping user-defined variable names to USGS parameter codes.
- `start` / `end` — start and end datetimes for the request.
- `service` — `'iv'` or `'dv'`.

**Returns:**

- A dictionary mapping variable name to a Pandas Series of that variable's data, indexed by timestamp.
