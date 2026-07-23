---
icon: lucide/cloud-sun
---

# UVM FEMC Meteorological Observations (Colchester Reef)

Downloads observed meteorological data from the UVM [Forest Ecosystem Monitoring Cooperative](https://www.uvm.edu/femc/) (FEMC) Colchester Reef station on Lake Champlain. The underlying data is split across three sources with different formats and update cadences — a large historical "v1" archive, a "v2" archive starting 2024-07-26, and a rolling "latest" feed covering roughly the last 21 days — which this module stitches together transparently.

## `get_data()`

Downloads and processes FEMC observational data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` — the start date for which to grab data.
- `end_date` — the end date for which to grab data.
- `locations` — a dictionary of location name to station ID. Defaults to `{'CR':'ColReefQAQC'}` (Colchester Reef is currently the only supported station).
- `variables` — a dictionary of variables to get; keys are user-defined names, values must be one of the FEMC abbreviations `T2` (air temp), `SWDOWN` (shortwave radiation), `RH2` (relative humidity), `WSPEED` (wind speed), `WDIR` (wind direction).

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series with units attached.

## `load_femc_data()`

Loads and merges the raw v1, v2, and latest datasets for the requested date range, choosing only the sources actually needed to cover it.

**Parameters:**

- `start_date` — the start date for which to load data.
- `end_date` — the end date for which to load data (cannot be in the future).

**Returns:**

- A Pandas DataFrame of raw meteorological data spanning the requested range, indexed by timestamp and sorted chronologically. The v1 dataset is loaded lazily via Dask due to its size; v2 and the latest feed are loaded with pandas.
