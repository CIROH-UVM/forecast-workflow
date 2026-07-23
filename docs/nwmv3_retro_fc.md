---
icon: lucide/droplet
---

# NWM v3.0 Retrospective Streamflow

Downloads channel-routing output (CHRTOUT variables — streamflow, velocity, lateral/bucket runoff) from version 3.0 of the National Water Model's single long-term [retrospective run](https://registry.opendata.aws/nwm-archive/), reading directly from the public Zarr store on AWS. Unlike the operational forecast products, NWM retrospective runs are produced only once per model version, so there is no forecast cycle or reference date to specify.

## `get_data()`

Downloads and processes NWM v3.0 retrospective CHRTOUT data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` / `end_date` — the date range to retrieve data for.
- `locations` — a dictionary mapping reach name to NWM `feature_id` (reach ID).
- `variables` — a dictionary mapping user-defined variable names to CHRTOUT variable names. Defaults to `{'streamflow':'streamflow'}`. Other available variables include `qBtmVertRunoff`, `qBucket`, `qSfcLatRunoff`, `q_lateral`, and `velocity`.

**Returns:**

- A nested dict where 1st-level keys are user-provided reach names and 2nd-level keys are variable names, with values as Pandas Series named `"{var} ({unit})"`. Series are end-date exclusive.
