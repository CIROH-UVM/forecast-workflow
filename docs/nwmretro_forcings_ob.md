---
icon: lucide/cloud-rain
---

# NWM v2.1 Retrospective Forcings

Downloads the meteorological forcing inputs (LDASIN files) used to drive version 2.1 of the National Water Model's retrospective run, from the [CIROH-hosted Zarr mirror on AWS](https://ciroh-nwm-zarr-retrospective-data-copy.s3.amazonaws.com/index.html). Unlike the v3.0 store, the v2.1 archive is exposed as per-timestep JSON references rather than a native Zarr dataset, so this module uses `kerchunk`'s `MultiZarrToZarr` to build a virtual, concatenated Zarr dataset on the fly. The dataset's native Lambert Conformal Conic grid is reprojected to lat/lon in order to locate each requested point.

## `get_data()`

Downloads and processes NWM v2.1 retrospective forcings data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` — the start date for which to grab data.
- `end_date` — the end date for which to grab data.
- `locations` — a dictionary mapping location name to a `(lat, lon)` coordinate tuple.
- `variables` — a dictionary mapping user-defined variable names to dataset variable names, drawn from the LDASIN forcing fields.

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series.

!!! note "Related module"
    For the newer, natively-Zarr NWM v3.0 retrospective forcings, see [`nwmretro_forcings3_ob`](nwmretro_forcings3_ob.md) — prefer that module unless v2.1-specific data is required.
