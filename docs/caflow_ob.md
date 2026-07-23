---
icon: lucide/droplet
---

# Canadian Hydrometric Streamflow (CEHQ)

Downloads observed streamflow data for Quebec tributary gauges (e.g. Pike and Rock, tributaries of Missisquoi Bay) from the [CEHQ](https://www.cehq.gouv.qc.ca/) (Centre d'expertise hydrique du Québec) historical data archive. Supports both instantaneous (15-minute) and daily-averaged data.

## `get_data()`

Downloads and processes Canadian observational hydrology data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` — the start date for which to grab data.
- `end_date` — the end date for which to grab data.
- `locations` — a dictionary of station name to CEHQ station ID. Defaults to `{'Pike':'030424', 'Rock':'030425'}`.
- `variables` — a dictionary of variables to download; the only currently supported value is `'Débit (m³/s)'` (streamflow). Defaults to `{'streamflow':'Débit (m³/s)'}`.
- `service` — which frequency of data to get: `'iv'` for instantaneous (15-minute, default) or `'dv'` for daily. See the [station page](https://www.cehq.gouv.qc.ca/hydrometrie/historique_donnees/fiche_station.asp?NoStation=030425) for details.

**Returns:**

- A nested dict where 1st-level keys are user-provided station names and 2nd-level keys are variable names, with values as Pandas Series in m³/s.

## `get_instantaneous()`

Downloads and concatenates instantaneous (15-minute frequency) streamflow records for a single station across one or more years.

**Parameters:**

- `id` — CEHQ station ID.
- `yearlist` — years to download instantaneous data for. Years with a failed request are skipped.

**Returns:**

- A Pandas DataFrame of concatenated instantaneous streamflow records, indexed by date/time.

## `get_daily()`

Downloads daily streamflow records for a single station.

**Parameters:**

- `id` — CEHQ station ID.

**Returns:**

- A Pandas DataFrame of daily streamflow records, indexed by date.
