---
icon: lucide/cloud-sun
---

# NOAA Local Climatological Data (LCD)

Downloads and cleans hourly station weather observations (temperature, precipitation, wind, relative humidity, sky cover) from NOAA's [Local Climatological Data](https://www.ncei.noaa.gov/access/services/support/v3/datasets.json) (LCD) API, e.g. from the Burlington International Airport station. The raw API response requires substantial cleanup — duplicate report handling, suspect/missing value markers, trace-precipitation encoding — which this module handles internally.

## `get_data()`

Downloads and processes LCD data, returning a nested dictionary of Pandas Series for each variable, for each location.

**Parameters:**

- `start_date` — the start date for which to grab data.
- `end_date` — the end date for which to grab data.
- `locations` — a dictionary of location name to LCD station ID.
- `variables` — a dictionary of variables to download; keys are user-defined names, values are LCD field names. Only the fields in `standard_var_units`/`metric_var_units` (`HourlySkyConditions`, `HourlyPrecipitation`, `HourlyDryBulbTemperature`, `HourlyRelativeHumidity`, `HourlyWindSpeed`, `HourlyWindDirection`) have been tested.
- `units` — unit system for the request: `'standard'` (US units, default) or `'metric'`.

**Returns:**

- A nested dict where 1st-level keys are user-provided location names and 2nd-level keys are variable names, with values as Pandas Series named `"{user_name} ({units})"`.

## `lcdRequest()`

Sends the raw request to the NOAA LCD API for a single station, variable list, date range, and unit system, retrying until a valid JSON response is received.

**Parameters:**

- `startDate` / `endDate` — request date range (time-of-day is ignored; the request always covers full UTC days).
- `var_list` — LCD-specific variable names to request.
- `station_id` — LCD station ID.
- `units` — `'standard'` or `'metric'`.

**Returns:**

- A Pandas DataFrame of the raw API response, one row per report.

??? info "Internal cleaning helpers"
    These handle the LCD API's quirks and are not intended to be called directly:

    - `clean_raw_df()` — deduplicates report timestamps (preferring `FM-16` > `FM-15` > `FM-12` > `SOD` reports) and indexes by UTC time.
    - `scrubSpecialChars()` — strips `'s'` (suspect) and `'*'` (missing) indicator characters and casts to float.
    - `process_clouds()` / `splitsky()` / `sky2prop()` — parse raw sky-condition strings (e.g. `'SCT:04 015'`) into a fractional sky-cover value.
    - `process_rain()` / `leavenotrace()` — convert trace ("T") precipitation readings to 0.00 and flag malformed values as NaN. `process_rain()` is not currently called by `get_data()`, which handles precipitation cleanup inline via `leavenotrace()` directly.
