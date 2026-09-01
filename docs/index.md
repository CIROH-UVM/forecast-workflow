---
icon: lucide/waves
---

# Welcome

This site documents **forecast-workflow**, the [CIROH](https://ciroh.ua.edu/) @ [UVM](https://www.uvm.edu/) pipeline for forecasting cyanobacteria harmful algal blooms (cyanoHABs) in Lake Champlain. The workflow acquires meteorological and hydrological data, computes tributary nutrient loads via concentration–discharge (CQ) relationships, drives the AEM3D 3D hydrodynamic and water quality model, and derives bloom metrics from the resulting chlorophyll-a output.

The source code lives at [CIROH-UVM/forecast-workflow](https://github.com/CIROH-UVM/forecast-workflow). This site documents the parts of the repository most useful on their own: the Python data acquisition modules and the R bloom metrics functions.

## Data acquisition

The modules under `data/` wrap upstream data providers behind a common interface. Unless noted otherwise, each module's `get_data()` returns a nested dictionary of Pandas Series, keyed first by location and then by variable:

```python
{location_name: {variable_name: pd.Series}}
```

The pages are grouped by the kind of data the module provides.

### Observations

Measured, monitored, and reanalysis data — records of what actually happened, used for model spin-up, hindcasting, and validation.

| Page | Data |
| --- | --- |
| [USGS Streamflow Observations](usgs_ob.md) | Streamflow and gage height from USGS NWIS streamgages |
| [Canadian Hydrometric Streamflow (CEHQ)](caflow_ob.md) | Streamflow for Quebec tributary gauges from the CEHQ archive |
| [NOAA Local Climatological Data (LCD)](lcd_ob.md) | Hourly station weather observations, e.g. Burlington International Airport |
| [UVM FEMC Meteorological Observations](femc_ob.md) | Weather observations from the Colchester Reef station on Lake Champlain |
| [NOAA AORC Meteorological Reanalysis](aorc_ob.md) | Gridded hourly reanalysis meteorology from the AORC v1.1 Zarr store |
| [NWM v2.1 Retrospective Forcings](nwmretro_forcings_ob.md) | Meteorological forcings behind the National Water Model v2.1 retrospective run |
| [NWM v3.0 Retrospective Forcings](nwmretro_forcings3_ob.md) | Meteorological forcings behind the National Water Model v3.0 retrospective run |
| [Sentinel-3 Cyanobacteria Index Imagery](sentinel3_ob.md) | Satellite cyanobacteria index rasters from NASA's Ocean Color CyAN portal |

### Forecasts

Operational forecast products and long-term model runs used to drive the lake model forward in time.

| Page | Data |
| --- | --- |
| [National Water Model Streamflow](nwm_fc.md) | NWM streamflow forecasts from Google Cloud Storage or NOMADS |
| [National Water Model Forecast Forcings](nwm_forcings_fc.md) | Meteorology driving NWM medium- and short-range forecasts |
| [NWM v3.0 Retrospective Streamflow](nwmv3_retro_fc.md) | Channel-routing output from the NWM v3.0 retrospective run |
| [NOAA GFS Forecast (THREDDS)](gfs_fc_thredds.md) | GFS 0.25° forecasts via NetCDF Subset Service — the recommended GFS module |
| [NOAA GFS Forecast (NOMADS)](gfs_fc.md) | GFS forecasts via the NOMADS GRIB filter — deprecated in favor of the THREDDS module |
| [NOAA CFS 9-Month Forecast](cfs_fc.md) | Long-range CFSv2 forecasts from the NCEI THREDDS server |

## Bloom metrics

The [Metrics](metrics/index.md) section documents `atomic_metrics.R`, the R functions that turn chlorophyll-a raster stacks from AEM3D into three bloom metrics for a lake segment:

- [Extent](metrics/extent.md) — what fraction of the segment is blooming each day
- [Incidence](metrics/incidence.md) — whether each day counts as a bloom event
- [Duration](metrics/duration.md) — how many bloom days fall within a rolling window
- [Subdomains](metrics/create-subdomain.md) — clipping a raster stack to a lake segment before computing anything

Start with the [Metrics overview](metrics/index.md), which explains how the three metrics layer on a single computation and how to compute all of them without repeating work.

## Contributing to these docs

The site is built with [Zensical](https://zensical.org/) from the Markdown files under `docs/` in the repository. The [Markdown in 5min](markdown.md) page is a quick syntax refresher for authoring new pages.
