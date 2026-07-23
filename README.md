# Welcome to CIROH @ UVM forecast-workflow Repository

This repository was developed by researchers at CIROH @ UVM to run multi-year, multi-scenario cyanobacteria harmful algal bloom (cyanoHAB) forecasts for Lake Champlain. It orchestrates the acquisition and preprocessing of meteorological and hydrological data, computes nutrient loading via CQ relationships, and drives the AEM3D 3D hydrodynamic and water quality model to produce forecasts of lake conditions.

## Data acquisition modules

The most broadly useful part of this repository to other researchers is the collection of data-acquisition modules under [`data/`](data/). Each module wraps the details of a specific upstream data source behind a common interface, returning ready-to-use time series for locations of interest. Sources currently supported include:

- **Streamflow**: NOAA's National Water Model (NWM) forecasts and retrospective runs, USGS streamgages, and CEHQ (Quebec) streamgages
- **Meteorology**: NOAA GFS and CFS forecasts, NOAA Local Climatological Data (LCD) observations, NWM forecast forcings, and AORC gridded meteorology
- **Water quality**: FEMC tributary monitoring data and Sentinel-3 satellite imagery

If you're building a similar hydrology or water-quality forecasting pipeline and need working code to pull one of these data sources, these modules are a good place to start.

For full documentation of each module see the [documentation site](https://ciroh-uvm.github.io/forecast-workflow/).
