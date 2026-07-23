---
icon: lucide/satellite
---

# Sentinel-3 Cyanobacteria Index (CI) Imagery

Downloads and processes Sentinel-3 OLCI Cyanobacterial Index (CI) satellite tiles from NASA's [Ocean Color CyAN](https://oceancolor.gsfc.nasa.gov/about/projects/cyan/) portal for a given area/tile ID, producing a time-stacked xarray Dataset of Digital Number (DN) rasters. Unlike the other acquisition modules, this one returns gridded imagery rather than a per-location Series dictionary, and its pipeline (download → reproject to WGS84 → crop to an area of interest → optional DN→CI conversion) is implemented as a stateful class, `CIFilesDownloadProcess`, which `get_data()` orchestrates.

## `get_data()`

Runs the full CI tile pipeline for a date range and area, returning the resulting raster time series.

**Parameters:**

- `start_date` / `end_date` — the date range to download CI tiles for.
- `appkey` — NASA Ocean Color app key for authentication (see the [CyAN docs](https://oceancolor.gsfc.nasa.gov/about/projects/cyan/) for how to obtain one).
- `cropbox` — `(lat, lon, lat, lon)` corners defining the area to crop each GeoTIFF to.
- `output_dir` — directory where output files (and intermediate temp files) are stored.
- `areaid` — string designating the tile/area ID to download, e.g. `"8_2"` for the Champlain Valley. Defaults to `"8_2"`.
- `convert_to_ci` — whether to convert Digital Number values to Cyanobacterial Index values in addition to returning the raw DN rasters. Defaults to `False` — as of 2025-07-30 this conversion has produced problematic results and is not currently recommended.
- `remove_temp` — whether to delete intermediate temp files after processing. Defaults to `True`.

**Returns:**

- An xarray Dataset with a `DN` variable stacked along a `time` dimension (one slice per downloaded tile), in the `EPSG:4326` (WGS84) CRS.

!!! note "Unused parameters"
    `locations`, `variables`, and `service` are accepted but not currently used by `get_data()` — the area downloaded is controlled entirely by `areaid` and `cropbox`.

??? info "Internal pipeline (`CIFilesDownloadProcess`)"
    `get_data()` instantiates this class and calls its methods in sequence; they are not intended to be used standalone:

    - `download_urls()` — queries the CyAN file-search API for tile URLs matching the date range and area, saving them to a text file.
    - `download_tiles()` — invokes `data/ci_data_download/ci-download.py` as a subprocess to download the listed tiles.
    - `convert_projection()` / `reproject_to_wgs84()` — reprojects downloaded GeoTIFFs to `EPSG:4326`.
    - `crop_aoi()` / `crop_using_gdal()` — crops reprojected tiles to `cropbox`.
    - `convert_dn_to_ci()` / `dn_to_ci()` — optional DN→CI conversion.
    - `remove_temp_dir()` — cleans up intermediate files.
    - `tif_to_ds()` — loads the cropped GeoTIFFs into a single time-stacked xarray Dataset.

## `plot()`

Plots a single time slice (2D array) of Sentinel-3 DN data, with an optional basemap.

**Parameters:**

- `da` — a 2D DN DataArray (a single time slice) to plot.
- `cmap` — colormap to use. If omitted, a custom colormap is built where DN value 254 (ground) renders white and 255 (no-data/cloud) renders black.
- `base` — whether to overlay a basemap using the DataArray's CRS. Defaults to `True`.
- `title` — plot title. If omitted, no title is set.

**Returns:**

- The matplotlib `Axes` object for the plot.
