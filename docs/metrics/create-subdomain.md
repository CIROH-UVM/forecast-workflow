---
icon: lucide/scissors
---

# Subdomains

Bloom metrics are only meaningful for a defined region of the lake. Missisquoi Bay, St. Albans Bay, and the Inland Sea bloom on different schedules and at different intensities, so averaging across all of Lake Champlain would wash out the signal each metric is meant to capture.

`create_subdomain()` clips a raster stack to one lake segment before any metric is computed. Every analysis begins with it.

## `create_subdomain()` { #create_subdomain }

```r
create_subdomain(r, shapefile)
```

`r`

:   A `SpatRaster` — a single layer or a full stack. Stacks are clipped layer by layer, so a whole season can be handled in one call.

`shapefile`

:   A `SpatVector` defining the lake segment, loaded with `terra::vect()`.

**Returns** a new `SpatRaster`, masked and cropped to the segment. The input is not modified.

```r
StAlbans <- vect(file.path(working_dir, "Detailed_Lake_Segments/StAlbans_prj.shp"))

baseline  <- rast(file.path(working_dir, "CYANO_DailyAvg_MaySep_2019.nc"))
subdomain <- create_subdomain(baseline, StAlbans)
```

It prints a progress line naming the segment as it works, which is useful when looping over several subdomains:

```
Masking and cropping raster for lake segment: StAlbans
```

The name is read from the shapefile's `SEG_NAME` attribute. A `SpatVector` without that field still clips correctly — the message simply comes out blank.

### Mask, then crop

The two operations do different jobs and the order matters:

1. **`mask()`** sets cells *outside* the polygon to `NA`, but leaves the raster's extent unchanged.
2. **`crop()`** then shrinks the extent to the polygon's bounding box.

Masking first means the cells that survive cropping are already `NA` outside the segment boundary. Cropping first would trim to the bounding box while leaving the corners — the parts of the box outside the irregular polygon — holding real values, which would then be counted as lake area.

This connects directly to how [`extent()`](extent.md#extent) works: its denominator is the count of non-`NA` cells. Masking is what makes that denominator equal the true area of the lake segment rather than the area of a rectangle drawn around it.

### CRS must match

`mask()` warns rather than errors when the raster and vector use different coordinate reference systems:

```
Warning message:
[mask] CRS do not match
```

The result will be wrong — likely empty or nonsensical — so treat this warning as fatal. The `_prj` suffix on the project's shapefiles (`StAlbans_prj.shp`) marks them as already projected to match the model output. Reproject explicitly if you bring in a segment definition from elsewhere:

```r
shapefile <- project(shapefile, crs(r))
```

## Available segments

The project's segment shapefiles live under `Detailed_Lake_Segments/` in the preprocessing results directory:

| Segment | Shapefile |
| --- | --- |
| Missisquoi Bay | `Missisquoi_prj.shp` |
| St. Albans Bay | `StAlbans_prj.shp` |
| Inland Sea (Northeast Arm) | `NEarm_prj.shp` |

Note that the Inland Sea is referred to as the Northeast Arm in filenames but as `InlandSea` in analysis-script variables — the same region under two names.

```r
shapefiles_list <- list(
  Missisquoi = vect(file.path(seg_dir, "Missisquoi_prj.shp")),
  StAlbans   = vect(file.path(seg_dir, "StAlbans_prj.shp")),
  InlandSea  = vect(file.path(seg_dir, "NEarm_prj.shp"))
)
```

### Checking the result

Masking and cropping are easy to get wrong in ways that are invisible in summary statistics but obvious on a map. Plot before and after:

```r
plot(baseline$CYANO_1)     # full lake
plot(subdomain$CYANO_1)    # clipped to the segment
```
