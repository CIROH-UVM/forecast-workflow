---
icon: lucide/chart-line
---

# Incidence

**Incidence** answers a yes-or-no question: was there a bloom on this day? A day counts as a bloom day when its [extent](extent.md) reaches `areal_threshold` — that is, when a large enough share of the subdomain is above the chlorophyll-a threshold.

The output is a `1`/`0` vector, one value per day, which [duration](duration.md) then counts within rolling windows.

## `incidence()` { #incidence }

```r
incidence(subdomain_layer, q_threshold, areal_threshold)
```

Computes incidence for a **single day**.

`subdomain_layer`

:   A 2D matrix or array of chlorophyll-a values in µg/L.

`q_threshold`

:   Chlorophyll-a threshold in µg/L, passed through to [`extent()`](extent.md#extent).

`areal_threshold`

:   Fraction of the subdomain that must be blooming for the day to count. `0.20` is the project default.

**Returns** `1` if the day is a bloom day, `0` otherwise.

```r
incidence(mockarr[, , 1], q_threshold = 15, areal_threshold = 0.20)
#> 1
```

### The two thresholds point in opposite directions

This trips people up, so it is worth stating plainly:

| Threshold | Comparison | A value exactly on the boundary is… |
| --- | --- | --- |
| `q_threshold` (per cell) | strictly greater than | **not** a bloom |
| `areal_threshold` (per day) | greater than or equal | **is** a bloom |

So a subdomain with extent of exactly `0.20` at `areal_threshold = 0.20` returns `1`:

```r
lay <- matrix(c(rep(30, 2), rep(1, 8)), nrow = 2, ncol = 5)
extent(lay, 15)[1]                 #> 0.2
incidence(lay, 15, 0.20)           #> 1
```

The asymmetry is inherited from the underlying implementations — `terra::classify(right = TRUE)` for cells, a plain `>=` for days — rather than chosen deliberately, but it is stable and worth relying on.

## `incidence_wrap()` { #incidence_wrap }

```r
incidence_wrap(subdomain_arr, q_threshold, areal_threshold)
```

Applies `incidence()` across every layer of a **3D array**.

**Returns** an integer vector of `1`s and `0`s, one per layer.

```r
incidences <- incidence_wrap(mockarr, q_threshold = 15, areal_threshold = 0.20)
#> [1] 1 1 0 0 1 1 0 1 1 1 ...
```

Works on rasters once cast to an array:

```r
incidence_wrap(as.array(mockrast), q_threshold = 15, areal_threshold = 0.20)
```

!!! warning "This recomputes extent"

    `incidence_wrap()` calls `extent()` internally for every layer. If you have already computed extent for this stack, use `incidence_from_extent()` below instead and skip the duplicated work.

## `incidence_from_extent()` { #incidence_from_extent }

```r
incidence_from_extent(extent_fractions, areal_threshold)
```

Derives incidence from extent results you already have. **This is the efficient path** and the one used in production analysis scripts.

`extent_fractions`

:   A numeric vector of extent fractions — row 1 of an [`extent_wrap()`](extent.md#extent_wrap) result.

`areal_threshold`

:   Fraction of the subdomain that must be blooming for the day to count.

**Returns** an integer vector of `1`s and `0`s.

```r
extents <- extent_wrap(mockarr, q_threshold = 15)
incidences <- incidence_from_extent(extents[1, ], areal_threshold = 0.20)
```

The whole function is a single vectorized comparison, `as.integer(extent_fractions >= areal_threshold)`, so it is effectively free next to the raster work `extent_wrap()` already did.

Note that it takes **no** `q_threshold` — that choice was already baked into the extent fractions you are passing in. If you want to vary `q_threshold`, you must recompute extent.

!!! tip "All three paths agree"

    The example script verifies this explicitly, which is a useful regression check after editing any of these functions:

    ```r
    all(incidences_ext == incidences_wrp)        #> TRUE
    all(incidences_wrp == incidences_wrp_rast)   #> TRUE
    ```
