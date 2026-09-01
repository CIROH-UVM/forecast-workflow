---
icon: lucide/timer
---

# Duration

**Duration** measures bloom persistence: within a rolling window of days, how many were bloom days? It counts [incidence](incidence.md) values inside each window and converts that count into elapsed time.

Note that duration counts bloom days *within* a window — it does not measure the length of an unbroken bloom episode. A window containing five scattered bloom days and a window containing five consecutive bloom days both report 5 days.

## `duration_from_incidence()` { #duration_from_incidence }

```r
duration_from_incidence(incidence_vect, time_index, window_size, step_size)
```

The primary function. Derives duration from an incidence vector you already have.

`incidence_vect`

:   Integer vector of `1`s and `0`s from any of the [incidence](incidence.md) functions.

`time_index`

:   A `POSIXct` vector of timestamps, one per element of `incidence_vect`. Must be the **same length** — enforced by `stopifnot()`.

`window_size`

:   Window length, in number of layers (days).

`step_size`

:   How far the window advances between evaluations, in layers. Set it equal to `window_size` for non-overlapping windows, or smaller for overlapping ones.

**Returns** an **unnamed list of two vectors**, accessed positionally:

```r
dur <- duration_from_incidence(incidences, time(mockrast),
                               window_size = 7, step_size = 7)

window_ends <- dur[[1]]   # POSIXct, one per window
bloom_time  <- dur[[2]]   # lubridate Duration, one per window
```

!!! warning "`dur[[1]]` holds window **end** timestamps, not centers"

    Some code and comments in this repository label this vector as window centers — `window_centers` in `atomic_metrics_example.R`, and a `Window_Center` column in `rb_metrics_example.R`. Those labels are wrong. The function stores `time_index[right_idx]`, the **last timestamp included in each window**:

    ```r
    inc <- c(1,1,0,0,1, 1,0,1,1,1)
    tt  <- seq(as.POSIXct("2019-05-01", tz = "UTC"), by = "day", length.out = 10)
    d   <- duration_from_incidence(inc, tt, window_size = 5, step_size = 5)

    format(d[[1]], "%m-%d")
    #> "05-05" "05-10"
    ```

    The first window spans 05-01 through 05-05 and reports `05-05`. Plot or join on these values as right-edge labels; treating them as centers shifts every window by half its width.

### Return units

The second element is a [`lubridate`](https://lubridate.tidyverse.org/) `Duration` object, which prints and stores in **seconds**. Convert to days explicitly:

```r
as.numeric(dur[[2]], "days")
#> [1] 3 4
```

The conversion factor comes from the spacing of the first two timestamps, `time_index[2] - time_index[1]`, which assumes an **evenly spaced** time index. With daily rasters this is exactly right. With irregular data — satellite scenes with gaps, for example — the first interval silently sets the scale for every window, so gap-fill to a regular grid first.

### Windows that do not divide evenly are truncated

Window start positions are `seq(1, length(incidence_vect) - window_size + 1, by = step_size)`. Any trailing days that cannot fill a complete window are dropped:

```r
# 10 days, 4-day windows, 4-day steps
duration_from_incidence(rep(1, 10), tt, window_size = 4, step_size = 4)[[1]]
#> covers days 1–4 and 5–8; days 9 and 10 are dropped
```

!!! danger "`window_size` must not exceed the number of layers"

    If `window_size > length(incidence_vect)`, the `seq()` call receives a negative endpoint with a positive step and R raises a cryptic error:

    ```
    Error in seq.default(...) : wrong sign in 'by' argument
    ```

    This is easy to hit when a bloom season is shorter than the window — a 30-day window over a partial-season stack, for instance. Check `nlyr(r) >= window_size` before calling.

## `duration()` { #duration }

```r
duration(subdomain_rast, q_threshold, areal_threshold, window_size, step_size)
```

Convenience wrapper that goes from raster to duration in one call, computing extent and incidence internally.

`subdomain_rast`

:   A `SpatRaster` **with a time index set**. This is the one metric function that takes a raster rather than an array, because it needs `terra::time()` to build the window timestamps.

The remaining arguments match the functions it delegates to.

**Returns** the same two-element list as `duration_from_incidence()`.

```r
dur <- duration(subdomain_rast = mockrast,
                q_threshold = 15, areal_threshold = 0.20,
                window_size = 7, step_size = 7)
```

Setting a time index on a raster that lacks one:

```r
time(mockrast) <- seq(as.POSIXct("2019-05-28", tz = "UTC"),
                      by = "day", length.out = nlyr(mockrast))
```

!!! warning "Recomputes everything"

    `duration()` calls `incidence_wrap()`, which calls `extent()` for every layer. In a pipeline that already reports extent and incidence, this triples the raster work. Prefer:

    ```r
    extents    <- extent_wrap(as.array(r), q_threshold = 15)
    incidences <- incidence_from_extent(extents[1, ], areal_threshold = 0.20)
    durations  <- duration_from_incidence(incidences, as.POSIXct(time(r)),
                                          window_size = 30, step_size = 15)
    ```

    Reach for `duration()` when duration is the *only* metric you need.
