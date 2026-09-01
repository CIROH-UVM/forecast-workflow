---
icon: lucide/ruler
---

# Extent

**Extent** is the fraction of a lake subdomain whose chlorophyll-a concentration exceeds a bloom threshold on a given day — the spatial footprint of a bloom.

It is the foundation of the other two metrics: [incidence](incidence.md) thresholds it, and [duration](duration.md) counts the result.

## `extent()` { #extent }

```r
extent(subdomain_layer, q_threshold)
```

Computes extent for a **single day**.

`subdomain_layer`

:   A 2D matrix or array of chlorophyll-a values in µg/L, typically one layer of a clipped raster stack. `NA` cells (outside the lake segment) are excluded from the calculation.

`q_threshold`

:   Chlorophyll-a threshold in µg/L. Cells strictly above this count as blooming.

**Returns** an unnamed numeric vector of length 2:

1. the fraction of non-`NA` cells that are blooming
2. the raw count of blooming cells

```r
e <- extent(mockarr[, , 1], q_threshold = 15)
e[1]   # 0.3421053 — fraction
e[2]   # 26 — cell count
```

The denominator is the number of non-`NA` cells, so it adapts to each subdomain's shape automatically. Because different layers of the same stack can have different `NA` patterns (cloud-masked satellite scenes, for example), the denominator is recomputed per layer rather than fixed across the stack.

### Threshold boundaries

Classification uses `terra::classify(..., right = TRUE)`, which makes intervals right-closed and left-open — `(lower, upper]`. The practical consequence:

| Chlorophyll-a | Classified as |
| --- | --- |
| `0` | Not bloom |
| below `q_threshold` | Not bloom |
| exactly `q_threshold` | **Not bloom** |
| above `q_threshold`, up to `200` | Bloom |

A cell sitting exactly on the threshold is *not* counted. Older analysis scripts in this repository instead offset the breakpoint (`q_threshold + 0.0001`) to achieve the same intent, so the two approaches can disagree on exact ties. With floating-point satellite and model output, exact ties are rare — but they are common in synthetic test data, which is worth remembering when writing tests.

!!! danger "Values above 200 µg/L are silently dropped from the bloom count"

    The reclassification matrix is hardcoded with an upper bound of `200`. Values **above** `200` fall outside every interval, and `terra::classify()` passes unmatched values through *unchanged* rather than setting them to `NA`. Such cells therefore land in the denominator but are never counted as blooming, which **understates** extent — badly, in exactly the situations that matter most.

    A subdomain where 40 % of cells are in a severe bloom at 350 µg/L reports:

    ```r
    lay <- matrix(1, 10, 10); lay[1:40] <- 350
    extent(lay, q_threshold = 15)[1]
    #> 0        <- should be 0.40

    incidence(lay, q_threshold = 15, areal_threshold = 0.20)
    #> 0        <- should be 1
    ```

    The same bloom at 150 µg/L correctly reports `0.40`. Negative values behave the same way.

    Whether this matters depends on your data: AEM3D `CYANO` output and satellite-derived chlorophyll for Lake Champlain generally sit well below 200 µg/L, so most existing results are unaffected. But the failure is silent, so it is worth checking `max(values(r), na.rm = TRUE)` on a new data source before trusting the output. Replacing the `200` bound with `Inf` removes the ceiling entirely and yields the correct `0.40`.

## `extent_wrap()` { #extent_wrap }

```r
extent_wrap(subdomain_arr, q_threshold)
```

Applies `extent()` across every layer of a **3D array**.

`subdomain_arr`

:   A 3D array with dimensions `(x, y, time)`.

**Returns** a `2 × N` matrix, where `N` is the number of layers:

- **row 1** — extent fractions, one per day
- **row 2** — blooming cell counts, one per day

```r
extents <- extent_wrap(mockarr, q_threshold = 15)

extents[1, ]   # fractions for every day
extents[2, ]   # cell counts for every day
```

!!! tip "Row 1 is what feeds the rest of the chain"

    `incidence_from_extent()` expects the fraction row specifically:

    ```r
    incidence_from_extent(extents[1, ], areal_threshold = 0.20)
    ```

    Passing the whole matrix, or row 2, will not error — it will quietly produce nonsense, since counts are compared against a fraction.

Converting from a raster is a one-liner:

```r
extent_wrap(as.array(my_raster), q_threshold = 15)
```
