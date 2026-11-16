# Ordination Wrapper Function

Runs any of ecodive's ordination methods by name, as
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
does for beta diversity metrics.

## Usage

``` r
ordination(x, method, ...)
```

## Arguments

- x:

  The data to ordinate: for `"pcoa"` and `"nmds"`, a `dist` object such
  as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  or a symmetric matrix of distances. See
  [`list_ordinations()`](https://cmmr.github.io/ecodive/reference/ordinations.md)
  for what each method takes.

- method:

  The name of an ordination method: `"nmds"` or `"pcoa"`. Flexible
  matching is supported (see below). Programmatic access via
  [`list_ordinations()`](https://cmmr.github.io/ecodive/reference/ordinations.md).

- ...:

  Further arguments for the method's function, such as `k` or `cpus`.
  Each method keeps its own defaults; see
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) and
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md).

## Value

The method's result, as described in "Ordination result" below.
`ordination(x, 'pcoa', ...)` is the same as `pcoa(x, ...)`.

## Details

**List of Ordination Methods**

|                        |                                            |
|------------------------|--------------------------------------------|
| Option / Function Name | Method Name                                |
| `nmds`                 | Non-metric Multidimensional Scaling (NMDS) |
| `pcoa`                 | Principal Coordinates Analysis (PCoA)      |

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters.
E.g. `method = 'PCoA'` or `'principal coord'` selects `pcoa`. Common
other names also work: `'MDS'`, `'PCO'`, `'cmdscale'`, or
`'classical MDS'` for `pcoa`, and `'non-metric MDS'`, `'metaMDS'`, or
`'monoMDS'` for `nmds`.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the list of methods.

## Ordination result

Every ordination method
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md), and
`ordination()`) returns the same kind of object: a plain numeric matrix
with one row per sample (named after the samples) and one column per
axis (named `"PCoA1"`, `"NMDS1"`, etc.), so that it can go straight into
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) or a data
frame. Its attributes always begin with these three, in this order:

- `method`: The method's ID, e.g. `"pcoa"` or `"nmds"`.

- `eig`: Each axis's eigenvalue, or `NA` for methods without them.

- `prop`: The proportion of the total variance that each axis explains,
  or `NA` where that isn't defined.

Attributes particular to the method follow them.

## See also

[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`list_ordinations()`](https://cmmr.github.io/ecodive/reference/ordinations.md)

## Examples

``` r
    d <- bray(ex_counts)

    ordination(d, 'pcoa')
#>             PCoA1       PCoA2        PCoA3
#> Saliva -0.3891250 -0.03721713  0.296748867
#> Gums   -0.3951099 -0.05550727 -0.293514268
#> Nose    0.3179676  0.53707433 -0.006811864
#> Stool   0.4662673 -0.44434993  0.003577264
#> attr(,"method")
#> [1] "pcoa"
#> attr(,"eig")
#> [1] 0.6260388 0.4903619 0.1742697
#> attr(,"prop")
#> [1] 0.4850493 0.3799280 0.1350226
#> attr(,"trace")
#> [1] 1.29067

    ordination(d, 'NMDS', k = 2)
#> Warning: stress is (nearly) zero: there may be too few samples
#>             NMDS1       NMDS2
#> Saliva -0.3891463 -0.03721917
#> Gums   -0.3951316 -0.05551031
#> Nose    0.3179851  0.53710375
#> Stool   0.4662929 -0.44437427
#> attr(,"method")
#> [1] "nmds"
#> attr(,"eig")
#> [1] NA NA
#> attr(,"prop")
#> [1] NA NA
#> attr(,"stress")
#> [1] 0
#> attr(,"converged")
#> [1] FALSE
#> attr(,"tries")
#> [1] 20
#> attr(,"best_try")
#> [1] 0
```
