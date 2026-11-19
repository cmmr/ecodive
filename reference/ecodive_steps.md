# How a result was made

ecodive's results record each step that made them in one attribute per
step, so that the record travels with the result through later steps:

## Details

- `dist`: how a distance matrix was computed, set by
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  its metric functions (such as
  [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md)),
  [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md).

- `ord`: how an ordination was computed, set by
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
  [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md)
  and every other ordination.

Each record is a small list of class `ecodive_step`, which
[`str()`](https://rdrr.io/r/utils/str.html) and
[`print()`](https://rdrr.io/r/base/print.html) show on one line. Its
first element, `name`, is the step's display name (such as
`"Bray-Curtis Dissimilarity"` or `"Principal Coordinates Analysis"`), so
that any package can label it without ecodive's tables; the rest are the
settings it was made with. Read them as a list:
`attr(x, "dist")$metric`.

The records accumulate as results pass along: an ordination of a
distance matrix keeps the matrix's `dist` record beside its own `ord`
record, and the permutation tests
([`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
...) keep their distance matrix's `dist`, and
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) its
ordination's `dist` and `ord`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) titles its
plots with them, such as "Principal Coordinates Analysis of Bray-Curtis
Dissimilarity".

## The `dist` record

- `name`: The metric's display name, as
  [`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md)
  gives it.

- `metric`: Its ID, such as `"bray"`.

- `norm`: The normalization applied before computing it (by name, such
  as `"percent"`), whether chosen with `norm` or fixed by the metric
  (such as
  [`hellinger()`](https://cmmr.github.io/ecodive/reference/hellinger.md)'s).
  Absent for the UniFrac metrics.

- `pseudocount`: The pseudocount used (`0` for none), for the metrics
  that take one. `denominator`, for
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md).

- `power`, `alpha`: For
  [`minkowski()`](https://cmmr.github.io/ecodive/reference/minkowski.md)
  and
  [`generalized_unifrac()`](https://cmmr.github.io/ecodive/reference/generalized_unifrac.md).

- `clust`, `order`: With the `clust` argument: the clustering method,
  and each sample's position in the input.

## The `ord` record

- `name`: The method's display name, as
  [`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  and
  [`list_count_ords()`](https://cmmr.github.io/ecodive/reference/count_ords.md)
  give it.

- `method`: Its ID, such as `"pcoa"` (also the ordination's `method`
  attribute).

- The other arguments it was called with that are single values, such as
  `n_axes`, `correction`, `seed`, `norm` or `formula` (as text). Data
  arguments (`dm`, `counts`, `data`, ...) and `cpus`, which never
  changes the result, are left out.

## Examples

``` r
    ord <- pcoa(bray(iris[, 1:4], norm = 'percent'))
    str(ord)
#>  'ecodive_ord' num [1:150, 1:3] 0.177 0.164 0.174 0.153 0.176 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:3] "PCoA1" "PCoA2" "PCoA3"
#>  - attr(*, "method")= chr "pcoa"
#>  - attr(*, "eig")= num [1:3] 2.1337 0.0263 0.0153
#>  - attr(*, "prop")= num [1:3] 0.98139 0.0121 0.00704
#>  - attr(*, "trace")= num 2.17
#>  - attr(*, "correction")= chr "none"
#>  - attr(*, "constant")= num 0
#>  - attr(*, "dist")= Bray-Curtis Dissimilarity (metric = bray, norm = percent, pseudocount = 0)
#>  - attr(*, "ord")= Principal Coordinates Analysis (method = pcoa, n_axes = 3, correction = none)
    attr(ord, 'dist')
#> Bray-Curtis Dissimilarity (metric = bray, norm = percent, pseudocount = 0)
    attr(ord, 'dist')$name
#> [1] "Bray-Curtis Dissimilarity"
    attr(ord, 'ord')$n_axes
#> [1] 3
```
