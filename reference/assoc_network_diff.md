# Differential Associations

Tests, for every pair of features, whether their correlation differs
between two groups of samples, with Fisher's \\z\\ test, as NetCoMi's
`diffnet(diffMethod = "fisherTest")` does.

## Usage

``` r
assoc_network_diff(x, y, alpha = 0.05, p_adjust = "BH", n_samples = NULL)
```

## Arguments

- x, y:

  The two groups' correlations between the same features: two results of
  the same method, one whose values are Pearson-type correlations
  (today,
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)), or
  two symmetric matrices of Pearson correlations with 1 on the diagonal
  (with `n_samples`). Features are matched by name.

- alpha:

  The adjusted p-value at most which a difference is called significant
  (the `reject` column). Default `0.05`.

- p_adjust:

  How to adjust the p-values for testing every pair: one of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html).
  Default `'BH'`.

- n_samples:

  The numbers of samples that `x` and `y` were estimated from: one
  number for both, or two. Needed for matrices; `NULL` (the default)
  takes them from the results' `n_samples` attributes.

## Value

A data frame of class `c("ecodive_network_diff", "data.frame")`, one row
per pair of features, with columns:

- `feature1`, `feature2`: the pair, `feature1` first in `x`'s input
  table.

- `association_x`, `association_y`: their correlation in each group.

- `z`: Fisher's \\z\\ statistic.

- `p_value`, `p_adjusted`: its two-sided p-value, and that adjusted by
  `p_adjust` over all pairs.

- `reject`: whether `p_adjusted` is at most `alpha`.

Attributes: `method` (the method that made `x` and `y`, when they are
results), `alpha`, `p_adjust`, and `display` (the title; see
[ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md)).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the pairs
with `reject` as a network
([`plot.ecodive_network_diff()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_network.md)).

## Details

Each correlation is transformed with Fisher's \\z = \mathrm{atanh}(r)\\,
and the difference is divided by its standard error:

\$\$z = \frac{\mathrm{atanh}\\r_x - \mathrm{atanh}\\r_y}{\sqrt{1/(n_x -
3) + 1/(n_y - 3)}}\$\$

which is standard normal when the two population correlations are equal.
Pairs with \\\|r\| \ge 1\\ in either group, or a missing correlation,
get `NA` (with a warning for \\\|r\| \ge 1\\).

For [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)
this is an approximation: its correlations are medians over Dirichlet
draws of an estimate, not sample Pearson correlations, so their sampling
variance isn't exactly \\1/(n - 3)\\.
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) isn't
accepted (a concordance, not a correlation); neither are partial
correlations, whose test needs other degrees of freedom.

**Names in NetCoMi.** `p_adjust` is `adjust`; `alpha` is `alpha`.
NetCoMi's default adjustment is the local false discovery rate; here
only [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)'s
methods are offered. NetCoMi's resampling-based tests (`"permute"`,
`"discordant"`) aren't offered.

## See also

[`assoc_network()`](https://cmmr.github.io/ecodive/reference/assoc_network.md),
[`assoc_network_compare()`](https://cmmr.github.io/ecodive/reference/assoc_network_compare.md),
[`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md).

## Examples

``` r
    # Two groups of 60 samples; taxa 1 and 2 correlated only in the first
    set.seed(1)
    make <- function (linked) {
      basis <- matrix(rnorm(60 * 10), 60, 10)
      if (linked) basis[, 2] <- basis[, 1] + rnorm(60, sd = 0.3)
      matrix(rpois(600, exp(basis + 4)), 60, 10,
             dimnames = list(NULL, paste0('OTU', 1:10)))
    }
    rx <- sparcc(make(TRUE),  n_draws = 20)
    ry <- sparcc(make(FALSE), n_draws = 20)

    d <- assoc_network_diff(rx, ry)
    d
#> # A data.frame: 45 × 8
#> # Differential Associations of SparCC Correlation
#>    feature1 feature2 association_x association_y       z  p_value p_adjusted
#>    <chr>    <chr>            <dbl>         <dbl>   <dbl>    <dbl>      <dbl>
#>  1 OTU1     OTU2           0.862         -0.0492  7.21   5.42e-13   2.44e-11
#>  2 OTU1     OTU3          -0.076          0.212  -1.55       0.12      0.491
#>  3 OTU1     OTU4          -0.006         -0.0542  0.258     0.797      0.974
#>  4 OTU1     OTU5           0.021         -0.0325  0.285     0.775      0.974
#>  5 OTU1     OTU6           0.0552         0.202  -0.799     0.424      0.974
#>  6 OTU1     OTU7           0.0389        -0.282   1.75     0.0794      0.408
#>  7 OTU1     OTU8          -0.0152         0.0537 -0.368     0.713      0.974
#>  8 OTU1     OTU9          -0.0461        -0.0338 -0.0656    0.948      0.974
#>  9 OTU1     OTU10          0.0141         0.0413 -0.145     0.885      0.974
#> 10 OTU2     OTU3           0.00872        0.0493 -0.217     0.828      0.974
#> # … 35 more rows
#> # … 1 more column: reject <lgl>
#> # Use print(x, n = ...) to see more rows
#> # Attributes: method, alpha, p_adjust
    d[d$reject, ]
#> # A data.frame: 2 × 8
#> # Differential Associations of SparCC Correlation
#>    feature1 feature2 association_x association_y     z  p_value p_adjusted
#>    <chr>    <chr>            <dbl>         <dbl> <dbl>    <dbl>      <dbl>
#> 1  OTU1     OTU2             0.862       -0.0492  7.21 5.42e-13   2.44e-11
#> 43 OTU8     OTU9             0.312       -0.371   3.8  0.000146    0.00328
#> # … 1 more column: reject <lgl>
#> # Attributes: method, alpha, p_adjust
```
