# Strongly Associated Pairs of Features

Lists the pairs of features whose correlation is at least `threshold` in
absolute value, and optionally significant, from a
[`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md) or
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) result:
one row per pair, strongest first.

## Usage

``` r
assoc_pairs(x, threshold = 0.3, alpha = NULL, p_adjust = "BH")
```

## Arguments

- x:

  Correlations between features: a result of
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md) or
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md), or a
  symmetric matrix of correlations with 1 on the diagonal (such as from
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)).

- threshold:

  Keep the pairs whose correlation is at least this in absolute value
  (\\\|r\| \ge\\ `threshold`), both positive and negative ones. Default
  `0.3`, the SparCC paper's cutoff for a network. `0` keeps every pair.

- alpha:

  Also keep only the pairs whose adjusted p-value is at most this.
  `NULL` (the default) for no p-value filter. Needs p-values: a
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)
  result made with `permutations`.

- p_adjust:

  How to adjust the p-values for testing every pair: one of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html).
  Default `'BH'` (Benjamini & Hochberg's false discovery rate). The
  adjustment is over all the pairs of `x`, not just those kept. Used
  only when `x` has p-values.

## Value

A data frame of class `c("ecodive_assoc_pairs", "data.frame")`, one row
per pair kept, ordered from the largest \\\|r\|\\ to the smallest, with
columns:

- `from`, `to`: the two features' names (their positions, as text, when
  the features have no names), `from` being the one that comes first in
  the input table.

- `association`: their correlation, with its sign.

- `p_value`, `p_adjusted`: when `x` has p-values, the pair's p-value and
  its adjustment by `p_adjust` over all pairs.

Attributes: `method` (the method that made `x`, such as `"sparcc"`;
absent for a matrix), `threshold`, and, when used, `alpha` and
`p_adjust`.

With no pair kept, the data frame has no rows.

## Details

Pairs whose correlation is `NA` or `NaN` are never kept. `threshold` and
`alpha` both apply when both are given; use `threshold = 0` to filter on
p-values alone. The pairs kept are the edges of a correlation network
with this cutoff.

**Enough permutations.** A permutation p-value is at least \\1 / (B +
1)\\ for \\B\\ null datasets, and adjusting it for testing every pair
multiplies it by up to the number of pairs \\m\\. So with too few
permutations no pair can pass `alpha`, however strong. With
Benjamini-Hochberg, \\k\\ pairs at the smallest p-value pass only when
\\m / (k (B + 1)) \le\\ `alpha`: 20 features (190 pairs) and one strong
pair need about 3,800.

Only correlations are accepted.
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md) and
[`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) are
dissimilarities (smaller means more proportional), so "above a cutoff"
doesn't apply to them: use
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md), their
correlation-scale counterpart.

## See also

[`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md),
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Examples

``` r
    # 60 samples of 20 taxa, with five correlated pairs
    set.seed(1)
    basis <- matrix(rnorm(60 * 20), 60, 20)
    for (i in c(1, 3, 5, 7, 9))
      basis[, i + 1] <- basis[, i] + rnorm(60, sd = 0.5)
    counts <- matrix(
      rpois(1200, exp(basis + 4)), 60, 20,
      dimnames = list(NULL, paste0('OTU', 1:20)) )

    r <- sparcc(counts, n_draws = 20)
    assoc_pairs(r)
#> # A data.frame: 11 × 3
#> # SparCC Correlation: pairs with |r| >= 0.3
#>    from  to    association
#>    <chr> <chr>       <dbl>
#>  1 OTU5  OTU6        0.885
#>  2 OTU9  OTU10       0.883
#>  3 OTU7  OTU8        0.871
#>  4 OTU3  OTU4        0.837
#>  5 OTU1  OTU2        0.725
#>  6 OTU14 OTU18      -0.45
#>  7 OTU2  OTU16       0.38
#>  8 OTU11 OTU13       0.365
#>  9 OTU4  OTU20       0.362
#> 10 OTU14 OTU16       0.355
#> 11 OTU3  OTU20       0.339
#> # Attributes: method, threshold

    # Only the strongest
    assoc_pairs(r, threshold = 0.5)
#> # A data.frame: 5 × 3
#> # SparCC Correlation: pairs with |r| >= 0.5
#>   from  to    association
#>   <chr> <chr>       <dbl>
#> 1 OTU5  OTU6        0.885
#> 2 OTU9  OTU10       0.883
#> 3 OTU7  OTU8        0.871
#> 4 OTU3  OTU4        0.837
#> 5 OTU1  OTU2        0.725
#> # Attributes: method, threshold

    # Significant pairs, from 999 null datasets
    r <- sparcc(counts, n_draws = 10, permutations = 999)
    assoc_pairs(r, alpha = 0.05)
#> # A data.frame: 5 × 5
#> # SparCC Correlation: pairs with |r| >= 0.3 and BH-adjusted p <= 0.05
#>   from  to    association p_value p_adjusted
#>   <chr> <chr>       <dbl>   <dbl>      <dbl>
#> 1 OTU5  OTU6        0.893   0.001      0.038
#> 2 OTU9  OTU10       0.883   0.001      0.038
#> 3 OTU7  OTU8        0.874   0.001      0.038
#> 4 OTU3  OTU4        0.839   0.001      0.038
#> 5 OTU1  OTU2        0.729   0.001      0.038
#> # Attributes: method, threshold, alpha, p_adjust

    # Correlations from elsewhere
    assoc_pairs(cor(log(counts + 1)), threshold = 0.5)
#> # A data.frame: 5 × 3
#> # Correlations: pairs with |r| >= 0.5
#>   from  to    association
#>   <chr> <chr>       <dbl>
#> 1 OTU5  OTU6        0.915
#> 2 OTU9  OTU10       0.911
#> 3 OTU7  OTU8        0.888
#> 4 OTU3  OTU4        0.853
#> 5 OTU1  OTU2        0.767
#> # Attributes: threshold
```
