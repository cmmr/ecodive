# Associations Between Features

Runs any of ecodive's measures of association between the features of a
count table by name, as
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
does for ordinations. Every one returns a `dist` over the features.

## Usage

``` r
feature_assoc(counts, method, ...)
```

## Arguments

- counts:

  A numeric matrix of counts (samples \\\times\\ features), or any of
  the objects listed under **Input Types**.

- method:

  The name of an association measure: `"vlr"`, `"rho_p"`, `"phi_s"` or
  `"sparcc"`. Flexible matching is supported (see below). Programmatic
  access via
  [`list_feature_assocs()`](https://cmmr.github.io/ecodive/reference/feature_assocs.md).

- ...:

  Further arguments for the method's function, such as `margin`,
  `pseudocount`, `clust` and `cpus`,
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md)'s and
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md)'s
  `denominator`, or
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)'s
  `n_draws`, `permutations` and `seed`. Each method keeps its own
  defaults. `seed` is ignored by the methods that have no random parts
  (all but
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)).
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md) takes
  no `pseudocount`: it handles zeros itself.

## Value

The method's result: a `dist` over the features, of class
`c("ecodive_dist", "dist")`, with attribute `method` set to the method's
ID. `feature_assoc(counts, 'sparcc', ...)` is the same as
`sparcc(counts, ...)`.

## Details

**List of Association Measures**

|  |  |  |  |
|----|----|----|----|
| Option / Function Name | Measure | Kind | Range |
| `vlr` | Variance of log ratios | dissimilarity | \\\[0, \infty)\\ |
| `rho_p` | Proportionality \\\rho_p\\ | similarity | \\\[-1, 1\]\\ |
| `phi_s` | Symmetric proportionality \\\phi_s\\ | dissimilarity | \\\[0, \infty)\\ |
| `sparcc` | SparCC correlation | similarity | \\\[-1, 1\]\\ |

The similarities
([`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
[`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)) can be
made into networks with
[`assoc_network()`](https://cmmr.github.io/ecodive/reference/assoc_network.md),
or listed with
[`assoc_pairs()`](https://cmmr.github.io/ecodive/reference/assoc_pairs.md).

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters.
Other names also work: `'rho'` or `'propr'` for `rho_p`; `'phi'` or
`'phs'` for `phi_s`; `'FastSpar'` for `sparcc`; `'log ratio variance'`
for `vlr`. `'proportionality'` matches both `rho_p` and `phi_s`, and so
is an error.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the list of methods.

## Input Types

The `counts` parameter is designed to accept a simple numeric matrix,
but seamlessly supports objects from the following biological data
packages:

- `phyloseq`

- `rbiom`

- `SummarizedExperiment`

- `TreeSummarizedExperiment`

For large datasets, standard matrix operations may be slow. See
[`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
for details on using optimized formats (e.g. sparse matrices) and
parallel processing.

## See also

[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md),
[`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md),
[`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md),
[`list_feature_assocs()`](https://cmmr.github.io/ecodive/reference/feature_assocs.md);
[`assoc_network()`](https://cmmr.github.io/ecodive/reference/assoc_network.md).

## Examples

``` r
    r <- feature_assoc(ex_counts, 'rho_p', pseudocount = 1)
    r
#> # A dist: 6 features (lower triangle)
#> # Proportionality (rho_p)
#>                   Streptococcus Bacteroides Corynebacterium Haemophilus
#> Bacteroides              -0.354
#> Corynebacterium          -0.753     -0.284
#> Haemophilus               0.887     -0.0663          -0.887
#> Propionibacterium        -0.491     -0.463            0.915      -0.668
#> Staphylococcus           -0.664     -0.253            0.955      -0.833
#> # … 1 more column: Propionibacterium
#> # Attributes: method

    # The same arguments for every log-ratio measure.
    lapply(c('vlr', 'rho_p', 'phi_s'), function (m)
      range(feature_assoc(ex_counts, m, pseudocount = 1)) )
#> [[1]]
#> [1]  0.3212435 23.7749710
#> 
#> [[2]]
#> [1] -0.8874735  0.9549243
#> 
#> [[3]]
#> [1]  0.0230575 16.7735862
#> 
```
