# Symmetric Proportionality

The symmetric proportionality measure \\\phi_s\\ of Quinn et al. (2017),
propr's `"phs"`: the variance of the log ratio of two features relative
to that of their log product. `0` means perfectly proportional.

## Usage

``` r
phi_s(
  counts,
  margin = 1L,
  denominator = NULL,
  pseudocount = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Counts,
  proportions, and the output of
  [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
  are all accepted: each sample's scale cancels out of the log ratios.
  Pairs are taken between the **features** (the columns when
  `margin = 1`), and the variance is over the samples.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- denominator:

  The reference that each feature's log abundance is taken relative to,
  for the variances \\\mathrm{var}(z_i)\\. `NULL` (the default) uses the
  centered log ratio (clr): each sample's geometric mean of all
  features, as propr's default `ivar = "clr"`. Otherwise features, as
  names, (1-based) indices, or a logical vector with one value per
  feature: their geometric mean. One feature gives the additive log
  ratio (alr, as `alr(denominator =)`), and that feature is then left
  out of the result (its alr is 0 in every sample); with several, every
  feature is kept.

- pseudocount:

  Added to every value before taking logs. `NULL` (the default) adds
  nothing when `counts` has no zeros, and otherwise uses half the
  smallest non-zero value, with a warning. `0` adds nothing, and is an
  error when `counts` has zeros. See **Zeros**.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A `dist` object over the features, as
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) returns,
with `method` `'phi_s'`. It is a dissimilarity, so it can be given as
`dm` to [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md) or
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

## Details

With \\z\\ the log-ratio-transformed abundances (clr or alr),
\$\$\phi\_{s,ij} = \frac{\mathrm{var}(z_i - z_j)}{\mathrm{var}(z_i +
z_j)} = \frac{1 - \rho\_{p,ij}}{1 + \rho\_{p,ij}}\$\$ where
\\\mathrm{var}(z_i - z_j)\\ is
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md) and
\\\rho_p\\ is
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md). It lies
in \\\[0, \infty)\\: 0 for proportional features, 1 for unrelated ones,
and infinite for reciprocal ones (\\z_i = -z_j\\). It is computed as
\\\mathrm{vlr}\_{ij} / (2\mathrm{var}(z_i) + 2\mathrm{var}(z_j) -
\mathrm{vlr}\_{ij})\\, since \\\mathrm{var}(a + b) = 2\mathrm{var}(a) +
2\mathrm{var}(b) - \mathrm{var}(a - b)\\. Unlike
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md), it depends
on the `denominator`. Near-reciprocal pairs (such as the two features of
a two-feature `denominator`, whose \\z\\ are exact opposites) have
\\\mathrm{var}(z_i + z_j)\\ near 0, where rounding dominates: their
value is `Inf` or merely very large.

Lovell et al.'s (2015) original \\\phi = \mathrm{vlr}\_{ij} /
\mathrm{var}(z_i)\\ is not symmetric (\\\phi\_{ij} \ne \phi\_{ji}\\), so
it can't be a `dist`; Quinn et al. (2017) replaced it with \\\phi_s\\.
This is the measure that propr computes with
`propr(counts, metric = "phs")` (see **Zeros** for its zero handling).

## Zeros

Logs of zero are undefined, so zeros need replacing first. With
`pseudocount = NULL` (the default), no pseudocount is used when `counts`
has no zeros; otherwise half the smallest non-zero value is added to
every value, with a warning. The recommended input is a table with its
zeros replaced by
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
which has none and is used as it is.

propr (Quinn et al. 2017) replaces only the zeros, by 1, which is not
`pseudocount = 1` (that adds 1 to every value). To reproduce it, write
`x[x == 0] <- 1` and then call these functions on `x`.

Pseudocounts are not neutral: in sparse tables the many zeros, all given
the same value, can dominate the variances. Filtering out rare features
first helps.

## Size

The result has \\D(D - 1)/2\\ values: 20,000 features give 2e8 pairs,
1.6 GB. Filter out rare features first for large tables. The table
itself stays sparse.

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

## References

Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., &
Bähler, J. (2015). Proportionality: a valid alternative to correlation
for relative data. *PLoS Computational Biology*, 11, e1004075.

Quinn, T. P., Richardson, M. F., Lovell, D., & Crowley, T. M. (2017).
propr: an R-package for identifying proportionally abundant features
using compositional data analysis. *Scientific Reports*, 7, 16252.

## See also

[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md),
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md).

## Examples

``` r
    x <- replace_zeros(ex_counts, 'mult')
    phi_s(x)
#>                   Streptococcus  Bacteroides Corynebacterium  Haemophilus
#> Bacteroides          3.03059414                                          
#> Corynebacterium      0.84332117   4.65087656                             
#> Haemophilus          0.07245552   1.54794798      1.93602839             
#> Propionibacterium   29.69577708   0.53014503      0.66220067 118.51898728
#> Staphylococcus       2.61515860   2.40122280      0.21092244   3.30740984
#>                   Propionibacterium
#> Bacteroides                        
#> Corynebacterium                    
#> Haemophilus                        
#> Propionibacterium                  
#> Staphylococcus           0.30095796

    # The same as (1 - rho_p) / (1 + rho_p)
    rho <- rho_p(x)
    all.equal(c(phi_s(x)), c((1 - rho) / (1 + rho)))
#> [1] TRUE
```
