# Abundance-based Coverage Estimator (ACE)

A non-parametric estimator of species richness that separates features
into abundant and rare groups.

## Usage

``` r
ace(counts, cutoff = 10L, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Also
  supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. See
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
  for optimizing large datasets.

- cutoff:

  The maximum number of observations to consider "rare". Default: `10`.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a special object
  class (e.g. `phyloseq`). Default: `1L`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

The ACE metric separates features into "abundant" and "rare" groups
based on a cutoff (usually 10 counts). It assumes that the presence of
abundant species is certain, while the true number of rare species must
be estimated.

**Equations:**

\$\$C\_{ace} = 1 - \frac{F_1}{X\_{rare}}\$\$

\$\$\gamma\_{ace}^2 = \max\left\[\frac{F\_{rare}
\sum\_{i=1}^{r}i(i-1)F_i}{C\_{ace}X\_{rare}(X\_{rare} - 1)} - 1,
0\right\]\$\$

\$\$D\_{ace} = F\_{abund} + \frac{F\_{rare}}{C\_{ace}} +
\frac{F_1}{C\_{ace}}\gamma\_{ace}^2\$\$

Where:

- \\r\\ : Rare cutoff (default 10). Features with \\\le r\\ counts are
  considered rare.

- \\F_i\\ : Number of features with exactly \\i\\ counts.

- \\F_1\\ : Number of features where \\X_i = 1\\ (singletons).

- \\F\_{rare}\\ : Number of rare features where \\X_i \le r\\.

- \\F\_{abund}\\ : Number of abundant features where \\X_i \> r\\.

- \\X\_{rare}\\ : Total counts belonging to rare features.

- \\C\_{ace}\\ : The sample abundance coverage estimator.

- \\\gamma\_{ace}^2\\ : The estimated coefficient of variation.

**Parameter: cutoff** The integer threshold distinguishing rare from
abundant species. Standard practice is to use 10.

## References

Chao, A., & Lee, S. M. (1992). Estimating the number of classes via
sample coverage. *Journal of the American Statistical Association*,
87(417), 210-217.
[doi:10.1080/01621459.1992.10475194](https://doi.org/10.1080/01621459.1992.10475194)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Richness metrics:
[`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md),
[`margalef()`](https://cmmr.github.io/ecodive/reference/margalef.md),
[`menhinick()`](https://cmmr.github.io/ecodive/reference/menhinick.md),
[`observed()`](https://cmmr.github.io/ecodive/reference/observed.md),
[`squares()`](https://cmmr.github.io/ecodive/reference/squares.md)

## Examples

``` r
    ace(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>    5.0    8.9    6.0    NaN 
```
