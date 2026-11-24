# Chao1 Richness Estimator

A non-parametric estimator of the lower bound of species richness.

## Usage

``` r
chao1(counts, bias_corrected = TRUE, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- bias_corrected:

  For Chao1 only: use the bias-corrected form (the default, as
  scikit-bio and QIIME 2), or the classic form when a sample has both
  singletons and doubletons (`FALSE`). See `chao1()`.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Details

The Chao1 estimator uses the ratio of singletons to doubletons to
estimate the number of missing species. By default it is the
bias-corrected form (Chao 1987), as in scikit-bio and QIIME 2: \$\$n +
\frac{F_1 (F_1 - 1)}{2 (F_2 + 1)}\$\$

With `bias_corrected = FALSE`, it is the classic form (Chao 1984):
\$\$n + \frac{(F_1)^2}{2 F_2}\$\$ which is infinite when there are no
doubletons. So, as scikit-bio's `chao1(bias_corrected = False)`, the
classic form is used only when \\F_1 \> 0\\ and \\F_2 \> 0\\, and the
bias-corrected form otherwise.

Where:

- \\n\\ : The number of observed features.

- \\F_1\\ : Number of features observed once (singletons).

- \\F_2\\ : Number of features observed twice (doubletons).

`vegan::estimateR()` and iNEXT multiply the bias-corrected term by
\\(N - 1)/N\\, where \\N\\ is the number of observations, so their
estimates are slightly smaller than these for shallow samples with
singletons.

ecodive before version 2.4.0 used the classic form without the fallback,
which gave `Inf` without doubletons and `NaN` without singletons and
doubletons.

**Base R Equivalent:**

    x  <- ex_counts[1,]
    F1 <- sum(x == 1)
    F2 <- sum(x == 2)
    sum(x > 0) + F1 * (F1 - 1) / (2 * (F2 + 1))

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

Chao, A. (1984). Nonparametric estimation of the number of classes in a
population. *Scandinavian Journal of Statistics*, 11, 265-270.

Chao, A. (1987). Estimating the population size for capture-recapture
data with unequal catchability. *Biometrics*, 43(4), 783-791.
[doi:10.2307/2531532](https://doi.org/10.2307/2531532)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Richness metrics:
[`ace()`](https://cmmr.github.io/ecodive/reference/ace.md),
[`margalef()`](https://cmmr.github.io/ecodive/reference/margalef.md),
[`menhinick()`](https://cmmr.github.io/ecodive/reference/menhinick.md),
[`observed()`](https://cmmr.github.io/ecodive/reference/observed.md),
[`squares()`](https://cmmr.github.io/ecodive/reference/squares.md)

## Examples

``` r
    chao1(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>      4      6      6     11 
    
    # The classic form
    chao1(ex_counts, bias_corrected = FALSE)
#> Saliva   Gums   Nose  Stool 
#>    4.5    6.0    6.0   11.0 
```
