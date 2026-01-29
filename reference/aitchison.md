# Aitchison distance

Calculates the Euclidean distance between centered log-ratio (CLR)
transformed abundances.

## Usage

``` r
aitchison(
  counts,
  pseudocount = NULL,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Also
  supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. See
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
  for optimizing large datasets.

- pseudocount:

  The value to add to all counts in `counts` to prevent taking `log(0)`
  for unobserved features. The default, `NULL`, selects the smallest
  non-zero value in `counts`.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a special object
  class (e.g. `phyloseq`). Default: `1L`

- pairs:

  Which combinations of samples should distances be calculated for? The
  default value (`NULL`) calculates all-vs-all. Provide a numeric or
  logical vector specifying positions in the distance matrix to
  calculate. See examples.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

The Aitchison distance is defined as: \$\$\sqrt{\sum\_{i=1}^{n}
\[(\ln{X_i} - X_L) - (\ln{Y_i} - Y_L)\]^2}\$\$

Where:

- \\X_i\\, \\Y_i\\ : Absolute counts for the \\i\\-th feature.

- \\X_L\\, \\Y_L\\ : Mean log of abundances. \\X_L =
  \frac{1}{n}\sum\_{i=1}^{n} \ln{X_i}\\.

- \\n\\ : The number of features.

**Parameter: pseudocount** Because the formula uses logarithms, zeros in
the data must be handled. The `pseudocount` argument adds a small value
to all counts prior to calculation.

Base R Equivalent:

    x <- log((x + pseudocount) / exp(mean(log(x + pseudocount))))
    y <- log((y + pseudocount) / exp(mean(log(y + pseudocount))))
    sqrt(sum((x-y)^2)) # Euclidean distance

## References

Aitchison, J. (1986). The statistical analysis of compositional data.
Chapman and Hall.
[doi:10.1007/978-94-009-4109-3](https://doi.org/10.1007/978-94-009-4109-3)

Aitchison, J. (1982). The statistical analysis of compositional data.
*Journal of the Royal Statistical Society: Series B (Methodological)*,
44(2), 139-160.
[doi:10.1111/j.2517-6161.1982.tb01195.x](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

Other Abundance metrics:
[`bhattacharyya()`](https://cmmr.github.io/ecodive/reference/bhattacharyya.md),
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md),
[`canberra()`](https://cmmr.github.io/ecodive/reference/canberra.md),
[`chebyshev()`](https://cmmr.github.io/ecodive/reference/chebyshev.md),
[`chord()`](https://cmmr.github.io/ecodive/reference/chord.md),
[`clark()`](https://cmmr.github.io/ecodive/reference/clark.md),
[`divergence()`](https://cmmr.github.io/ecodive/reference/divergence.md),
[`euclidean()`](https://cmmr.github.io/ecodive/reference/euclidean.md),
[`gower()`](https://cmmr.github.io/ecodive/reference/gower.md),
[`hellinger()`](https://cmmr.github.io/ecodive/reference/hellinger.md),
[`horn()`](https://cmmr.github.io/ecodive/reference/horn.md),
[`jensen()`](https://cmmr.github.io/ecodive/reference/jensen.md),
[`jsd()`](https://cmmr.github.io/ecodive/reference/jsd.md),
[`lorentzian()`](https://cmmr.github.io/ecodive/reference/lorentzian.md),
[`manhattan()`](https://cmmr.github.io/ecodive/reference/manhattan.md),
[`matusita()`](https://cmmr.github.io/ecodive/reference/matusita.md),
[`minkowski()`](https://cmmr.github.io/ecodive/reference/minkowski.md),
[`morisita()`](https://cmmr.github.io/ecodive/reference/morisita.md),
[`motyka()`](https://cmmr.github.io/ecodive/reference/motyka.md),
[`psym_chisq()`](https://cmmr.github.io/ecodive/reference/psym_chisq.md),
[`soergel()`](https://cmmr.github.io/ecodive/reference/soergel.md),
[`squared_chisq()`](https://cmmr.github.io/ecodive/reference/squared_chisq.md),
[`squared_chord()`](https://cmmr.github.io/ecodive/reference/squared_chord.md),
[`squared_euclidean()`](https://cmmr.github.io/ecodive/reference/squared_euclidean.md),
[`topsoe()`](https://cmmr.github.io/ecodive/reference/topsoe.md),
[`wave_hedges()`](https://cmmr.github.io/ecodive/reference/wave_hedges.md)

## Examples

``` r
    aitchison(ex_counts, pseudocount = 1)
#>         Saliva     Gums     Nose
#> Gums  1.748405                  
#> Nose  9.710741 9.862353         
#> Stool 8.245648 8.372391 9.409111
```
