# Minkowski distance

A generalized metric that includes Euclidean and Manhattan distance as
special cases.

## Usage

``` r
minkowski(
  counts,
  norm = "percent",
  power = 1.5,
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

- norm:

  Normalize the incoming counts. Options are: \* `'percent'`: Relative
  abundance (sample abundances sum to 1). \* `'binary'`: Unweighted
  presence/absence (each count is either 0 or 1). \* `'clr'`: Centered
  log ratio. \* `'none'`: No transformation.

         Default: `'percent'`, which is the expected input for these formulas.

- power:

  Scaling factor for the magnitude of differences between communities
  (\\p\\). Default: `1.5`

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

The Minkowski distance is defined as: \$\$\sqrt\[p\]{\sum\_{i=1}^{n}
(P_i - Q_i)^p}\$\$

Where:

- \\P_i\\, \\Q_i\\ : Proportional abundances of the \\i\\-th feature.

- \\n\\ : The number of features.

- \\p\\ : The geometry of the space (power parameter).

**Parameter: power** The `power` parameter (default 1.5) determines the
value of \\p\\ in the equation.

Base R Equivalent:

    p <- 1.5
    sum(abs(x - y)^p) ^ (1/p)

## References

Deza, M. M., & Deza, E. (2009). Encyclopedia of distances. Springer.

Minkowski, H. (1896). *Geometrie der Zahlen*. Teubner.

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

Other Abundance metrics:
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
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
    minkowski(ex_counts, power = 2) # Equivalent to Euclidean
#>          Saliva      Gums      Nose
#> Gums  0.6003499                    
#> Nose  0.9099609 1.0628029          
#> Stool 1.2103914 1.3362601 1.1575500
```
