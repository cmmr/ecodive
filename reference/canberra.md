# Canberra distance

A weighted version of the Manhattan distance, sensitive to differences
when both values are small.

## Usage

``` r
canberra(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Also
  supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. See
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
  for optimizing large datasets.

- norm:

  Normalize the incoming counts. Options are:

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'binary'`: Unweighted presence/absence (each count is either 0 or
    1).

  - `'clr'`: Centered log ratio.

  - `'none'`: No transformation.

  Default: `'percent'`, which is the expected input for these formulas.

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

The Canberra distance is defined as: \$\$\sum\_{i=1}^{n} \frac{\|P_i -
Q_i\|}{P_i + Q_i}\$\$

Where:

- \\P_i\\, \\Q_i\\ : Proportional abundances of the \\i\\-th feature.

- \\n\\ : The number of features.

Base R Equivalent:

    sum(abs(x-y) / (x+y))

## References

Lance, G. N., & Williams, W. T. (1966). Computer programs for
hierarchical polythetic classification ("similarity analyses"). *The
Computer Journal*, 9(1), 60-64.
[doi:10.1093/comjnl/9.1.60](https://doi.org/10.1093/comjnl/9.1.60)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

Other Abundance metrics:
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`bhattacharyya()`](https://cmmr.github.io/ecodive/reference/bhattacharyya.md),
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md),
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
    canberra(ex_counts)
#>         Saliva     Gums     Nose
#> Gums  2.558846                  
#> Nose  5.371936 5.285004         
#> Stool 5.975282 5.135293 4.934279
```
