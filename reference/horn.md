# Horn-Morisita dissimilarity

A similarity index based on Simpson's diversity index, suitable for
abundance data.

## Usage

``` r
horn(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())
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

The Horn-Morisita dissimilarity is defined as: \$\$1 -
\frac{2\sum\_{i=1}^{n}P\_{i}Q\_{i}}{\sum\_{i=1}^{n}P_i^2 +
\sum\_{i=1}^{n}Q_i^2}\$\$

Where:

- \\P_i\\, \\Q_i\\ : Proportional abundances of the \\i\\-th feature.

- \\n\\ : The number of features.

Base R Equivalent:

    z <- sum(x^2) / sum(x)^2 + sum(y^2) / sum(y)^2
    1 - ((2 * sum(x * y)) / (z * sum(x) * sum(y)))

## References

Horn, H. S. (1966). Measurement of "overlap" in comparative ecological
studies. *The American Naturalist*, 100(914), 419-424.
[doi:10.1086/282436](https://doi.org/10.1086/282436)

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
    horn(ex_counts)
#>          Saliva      Gums      Nose
#> Gums  0.2765026                    
#> Nose  0.9718747 0.9654300          
#> Stool 0.9900374 0.9932115 0.9952692
```
