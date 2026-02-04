# Minkowski distance

A generalized metric that includes Euclidean and Manhattan distance as
special cases.

## Usage

``` r
minkowski(
  counts,
  norm = "none",
  power = 1.5,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- norm:

  Normalize the incoming counts. Options are:

  - `'none'`: No transformation.

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'binary'`: Unweighted presence/absence (each count is either 0 or
    1).

  - `'clr'`: Centered log ratio.

  Default: `'none'`.

- power:

  Scaling factor for the magnitude of differences between communities
  (\\p\\). Default: `1.5`

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

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
(X_i - Y_i)^p}\$\$

Where:

- \\X_i\\, \\Y_i\\ : Absolute abundances of the \\i\\-th feature.

- \\n\\ : The number of features.

- \\p\\ : The geometry of the space (power parameter).

**Parameter: power**

The `power` parameter (default 1.5) determines the value of \\p\\ in the
equation.

**Special Cases**

- **Manhattan distance**: When \\p = 1\\, the formula reduces to the sum
  of absolute differences.

- **Euclidean distance**: When \\p = 2\\, the formula reduces to the
  standard straight-line distance.

- **Chebyshev distance**: When \\p \to \infty\\, the formula reduces to
  the maximum absolute difference.

Base R Equivalent:

    p <- 1.5
    sum(abs(x - y)^p) ^ (1/p)

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
#> Gums   637.8205                    
#> Nose   646.1300  983.5644          
#> Stool  654.8633 1001.5543  858.2296
```
