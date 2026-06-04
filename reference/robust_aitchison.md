# Robust Aitchison distance

Calculates the pairwise Robust Aitchison distance for compositional
data. This method is specifically engineered for sparse datasets - such
as microbiome OTU/ASV tables - by calculating distances based only on
observed positive abundances, avoiding the need for pseudo-counts.

## Usage

``` r
robust_aitchison(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

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

The Robust Aitchison distance is defined as: \$\$\sqrt{\sum\_{i=1}^{n}
(X^\*\_i - Y^\*\_i)^2}\$\$

Where:

- \\X^\*\_i\\, \\Y^\*\_i\\ : The rclr-transformed counts for the
  \\i\\-th feature. For a given sample \\X\\, \\X^\*\_i = \ln(X_i) -
  X_L\\ if \\X_i \> 0\\, and \\0\\ otherwise.

- \\X_L\\, \\Y_L\\ : Mean log of strictly positive abundances. \\X_L =
  \frac{1}{\|P_X\|}\sum\_{j \in P_X} \ln{X_j}\\, where \\P_X\\ is the
  set of indices where \\X \> 0\\.

- \\\|P_X\|\\, \\\|P_Y\|\\ : The number of strictly positive features in
  the respective samples.

- \\n\\ : The total number of features.

Base R Equivalent:

    x <- ifelse(x > 0, log(x) - mean(log(x[x > 0])), 0)
    y <- ifelse(y > 0, log(y) - mean(log(y[y > 0])), 0)
    sqrt(sum((x-y)^2)) # Euclidean distance

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

Martino, C., Morton, J. T., Marotz, C. A., Thompson, L. R., Tripathi,
A., Knight, R., and Zengler, K. (2019). A novel sparse compositional
technique reveals microbial perturbations. *mSystems*, 4(1), e00016-19.
[doi:10.1128/mSystems.00016-19](https://doi.org/10.1128/mSystems.00016-19)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)

[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
`vignette('bdiv_guide')`

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
    robust_aitchison(ex_counts)
#>         Saliva     Gums     Nose
#> Gums  3.282546                  
#> Nose  8.345744 9.711832         
#> Stool 9.460092 9.421025 9.830042
```
