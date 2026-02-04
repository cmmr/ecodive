# Clark's divergence distance

Also known as the coefficient of divergence.

## Usage

``` r
clark(counts, norm = "none", margin = 1L, pairs = NULL, cpus = n_cpus())
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

Clark's divergence distance is defined as:
\$\$\sqrt{\sum\_{i=1}^{n}\left(\frac{X_i - Y_i}{X_i +
Y_i}\right)^{2}}\$\$

Where:

- \\X_i\\, \\Y_i\\ : Absolute abundances of the \\i\\-th feature.

- \\n\\ : The number of features.

Base R Equivalent:

    sqrt(sum((abs(x - y) / (x + y)) ^ 2))

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

Clark, P. J. (1952). An extension of the coefficient of divergence for
use with multiple characters. *Copeia*, 1952(2), 61-64.
[doi:10.2307/1438598](https://doi.org/10.2307/1438598)

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
    clark(ex_counts)
#>         Saliva     Gums     Nose
#> Gums  1.291899                  
#> Nose  2.125944 2.210332         
#> Stool 2.437340 2.219050 2.215163
```
