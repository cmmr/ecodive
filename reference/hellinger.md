# Hellinger distance

A distance metric related to the Bhattacharyya distance, often used for
community data with many zeros.

## Usage

``` r
hellinger(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
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

The Hellinger distance is defined as:
\$\$\sqrt{\sum\_{i=1}^{n}(\sqrt{P_i} - \sqrt{Q_i})^{2}}\$\$

Where:

- \\P_i\\, \\Q_i\\ : Proportional abundances of the \\i\\-th feature.

- \\n\\ : The number of features.

Base R Equivalent:

    p <- x / sum(x)
    q <- y / sum(y)
    sqrt(sum((sqrt(p) - sqrt(q)) ^ 2))

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

Rao, C. R. (1995). A review of canonical coordinates and an alternative
to correspondence analysis using Hellinger distance. *Qüestiió*, 19,
23-63.

Hellinger, E. (1909). Neue Begründung der Theorie quadratischer Formen
von unendlichvielen Veränderlichen. *Journal für die reine und
angewandte Mathematik*, 136, 210–271.
[doi:10.1515/crll.1909.136.210](https://doi.org/10.1515/crll.1909.136.210)

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
    hellinger(ex_counts)
#>          Saliva      Gums      Nose
#> Gums  0.4867106                    
#> Nose  1.2935044 1.2732200          
#> Stool 1.3170808 1.3273191 1.3417468
```
