# Aitchison distance

Calculates the Euclidean distance between centered log-ratio (CLR)
transformed abundances.

## Usage

``` r
aitchison(
  counts,
  margin = 1L,
  pseudocount = NULL,
  pairs = NULL,
  cpus = n_cpus()
)
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

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'`. Ignored for
  other normalization methods. See **Pseudocount** section.

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

Base R Equivalent:

    x <- log((x + pseudocount) / exp(mean(log(x + pseudocount))))
    y <- log((y + pseudocount) / exp(mean(log(y + pseudocount))))
    sqrt(sum((x-y)^2)) # Euclidean distance

## Pseudocount

Zeros are undefined in the Aitchison (CLR) transformation. If
`pseudocount` is `NULL` (the default) and zeros are detected, the
function uses half the minimum non-zero value (`min(x[x>0]) / 2`) and
issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`).

**Why this matters:** The choice of pseudocount is not neutral; it acts
as a weighting factor that can significantly distort downstream results,
especially for sparse datasets. See Gloor et al. (2017) and Kaul et al.
(2017) for open-access discussions on the mathematical implications, or
Costea et al. (2014) for the impact on community clustering.

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

Aitchison, J. (1986). The statistical analysis of compositional data.
Chapman and Hall.
[doi:10.1007/978-94-009-4109-3](https://doi.org/10.1007/978-94-009-4109-3)

Aitchison, J. (1982). The statistical analysis of compositional data.
*Journal of the Royal Statistical Society: Series B (Methodological)*,
44(2), 139-160.
[doi:10.1111/j.2517-6161.1982.tb01195.x](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x)

Costea, P. I., Zeller, G., Sunagawa, S., & Bork, P. (2014). A fair
comparison. *Nature Methods*, 11(4), 359.
[doi:10.1038/nmeth.2897](https://doi.org/10.1038/nmeth.2897)

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
(2017). Microbiome datasets are compositional: and this is not optional.
*Frontiers in Microbiology*, 8, 2224.
[doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)

Kaul, A., Mandal, S., Davidov, O., & Peddada, S. D. (2017). Analysis of
microbiome data in the presence of excess zeros. *Frontiers in
Microbiology*, 8, 2114.
[doi:10.3389/fmicb.2017.02114](https://doi.org/10.3389/fmicb.2017.02114)

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
