# Fisher's Alpha

A parametric diversity index assuming species abundances follow a
log-series distribution.

## Usage

``` r
fisher(counts, digits = 3L, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Also
  supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. See
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
  for optimizing large datasets.

- digits:

  Precision of the returned values, in number of decimal places. E.g.
  the default `digits=3` could return `6.392`.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a special object
  class (e.g. `phyloseq`). Default: `1L`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

Fisher's Alpha (\\\alpha\\) is the parameter in the equation:
\$\$\frac{n}{\alpha} = \ln{\left(1 + \frac{X_T}{\alpha}\right)}\$\$

Where:

- \\n\\ : The number of features.

- \\X_T\\ : Total of all counts (sequencing depth).

The value of \\\alpha\\ is solved for iteratively.

**Parameter: digits** The precision (number of decimal places) to use
when solving the equation.

## References

Fisher, R. A., Corbet, A. S., & Williams, C. B. (1943). The relation
between the number of species and the number of individuals in a random
sample of an animal population. *Journal of Animal Ecology*, 12, 42-58.
[doi:10.2307/1411](https://doi.org/10.2307/1411)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Diversity metrics:
[`brillouin()`](https://cmmr.github.io/ecodive/reference/brillouin.md),
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md),
[`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md)

## Examples

``` r
    fisher(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>  0.635  0.700  0.847  0.744 
```
