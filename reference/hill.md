# Hill Numbers

The effective number of species of order \\q\\: how many equally
abundant features would give the same diversity.

## Usage

``` r
hill(counts, q = 2, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- q:

  The order of a Hill number, or of Rényi's or Tsallis's entropy: a
  single number of at least 0, where `Inf` is allowed. Only used by
  `hill()`,
  [`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) and
  [`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md).
  Default: `2`, as scikit-bio (which calls it `order`).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Details

The Hill number of order \\q\\ is defined as: \$\${}^qD = \left(\sum\_{i
= 1}^{n} P_i^q\right)^{1/(1-q)}\$\$

with its limits \\{}^1D = e^H\\ (\\H\\: Shannon's entropy) and
\\{}^\infty D = 1 / \max(P_i)\\.

Where:

- \\n\\ : The number of features.

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

The order sets how much weight abundant features get. The Hill numbers
are a profile that unifies the common indices:

- \\q = 0\\: richness,
  [`observed()`](https://cmmr.github.io/ecodive/reference/observed.md).

- \\q = 1\\: the exponential of Shannon's entropy, `exp(shannon())`,
  which weighs every feature by its frequency (Jost 2006).

- \\q = 2\\: the inverse Simpson index,
  [`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md).

- \\q = \infty\\: the inverse Berger-Parker index, `1 / berger()`.

The default order, 2, is scikit-bio's, which makes `hill()` equal to
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md);
give `q` for any other order.
[`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) is
`log(hill())`, and
[`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md) is
another function of it.

Observed Hill numbers grow with sequencing depth, most for small \\q\\;
compare samples at equal depth, for instance with
[`alpha_div_rarefied()`](https://cmmr.github.io/ecodive/reference/rarefied_diversity.md).

**Base R Equivalent:**

    x <- ex_counts[1,]
    q <- 2
    p <- x[x > 0] / sum(x)
    sum(p ** q) ** (1 / (1 - q))

## Same names in other packages

scikit-bio's `hill()` is the same (its `order` is `q`).
`vegan::renyi(x, scales, hill = TRUE)` gives a profile of Hill numbers
over several orders.

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

Hill, M. O. (1973). Diversity and evenness: a unifying notation and its
consequences. *Ecology*, 54, 427-432.
[doi:10.2307/1934352](https://doi.org/10.2307/1934352)

Jost, L. (2006). Entropy and diversity. *Oikos*, 113, 363-375.
[doi:10.1111/j.2006.0030-1299.14714.x](https://doi.org/10.1111/j.2006.0030-1299.14714.x)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md),
[`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Diversity metrics:
[`brillouin()`](https://cmmr.github.io/ecodive/reference/brillouin.md),
[`fisher()`](https://cmmr.github.io/ecodive/reference/fisher.md),
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md),
[`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md),
[`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md),
[`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md)

## Examples

``` r
    hill(ex_counts)
#>   Saliva     Gums     Nose    Stool 
#> 2.029446 1.233425 2.783607 1.013125 
    
    # Richness, exp(Shannon), inverse Simpson
    sapply(c(0, 1, 2), function (q) hill(ex_counts, q = q))
#>        [,1]     [,2]     [,3]
#> Saliva    4 2.098450 2.029446
#> Gums      5 1.443173 1.233425
#> Nose      6 3.133746 2.783607
#> Stool     5 1.049432 1.013125
```
