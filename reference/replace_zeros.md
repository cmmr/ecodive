# Zero Replacement

Replaces the zeros of a count table with small positive values, so that
log ratios can be taken without a pseudocount. Runs any of ecodive's
zero replacement methods by name, as
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
does for ordinations.

## Usage

``` r
replace_zeros(counts, method, ...)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Counts
  or proportions.

- method:

  The name of a zero replacement method: `"mult"`. Flexible matching is
  supported (see below). Programmatic access via
  [`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md).
  Required: there is no default, because the choice of method affects
  every downstream log ratio.

- ...:

  Further arguments for the method's function, such as `margin` and
  [`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md)'s
  `delta`. Each method keeps its own defaults.

## Value

The method's result: a dense numeric `matrix` of proportions with no
zeros, in the orientation of `counts` (samples \\\times\\ features for
`margin = 1`, features \\\times\\ samples for `margin = 2`), with its
names. Each sample sums to 1. Attributes `method` (the method's id) and
the method's parameters, ending with `n_replaced`, the number of zeros
replaced. `replace_zeros(counts, 'mult', ...)` is the same as
`replace_zeros_mult(counts, ...)`.

## Details

Logs of zero are undefined, so every log-ratio method, such as
`norm = 'clr'`,
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md) and
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md), needs the
zeros dealt with. A `pseudocount` added to every value changes the
ratios between the non-zero parts, and acts as an arbitrary weighting. A
zero replacement instead gives each zero a small value and scales the
sample's non-zero parts down to make room for it, keeping their ratios.

The result has no zeros, so it can be given to any of those functions
without a pseudocount (and without the warning that `pseudocount = NULL`
gives for zeros). The result is in proportions: every log-ratio method
is scale invariant, so nothing is lost.

Conclusions that depend on how zeros were treated are fragile. Repeating
an analysis with each method of `list_zero_replacements('id')` is a
simple check.

**List of Zero Replacement Methods**

|  |  |
|----|----|
| Option / Function Name | Method Name |
| `mult` | Multiplicative replacement, with a fixed value \\\delta\\ |

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters.
E.g. `method = 'MULT'`, `'multiplicative'` or
`'multiplicative replacement'` selects `mult`. Common other names also
work: `'multRepl'` (from zCompositions) and `'multi_replace'` (from
scikit-bio) for `mult`.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the list of methods.

**Same names in other packages**

None known. zCompositions' `multRepl()` and `cmultRepl()`, and
scikit-bio's `multi_replace()`, do the same job under other names.

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

Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V.
(2003). Dealing with zeros and missing values in compositional data sets
using nonparametric imputation. *Mathematical Geology*, 35, 253–278.
[doi:10.1023/A:1023866030544](https://doi.org/10.1023/A%3A1023866030544)

Lubbe, S., Filzmoser, P., & Templ, M. (2021). Comparison of zero
replacement strategies for compositional data with large numbers of
zeros. *Chemometrics and Intelligent Laboratory Systems*, 210, 104248.
[doi:10.1016/j.chemolab.2021.104248](https://doi.org/10.1016/j.chemolab.2021.104248)

## See also

[`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md),
[`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md);
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) and
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
(`'clr'`), which take the result.

## Examples

``` r
    x <- replace_zeros(ex_counts, 'mult')
    x[1:3, 1:4]
#>        Streptococcus Bacteroides Corynebacterium Haemophilus
#> Saliva    0.44347826 0.005475040      0.02777778 0.492753623
#> Gums      0.87017181 0.004389265      0.02777778 0.095466516
#> Nose      0.02176063 0.001978239      0.49258160 0.001978239
    attr(x, 'n_replaced')
#> [1] 4

    # No zeros left, so no pseudocount is needed.
    aitchison(x)
#> # A dist: 4 samples (lower triangle)
#> # Aitchison Distance
#>       Saliva  Gums Nose
#> Gums    3.12
#> Nose    8.58  9.54
#> Stool   9.89 10.1  9.63
#> # Attributes: method

    # The same analysis with every method.
    lapply(list_zero_replacements('id'), function (m) {
      aitchison(replace_zeros(ex_counts, m))
    })
#> [[1]]
#> # A dist: 4 samples (lower triangle)
#> # Aitchison Distance
#>       Saliva  Gums Nose
#> Gums    3.12
#> Nose    8.58  9.54
#> Stool   9.89 10.1  9.63
#> # Attributes: method
#> 
```
