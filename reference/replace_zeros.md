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

  The name of a zero replacement method: `"mult"`, `"czm"`, `"gbm"`,
  `"sq"` or `"bl"`. Flexible matching is supported (see below).
  Programmatic access via
  [`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md).
  Required: there is no default, because the choice of method affects
  every downstream log ratio. `"gbm"` is the common choice for count
  tables (it is zCompositions' default).

- ...:

  Further arguments for the method's function, such as `margin`,
  [`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md)'s
  `delta`, or `frac`, `threshold` and `adjust` (see
  [`replace_zeros_gbm()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md)).
  Each method keeps its own defaults.

## Value

The method's result: a dense numeric `matrix` of proportions with no
zeros, in the orientation of `counts` (samples \\\times\\ features for
`margin = 1`, features \\\times\\ samples for `margin = 2`), with its
names. Each sample sums to 1. Attributes `method` (the method's id),
then the method's results, ending with `n_replaced`, the number of zeros
replaced. `replace_zeros(counts, 'gbm', ...)` is the same as
`replace_zeros_gbm(counts, ...)`.

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

`mult` gives every zero the same proportion, whatever the sample's
depth. The other methods are for counts: `czm` imputes a fraction of one
read in each sample, and the Bayesian-multiplicative methods (`gbm`,
`sq`, `bl`) also use how abundant the feature is in the other samples.
`gbm` is the common default for count tables.

Conclusions that depend on how zeros were treated are fragile. Repeating
an analysis with each method of `list_zero_replacements('id')` is a
simple check.

**List of Zero Replacement Methods**

|  |  |
|----|----|
| Option / Function Name | Method Name |
| `mult` | Multiplicative replacement, with a fixed value \\\delta\\ ([`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md)) |
| `czm` | Count zero multiplicative ([`replace_zeros_czm()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md)) |
| `gbm` | Geometric Bayesian-multiplicative ([`replace_zeros_gbm()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md)) |
| `sq` | Square-root Bayesian-multiplicative ([`replace_zeros_sq()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md)) |
| `bl` | Bayes-Laplace Bayesian-multiplicative ([`replace_zeros_bl()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md)) |

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters.
E.g. `method = 'MULT'`, `'multiplicative'` or
`'multiplicative replacement'` selects `mult`, and `'GBM'` or
`'geometric Bayesian-multiplicative'` selects `gbm`. Common other names
also work: `'multRepl'` (from zCompositions) and `'multi_replace'` (from
scikit-bio) for `mult`; `'cmultRepl'` (zCompositions' function, whose
default method is GBM) for `gbm`; `'square root'` for `sq`; and
`'Laplace'` for `bl`.

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

Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., &
Palarea-Albaladejo, J. (2015). Bayesian-multiplicative treatment of
count zeros in compositional data sets. *Statistical Modelling*, 15,
134–158.
[doi:10.1177/1471082X14535524](https://doi.org/10.1177/1471082X14535524)

Lubbe, S., Filzmoser, P., & Templ, M. (2021). Comparison of zero
replacement strategies for compositional data with large numbers of
zeros. *Chemometrics and Intelligent Laboratory Systems*, 210, 104248.
[doi:10.1016/j.chemolab.2021.104248](https://doi.org/10.1016/j.chemolab.2021.104248)

## See also

[`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md),
[`replace_zeros_gbm()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md),
[`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md);
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) and
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
(`'clr'`), which take the result.

## Examples

``` r
    x <- replace_zeros(ex_counts, 'gbm')
    x[1:3, 1:4]
#>        Streptococcus Bacteroides Corynebacterium Haemophilus
#> Saliva    0.46872444 0.005786721     0.001056911 0.520804933
#> Gums      0.89408789 0.004509901     0.001056911 0.098090349
#> Nose      0.02176063 0.001978239     0.492581602 0.001978239
    attr(x, 'n_replaced')
#> [1] 4

    # No zeros left, so no pseudocount is needed.
    aitchison(x)
#> # A dist: 4 samples (lower triangle)
#> # Aitchison Distance
#>       Saliva  Gums Nose
#> Gums    1.95
#> Nose   11.2  10.6
#> Stool   9.36  9.09 10.4
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
#> [[2]]
#> # A dist: 4 samples (lower triangle)
#> # Aitchison Distance
#>       Saliva  Gums Nose
#> Gums    1.9
#> Nose   11.1  11.1
#> Stool   9.35  9.28 10.5
#> # Attributes: method
#> 
#> [[3]]
#> # A dist: 4 samples (lower triangle)
#> # Aitchison Distance
#>       Saliva  Gums Nose
#> Gums    1.95
#> Nose   11.2  10.6
#> Stool   9.36  9.09 10.4
#> # Attributes: method
#> 
#> [[4]]
#> # A dist: 4 samples (lower triangle)
#> # Aitchison Distance
#>       Saliva  Gums Nose
#> Gums    1.95
#> Nose   11.2  10.6
#> Stool   9.36  9.09 10.4
#> # Attributes: method
#> 
#> [[5]]
#> # A dist: 4 samples (lower triangle)
#> # Aitchison Distance
#>       Saliva  Gums Nose
#> Gums    1.95
#> Nose   11.2  10.6
#> Stool   9.36  9.11 10.3
#> # Attributes: method
#> 
```
