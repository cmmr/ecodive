# Bayesian-Multiplicative and Count Zero Replacement

Replaces the zeros of a count table with small proportions that depend
on each sample's depth, and scales each sample's non-zero proportions
down by the same factor, keeping their ratios. The methods of
zCompositions' `cmultRepl()` (Martín-Fernández et al. 2015): geometric
Bayesian-multiplicative (`gbm`, the usual choice), square-root
Bayesian-multiplicative (`sq`), Bayes-Laplace Bayesian-multiplicative
(`bl`), and count zero multiplicative (`czm`).

## Usage

``` r
replace_zeros_gbm(counts, margin = 1L, frac = 0.65, adjust = TRUE)

replace_zeros_sq(counts, margin = 1L, frac = 0.65, adjust = TRUE)

replace_zeros_bl(counts, margin = 1L, frac = 0.65, adjust = TRUE)

replace_zeros_czm(
  counts,
  margin = 1L,
  frac = 0.65,
  threshold = 0.5,
  adjust = TRUE
)
```

## Arguments

- counts:

  A numeric matrix of non-negative counts (samples \\\times\\ features),
  or any of the objects listed under **Input Types**. Each sample's
  total is its number of reads, so give counts, not proportions (other
  values give a warning).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- frac:

  For `czm`, the fraction of the upper limit (`threshold` reads) that a
  zero becomes; and, with `adjust`, the fraction of a feature's smallest
  observed proportion that an imputed value above it becomes. Between 0
  and 1 (exclusive). Default: `0.65`.

- adjust:

  If `TRUE` (the default), an imputed proportion above the smallest
  proportion observed for that feature in any sample is replaced by
  `frac` times that proportion, so that no zero is imputed above what
  was actually seen.

- threshold:

  For `czm`, the upper limit of a zero, in reads: a zero in a sample of
  \\n\\ reads is below `threshold / n`. Above 0 and at most 1. Default:
  `0.5`.

## Value

A dense numeric `matrix` of proportions with no zeros, in the
orientation of `counts` (samples \\\times\\ features for `margin = 1`,
features \\\times\\ samples for `margin = 2`, and for a single sample
given as a vector), with its names. Each sample sums to 1. Attributes:

- `method`: `'gbm'`, `'sq'`, `'bl'` or `'czm'`.

- `n_adjusted`: the number of imputed values that `adjust` lowered.

- `n_replaced`: the number of zeros replaced.

## Details

For sample \\i\\ with counts \\c\_{ij}\\, depth \\n_i = \sum_j
c\_{ij}\\, proportions \\x\_{ij} = c\_{ij} / n_i\\ and \\D\\ features,
each zero becomes a proportion \\r\_{ij}\\, and each non-zero proportion
is scaled so that the sample sums to 1: \$\$r\_{ij} = x\_{ij} \Bigl(1 -
\sum\_{k:\\ c\_{ik} = 0} r\_{ik}\Bigr) \quad \text{if } c\_{ij} \>
0.\$\$

**Bayesian-multiplicative methods** (`gbm`, `sq`, `bl`). A zero is the
posterior estimate of a multinomial probability under a Dirichlet prior
with expectation \\t\_{ij}\\ and strength \\s_i\\: \$\$r\_{ij} =
\frac{s_i\\ t\_{ij}}{n_i + s_i} \quad \text{if } c\_{ij} = 0.\$\$ The
prior expectation is the feature's pooled proportion in the other
samples, \$\$t\_{ij} = \frac{\sum\_{k \ne i} c\_{kj}}{\sum\_{k \ne i}
n_k},\$\$ so a feature that is abundant elsewhere gets a larger value.
The methods differ in the strength:

|  |  |
|----|----|
| Method | \\s_i\\ |
| `gbm` | \\1 / g_i\\, with \\g_i\\ the geometric mean of \\t\_{i1}, \ldots, t\_{iD}\\ |
| `sq` | \\\sqrt{n_i}\\ |
| `bl` | \\D\\ |

**Count zero multiplicative** (`czm`). Every zero of sample \\i\\
becomes \\r\_{ij} = frac \cdot threshold / n_i\\: by default 65% of half
a read. The other samples are used only by `adjust`.

**Adjustment.** With `adjust = TRUE`, an \\r\_{ij}\\ above \\m_j =
\min\_{k:\\ c\_{kj} \> 0} x\_{kj}\\, the smallest proportion of feature
\\j\\ observed in any sample, becomes \\frac \cdot m_j\\ (before the
non-zero parts are scaled). This happens mostly in shallow samples.

**Requirements.** No sample may be all zeros. The Bayesian methods need
at least two samples, and no feature may be zero in every sample (its
prior would be 0). `gbm` also needs every feature to be non-zero in at
least two samples: a feature seen in a single sample makes that sample's
\\t\_{ij}\\ zero, its geometric mean zero and its prior strength
infinite. Filter such features out first (e.g.
`counts[, colSums(counts > 0) >= 2]`), or use `sq` or `bl`, which handle
them. `czm` works on a single sample, and replaces the zeros of a
feature that is zero everywhere like any other.

The proportion of zeros is not checked: features that are zeros in most
samples, as are most features of a microbiome table, are imputed mostly
from their prior. zCompositions warns about (or drops) rows and columns
that are at least 80% zeros; ecodive never drops any.

**Compared with zCompositions.** The results equal
`cmultRepl(counts, method = ..., output = "prop", z.delete = FALSE)` of
zCompositions 1.6.2, to within about \\10^{-15}\\ (relative). Unlike
`cmultRepl()`, a table without zeros gives its proportions instead of an
error, `czm` accepts a single sample and features that are zero
everywhere, and `czm` refuses zeros whose imputed values would leave the
sample's non-zero parts at zero or below.

The result feeds straight into
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) and any
function with `norm = 'clr'`, with no pseudocount.

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

Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., &
Palarea-Albaladejo, J. (2015). Bayesian-multiplicative treatment of
count zeros in compositional data sets. *Statistical Modelling*, 15,
134–158.
[doi:10.1177/1471082X14535524](https://doi.org/10.1177/1471082X14535524)

Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015). zCompositions
— R package for multivariate imputation of left-censored data under a
compositional approach. *Chemometrics and Intelligent Laboratory
Systems*, 143, 85–96.
[doi:10.1016/j.chemolab.2015.02.019](https://doi.org/10.1016/j.chemolab.2015.02.019)

## See also

[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
[`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md),
[`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md)

## Examples

``` r
    # Geometric Bayesian-multiplicative replacement.
    x <- replace_zeros_gbm(ex_counts)
    x
#>        Streptococcus Bacteroides Corynebacterium Haemophilus Propionibacterium
#> Saliva   0.468724440 0.005786721     0.001056911 0.520804933      0.0028933607
#> Gums     0.894087889 0.004509901     0.001056911 0.098090349      0.0011274753
#> Nose     0.021760633 0.001978239     0.492581602 0.001978239      0.2482690406
#> Stool    0.001624823 0.992767072     0.001624823 0.001624823      0.0007336343
#>        Staphylococcus
#> Saliva   0.0007336343
#> Gums     0.0011274753
#> Nose     0.2334322453
#> Stool    0.0016248234
#> attr(,"method")
#> [1] "gbm"
#> attr(,"n_adjusted")
#> [1] 4
#> attr(,"n_replaced")
#> [1] 4

    # The ratios between non-zero parts are kept.
    nz <- ex_counts[1, ] > 0
    all.equal(x[1, nz] / sum(x[1, nz]), ex_counts[1, nz] / sum(ex_counts[1, nz]))
#> [1] TRUE

    # The other methods.
    replace_zeros_sq(ex_counts)
#>        Streptococcus Bacteroides Corynebacterium Haemophilus Propionibacterium
#> Saliva   0.468724440 0.005786721     0.001056911 0.520804933      0.0028933607
#> Gums     0.894087889 0.004509901     0.001056911 0.098090349      0.0011274753
#> Nose     0.021760633 0.001978239     0.492581602 0.001978239      0.2482690406
#> Stool    0.001624823 0.992767072     0.001624823 0.001624823      0.0007336343
#>        Staphylococcus
#> Saliva   0.0007336343
#> Gums     0.0011274753
#> Nose     0.2334322453
#> Stool    0.0016248234
#> attr(,"method")
#> [1] "sq"
#> attr(,"n_adjusted")
#> [1] 4
#> attr(,"n_replaced")
#> [1] 4
    replace_zeros_bl(ex_counts)
#>        Streptococcus Bacteroides Corynebacterium Haemophilus Propionibacterium
#> Saliva   0.468724440 0.005786721     0.001056911 0.520804933       0.002893361
#> Gums     0.894087889 0.004509901     0.001056911 0.098090349       0.001127475
#> Nose     0.021760633 0.001978239     0.492581602 0.001978239       0.248269041
#> Stool    0.001624243 0.992412730     0.001624243 0.001624243       0.001090296
#>        Staphylococcus
#> Saliva   0.0007336343
#> Gums     0.0011274753
#> Nose     0.2334322453
#> Stool    0.0016242434
#> attr(,"method")
#> [1] "bl"
#> attr(,"n_adjusted")
#> [1] 3
#> attr(,"n_replaced")
#> [1] 4
    replace_zeros_czm(ex_counts)
#>        Streptococcus Bacteroides Corynebacterium Haemophilus Propionibacterium
#> Saliva   0.468680529 0.005786179    0.0009420290 0.520756144      0.0028930897
#> Gums     0.894705546 0.004513017    0.0003668172 0.098158112      0.0011282542
#> Nose     0.021760633 0.001978239    0.4925816024 0.001978239      0.2482690406
#> Stool    0.001625157 0.992970917    0.0016251570 0.001625157      0.0005284553
#>        Staphylococcus
#> Saliva    0.000942029
#> Gums      0.001128254
#> Nose      0.233432245
#> Stool     0.001625157
#> attr(,"method")
#> [1] "czm"
#> attr(,"n_adjusted")
#> [1] 0
#> attr(,"n_replaced")
#> [1] 4

    # Without the adjustment.
    attr(replace_zeros_czm(ex_counts), 'n_adjusted')
#> [1] 0
    attr(replace_zeros_czm(ex_counts, adjust = FALSE), 'n_adjusted')
#> [1] 0
```
