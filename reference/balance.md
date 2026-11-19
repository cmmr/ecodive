# Isometric Log-Ratio Balance

One value per sample: the balance between two groups of features, the
scaled log ratio of their geometric means. This is an isometric
log-ratio balance (Egozcue & Pawlowsky-Glahn 2005), not a balanced
design.

## Usage

``` r
balance(counts, numerator, denominator, margin = 1L, pseudocount = NULL)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Counts,
  proportions, and the output of a zero replacement are all accepted:
  each sample's scale cancels out of the ratio.

- numerator, denominator:

  The features on the top and the bottom of the ratio, each given as
  feature names, (1-based) feature indices, or a logical vector with one
  value per feature. Both are required, each must select at least one
  feature, and they must not share any. The two may use different forms
  (names for one, indices for the other).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pseudocount:

  Added to every value of `counts` before taking logs. `NULL` (the
  default) adds nothing when `counts` has no zeros, and otherwise uses
  half the smallest non-zero value of the whole table, with a warning,
  as [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) does.
  `0` adds nothing, and gives `NA` for samples with a zero in a selected
  feature. See **Pseudocount**.

## Value

A numeric vector with one value per sample, in the order of `counts` and
named after the samples (unnamed if they have no names), whatever the
`margin`. Larger values mean relatively more of the numerator.
Attributes:

- `method`: `'balance'`.

- `numerator`, `denominator`: the selected features' names (their
  indices when the features have no names), in the order of `counts`.

- `pseudocount`: the pseudocount used (`0` for none).

## Details

For a sample \\x\\, numerator features \\N\\ (\\r\\ of them) and
denominator features \\D\\ (\\s\\ of them), \$\$\mathrm{balance}(x) =
\sqrt{\frac{rs}{r + s}}\\\log\frac{g(x_N + c)}{g(x_D + c)}\$\$ where
\\g\\ is the geometric mean and \\c\\ the pseudocount (natural logs).
This is the inner product of \\\mathrm{clr}(x)\\ with the basis vector
that has \\\sqrt{s / (r (r + s))}\\ on the numerator features and
\\-\sqrt{r / (s (r + s))}\\ on the denominator ones: one row of
scikit-bio's `sbp_basis()`. So `balance()` equals, to the last bit, the
coordinate of [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md)
whose row of a sequential binary partition has `1` on `numerator` and
`-1` on `denominator` (with the same `pseudocount`), without building
the full basis. Features in neither group don't affect the value.

The factor \\\sqrt{rs/(r+s)}\\ puts balances of different sizes on the
same scale: each is a coordinate of an orthonormal basis, so their
variances can be compared directly. Swapping `numerator` and
`denominator` flips the sign. Like every log ratio, a balance is
unchanged when a sample is multiplied by a constant (with no
pseudocount). See
[`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
for when to prefer the ratio of sums.

Balances are what selbal (Rivera-Pinto et al. 2018) selects, gneiss and
PhILR build from trees, and coda4microbiome generalizes to
log-contrasts.

**Same names in other packages**: `compositions::balance()` computes the
same balances, with the same sign, from a formula (`~ (A*B)/C`);
whichever package was attached last masks the other's. Write
`ecodive::balance()` to be sure.

## Pseudocount

A zero in any numerator or denominator feature makes its geometric mean
0, so zeros need replacing first. With `pseudocount = NULL` (the
default), no pseudocount is used when `counts` has no zeros; otherwise
half the smallest non-zero value of the **whole table** is added to
every value, with a warning, exactly as
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) does, so a
balance agrees with the matching
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) coordinate.
To avoid a pseudocount, replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md):
its result has none. With `pseudocount = 0`, samples with a zero in a
selected feature get `NA`, with one warning naming them.

The choice of pseudocount is not neutral: it can distort the results,
especially for sparse tables. See
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

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

Egozcue, J. J., & Pawlowsky-Glahn, V. (2005). Groups of parts and their
balances in compositional data analysis. *Mathematical Geology*, 37,
795–828.

Rivera-Pinto, J., Egozcue, J. J., Pawlowsky-Glahn, V., et al. (2018).
Balances: a new perspective for microbiome analysis. *mSystems*, 3,
e00053-18.

## See also

[`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
for the ratio of sums;
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) for all the
balances of a sequential binary partition at once;
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md);
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
[`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) for
finding features that move together.

## Examples

``` r
    # Replace zeros first, then take the balance
    x <- replace_zeros(ex_counts, 'mult')
    b <- balance(
      x,
      numerator   = c('Streptococcus', 'Haemophilus'),
      denominator = c('Propionibacterium', 'Staphylococcus', 'Corynebacterium') )
    b
#>    Saliva      Gums      Nose     Stool 
#>  3.938659  4.922648 -4.207824 -1.046615 
#> attr(,"method")
#> [1] "balance"
#> attr(,"numerator")
#> [1] "Streptococcus" "Haemophilus"  
#> attr(,"denominator")
#> [1] "Corynebacterium"   "Propionibacterium" "Staphylococcus"   
#> attr(,"pseudocount")
#> [1] 0

    # The first coordinate of ilr() with a partition that starts with it
    sbp <- rbind(
      c(1,  0, -1,  1, -1, -1),
      c(1,  0,  0, -1,  0,  0),
      c(0,  0,  1,  0, -1, -1),
      c(0,  0,  0,  0,  1, -1),
      c(1, -1,  1,  1,  1,  1) )
    all.equal(unname(b), unname(ilr(x, basis = sbp)[, 1]))
#> [1] "Attributes: < Modes: list, NULL >"                   
#> [2] "Attributes: < Lengths: 4, 0 >"                       
#> [3] "Attributes: < names for target but not for current >"
#> [4] "Attributes: < current is not list-like >"            
```
