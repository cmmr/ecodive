# Log Ratio of Summed Abundances

One value per sample: the log of the summed abundance of the `numerator`
features over that of the `denominator` features, as Qurro plots it.

## Usage

``` r
log_ratio(counts, numerator, denominator, margin = 1L, pseudocount = 0)
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

  Added to every selected value before summing, so \\r c\\ is added to a
  numerator of \\r\\ features. `0` (the default) adds nothing; `NULL`
  applies the rule of the other log-ratio functions (half the smallest
  non-zero value of the whole table, with a warning, when it has zeros).
  See **Zeros**.

## Value

A numeric vector with one value per sample, in the order of `counts` and
named after the samples (unnamed if they have no names), whatever the
`margin`. Larger values mean relatively more of the numerator.
Attributes:

- `method`: `'log_ratio'`.

- `numerator`, `denominator`: the selected features' names (their
  indices when the features have no names), in the order of `counts`.

- `pseudocount`: the pseudocount used (`0` for none).

The vector drops straight into metadata (`data$ratio <- log_ratio(...)`)
and into any test or fit that takes a per-sample variable, such as
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) or
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md).

## Details

For a sample \\x\\, numerator features \\N\\ and denominator features
\\D\\, the log ratio is \$\$\log\frac{\sum\_{i \in N} (x_i +
c)}{\sum\_{j \in D} (x_j + c)}\$\$ with \\c\\ the pseudocount (natural
logs). With one feature on each side, this is the log ratio of two
features, and equals the matching
[`balance()`](https://cmmr.github.io/ecodive/reference/balance.md)
divided by \\\sqrt{1/2}\\.

Ratios between features can be interpreted without knowing each sample's
total load: if the ratio of A to B goes up, it went up whatever happened
to the total (the "reference frames" of Morton et al. 2019). Qurro
(Fedarko et al. 2020) is a tool for choosing the numerator and
denominator interactively from feature rankings; this function computes
the same log ratio (`log(Num_Sum / Denom_Sum)` in Qurro's `qarcoal`) for
use in R.

**Log ratio or balance?** Sums are robust to zeros in single features
and easy to explain ("the log of the ratio of total Firmicutes to total
Bacteroidetes"), but a sum is dominated by its most abundant features.
[`balance()`](https://cmmr.github.io/ecodive/reference/balance.md) uses
geometric means instead: every feature counts equally, it is a
coordinate of the isometric log-ratio geometry the CoDa theory is built
on, and its scale makes balances of different sizes comparable; but a
zero in any one of its features needs replacing first.

## Zeros

Zeros in single features are harmless as long as each sum is positive,
so no pseudocount is needed in almost every case, and none is used by
default (unlike the other log-ratio functions, whose default is `NULL`).
A sample whose numerator or denominator sums to 0 has no defined ratio:
its value is `NA`, with one warning naming such samples. (Qurro drops
those samples; ecodive keeps their place, so the result lines up with
`counts`.) A positive `pseudocount` gives them a value.

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

Morton, J. T., Marotz, C., Washburne, A., et al. (2019). Establishing
microbial composition measurement standards with reference frames.
*Nature Communications*, 10, 2719.

Fedarko, M. W., Martino, C., Morton, J. T., et al. (2020). Visualizing
'omic feature rankings and log-ratios using Qurro. *NAR Genomics and
Bioinformatics*, 2, lqaa023.

## See also

[`balance()`](https://cmmr.github.io/ecodive/reference/balance.md) for
the ratio of geometric means;
[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md) for every
feature against one;
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
[`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) for
finding features that move together.

## Examples

``` r
    # Oral streptococci and Haemophilus against skin taxa
    r <- log_ratio(
      ex_counts,
      numerator   = c('Streptococcus', 'Haemophilus'),
      denominator = c('Propionibacterium', 'Staphylococcus', 'Corynebacterium') )
    r
#>    Saliva      Gums      Nose     Stool 
#>  5.834811  6.086775 -3.714588  0.000000 
#> attr(,"method")
#> [1] "log_ratio"
#> attr(,"numerator")
#> [1] "Streptococcus" "Haemophilus"  
#> attr(,"denominator")
#> [1] "Corynebacterium"   "Propionibacterium" "Staphylococcus"   
#> attr(,"pseudocount")
#> [1] 0

    # The same, by index; zero sums give NA
    log_ratio(ex_counts, c(1, 4), c(3, 5, 6))
#>    Saliva      Gums      Nose     Stool 
#>  5.834811  6.086775 -3.714588  0.000000 
#> attr(,"method")
#> [1] "log_ratio"
#> attr(,"numerator")
#> [1] "Streptococcus" "Haemophilus"  
#> attr(,"denominator")
#> [1] "Corynebacterium"   "Propionibacterium" "Staphylococcus"   
#> attr(,"pseudocount")
#> [1] 0
    suppressWarnings(log_ratio(ex_counts, 2, 6))
#>    Saliva      Gums      Nose     Stool 
#>        NA  1.386294 -4.770685  6.415097 
#> attr(,"method")
#> [1] "log_ratio"
#> attr(,"numerator")
#> [1] "Bacteroides"
#> attr(,"denominator")
#> [1] "Staphylococcus"
#> attr(,"pseudocount")
#> [1] 0
```
