# Print ecodive's results

ecodive's tables, ordinations and distance matrices print a preview, in
the manner of the tibble package, rather than every value: a header
saying what the object is (a data frame, matrix or `dist`) and how big,
and its title (see
[ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md));
the first rows, with each column's type under its name; as many columns
as fit the console; and a footer counting the rows and columns left out
and naming the object's other attributes. The object itself is
unchanged: `as.data.frame(x)`, `as.matrix(x)` or `unclass(x)` print it
in full.

## Usage

``` r
# S3 method for class 'ecodive_types'
print(x, n = NULL, width = NULL, digits = 3L, ...)

# S3 method for class 'ecodive_ord'
print(x, n = NULL, width = NULL, digits = 3L, ...)

# S3 method for class 'ecodive_dist'
print(x, n = NULL, width = NULL, digits = 3L, ...)
```

## Arguments

- x:

  A result of class `ecodive_ord` (an ordination), `ecodive_dist` (a
  distance matrix), or one of ecodive's data frames: from
  [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md),
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc_global()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`ancombc_pairwise()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`ancombc_dunnett()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`ancombc_trend()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md),
  [`da_sensitivity()`](https://cmmr.github.io/ecodive/reference/da_sensitivity.md),
  [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md),
  [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md),
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md),
  [`assoc_pairs()`](https://cmmr.github.io/ecodive/reference/assoc_pairs.md)
  or
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md).
  A community-types result (see
  [ecodive_types](https://cmmr.github.io/ecodive/reference/ecodive_types.md))
  prints its types and how their number was chosen before its first
  samples.

- n:

  The number of rows to print (samples, for a distance matrix). `NULL`,
  the default, prints every row of up to 20, and the first 10 of more;
  `Inf` prints them all.

- width:

  The width of the output, in characters: columns that don't fit are
  left out and named in the footer. Default `getOption("width")`.

- digits:

  The significant digits of each number. Default `3`.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
    # An ordination: the first rows of its coordinates.
    ord <- pcoa(bray(iris[, 1:4]))
    ord
#> # A matrix: 150 samples × 3 axes
#> # Principal Coordinates Analysis of Bray-Curtis Dissimilarity
#> # % of variance: 89.0, 9.3, 2.9
#>    PCoA1   PCoA2    PCoA3
#>  1 0.169 -0.0259 -0.00505
#>  2 0.171  0.0235  0.0308
#>  3 0.187 -0.0032  0.0322
#>  4 0.18   0.0113  0.0372
#>  5 0.176 -0.03   -0.00613
#>  6 0.136 -0.0398 -0.0431
#>  7 0.186 -0.0182  0.0156
#>  8 0.167 -0.0165  0.00287
#>  9 0.194  0.0406  0.0413
#> 10 0.171  0.011   0.0286
#> # … 140 more samples
#> # Use print(x, n = ...) to see more samples
#> # Attributes: method, eig, prop, trace, correction, constant

    # A distance matrix: the top-left corner of its lower triangle.
    bray(iris[, 1:4])
#> # A dist: 150 samples (lower triangle)
#> # Bray-Curtis Dissimilarity
#>         1      2      3      4      5      6      7      8      9     10
#> 2  0.0355
#> 3  0.0408 0.0265
#> 4  0.051  0.0265 0.0213
#> 5  0.0098 0.0355 0.0408 0.051
#> 6  0.0556 0.0909 0.0962 0.0962 0.0556
#> 7  0.0352 0.0417 0.0262 0.0262 0.0352 0.0806
#> 8  0.0148 0.0306 0.0359 0.0359 0.0148 0.0605 0.0303
#> 9  0.0681 0.0326 0.0383 0.0273 0.0681 0.123  0.043  0.0632
#> 10 0.0404 0.0157 0.0316 0.0211 0.0404 0.0857 0.0466 0.0254 0.0486
#> 11 0.0286 0.064  0.0693 0.0693 0.0286 0.027  0.0634 0.0335 0.0964 0.0588
#> # … 139 more samples
#> # Use print(x, n = ...) to see more samples
#> # Attributes: method

    # The values as a plain matrix.
    head(unclass(ord))
#>          PCoA1       PCoA2        PCoA3
#> [1,] 0.1692812 -0.02591990 -0.005053611
#> [2,] 0.1706354  0.02346821  0.030848760
#> [3,] 0.1870090 -0.00320247  0.032221661
#> [4,] 0.1801593  0.01132449  0.037211653
#> [5,] 0.1763516 -0.02998591 -0.006128668
#> [6,] 0.1363264 -0.03979963 -0.043074846
```
