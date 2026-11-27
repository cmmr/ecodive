# Plot zero patterns

Draws a
[`zero_pattern()`](https://cmmr.github.io/ecodive/reference/zero_pattern.md)
result: the most frequent patterns as rows, the features that are
sometimes zero as columns, zeros dark; with each feature's percentage of
zeros above, and each pattern's percentage of the samples on the right.

## Usage

``` r
# S3 method for class 'ecodive_zero_pattern'
plot(x, n_patterns = 30, n_features = 50, ...)
```

## Arguments

- x:

  A result of
  [`zero_pattern()`](https://cmmr.github.io/ecodive/reference/zero_pattern.md).

- n_patterns:

  The number of patterns to draw, the most frequent first. Default `30`.
  `Inf` or `TRUE` for all.

- n_features:

  The largest number of features to draw. Default `50`. `Inf` or `TRUE`
  for all.

- ...:

  Passed on to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html); `main`
  and `sub` replace the title and subtitle.

## Value

`x`, invisibly.

## Details

**Features drawn.** Only features that are zero in some samples but not
all: those never zero and those always zero say nothing about patterns,
and are counted in the subtitle instead. When more than `n_features`
remain, those whose proportion of zeros is closest to one half are
drawn, as the most informative. Columns are sorted by their proportion
of zeros, the fewest on the left.

**Patterns drawn.** The `n_patterns` most frequent, top to bottom.
Patterns are defined on all the features, so two rows can look the same
when they differ only in features that aren't drawn; the subtitle says
how many are hidden.

When nearly every sample has its own pattern, as in most microbiome
tables, the top bars (each feature's zeros) are the useful part.

## See also

[`zero_pattern()`](https://cmmr.github.io/ecodive/reference/zero_pattern.md),
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
[`structural_zeros()`](https://cmmr.github.io/ecodive/reference/structural_zeros.md).

## Examples

``` r
    plot(zero_pattern(ex_counts), n_features = 20)

```
