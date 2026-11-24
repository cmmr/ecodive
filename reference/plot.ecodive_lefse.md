# Plot a LEfSe result

LEfSe's bar chart of LDA scores: one horizontal bar per reported feature
(`reject`), its length the feature's LDA score, coloured by the class it
is enriched in.

## Usage

``` r
# S3 method for class 'ecodive_lefse'
plot(x, n_features = Inf, ...)
```

## Arguments

- x:

  The result of
  [`lefse()`](https://cmmr.github.io/ecodive/reference/lefse.md).

- n_features:

  How many reported features to show, those with the largest scores.
  Default `Inf`: all of them.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `sub` or `xlim`.

## Value

`x`, invisibly.

## Details

As LEfSe draws it: with two classes, the first class's bars point left
(negative) and the second's right, each sorted by score; with more, the
bars are grouped by class, in the classes' order, and sorted by score
within each. The dashed lines are at `lda_threshold`. The features that
passed both tests but scored at or below the threshold are not shown.

## Examples

``` r
    set.seed(1)
    x <- matrix(rpois(40 * 12, 50), 40, 12, dimnames = list(NULL, paste0('F', 1:12)))
    x[21:40, 1:4] <- x[21:40, 1:4] * 4
    x[1:20, 5:6]  <- x[1:20, 5:6] * 3
    plot(lefse(x, rep(c('a', 'b'), each = 20), cpus = 1))
```
