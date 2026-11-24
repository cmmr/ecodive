# Plot a Neutral Community Model Fit

The occurrence frequency of each feature against its mean relative
abundance (log scale), with the neutral model's fitted curve and its
Wilson band, as Burns et al. (2016) and neufit draw it. Points are
coloured by partition: above the band, within it, or below it.

## Usage

``` r
# S3 method for class 'ecodive_neutral_model'
plot(
  x,
  type = c("neutral", "binomial"),
  group = NULL,
  legend = "bottomright",
  ...
)
```

## Arguments

- x:

  A
  [`neutral_model()`](https://cmmr.github.io/ecodive/reference/neutral_model.md)
  result.

- type:

  `'neutral'` (the default): the neutral model's curve and band.
  `'binomial'`: the binomial (random sampling) model's curve too, \\1 -
  (1 - p)^N\\, dashed.

- group:

  With `groups`, which groups to draw (names or positions), one panel
  each. Default: `NULL`, all of them.

- legend:

  Where the legend goes, as a
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html) keyword
  (`"topleft"`, ...), or `NULL` for none.

- ...:

  Passed on to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `pch`, `cex` or `xlim`.

## Value

`x`, invisibly.

## Details

The curve and the band are drawn from the result's attributes
(`n_individuals`, `m`, `n_samples` and `conf_level`) over a grid of
abundances from a tenth of the smallest to 1, as neufit does.

## See also

[`neutral_model()`](https://cmmr.github.io/ecodive/reference/neutral_model.md)

## Examples

``` r
    set.seed(1)
    p <- rgamma(60, 0.3); p <- p / sum(p)
    counts <- t(replicate(30, {
      x <- rgamma(60, 1000 * 0.05 * p)
      as.vector(rmultinom(1, 1000, x / sum(x)))
    }))

    res <- neutral_model(counts)
    plot(res)

    plot(res, 'binomial', legend = 'topleft')

```
