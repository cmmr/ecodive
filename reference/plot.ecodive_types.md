# Plot community types

Draws a result of
[`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md)
with base graphics.

## Usage

``` r
# S3 method for class 'ecodive_types'
plot(x, type = c("profiles", "selection"), n_features = 10, ...)
```

## Arguments

- x:

  A result of
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md).

- type:

  `"profiles"` or `"selection"`.

- n_features:

  The number of features to show in `"profiles"`. Default `10`.

- ...:

  Passed on to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html) or
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `ylab` or `col`.

## Value

`x`, invisibly.

## Details

- `"profiles"` (the default): stacked bars of each type's mean
  composition, with the `n_features` features that reach the largest
  mean proportion in any type, and the rest pooled as "Other". The bars
  are labelled with the type and its number of samples.

- `"selection"`: the criterion that chose the number of types (for
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md),
  the Calinski-Harabasz index) against the number of types, with the one
  chosen circled.

The types keep the order and colours of `levels(x$type)` in every plot.

## Examples

``` r
    et <- enterotype(iris[, 1:4], n_types = 2:5)
    plot(et)

    plot(et, 'selection')
```
