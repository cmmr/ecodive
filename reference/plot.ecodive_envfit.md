# Plot fitted variables

Draws the variables that
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) fitted
onto an ordination: an arrow for each numeric variable (or feature),
pointing the way it increases, with a length proportional to
\\\sqrt{r^2}\\, and a label at the centroid of each factor level. With
`add = TRUE`, they go on top of an ordination plot, such as `plot(ord)`.

## Usage

``` r
# S3 method for class 'ecodive_envfit'
plot(x, add = FALSE, axes = c(1, 2), p_max = 1, col = "blue3", ...)
```

## Arguments

- x:

  The result of
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md).

- add:

  Whether to add to the current plot, such as `plot(ord)` of the
  ordination that `x` was fitted to (`TRUE`), or start a new one
  (`FALSE`, the default).

- axes:

  Which two of `x`'s axes to draw, by number among them. Default
  `c(1, 2)`.

- p_max:

  Draw only the variables with a `p_value` at most `p_max`. Default `1`:
  all of them (including those with no p-value).

- col:

  The colour of the arrows and labels.

- ...:

  With `add = FALSE`, passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `xlim` or `asp`.

## Value

`x`, invisibly.

## Details

The arrows are scaled to fit the plot: with `add = TRUE`, the longest
reaches 80% of the way from the origin to the nearest edge of the plot;
with `add = FALSE`, 80% of the way to the furthest centroid (or they
keep their own lengths, \\\sqrt{r^2}\\ at most 1, when there are no
centroids). Only their directions and relative lengths have meaning, as
in any biplot. The centroids are in the ordination's units, so they are
drawn where they are. Factor levels are labelled with the level, or with
the row name (the variable and the level pasted together) when two
variables share a level's name.

## Examples

``` r
    ord <- pcoa(bray(iris[, 1:4]))
    fit <- envfit(ord, iris, permutations = 99)

    plot(ord, 'Species', iris)
    plot(fit, add = TRUE)


    # On their own.
    plot(fit)
```
