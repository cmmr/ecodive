# Plot a SIMPER result

A dot chart of the taxa that contribute most to the Bray-Curtis
dissimilarity between each pair of groups, from
[`simper()`](https://cmmr.github.io/ecodive/reference/simper.md): each
taxon's average contribution, largest first, in one panel per pair.
Filled points mark the taxa whose adjusted permutation p-value is at
most `alpha`.

## Usage

``` r
# S3 method for class 'ecodive_simper'
plot(x, n_taxa = 10, pairs = NULL, alpha = 0.05, ...)
```

## Arguments

- x:

  The result of
  [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md).

- n_taxa:

  How many taxa to show for each pair, those with the largest average
  contributions. Default `10`; `Inf` for all.

- pairs:

  Which pairs of groups to show, by name (such as `"setosa-virginica"`).
  Default: all of them.

- alpha:

  Taxa with a `p_adjusted` at most `alpha` are drawn with filled points.
  Default `0.05`.

- ...:

  Passed on to
  [`graphics::dotchart()`](https://rdrr.io/r/graphics/dotchart.html),
  such as `main`, `xlab`, `cex` or `col` (its `color`); and `sub`, to
  [`graphics::title()`](https://rdrr.io/r/graphics/title.html).

## Value

`x`, invisibly.

## Details

Each pair's label gives its mean between-group Bray-Curtis dissimilarity
(the `overall` attribute), which the contributions of all its taxa add
up to. Without permutations (`permutations = 0`), all points are open.

## Examples

``` r
    res <- simper(round(iris[, 1:4] * 10), 'Species', iris, permutations = 99)
    plot(res)

    plot(res, pairs = 'setosa-virginica')
```
