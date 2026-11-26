# Plot a multiple regression on distance matrices

Draws each term's standardized coefficient (`std_estimate`), which can
be compared between terms in different units: filled for the terms whose
p-value is at most `alpha`, open for the others, with a dashed line at
0. The subtitle gives the model's \\R^2\\ and p-value.

## Usage

``` r
# S3 method for class 'ecodive_mantel_regression'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  A
  [`mantel_regression()`](https://cmmr.github.io/ecodive/reference/mantel_regression.md)
  result.

- alpha:

  The significance level of the filled points, compared with each term's
  `p_value`. Default: `0.05`.

- ...:

  Passed on to
  [`graphics::dotchart()`](https://rdrr.io/r/graphics/dotchart.html),
  such as `main`, `xlab`, `col`, `xlim`, or `sub`.

## Value

`x`, invisibly.

## Examples

``` r
    set.seed(1)
    meta <- data.frame(ph = runif(20, 6, 8), temp = runif(20, 10, 20))
    counts <- matrix(rpois(20 * 10, 5), 20, 10)
    counts[, 1:3] <- counts[, 1:3] + round(4 * (meta$ph - 6))

    res <- mantel_regression(bray(counts), data = meta, formula = ~ ph + temp)
    plot(res)

```
