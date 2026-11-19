# Plot a permutation test

Plots the result of one of ecodive's permutation tests. By default,
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
[`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md),
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md) and
[`mantel_partial()`](https://cmmr.github.io/ecodive/reference/mantel_partial.md)
results show the permutation distribution of their statistic, with the
observed value marked;
[`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) shows each
group's mean within-group distance;
[`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
each group's dispersion; and
[`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md)
the Procrustes fit. `type = "null"` shows the permutation distribution
for any of them. Nothing is recomputed: the plot is made from the stored
results.

## Usage

``` r
# S3 method for class 'ecodive_test'
plot(x, type = "null", ...)

# S3 method for class 'ecodive_mrpp'
plot(x, type = c("groups", "null"), ...)

# S3 method for class 'ecodive_permdisp'
plot(x, type = c("groups", "null"), ...)

# S3 method for class 'ecodive_procrustes_test'
plot(x, type = c("fit", "null"), ...)
```

## Arguments

- x:

  A permutation test result, of class `ecodive_test`: from
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
  [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md),
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md),
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md),
  [`mantel_partial()`](https://cmmr.github.io/ecodive/reference/mantel_partial.md),
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md),
  [`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md),
  or
  [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md).

- type:

  What to draw: `"null"` for the permutation distribution of the
  statistic; or for
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) and
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md),
  `"groups"` (their default) for a bar for each group; or for
  [`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md),
  `"fit"` (its default) for the Procrustes rotation, drawn by
  [`plot.ecodive_ord()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_ord.md).

- ...:

  Passed on to the plotting function:
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html) for
  `"null"`,
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html) for
  `"groups"`, and
  [`plot.ecodive_ord()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_ord.md)
  for `"fit"`. Such as `main`, `sub`, `xlab`, `col` or `breaks`; they
  replace the defaults made from `x`.

## Value

`x`, invisibly.

## Details

The title is the test's name, and the subtitle its statistic, p-value
and standardized effect size \\z\\, as printed.

- `"null"`: a histogram of the statistic in each permutation (the
  `f_perm`, `r_perm` or `delta_perm` attribute), with a line at the
  observed statistic. The p-value is the proportion of permutations,
  with the observed one, at least as extreme as the observed statistic.
  It needs a test run with `permutations` above 0.

- `"groups"` for
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md): the mean
  within-group distance of each group (`class_delta`), with lines at the
  observed \\\delta\\ (their weighted mean) and the \\\delta\\ expected
  by chance (`expected_delta`, the mean of all distances).

- `"groups"` for
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md):
  each group's mean distance to its centroid or median (`dispersions`).

- `"fit"` for
  [`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md):
  the symmetric Procrustes rotation it tested (its `procrustes` field),
  with an arrow from each sample to the target ordination.

## Examples

``` r
    dm <- bray(iris[, 1:4])

    # The pseudo-F of each permutation, and the observed F.
    res <- permanova(dm, 'Species', iris, permutations = 199)
    plot(res)


    # MRPP: each species' mean within-group distance.
    res <- mrpp(dm, 'Species', iris, permutations = 199)
    plot(res)

    plot(res, type = 'null')


    # PERMDISP: the dispersion of each species.
    plot(permdisp(dm, 'Species', iris, permutations = 199))
#> Warning: Some squared distances to group centres are negative, and were set to 0.


    # PROTEST: the Procrustes fit of two ordinations.
    res <- procrustes_test(pcoa(dm, n_axes = 2), pca(iris[, 1:4], n_axes = 2))
    plot(res)
```
