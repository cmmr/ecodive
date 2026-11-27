# Plot interpolated and extrapolated Hill numbers

Draws the curves of
[`hill_extrapolation()`](https://cmmr.github.io/ecodive/reference/hill_extrapolation.md):
one per sample, rarefaction as a solid line, extrapolation dashed, each
sample's own size as a point, with interval bands when `x` has them.

## Usage

``` r
# S3 method for class 'ecodive_hill_extrapolation'
plot(
  x,
  type = "depth",
  q = NULL,
  samples = NULL,
  legend = "bottomright",
  col = NULL,
  ...
)
```

## Arguments

- x:

  The result of
  [`hill_extrapolation()`](https://cmmr.github.io/ecodive/reference/hill_extrapolation.md).

- type:

  What to draw: `"depth"` (the default; `"size"` too), the Hill number
  against depth (iNEXT's type 1); `"completeness"`, sample coverage
  against depth (type 2); or `"coverage"`, the Hill number against
  coverage (type 3). A prefix is enough.

- q:

  Which order to draw, one of `x`'s. Default: the first.

- samples:

  Which samples to draw, by name or position. Default: all.

- legend:

  Where to put the legend of samples (drawn for up to 10): a keyword for
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html), or
  `NULL` for none. Default: `"bottomright"`.

- col:

  The samples' colours, one per sample drawn. Default:
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)'s
  `"Dark 3"` palette.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `ylab`, `xlim`, `ylim` or `log`.

## Value

`x`, invisibly.

## Details

The titles are "Sample-size-based Rarefaction and Extrapolation",
"Sample Completeness Curves" and "Coverage-based Rarefaction and
Extrapolation"; the y axis is "Species richness" (\\q = 0\\), "Shannon
diversity" (\\q = 1\\) or "Simpson diversity" (\\q = 2\\), as in Chao et
al. (2014). The subtitle names the bands' confidence level.

## References

Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H.,
Colwell, R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation
with Hill numbers: a framework for sampling and estimation in species
diversity studies. *Ecological Monographs*, 84, 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1)

## Examples

``` r
    he <- hill_extrapolation(ex_counts, n_boot = 10)
    plot(he)

    plot(he, q = 2, type = 'coverage')

    plot(he, type = 'completeness', samples = c('Saliva', 'Stool'))
```
