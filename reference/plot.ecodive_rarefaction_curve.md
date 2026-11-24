# Plot rarefaction curves

Draws the curves from
[`rarefaction_curve()`](https://cmmr.github.io/ecodive/reference/rarefaction_curve.md):
one line per sample, of a metric against depth (or coverage), or of
coverage against depth.

## Usage

``` r
# S3 method for class 'ecodive_rarefaction_curve'
plot(
  x,
  metric = NULL,
  groups = NULL,
  data = NULL,
  type = "depth",
  depth = NULL,
  legend = "bottomright",
  col = NULL,
  ...
)
```

## Arguments

- x:

  The result of
  [`rarefaction_curve()`](https://cmmr.github.io/ecodive/reference/rarefaction_curve.md).

- metric:

  Which of `x`'s metrics to draw, by ID or name (matched as
  [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
  matches them). Default: the first in `x`.

- groups:

  Colours the samples by group, with a legend, and draws each group's
  median curve as a thicker line: the name of a column of `data`, or a
  vector with one value per sample (all the samples in the `totals`
  attribute, matched by name when both are named, as for
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)).

- data:

  A data frame of sample metadata, with one row per sample (matched by
  row name, else in order), when `groups` names one of its columns.

- type:

  What to draw: `"depth"` (the default), the metric against depth;
  `"coverage"`, the metric against expected coverage (coverage-based
  rarefaction curves); or `"completeness"`, expected coverage against
  depth (sample completeness curves, Chao et al. 2014, Fig. 4). A prefix
  is enough.

- depth:

  A depth to mark with a vertical line (for `type = "depth"` or
  `"completeness"`), such as one you plan to rarefy to. It is labelled
  with how many of the samples have at least that many observations.
  Default: `NULL` (none).

- legend:

  Where to put the legend of `groups`: a keyword for
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html), such
  as `"bottomright"` (the default), or `NULL` for none.

- col:

  The lines' colours: one colour for every sample, or with `groups`, one
  per group (in order, or named by group). Defaults to grey without
  groups, and to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)'s
  `"Dark 3"` palette with them.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `sub`, `xlab`, `ylab`, `xlim`, `ylim` or `log`. They
  replace the defaults made from `x`.

## Value

`x`, invisibly.

## Details

The title is "Rarefaction Curves of" and the metric's name (see
[`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md)),
or "Sample Completeness Curves"; the subtitle says how the values were
made: "Expected values" for a closed form, or, for instance, "Mean of 10
rarefactions". A group's median curve is the median of its samples'
values at each depth that at least two of them share (or every depth,
for a group of one); past a sample's library size, the sample has no
value, so the median is over fewer samples (as QIIME 2's group medians).

## References

Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H.,
Colwell, R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation
with Hill numbers: a framework for sampling and estimation in species
diversity studies. *Ecological Monographs*, 84, 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1)

## Examples

``` r
    rc <- rarefaction_curve(ex_counts, c('observed', 'shannon'), times = 5)

    plot(rc, depth = 1000)

    plot(rc, 'shannon', groups = c('A', 'A', 'B', 'C'))

    plot(rc, type = 'coverage')

    plot(rc, type = 'completeness', main = 'How complete is each sample?')
```
