# Plot source contributions

The usual stacked bar chart of a source-tracking result: one bar per
sink, divided into its estimated proportions from each source, with the
unknown source last, in grey.

## Usage

``` r
# S3 method for class 'ecodive_source'
plot(x, horiz = FALSE, legend = "right", ...)
```

## Arguments

- x:

  The result of
  [`sourcetracker()`](https://cmmr.github.io/ecodive/reference/sourcetracker.md).

- horiz:

  Draw the bars horizontally. Default: `FALSE`.

- legend:

  Where to put the legend of sources: a position keyword for
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html) (such
  as `"topright"`), or `NULL` for none. Default: `"right"`, outside the
  plot.

- ...:

  Passed on to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html), such
  as `main`, `col`, `names.arg`, `las` or `border`.

## Value

`x`, invisibly.

## Examples

``` r
    counts <- rbind(
      a1 = c(90, 10, 0, 0), a2 = c(80, 20, 0, 0),
      b1 = c(0, 0, 70, 30), b2 = c(0, 5, 60, 35),
      s1 = c(45, 10, 30, 15), s2 = c(80, 15, 5, 0) )
    st <- sourcetracker(
      counts, c('a', 'a', 'b', 'b', NA, NA), depth = NULL, n_burnin = 10 )
    plot(st)
```
