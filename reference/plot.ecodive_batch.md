# Plot a batch correction

MMUPHin's diagnostic figure for a corrected table, in base graphics: how
far the empirical Bayes step shrank each batch effect, and each
feature's batch means before and after the correction.

## Usage

``` r
# S3 method for class 'ecodive_batch'
plot(x, type = c("shrinkage", "means"), legend = "topleft", ...)
```

## Arguments

- x:

  The result of
  [`batch_mmuphin()`](https://cmmr.github.io/ecodive/reference/batch_mmuphin.md)
  or
  [`batch_combat()`](https://cmmr.github.io/ecodive/reference/batch_combat.md).

- type:

  `"shrinkage"` (the default): each feature's location effect in each
  batch, \\\hat\gamma\\ (x), against its shrunken value \\\gamma^\*\\
  (y), coloured by batch, with the 1:1 line. Points below the line in
  the upper right (and above it in the lower left) were pulled towards
  the batch's average effect. `"means"`: for each feature (one vertical
  line of points per feature), its mean in each batch (y) against its
  mean over all samples (x), before the correction (black) and after
  (red, just to the right), with the 1:1 line. A batch effect spreads a
  feature's batch means away from the line; the correction draws them
  together.
  [`batch_mmuphin()`](https://cmmr.github.io/ecodive/reference/batch_mmuphin.md)'s
  means are of relative abundances;
  [`batch_combat()`](https://cmmr.github.io/ecodive/reference/batch_combat.md)'s
  are of the values it corrected.

- legend:

  Where to put the legend: a position keyword for
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html) (such
  as `"topleft"`), or `NULL` for none. Default: `"topleft"`.

- ...:

  Passed on to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `ylab`, `log`, `pch`, `cex`, `col`
  (`"shrinkage"`: one colour per batch; `"means"`: before and after).

## Value

`x`, invisibly.

## Details

The plot reads only the attributes stored with the result, so a table
taken apart with `[` can't be plotted. The real check of a correction is
an ordination coloured by batch, and the PERMANOVA \\R^2\\ of batch,
before and after (see the examples of
[`batch_mmuphin()`](https://cmmr.github.io/ecodive/reference/batch_mmuphin.md)).

## See also

[`batch_mmuphin()`](https://cmmr.github.io/ecodive/reference/batch_mmuphin.md),
[`batch_combat()`](https://cmmr.github.io/ecodive/reference/batch_combat.md)

## Examples

``` r
    set.seed(1)
    counts <- matrix(rpois(20 * 30, 20), 20, 30)
    batch  <- rep(c('run1', 'run2'), each = 10)
    counts[11:20, 1:10] <- counts[11:20, 1:10] * 3

    x <- batch_mmuphin(counts, batch)
    plot(x)

    plot(x, type = 'means')
```
