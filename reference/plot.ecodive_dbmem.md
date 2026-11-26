# Plot Moran's eigenvector maps

Each MEM's Moran's I, or maps of the MEMs over the samples' coordinates.

## Usage

``` r
# S3 method for class 'ecodive_dbmem'
plot(x, coords = NULL, mems = NULL, ...)
```

## Arguments

- x:

  The result of
  [`dbmem()`](https://cmmr.github.io/ecodive/reference/dbmem.md).

- coords:

  `NULL` (the default) for a bar chart of each MEM's Moran's I
  (attribute `moran_i`). Or the samples' coordinates, in the order of
  `x`'s rows: a matrix or data frame whose first two numeric columns are
  drawn as x and y, to map MEMs `mems`. Each sample is a square, filled
  for positive values and open for negative ones, sized by its absolute
  value.

- mems:

  Which MEMs to map with `coords`: column numbers or names of `x`.
  Default: the first four (or fewer).

- ...:

  Passed on to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html) (for
  the bar chart) or
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)
  (for each map), such as `main`, `col` or `cex`.

## Value

`x`, invisibly.

## Details

The bar chart reads only attribute `moran_i`, so columns taken with `[`
can't be drawn as one. Maps of several MEMs share the device, in a grid
of panels (restored afterwards).

## See also

[`dbmem()`](https://cmmr.github.io/ecodive/reference/dbmem.md)

## Examples

``` r
    g   <- expand.grid(x = 1:10, y = 1:6)
    mem <- dbmem(g)
    plot(mem)

    plot(mem, g, mems = 1:6)
```
