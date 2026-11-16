# documentation

documentation

## Ordination result

Every ordination method
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md), and
[`ordination()`](https://cmmr.github.io/ecodive/reference/ordination.md))
returns the same kind of object: a plain numeric matrix with one row per
sample (named after the samples) and one column per axis (named
`"PCoA1"`, `"NMDS1"`, etc.), so that it can go straight into
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) or a data
frame. Its attributes always begin with these three, in this order:

- `method`: The method's ID, e.g. `"pcoa"` or `"nmds"`.

- `eig`: Each axis's eigenvalue, or `NA` for methods without them.

- `prop`: The proportion of the total variance that each axis explains,
  or `NA` where that isn't defined.

Attributes particular to the method follow them.
