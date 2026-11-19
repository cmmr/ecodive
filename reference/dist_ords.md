# Find and Browse Distance-based Ordination Methods

Programmatic access to the list of ordination methods that
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md)
runs, and their associated functions.

## Usage

``` r
list_dist_ords(
  val = c("data.frame", "list", "func", "id", "name"),
  nm = c(NA, "id", "name")
)

match_dist_ord(method)
```

## Arguments

- val:

  Sets the return value for this function call. See "Value" section
  below. Default: `"data.frame"`

- nm:

  What value to use for the names of the returned object. Default is
  `"id"` when `val` is `"list"` or `"func"`, otherwise the default is
  `NA` (no name).

- method:

  The name of an ordination method to search for. Supports partial
  matching. All non-alpha characters are ignored.

## Value

**`match_dist_ord()`**

A `list` with the following elements.

- `name` : Method name, e.g. `"Principal Coordinates Analysis"`

- `id` : Method ID - also the name of the function, e.g. `"pcoa"`

- `func` : The function for this method, e.g.
  [`ecodive::pcoa`](https://cmmr.github.io/ecodive/reference/pcoa.md)

- `params` : Formal args for `func`, e.g. `c("dm", "n_axes", "cpus")`

**`list_dist_ords()`**

The returned object's type and values are controlled with the `val` and
`nm` arguments.

- `val = "data.frame"` : The data.frame from which the below options are
  sourced.

- `val = "list"` : A list of objects as returned by `match_dist_ord()`
  (above).

- `val = "func"` : A list of functions.

- `val = "id"` : A character vector of method IDs.

- `val = "name"` : A character vector of method names.

If `nm` is set, then the names of the vector or list will be the method
ID (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names
will be applied to the
[`rownames()`](https://rdrr.io/r/base/colnames.html) property of the
`data.frame`.

## Examples

``` r

    # A data.frame of all available methods.
    list_dist_ords()
#>                                                       name       id
#> 1 Discriminant Canonical Analysis of Principal Coordinates    capda
#> 2            Constrained Analysis of Principal Coordinates capscale
#> 3                       Distance-based Redundancy Analysis    dbrda
#> 4                                Isometric Feature Mapping   isomap
#> 5                          Metric Multidimensional Scaling      mds
#> 6                      Non-metric Multidimensional Scaling     nmds
#> 7                           Principal Coordinates Analysis     pcoa
#> 8              t-distributed Stochastic Neighbor Embedding     tsne
#> 9            Uniform Manifold Approximation and Projection     umap

    # Find the method for 'MDS'.
    m <- match_dist_ord('MDS')
    m$id
#> [1] "mds"
```
