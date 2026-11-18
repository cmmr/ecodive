# Find and Browse Distance-based Ordination Methods

Programmatic access to the list of ordination methods that
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md)
runs, and their associated functions.

## Usage

``` r
list_dist_ords(
  val = c("data.frame", "list", "func", "id", "name", "input"),
  nm = c(NA, "id", "name"),
  input = NULL
)

match_dist_ord(method, input = NULL)
```

## Arguments

- val:

  Sets the return value for this function call. See "Value" section
  below. Default: `"data.frame"`

- nm:

  What value to use for the names of the returned object. Default is
  `"id"` when `val` is `"list"` or `"func"`, otherwise the default is
  `NA` (no name).

- input:

  Filter by the kind of data the method takes: `"dist"` for distances.
  Default: `NULL` (no filtering).

- method:

  The name of an ordination method to search for. Supports partial
  matching. All non-alpha characters are ignored.

## Value

**`match_dist_ord()`**

A `list` with the following elements.

- `name` : Method name, e.g. `"Principal Coordinates Analysis"`

- `id` : Method ID - also the name of the function, e.g. `"pcoa"`

- `input` : The kind of data the method takes, e.g. `"dist"`.

- `func` : The function for this method, e.g.
  [`ecodive::pcoa`](https://cmmr.github.io/ecodive/reference/pcoa.md)

- `params` : Formal args for `func`, e.g. `c("dm", "k", "cpus")`

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

- `val = "input"` : A character vector of the kinds of data each method
  takes.

If `nm` is set, then the names of the vector or list will be the method
ID (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names
will be applied to the
[`rownames()`](https://rdrr.io/r/base/colnames.html) property of the
`data.frame`.

## Examples

``` r

    # A data.frame of all available methods.
    list_dist_ords()
#>                                            name       id input
#> 1 Constrained Analysis of Principal Coordinates capscale  dist
#> 2            Distance-based Redundancy Analysis    dbrda  dist
#> 3                     Isometric Feature Mapping   isomap  dist
#> 4               Metric Multidimensional Scaling      mds  dist
#> 5           Non-metric Multidimensional Scaling     nmds  dist
#> 6                Principal Coordinates Analysis     pcoa  dist
#> 7   t-distributed Stochastic Neighbor Embedding     tsne  dist
#> 8 Uniform Manifold Approximation and Projection     umap  dist

    # Find the method for 'MDS'.
    m <- match_dist_ord('MDS')
    m$id
#> [1] "mds"
```
