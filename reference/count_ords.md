# Find and Browse Count Table Ordination Methods

Programmatic access to the list of ordination methods that
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
runs, and their associated functions.

## Usage

``` r
list_count_ords(
  val = c("data.frame", "list", "func", "id", "name"),
  nm = c(NA, "id", "name")
)

match_count_ord(method)
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

**`match_count_ord()`**

A `list` with the following elements.

- `name` : Method name, e.g. `"Principal Component Analysis"`

- `id` : Method ID - also the name of the function, e.g. `"pca"`

- `func` : The function for this method, e.g.
  [`ecodive::pca`](https://cmmr.github.io/ecodive/reference/pca.md)

- `params` : Formal args for `func`, e.g. `c("counts", "margin", ...)`

**`list_count_ords()`**

The returned object's type and values are controlled with the `val` and
`nm` arguments.

- `val = "data.frame"` : The data.frame from which the below options are
  sourced.

- `val = "list"` : A list of objects as returned by `match_count_ord()`
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
    list_count_ords()
#>                                                     name     id
#> 1                                Correspondence Analysis     ca
#> 2                      Canonical Correspondence Analysis    cca
#> 3                           Linear Discriminant Analysis    lda
#> 4 Orthogonal Partial Least Squares Discriminant Analysis oplsda
#> 5                           Principal Component Analysis    pca
#> 6            Partial Least Squares Discriminant Analysis  plsda
#> 7                                    Redundancy Analysis    rda
#> 8                                   Robust Aitchison PCA   rpca

    # Find the method for 'prcomp'.
    m <- match_count_ord('prcomp')
    m$id
#> [1] "pca"
```
