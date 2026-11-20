# Find and Browse Association Measures

Programmatic access to the list of association measures that
[`feature_assoc()`](https://cmmr.github.io/ecodive/reference/feature_assoc.md)
runs, and their associated functions.

## Usage

``` r
list_feature_assocs(
  val = c("data.frame", "list", "func", "id", "name"),
  nm = c(NA, "id", "name")
)

match_feature_assoc(method)
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

  The name of an association measure to search for. Supports partial
  matching. All non-alpha characters are ignored.

## Value

**`match_feature_assoc()`**

A `list` with the following elements.

- `name` : Measure name, e.g. `"SparCC Correlation"`

- `id` : Measure ID - also the name of the function, e.g. `"sparcc"`

- `kind` : `"similarity"` (larger means more associated) or
  `"dissimilarity"`.

- `lower`, `upper` : The range of its values.

- `fisher_z` : Whether its values are Pearson-type correlations, which
  [`assoc_network_diff()`](https://cmmr.github.io/ecodive/reference/assoc_network_diff.md)
  can compare with Fisher's \\z\\.

- `sparse` : Whether its zeros mean "no edge" (a selected graph).

- `strict` : Whether its own cutoff keeps \\\|r\| \>\\ the cutoff,
  rather than \\\ge\\.

- `func` : The function for this measure, e.g.
  [`ecodive::sparcc`](https://cmmr.github.io/ecodive/reference/sparcc.md)

- `params` : Formal args for `func`, e.g. `c("counts", "margin", ...)`

**`list_feature_assocs()`**

The returned object's type and values are controlled with the `val` and
`nm` arguments.

- `val = "data.frame"` : The data.frame from which the below options are
  sourced.

- `val = "list"` : A list of objects as returned by
  `match_feature_assoc()` (above).

- `val = "func"` : A list of functions.

- `val = "id"` : A character vector of measure IDs.

- `val = "name"` : A character vector of measure names.

If `nm` is set, then the names of the vector or list will be the measure
ID (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names
will be applied to the
[`rownames()`](https://rdrr.io/r/base/colnames.html) property of the
`data.frame`.

## Examples

``` r

    # A data.frame of all available measures.
    list_feature_assocs()
#>                      name     id          kind lower upper fisher_z sparse
#> 1  Variance of Log Ratios    vlr dissimilarity     0   Inf    FALSE  FALSE
#> 2 Proportionality (rho_p)  rho_p    similarity    -1     1    FALSE  FALSE
#> 3 Proportionality (phi_s)  phi_s dissimilarity     0   Inf    FALSE  FALSE
#> 4      SparCC Correlation sparcc    similarity    -1     1     TRUE  FALSE
#>   strict
#> 1  FALSE
#> 2  FALSE
#> 3  FALSE
#> 4  FALSE

    # Find the measure for 'propr'.
    m <- match_feature_assoc('propr')
    m$id
#> [1] "rho_p"
```
