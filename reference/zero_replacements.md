# Find and Browse Zero Replacement Methods

Programmatic access to the list of zero replacement methods that
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
runs, and their associated functions.

## Usage

``` r
list_zero_replacements(
  val = c("data.frame", "list", "func", "id", "name"),
  nm = c(NA, "id", "name")
)

match_zero_replacement(method)
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

  The name of a zero replacement method to search for. Supports partial
  matching. All non-alpha characters are ignored.

## Value

**`match_zero_replacement()`**

A `list` with the following elements.

- `name` : Method name, e.g. `"Multiplicative Replacement"`

- `id` : Method ID, e.g. `"mult"`

- `func` : The function for this method, e.g.
  [`ecodive::replace_zeros_mult`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md)

- `params` : Formal args for `func`, e.g.
  `c("counts", "margin", "delta")`

**`list_zero_replacements()`**

The returned object's type and values are controlled with the `val` and
`nm` arguments.

- `val = "data.frame"` : The data.frame from which the below options are
  sourced.

- `val = "list"` : A list of objects as returned by
  `match_zero_replacement()` (above).

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
    list_zero_replacements()
#>                         name   id
#> 1 Multiplicative Replacement mult

    # Find the method for zCompositions' 'multRepl'.
    m <- match_zero_replacement('multRepl')
    m$id
#> [1] "mult"
```
