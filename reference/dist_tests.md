# Find and Browse Distance Tests

Programmatic access to the list of available statistical tests of a
distance matrix, and their associated functions.

## Usage

``` r
list_dist_tests(
  val = c("data.frame", "list", "func", "id", "name", "groups"),
  nm = c(NA, "id", "name"),
  groups = NULL
)

match_dist_test(test, groups = NULL)
```

## Arguments

- val:

  Sets the return value for this function call. See "Value" section
  below. Default: `"data.frame"`

- nm:

  What value to use for the names of the returned object. Default is
  `"id"` when `val` is `"list"` or `"func"`, otherwise the default is
  `NA` (no name).

- groups:

  Filter by the kind of variable the test takes: `"categorical"` or
  `"numeric"`. Default: `NULL` (no filtering).

- test:

  The name of a test to search for. Supports partial matching. All
  non-alpha characters are ignored.

## Value

**`match_dist_test()`**

A `list` with the following elements.

- `name` : Test name, e.g. `"Multi-Response Permutation Procedure"`

- `id` : Test ID - also the name of the function, e.g. `"mrpp"`

- `groups` : The kinds of variable the test takes, e.g. `"categorical"`,
  or `"categorical numeric"`.

- `func` : The function for this test, e.g.
  [`ecodive::mrpp`](https://cmmr.github.io/ecodive/reference/mrpp.md)

- `params` : Formal args for `func`, e.g.
  `c("dm", "groups", "data", ...)`

**`list_dist_tests()`**

The returned object's type and values are controlled with the `val` and
`nm` arguments.

- `val = "data.frame"` : The data.frame from which the below options are
  sourced.

- `val = "list"` : A list of objects as returned by `match_dist_test()`
  (above).

- `val = "func"` : A list of functions.

- `val = "id"` : A character vector of test IDs.

- `val = "name"` : A character vector of test names.

- `val = "groups"` : A character vector of the kinds of variable each
  test takes.

If `nm` is set, then the names of the vector or list will be the test ID
(`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names
will be applied to the
[`rownames()`](https://rdrr.io/r/base/colnames.html) property of the
`data.frame`.

## Examples

``` r

    # A data.frame of all available tests.
    list_dist_tests()
#>                                              name        id              groups
#> 1                        Analysis of Similarities    anosim         categorical
#> 2                                     Mantel Test    mantel             numeric
#> 3            Multi-Response Permutation Procedure      mrpp         categorical
#> 4 Permutational Multivariate Analysis of Variance permanova categorical numeric
#> 5         Homogeneity of Multivariate Dispersions  permdisp         categorical

    # The tests of a numeric variable.
    list_dist_tests(val = 'id', groups = 'numeric')
#> [1] "mantel"    "permanova"

    # Find the test for 'adonis2'.
    m <- match_dist_test('adonis2')
    m$id
#> [1] "permanova"
```
