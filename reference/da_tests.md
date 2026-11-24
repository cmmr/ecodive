# Find and Browse Differential Abundance Methods

Programmatic access to the list of differential-abundance methods that
[`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md) runs,
and their associated functions.

## Usage

``` r
list_da_tests(
  val = c("data.frame", "list", "func", "id", "name"),
  nm = c(NA, "id", "name")
)

match_da_test(method)
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

  The name of a DA method to search for. Supports partial matching. All
  non-alpha characters are ignored.

## Value

**`match_da_test()`**

A `list` with the following elements.

- `name` : Method name, e.g. `"Welch t-test"`

- `id` : Method ID - also the name of the function, e.g. `"da_ttest"`

- `design` : `"groups"` (one grouping, given as `groups`) or `"formula"`
  (a model, given as `formula` and/or `groups`)

- `effect` : The result's column of effect sizes, e.g. `"log2_fc"`, or
  `NA` (ANCOM gives none)

- `base` : The base of the effect sizes' logs: `"e"`, `"2"`, or `"norm"`
  (natural logs for a log-ratio `norm`, no logs otherwise); `NA` without
  effect sizes

- `zero_free_ok` : Whether the method may be given a table with its
  zeros replaced
  ([`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)):
  `FALSE` for the methods that model counts

- `func` : The function for this method, e.g.
  [`ecodive::da_ttest`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)

- `params` : Formal args for `func`, e.g. `c("counts", "groups", ...)`

**`list_da_tests()`**

The returned object's type and values are controlled with the `val` and
`nm` arguments.

- `val = "data.frame"` : The data.frame from which the below options are
  sourced.

- `val = "list"` : A list of objects as returned by `match_da_test()`
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
    list_da_tests()
#>                                     name            id  design      effect base
#> 1 Analysis of Composition of Microbiomes         ancom  groups        <NA> <NA>
#> 2                               ANCOM-BC       ancombc formula      log_fc    e
#> 3                              ANCOM-BC2      ancombc2 formula      log_fc    e
#> 4           Dirichlet-Multinomial t-test dirmult_ttest  groups     log2_fc    2
#> 5                 Wilcoxon Rank-Sum Test   da_wilcoxon  groups median_diff norm
#> 6                           Welch t-test      da_ttest  groups        diff norm
#> 7                           Linear Model         da_lm formula    estimate norm
#> 8                                  LinDA         linda formula     log2_fc    2
#> 9                                  LEfSe         lefse  groups   lda_score <NA>
#>   zero_free_ok
#> 1         TRUE
#> 2         TRUE
#> 3         TRUE
#> 4        FALSE
#> 5         TRUE
#> 6         TRUE
#> 7         TRUE
#> 8         TRUE
#> 9         TRUE

    # Find the method for 'Mann-Whitney'.
    m <- match_da_test('Mann-Whitney')
    m$id
#> [1] "da_wilcoxon"
```
