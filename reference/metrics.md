# Find and Browse Available Metrics

Programmatic access to the lists of available metrics, and their
associated functions.

## Usage

``` r
list_metrics(
  div = c(NA, "alpha", "beta"),
  val = c("data.frame", "list", "func", "id", "name", "div", "phylo", "weighted",
    "int_only", "true_metric"),
  nm = c(NA, "id", "name"),
  phylo = NULL,
  weighted = NULL,
  int_only = NULL,
  true_metric = NULL
)

match_metric(
  metric,
  div = NULL,
  phylo = NULL,
  weighted = NULL,
  int_only = NULL,
  true_metric = NULL
)
```

## Arguments

- div:

  Filter by diversity type. One of `"alpha"` or `"beta"`. Default: `NA`
  (no filtering).

- val:

  Sets the return value for this function call. See "Value" section
  below. Default: `"data.frame"`

- nm:

  What value to use for the names of the returned object. Default is
  `"id"` when `val` is `"list"` or `"func"`, otherwise the default is
  `NA` (no name).

- phylo:

  Filter by whether a phylogenetic tree is required. `TRUE` returns only
  phylogenetic metrics. `FALSE` returns only non-phylogenetic metrics.
  Default: `NULL` (no filtering).

- weighted:

  Filter by whether relative abundance is used. `TRUE` returns
  quantitative metrics. `FALSE` returns qualitative (presence/absence)
  metrics. Default: `NULL` (no filtering).

- int_only:

  Filter by whether integer counts are required. `TRUE` returns metrics
  requiring integers (e.g. richness estimators). `FALSE` returns metrics
  that accept proportions. Default: `NULL` (no filtering).

- true_metric:

  Filter by whether the metric satisfies the triangle inequality. `TRUE`
  returns proper distance metrics. `FALSE` returns dissimilarities.
  Default: `NULL` (no filtering).

- metric:

  The name of an alpha/beta diversity metric to search for. Supports
  partial matching. All non-alpha characters are ignored.

## Value

**`match_metric()`**

A `list` with the following elements.

- `name` : Metric name, e.g. `"Faith's Phylogenetic Diversity"`

- `id` : Metric ID - also the name of the function, e.g. `"faith"`

- `div` : Either `"alpha"` or `"beta"`.

- `phylo` : `TRUE` if metric requires a phylogenetic tree; `FALSE`
  otherwise.

- `weighted` : `TRUE` if metric takes relative abundance into account;
  `FALSE` if it only uses presence/absence.

- `int_only` : `TRUE` if metric requires integer counts; `FALSE`
  otherwise.

- `true_metric` : `TRUE` if metric is a true metric and satisfies the
  triangle inequality; `FALSE` if it is a non-metric dissimilarity; `NA`
  for alpha diversity metrics.

- `func` : The function for this metric, e.g.
  [`ecodive::faith`](https://cmmr.github.io/ecodive/reference/faith.md)

- `params` : Formal args for `func`, e.g.
  `c("counts", "norm", "tree", "cpus")`

**`list_metrics()`**

The returned object's type and values are controlled with the `val` and
`nm` arguments.

- `val = "data.frame"` : The data.frame from which the below options are
  sourced.

- `val = "list"` : A list of objects as returned by `match_metric()`
  (above).

- `val = "func"` : A list of functions.

- `val = "id"` : A character vector of metric IDs.

- `val = "name"` : A character vector of metric names.

- `val = "div"` : A character vector `"alpha"` and/or `"beta"`.

- `val = "phylo"` : A logical vector indicating which metrics require a
  tree.

- `val = "weighted"` : A logical vector indicating which metrics take
  relative abundance into account (as opposed to just presence/absence).

- `val = "int_only"` : A logical vector indicating which metrics require
  integer counts.

- `val = "true_metric"` : A logical vector indicating which metrics are
  true metrics and satisfy the triangle inequality, which work better
  for ordinations such as PCoA.

If `nm` is set, then the names of the vector or list will be the metric
ID (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names
will be applied to the
[`rownames()`](https://rdrr.io/r/base/colnames.html) property of the
`data.table`.

## Examples

``` r
    # A data.frame of all available metrics.
    head(list_metrics())
#>                                       name            id phylo weighted
#> 1 Abundance-based Coverage Estimator (ACE)           ace FALSE     TRUE
#> 2                       Aitchison Distance     aitchison FALSE     TRUE
#> 3                      Berger-Parker Index        berger FALSE     TRUE
#> 4                   Bhattacharyya Distance bhattacharyya FALSE     TRUE
#> 5                Bray-Curtis Dissimilarity          bray FALSE     TRUE
#> 6                          Brillouin Index     brillouin FALSE     TRUE
#>   int_only true_metric   div
#> 1     TRUE          NA alpha
#> 2    FALSE        TRUE  beta
#> 3    FALSE          NA alpha
#> 4    FALSE        TRUE  beta
#> 5    FALSE       FALSE  beta
#> 6     TRUE          NA alpha
    
    # All alpha diversity function names.
    list_metrics('alpha', val = 'id')
#>  [1] "ace"         "berger"      "brillouin"   "chao1"       "faith"      
#>  [6] "fisher"      "simpson"     "inv_simpson" "margalef"    "mcintosh"   
#> [11] "menhinick"   "observed"    "shannon"     "squares"    
    
    # Try to find a metric named 'otus'.
    m <- match_metric('otus')
    
    # The result is a list that includes the function.
    str(m)
#> List of 9
#>  $ name       : chr "Observed Features"
#>  $ id         : chr "observed"
#>  $ phylo      : logi FALSE
#>  $ weighted   : logi FALSE
#>  $ int_only   : logi FALSE
#>  $ true_metric: logi NA
#>  $ div        : chr "alpha"
#>  $ func       :function (counts, margin = 1L, cpus = n_cpus())  
#>  $ params     : chr [1:3] "counts" "margin" "cpus"
```
