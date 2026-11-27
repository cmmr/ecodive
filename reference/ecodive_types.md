# Community types

[`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md),
[`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md) and
[`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md) assign each
sample to one of a few community types, and return a data frame of class
`c("ecodive_<method>", "ecodive_types", "data.frame")`, with one row per
sample, in the input's order, and the sample names as row names.

## Columns

- `sample`: the sample's name (or its position, for unnamed samples).

- `type`: the type assigned to it, a factor whose levels are the types
  in order: `"1"`, `"2"`, ... for
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md)
  and [`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md); the
  sub-CSTs `"I-A"` to `"V"` for
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md).

- Columns of the method: for
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md),
  `silhouette`, the sample's silhouette width; for
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md),
  `cst` (the higher-order community state type), `score` (the similarity
  to the assigned centroid) and `matched`; for
  [`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md),
  `probability`, the posterior probability of the assigned type.

## Attributes

- `types`: a data frame, one row per type: `type`, `n_samples`,
  `top_feature` (the feature of largest mean proportion, Arumugam et
  al.'s "driver"; `NA` without abundances), and columns of the method,
  such as
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md)'s
  `medoid` and `silhouette`,
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)'s
  `cst` and reference scores, or
  [`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md)'s `weight`
  and `precision`.

- `profiles`: a types \\\times\\ features matrix of each type's mean
  proportions (`NULL` without abundances); for
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md),
  the centroids; for
  [`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md), the
  Dirichlet means.

- `selection`: a data frame, one row per number of types tried, with the
  criterion that chose among them, and `selected`, `TRUE` for the one
  chosen. Absent (`NULL`) for
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md),
  whose centroids fix the types.

- Attributes of the method, such as
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)'s
  `scores` (every sample's similarity to every centroid), or
  [`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md)'s
  `posterior` and `alpha`.

- `method`: the function that made the result, such as `"enterotype"`.

## Using the types

`res$type` is a plain factor, so it can be given as `groups` to any
test, such as
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
or to an ordination's
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), to colour the
samples by type.

## See also

[`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md),
[`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md),
[`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md),
[`plot.ecodive_types()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_types.md)

## Examples

``` r
    et <- enterotype(iris[, 1:4], n_types = 2:5)
    et
#> # Enterotypes of Jensen-Shannon Distance
#> # 150 samples in 2 types, chosen by the Calinski-Harabasz index from 2 to 5
#> #
#> # Types:
#>   type n_samples top_feature  medoid silhouette
#> 1 1           50 Sepal.Length 39          0.825
#> 2 2          100 Sepal.Length 78          0.821
#> #
#> # Selection:
#>   n_types   ch silhouette objective selected
#> 1       2 2000      0.822    0.0268     TRUE
#> 2       3 1630      0.565    0.0208    FALSE
#> 3       4 1310      0.549    0.0191    FALSE
#> 4       5  988      0.405    0.0175    FALSE
#> #
#> # Samples: a data.frame: 150 × 3
#>    sample type  silhouette
#>    <chr>  <fct>      <dbl>
#> 1  1      1          0.877
#> 2  2      1          0.85
#> 3  3      1          0.88
#> 4  4      1          0.855
#> 5  5      1          0.871
#> 6  6      1          0.808
#> 7  7      1          0.836
#> 8  8      1          0.877
#> 9  9      1          0.858
#> 10 10     1          0.795
#> # … 140 more rows
#> # Use print(x, n = ...) to see more rows
#> # Attributes: types, profiles, selection, objective, ch, method
    table(et$type)
#> 
#>   1   2 
#>  50 100 
    attr(et, 'types')
#>   type n_samples  top_feature medoid silhouette
#> 1    1        50 Sepal.Length     39  0.8245819
#> 2    2       100 Sepal.Length     78  0.8207907
    attr(et, 'selection')
#>   n_types        ch silhouette  objective selected
#> 1       2 2000.8520  0.8220544 0.02681756     TRUE
#> 2       3 1630.9986  0.5647328 0.02078605    FALSE
#> 3       4 1311.9767  0.5491302 0.01906313    FALSE
#> 4       5  988.3541  0.4046429 0.01748311    FALSE

    # The types as groups.
    permanova(jensen(iris[, 1:4]), et$type)
#> Permutational multivariate analysis of variance (PERMANOVA) of Jensen-Shannon Distance
#> 
#> F =  1724, p = 0.001, z =  1274
#> Permutations: 999 (free)
#> 
#> Model R2 = 0.921. Model tested as a whole:
#> 
#>           Df SumOfSqs    R2    F Pr(>F)    
#> Model      1   1.4481 0.921 1724  0.001 ***
#> Residual 148   0.1243 0.079                
#> Total    149   1.5724 1.000                
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
