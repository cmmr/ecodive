# Enterotypes

Clusters samples into community types as Arumugam et al. (2011) did for
the human gut: partitioning around medoids (PAM) on Jensen-Shannon
distances, with the number of types chosen by the Calinski-Harabasz (CH)
index, and each sample's fit given by its silhouette.

## Usage

``` r
enterotype(
  x,
  n_types = 2:10,
  ch = c("medoids", "centroids"),
  counts = NULL,
  margin = 1L,
  cpus = n_cpus()
)
```

## Arguments

- x:

  A count table (or proportions), from which
  [`jensen()`](https://cmmr.github.io/ecodive/reference/jensen.md)
  distances are computed; or a `dist` object, such as one from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  used as it is, so that other distances can be tried. Samples are rows
  unless `margin = 2`.

- n_types:

  The numbers of types to try (each from 2 to the number of samples
  minus 1); a single value fixes it. Default `2:10`, cut to `2:(n - 1)`
  for fewer than 11 samples.

- ch:

  Which form of the Calinski-Harabasz index chooses the number of types:
  `"medoids"` (the default), as
  `clusterSim::index.G1(centrotypes = "medoids")` computes it in
  Arumugam et al.'s published script, or `"centroids"`, the pseudo-F of
  the distances. See Details.

- counts:

  Only with distances as `x`: an optional count table of the same
  samples (matched by name, or else by position), from which each type's
  mean composition (`profiles`) and `top_feature` are computed. Without
  it, a `dist` gives no profiles.

- margin:

  If your samples are in the table's rows, set to `1L` (the default). If
  they are in columns, set to `2L`.

- cpus:

  How many parallel processing threads to use. Results are the same for
  any number. Default: all available cores.

## Value

A data frame of class
`c("ecodive_enterotype", "ecodive_types", "data.frame")`, one row per
sample (see
[ecodive_types](https://cmmr.github.io/ecodive/reference/ecodive_types.md)),
with columns `sample`, `type` (a factor) and `silhouette`. Its
attributes:

- `types`: a data frame, one row per type, with `type`, `n_samples`,
  `top_feature` (the feature of largest mean proportion), `medoid` (the
  medoid sample's name, or its position) and `silhouette` (the mean of
  the type's samples).

- `profiles`: a types \\\times\\ features matrix of mean proportions, or
  `NULL` for distances without `counts`.

- `selection`: a data frame, one row per number of types tried, with
  `n_types`, `ch`, `silhouette` (the mean), `objective` (PAM's mean
  distance to the nearest medoid) and `selected`.

- `objective`: PAM's mean distance from each sample to its nearest
  medoid for the chosen number of types, after BUILD and after SWAP, as
  [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html) reports
  it.

- `ch`: how the CH index was computed, `"medoids"` or `"centroids"`.

- `method`: `"enterotype"`.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
types' profiles or the selection (see
[`plot.ecodive_types()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_types.md));
its title is "Enterotypes of" the distance's title.

## Details

**PAM** (Kaufman & Rousseeuw 1990, ch. 2) chooses \\k\\ samples as
medoids, and assigns every sample to its nearest medoid. BUILD picks the
medoids one at a time: first the sample with the smallest total distance
to all others, then each time the sample that most reduces the total
distance to the nearest medoid. SWAP then makes the best exchange of a
medoid with a non-medoid, as long as one lowers the total. Ties go to
the Types are numbered in the order of their medoids' positions in the
input. Where samples tie, PAM breaks the tie as
[`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html) does in
nearly every case tested: in BUILD, the highest sample index; in SWAP,
the first swap found; and a sample equally near two medoids joins the
lower-indexed one.

**The Calinski-Harabasz index** compares the spread between types with
the spread within them; the number of types with the largest index is
chosen (the smallest, on ties). Both forms are \$\$CH = \frac{(T -
W)/(k - 1)}{W/(n - k)}\$\$ for \\n\\ samples in \\k\\ types, with \\T\\
the total spread and \\W\\ the spread within types.

- `ch = "medoids"` follows
  `clusterSim::index.G1(x, cl, d, centrotypes = "medoids")`, which
  Arumugam et al.'s script calls with the Jensen-Shannon distances as
  `d` and the samples' relative abundances as `x`. The medoid \\m_c\\ of
  each type, and \\M\\ of all samples, is the member with the smallest
  total distance to the others. The spreads are then measured by squared
  Euclidean distances \\e\\ between the samples' proportions: \\T =
  \sum_i e(i, M)^2\\ and \\W = \sum_i e(i, m\_{c(i)})^2\\. For distances
  given without `counts`, \\e\\ is the distance itself; this equals
  `index.G1()` when the distances are Euclidean distances of its `x`.

- `ch = "centroids"` is the pseudo-F of the distances (the F statistic
  of
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)):
  \\T = \frac1n\sum\_{i\<j}d\_{ij}^2\\ and \\W = \sum_c
  \frac{1}{n_c}\sum\_{i\<j \in c} d\_{ij}^2\\, for types of \\n_c\\
  samples. For Euclidean distances, this is exactly Calinski &
  Harabasz's (1974) index. It uses only the distances.

**The silhouette** of sample \\i\\ (Rousseeuw 1987) is \\s_i = (b_i -
a_i) / \max(a_i, b_i)\\, where \\a_i\\ is its mean distance to the rest
of its type and \\b_i\\ its smallest mean distance to another type; 0
for a sample alone in its type. Values near 1 mean a sample sits well
within its type; values near 0 or below, that it lies between types.
Arumugam et al. checked their clusters with the mean silhouette, which
`selection` reports for every number of types.

**Distances.** The paper used the Jensen-Shannon *distance*, the square
root of the divergence:
[`jensen()`](https://cmmr.github.io/ecodive/reference/jensen.md), not
[`jsd()`](https://cmmr.github.io/ecodive/reference/jsd.md). The base of
the logarithm doesn't matter, since PAM, the CH index and silhouettes
are unchanged when every distance is scaled. Zeros need no pseudocount
(\\0 \log 0 = 0\\); some published scripts add a small one (such as
\\10^{-6}\\), which changes the distances, and so the results, slightly.
Arumugam et al. used genus abundances, so aggregate features to genus
first to follow them.

**A caution.** Costea et al. (2018), with the authors of the original
paper, describe enterotypes as a coarse stratification of what is a
continuous landscape of community composition, rather than discrete
types. Check the silhouettes, and consider ordinations (e.g.
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)) alongside
them.

## Differences from cluster and clusterSim

- Medoids, assignments, objectives and silhouettes match
  [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html) (with its
  default `pamonce = FALSE`) and
  [`cluster::silhouette()`](https://rdrr.io/pkg/cluster/man/silhouette.html)
  to rounding, and the CH index matches `clusterSim::index.G1()`, as
  described above.

- [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html) numbers
  its clusters by their first appearance in the data; `enterotype()`
  numbers types by their medoids' positions. The partitions are the
  same.

- Exact ties can still be broken differently:
  [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html) breaks a
  few ties in its first BUILD step otherwise, which can lead to another,
  equally good, set of medoids.

## Same names in other packages

None known.

## References

Arumugam, M., Raes, J., Pelletier, E., et al. (2011). Enterotypes of the
human gut microbiome. *Nature*, 473, 174-180.
[doi:10.1038/nature09944](https://doi.org/10.1038/nature09944)

Costea, P. I., Hildebrand, F., Arumugam, M., et al. (2018). Enterotypes
in the landscape of gut microbial community composition. *Nature
Microbiology*, 3, 8-16.
[doi:10.1038/s41564-017-0072-8](https://doi.org/10.1038/s41564-017-0072-8)

Kaufman, L., & Rousseeuw, P. J. (1990). *Finding Groups in Data: An
Introduction to Cluster Analysis*. Wiley.
[doi:10.1002/9780470316801](https://doi.org/10.1002/9780470316801)

Calinski, T., & Harabasz, J. (1974). A dendrite method for cluster
analysis. *Communications in Statistics*, 3, 1-27.
[doi:10.1080/03610927408827101](https://doi.org/10.1080/03610927408827101)

Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the
interpretation and validation of cluster analysis. *Journal of
Computational and Applied Mathematics*, 20, 53-65.
[doi:10.1016/0377-0427(87)90125-7](https://doi.org/10.1016/0377-0427%2887%2990125-7)

## See also

[ecodive_types](https://cmmr.github.io/ecodive/reference/ecodive_types.md),
[`plot.ecodive_types()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_types.md),
[`jensen()`](https://cmmr.github.io/ecodive/reference/jensen.md)

## Examples

``` r
    # Types of iris flowers, from the proportions of their four
    # measurements (a stand-in for a table of genus abundances).
    et <- enterotype(iris[, 1:4])
    et
#> # Enterotypes of Jensen-Shannon Distance
#> # 150 samples in 2 types, chosen by the Calinski-Harabasz index from 2 to 10
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
#> 5       6  895      0.346    0.0163    FALSE
#> 6       7  796      0.3      0.0155    FALSE
#> 7       8  805      0.263    0.0146    FALSE
#> 8       9  739      0.261    0.014     FALSE
#> 9      10  647      0.266    0.0135    FALSE
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
    table(et$type, iris$Species)
#>    
#>     setosa versicolor virginica
#>   1     50          0         0
#>   2      0         50        50

    # The types as groups in a test.
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

    # Another distance, with the abundances for the profiles.
    enterotype(bray(iris[, 1:4]), n_types = 3, counts = iris[, 1:4])
#> # Enterotypes of Bray-Curtis Dissimilarity
#> # 150 samples in 3 types
#> #
#> # Types:
#>   type n_samples top_feature  medoid silhouette
#> 1 1           50 Sepal.Length 8           0.759
#> 2 2           54 Sepal.Length 56          0.369
#> 3 3           46 Sepal.Length 113         0.45
#> #
#> # Samples: a data.frame: 150 × 3
#>    sample type  silhouette
#>    <chr>  <fct>      <dbl>
#> 1  1      1          0.829
#> 2  2      1          0.783
#> 3  3      1          0.795
#> 4  4      1          0.776
#> 5  5      1          0.829
#> 6  6      1          0.677
#> 7  7      1          0.8
#> 8  8      1          0.832
#> 9  9      1          0.712
#> 10 10     1          0.784
#> # … 140 more rows
#> # Use print(x, n = ...) to see more rows
#> # Attributes: types, profiles, selection, objective, ch, method
```
