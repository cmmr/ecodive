# Redundancy Analysis (RDA)

Constrained ordination of samples from their feature abundances: the
axes of a (normalized) count table that sample metadata explain,
followed by principal components of what they leave unexplained. Gives
the same results as `vegan::rda()`, computed straight from sparse input
without forming any \\n \times n\\ or dense table, so it takes far less
time and memory for large datasets. With `norm = 'hellinger'` or
`norm = 'chord'` it is transformation-based RDA (tb-RDA); with
`Condition()` terms, partial RDA.

## Usage

``` r
rda(
  counts,
  groups,
  data = NULL,
  formula = NULL,
  margin = 1L,
  norm = "none",
  pseudocount = NULL,
  scale = FALSE,
  n_axes = 3L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of feature abundances (samples \\\times\\ features),
  or any of the table types listed under "Input Types". Alternatively,
  as in `vegan::rda()`, a formula with the table on its left-hand side
  and the constraints on its right, such as `counts ~ treatment`, in
  place of `counts` and `formula`; the metadata may then be given
  second, as in `rda(counts ~ treatment, meta)`.

- groups:

  The one variable that constrains the ordination: a vector or factor
  with one value per sample, matched to the samples in `counts` by name
  when both are named, else in the same order. Or, when `data` is given,
  the name of one of its columns (or a vector with one value per row of
  `data`). Categories (factor, character, or logical) or a numeric
  covariate. The same as `formula = ~ groups`. Give `groups` or
  `formula`, or neither for an unconstrained RDA, which is a PCA.

- data:

  Sample metadata: a data frame with one row per sample, holding the
  columns named in `groups` or `formula`. Rows are matched to the
  samples in `counts` by row name, and rows for other samples are
  ignored; without row names, they must be in the same order as the
  samples. Default: `NULL`.

- formula:

  Instead of `groups`, for models of several variables: a one-sided
  formula naming the variables in `data` that constrain the ordination,
  such as `~ treatment + ph`. Factors, numeric covariates, and
  interactions are all allowed. Wrap variables in `Condition()` to
  remove their effect first (partial RDA), as in `vegan::rda()`:
  `~ treatment + Condition(site)`. With only `Condition()` terms, the
  result is the PCA of the table after removing their effect. A
  left-hand side, if given, is ignored. Default: `NULL`.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- norm:

  Normalize the incoming counts before the PCA. Options are:

  - `'none'`: No transformation.

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'chord'`: Each sample scaled to unit length (Euclidean norm 1).

  - `'hellinger'`: Square root of relative abundance (the Hellinger
    transformation).

  - `'binary'`: Presence/absence (each count is either 0 or 1).

  - `'clr'`: Centered log ratio, giving an Aitchison PCA.

  - `'rclr'`: Robust centered log ratio. Zeros stay 0.

  Default: `'none'`.

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'`. Ignored for
  other normalization methods. See **Pseudocount** section.

- scale:

  Scale each feature to unit variance (after `norm`, before the
  analysis), as in `vegan::rda(scale = TRUE)` and
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s `scale`,
  so that every feature has equal weight. Features with constant values
  are left out. Default: `FALSE`

- n_axes:

  The number of axes to return. Default: `3`. Reduced, with a warning,
  to the number of axes with positive eigenvalues, which is at most
  \\n - 1\\ for \\n\\ samples.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A numeric matrix with one row per sample and one column per axis, as
described in "Ordination result" below. Constrained axes come first,
named `"RDA1"`, `"RDA2"`, etc. If the constraints give fewer than
`n_axes` of them, unconstrained axes follow, named `"PC1"`, `"PC2"`,
etc. The values are `vegan`'s site scores ("wa" scores, scaling 2): what
`vegan::scores(x, display = "sites")` returns and what
`phyloseq::plot_ordination()` plots. Its attributes are those of
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md), with the
same meanings, plus `loadings`:

- `method`: `"rda"`.

- `eig`: The variance of each axis.

- `prop`: The proportion of the variance that each axis explains, not
  counting the conditional variance:
  `eig / (trace - inertia[["conditional"]])`. `vegan`'s "Proportion
  Explained".

- `trace`: The total variance: the sum of all features' variances (after
  `norm` and `scale`), as
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s `trace`.

- `inertia`: The total variance split into three parts, which add up to
  `trace`: a named vector of

  - `conditional`: explained by the `Condition()` variables (`0` if
    there are none);

  - `constrained`: explained by the constraints, after the conditions
    (the sum of all constrained eigenvalues);

  - `unconstrained`: the rest.

- `constrained_eig`: Every non-zero constrained eigenvalue.

- `r_squared`: The proportion of the total variance that the constraints
  explain, `inertia[["constrained"]] / trace`, as
  `vegan::RsquareAdj(x)$r.squared`.

- `adj_r_squared`: `r_squared` adjusted for the number of constraint
  columns, as `vegan::RsquareAdj(x)$adj.r.squared`; see
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) for the
  formula, with and without conditions. `NA` when the model leaves no
  residual degrees of freedom.

- `biplot`: Biplot scores of the constraints: one row per column of the
  model matrix and one column per axis (`0` for unconstrained axes).
  `vegan::scores(x, display = "bp")`.

- `centroids`: The centroid of each factor level on each axis, as
  `vegan::scores(x, display = "cn")`, or `NULL` if no variable is a
  factor.

- `loadings`: The species scores: one row per feature and one column per
  axis, `vegan::scores(x, display = "species")` (scaling 2). Unlike
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s
  `loadings`, these are not unit-length: each axis's are scaled by the
  square root of its share of the variance, so that the sites and
  species form a correlation biplot, as in `vegan`.

## Details

**Method**

With \\Z\\ the normalized table (samples in rows), \\J\\ the centering
matrix, and \\D\\ the diagonal matrix of the features' standard
deviations (or the identity when `scale = FALSE`), RDA (Rao 1964;
Legendre & Legendre 2012, §11.1) analyzes the centered table \\Y =
JZD^{-1} / \sqrt{n - 1}\\. With \\H\\ the projection onto the (centered)
constraints, after removing any conditions (Borcard et al. 1992), the
constrained axes are the principal components of the fitted values
\\HY\\, and the unconstrained axes those of the residuals \\(I - H)Y\\.
Each axis's eigenvalue is a variance, so they add up to the total
variance of \\Y\\.

`rda()` shares its solver with
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md): since
\\YY^T\\ is Gower's centered matrix of the Euclidean distances between
the rows of \\Y\\, the constrained axes come from a small eigenproblem
with one row per constraint column, and the unconstrained axes, when
needed, from block Lanczos iteration. Both need only products of the
table with a few vectors at a time, computed in parallel from the
non-zero values alone, as in
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md), so a sparse
table is never made dense.

**Transformation-based RDA**

RDA of raw counts measures differences between samples by Euclidean
distances, which suit species abundance data poorly (Legendre &
Gallagher 2001). The Hellinger and chord transformations fix that:
`rda(counts, ..., norm = 'hellinger')` is an RDA that preserves
Hellinger distances between samples, and `norm = 'chord'` chord
distances. `norm = 'clr'` (with a `pseudocount`) preserves Aitchison
distances.

**Relation to dbrda(), pca() and permanova()**

- `rda(counts, ..., norm = 'hellinger')` has the same axes, `prop`,
  \\R^2\\ and biplot scores as
  `dbrda(euclidean(counts, norm = 'hellinger'), ...)`: db-RDA of
  Euclidean distances is RDA.
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)
  analyzes \\YY^T\\ without dividing by \\n - 1\\, so its `trace`,
  `eig`, `inertia`, and `constrained_eig` are `rda()`'s times \\n - 1\\,
  and its scores `rda()`'s times \\(n - 1)^{1/4}\\. When any distance is
  4 or more (see
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)'s
  "Scale"), as Euclidean distances between raw counts often are,
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) divides
  by \\n - 1\\ too, and all of these are the same as `rda()`'s.

- With no constraints, `rda(counts)` is a PCA: its `eig`, `prop`, and
  `trace` are those of
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) with the
  same `norm` and `scale`. Its scores point the same way as
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s, but are
  in `vegan`'s scaling 2: each axis's scores are
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s divided
  by \\\sqrt{(n - 1) \cdot eig}\\ and multiplied by \\\sqrt\[4\]{(n - 1)
  \cdot trace}\\.

- The permutation test of an RDA (`vegan::anova.cca()`) is PERMANOVA of
  the Euclidean distances between the normalized samples: e.g.
  `permanova(euclidean(counts, norm = 'hellinger'), 'treatment', meta)`.
  Its \\R^2\\ is `rda()`'s `r_squared`, and its pseudo-F and p-value
  those of `anova(vegan::rda(...))` for the same permutations.

[`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md),
[`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md)
and `principal_response_curves()`, which are planned for ecodive, will
build on `rda()`.

**Matching vegan**

`rda()` gives the same scores, eigenvalues, inertia, \\R^2\\, biplot
scores, centroids and species scores as `vegan::rda()` on the normalized
table, up to the signs of the axes. Columns of the model matrix that are
linear combinations of earlier ones are dropped, with a message, as in
`vegan`. With `scale = TRUE`, features with constant values are left out
(their loadings are 0), where `vegan` would stop with an error.

**Same names in other packages**

`vegan::rda()` shares its name with ecodive's `rda()`, and whichever
package was attached last masks the other's. Write `ecodive::rda()`, or
use
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md),
which always runs ecodive's, to be sure.

**Axis signs and reproducibility**

An eigenvector's sign is arbitrary, so different software can return
mirror images of the same axis. `rda()` flips each axis so that its
largest score in absolute value is positive, and flips the biplot
scores, centroids and loadings with it. Results are identical from run
to run and for any number of `cpus`.

## Ordination result

Every ordination function
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), and so
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md))
returns the same kind of object: a numeric matrix with one row per
sample and one column per axis, so that it can go straight into
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) or a data
frame. Rows are named after the samples, when the input has names.
Columns are named after the method: `"PCoA1"`, `"PC1"`, `"RPCA1"`,
`"dbRDA1"` or `"CAP1"` (then `"MDS1"` for unconstrained axes), `"MDS1"`,
`"NMDS1"`, `"Isomap1"`, `"TSNE1"`, `"UMAP1"`, and so on.

Its attributes always begin with these three, in this order, each with
one value per column:

- `method`: The method's ID, e.g. `"pcoa"` or `"nmds"`.

- `eig`: Each axis's eigenvalue, on the method's own scale:
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) and
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) give
  the eigenvalues of Gower's centered matrix (sums of squares);
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) give
  variances, which are
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)'s divided
  by \\n - 1\\ for \\n\\ samples; and
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) and
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  follow `vegan`, which also divides by \\n - 1\\ when any distance is 4
  or more. `NA` for
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), which
  have no eigenvalues.

- `prop`: The proportion of the total variance (or inertia) that each
  axis explains. Unlike `eig`, it is comparable between methods, so use
  it to label axes. `NA` where it isn't defined.

Attributes particular to the method follow them. The matrix has class
`"ecodive_ord"` (as well as `"matrix"`), only so that
[`print()`](https://rdrr.io/r/base/print.html) lists the other
attributes by name rather than in full. Subsetting with `[`, or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), gives a
plain matrix or data frame without them.

## Input Types

The `counts` parameter is designed to accept a simple numeric matrix,
but seamlessly supports objects from the following biological data
packages:

- `phyloseq`

- `rbiom`

- `SummarizedExperiment`

- `TreeSummarizedExperiment`

For large datasets, standard matrix operations may be slow. See
[`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
for details on using optimized formats (e.g. sparse matrices) and
parallel processing.

## Pseudocount

The `pseudocount` parameter is only relevant when `norm = 'clr'`.

Zeros are undefined in the centered log-ratio (CLR) transformation. If
`norm = 'clr'`, `pseudocount` is `NULL` (the default), and zeros are
detected, the function uses half the minimum non-zero value
(`min(x[x>0]) / 2`) and issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`), or
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
whose result needs no pseudocount.

**Why this matters:** The choice of pseudocount is not neutral; it acts
as a weighting factor that can significantly distort downstream results,
especially for sparse datasets. See Gloor et al. (2017) and Kaul et al.
(2017) for open-access discussions on the mathematical implications, or
Costea et al. (2014) for the impact on community clustering.

See [`aitchison`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## References

Rao, C. R. (1964). The use and interpretation of principal component
analysis in applied research. *Sankhya A*, 26(4), 329-358.

Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful
transformations for ordination of species data. *Oecologia*, 129(2),
271-280.
[doi:10.1007/s004420100716](https://doi.org/10.1007/s004420100716)

Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the
spatial component of ecological variation. *Ecology*, 73(3), 1045-1055.
[doi:10.2307/1940179](https://doi.org/10.2307/1940179)

Legendre, P., Oksanen, J., & ter Braak, C. J. F. (2011). Testing the
significance of canonical axes in redundancy analysis. *Methods in
Ecology and Evolution*, 2(3), 269-277.
[doi:10.1111/j.2041-210X.2010.00078.x](https://doi.org/10.1111/j.2041-210X.2010.00078.x)

Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.).
Elsevier.

## See also

[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)

## Examples

``` r
    # tb-RDA: RDA of Hellinger-transformed abundances.
    counts <- round(iris[, 1:4] * 10)
    ord    <- rda(counts, 'Species', iris, norm = 'hellinger')
    plot(ord, col = iris$Species)

    attr(ord, 'r_squared')
#> [1] 0.9480809
    attr(ord, 'loadings')
#>                    RDA1        RDA2         PC1
#> Sepal.Length  0.3150483  0.04442764 -0.07668916
#> Sepal.Width   0.6627525 -0.01599961  0.03400069
#> Petal.Length -0.7703689  0.03375146 -0.01565281
#> Petal.Width  -0.7519133 -0.03006730  0.19051616

    # Chord-transformed, in vegan's style.
    ord <- rda(counts ~ Species, iris, norm = 'chord', n_axes = 2)
    str(ord)
#>  'ecodive_ord' num [1:150, 1:2] 0.193 0.174 0.191 0.168 0.197 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "RDA1" "RDA2"
#>  - attr(*, "method")= chr "rda"
#>  - attr(*, "eig")= num [1:2] 0.042308 0.000239
#>  - attr(*, "prop")= num [1:2] 0.94435 0.00534
#>  - attr(*, "trace")= num 0.0448
#>  - attr(*, "inertia")= Named num [1:3] 0 0.04255 0.00225
#>   ..- attr(*, "names")= chr [1:3] "conditional" "constrained" "unconstrained"
#>  - attr(*, "constrained_eig")= num [1:2] 0.042308 0.000239
#>  - attr(*, "r_squared")= num 0.95
#>  - attr(*, "adj_r_squared")= num 0.949
#>  - attr(*, "biplot")= num [1:2, 1:2] -0.345 -0.641 0.939 -0.768
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:2] "Speciesversicolor" "Speciesvirginica"
#>   .. ..$ : chr [1:2] "RDA1" "RDA2"
#>  - attr(*, "centroids")= num [1:3, 1:2] 0.1829 -0.064 -0.1189 -0.0317 0.1742 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "Speciessetosa" "Speciesversicolor" "Speciesvirginica"
#>   .. ..$ : chr [1:2] "RDA1" "RDA2"
#>  - attr(*, "loadings")= num [1:4, 1:2] 0.2875 0.7705 -1.2005 -0.5678 0.0893 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
#>   .. ..$ : chr [1:2] "RDA1" "RDA2"

    # Its permutation test: PERMANOVA of the Hellinger distances.
    permanova(euclidean(counts, norm = 'hellinger'), 'Species', iris, permutations = 99)
#> Permutational multivariate analysis of variance (PERMANOVA)
#> 
#> F =  1342, p = 0.01, z =  1215
#> Permutations: 99 (free)
#> 
#> Model R2 = 0.9481. Model tested as a whole:
#> 
#>           Df SumOfSqs     R2    F Pr(>F)   
#> Model      2    3.054 0.9481 1342   0.01 **
#> Residual 147    0.167 0.0519               
#> Total    149    3.221 1.0000               
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

    # Without constraints, a PCA.
    head(rda(counts))
#>            PC1        PC2         PC3
#> [1,] -1.727729  0.8581736 -0.13211676
#> [2,] -1.747050 -0.4755764 -0.99609633
#> [3,] -1.859597 -0.3894579  0.08471927
#> [4,] -1.767134 -0.8552227  0.14936586
#> [5,] -1.756432  0.8779415  0.42633175
#> [6,] -1.468153  1.9918464  0.79832646
```
