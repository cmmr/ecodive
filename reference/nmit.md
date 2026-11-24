# Nonparametric Microbial Interdependence Test (NMIT) Distances

Distances between subjects in how their features vary together over time
(Zhang et al. 2017), as q2-longitudinal's `nmit` action computes them:
for each subject, the correlation of every pair of features over that
subject's samples (its time points); then, for each pair of subjects, a
matrix norm of the difference of their correlation matrices. NMIT's test
is a PERMANOVA of these distances against a grouping of the subjects:
run
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
on the result.

## Usage

``` r
nmit(
  counts,
  subjects,
  data = NULL,
  method = "kendall",
  matrix_norm = "frobenius",
  margin = 1L,
  norm = "percent",
  pseudocount = NULL,
  clust = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A count table: samples \\\times\\ features, or any of the types listed
  under **Input Types**. Every sample of a subject enters its
  correlations, in any order.

- subjects:

  The subject (individual) of each sample: one value per sample, matched
  to the samples by name when both are named, else in the samples'
  order; or, with `data`, the name of one of its columns (q2's
  `individual_id_column`). No `NA`s.

- data:

  Optional sample metadata: a data frame with one row per sample,
  holding the `subjects` column. Rows are matched to the samples in
  `counts` by row name, and rows for other samples are ignored; without
  row names, they must be in the samples' order. Default: `NULL`.

- method:

  The correlation over time: `'kendall'` (Kendall's \\\tau_b\\, the
  default, as q2's), `'pearson'` or `'spearman'`.

- matrix_norm:

  The norm of the difference of two subjects' correlation matrices:
  `'frobenius'` (the default, as q2's) or `'nuclear'` (the sum of its
  singular values).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- norm:

  How to normalize `counts` first: `'percent'` (relative abundances, the
  default, as q2 expects), `'none'`, `'clr'`, or any other `norm` of
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md).

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

- clust:

  Whether and how to order the subjects by clustering, as for
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md):
  `NULL` (the default) or `FALSE` to keep their order, or one of
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)'s linkage
  methods (such as `"average"`), or `"pcoa"`. The `order` attribute then
  gives each subject's position before clustering.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A `dist` object over the **subjects** (class
`c("ecodive_dist", "dist")`), labelled by subject and in their order (a
factor's levels, or sorted), with `method` `"nmit"` and title "NMIT
Distance". It can be given as `dm` to
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) or
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

## Details

**Correlations** are computed as pandas' `DataFrame.corr()` (which q2
calls) computes them. A pair involving a feature that doesn't vary
within the subject has no correlation, and counts as 0 (q2's
`fillna(0)`). On the diagonal, pandas puts 1 for Kendall's \\\tau_b\\
whatever the feature, but for Pearson and Spearman a constant feature's
own correlation is also undefined, so 0: subjects differ there too, and
that adds to their distance.

**Subjects with fewer than 3 samples** have correlations of \\\pm 1\\ or
0 only (one sample: none at all), which say little about how their
features vary together. q2 computes them all the same; so does `nmit()`,
with a warning naming those subjects. Leave them out first if you can.

**Distances**: the Frobenius norm of the difference \\\Delta\\ is
\\\sqrt{2\sum\_{a\<b} \Delta\_{ab}^2 + \sum_a \Delta\_{aa}^2}\\,
computed from the stored lower triangles. The nuclear norm is the sum of
the absolute eigenvalues of \\\Delta\\ (symmetric, so these are its
singular values): an eigendecomposition of a \\D \times D\\ matrix for
each pair of subjects, \\O(S^2 D^3)\\ in all for \\S\\ subjects and
\\D\\ features, which for hundreds of features takes minutes.

**Memory**: every subject's correlations are kept, \\S D (D + 1) / 2\\
numbers (50 subjects and 500 features: 50 MB). `nmit()` warns above 500
MB; as q2's documentation advises, filter out rare features first.

**The test**: give
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
the distances and one group per subject, matched by subject name, for
instance from the sample metadata:


      d    <- nmit(counts, 'subject', meta)
      subj <- unique(meta[, c('subject', 'group')])
      rownames(subj) <- subj$subject
      permanova(d, 'group', subj)

Its p-values equal those of other PERMANOVAs of the same distances, such
as QIIME 2's `beta-group-significance`, given the same permutations.

**Differences from q2-longitudinal**: q2 labels each subject by the ID
of its first sample, in the order the subjects first appear; ecodive
labels them by subject. Distances are the same to rounding (about
\\10^{-15}\\).

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

## References

Zhang, Y., Han, S. W., Cox, L. M., & Li, H. (2017). A multivariate
distance-based analytic framework for microbial interdependence
association test in longitudinal study. *Genetic Epidemiology*, 41,
769–778. [doi:10.1002/gepi.22065](https://doi.org/10.1002/gepi.22065)

Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
longitudinal and paired-sample analyses of microbiome data. *mSystems*,
3, e00219-18.
[doi:10.1128/mSystems.00219-18](https://doi.org/10.1128/mSystems.00219-18)

## See also

Other longitudinal:
[`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md),
[`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md),
[`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)

## Examples

``` r
    # Six subjects, each sampled 5 times; two diets.
    meta <- data.frame(
      subject = rep(sprintf('m%d', 1:6), each = 5),
      diet    = rep(c('x', 'y'), each = 15) )
    counts <- matrix(rpois(30 * 6, 30), 30, dimnames = list(NULL, letters[1:6]))

    d <- nmit(counts, 'subject', meta)
    d
#> # A dist: 6 subjects (lower triangle)
#> # NMIT Distance
#>      m1   m2   m3   m4   m5
#> m2 3.62
#> m3 3.39 2.04
#> m4 2.94 3.05 3.67
#> m5 3.97 3.01 3.26 2.73
#> m6 3.76 3.21 3.16 3.41 3.05
#> # Attributes: method

    # NMIT's test: a PERMANOVA of the subjects.
    subj <- unique(meta)
    rownames(subj) <- subj$subject
    permanova(d, 'diet', subj, permutations = 99)
#> Permutational multivariate analysis of variance (PERMANOVA) of NMIT Distance
#> 
#> F = 1.546, p = 0.22, z = 1.207
#> Permutations: 99 (free)
#> 
#> Model R2 = 0.2788. Model tested as a whole:
#> 
#>          Df SumOfSqs     R2     F Pr(>F)
#> Model     1     7.36 0.2788 1.546   0.22
#> Residual  4    19.04 0.7212             
#> Total     5    26.40 1.0000             
```
