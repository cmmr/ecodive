# Similarity Percentages (SIMPER)

Breaks down the mean Bray-Curtis dissimilarity between each pair of
groups into the contribution of each taxon, to rank the taxa that differ
most between the groups (Clarke 1993). Permutation tests give each
taxon's contribution a p-value. Gives the same results as
`vegan::simper()` and [`summary()`](https://rdrr.io/r/base/summary.html)
of it, but runs in parallel C on sparse tables. Its p-values count ties
with the observed contribution, as Clarke (1993) does, so they are never
smaller than vegan's, which are too small for small groups (see
**Permutations**).

## Usage

``` r
simper(
  counts,
  groups,
  data = NULL,
  margin = 1L,
  norm = "none",
  p_adjust = "holm",
  permutations = 999L,
  strata = NULL,
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A feature table, such as taxa abundances. Accepts the same objects as
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md):
  a numeric matrix (samples in rows, or in columns with `margin = 2`), a
  sparse matrix from Matrix or slam, or a `phyloseq`, `rbiom`, or
  `SummarizedExperiment` object. Sparse tables are never made dense.
  Values must be non-negative, and every sample must have a non-zero
  total.

- groups:

  Each sample's group: a vector or factor with one value per sample,
  matched to the samples in `counts` by name when both are named, else
  in the same order. Or, when `data` is given, the name of one of its
  columns (or a vector with one value per row of `data`). Every distinct
  value is a group, numbers included; there must be at least two. A
  factor's unused levels are dropped.

- data:

  Sample metadata: a data frame with one row per sample, holding the
  columns named in `groups` or `strata`. Rows are matched to the samples
  in `counts` by row name, and rows for other samples are ignored;
  without row names, they must be in the same order as the samples.
  Default: `NULL`.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- norm:

  Normalize `counts` first: `'none'` (as `vegan::simper()`, which uses
  the values as given), `'percent'`, `'chord'`, `'hellinger'`, or
  `'binary'`, as in
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).
  `'clr'` and `'rclr'` are not allowed: they give negative values, which
  Bray-Curtis can't use. Default: `'none'`.

- p_adjust:

  How to adjust each pair's p-values for testing many taxa: one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html)
  (`'holm'`, `'hochberg'`, `'hommel'`, `'bonferroni'`, `'BH'`, `'BY'`,
  `'fdr'`, or `'none'`), matched as
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) matches
  them. The taxa of each pair are adjusted together, leaving out those
  with a `p_value` of `NA`. Default: `'holm'`.

- permutations:

  The number of permutations. `0` skips the tests, and gives `NA`
  p-values. Default: `999`. Alternatively, a matrix with one permutation
  of `1:n` per row, as from `permute::shuffleSet()`, to use instead of
  random ones: it is cut down to each pair's samples (see
  **Permutations**).

- strata:

  Permute samples only within these groups (blocks), for designs such as
  repeated measures. With `data`, the name of one of its columns, or a
  vector with one value per row of `data`. Without `data`, a vector with
  one value per sample, matched to the samples by name when it has
  names, otherwise in the same order. Default: `NULL` (free
  permutation).

- seed:

  An integer seed for the random permutations. The same seed gives the
  same p-values. Default: `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A data frame with one row per pair of groups and taxon, the pairs in the
order of the groups' factor levels (sorted values, for other vectors),
and the taxa of each pair in decreasing order of `average`. Its columns
are:

- `pair`: The pair of groups, as `"a-b"`.

- `taxon`: The taxon's name (its column number, as text, when `counts`
  has no taxa names).

- `average`: The taxon's mean contribution to the Bray-Curtis
  dissimilarity between the two groups' samples. A pair's `average`s sum
  to its mean between-group Bray-Curtis dissimilarity.

- `sd`: The standard deviation of the taxon's contribution over the
  pairs of samples (`NA` when the groups have one sample each).

- `ratio`: `average / sd`. Taxa with a high ratio contribute
  consistently, not just a lot in a few pairs of samples.

- `mean_a`, `mean_b`: The taxon's mean abundance (after `norm`) in the
  first and second group of the pair.

- `cumsum`: The cumulative proportion of the pair's mean dissimilarity
  that this taxon and those above it account for.

- `p_value`: The permutation p-value: the proportion of permutations
  giving the taxon a larger `average` (see **Permutations**). `NA` with
  `permutations = 0`, and for taxa whose `average` is 0 (those absent
  from both groups, or with the same abundance in every sample).

- `p_adjusted`: `p_value` adjusted for multiple testing with `p_adjust`,
  over the taxa of the pair.

The data frame has class `"ecodive_simper"`, only to mark what made it.
Attributes: `overall`, each pair's mean between-group Bray-Curtis
dissimilarity (named by pair); `n_perm`, the number of permutations;
`p_adjust`, the adjustment method; and `permutation`, `"free"`,
`"within strata"` or `"supplied matrix"`.

## Details

**Method**

For samples \\i\\ and \\j\\ with totals \\T_i\\ and \\T_j\\, the
Bray-Curtis dissimilarity \\\sum_k \|x\_{ik} - x\_{jk}\| / (T_i + T_j)\\
is a sum over the taxa, so taxon \\k\\ contributes \\c\_{ijk} =
\|x\_{ik} - x\_{jk}\| / (T_i + T_j)\\. For two groups \\a\\ and \\b\\, a
taxon's `average` is the mean of \\c\_{ijk}\\ over the \\n_a n_b\\ pairs
of a sample from \\a\\ and one from \\b\\, and `sd` their standard
deviation (with denominator \\n_a n_b - 1\\). SIMPER is defined for
Bray-Curtis only.

**Permutations**

A taxon's p-value (Clarke 1993) comes from shuffling the group labels of
the pair's samples (only within `strata`, when given), and recomputing
its `average`. The p-value is \\(1 + \\\\average\_{perm} \ge average\\)
/ (1 + permutations)\\, Clarke's (1993) "at least as large" rule, with
permuted averages within a relative \\10^{-7}\\ of the observed one
counted as ties (for rounding). Permutations that only swap samples
within groups reproduce the observed `average` exactly, so they count.
`vegan::simper()` counts only permuted averages that are strictly
larger, which makes its p-values too small when groups are small: with
two samples per group, a third of all permutations reproduce the
observed `average`, and a taxon whose observed `average` is the largest
possible gets \\p \approx 1 / (1 + permutations)\\ from vegan, where its
p-value is about 1/3. So ecodive's p-values are never smaller than
vegan's on the same permutations, and are equal for taxa that no
permutation ties. Each pair has random permutations of its own, drawn
from `seed` and the pair's position, as the pairwise tests of
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
and the others have. So `p_value`s are identical for any number of
`cpus`, but differ by sampling error from `simper()` run by hand on the
pair's samples with the same `seed`. A supplied `permutations` matrix is
cut down to each pair's samples by following each permutation's cycles,
as in the pairwise tests (see
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)'s
**Pairwise comparisons**). Give `vegan::simper()` those permutations, on
the pair's samples, for the same permutations.

The adjusted p-values (`p_adjusted`) are an addition: `vegan::simper()`
doesn't adjust. A study of a few groups and hundreds of taxa tests many
hypotheses, so look at `p_adjusted`, not `p_value`.

**Interpretation: location and dispersion**

Warton, Wright & Wang (2012) show that SIMPER confounds differences in
mean abundance (location) with differences in variability (dispersion):
the most variable taxa get the largest contributions even when their
mean abundance is the same in both groups, so a high `average` doesn't
mean a taxon differs between the groups. Check the groups' dispersions
with
[`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md),
compare `mean_a` with `mean_b`, and test the differences of interest
taxon by taxon with models of abundance (such as negative binomial
GLMs), rather than reading the ranking as evidence. The p-values test a
change in a taxon's contribution, which a change in its variance alone
can also give.

**Speed and memory**

A taxon's contribution is 0 for pairs of samples that both lack it, so
each taxon costs about \\n/2\\ operations per sample that has it, and
sparse tables are never made dense. The observed sums run in parallel
over taxa, and the permutations in parallel over permutations; sums are
always in the same order, so results are identical for any `cpus`.

**Differences from vegan::simper()**

- The result is one data frame (as
  [`summary()`](https://rdrr.io/r/base/summary.html) of vegan's result,
  with its pairs stacked), with the pair and taxon as columns, and
  vegan's `ava`, `avb`, and `cusum` named `mean_a`, `mean_b`, and
  `cumsum`.

- Pairs are named `"a-b"` and come in the order of the groups' factor
  levels; vegan names them `"a_b"`, in the order the groups first
  appear.

- vegan permutes all the samples' groups for every pair (a supplied
  matrix covers all of them), so its p-values for a pair also depend on
  the other groups' samples. ecodive permutes only the pair's samples,
  as its pairwise tests do.

- p-values count permutations whose `average` ties the observed one, as
  Clarke (1993) does; vegan counts only larger ones (see
  **Permutations**).

- `p_adjusted`, `norm`, and `strata` given directly.

## Same names in other packages

`vegan::simper()` has the same name: whichever package was attached last
masks the other. Call `ecodive::simper()` or `vegan::simper()` to be
sure.

## References

Clarke, K. R. (1993). Non-parametric multivariate analyses of changes in
community structure. *Australian Journal of Ecology*, 18, 117-143.
[doi:10.1111/j.1442-9993.1993.tb00438.x](https://doi.org/10.1111/j.1442-9993.1993.tb00438.x)

Warton, D. I., Wright, S. T., & Wang, Y. (2012). Distance-based
multivariate analyses confound location and dispersion effects. *Methods
in Ecology and Evolution*, 3, 89-101.
[doi:10.1111/j.2041-210X.2011.00127.x](https://doi.org/10.1111/j.2041-210X.2011.00127.x)

## See also

[`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md),
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md)

## Examples

``` r
    # Oral (Saliva, Gums) versus other (Nose, Stool) samples.
    grp <- c('oral', 'oral', 'other', 'other')
    # Four samples are far too few for permutation p-values.
    res <- simper(ex_counts, grp, permutations = 0, cpus = 1)
    res
#>         pair             taxon    average         sd     ratio mean_a mean_b
#> 1 oral-other     Streptococcus 0.30125815 0.19946484 1.5103321  477.5   11.5
#> 2 oral-other       Bacteroides 0.25995659 0.31393206 0.8280664    3.0  306.5
#> 3 oral-other   Corynebacterium 0.15787107 0.18628159 0.8474862    0.0  249.5
#> 4 oral-other       Haemophilus 0.10495737 0.06639368 1.5808337  133.5    1.5
#> 5 oral-other Propionibacterium 0.07946518 0.09327600 0.8519359    1.0  125.5
#> 6 oral-other    Staphylococcus 0.07474069 0.08811542 0.8482136    0.5  118.5
#>      cumsum p_value p_adjusted
#> 1 0.3079565      NA         NA
#> 2 0.5736931      NA         NA
#> 3 0.7350744      NA         NA
#> 4 0.8423654      NA         NA
#> 5 0.9235975      NA         NA
#> 6 1.0000000      NA         NA

    # The averages sum to the mean between-group Bray-Curtis dissimilarity.
    sum(res$average)
#> [1] 0.9782491
    mean(as.matrix(bray(ex_counts))[1:2, 3:4])
#> [1] 0.9782491
```
