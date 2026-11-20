# `cross_assoc()`: associations between the features of two tables

**Priority:** 39. **Tier:** formula, with SciPy's p-value formulas as the
port (`scipy.stats.pearsonr`, `spearmanr`, `kendalltau`; SciPy 1.x, **BSD-3**,
`LICENSE.txt` "Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy
Developers", checked 2026-11-20). `stats::cor()` / `cor.test()` (R, GPL-2),
Hmisc 5.3-1 `rcorr()` (**GPL (>= 2)**, CRAN `DESCRIPTION` checked
2026-11-20), psych 2.6.9 `corr.test()` (**GPL (>= 2)**, checked 2026-11-20)
and microeco 2.4.3 `trans_env$cal_cor()` (**GPL-3**, checked 2026-11-20) are
output-only: their docs and outputs may be used, their R source may not.
Read [README.md](README.md) first, then [sparcc.md](sparcc.md) and
[assoc_network.md](assoc_network.md) for the one-table association
conventions this plan extends. [halla.md](halla.md) builds on this item.


## What it is

The "correlation heatmap between genera and environmental factors" (or
metabolites, host genes, clinical variables): for every feature of one table
and every variable of another, measured on the same samples, a correlation,
its p-value and an FDR-adjusted p-value, drawn as a heatmap with stars on the
significant cells. It is one of the most common figures in microbiome papers:
"Spearman" AND "heatmap" AND ("environmental factors" OR "metabolites") has
**7,120** Europe PMC microbiome full-text mentions 2022–26 (retrieved
2026-11-20). Users now make it with `Hmisc::rcorr()` or `psych::corr.test()`
plus pheatmap, or microeco's `trans_env$cal_cor()`, all of which compute the
full (p + q) × (p + q) matrix when only the p × q block is wanted, and none of
which is threaded.

ecodive's `feature_assoc()` (`FEATURE_ASSOCS`), `rho_p()`, `sparcc()` and
`assoc_pairs()` all work **within one table** and return a `dist` over its
features. Two tables give a rectangular p × q result, so this is a new
function, not a new `FEATURE_ASSOCS` row (see "Names").


## Names

- **`cross_assoc()`** (recommended): associations *across* two tables. It
  parallels `feature_assoc()` (`<what>_assoc`), reads as what it does, and
  no registry has an alias starting `cross` (check with README's snippet,
  adding `FEATURE_ASSOCS` and `DA_TESTS`). Alternatives: `assoc_cross()`
  (sorts with `assoc_pairs()` / `assoc_network()`, rule 2), `cor_tables()`.
  **Ask the user.**
- **Not a `feature_assoc()` member.** `FEATURE_ASSOCS` rows mean "a `dist`
  over the features of one table": `dist_clust_order()`,
  `plot.ecodive_dist()` and `assoc_network()` rely on that. A p × q table
  has no diagonal, no packed triangle and no clustering of one feature set.
  It takes two inputs, so per README it is a plain function with no
  wrapper, free of the wrapper-argument constraints.
- **Same names in other packages:** none known.


## Design

### Signature

```r
cross_assoc(x, y, data = NULL, method = 'spearman', margin = 1L,
            norm = 'percent', norm_y = 'none', pseudocount = NULL,
            given = NULL, permutations = 0L, strata = NULL,
            alpha = 0.05, p_adjust = 'BH', min_samples = 3L,
            seed = 0L, cpus = n_cpus())
```

- **`x`:** a count table (any type `validate_counts()` reads; `margin` as
  everywhere). It may also be any numeric table (metabolites, host gene
  expression) with `norm = 'none'`.
- **`y`:** the second table, **samples in rows**: a numeric matrix, a data
  frame of numeric columns (environmental or clinical variables), or a
  count table. A character vector of column names of `data` selects those
  columns. Logical columns become 0/1. Factor or character columns are an
  error: "`y` column 'site' is not numeric; use da_test() for groups". (A
  two-level factor as 0/1 gives the point-biserial correlation; allow it
  only if the user asks.) Samples are matched to `x` by name when both are
  named (`match_rows()`, as `envfit()` does), else by position.
- **`method`:** `'spearman'` (default, as Hmisc, microeco and HAllA),
  `'pearson'`, `'kendall'` (tau-b), or `'rho_p'` (below). The glossary's
  meaning of `method` (a correlation coefficient, as `mantel()`).
- **`norm`, `norm_y`, `pseudocount`:** how each table is transformed first
  (`normalize()`'s norms: `'none'`, `'percent'`, `'clr'`, `'rclr'`,
  `'log'`, ...). Default `'percent'` for `x` (relative abundances, what
  papers correlate) and `'none'` for `y`. Spearman and Kendall are
  invariant to monotone *per-feature* transforms but not to per-sample ones,
  so `'percent'` versus raw counts matters; the docs say so. `norm_y` is a
  new name (add it to parameter_names.md); **ask the user** whether they
  prefer `y_norm`.
- **`method = 'rho_p'`** (the "compositional side"): needs both tables to be
  count tables (or compositions). Each is clr-transformed **within its own
  table** (its own geometric mean; `pseudocount` as `rho_p()`), then
  \eqn{\rho_p(i,j) = 1 - \mathrm{var}(z^x_i - z^y_j) /
  (\mathrm{var}(z^x_i) + \mathrm{var}(z^y_j))}, Lin's concordance of the
  two clr values (`rho_p()`'s formula, see sparcc.md "Why it is not
  rho_p()"). It is not in propr, which works within one composition: the
  two clr references differ, so this measures concordance of each part
  relative to its own table's centre. **Ask the user** whether to include
  it; `norm = 'clr', norm_y = 'clr', method = 'pearson'` is the more
  standard compositional choice and needs nothing new. No p-value
  formula exists for `rho_p`: p-values only with `permutations`.
- **`given`:** covariates to partial out (a one-sided formula over `data`,
  or a numeric data frame), giving **partial** correlations: both tables'
  values (ranks, for Spearman) are residualised on `cbind(1, given)` by one
  QR (`condition_basis()`), and Pearson's r of the residuals is reported,
  with \eqn{n - 2 - k} degrees of freedom (k = columns of `given`). This is
  `ppcor::pcor.test()`'s definition (ppcor 1.1, **GPL-2**, CRAN
  `DESCRIPTION` checked 2026-11-20: output-only). Not available for
  Kendall (error).
- **`permutations`, `strata`, `seed`:** optional permutation p-values: the
  rows of `y` are permuted (`perm_matrix()`, `src/permute.c`, so a supplied
  matrix reproduces another tool's permutations exactly), the **same**
  permutation for every pair, which keeps each table's correlation
  structure. With `given`, the residualised `y` is permuted (ter Braak's
  approach). p = (1 + #{|r_b| ≥ |r| − √ε})/(B + 1) via `perm_p()`.
- **`alpha`, `p_adjust`:** adjustment over **all** p × q pairs
  (`stats::p.adjust()`, default `'BH'`, as Hmisc users do by hand, microeco's
  default `p_adjust_method = "fdr"`, and HAllA); `alpha` makes the `reject`
  column and the plot's default marks.
- **`min_samples`:** pairs with fewer complete samples get `NA`.

### Missing values

`x` (counts) has none. `y` (metadata) often does. Use **pairwise complete
observations**, as `rcorr()` and `cor(use = 'pairwise.complete.obs')` do:
pair (i, j) uses the samples where `y[, j]` is not `NA` (and `given` is
complete). Ranks of `x[, i]` must then be recomputed on that subset. Group
the columns of `y` by their missingness pattern (hash of the logical
vector), and rank `x` once **per pattern**; a typical metadata table has a
handful of patterns. Record `n` per pair. `given` with `NA`: drop those
samples from everything (warn with the count).

### Algorithm

1. Validate, match samples, transform each table (`normalize()` in C).
2. **Ranks once per column** (Spearman): average ranks for ties, as
   `rank()` and `scipy.stats.rankdata`. Then each column is centred and
   scaled to unit norm, so \eqn{r_{ij} = \sum_s \tilde x_{si} \tilde y_{sj}},
   a product of a p × n and an n × q matrix. Pearson: the same on the
   values.
3. **Sparse `x`** (Spearman, Pearson): every zero of a column has the same
   value (rank \eqn{c_i}, the average rank of the zeros), and \eqn{\sum_s
   \tilde y_{sj} = 0}, so \eqn{r_{ij} = \sum_{s: x_{si} > 0} (\tilde x_{si}
   - \tilde c_i) \tilde y_{sj}}: a sum over the non-zeros only. Never
   densify `x`; cost O(nnz × q).
4. **Kendall tau-b:** per pair, Knight's (1966) O(n log n) algorithm on `y`
   values sorted by `x` order: `kendall_sort_count()` and
   `kendall_tied_sorted()` already exist (`src/mantel.c`); move them to a
   shared header if needed. Each `x` column's order and tie counts are
   computed once.
5. **p-values** (in R, vectorised, after C returns `r`, `n` and for Kendall
   the tie sums):
   - Pearson, Spearman: \eqn{t = r\sqrt{(n-2-k)/(1-r^2)}}, two-sided
     Student t with \eqn{n - 2 - k} df: `scipy.stats.pearsonr`'s and
     `spearmanr`'s formula, and Hmisc's and psych's. **`cor.test(method =
     'spearman')` differs**: for n < 1290 without ties it uses the exact
     AS 89 distribution (`exact = NULL`). Document; fixtures from
     `cor.test()` use `exact = FALSE`. |r| = 1 gives p = 0.
   - Kendall: normal approximation with the tie-corrected variance (Kendall
     1970), as `scipy.stats.kendalltau(method = 'asymptotic')` (port that
     formula, BSD-3). `cor.test(method = 'kendall')` is exact for n < 50
     without ties: document, and use `exact = FALSE` in fixtures.
   - `rho_p`: only by permutation.
6. **Adjust** over all pairs with `p.adjust()`; `reject = p_adjusted <=
   alpha`.

### Threads, memory, determinism

- Phase 1: one task per column (rank, centre, scale): `x` columns per
  missingness pattern, then `y` columns.
- Phase 2: tiles of 64 `x` features × all `y` columns; each dot product
  sums over samples in sample order, so results are `identical()` for any
  `cpus`. No BLAS: `dgemm`'s summation order depends on the library and its
  threads (see `blas_threads()`).
- Permutations: one task per (tile, block of permutations); each task
  applies its permutations to the standardised `y` and keeps per-pair hit
  counts (integers), reduced in task order.
- Memory: the standardised `y` (n × q doubles), `x` sparse or n × p
  doubles, and the p × q result. 1,000 genera × 500 metabolites is 500,000
  pairs: 4 MB per numeric column of the result.

### Result

A data frame of class `c('ecodive_cross_assoc', 'data.frame')`, one row per
pair in `x`-major order (all of `y` for the first `x` feature, then the
next), with columns:

- `x`, `y`: the feature and variable names (positions as text when
  unnamed);
- `association`: the correlation (named as in `assoc_pairs()`);
- `statistic`: t (Pearson, Spearman) or z (Kendall); absent for `rho_p`;
- `p_value`, `p_adjusted`, `reject`;
- `n`: the samples used (only when it varies between pairs).

Attributes (README rule 1 and 2): `method` (`'spearman'`, ...), `alpha`,
`p_adjust`, `n_samples`, `given` (the names of the partialled covariates,
since they qualify the coefficient), `n_perm` when permuted, and
`x_features`, `y_features` (each table's feature names in input order, so
`as.matrix()` and `plot()` can lay out the p × q grid without guessing).
Print with `print_table_result` and a `print_df_label()` entry, e.g.
"Spearman correlations: 120 taxa x 12 variables, 37 pairs with BH-adjusted
p <= 0.05".

`as.matrix.ecodive_cross_assoc(x, what = 'association')` returns the p × q
matrix of any column (`'p_adjusted'`, ...), rows `x_features`, columns
`y_features`.


## Reference implementations (all output-only except SciPy)

- `stats::cor()`, `cor.test(exact = FALSE)`: Pearson, Spearman, Kendall
  values and p-values.
- `Hmisc::rcorr(type = 'spearman')`: values, pairwise n and p (on the full
  combined matrix; take the p × q block). Check from its outputs how
  precisely it reports very small p-values, and compare it only where
  p > 1e-8.
- `psych::corr.test(x, y, method = 'spearman', adjust = 'BH')`: the
  two-table route; its adjusted p-values are over the p × q block.
- microeco `trans_env$cal_cor(use_data = 'genus', cor_method = 'spearman')`:
  the microbiome front end; check its `p_adjust_type` grouping (it can
  adjust per environmental variable rather than overall; document if so).
- `ppcor::pcor.test(method = 'spearman')`: partial correlations (GPL-2).
- **SciPy** (BSD-3) `scipy.stats.spearmanr`, `pearsonr`,
  `kendalltau(method = 'asymptotic')`: values and p-values, in
  `data-raw/expected/cross_assoc.py`.
- Fixtures: `data-raw/expected/cross_assoc.r` (40 samples × 25 taxa from
  `ex_counts`-like data, × 6 environmental variables with some `NA`, a
  constant column, ties), stored in `tests/tinytests/expected/cross_assoc.rds`.
  Match values at 1e-12 and p-values at 1e-10.


## Plotting

`plot.ecodive_cross_assoc(x, n_x = 30, n_y = 30, clust = 'average',
stars = c(0.05, 0.01, 0.001), col = NULL, key = TRUE, ...)`: the figure.

- Rows are `x` features, columns `y` variables. With more than `n_x` (or
  `n_y`) features, keep those with the smallest `p_adjusted` (ties by
  largest |r|), and say "30 of 120 taxa shown" in the subtitle.
- Order: `clust = 'average'` clusters rows and columns by average linkage on
  \eqn{1 - r} between their profiles in the shown block (`stats::hclust()`,
  computed from the result, never from the data); `'none'` keeps input
  order.
- Colours: the diverging `"Blue-Red 3"` scale from −1 to 1, as
  `plot.ecodive_dist()` draws similarities; the key via
  `dist_plot_key()`.
- Marks: `*`, `**`, `***` for `p_adjusted` below each of `stars` (the
  convention of these figures); `stars = NULL` marks `reject` with a dot.
  The subtitle says which p-values were starred ("BH-adjusted").
- Title: the display name of `method` from a small internal table
  (`CROSS_ASSOCS`: `spearman` "Spearman Correlation", `pearson`, `kendall`
  "Kendall's tau-b", `rho_p` "Proportionality (rho_p)"; "Partial ..." when
  `given` is set). Count tables have no titles, so there is no `display`
  attribute.
- Axis labels from `xlab`/`ylab` (default none; feature names label the
  rows and columns). Returns `x` invisibly; `par()` restored.

A bipartite network of the significant pairs (`assoc_network()` for two
tables) is a possible follow-up, not part of this item.


## Tests

- Fixtures above: Pearson, Spearman, Kendall values and p-values; partial
  Spearman; pairwise-complete `n` with `NA`s in `y`.
- The sparse path equals the dense path (`identical()` is not expected for
  sums in a different order: 1e-14), and dense input with zeros equals
  `cor(method = 'spearman')`.
- Invariances: per-feature monotone transforms of `y` don't change
  Spearman or Kendall; swapping `x` and `y` transposes the result.
- Permutations: a supplied permutation matrix gives p-values equal to a
  base-R loop over the same permutations; p in `[1/(B+1), 1]`.
- `identical()` for `cpus` 1, 2, 4; same `seed` same result.
- Edge cases: n = 3; a constant column (`NA` correlation, as `cor()`,
  with a warning); all-`NA` column; `y` with a factor column (error);
  sample names that don't match (error naming them); `rho_p` with a
  non-count `y` (error); Kendall with `given` (error).
- One smoke-test line (10 samples × 5 taxa × 3 variables).


## Docs

- Rd page in a new pkgdown subsection "Associations between tables" (with
  `halla()` later), linked from `feature_assoc()`, `rho_p()`, `sparcc()`.
- Explain: what the p-values assume (no multiple-testing structure; BH over
  all pairs); why relative abundances are compositional and that
  `norm = 'clr'` is the safer choice for taxon-taxon or taxon-metabolite
  pairs; that Spearman on relative abundances with many zeros is dominated
  by presence/absence.
- "Differences from cor.test(), Hmisc, psych, microeco": t-approximation
  p-values for Spearman (cor.test's exact AS 89); Kendall asymptotic; how
  NAs are handled; adjustment over all pairs.
- NEWS; `norm_y` in parameter_names.md.


## References

- Spearman, C. (1904). The proof and measurement of association between
  two things. *American Journal of Psychology*, 15, 72–101.
  doi:10.2307/1412159
- Kendall, M. G. (1970). *Rank Correlation Methods* (4th ed.). Griffin.
- Knight, W. R. (1966). A computer method for calculating Kendall's tau
  with ungrouped data. *JASA*, 61, 436–439.
  doi:10.1080/01621459.1966.10480879
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
  rate. *JRSS B*, 57, 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Lin, L. I.-K. (1989). A concordance correlation coefficient to evaluate
  reproducibility. *Biometrics*, 45, 255–268. doi:10.2307/2532051
- Erb, I., & Notredame, C. (2016). How should we measure proportionality
  on relative gene expression data? *Theory in Biosciences*, 135, 21–36.
  doi:10.1007/s12064-015-0220-8
- ter Braak, C. J. F. (1992). Permutation versus bootstrap significance
  tests in multiple regression and ANOVA. In *Bootstrapping and Related
  Techniques* (pp. 79–86). Springer.
- Virtanen, P., et al. (2020). SciPy 1.0. *Nature Methods*, 17, 261–272.
  doi:10.1038/s41592-019-0686-2
