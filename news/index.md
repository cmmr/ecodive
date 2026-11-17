# Changelog

## ecodive 2.4.0

- New [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) runs
  a principal coordinates analysis (classical multidimensional scaling)
  on a `dist` object, such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).
  It gives the same coordinates and eigenvalues as
  [`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html), `ape::pcoa()`,
  and `vegan::wcmdscale()` (up to each axis’s sign), but computes only
  the first few axes (`k = 3` by default), in parallel C, without
  forming an n-by-n matrix. Its results are identical for any number of
  `cpus`. `correction = "lingoes"` or `"cailliez"` adds the smallest
  constant to the squared distances, or to the distances, that leaves no
  negative eigenvalue, giving the same constants, eigenvalues and
  coordinates as `ape::pcoa()`, `vegan::wcmdscale(add = )` and
  `cmdscale(add = TRUE)`, again without an n-by-n (or 2n-by-2n) matrix:
  with 2,000 samples, `"cailliez"` takes 20 seconds where
  `cmdscale(add = TRUE)` takes 2 minutes. Every result has `correction`
  and `constant` attributes.

- New
  [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md)
  projects the features of a table (taxa counts, or environmental
  variables) onto the axes of a
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) result,
  for drawing arrows on a PCoA plot. It computes the same scores as
  scikit-bio’s
  [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md),
  in parallel C, keeps sparse tables sparse, and takes every table type
  and `norm` that
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) does.
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) gives
  the same directions with permutation tests.

- New
  [`ordination()`](https://cmmr.github.io/ecodive/reference/ordination.md)
  runs any ordination method by name, as
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  does for beta diversity metrics, with the method’s own arguments.
  [`list_ordinations()`](https://cmmr.github.io/ecodive/reference/ordinations.md)
  and
  [`match_ordination()`](https://cmmr.github.io/ecodive/reference/ordinations.md)
  list and look up the methods. Every ordination method returns the same
  kind of result: a sample-by-axis matrix (of class `"ecodive_ord"`,
  which prints its attributes’ names rather than their contents) whose
  attributes start with `method`, `eig`, and `prop`. Every method
  rejects negative distances (beyond rounding error), and `seed` can be
  given to any of them through
  [`ordination()`](https://cmmr.github.io/ecodive/reference/ordination.md).

- New [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) runs
  non-metric multidimensional scaling with the SMACOF algorithm of
  scikit-learn’s non-metric `MDS()`, ported to parallel C. It starts
  from the PCoA solution and from `n_starts` random configurations,
  keeps the one with the lowest Kruskal stress-1, and centers it,
  rotates it to principal axes, and scales it to the distances. Tied
  distances use Kruskal’s primary approach (as in `vegan::metaMDS()`),
  or scikit-learn’s secondary one with `ties = "secondary"`. The
  distances are sorted once and shared by all starts; nothing n-by-n is
  formed. Each start runs at most `max_iters` iterations, or stops when
  the stress improves by less than `tol`. Random starts come from a
  `seed` argument, and results are identical for any number of `cpus`.
  The result’s attributes are `stress`, `converged`, `n_iters`,
  `n_starts`, `best_start`, `n_best`, and `ties`.

- New [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) runs
  metric multidimensional scaling with the SMACOF algorithm of
  scikit-learn’s metric `MDS()`, fitting the distances themselves rather
  than their ranks. It shares
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md)’s engine,
  starts, arguments (without `ties`) and result format, and has no limit
  on the number of samples. `ordination(dm, "mds")` and `"metric MDS"`
  select it; `"classical MDS"`, `"pco"` and `"cmdscale"` select
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md).

- New [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md)
  runs an Isomap (isometric feature mapping) ordination: it links each
  sample to its `n_neighbors` nearest neighbors (or to those within a
  `radius`), replaces each distance with the shortest path through that
  graph, and runs
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) on those
  geodesic distances, which keeps long ecological gradients straight
  where [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)
  bends them into an arch. A C port of scikit-learn’s
  `Isomap(metric = "precomputed")`, it gives the same geodesic distances
  and embedding to rounding, breaks ties by sample order, and joins
  disconnected neighbor graphs at their closest samples, with a warning
  (`vegan::isomap()` stops instead; note that its `k` is the number of
  neighbors, while
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md)’s `k`
  is the number of axes). The neighbor and shortest-path searches run in
  parallel without an n-by-n matrix, and results are identical for any
  number of `cpus`.

- New [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md) runs
  t-distributed stochastic neighbor embedding (t-SNE) on a `dist`
  object, as a port of scikit-learn’s `TSNE` to parallel C.
  `gradient = "auto"` uses the exact gradient for up to 500 samples (or
  more than 3 dimensions) and the Barnes-Hut approximation (a quadtree
  in 2 dimensions, an octree in 3) above that. It finds nearest
  neighbors exactly from the distances and never forms an n-by-n matrix.
  It starts from the PCoA solution by default, sets the learning rate
  with scikit-learn’s rule (`learning_rate = "auto"`), stops early once
  the layout stops improving (at most `max_iters` iterations; the number
  run is the `n_iters` attribute), and lowers `perplexity` for small
  datasets. Its layouts match scikit-learn’s and `Rtsne::Rtsne()`’s in
  quality. Random starts come from a `seed` argument, and results are
  identical for any number of `cpus`.

- New [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md) runs
  a uniform manifold approximation and projection (UMAP) on a `dist`
  object, following `uwot::umap()` and its defaults: exact nearest
  neighbours, smooth-kNN calibration and fuzzy union (identical to
  uwot’s), a spectral start from the normalized graph Laplacian, and the
  layout optimization of uwot’s batch mode (`batch = TRUE`: per-epoch
  gradients, then Adam), which gives uwot’s coordinates from the same
  start and seeds. Every stage runs in parallel C, the neighbours
  straight from the dist triangle without an n-by-n matrix. Random
  numbers come from a `seed` argument, and results are identical for any
  number of `cpus`. `min_dist` defaults to the canonical 0.1
  (umap-learn’s), not uwot’s 0.01.

- New
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  runs a permutational multivariate analysis of variance on a `dist`
  object and sample metadata. Like ecodive’s other tests, it takes
  `groups`: one value per sample (categories or a numeric covariate), or
  a column name of `data`. For several variables, give `data` and a
  `formula` instead: factors, covariates, interactions, and
  `Condition()` terms. It gives the same table as `vegan::adonis2()`
  (Df, sums of squares, R2, F, and p-values), testing the model as a
  whole or each term sequentially, marginally, or by contrast, with
  optional `strata`. Given the same permutation matrix, its permuted F
  statistics and p-values are identical to `adonis2()`’s. It returns a
  list of class `"ecodive_test"`, like
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md): the
  model’s F, p-value, z-score (as `summary(vegan::permustats())`) and
  number of permutations come first, whatever `by` is, then its R2 and
  `adonis2()`’s table, as `table`. The sums of squares come straight
  from the dist triangle, without the n-by-n matrices `adonis2()` builds
  for each permutation, and permutations run in parallel. On one thread
  it is 25 to 50 times faster than `adonis2()` with 200 samples, and 100
  to 300 times faster with 1,000. Permutations come from a `seed`
  argument, and results are identical for any number of `cpus`. Like
  `adonis2()`, it also takes the distances on the left-hand side of a
  formula: `permanova(dm ~ treatment + ph, data)`.

- New [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) runs
  a multi-response permutation procedure on a `dist` object and each
  sample’s group, as `vegan::mrpp()`: the observed and expected weighted
  mean within-group distance (delta), the chance-corrected within-group
  agreement A, each group’s mean distance, and a permutation p-value,
  with group weights `'n'`, `'n-1'` or `'n(n-1)'` (vegan’s `weight.type`
  1 to 3) and optional `strata`. It also gives the standardized effect
  size z of `summary(vegan::permustats())`. Given the same permutation
  matrix, its permuted deltas, p-values and z are those of
  `vegan::mrpp()`. Each permutation costs one comparison per pair of
  samples, read straight from the dist triangle; with 500 samples and
  999 permutations it is about 60 times faster than `vegan::mrpp()` on
  one thread. Results are identical for any number of `cpus`. Like
  ecodive’s other permutation tests, it returns a list of class
  `"ecodive_test"` whose `statistic`, `p_value`, `z` and `n_perm` come
  first, and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) turns
  it into one row with the same columns for every test.

- New [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md)
  runs an analysis of similarities (ANOSIM) on a `dist` object and each
  sample’s group, given as a vector or as a column of sample metadata,
  with optional `strata`. Given the same permutation matrix, its R
  statistic, permuted R values and p-value match `vegan::anosim()`’s. It
  ranks the distances once, with a radix sort, and then each permutation
  needs only one pass over the ranks, without the n-by-n matrices
  `vegan::anosim()` builds; permutations run in parallel. With 300
  samples and 999 permutations it takes under a tenth of a second,
  against about seven for `vegan::anosim()`. Results are identical for
  any number of `cpus`.

- New
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
  tests the homogeneity of multivariate dispersions (PERMDISP), as
  `vegan::betadisper()` followed by `vegan::permutest()`: each sample’s
  distance to its group’s centroid (the default) or spatial median, each
  group’s mean distance, and the ANOVA F with a permutation p-value,
  plus optional pairwise t tests and small-group bias adjustment.
  Distances to centroids come from the distances within groups alone,
  with no eigendecomposition: they match `betadisper()`’s to rounding
  error, and are about 100 times faster with 1,000 samples.
  `betadisper()` defaults to spatial medians instead: pass
  `type = "median"` to compare with it.
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
  solves them exactly, where `betadisper()`’s optimizer stops at about
  1e-5. Given the same permutation matrix, its p-values are identical to
  `permutest()`’s. Permutations come from a `seed` argument, can be
  restricted to `strata`, and results are identical for any number of
  `cpus`.

- New [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)
  runs a Mantel test of the correlation between two distance matrices
  (Pearson, Spearman, or Kendall), or a partial Mantel test that removes
  the effect of a third (`given`). It gives the same statistic as
  `vegan::mantel()` and `vegan::mantel.partial()`, and, given the same
  permutation matrix, the same permuted statistics and p-values. The
  second matrix may also be a numeric vector, one value per sample (such
  as age), and samples are matched by name. It adds lower-tail and
  two-sided tests (`alternative`). Each permutation is one pass over the
  dist triangle, with no n-by-n matrix (Kendall’s tau-b by Knight’s
  merge sort rather than [`cor()`](https://rdrr.io/r/stats/cor.html)’s
  quadratic count), and permutations run in parallel. Permutations come
  from a `seed` argument, and results are identical for any number of
  `cpus`.

- New `dist_test(dm, test, groups, data, ...)` runs any of these tests
  by name, as
  [`ordination()`](https://cmmr.github.io/ecodive/reference/ordination.md)
  does for ordination methods:
  `dist_test(dm, "mrpp", "treatment", meta)` is
  `mrpp(dm, "treatment", meta)`, and other arguments (such as
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)’s
  `method`) pass through. It checks that the variable is of a kind the
  test takes: categorical groups for `anosim`, `mrpp` and `permdisp`,
  numeric values for `mantel` (passed as its `y`), and either for
  `permanova`.
  [`list_dist_tests()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
  and
  [`match_dist_test()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
  list and look up the tests, with forgiving name matching (`"adonis2"`
  finds `permanova`, `"betadisper"` finds `permdisp`).

- One rule aligns every per-sample argument (`groups`, `strata`, and
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)’s `y`
  and `given`) in
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md),
  [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md),
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md),
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md),
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) and
  [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md):
  with `data`, a column name or one value per row of `data`; without it,
  one value per sample, matched by name when named, else by position.
  The same checks and messages apply everywhere. As with
  `vegan::adonis2()`,
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  needs `groups` or a `formula` (use `formula = ~ .` to test every
  column of `data`). `permutations` is limited to 1e8.

- New [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)
  runs a distance-based redundancy analysis (db-RDA): the axes of a
  `dist` object that sample metadata explain, given as a formula such as
  `~ treatment + Condition(site)`. It gives the same eigenvalues,
  inertia, site scores, biplot scores and centroids as `vegan::dbrda()`
  (up to each axis’s sign; not `vegan::capscale()`, which differs for
  non-Euclidean distances), without forming an n-by-n matrix: the
  constrained axes come from a small eigenproblem, and unconstrained
  axes (`MDS1`, …) are only computed when the constraints give fewer
  than `k` axes. Its results are identical for any number of `cpus`. As
  in `vegan`, `dbrda(dm ~ treatment, data)` works too, and a formula of
  only `Condition()` terms gives the unconstrained ordination after
  removing their effect.

- [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) takes
  the same `groups`, `data` and `formula` as
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md):
  one variable as `groups` (a vector or a column of `data`), or `data`
  and a `formula`, or `dm ~ terms`, so
  `ordination(dm, "dbrda", "treatment", meta)` works. Its results gain
  `r_squared` and `adj_r_squared`, as `vegan::RsquareAdj()`, including
  for partial models.

- New
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  runs a constrained analysis of principal coordinates. It gives the
  same results as `vegan::capscale()` (and
  `phyloseq::ordinate(method = "CAP")`), with any imaginary inertia from
  negative eigenvalues reported in `inertia`. Euclidean distances take
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)’s fast
  path; otherwise a threaded eigendecomposition makes it about five
  times faster than `vegan::capscale()` on one thread. Also available as
  `ordination(dm, "capscale", ...)`.

- New [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md)
  fits sample metadata (`data`) and feature tables (`counts`, such as
  taxa abundances) onto any ordination, for biplots: an arrow and
  r-squared for each numeric variable or feature, the centroid of each
  factor level, and permutation p-values. It gives the same arrows,
  r-squared and centroids as `vegan::envfit()`, and the same p-values
  for the same permutation matrix. Each permutation costs a few
  multiply-adds per non-zero value of the table, with no refitting, and
  sparse tables are never made dense, so thousands of taxa with 999
  permutations take well under a second (minutes in `vegan`). It fits
  every axis of the ordination by default (`vegan` fits the first two),
  and its results are identical for any number of `cpus`.

- New
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
  returns a count table after any of the transformations that
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) apply
  with their `norm` argument (`'percent'`, `'chord'`, `'binary'`,
  `'clr'`, `'rclr'`), with exactly the values they use. It keeps the
  table’s orientation, dimnames and form: sparse Matrix and slam input
  stays sparse. `norm = 'clr'` gives a dense matrix, since every zero
  gets a value; the functions with a `norm` argument avoid that copy.

- New [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) runs a
  principal component analysis of a count table, after any of
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)‘s
  normalizations (e.g. `norm = 'clr'` for an Aitchison PCA), optionally
  scaling features to unit variance (`scale = TRUE`). It gives the same
  scores, loadings, and variances as
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) on the normalized
  table (up to each axis’s sign), but computes only the first few axes,
  in parallel C, straight from sparse input without making it dense. Its
  scores equal
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)’s on the
  matching Euclidean or Aitchison distances; its `loadings` attribute
  holds the features’ weights. Results are identical for any number of
  `cpus`.

- New [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) runs
  a robust Aitchison PCA (Martino et al. 2019) straight from a count
  table: a robust centered log-ratio transform that treats zeros as
  missing, OptSpace matrix completion, and a PCA of the completed
  matrix. It returns sample scores, feature loadings, and the fit’s
  progress in the shared ordination format (with `eig` as variances, on
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)’s scale),
  and optionally gemelli’s distances between samples as a `dist` object
  (`distances = TRUE`). A C port of `gemelli`’s
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md), it gives
  the same results to within rounding error, but works only with the
  non-zero counts, so sparse tables stay sparse and nothing of size
  samples by features is formed. It takes every count-table class that
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) takes, and
  its results are identical for any number of `cpus`.

- New
  [`blas_threads()`](https://cmmr.github.io/ecodive/reference/blas_threads.md)
  gets or sets how many threads R’s BLAS/LAPACK library uses, while R is
  running. It works with OpenBLAS, Intel MKL, BLIS, and FlexiBLAS,
  returns `NA` where the count can’t be read (e.g. R’s reference BLAS or
  Apple’s Accelerate), and never raises an error, warning, or message.
  Setting is limited to what is safe: the count is never raised above
  the library’s own, and OpenMP’s thread count is left as it is. No
  other ecodive function calls it.

- Sparse counts with samples in rows (a `dgCMatrix` with `margin = 1`),
  or as triplets (`dgTMatrix`, slam’s `simple_triplet_matrix`), are now
  imported with counting sorts in linear time, instead of a quicksort: 4
  s became 0.06 s for a 5000 x 4000 `dgCMatrix` with 2.5 million
  non-zero values, and 26 s became 0.3 s for 9.4 million. Every alpha
  and beta diversity function benefits, and results are unchanged.
  Out-of-range indices in these inputs are now an error.

- Sparse counts: explicitly stored zeros are now ignored, and repeated
  entries of a `dgTMatrix` summed, as the Matrix package defines them,
  so every function gives the same results as for the dense matrix
  (before, `norm = 'binary'` counted a stored zero as present).
  Malformed sparse objects (wrong index types or lengths, bad column
  pointers) are now an error rather than read out of bounds, and tables
  with more than ~268 million non-zero values no longer overflow an
  internal size.

- Threads now come from a pool that lasts the whole session, and jobs
  are shared out to whichever threads are free. On a machine busy with
  other processes, multithreaded calls are no longer much slower than
  single-threaded ones
  ([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) of 2,000
  samples: 0.8 s with 6 threads, where restarting threads for every step
  took 16 s), and small jobs no longer start threads at all. Results are
  unchanged, and still the same for any `cpus`.

- Every function’s `cpus` is now capped at
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), the
  cores available, so that asking for more threads than that costs
  nothing. Results are the same for any `cpus`. Without the parallelly
  package,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md) now
  returns at most 2 under `R CMD check --as-cran`, as parallelly’s does.

- Fixed `NA` in an integer or logical count table being read as
  -2147483648: every alpha and beta diversity function now treats it as
  `NA`, as it does in a double table. Before,
  [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) returned
  -1 for every pair with that sample.

- Fixed the `pairs` argument of the beta diversity functions, which
  counted the columns of `counts` as samples even when samples are rows
  (`margin = 1`, the default). Out-of-range pairs could crash R, and a
  table with one feature was rejected as having fewer than two samples.

## ecodive 2.3.2

- Fixed
  [`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
  and every beta metric called with `norm = 'clr'`, returning `NaN` (or,
  for
  [`chebyshev()`](https://cmmr.github.io/ecodive/reference/chebyshev.md),
  wrong values) for every sample pair. Two bugs were involved:

  - Data with no zeros (pseudocount 0) hit `log(0) * 0`, which is `NaN`.
  - The pseudocount was truncated to an integer in C, so values below 1,
    including the automatic `min(counts[counts > 0]) / 2`, became 0, and
    values like 1.5 were silently rounded down.

- [`gower()`](https://cmmr.github.io/ecodive/reference/gower.md) with
  `norm = 'clr'` now computes each feature’s range, and the distance for
  features absent from both samples, from that feature’s CLR values
  instead of treating absent features as 0.
  [`gower()`](https://cmmr.github.io/ecodive/reference/gower.md) with
  negative inputs (e.g. `norm = 'rclr'`) now also uses the full
  min-to-max range.

## ecodive 2.3.1

CRAN release: 2026-11-05

- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)
  now handles quoted labels per the Newick standard:

  - `''` inside a quoted label is decoded as one literal apostrophe, so
    `'it''s'` is read as `it's`.
  - Square brackets inside quotes are literal characters rather than
    comments, so labels like `'o__[Saprospirales]'` are no longer
    truncated.
  - Newlines and whitespace inside quotes are preserved.
  - Comments (`[...]`) outside of quotes may contain apostrophes.

- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)
  raises an informative error, with the position of the problem, for
  malformed input such as an unclosed quote, comment, or parenthesis.
  Previously these could return a garbled tree.

- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)’s
  parser is now linear-time and non-recursive, making it about 4x faster
  on a typical 100,000-tip tree and much faster on very deep trees.

## ecodive 2.3.0

- Added `ecodive` JOSS citation information.
- New
  [`robust_aitchison()`](https://cmmr.github.io/ecodive/reference/robust_aitchison.md)
  beta diversity function.
- Switch from `testthat` to dependency-free `tinytest` for test suite.
- `read_tree` fixes for edge cases.

## ecodive 2.2.6

CRAN release: 2026-04-14

- When using CLR normalization:

  - Warns when zeros are detected and `pseudocount` is `NULL`.
  - Then continues with `pseudocount <- min(counts[counts>0]) / 2`.

### Breaking Changes

- `norm` now defaults to `"none"` instead of `"percent"`.
- `norm` parameter was removed from some beta diversity functions.

## ecodive 2.2.4

### Breaking Changes

- The [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md)
  function has been simplified:

  - `drop = FALSE` now returns samples with original counts when sum \<
    depth.
  - `drop = TRUE` no longer drops zero-sum features - only zero-sum
    samples.
  - `depth` can no longer be a fraction (only positive integers now).
  - `n_samples` parameter has been removed.
  - `warn` parameter has been added.

## ecodive 2.2.3

- Coerce counts to double before rarefaction.

## ecodive 2.2.2

CRAN release: 2026-01-16

- Drop explicit zeros from sparse matrices after rarefaction.

## ecodive 2.2.1

CRAN release: 2025-12-12

- Corrected CRAN `rchk` issue.

## ecodive 2.2.0

CRAN release: 2025-12-03

- Dedicated rarefaction backends for dense and sparse matrices.
- Auto-conversion to col-major sparse matrix for distance/diversity
  calculations.

## ecodive 2.1.0

CRAN release: 2025-09-26

### Breaking Changes

- `rescale` parameter is superseded by `norm` parameter.

### Fixes

- [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
  and
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  will only forward applicable parameters.
- No longer crashes when `pairs = integer(0)`.
- `pseudocount = NULL` correctly selects smallest non-zero count.

## ecodive 2.0.0

CRAN release: 2025-09-21

### Breaking Changes

- Input matrix is now samples as rows.
- Removed `weighted` parameter.

### Added New Diversity Methods

- Updated Alpha Options:
  `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`

- Updated Beta Options:
  `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`

## ecodive 1.0.0

CRAN release: 2025-08-22

- This is the first release of ecodive.
