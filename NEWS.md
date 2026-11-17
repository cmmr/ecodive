# ecodive 2.4.0

* New `pcoa()` runs a principal coordinates analysis (classical
  multidimensional scaling) on a `dist` object, such as from `beta_div()`. It
  gives the same coordinates and eigenvalues as `cmdscale()`, `ape::pcoa()`,
  and `vegan::wcmdscale()` (up to each axis's sign), but computes only the
  first few axes (`k = 3` by default), in parallel C, without forming an
  n-by-n matrix. Its results are identical for any number of `cpus`.

* New `ordination()` runs any ordination method by name, as `beta_div()`
  does for beta diversity metrics, with the method's own arguments.
  `list_ordinations()` and `match_ordination()` list and look up the methods.
  Every ordination method returns the same kind of result: a sample-by-axis
  matrix (of class `"ecodive_ord"`, which prints its attributes' names
  rather than their contents) whose attributes start with `method`, `eig`,
  and `prop`. Every method rejects negative distances (beyond rounding
  error), and `seed` can be given to any of them through `ordination()`.

* New `nmds()` runs a non-metric multidimensional scaling, following
  `vegan::metaMDS()` with its `monoMDS()` engine: a start from the PCoA
  solution, then random starts until the best solution is found again, then
  centering, rotation to principal components, and scaling. Its C port of
  `monoMDS` gives the same result as `vegan::monoMDS()` from the same starting
  configuration. Random starts run in parallel, one per thread, and share one
  sorted copy of the distances, so memory use is a fraction of `metaMDS()`'s.
  Random starts come from a `seed` argument, and results are identical for
  any number of `cpus`. Each start runs at most `max_iter` iterations.

* New `tsne()` runs t-distributed stochastic neighbor embedding (t-SNE) on a
  `dist` object, exactly (`theta = 0`) or with the Barnes-Hut approximation
  (`theta > 0`, the default). It follows `Rtsne::Rtsne()` given distances,
  and from the same starting coordinates gives exactly the same Barnes-Hut
  coordinates (up to each axis's sign) and KL divergence. The number of
  iterations is `n_iter`. It works straight from the distances, without an
  n-by-n matrix (Barnes-Hut keeps only each sample's nearest neighbors), and
  runs the neighbor search, perplexity calibration, and gradients in
  parallel C. As scikit-learn, openTSNE and FIt-SNE do, it starts from PCoA
  coordinates and sets the learning rate from the number of samples by
  default (`init = "random", learning_rate = 200` gives Rtsne's procedure).
  Random starts come from a `seed` argument, and results are identical for
  any number of `cpus`.

* New `umap()` runs a uniform manifold approximation and projection (UMAP)
  on a `dist` object, following `uwot::umap()` and its defaults: exact
  nearest neighbours, smooth-kNN calibration and fuzzy union (identical to
  uwot's), a spectral start from the normalized graph Laplacian, and the
  layout optimization of uwot's batch mode (`batch = TRUE`: per-epoch
  gradients, then Adam), which gives uwot's coordinates from the same start
  and seeds. Every stage runs in parallel C, the neighbours straight from
  the dist triangle without an n-by-n matrix. Random numbers come from a
  `seed` argument, and results are identical for any number of `cpus`.
  `min_dist` defaults to the canonical 0.1 (umap-learn's), not uwot's 0.01.

* New `permanova()` runs a permutational multivariate analysis of variance
  on a `dist` object and sample metadata. Like ecodive's other tests, it
  takes `groups`: one value per sample (categories or a numeric covariate),
  or a column name of `data`. For several variables, give `data` and a
  `formula` instead: factors, covariates, interactions, and `Condition()`
  terms. It gives the same table as `vegan::adonis2()` (Df, sums of squares,
  R2, F, and p-values), testing the model as a whole or each term
  sequentially, marginally, or by contrast, with optional `strata`. Given the
  same permutation matrix, its permuted F statistics and p-values are
  identical to `adonis2()`'s. It returns a list of class `"ecodive_test"`,
  like `mrpp()`: the model's F, p-value, z-score (as
  `summary(vegan::permustats())`) and number of permutations come first,
  whatever `by` is, then its R2 and `adonis2()`'s table, as `table`. The
  sums of squares come straight from the dist triangle, without the n-by-n
  matrices `adonis2()` builds for each permutation, and permutations run in
  parallel. On one thread it is 25 to 50 times faster than `adonis2()` with
  200 samples, and 100 to 300 times faster with 1,000.
  Permutations come from a `seed` argument, and results are identical for
  any number of `cpus`. Like `adonis2()`, it also takes the distances on the
  left-hand side of a formula: `permanova(dm ~ treatment + ph, data)`.

* New `mrpp()` runs a multi-response permutation procedure on a `dist`
  object and each sample's group, as `vegan::mrpp()`: the observed and
  expected weighted mean within-group distance (delta), the chance-corrected
  within-group agreement A, each group's mean distance, and a permutation
  p-value, with group weights `'n'`, `'n-1'` or `'n(n-1)'` (vegan's
  `weight.type` 1 to 3) and optional `strata`. It also gives the
  standardized effect size z of `summary(vegan::permustats())`. Given the
  same permutation matrix, its permuted deltas, p-values and z are those of
  `vegan::mrpp()`. Each permutation costs one comparison per pair of
  samples, read straight from the dist triangle; with 500 samples and 999
  permutations it is about 60 times faster than `vegan::mrpp()` on one
  thread. Results are identical for any number of `cpus`. Like ecodive's
  other permutation tests, it returns a list of class `"ecodive_test"`
  whose `statistic`, `p_value`, `z` and `n_perm` come first, and
  `as.data.frame()` turns it into one row with the same columns for every
  test.

* New `anosim()` runs an analysis of similarities (ANOSIM) on a `dist`
  object and each sample's group, given as a vector or as a column of sample
  metadata, with optional `strata`. Given the same permutation matrix, its R
  statistic, permuted R values and p-value match `vegan::anosim()`'s. It
  ranks the distances once, with a radix sort, and then each permutation
  needs only one pass over the ranks, without the n-by-n matrices
  `vegan::anosim()` builds; permutations run in parallel. With 300 samples
  and 999 permutations it takes under a tenth of a second, against about
  seven for `vegan::anosim()`. Results are identical for any number of `cpus`.
* New `permdisp()` tests the homogeneity of multivariate dispersions
  (PERMDISP), as `vegan::betadisper()` followed by `vegan::permutest()`: each
  sample's distance to its group's centroid (the default) or spatial
  median, each group's mean distance, and the ANOVA F with a permutation
  p-value, plus optional pairwise t tests and small-group bias adjustment.
  Distances to centroids come from the distances within groups alone, with
  no eigendecomposition: they match `betadisper()`'s to rounding error, and
  are about 100 times faster with 1,000 samples. `betadisper()` defaults to
  spatial medians instead: pass `type = "median"` to compare with it.
  `permdisp()` solves them exactly, where `betadisper()`'s optimizer stops
  at about 1e-5. Given the same permutation matrix, its p-values are
  identical to `permutest()`'s. Permutations come from a `seed` argument,
  can be restricted to `strata`, and results are identical for any number
  of `cpus`.

* New `mantel()` runs a Mantel test of the correlation between two distance
  matrices (Pearson, Spearman, or Kendall), or a partial Mantel test that
  removes the effect of a third (`given`). It gives the same statistic as
  `vegan::mantel()` and `vegan::mantel.partial()`, and, given the same
  permutation matrix, the same permuted statistics and p-values. The second
  matrix may also be a numeric vector, one value per sample (such as age),
  and samples are matched by name. It adds lower-tail and two-sided tests
  (`alternative`). Each permutation is one pass over the dist triangle, with
  no n-by-n matrix (Kendall's tau-b by Knight's merge sort rather than
  `cor()`'s quadratic count), and permutations run in parallel. Permutations
  come from a `seed` argument, and results are identical for any number of
  `cpus`.

* New `dbrda()` runs a distance-based redundancy analysis (db-RDA): the axes
  of a `dist` object that sample metadata explain, given as a formula such as
  `~ treatment + Condition(site)`. It gives the same eigenvalues, inertia,
  site scores, biplot scores and centroids as `vegan::dbrda()` (up to each
  axis's sign; not `vegan::capscale()`, which differs for non-Euclidean
  distances), without forming an n-by-n matrix: the constrained axes come from
  a small eigenproblem, and unconstrained axes (`MDS1`, ...) are only computed
  when the constraints give fewer than `k` axes. Its results are identical for
  any number of `cpus`. As in `vegan`, `dbrda(dm ~ treatment, data)` works
  too, and a formula of only `Condition()` terms gives the unconstrained
  ordination after removing their effect.

* New `envfit()` fits sample metadata (`data`) and feature tables
  (`counts`, such as taxa abundances) onto any ordination, for biplots: an
  arrow and r-squared for each numeric variable or feature, the centroid of
  each factor level, and permutation p-values. It gives the same arrows,
  r-squared and centroids as `vegan::envfit()`, and the same p-values for the
  same permutation matrix. Each permutation costs a few multiply-adds per
  non-zero value of the table, with no refitting, and sparse tables are never
  made dense, so thousands of taxa with 999 permutations take well under a
  second (minutes in `vegan`). It fits every axis of the ordination by
  default (`vegan` fits the first two), and its results are identical for any
  number of `cpus`.

* New `normalize()` returns a count table after any of the transformations
  that `beta_div()`, `pca()` and `envfit()` apply with their `norm` argument
  (`'percent'`, `'chord'`, `'binary'`, `'clr'`, `'rclr'`), with exactly the
  values they use. It keeps the table's orientation, dimnames and form:
  sparse Matrix and slam input stays sparse. `norm = 'clr'` gives a
  dense matrix, since every zero gets a value; the functions with a `norm`
  argument avoid that copy.

* New `pca()` runs a principal component analysis of a count table, after
  any of `beta_div()`'s normalizations (e.g. `norm = 'clr'` for an Aitchison
  PCA), optionally scaling features to unit variance (`scale = TRUE`). It
  gives the same scores, loadings, and variances as `prcomp()` on the
  normalized table (up to each axis's sign), but computes only the first few
  axes, in parallel C, straight from sparse input without making it dense.
  Its scores equal `pcoa()`'s on the matching Euclidean or Aitchison
  distances; its `loadings` attribute holds the features' weights. Results
  are identical for any number of `cpus`.

* New `blas_threads()` gets or sets how many threads R's BLAS/LAPACK library
  uses, while R is running. It works with OpenBLAS, Intel MKL, BLIS, and
  FlexiBLAS, returns `NA` where the count can't be read (e.g. R's reference
  BLAS or Apple's Accelerate), and never raises an error, warning, or
  message. Setting is limited to what is safe: the count is never raised
  above the library's own, and OpenMP's thread count is left as it is. No
  other ecodive function calls it.

* Sparse counts with samples in rows (a `dgCMatrix` with `margin = 1`), or as
  triplets (`dgTMatrix`, slam's `simple_triplet_matrix`), are now imported
  with counting sorts in linear time, instead of a quicksort: 4 s became
  0.06 s for a 5000 x 4000 `dgCMatrix` with 2.5 million non-zero values, and
  26 s became 0.3 s for 9.4 million. Every alpha and beta diversity function
  benefits, and results are unchanged. Out-of-range indices in these inputs
  are now an error.

* Sparse counts: explicitly stored zeros are now ignored, and repeated
  entries of a `dgTMatrix` summed, as the Matrix package defines them, so
  every function gives the same results as for the dense matrix (before,
  `norm = 'binary'` counted a stored zero as present). Malformed sparse
  objects (wrong index types or lengths, bad column pointers) are now an
  error rather than read out of bounds, and tables with more than ~268
  million non-zero values no longer overflow an internal size.

* Threads now come from a pool that lasts the whole session, and jobs are
  shared out to whichever threads are free. On a machine busy with other
  processes, multithreaded calls are no longer much slower than
  single-threaded ones (`pcoa()` of 2,000 samples: 0.8 s with 6 threads,
  where restarting threads for every step took 16 s), and small jobs no
  longer start threads at all. Results are unchanged, and still the same
  for any `cpus`.

* Every function's `cpus` is now capped at `n_cpus()`, the cores available,
  so that asking for more threads than that costs nothing. Results are the
  same for any `cpus`. Without the parallelly package, `n_cpus()` now
  returns at most 2 under `R CMD check --as-cran`, as parallelly's does.

* Fixed `NA` in an integer or logical count table being read as
  -2147483648: every alpha and beta diversity function now treats it as
  `NA`, as it does in a double table. Before, `bray()` returned -1 for every
  pair with that sample.

* Fixed the `pairs` argument of the beta diversity functions, which counted
  the columns of `counts` as samples even when samples are rows
  (`margin = 1`, the default). Out-of-range pairs could crash R, and a table
  with one feature was rejected as having fewer than two samples.



# ecodive 2.3.2

* Fixed `aitchison()`, and every beta metric called with `norm = 'clr'`,
  returning `NaN` (or, for `chebyshev()`, wrong values) for every sample
  pair. Two bugs were involved:

    - Data with no zeros (pseudocount 0) hit `log(0) * 0`, which is `NaN`.
    - The pseudocount was truncated to an integer in C, so values below 1,
      including the automatic `min(counts[counts > 0]) / 2`, became 0, and
      values like 1.5 were silently rounded down.

* `gower()` with `norm = 'clr'` now computes each feature's range, and the
  distance for features absent from both samples, from that feature's CLR
  values instead of treating absent features as 0. `gower()` with negative
  inputs (e.g. `norm = 'rclr'`) now also uses the full min-to-max range.



# ecodive 2.3.1

* `read_tree()` now handles quoted labels per the Newick standard:

    - `''` inside a quoted label is decoded as one literal apostrophe, so
      `'it''s'` is read as `it's`.
    - Square brackets inside quotes are literal characters rather than
      comments, so labels like `'o__[Saprospirales]'` are no longer truncated.
    - Newlines and whitespace inside quotes are preserved.
    - Comments (`[...]`) outside of quotes may contain apostrophes.

* `read_tree()` raises an informative error, with the position of the
  problem, for malformed input such as an unclosed quote, comment, or
  parenthesis. Previously these could return a garbled tree.

* `read_tree()`'s parser is now linear-time and non-recursive, making it
  about 4x faster on a typical 100,000-tip tree and much faster on very deep
  trees.



# ecodive 2.3.0

* Added `ecodive` JOSS citation information.
* New `robust_aitchison()` beta diversity function.
* Switch from `testthat` to dependency-free `tinytest` for test suite.
* `read_tree` fixes for edge cases.



# ecodive 2.2.6

* When using CLR normalization:

    - Warns when zeros are detected and `pseudocount` is `NULL`.
    - Then continues with `pseudocount <- min(counts[counts>0]) / 2`.

## Breaking Changes

* `norm` now defaults to `"none"` instead of `"percent"`.
* `norm` parameter was removed from some beta diversity functions.



# ecodive 2.2.4

## Breaking Changes

* The `rarefy()` function has been simplified:

    - `drop = FALSE` now returns samples with original counts when sum < depth.
    - `drop = TRUE` no longer drops zero-sum features - only zero-sum samples.
    - `depth` can no longer be a fraction (only positive integers now).
    - `n_samples` parameter has been removed.
    - `warn` parameter has been added.



# ecodive 2.2.3

* Coerce counts to double before rarefaction.



# ecodive 2.2.2

* Drop explicit zeros from sparse matrices after rarefaction.



# ecodive 2.2.1

* Corrected CRAN `rchk` issue.



# ecodive 2.2.0

* Dedicated rarefaction backends for dense and sparse matrices.
* Auto-conversion to col-major sparse matrix for distance/diversity calculations.



# ecodive 2.1.0

## Breaking Changes

* `rescale` parameter is superseded by `norm` parameter.


## Fixes

* `alpha_div()` and `beta_div()` will only forward applicable parameters.
* No longer crashes when `pairs = integer(0)`.
* `pseudocount = NULL` correctly selects smallest non-zero count.



# ecodive 2.0.0

## Breaking Changes

* Input matrix is now samples as rows.
* Removed `weighted` parameter.

## Added New Diversity Methods

* Updated Alpha Options: `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`

* Updated Beta Options: `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`



# ecodive 1.0.0

* This is the first release of ecodive.
