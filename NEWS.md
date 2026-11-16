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
  matrix whose attributes start with `method`, `eig`, and `prop`.

* New `nmds()` runs a non-metric multidimensional scaling, following
  `vegan::metaMDS()` with its `monoMDS()` engine: a start from the PCoA
  solution, then random starts until the best solution is found again, then
  centering, rotation to principal components, and scaling. Its C port of
  `monoMDS` gives the same result as `vegan::monoMDS()` from the same starting
  configuration. Random starts run in parallel, one per thread, and share one
  sorted copy of the distances, so memory use is a fraction of `metaMDS()`'s.
  Random starts come from a `seed` argument, and results are identical for
  any number of `cpus`.

* New `blas_threads()` gets or sets how many threads R's BLAS/LAPACK library
  uses, while R is running. It works with OpenBLAS, Intel MKL, BLIS, and
  FlexiBLAS, returns `NA` where the count can't be read (e.g. R's reference
  BLAS or Apple's Accelerate), and never raises an error, warning, or
  message. Setting is limited to what is safe: the count is never raised
  above the library's own, and OpenMP's thread count is left as it is. No
  other ecodive function calls it.



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
