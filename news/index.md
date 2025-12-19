# Changelog

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
