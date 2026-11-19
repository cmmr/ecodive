# `coinertia()`: co-inertia analysis, with the RV coefficient and its test

**Priority:** 38. **Tier:** clean room for the weighted analysis; port the
RV coefficient from hoggorm (BSD-2). ade4 (GPL) is output-only. Read
[README.md](README.md) first.


## What it is

Co-inertia analysis (Dolédec & Chessel 1994; Dray, Chessel & Thioulouse
2003) finds paired axes in two tables of the same samples: e.g. taxa and
metabolites, or community and environment. The axes maximise the
covariance between the two sides' scores. The RV coefficient (Robert &
Escoufier 1976) measures overall co-structure, and the RV test permutes
rows to test it. It's used for multi-omics integration.


## Design

- **Signature:** `coinertia(x, y, x_method = 'pca', y_method = 'pca', n_axes = 3L, permutations = 999L, seed = 0L, cpus = n_cpus())`.
  - `x` and `y` are two tables with the same samples (matched by name),
    or ecodive ordination results.
  - `*_method` is the analysis applied to each side: `'pca'` (centred),
    `'pca_scaled'`, `'ca'` (once `ca()` exists), or `'pcoa'` (for a
    `dist`).
  - Each method brings its own row and column weights (ade4's "dudi"
    weights).
- **Result:** a list of class `ecodive_coinertia`, with
  - `x_scores` and `y_scores` (`ecodive_ord` matrices, columns `CoIA1..`);
  - `eig`, `prop`;
  - `rv` (the RV coefficient);
  - `test`, an `ecodive_test` with `statistic = c(RV = ...)`, `p_value`,
    `z`, `n_perm` and attribute `rv_perm`;
  - the `x_loadings` and `y_loadings`.

  It's a plain function, since it takes two inputs.
- **Masking:** `ade4::coinertia`.


## Reference implementations

- **RV coefficient, port from:** hoggorm `RVcoeff` and `RV2coeff`, BSD-2
  ([LICENSE](https://raw.githubusercontent.com/olivertomic/hoggorm/master/LICENSE),
  © 2016 Oliver Tomic). RV = tr(XXᵀYYᵀ) / √(tr((XXᵀ)²) tr((YYᵀ)²)).
- **Uniform weights:** co-inertia of two centred PCAs is the SVD of XᵀY,
  i.e. Tucker's inter-battery analysis. That's scikit-learn's `PLSSVD`
  (BSD-3) or linfa's `PlsSvd`, so port that case from scikit-learn.
- **Weighted (dudi) cases:** clean room, from Dray et al. (2003). Take the
  weighted SVD of D_x^{1/2} Xᵀ D_r Y D_y^{1/2}, with row weights D_r and
  column weights D_x and D_y from each analysis.
- **Output-only fixtures (GPL):**
  - ade4 1.7-24 `coinertia(dudi.pca(X), dudi.pca(Y))`, with `$eig`,
    `$RV`, `$li`, `$co`, `$l1`, `$c1`;
  - `RV.rtest()` with a supplied permutation scheme if possible
    (otherwise compare statistically);
  - FactoMineR `coeffRV()` as an RV cross-check.


## Tests

- RV matches hoggorm (a Python fixture) and ade4 at 1e-12.
- Co-inertia eigenvalues and scores match ade4 at 1e-10 (sign-matched),
  for PCA–PCA, PCA–CA (once `ca()` exists) and PCoA–PCA.
- The uniform-weight case equals scikit-learn's `PLSSVD`.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

An Rd page in a "Comparing ordinations" group with `procrustes()`, and
NEWS.


## References

- Dolédec, S., & Chessel, D. (1994). Co-inertia analysis: an alternative
  method for studying species–environment relationships. *Freshwater
  Biology*, 31, 277–294.
- Dray, S., Chessel, D., & Thioulouse, J. (2003). Co-inertia analysis and
  the linking of ecological data tables. *Ecology*, 84, 3078–3089.
- Robert, P., & Escoufier, Y. (1976). A unifying tool for linear
  multivariate statistical methods: the RV-coefficient. *Applied
  Statistics*, 25, 257–265.
