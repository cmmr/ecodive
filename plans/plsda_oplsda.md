# `plsda()` and `oplsda()`: partial least squares discriminant analysis, and its orthogonal form

**Priority:** 7 (`plsda`) and 11 (`oplsda`). **Tier:** port from
scikit-learn (BSD-3) and pyopls (MIT). The sparse variant is a separate,
clean-room plan: [splsda.md](splsda.md). Read [README.md](README.md) first.

Following the one-function-per-method rule, PLS-DA, OPLS-DA and sparse
PLS-DA are **three functions**.


## What they are

- **PLS-DA** is a supervised ordination. It finds the components of a
  (usually CLR-transformed) feature table that best separate known groups:
  PLS regression on a dummy-coded group matrix, solved by NIPALS. Common in
  microbiome and metabolomics work through mixOmics `plsda()`.
- **OPLS-DA** (Trygg & Wold 2002) splits the variation into one
  predictive component and orthogonal ones. It predicts the same as PLS-DA
  with one predictive component; it's mainly a rotation that's easier to
  interpret. Popular in metabolomics (ropls).


## Design

- **`plsda(counts, groups, data = NULL, margin = 1L, norm = 'none', pseudocount = NULL, scale = TRUE, n_axes = 2L, cpus = n_cpus())`.**
  - `groups` is a categorical vector or a column name of `data`, handled
    by `test_groups()`-style code.
  - `scale = TRUE` standardises features, as scikit-learn and mixOmics do
    by default. Document that `norm = 'clr'` is the usual choice for
    microbiome counts.
  - **Result:** an `ecodive_ord` matrix of sample scores, columns
    `PLSDA1..`.
    - Attributes: `method`, `eig` (`NA`; or the variance explained in X —
      decide and document), `prop` (the proportion of X variance per
      component).
    - Then `loadings` (X loadings), `weights`, `y_loadings`, `vip`
      (variable importance in projection), `groups`, and `centroids`
      (group means of the scores).
  - **Prediction:** probably none in v1 (ecodive has no `predict()`
    methods). If it's added, mirror mixOmics' `max.dist`,
    `centroids.dist` and `mahalanobis.dist` rules, as a separate
    `predict.ecodive_plsda()`.
- **`oplsda(counts, groups, data = NULL, margin = 1L, norm = 'none', pseudocount = NULL, scale = TRUE, n_orthogonal = 1L, cpus = n_cpus())`.**
  - Columns `p1` and `o1..`, or `OPLS1` and `ORTH1..` (pick a style that
    is consistent with ecodive's others).
  - Attributes: R2X, R2Y, the orthogonal loadings, and VIP.
  - **Two groups only** at first, as pyopls and most OPLS-DA tools have
    (a single-column Y). Error for more, with a hint to use `plsda()`.
- **Wrapper:** both join `COUNT_ORDS` (they ordinate a count table), with
  aliases `plsda`, `partialleastsquaresda`,
  `partialleastsquaresdiscriminantanalysis`, `oplsda`, `opls`,
  `orthogonalplsda`. `count_ord()` passes `groups` and `data` through
  `...`, as `dist_ord()` does for `dbrda`.
  - **Check:** `oplsda` against `opls`, and `plsda` against `pls`; avoid a
    bare `pls` alias, which would match both.
- **Masking:** `mixOmics::plsda`, `caret::plsda`, `ropls::opls`
  (different name).


## Reference implementations

- **`plsda`, port from:** scikit-learn 1.9.1 `PLSRegression` (NIPALS),
  BSD-3, [COPYING](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING),
  in `sklearn/cross_decomposition/_pls.py`.
  - Run it on a dummy-coded Y with `scale = True`.
  - scikit-learn has no classification rule and doesn't compute VIP. Write
    VIP from Wold et al. (2001) or Mehmood et al. (2012): it's a formula.
- **Alternatives:**
  - linfa-pls (Rust, MIT/Apache, 0.8.1);
  - hoggorm `nipalsPLS2` (BSD-2, last updated 2021).
- **`oplsda`, port from:** pyopls (MIT,
  [LICENSE](https://raw.githubusercontent.com/BiRG/pyopls/master/LICENSE),
  © 2019 Wright State University). Unmaintained since 2020; single-column
  Y only. Its R2X is defined differently from libPLS: document which
  definition ecodive uses.
  - Alternative: pypls (Apache-2.0, Cython NIPALS, binary classes only),
    which must keep its Apache notices if ported.
- **Output-only fixtures:**
  - mixOmics 6.36 `plsda()` (GPL): scores, loadings, explained variance;
  - ropls 1.44 `opls(orthoI = 1)` (CeCILL, copyleft): fixtures only.

  Expect sign flips per component, and check how each scales X before
  comparing.
- **Python fixture:** scikit-learn `PLSRegression` gives the exact NIPALS
  path. Add `data-raw/expected/plsda.py`, which asserts scikit-learn's
  version, as `nmds.py` does.


## Tests

- `plsda()` scores, loadings and weights match scikit-learn's
  `PLSRegression(n_components = n_axes, scale = True)` on the dummy matrix at
  1e-10, sign-matched; the explained variance matches mixOmics.
- VIP matches a direct implementation of the formula.
- `oplsda()` matches pyopls (Python fixture) at 1e-10, and ropls's R2Y on
  a shared dataset.
- Two groups and many groups (`plsda` only); unequal group sizes; a
  constant feature (with `scale = TRUE`: error or drop? decide); CLR via
  `norm`; sparse input.
- Results are identical for any `cpus`.
- `count_ord(counts, 'plsda', 'grp', meta)` equals `plsda(...)`.
- A smoke-test line for each function.


## Docs

- Rd pages; add them to the Ordination section of `pkgdown/_pkgdown.yml`,
  in a "Supervised" group with `lda()` and `capda()`.
- Vignette: a warning that supervised ordinations always separate groups,
  so validate them with a permutation test (`permanova()`) or
  cross-validation.
- NEWS.


## References

- Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic
  tool of chemometrics. *Chemometrics and Intelligent Laboratory Systems*,
  58, 109–130.
- Barker, M., & Rayens, W. (2003). Partial least squares for
  discrimination. *Journal of Chemometrics*, 17, 166–173.
- Trygg, J., & Wold, S. (2002). Orthogonal projections to latent
  structures (O-PLS). *Journal of Chemometrics*, 16, 119–128.
- Lê Cao, K.-A., Boitard, S., & Besse, P. (2011). Sparse PLS discriminant
  analysis. *BMC Bioinformatics*, 12, 253. (See [splsda.md](splsda.md).)
