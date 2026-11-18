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


## Decisions made while implementing (plsda)

Item 7 (`plsda()` only; `oplsda()` is item 11 and still open).

Asked of the user at the start (2026-11-18), and answered:

1. **Weights are computed exactly, not by iterating NIPALS.** Each
   component's weights are the top left singular vector of the small
   p × g matrix M = X_aᵀY_a, found from the top eigenvector c of the
   g × g matrix MᵀM (`lz_sym_eigen()`), as w = Mc / |Mc|. This is the
   point NIPALS converges to, and scikit-learn's own `algorithm = "svd"`
   path (`_get_first_singular_vectors_svd()`). scikit-learn's
   `PLSRegression` stops NIPALS when the squared change in w is below
   `tol = 1e-6`, so with three or more groups its default output differs
   from the converged one by about 0.1% (0.09–0.19% on the fixture data);
   with two groups (rank-1 centred Y) NIPALS converges at once. There are
   no `tol` or `max_iters` arguments. Fixtures use scikit-learn with
   `tol = 1e-28` (the converged path) and mixOmics with `tol = 1e-15`;
   the two-group case is also checked against scikit-learn's defaults.
2. **`eig` is the variance of X each component explains,**
   |t|²|p|² / (n − 1), so that `prop = eig / trace` (a `trace` attribute
   was added, as in `pca()` and `rda()`). This equals mixOmics'
   `prop_expl_var$X` times the trace. Because the scores are orthogonal,
   it's the same whether X or the deflated X_a is used, and the values
   add up over components. They need not decrease.
3. **No prediction in v1.** No `predict()` method and no new S3 class.
4. **Constant features under `scale = TRUE` are dropped silently:**
   weights, loadings and VIP of 0, and `trace` counts only the others.
   This is what scikit-learn's output amounts to (it sets a zero std to 1,
   leaving a zero column), and what `pca()` and `rda()` do (through
   `pca_table_init()`). mixOmics was not run on that case.

Made while implementing:

5. **Port scope and licensing.** Ported from scikit-learn 1.9.1
   `sklearn/cross_decomposition/_pls.py` (BSD-3): `_PLS.fit()` (regression
   deflation, scores, loadings, the near-zero-column mask on Y),
   `_center_scale_xy()` (ddof = 1, X and Y both scaled by `scale`),
   `_get_first_singular_vectors_svd()`. The BSD-3 text is in
   `src/plsda.c`'s header, `src/plsda.c` and `R/plsda.r` were added to the
   existing scikit-learn stanza of `inst/COPYRIGHTS`, and the existing
   scikit-learn `Authors@R` entry's comment now mentions PLS-DA (no new
   person).
6. **Matrix-free, sparse.** X is never centred, scaled or deflated in
   memory. `pca_table_init()` sets up A = JZD⁻¹ (any table type, every
   `norm`, CLR's dense-zero trick); a new `pca_mult()` in `src/pca.c`
   gives AV next to the existing `pca_cross()` (AᵀU). Deflation is
   implicit: X_a v = Av − T(Pᵀv), X_aᵀu = Aᵀu − P(Tᵀu). Y (n × g) is
   dense and deflated in place. 1,000 samples × 5,000 taxa, 3 components:
   0.06 s (0.02 s from a sparse Matrix) vs 1.2 s for `mixOmics::plsda()`.
7. **A reusable core for `oplsda()`** (declared in `src/ecodive.h`):
   `pls_x_t` (the deflated table: T, P, k), `pls_x_mult()` (X_a v),
   `pls_xt_mult()` (X_aᵀU), `pls_x_deflate()` (append a (t, p) pair) and
   `pls_weights()` (w and the top singular value of X_aᵀY; for one Y
   column, w = Xᵀy / |Xᵀy|, which is NIPALS' single step). OPLS's
   orthogonal components need exactly these: w_o from p and w, t_o = X_a
   w_o, p_o = X_aᵀt_o / |t_o|², then deflate by (t_o, p_o).
8. **Signs:** each component is flipped so that its largest |score| is
   positive (`lz_sign_index()`, ecodive's rule for every ordination), with
   w, p and q flipped together; scikit-learn instead makes w's largest
   element positive. Tests sign-match.
9. **Early stop:** fitting stops, with a warning ("only k of the first m
   components could be fitted"), when M's top singular value is at most
   √ε ≈ 1.5e-8 times the first component's (e.g., once the components
   span X, or every feature is constant: 0 components). `n_axes` is
   first capped silently at min(n − 1, p), as `pca()` does.
10. **VIP:** the textbook formula (Wold et al. 2001; Mehmood et al. 2012)
    with SSY_b = |t_b|²|q_b|² of Y as analysed, returned as a p × k
    matrix whose column a uses the first a components (mixOmics'
    `vip()` layout). P counts every feature, constant ones included.
    mixOmics' `vip()` weights components by Σ cor(y_j, t_b)², which is
    the same when Y is scaled (the default); with `scale = FALSE` they
    differ (1.3% on the Hellinger fixture). The test checks both: ours
    against the formula and against mixOmics for `scale = TRUE`, and
    mixOmics' correlation form against its fixture for `scale = FALSE`.
11. **Attributes, in order:** `method`, `eig`, `prop`, `trace`,
    `loadings` (P, scikit-learn's `x_loadings_`), `weights` (W,
    `x_weights_`, mixOmics' `loadings$X`), `y_loadings` (Q, g × k, for the
    centred/scaled dummy matrix), `vip`, `groups` (a factor, named by
    sample), `centroids` (g × k group means of the scores). Columns are
    `PLSDA1..`. No `y_scores`, rotations or `x_rotations` (no prediction).
12. **Groups:** through `test_groups()`; every distinct value is a group
    (numbers too), factor levels keep their order and unused ones are
    dropped, dummy columns follow the levels; at least two groups are
    required; NAs are an error (as for every test). A singleton group is
    fine. `scale` scales the dummy columns too, as scikit-learn and
    mixOmics do.
13. **Registry:** `plsda` sits between `pca` and `rda` in `COUNT_ORDS`,
    with aliases `partialleastsquaresda`,
    `partialleastsquaresdiscriminantanalysis` and
    `plsdiscriminantanalysis`; no bare `pls`. `'pls'` and `'PLS-DA'`
    reach `plsda` by prefix for now; once `oplsda` exists, `'pls'` still
    selects only `plsda` (prefixes are matched from the start). The
    single letter `'p'` is now ambiguous (`pca`, `plsda`); it used to
    select `pca`. No alias collides with the other three registries.
    `count_ord(m, 'plsda', groups, data)` passes `groups`/`data`
    through `...`; the `count_ord()` example skips `plsda` in its
    "same arguments for every method" loop.
14. **Fixtures:** `data-raw/expected/plsda.r` (inputs, mixOmics, folding)
    and `plsda.py` (scikit-learn, asserts 1.9.1), CSVs in
    `data-raw/expected/plsda/`. mixOmics is **6.26.0** (Bioconductor
    3.18, the newest for the container's R 4.3.3), not the plan's 6.36.
    Cases: dune by Management (4 unequal groups), two groups, constant and
    all-zero features added, Hellinger with `scale = FALSE`, CLR by
    Moisture with 4 components, and a 110 × 40 sparse table with 3 groups
    and 5 components (also run as Matrix, TsparseMatrix, slam and
    `margin = 2`). Largest relative differences: scikit-learn 2.2e-14,
    mixOmics 3.0e-15 (scores, weights, prop, VIP).
15. **Docs:** a "Supervised" subtitle under Ordination in
    `pkgdown/_pkgdown.yml` (for `lda()`, `capda()` later); a "Supervised
    ordination" section in `vignette("ordination")` with the warning to
    validate with `permanova()` or cross-validation; rows in the README
    and `vignette("ecodive")` tables; NEWS under 2.4.0. The shared
    "Ordination result" section's description of `eig` was not changed
    (its method list was already incomplete for `rda()`/`cca()`/`ca()`);
    `plsda()`'s own docs define its `eig`. The benchmark article's
    reference table was not updated (as for `rda()`/`cca()`).
