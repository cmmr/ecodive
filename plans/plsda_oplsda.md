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


## Decisions made while implementing (oplsda)

Item 11 (`oplsda()` only; sparse PLS-DA is [splsda.md](splsda.md)).

Asked of the user at the start (2026-11-19), and answered:

1. **Columns are `p1`, `o1`, `o2`, ...** (ropls's and SIMCA's labels), not
   `OPLS1`/`ORTH1`, so `plot(ord)` is the usual OPLS-DA score plot.
2. **R²X and R²Y as ropls (and SIMCA) report them.** `prop` is each
   component's R²X, |t|²|p|² over the total sum of squares of the scaled
   table (ropls's `modelDF$R2X`); `eig` = |t|²|p|² / (n − 1) and `trace`
   as in `plsda()`; `r2x` = `sum(prop)` (`R2X(cum)`); `r2y` = the R² of
   the fitted values t_p q, which is ropls's `R2Y(cum)` and pyopls's
   `OPLSValidator.r_squared_Y_` (the same in both). pyopls's own R²X
   (`OPLS.score()`, the share of X left after filtering, which divides a
   scaled residual by the unscaled X when `scale = True`) is not reported;
   the docs say it would be `1 - sum(prop[-1])`.
3. **Two groups only.** More is an error: "`groups` has 3 different
   values, but oplsda() takes exactly two; use plsda() for more."
4. **No Q² or cross-validation, and no orthogonal VIP** (the user asked for
   no extras). The predictive VIP stays, since the plan lists VIP and it
   is what OPLS-DA users rank features by.

Made while implementing:

5. **Port scope and licensing.** Ported from pyopls 20.3.post1 (MIT, ©
   2019 Wright State University), `OPLS.fit()` and `_center_scale_xy()` in
   `pyopls/opls.py`; the MIT text is in `src/oplsda.c`'s header, a new
   `Files: src/oplsda.c R/oplsda.r` stanza is in `inst/COPYRIGHTS`, and
   "Wright State University" (ctb, cph) is in `Authors@R`.
6. **The predictive component is Trygg & Wold's:** t_p = X_k w,
   p_p = X_kᵀt_p / |t_p|², q = yᵀt_p / |t_p|² on the filtered table X_k,
   with w = Xᵀy / |Xᵀy| fixed (X_kᵀy = Xᵀy). pyopls only filters; its
   README fits `PLSRegression(1)` with the default `scale = True` to the
   filtered table, which rescales the filtered features. The fixture uses
   `PLSRegression(1, scale = False)` instead, which is Trygg & Wold's (and
   ropls's) predictive component.
7. **Reuses plsda.c's core** (`pls_x_t`, `pls_x_mult()`, `pls_xt_mult()`,
   `pls_x_deflate()`, `pls_weights()` with g = 1) and `pca_table_init()`:
   matrix-free, sparse input never densified, every `norm`, constant
   features dropped under `scale = TRUE` (zero weights, loadings and VIP;
   P in VIP counts every feature, as in `plsda()`). No changes to
   plsda.c. 1,000 samples × 5,000 taxa, 2 orthogonal components: 0.13 s
   dense, 0.06 s from a sparse Matrix, vs 8.2 s for `ropls::opls()`
   (whose time includes its 7-fold cross-validation).
8. **y** is the 0/1 indicator of the second group (factor level order),
   centred and, with `scale = TRUE`, scaled (ddof = 1), as pyopls's
   `_center_scale_xy()` does. Swapping the groups changes only the sign of
   `y_loading`.
9. **Signs:** each component, predictive and orthogonal, is flipped so its
   largest |score| is positive (ecodive's rule), with its weights and
   loadings (and q for p1). ropls and pyopls are sign-matched in tests.
10. **`n_orthogonal`** (default 1) may be 0 (a one-component `plsda()`,
    tested to match). It is capped silently at min(n − 2, p − 1), so that
    a predictive component is always left. Orthogonal fitting stops with a
    warning ("only k of the first m orthogonal components could be fitted")
    when |w_o| before normalizing is at most √ε |p|. No `n_axes` argument:
    `count_ord(m, 'oplsda', g, n_orthogonal = 2)`; the `count_ord()`
    example and test loop treat it like `lda`/`plsda`.
11. **No covariance with the groups is an error** ("the table has no
    covariance with `groups` (for instance, every feature is constant)"),
    when |Xᵀy| ≤ √ε |X|_F |y|, rather than `plsda()`'s 0-column result with
    a warning, because an OPLS-DA without a predictive component has
    nothing to report.
12. **VIP** (predictive only) is √P |p_pj| / |p_p| (Galindo-Prieto et al.
    2014 with one predictive component), a named vector; it matches ropls's
    `vipVn` to 4e-14. (ropls's VIP uses the loadings p, not the weights w.)
13. **Attributes, in order:** `method`, `eig`, `prop`, `trace`, `r2x`,
    `r2y`, `loadings` (p × (1 + k): p_p, then the p_o), `weights` (w, then
    the w_o), `y_loading` (q, a number), `vip`, `groups` (two-level factor,
    named by sample), `centroids` (2 × (1 + k)).
14. **Registry:** `oplsda` sits between `lda` and `pca` in `COUNT_ORDS`,
    with aliases `opls`, `orthogonalplsda`,
    `orthogonalpartialleastsquaresda` and
    `orthogonalpartialleastsquaresdiscriminantanalysis`; `'pls'` still
    selects only `plsda`, and `'p'` stays ambiguous (`pca`, `plsda`). No
    alias collides with the other three registries.
15. **Fixtures:** `data-raw/expected/oplsda.r` (inputs, ropls, folding)
    and `oplsda.py` (pyopls 20.3.post1 + scikit-learn 1.9.1, both
    asserted), CSVs in `data-raw/expected/oplsda/`. ropls is **1.34.0**
    (Bioconductor 3.18, the newest for the container's R 4.3.3), not the
    plan's 1.44. ropls rounds R2X/R2Y to 3 digits, so the fixture also
    derives unrounded ones from its outputs (|t|²|p|² over its scaled
    table `xModelMN`; R² of its fitted values `yPreMN`). Cases: dune by
    management (natural vs managed) with 1 and 2 orthogonal components,
    the same with constant and all-zero features added, Hellinger with
    `scale = FALSE` (hayfield vs other), CLR by moisture (wet vs dry) with
    3, and a 110 × 40 sparse table with 4 (also run as Matrix,
    TsparseMatrix, slam and `margin = 2`). Largest relative differences:
    pyopls 6.8e-14 (orthogonal weights, CLR case), ropls 1.5e-13
    (weights, big case); R²X and R²Y 6.7e-15.
16. **Docs:** in the "Supervised" pkgdown group after `plsda`; a paragraph
    and example in `vignette("ordination")`'s "Supervised ordination"
    section; rows in the README and `vignette("ecodive")` tables;
    `count_ord()` docs; NEWS under 2.4.0. `oplsda()` docs say Q² is not
    computed and point to `permanova()` or cross-validation.
