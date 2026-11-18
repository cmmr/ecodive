# `rda()`: redundancy analysis, including transformation-based and partial RDA

**Priority:** 3. **Tier:** port the core from scikit-bio (BSD-3), and build
the partial and constrained machinery on ecodive's own `dbrda()` code.
Read [README.md](README.md) first.


## What it is

RDA is constrained PCA: the axes of a (transformed) count table that sample
metadata explain, followed by unconstrained PCA axes of the residuals.

- **tb-RDA** (transformation-based RDA) is RDA after a Hellinger or chord
  transformation (Legendre & Gallagher 2001). In ecodive that's just
  `norm = 'hellinger'` (added in ee945d5) or `norm = 'chord'`, not a
  separate method.
- **Partial RDA** removes `Condition()` covariates first.

It has high demand: vegan `rda`, phyloseq `ordinate(method = "RDA")`,
microeco, mia `runRDA`.


## Design

- **Signature:** `rda(counts, groups, data = NULL, formula = NULL, margin = 1L, norm = 'none', pseudocount = NULL, scale = FALSE, n_axes = 3L, cpus = n_cpus())`.
  - `groups`, `data` and `formula` work as in `dbrda()` (`model_design()`),
    including `Condition()`.
  - With no constraints, it is a PCA (vegan's `rda(X)`): keep that
    working, and its scores should equal `pca()`'s.
  - `scale = TRUE` standardises features, as vegan `rda(scale = TRUE)` and
    `pca(scale)` do.
- **Wrapper:** add `rda` to `COUNT_ORDS` in `R/count_ord.r`, with aliases
  `redundancyanalysis`, `tbrda`, `transformationbasedrda`.
  - Check that `tbrda` doesn't prefix-collide with anything.
  - Make `match_dist_ord('rda')` point to `count_ord()` (that is automatic
    once `rda` is in `COUNT_ORDS`).
- **Result:** an `ecodive_ord` matrix, with columns `RDA1..` and then
  `PC1..` for unconstrained axes, filling up to `n_axes`.
  - Site scores ("wa" or "lc" — decide which and document it; vegan's
    `scores()` default for sites is "wa").
  - Use scaling 2, matching `dbrda()`'s choice.
  - Attributes, in order: `method`, `eig`, `prop`; then `trace`, `inertia`
    (conditional, constrained, unconstrained), `constrained_eig`,
    `r_squared`, `adj_r_squared`, `biplot`, `centroids`, and **`loadings`**
    (species scores, which RDA has and dbrda doesn't).
  - Mirror `dbrda()`'s attribute names exactly where the meaning is the
    same.
- **Name collision:** `vegan::rda` masks or is masked by it. Document this
  in "Same names in other packages".


## Reference implementations

- **Port from:** scikit-bio 0.7.4, `skbio.stats.ordination.rda()`
  (`_redundancy_analysis.py`), BSD-3,
  [LICENSE](https://github.com/scikit-bio/scikit-bio/blob/main/LICENSE.txt).
  - **Caveat (verified in the source):** skbio's "eigenvalues" are
    **singular values `s`, not squared**, so its `proportion_explained`
    (`s / s.sum()`) is not a share of variance. **Don't copy that:**
    ecodive's `eig` must be variances, as vegan's are.
  - skbio has no `Condition()` / partial RDA and no permutation test. It
    centres X only, needs dummy-coded factors, and supports scalings 1 and
    2.
  - Its scaling constant `(Σ s²)^¼` matches vegan's.
- **Partial RDA and the constrained solver:** extend ecodive's own
  `src/dbrda.c` and `R/dbrda.r` (`dbrda_core`, `condition_basis`,
  `design_basis`), which already do `Condition()` with QR, rather than
  writing a second solver.
  - dbrda works on Gower's centred matrix G; RDA works on the centred
    (scaled) table Y, where G = YYᵀ.
  - Generalise the kernel to take either a table or G. Keep the table side
    sparse when `norm` keeps zeros.
  - Use the block Lanczos code (`src/lanczos.c`) for the leading
    constrained and residual axes, as `pca()` does.
- **Output-only fixtures (GPL):** vegan 2.7-x `rda()`, `scores()`,
  `RsquareAdj()` and `eigenvals()`, including `rda(Y ~ X + Condition(Z))`.
- **ecopy `rda` (MIT):** don't use it as a reference. It defaults to
  `scale_y = TRUE`, and its `anova` is wrong (it resamples with
  replacement and uses the wrong df).


## Tests

- Fixtures from vegan, generated with `data-raw/expected/rda.r`, mirroring
  `data-raw/expected/dbrda.r`.
  - Data: dune / dune.env, varespec / varechem, and a larger random table
    (n of about 110).
  - Variants: Hellinger, chord and none; `scale = TRUE`; factor and
    numeric constraints; `Condition()`; constraints with no rank; `n_axes`
    larger than the number of constrained axes.
  - Compare scores (sign-matched), `eig`, `prop`, `trace`, `inertia`,
    `constrained_eig`, R², adjusted R², biplot, centroids and loadings, at
    1e-10.
- With no constraints, scores equal `pca(norm = same)`'s up to sign.
- `rda(norm = 'hellinger')` with `scale = FALSE` has the same constrained
  eigenvalues as `dbrda(hellinger(counts))`. The scales differ by n − 1,
  as the `eig` docs explain, so check and document the relation.
- Sparse and dense input give the same results; results are identical for
  `cpus` 1, 2 and 4.
- Edge cases: a constant feature, a singleton level, n = 3.
- `count_ord(counts, 'rda', 'grp', meta)` equals `rda(counts, 'grp', meta)`.
- A smoke-test line.

**Permutation tests:** vegan's `anova.cca` on an RDA equals a PERMANOVA on
Euclidean distances of the transformed table. Document that
`permanova(euclidean(counts, norm = 'hellinger'), ...)` gives the test, and
check it against the fixture. Don't add a separate test function.


## Docs

- An Rd page, with a "Same names" section and the Hellinger and chord
  tb-RDA example.
- Add to the Ordination section of `pkgdown/_pkgdown.yml`.
- Update the "Constrained ordination" section of `vignettes/ordination.Rmd`.
- NEWS.
- Note in the docs that `varpart()`, `forward_select()` and `principal_response_curves()` will
  build on `rda()`.


## References

- Rao, C. R. (1964). The use and interpretation of principal component
  analysis in applied research. *Sankhyā A*, 26, 329–358.
- Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful
  transformations for ordination of species data. *Oecologia*, 129,
  271–280.
- Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the
  spatial component of ecological variation. *Ecology*, 73, 1045–1055.
- Legendre, P., Oksanen, J., & ter Braak, C. J. F. (2011). Testing the
  significance of canonical axes in redundancy analysis. *Methods in
  Ecology and Evolution*, 2, 269–277.
- Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
  §11.1.


## Decisions made while implementing

Asked of the user at the start (2026-11-18), and answered:

1. **Site scores are "wa" scores** (vegan's `scores()` default for sites,
   and what `dbrda()` / `capscale()` return), in scaling 2. No `lc`
   attribute.
2. **`loadings` are vegan's species scores in scaling 2**
   (`scores(display = "species")`), so that an `rda()` biplot matches
   vegan's. The docs say they differ from `pca()`'s unit-length loadings.
3. **Scaling 2 everywhere.** The plan's "scores equal `pca()`'s up to
   sign" conflicts with scaling 2, so the test checks that unconstrained
   `rda()` has `pca()`'s `eig`, `prop` and `trace`, and that its scores and
   loadings are `pca()`'s rescaled per axis (same directions).
4. **`rda(counts ~ terms, meta)` is accepted,** as `dbrda(dm ~ terms)` and
   `vegan::rda(Y ~ X)` are.

Made while implementing:

5. **No code was ported from scikit-bio.** `dbrda()`'s solver already does
   all of RDA's algebra (a QR basis of the centred design, `Condition()`,
   the small constrained eigenproblem, Lanczos for residual axes), since
   G = YYᵀ. So there is no scikit-bio licence header, `inst/COPYRIGHTS`
   stanza or `Authors@R` change. scikit-bio 0.7.4 is used as an
   output reference instead (`data-raw/expected/rda.py`), in scaling 2.
6. **The kernel was generalised, not duplicated.**
   - `src/dbrda.c`: `bop_t` gained a table mode (`g == NULL`: B = AAᵀ
     through a product callback). `C_dbrda()` was split into
     `cord_design()` (centring and Gram-Schmidt of the design) and
     `cord_solve()` (everything after the operator is set up).
     `C_rda()` calls both.
   - `src/pca.c`: the table setup of `C_pca()` became
     `pca_table_init()`, with `pca_gram_mult()` (AAᵀV) and `pca_cross()`
     (AᵀU), shared by `C_pca()` and `C_rda()`. Sparse tables stay sparse,
     including CLR's dense-zero trick.
   - On the R side, `constrained_ord()` was split into `ord_model()` (the
     model from `groups` / `data` / `formula`), `cord_design()`,
     `cord_empty()` and `cord_result()` (the shared result builder), used
     by `dbrda()`, `capscale()` and `rda()`.
   - `dbrda()`, `capscale()` and `pca()` results are unchanged: their tests
     pass as before.
7. **Unconstrained axes always use the sample-side Lanczos** (as
   `dbrda()` does), even when there are fewer features than samples. Each
   product costs O(nnz), so there was no need for `pca()`'s choice of side.
8. **With no `groups` and no `formula`, `rda()` is an unconstrained PCA.**
   An explicit `NULL` `groups` (a misspelled column) is still an error, as
   in `dbrda()`. `count_ord(counts, 'rda')` therefore runs a PCA.
9. **Scale:** eigenvalues and inertia are variances (G is multiplied by
   1/(n − 1), as in `vegan::rda()`); `dbrda()`'s "mean squared distances"
   rule doesn't apply. So, against `dbrda(euclidean(counts, norm =
   'hellinger'))`, `eig`, `trace`, `inertia` and `constrained_eig` differ
   by n − 1, and scores by (n − 1)^¼. `prop`, R² and biplot scores are the
   same. When a distance reaches 4, `dbrda()` divides by n − 1 too, and
   everything matches. Both cases are documented and tested.
10. **`scale = TRUE` with constant features:** they are left out (loading
    0, and `trace` counts only the other features), as `pca(scale = TRUE)`
    does, where `vegan::rda(scale = TRUE)` stops. The fixture gives vegan
    the table without them.
11. **`n_axes`** is reduced only to the positive eigenvalues (at most
    n − 1). It is not bounded by the number of features: one feature gives
    `RDA1` and `PC1` (the fitted values and the residuals).
12. **Attributes:** `dbrda()`'s, in its order, then `loadings`.
    `centroids` is absent when no constraint is a factor, as in `dbrda()`.
    Columns are `RDA1..`, then `PC1..`.
13. **Fixtures use vegan 2.7.8 and permute 0.10.0** (the plan asks for
    2.7-x). The apt build's vegan 2.6-4 gives the same RDA results, but
    `anova.cca()` permutes partial models differently from 2.7-x.
    `permanova()` matches 2.7.8's F and every permuted F exactly, for full
    and partial models. vegan 2.7.8 needs a `shuffleSet()` matrix rather
    than a plain one.
14. **Registry order:** `rda` sits between `pca` and `rpca` in
    `COUNT_ORDS` (alphabetical). So `match_count_ord('r')` is now
    ambiguous between `rda` and `rpca`; it used to select `rpca`. The
    aliases `redundancyanalysis`, `tbrda` and `transformationbasedrda`
    collide with nothing in the four registries.
15. **Permutation test:** documented as
    `permanova(euclidean(counts, norm = 'hellinger'), ...)`, and tested
    against `anova.cca()` on a tb-RDA, full and partial (same F and
    p-value for the same permutations). No separate test function.
16. **Docs:** NEWS (under 2.4.0, the version in development); README and
    the `ecodive` vignette's method tables; the `ordination` vignette's
    "Constrained ordination" section; `count_ord()`, `dist_ord()` and
    `print.ecodive_ord()` docs; and pkgdown's Ordination section.

Reference results:

- **vegan 2.7.8**, 29 cases (dune, varespec, a 110 × 40 table; Hellinger,
  chord, none, percent, binary, CLR; `scale`; factors, ordered factors,
  numeric constraints and interactions; `Condition()`; conditions only;
  unconstrained; aliased; saturated; constant features; `n_axes` beyond the
  constrained axes, n − 1 and the rank). Sites, species, `eig`, `prop`,
  `trace`, `inertia`, `constrained_eig`, R², adjusted R², biplot and
  centroids all match at the 1e-10 tolerance. The largest differences seen
  were about 2e-15 for constrained results and 7e-14 for the
  Lanczos-computed unconstrained axes.
- **scikit-bio 0.7.4**, scaling 2: all 19 axes' sites, species and biplot
  scores match to about 7e-15, and `eig` = s² / (n − 1) for its "eigvals"
  s.
- **Speed:** 3,000 samples × 2,000 taxa (5% non-zero), `~ g + x +
  Condition(z)`, 8 axes. `rda()` takes 0.46 s on 4 threads and
  `vegan::rda()` 63 s, with eigenvalues equal to 6e-18.
