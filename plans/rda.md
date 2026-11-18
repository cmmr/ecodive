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

- **Signature:** `rda(counts, groups, data = NULL, formula = NULL, margin = 1L, norm = 'none', pseudocount = NULL, scale = FALSE, k = 3L, cpus = n_cpus())`.
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
  `PC1..` for unconstrained axes, filling up to `k`.
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
    numeric constraints; `Condition()`; constraints with no rank; k larger
    than the number of constrained axes.
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
