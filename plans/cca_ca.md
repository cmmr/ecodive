# `cca()` and `ca()`: canonical correspondence analysis and correspondence analysis

**Priority:** 4 (`cca`) and 5 (`ca`). **Tier:** port from scikit-bio
(BSD-3). Build after [rda.md](rda.md), whose constrained solver these
reuse with chi-square weights. Read [README.md](README.md) first.

These are **two functions**, following the one-function-per-method rule:
CA is unconstrained, and CCA is constrained by metadata. They share one
weighted core.


## What they are

- **CA** (correspondence analysis) ordinates a table of non-negative
  abundances in the chi-square metric: an SVD of the matrix of
  standardised residuals, with row and column weights.
- **CCA** (canonical correspondence analysis; ter Braak 1986) is the
  constrained form: weighted regression of that matrix on the
  constraints, then SVD, plus unconstrained CA axes of the residuals.
  Partial CCA uses `Condition()`.

Both are high demand in ecology, and medium in microbiome work: phyloseq
`ordinate("CCA")`, microeco, mia `runCCA`, vegan `cca`.

**Name clash to document:** scikit-learn's `CCA`, MultivariateStats.jl's
`CCA` and linfa's `Cca` are canonical **correlation** analysis, a different
method. Always say "canonical correspondence analysis" in the docs.


## Design

- **`ca(counts, margin = 1L, n_axes = 3L, cpus = n_cpus())`.**
  - Raw non-negative abundances; no `norm`, since the chi-square metric is
    the method.
  - Error on negative values, and on all-zero rows or columns: drop them
    with a warning, as vegan does, or error — decide and document.
  - Columns `CA1..`.
  - Attributes: `method`, `eig`, `prop`, `trace` (total inertia), and
    `loadings` (species scores).
- **`cca(counts, groups, data = NULL, formula = NULL, margin = 1L, n_axes = 3L, cpus = n_cpus())`.**
  - `groups`, `data` and `formula` as in `dbrda()` / `rda()`, including
    `Condition()`.
  - Columns `CCA1..`, then `CA1..`.
  - Attributes as in `rda()`: `trace`, `inertia`, `constrained_eig`,
    `r_squared`, `adj_r_squared`, `biplot`, `centroids`, `loadings`.
  - Site scores: decide between "wa" and "lc", document it, and match
    vegan's default. Use scaling 2 (or match vegan's `scores()` default).
  - With no constraints, `cca()` equals `ca()`.
- **Wrapper:** both join `COUNT_ORDS`.
  - Aliases: `cca` takes `canonicalcorrespondenceanalysis` and
    `constrainedcorrespondenceanalysis`; `ca` takes
    `correspondenceanalysis`.
  - **`'CA'` and `'ca'` must select `ca`, not capscale.** The exact-alias
    matching added in ee945d5 handles this; add a test.
  - Avoid adding an alias `cap...` to these.
- **Masking:** `vegan::cca`, `ca::ca`, `ade4::dudi.coa` (different name).


## Reference implementations

- **Port from:** scikit-bio 0.7.4 (BSD-3).
  - `skbio.stats.ordination.ca()` (Legendre & Legendre §9.4, eqs 9.43a and
    9.44a). skbio says it matches vegan `cca(X)`.
  - `skbio.stats.ordination.cca()` (ter Braak 1986; L&L §11.2). skbio says
    it is "expected to give the same results as" vegan `cca`, **but**:
    - it doesn't drop collinear constraints (ecodive's
      `design_basis()` does);
    - it has no partial CCA (build it on ecodive's `condition_basis()`);
    - it supports scalings 1 and 2 only;
    - it has no permutation test.
- **ecopy `ca` and `cca` (MIT):** secondary, and unmaintained. ecopy's
  `anova` permutes raw rows rather than reduced-model residuals, so don't
  copy its test.
- **Core:** the weighted variant of the RDA solver. Apply row weights
  r = row sums / grand total and column weights c = column sums / grand
  total, form the residual matrix (P − rcᵀ) / √(rcᵀ), and do a weighted
  regression of it on the constraints (weighted QR, with weights r). Keep
  this matrix-free where possible, so sparse tables stay sparse.
- **Output-only fixtures (GPL):** vegan `cca()` (CA with `cca(X)`),
  `scores()`, `eigenvals()`, `RsquareAdj()`; ade4 `dudi.coa` as a second
  CA check.


## Tests

- Fixtures from vegan, written by `data-raw/expected/cca.r`.
  - Data: dune, dune.env, varespec, varechem, and a random sparse table.
  - Cover factor and numeric constraints, `Condition()`, collinear
    constraints, and `n_axes` larger than the number of axes.
  - Compare scores (sign-matched), eig, prop, inertia, R² and adjusted R²,
    biplot, centroids and loadings at 1e-10.
- `ca(x)` equals `cca(x)` with no constraints.
- Name matching: `match_count_ord('CA')$id == 'ca'`,
  `match_count_ord('CCA')$id == 'cca'`, and
  `match_dist_ord('CA')` points to `count_ord()`.
- Sparse and dense input agree; results are identical for any `cpus`.
- Edge cases: an empty row or column, a single non-zero cell, n = 3.
- A smoke-test line for each function.


## Docs

- Rd pages for both, linked to each other, saying "canonical
  correspondence".
- Add them to the Ordination section of `pkgdown/_pkgdown.yml`.
- In `vignettes/ordination.Rmd`, add both to the method table, with a note
  on choosing between unimodal methods (CA, CCA) and linear ones (PCA,
  RDA).
- NEWS.


## References

- ter Braak, C. J. F. (1986). Canonical correspondence analysis: a new
  eigenvector technique for multivariate direct gradient analysis.
  *Ecology*, 67, 1167–1179.
- Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
  §9.2 (CA) and §11.2 (CCA).
- Greenacre, M. (2007). *Correspondence Analysis in Practice* (2nd ed.).


## Decisions made while implementing (cca)

Item 4 (`cca()` only; `ca()` is item 5 and is left open here).

Asked of the user at the start (2026-11-18), and answered:

1. **`adj_r_squared` is estimated by permutation, as
   `vegan::RsquareAdj()` does for CCA** (Peres-Neto et al. 2006):
   `1 - (1 - R²) / (1 - mean permuted R²)`. `cca()` takes
   `permutations = 1000L` (vegan's default; a count or a `shuffleSet()`
   matrix) and `seed = 0L`, with PCG streams, so it is reproducible and
   identical for any `cpus`. A permutation moves the samples' rows of the
   table (with their weights) against the fixed constraints; with the same
   matrix, the value equals vegan 2.7.8's to about 2e-16. `NA` for
   `permutations = 0`, for partial CCA (see 9) and for saturated models
   (no residual inertia; vegan gives `numeric(0)` there).
2. **All-zero samples are an error**, naming them (as vegan stops: "all
   row sums must be >0"); the result must keep one row per sample.
   **All-zero features are left out** of the analysis silently, with
   loadings of 0, so the loadings keep one row per feature (vegan drops
   them from its species scores).
3. **Scores as in `rda()`:** site scores are "wa" scores in scaling 2
   (vegan's `scores()` default), `loadings` are vegan's species scores in
   scaling 2, and `cca(counts ~ terms, meta)` is accepted.

Made while implementing:

4. **No code was ported from scikit-bio.** As in item 3, ecodive's own
   constrained solver does all of the algebra once given the weights, so
   there is no licence header, `inst/COPYRIGHTS` stanza or `Authors@R`
   change. scikit-bio 0.7.4's `cca()` (scaling 2) is an output reference
   (`data-raw/expected/cca.py`).
5. **The weighted core.**
   - `src/cca.c` holds the chi-square table: `Q̄ = A − √r √cᵀ` with
     `A = D_r^{-1/2} P D_c^{-1/2}` as sparse as the counts (CSR and CSC
     copies of A's values), so products with `Q̄` and `Q̄Q̄ᵀ` cost
     O(nnz) plus a rank-one correction, parallel over features (`Q̄ᵀX`)
     and samples (`Q̄X`), deterministic for any `cpus`.
     `cca_table_init()`, `cca_gram_mult()` (a `lanczos_op_t`),
     `cca_cross()` (species scores), `cca_sqrt_r()` / `cca_sqrt_c()`,
     and `cca_perm_inertia()` (the permuted constrained inertia, parallel
     over permutations). Total inertia is `sum(A²) − 1`.
   - `src/dbrda.c`: `cord_design()` takes optional weights `w = √r`: each
     design column is centred with weights `r` and multiplied by `√r`
     (so it is orthogonal to `√r`, which takes the place of the ones
     vector). `cord_solve()` then deflates `√r` in the Lanczos run for the
     unconstrained axes (as the first basis column, with `deflate = 0`),
     and divides the site scores by `√r` before fixing signs. Everything
     else (conditional inertia, the small constrained eigenproblem, the
     biplot scores, which become weighted correlations) is unchanged.
     `C_cca()` sits beside `C_rda()`. `dbrda()`, `capscale()`, `rda()`
     and `pca()` results are unchanged (their tests pass as before).
   - R: `cord_result()` gained `scale` (vegan's constant is 1 for CCA),
     `r2` and `weights`; `dbrda_centroids()` gained `weights` (CCA
     centroids are weighted means, weights `r`, of the lc scores `u/√r`).
   - **For `ca()` (item 5):** `ca(counts)` is `C_cca()` with an empty
     design: `cca(counts)` already returns the CA (columns `CA1..`,
     `r_squared` 0, `adj_r_squared` 0), matching `vegan::cca(X)` (tested
     on dune and varespec). `ca()` can call the same code and keep only
     the attributes the plan lists for it.
6. **Unconstrained CCA is CA:** with no `groups` or `formula`, `cca()`
   returns CA axes (`CA1..`), as `vegan::cca(X)` does; an explicit `NULL`
   `groups` is still an error, as in `rda()`.
7. **R² without constraints is 0, and so is the adjusted R²** (as
   `rda()`), where vegan gives `numeric(0)`.
8. **Signs:** each axis is flipped so that its largest site score (after
   dividing by `√r`) is positive, as for the other methods.
9. **Partial CCA's adjusted R² is `NA`.** At the start, the plan's and
   my understanding was that vegan gives `NA` here, and the user approved
   that. vegan 2.7.8 in fact returns a number for partial models, but its
   documentation doesn't say how, and none of the candidates tried
   reproduced it from outputs (the mean permuted R² of `permutest()`'s
   statistics, of refits with the community or either set of variables
   permuted, Peres-Neto's fraction `[a] = adj(X+Z) − adj(Z)`, with either
   denominator). Without reading vegan's (GPL) source, `NA` is kept, and
   the fixtures store vegan's value for reference. `varpart()` (item 12)
   should compute CCA fractions from non-partial models, as Peres-Neto et
   al. (2006) do.
10. **Registry:** `cca` is first in `COUNT_ORDS`, with aliases
    `canonicalcorrespondenceanalysis` and
    `constrainedcorrespondenceanalysis`, which collide with nothing in the
    four registries. `match_dist_ord('cca')` points to `count_ord()`.
    **Until `ca()` exists, `match_count_ord('ca')` (and `'CA'`)
    prefix-matches `cca`** through the `canonical...` alias; item 5's
    exact id `ca` will take precedence and should be tested then.
    `match_dist_ord('ca')` still selects `capscale` by prefix.
11. **No permutation test** was added (not in the plan; CCA's test is not
    a PERMANOVA of any distance). The docs point to `vegan::anova.cca()`.
12. **No `norm` argument**, as the plan says; negative or non-finite values
    are errors. Non-integer abundances are fine.
13. **Fixtures use vegan 2.7.8 and permute 0.10.0**, from a scratch
    library; `RsquareAdj()` there needs a `shuffleSet()` object, not a plain
    matrix. Each case stores its 99-permutation matrix.
14. **Docs:** NEWS (2.4.0), README and the `ecodive` vignette's tables, the
    `ordination` vignette (method table, and a paragraph on unimodal versus
    linear methods with an example), `count_ord()`, `dist_ord()` and
    `print.ecodive_ord()` docs, and pkgdown's Ordination section. The
    shared "Ordination result" section (`R/documentation.r`) wasn't
    touched, as item 3 didn't touch it for `rda()`.

Reference results:

- **vegan 2.7.8**, 24 cases (dune, varespec, a 110 × 40 random sparse
  table; factor, ordered, numeric and interaction constraints;
  `Condition()` with numeric and factor conditions; conditions only;
  unconstrained CA; aliased and aliased partial; saturated; an all-zero
  feature; `n_axes` beyond the constrained axes, n − 1 and the rank).
  Sites, species, `eig`, `prop`, `trace`, `inertia`, `constrained_eig`,
  R², permutation adjusted R² (non-partial), biplot and centroids all
  match at the 1e-10 tolerance. Largest absolute differences over all
  cases: site scores 6e-13 and species scores 1e-13 (Lanczos-computed CA
  axes; constrained axes about 1e-14), centroids 2e-14, biplot 7e-15,
  trace 2e-15, eigenvalues 1e-15, adjusted R² 6e-16.
- **scikit-bio 0.7.4**, scaling 2, all 19 axes: sites, species and
  eigenvalues match at 1e-10. Its biplot scores are unweighted
  correlations of the weighted constraints, where vegan's are weighted
  correlations, so they are not compared.
- **Speed:** 3,000 samples × 2,000 taxa (5% non-zero),
  `~ g + x + Condition(z)`, 8 axes: `cca()` 0.53 s on 4 threads,
  `vegan::cca()` 66 s, eigenvalues equal to 2e-16. Without the condition
  and with 1,000 permutations for the adjusted R²: 0.9 s.
