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
