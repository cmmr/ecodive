# `lda()` and `capda()`: linear discriminant analysis, and discriminant canonical analysis of principal coordinates

**Priority:** 10 (`lda`) and 14 (`capda`). **Tier:**
- `lda()` is ported from scikit-learn (BSD-3).
- `capda()` is *derived*: `pcoa()` followed by `lda()`, with the
  number of axes chosen as Anderson & Willis (2003) describe.

Read [README.md](README.md) first.


## What they are

- **LDA** finds the axes that best separate known groups, relative to the
  variation within them. For microbiome tables it is singular when taxa
  outnumber samples, so it is usually run on a few PCoA or PCA axes (which
  is what CAP does).
- **Discriminant CAP** (Anderson & Willis 2003; "CAP" in PRIMER-PERMANOVA+)
  runs LDA on the first m PCoA axes of any distance matrix. It chooses m by
  leave-one-out classification success, and reports that success rate and
  a permutation test.

**Name:** vegan's `capscale()` ("constrained analysis of principal
coordinates") is db-RDA, not Anderson's CAP. ecodive keeps the `cap` alias
for `capscale()`, as already documented, so this is **`capda()`**.
**Don't give it the alias `cap`.** Give it the aliases `capda`,
`capdiscriminant` and `discriminantcap`, and check them against
`capscale`'s aliases.


## Design

- **`lda(counts, groups, data = NULL, margin = 1L, norm = 'none', pseudocount = NULL, k = NULL, cpus = n_cpus())`.**
  - `k` defaults to min(number of groups − 1, number of features).
  - **Result:** an `ecodive_ord` matrix of discriminant scores, columns
    `LD1..`.
    - `eig` holds the discriminant eigenvalues (between-group / within-group
      variance) and `prop` the proportion of trace (as MASS's
      "Proportion of trace").
    - Then `scaling` (the coefficients), `means` (group means), `prior`,
      `groups`.
  - Error when the within-group covariance is singular, suggesting
    `capda()` or fewer features.
  - Joins `COUNT_ORDS`, with aliases `lda`, `lineardiscriminantanalysis`.
- **`capda(dm, groups, data = NULL, n_axes = NULL, permutations = 999L, seed = 0L, cpus = n_cpus())`.**
  - `n_axes` is the number of PCoA axes (Anderson & Willis's *m*). The
    default (`NULL`) chooses it by leave-one-out classification success,
    as they do.
  - Don't name it `m`: R would partially match `m` to `dist_ord()`'s own
    `method` argument, so `dist_ord(dm, 'capda', m = 3)` would fail (see
    README, "Argument names").
  - **Result:** an `ecodive_ord` matrix, columns `CAP1..`.
    - Attributes: `eig` (the squared canonical correlations, δ²), `prop`,
      `n_axes`, `success` (the leave-one-out classification table and overall
      rate), `trace_stat`, `p_value` (permutation test of the trace
      statistic), `n_perm`.
  - Joins `ORDINATIONS` (it takes distances), passing `groups` and `data`
    through `dist_ord()`'s `...`, like `dbrda`.
- **Masking:** `MASS::lda` (MASS is a recommended package and often
  attached). BiodiversityR's function is `CAPdiscrim`, a different name.


## Reference implementations

- **`lda`, port from:** scikit-learn 1.9.1 `LinearDiscriminantAnalysis`
  (`sklearn/discriminant_analysis.py`), BSD-3, using the `'svd'` solver.
  - **Caveat:** scikit-learn's scaling of the discriminant scores differs
    from MASS `lda()`'s `scaling`, so compare up to scale and sign. Decide
    which scaling ecodive reports: MASS's is the R users' expectation. Not
    yet verified.
  - Alternative: MultivariateStats.jl `MulticlassLDA` (MIT, v0.10.5).
- **`capda`:** no permissive implementation. Build it from
  `pcoa()` and `lda()`, following Anderson & Willis (2003) and Anderson,
  Gorley & Clarke (2008, the PERMANOVA+ manual).
- **Output-only fixtures:**
  - MASS `lda()` (GPL-2|3): coefficients up to scale, scores, proportion
    of trace;
  - BiodiversityR 2.18 `CAPdiscrim()` (GPL-3): m, δ², classification
    success, p-value with a matched permutation scheme;
  - optionally the Fathom MATLAB toolbox `f_cap` (GPL-2).
- **Python fixture:** scikit-learn LDA (`data-raw/expected/lda.py`) for the
  exact SVD-solver path.


## Tests

- `lda()` scores match scikit-learn's at 1e-10 (sign- and scale-matched),
  and MASS's proportion of trace and group separation.
- `capda()` matches BiodiversityR's m, δ² and success table on dune
  or varespec, and its p-value with shared permutations.
- `capda(dm, n_axes = M)` equals `lda()` on `pcoa(dm, k = M)` scores.
- Edge cases: two groups (k = 1); a singleton group; singular input to
  `lda()` (a clear error).
- Results are identical for any `cpus`; a smoke-test line for each
  function.


## Docs

- Rd pages; put them in a "Supervised" group of the reference index with
  `plsda()`.
- Vignette: point to `capda()` for microbiome tables, and note that
  `capscale()` is a different method.
- NEWS.


## References

- Fisher, R. A. (1936). The use of multiple measurements in taxonomic
  problems. *Annals of Eugenics*, 7, 179–188.
- Anderson, M. J., & Willis, T. J. (2003). Canonical analysis of principal
  coordinates: a useful method of constrained ordination for ecology.
  *Ecology*, 84, 511–525.
- Anderson, M. J., Gorley, R. N., & Clarke, K. R. (2008). *PERMANOVA+ for
  PRIMER: Guide to Software and Statistical Methods*. PRIMER-E.
