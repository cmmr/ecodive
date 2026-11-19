# `lda()` and `capda()`: linear discriminant analysis, and discriminant canonical analysis of principal coordinates

**Priority:** 10 (`lda`) and 26 (`capda`). **Tier:**
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

- **`lda(counts, groups, data = NULL, margin = 1L, norm = 'none', pseudocount = NULL, n_axes = NULL, cpus = n_cpus())`.**
  - `n_axes = NULL` means min(3, number of groups − 1), the most
    discriminant axes there are, capped at 3 as for other ordinations
    with eigenvalues (see [parameter_names.md](parameter_names.md) §3.1).
    It can't exceed the number of features either.
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
- **`capda(dm, groups, data = NULL, n_pcoa_axes = NULL, n_axes = NULL, permutations = 999L, seed = 0L, cpus = n_cpus())`.**
  - `n_pcoa_axes` is the number of PCoA axes *used* (Anderson & Willis's
    *m*). The default (`NULL`) chooses it by leave-one-out classification
    success, as they do.
  - `n_axes` is the number of discriminant axes *returned*, as in every
    ordination. `NULL` means min(3, number of groups − 1), as for `lda()`.
  - Don't name it `m`: R would partially match `m` to `dist_ord()`'s own
    `method` argument, so `dist_ord(dm, 'capda', m = 3)` would fail (see
    README, "Argument names").
  - **Result:** an `ecodive_ord` matrix, columns `CAP1..`.
    - Attributes: `eig` (the squared canonical correlations, δ²), `prop`,
      `n_pcoa_axes`, `success` (the leave-one-out classification table and
      overall rate), `trace_stat`, `p_value` (permutation test of the trace
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
- `capda(dm, n_pcoa_axes = M)` equals `lda()` on `pcoa(dm, n_axes = M)`
  scores.
- Edge cases: two groups (one discriminant axis, so `n_axes = NULL`
  gives 1); a singleton group; singular input to
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


## Decisions made while implementing (lda)

Item 10 (`lda()` only; `capda()` is item 26 and still open).

Asked of the user at the start (2026-11-19), and answered:

1. **Scores are scaled as MASS scales them** (the user's choice, had the
   two differed). Checked: they do differ, by a constant. scikit-learn
   1.9.1's `_solve_svd()` divides the within-group scatter by n
   (`fac = 1 / n_samples`), MASS by n − g, so scikit-learn's `transform()`
   scores and `scalings_` are MASS's times √(n / (n − g)) (1.0102 for
   iris), and its between-group singular values² are MASS's `svd`² times
   n / (n − g). `explained_variance_ratio_` equals MASS's proportion of
   trace. The port divides by n − g; everything else follows scikit-learn.
   So the pooled within-group covariance of the scores is the identity.
2. **`eig` is MASS's `svd`²**: the between- over within-group mean square
   of each axis's scores (the "canonical F"), not the SS-ratio eigenvalues
   of W⁻¹B. `prop = eig / sum(eig over all min(g, p) axes)`.
3. **`prior = NULL` argument**: the group proportions by default, else one
   positive value per group, named by group (any order) or in level order;
   it must add up to 1 within 1e-5 (MASS errors, scikit-learn rescales with
   a warning; here it is an error), and is then rescaled to exactly 1.
   Zero priors are rejected (MASS allows them). Priors set the centre
   x̄ = Σ πₖmₖ and weight the between-group term, as in both references.
4. **No public prediction.** No `predict()` method and no new S3 class
   (as `plsda()`). Instead an internal, exact leave-one-out classifier
   (`lda_loo()` in R, `C_lda_loo` in C) for `capda()`, tested against
   `MASS::lda(CV = TRUE)`.

The plan's own call, not re-asked: **a singular within-group covariance is
an error** (scikit-learn and MASS warn and drop null directions).

Made while implementing:

5. **Port scope and licensing.** Ported from scikit-learn 1.9.1
   `sklearn/discriminant_analysis.py` (BSD-3):
   `LinearDiscriminantAnalysis._solve_svd()`, the priors handling in
   `fit()`, and `transform()`. The BSD-3 text is in `src/lda.c`'s header;
   `src/lda.c` and `R/lda.r` were added to the existing scikit-learn stanza
   of `inst/COPYRIGHTS`, and the scikit-learn `Authors@R` comment now
   mentions LDA (no new person). While copying the licence text, the last
   line of the BSD notice ("OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.") was found missing from `src/plsda.c`'s
   header (item 7) and was restored.
6. **Dense, small linear algebra in C, no LAPACK.** LDA needs p ≤ n − g, so
   the normalized table is made dense (n × p ≤ n²) through `C_normalize`
   (every table type, `margin`, every `norm`). The SVDs are a Householder
   QR of the scaled residuals Z (n × p, the reflections applied to the
   columns in parallel), then one-sided (Hestenes) Jacobi on the p × p R;
   and one-sided Jacobi on the g × p matrix of sphered weighted group
   means. Both give singular values to full relative accuracy, without
   forming W or B. Scores are computed in parallel over samples. Columns
   at rounding level (norm² ≤ ε²‖A‖²_F) are skipped by Jacobi (B has rank
   ≤ g − 1). 5,000 × 200 with 5 groups: 0.2 s vs 0.77 s for `MASS::lda()`.
7. **Rank and tolerances.** As scikit-learn: within-group rank = number of
   singular values of Z (each feature scaled by its within-group sd,
   divisor n; Z divided by √(n − g)) above `tol = 1e-4`; if below p, error.
   p > n − g errors before the table is made dense. Between-group axes are
   kept if their singular value is above 1e-4 times the first **and above
   1e-8** (the latter is new: groups with identical means give no axes,
   with a warning, "there are no discriminant axes", and a 0-column
   result, rather than an axis of rounding noise). n ≤ g is an error
   (scikit-learn's rule).
8. **The singular error message** suggests "a few axes of pcoa() or pca()".
   The plan says to suggest `capda()`, which doesn't exist yet;
   **item 26 (formerly 24, and 14 before that) should add `capda()` to `lda_singular()`'s message and to
   the "Singular tables" section of `lda()`'s docs**, and to the
   vignette's LDA paragraph (the plan's "point to `capda()` for microbiome
   tables").
9. **Signs:** each axis is flipped so that its largest |score| is
   positive (`ord_flip()`'s rule), with its column of `scaling`.
10. **`n_axes`:** `NULL` means min(3, number of axes); an explicit value is
    capped silently at the number of discriminant axes (≤ min(g − 1, p)).
11. **Attributes, in order:** `method`, `eig`, `prop`, `scaling` (p × k,
    columns `LD1..`), `means` (g × p), `prior` (named), `groups` (a factor,
    named by sample). `eig`, `prop` and `scaling` cover the returned axes
    only.
12. **Groups:** through `test_groups()`, as `plsda()`: every distinct value
    is a group, unused levels dropped, at least two groups, NAs an error.
    A singleton group is fine for the fit (it adds nothing to W).
13. **Leave-one-out core for `capda()`** (`lda_loo(x, grp, prior, cpus)`,
    in `R/lda.r`; C in `src/lda.c`). For a dense n × p `x`, it returns
    `class`, an n × p integer matrix of level codes whose column m is the
    LOO class using the **first m columns** (for every m = 1..p at once),
    and `posterior` (n × g) for all p columns. It is exact: the fit without
    sample i (group means and pooled covariance, divisor n − 1 − g),
    with the **same priors** (MASS's `CV = TRUE`; checked equal to
    refitting MASS without i with `prior =` the full-data priors, 2e-15).
    It spheres by the R of a QR of the residuals (nested: R⁻ᵀ is lower
    triangular, so the first m sphered coordinates use only the first m
    columns), and downdates by Sherman–Morrison, so every m costs O(g)
    more per sample: 5,000 × 200, all 200 values of m, 0.1 s on 4 threads
    (MASS takes 0.9 s for one m). Classification uses the full-rank
    Mahalanobis rule, which equals classifying on all the discriminant
    axes. NA where a fit is singular: from the first m whose R[m, m] is
    below 1e-4·√(n − g), where x_i alone spans a direction
    (1 − c·u'u ≤ √ε), when n − 1 − g < 1, and **for a sample alone in its
    group** (as MASS, which returns NA there). Ties go to the first group.
    `capda()` must decide how to count NAs in its success rate.
14. **Registry:** `lda` in `COUNT_ORDS` (between `cca` and `pca`), alias
    `lineardiscriminantanalysis`. `'l'` and `'LDA'` select it; no alias
    collides with the other registries. `count_ord(m, 'lda', groups,
    data)` passes `groups`/`data` through `...`; the `count_ord()` example
    skips `lda` in its loop, as `plsda`.
15. **Fixtures:** `data-raw/expected/lda.r` (inputs, MASS 7.3-60.0.1) and
    `lda.py` (scikit-learn, asserts 1.9.1), CSVs in
    `data-raw/expected/lda/`. Cases: iris; iris with priors (.2, .3, .5);
    iris two species (one axis); classical PCoA (cmdscale of Bray-Curtis,
    6 axes) of dune by Management (capda's input); Hellinger of 8 common
    dune species by Moisture; a 151 × 12 sparse Poisson table with 4
    unequal groups including a singleton (also run as Matrix,
    TsparseMatrix, slam and `margin = 2`). MASS: scores, scaling, svd²,
    prop, means, prior, and `CV = TRUE` classes for every m and posteriors.
    scikit-learn: scores and `scalings_` (÷ √(n / (n − g))),
    `explained_variance_ratio_`, `means_`, `xbar_`. Largest relative
    differences: MASS 1.1e-14 (scaling, iris with priors), 4.5e-15
    (scores); scikit-learn 8.5e-15 (scores), 8.4e-16 (prop); LOO
    posteriors 5.4e-15, classes identical.
16. **Docs:** `lda()` joins the "Supervised" subtitle in
    `pkgdown/_pkgdown.yml`; a paragraph and table row in
    `vignette("ordination")`; rows in the README and `vignette("ecodive")`
    tables; `count_ord()`'s docs; NEWS under 2.4.0. "Same names in other
    packages" covers `MASS::lda()`, latent Dirichlet allocation and LEfSe's
    LDA score. The benchmark article was not updated.
