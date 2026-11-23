# `ancombc()` and `ancombc2()`: ANCOM with bias correction

**Priority:** 22 (`ancombc()`, `ancombc2()`) and 23 (the post-hoc tests
and `structural_zeros()`). **Tier:** port (scikit-bio 0.7.4,
`skbio/stats/composition/_ancombc.py`, BSD-3, itself derived from the
Artistic-2.0 R package ANCOMBC; both notices apply, see "Licensing"). Read [README.md](README.md) first, and
[compositional.md](compositional.md) for the conventions the compositional
plans share.


## Licensing (decided by the user, 2026-11-19)

scikit-bio's file header says its implementation "is based on an analysis
of the source code of the R package ANCOMBC", which is **Artistic-2.0**,
and thanks ANCOMBC's author for advice; scikit-bio releases it as BSD-3.
The user has allowed Artistic-2.0 code in ecodive (README.md, "Licensing"),
so the port goes ahead with **both notices**, and it complies whichever
licence is held to apply:

- Every ported file (`src/ancombc.c`, `R/ancombc.r`, ...) carries
  scikit-bio's BSD-3 notice ("Copyright (c) 2013--, scikit-bio development
  team.") **and** ANCOMBC's Artistic-2.0 notice. Take ANCOMBC's copyright
  line from its DESCRIPTION / LICENSE (Bioconductor; Huang Lin is the
  author and maintainer), and include the Artistic-2.0 text.
- Each file header says how it differs from the original (Artistic-2.0
  section 4), e.g. "C/R port of ANCOM-BC from scikit-bio 0.7.4
  (`_ancombc.py`), itself derived from the R package ANCOMBC
  <version>; restructured for ecodive." User-visible differences go in the
  Rd's "Differences from ANCOMBC / scikit-bio" section.
- `inst/COPYRIGHTS`: one `Files:` stanza for these files with both
  licences. `Authors@R`: scikit-bio is already credited; add Huang Lin as
  `ctb`/`cph` ("ANCOM-BC, via scikit-bio's port").
- Because ANCOMBC is now portable under README.md's rule, its R source may
  be **read** to resolve differences between scikit-bio and ANCOMBC
  outputs (not only run for fixtures). Port from scikit-bio first: it's
  restructured for batching and already matched to ANCOMBC.
- Optional courtesy (the user's call): tell Huang Lin about the port; a
  permission to use MIT would let the Artistic notice go.


## What it is

ANCOM-BC (Lin & Peddada 2020) models each feature's log abundance with a
linear model, and estimates and removes the **sample-specific sampling
fraction** (the unknown ratio between a sample's sequenced and true
absolute abundance), which biases naive log fold changes. It gives each
feature × covariate a log fold change, standard error, test statistic and
p-value, with covariates via a formula. ANCOM-BC2 (Lin & Peddada 2024)
adds explicit sampling-fraction correction before refitting, variance
regularisation for rare features, and zeros handled as missing, plus
multigroup post-hoc tests. Both are among the most used
differential-abundance methods in current microbiome papers, and
Nearing et al. (2022) found ANCOM-family methods among the most
consistent.

They're separate named methods (README rule): `ancombc()` and
`ancombc2()`. scikit-bio's post-hoc analyses (global, pairwise, Dunnett,
trend) and structural-zero detection are separate named tests too:

| Function | scikit-bio | What |
|:---------|:-----------|:-----|
| `ancombc()` | `ancombc()` | ANCOM-BC primary results |
| `ancombc2()` | `ancombc2()` | ANCOM-BC2 primary results |
| `ancombc_global()` | `ANCOMBCResult.global_test()` | is a multi-level factor associated with the feature at all? |
| `ancombc_pairwise()` | `.pairwise_test()` | all pairs of levels, mdFDR-controlled |
| `ancombc_dunnett()` | `.dunnett_test()` | each level vs the reference |
| `ancombc_trend()` | `.trend_test()` | monotone / ordered trend across levels |
| `structural_zeros()` | `struc_zero()` | features systematically absent from a group |

Names: the post-hoc tests take an `ancombc()` **or** `ancombc2()` result,
so they share the `ancombc_` prefix (rule 2). `structural_zeros()` is
spelled out (rule 4: `struc_zero` reads oddly, and the test isn't
ANCOM-specific). **No wrapper** (per-feature tables; see compositional.md).


## Design

### Signatures

```r
ancombc (counts, formula, data, groups = NULL, margin = 1L,
         pseudocount = NULL, max_iters = 100L, tol = 1e-5, alpha = 0.05,
         p_adjust = 'holm', cpus = n_cpus())

ancombc2(counts, formula, data, groups = NULL, margin = 1L,
         pseudocount = 0, var_quantile = 0.05, max_iters = 100L,
         tol = 1e-5, alpha = 0.05, p_adjust = 'holm', cpus = n_cpus())

ancombc_global  (fit, ...)
ancombc_pairwise(fit, ...)
ancombc_dunnett (fit, ...)
ancombc_trend   (fit, ...)
structural_zeros(counts, groups, data = NULL, margin = 1L, neg_lb = FALSE)
```

- Post-hoc arguments (`...` above) are those of scikit-bio's
  `ANCOMBCResult` methods, renamed to the glossary (`alpha`, `p_adjust`
  default to the fit's; bootstrap counts as `n_<unit>`; `seed = 0L` where
  they draw at random). Read their signatures from the 0.7.4 source.
- `formula`, `data`: a one-sided formula (`~ status + age`) over `data`,
  through `model_design()` / `match_rows()` (R/design.r). scikit-bio uses
  patsy with an intercept and **treatment contrasts whose reference is the
  first category in sorted order**. For parity: character columns are
  sorted (C locale); factors keep their level order (the R way to choose a
  reference; document it). Interactions: support what `model_design()`
  supports and check patsy's column naming for fixtures.
- `groups`: the formula term used by the post-hoc tests (scikit-bio's
  `grouping`; ≥ 3 levels). Doesn't change the primary results.
- `pseudocount`:
  - `ancombc()`: scikit-bio requires strictly positive input and
    documents adding 1 (R's `pseudo = 1`) for zeros. ecodive uses the
    shared rule (compositional.md); fixtures pass an explicit value.
  - `ancombc2()`: scikit-bio's default 0, with zeros treated as
    **missing** (its "sparse" estimation path). Keep 0 as the default;
    don't apply the shared `NULL` rule.
- `max_iters`, `tol`: the bias-estimation EM (scikit-bio's `max_iter`,
  `tol`; glossary names).
- `var_quantile`, `alpha`, `p_adjust`: as scikit-bio.
- Not carried over: random effects (the R package's `rand_formula` for
  repeated measures; scikit-bio has no mixed model), scikit-bio's
  `aggregator` (feature aggregation
  before the final regression; users can aggregate first) and the R
  package's preprocessing (`prv_cut`, `lib_cut`, `tax_level`), which
  scikit-bio also leaves to the user; show them as examples.

### Algorithm

Port scikit-bio's `_ancombc_core()` and its helpers as they are:
`_transform_data`, `_estimate_params` (dense, and sparse for missing
entries), the batched least squares (`_lstsq_*`, `_apply_pinv`,
`_coef_estimable`), `_calc_var_cov`, the bias EM `_estimate_bias_em`,
`_sample_fractions` (ANCOM-BC2), `_adjust_variance`, `_calc_statistics`,
`_calc_pvalues`; then for the post-hoc tests `_global_stats`,
`_pairwise_test` / `_mdfdr_pairwise`, `_dunnett_test` / `_mdfdr_dunnett`,
`_trend_test` (which uses `scipy.optimize.minimize(method="SLSQP")` in
place of R's `quadprog::solve.QP`: a small dense QP; write a QP solver
or port a permissive one), and `struc_zero`. Read the 0.7.4 source
closely: it is ~3,900 lines including docs, with many numerical
safeguards that parity depends on.

- **C:** per-feature least squares and variances are independent: thread
  over feature blocks. The EM and the per-covariate summaries are small;
  keep them in R unless profiling says otherwise.
- **Random numbers:** the mdFDR bootstraps (pairwise, Dunnett) are random
  in scikit-bio (NumPy); use PCG with `seed` and match those by
  distribution only. Everything else is deterministic and must match at
  about 1e-10.

### Result (shared per-feature structure, compositional.md)

- `data.frame` of class `c('ecodive_ancombc', 'data.frame')` (also for
  `ancombc2()`, with attribute `method = 'ancombc2'`), one row per
  **feature × covariate** (feature-major, as scikit-bio's multi-index):
  `feature`, `covariate`, `log_fc` (natural log, as scikit-bio 0.7.4's
  `Log(FC)`), `se`, `w`, `p_value`, `p_adjusted`, `reject`.
- Attributes: `method`, `formula`, `groups`, `alpha`, `p_adjust`,
  `max_iters`, `tol`, `var_quantile`, `pseudocount`, and the fitted
  quantities the post-hoc tests need (scikit-bio's `_beta_hat`,
  `_var_hat`, `_vcov_hat`, `_dof`, `_estimable`, the design), kept as
  attributes so `fit` is self-contained.
- Post-hoc results: the same structure, with the comparison in place of
  `covariate` (e.g. `contrast` "b - a"), per scikit-bio's tables.
- `structural_zeros()`: a logical feature × group data frame, as
  scikit-bio's.


## Reference implementations

- **Port from:** scikit-bio 0.7.4 `_ancombc.py` (BSD-3, derived from
  Artistic-2.0 ANCOMBC; see "Licensing").
- **Fixtures:** `data-raw/expected/ancombc.py` with scikit-bio (one and
  two covariates, a categorical with 3 levels and a numeric, zeros
  handled by pseudocount for `ancombc()` and as missing for `ancombc2()`,
  each post-hoc test); and `data-raw/expected/ancombc.r` with the R
  package ANCOMBC (Artistic-2.0, portable; see "Licensing") on the same data, since
  scikit-bio says it matches it. Compare at 1e-8 (the EM is iterative).


## Tests

- Fixtures; feature × covariate order; reference levels; zeros as missing
  in `ancombc2()`; a feature with too few observations (non-estimable →
  NA rows, as scikit-bio).
- Scale invariance per sample (the sampling fraction absorbs it).
- Identical for `cpus` 1, 2 and 4; post-hoc bootstraps identical for the
  same `seed`.
- Errors: a `formula` term missing from `data`; `groups` not a term or
  with < 3 levels for post-hoc; NA in covariates; zeros in `ancombc()`
  with `pseudocount = 0`.
- Smoke-test lines.


## Docs

- Rd pages in "Compositional Data" → "Differential abundance".
- Explain the sampling-fraction bias in two sentences; the reference
  level rule; how to reproduce the R package's preprocessing.
- NEWS.


## References

- Lin, H., & Peddada, S. D. (2020). Analysis of compositions of
  microbiomes with bias correction. *Nature Communications*, 11, 3514.
- Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions
  of microbiomes with covariate adjustments and repeated measures.
  *Nature Methods*, 21, 83–91.
- Nearing, J. T., Douglas, G. M., Hayes, M. G., et al. (2022). Microbiome
  differential abundance methods produce different results across 38
  datasets. *Nature Communications*, 13, 342.


## Decisions made while implementing (item 22)

Asked of the user at the start (all four answered with the recommendation):

- **Covariate names:** R's `model.matrix()` names (`(Intercept)`,
  `statusmoderate`, `age`, `statusmoderate:age`, `flagTRUE`), in R's order
  (patsy puts categorical terms before numeric ones; R follows the formula).
  The design is `model.matrix()` with treatment contrasts forced for every
  factor (ordered ones too), so its columns are patsy's; any formula
  `model.matrix()` takes is allowed, including interactions and `- 1`.
  Intercept rows are kept, as both references report them. Fixtures are
  matched by translating patsy's names (`status[T.moderate]` →
  `statusmoderate`, `[T.True]` → `TRUE`, `Intercept` → `(Intercept)`).
- **Non-estimable coefficients:** as scikit-bio and ANCOMBC: `log_fc`, `se`,
  `w` NaN, `p_value` 1 (counted in the adjustment), `reject` FALSE.
- **The fit is always stored, in full** (attributes): `coef` (D × P,
  bias-corrected, unmasked, as scikit-bio's `_beta_hat`), `var` (D × P;
  `ancombc2()`'s adjusted variances, as `_var_hat`), `vcov` (D × P × P, every
  feature's full covariance, `var` on the diagonal; scikit-bio keeps only the
  `grouping` block), `dof` (NULL for `ancombc()`; a per-feature vector for
  `ancombc2()`, NaN where observed samples ≤ rank), `estimable` (NULL, or a
  logical D × P matrix, rows all TRUE or all FALSE), `design` (the n × P
  model matrix with `assign` and `contrasts`), `groups` (the term name or
  NULL), plus `bias` (P × 3: `delta_em`, `delta_wls`, `var_delta`),
  `n_iters` (EM iterations per covariate) and, for `ancombc2()`,
  `sample_fractions`. Memory 8 D P² bytes for `vcov`.
- **`ancombc()` zeros:** the shared rule (`pseudocount = NULL`: nothing
  added without zeros, else half the smallest non-zero value with a
  warning). An explicit `pseudocount = 0` is allowed for a zero-free table
  and an error otherwise ("`counts` has zeros, and `pseudocount = 0` leaves
  their logs undefined.").

Decided along the way:

- **Signatures** as planned: `ancombc(counts, formula, data, groups = NULL,
  margin, pseudocount = NULL, max_iters = 100L, tol = 1e-5, alpha = 0.05,
  p_adjust = 'holm', cpus)`; `ancombc2()` adds `var_quantile = 0.05` and has
  `pseudocount = 0`. `formula` and `data` are required. `formula` may be a
  string, with or without `~` (scikit-bio takes `'status + age'`); a
  left-hand side is dropped. `data` may be a vector (one variable,
  `group`), and is matched to the samples by row name (`match_rows()`).
  `ancombc2(pseudocount = NULL)` is an error (no shared rule there).
- **Design helper:** `ancombc_design()` in `R/ancombc.r`, not
  `model_design()` (which centres and drops the intercept). It reuses
  `match_rows()` and `model_design()`'s checks and messages (variables not in
  `data`, NA, Inf). Character columns → factors with `sort(method =
  'radix')` levels; logical → factor(FALSE, TRUE); factors keep their level
  order, with unused levels dropped and ordered-ness removed. A categorical
  variable with one value is an error. `n <= P` is an error ("needs more
  samples than design columns"); a rank-deficient design is a warning (its
  coefficients come from the pseudoinverse, as scikit-bio's).
- **`groups` checks** (scikit-bio's `_validate_grouping()`): a column of
  `data`, a term label of `formula`, at least 3 distinct values, and at least
  2 design columns. It changes nothing else.
- **p-values:** z (`2 pnorm(-|w|)`) for `ancombc()`, t with each feature's
  dof for `ancombc2()`; `stats::p.adjust()` per covariate over the
  non-NaN p-values (scikit-bio's rule; R's default would count NaNs in `n`),
  any `p.adjust.methods` via `pairwise_method()`. `reject = p_adjusted <=
  alpha` (scikit-bio's `<=`, unlike `ancom()`'s `<`).
- **Quantiles:** numpy's linear method, including its interpolation form
  (`b - d (1 - t)` for `t >= 0.5`), in C (EM) and in R (`np_quantile()`, for
  `var_quantile`).
- **C (`src/ancombc.c`), one entry point `C_ancombc`:** the log table (n × D
  doubles, NaN = missing), the initial fit, the EM, the sampling fractions and
  the final fit; R builds the design, its pseudoinverse and `(X'X)+ 1` with
  `svd()` (numpy's 1e-15 cutoff), and does the variance adjustment and the
  statistics. Features are split into blocks of 64; per-sample sums over
  features (θ, sampling fractions, the convergence sums) are kept per block
  and added in block order, so results are identical for any `cpus`
  (checked for 1, 2, 3, 4, 8) and for dense, `Matrix`, `slam` and `margin =
  2` input. The EM runs one covariate per task.
- **Per-feature masked-design SVD:** a one-sided Jacobi SVD written for
  ecodive (no LAPACK in C), with scikit-bio's rank rule, its compact operator
  `V S⁻² V' X'r`, and the full pseudoinverse for features whose retained
  condition number exceeds 1e4. Computed once and reused by ANCOM-BC2's
  final fit (same missing pattern).
- **R-compatible failed fits and rebasing** (`_r_fit_info()` /
  `_r_rebase_categorical()`, scikit-bio's `match_r = True`, its only mode):
  R passes, per factor with a main-effect term, its level codes, its columns,
  whether it is treatment-coded (L − 1 columns; full coding without an
  intercept is the "custom" path: only the ≥ 2 observed levels check), and
  whether it may be rebased (intercept present, not in an interaction).
- **EM minimiser:** a port of SciPy 1.18.1's `_minimize_neldermead()` for
  one variable with bounds [0, ∞) (clipping, the 1.05 / 0.00025 initial
  simplex, 200 evaluations, xatol = fatol = 1e-8, numpy's argsort order for
  ties and NaN), credited to SciPy (`inst/COPYRIGHTS`, `Authors@R`).
- **Licensing:** `src/ancombc.c` and `R/ancombc.r` carry scikit-bio's BSD-3
  text, ANCOMBC's Artistic-2.0 notice (copyright "Huang Lin": ANCOMBC's
  DESCRIPTION names him as author and gives no years) with the "modified
  version" statement, and SciPy's BSD-3 (C file). One `inst/COPYRIGHTS`
  stanza holds both licences in full (the Artistic-2.0 text from R's
  `share/licenses`); SciPy's stanza now also lists `src/ancombc.c`. Huang
  Lin is `ctb`/`cph` in `Authors@R`.
- **ANCOMBC reference:** 2.14.0 (Bioconductor 3.23) needs R ≥ 4.5, so the
  fixtures use **2.8.1** (Bioconductor 3.20), the newest that installs on
  R 4.3.3 (with CVXR 1.0-15 from the CRAN archive: CVXR 1.0-17 needs
  Matrix ≥ 1.7). Its `R/ancombc_bias_correct.R` (`.iter_mle()`,
  `.bias_em()`) is identical to 2.14.0's; `ancombc2()`'s flow around it was
  restructured in 2.14.0 but computes the same quantities. `ancombc()` there
  always adds 1, so it is compared only where scikit-bio's table had 1
  added. ANCOMBC drops a one-column `meta_data` to a vector; the fixture
  script adds a `sample_id` column.

### Reference results

- **scikit-bio 0.7.4** (NumPy 2.5.5, SciPy 1.18.1, pandas 3.1.1, patsy
  1.0.3), 24 cases (both docstring examples, a 42 × 60 three-site table with
  numeric and logical covariates, a 30 × 25 sparse table with features
  missing from the reference group, from all but one group, or observed once,
  a factor with its own level order plus a covariate in the tens of
  thousands, interactions, a no-intercept model, proportions; `p_adjust`
  holm / BH / bonferroni / BY / none, `var_quantile` 0 / 0.05 / 0.1 / 0.2,
  pseudocounts 0 / 0.1 / 0.5 / 1, `max_iters = 3`): largest differences
  `log_fc` 8e-9 (absolute), `se` 7.9e-10 (relative), p-values 3.1e-8
  (relative), the biases 8e-9, sampling fractions 4.3e-9, `coef` 8e-9, `var`
  and `vcov` 1.6e-9 (relative); `reject`, `dof` and `estimable` identical
  everywhere. The EM's iteration counts are equal except one covariate (93 vs
  92, at the tolerance's edge). Where the EM converges early the agreement is
  1e-12 to 1e-14; the 1e-8 cases are those where it runs to `max_iters`, and
  the difference is summation order (numpy's pairwise and BLAS sums) moving
  the Nelder-Mead path.
- **ANCOMBC 2.8.1**, 16 comparable cases: where its EM converges for every
  covariate (7 cases, tested), `log_fc` within 3e-6, `se` 3e-6 (relative),
  p-values 4e-5 (relative), `reject` identical: nloptr's Nelder-Mead
  (default tolerances) vs SciPy's. Larger differences, all explained:
  - **EM not converged** (a covariate reaching 100 iterations): `log_fc` up to
    0.02 (`doc3`; 4e-3 `sim_bc2`, 2e-3 `sparse_bc2`). Running ANCOMBC's
    `.iter_mle()` on `doc3` gives scikit-bio's initial fit to every printed
    digit, and its `.bias_em()` on those inputs gives scikit-bio's bias to
    1e-8–1e-13 for the converged covariates and 0.02 off for the one that
    isn't.
  - **`var_quantile = 0`:** ANCOMBC's `s0_perc = 0` adds the smallest
    variance (`quantile(x, 0)`), scikit-bio's `if var_quantile:` adds
    nothing (se up to 41% apart).
  - **Interactions with zeros as missing:** scikit-bio deliberately doesn't
    rebase factors in interactions; R's `lm()` does in effect (one feature
    each in `sim_bc2_ia` and `levels_bc2_ia`).
  - **A covariate in the tens of thousands** (`levels_bc2`, `days`):
    ANCOMBC's `MASS::ginv(X'X)` (tolerance 1.5e-8) drops 4 of the Gram
    matrix's 5 singular values (they span 1e-10), so its standard errors are
    wrong (74 decisions differ); scikit-bio's pseudoinverse of X (1e-15 on
    singular values spanning 1e-5) keeps them. Documented.
- **Scale invariance:** without missing values, rescaling samples changes
  the slopes by < 1e-6 (tested). With zeros as missing values, the
  alternating fit stops at ANCOMBC's 0.01 tolerance, so per-sample rescaling
  moves `log_fc` by up to 0.12 in the `sparse` table: scikit-bio moves by the
  same 0.12 (ecodive's rescaled results equal scikit-bio's rescaled ones).
  Tested with a common factor instead; documented.
- **Timing** (500 samples × 5,000 features, 2% zeros, 5 covariates, -O2):
  `ancombc()` 0.5 s on 1 core, 0.4 s on 4; `ancombc2()` with zeros as
  missing 1.0 s / 0.5 s. scikit-bio 0.7.4 (vectorised NumPy): 0.5 s and
  1.2 s. Memory: the n × D log table, 8 D P² bytes for `vcov`, and with
  missing values 8 D P² more for the per-feature SVDs.

### For item 23 (post-hoc tests)

The result carries everything scikit-bio's `ANCOMBCResult` keeps, for any
fit (with or without `groups`):

- `attr(fit, 'coef')` = `_beta_hat` (D × P, R covariate names as
  colnames), `attr(fit, 'var')` = `_var_hat`, `attr(fit, 'dof')` = `_dof`
  (NULL for ANCOM-BC, which uses z; a length-D vector for ANCOM-BC2),
  `attr(fit, 'estimable')` = `_estimable` (NULL = all), `attr(fit,
  'alpha')`, `attr(fit, 'p_adjust')`, `attr(fit, 'method')`.
- `_vcov_hat` is `attr(fit, 'vcov')[, gi, gi]`, where the grouping's
  columns are `gi <- which(attr(attr(fit, 'design'), 'assign') ==
  match(attr(fit, 'groups'), attr(terms(attr(fit, 'formula')),
  'term.labels')))` — scikit-bio's `_groups` (0-based there). Because the
  full covariance is stored, a post-hoc test could also accept a `groups`
  term the fit wasn't given (validate it as `ancombc_design()` does).
  `_posthoc_estimable()` is `apply(estimable[, gi], 1, all)`.
- `attr(fit, 'design')` is `_dmat` (the model matrix; its `assign`,
  `contrasts`, and the column names give the comparisons' labels).
- The tests compare `coef`, `var`, `vcov` (grouping block), `dof` and
  `estimable` with scikit-bio's `_beta_hat`, `_var_hat`, `_vcov_hat`,
  `_dof`, `_estimable` for every fixture case with a grouping (in
  `tests/tinytests/expected/ancombc.rds`, `$skbio[[case]]$coef` etc.), so
  item 23's fixtures can start from the same cases (`data-raw/expected/
  ancombc.py` builds them; the post-hoc calls can be added there).

### Not changed

- The draft article `vignettes/articles/compositional.Rmd` still lists
  `ancombc()` / `ancombc2()` as planned (its chunk is `eval = FALSE`; its
  call `ancombc2(counts, ~ group, data = ...)` matches the signature). Left
  for the article's own update, as items 19 and 20 did.


## Decisions made while implementing (item 23)

Item 23 was given to the item-22 session after item 22 was pushed (the user's
message). Asked of the user at the start (all four answered with the
recommendation):

- **Comparison labels:** level names, `'moderate - mild'` (the first more
  abundant when `log_fc > 0`); pairs as `'severe - moderate'`. Without an
  intercept (full coding) a level's own coefficient is labelled with the
  level alone (`'gut'`). scikit-bio: `status[T.severe]_status[T.moderate]`.
- **`structural_zeros()` returns a logical matrix**, features × groups with
  dimnames.
- **`groups` can be overridden** at post-hoc time (`groups =
  attr(fit, 'groups')` by default): any categorical main-effect term with ≥ 3
  observed levels and ≥ 2 design columns. A fit without `groups` needs one
  given ("`groups` is required").
- **`n_bootstraps = 100`** for Dunnett's and the trend test (scikit-bio's and
  ANCOMBC's default); the docs recommend about 10,000 for the trend test.

Decided along the way:

- **Signatures:** `ancombc_global(fit, groups, alpha, p_adjust, cpus)`,
  `ancombc_pairwise(...)` (same), `ancombc_dunnett(..., n_bootstraps =
  100L, seed = 0L, cpus)`, `ancombc_trend(..., patterns = NULL, nodes =
  NULL, n_bootstraps = 100L, seed = 0L, cpus)`, with `alpha` and `p_adjust`
  defaulting to the fit's; `structural_zeros(counts, groups, data = NULL,
  margin = 1L, neg_lb = FALSE)`. `n_bootstraps` follows the `n_<unit>` rule
  (`n_draws`, `n_starts`). `patterns` / `nodes` replace scikit-bio's
  `trend_contrast` / `trend_node` (`contrast` already means the two groups of
  `dirmult_ttest()`): a named list of constraint matrices (`C %*% beta >= 0`,
  one column per coefficient of `groups`) and 1-based node indices
  (default: the last). Default patterns `increasing` (`diag(k)` with -1 on
  the subdiagonal) and `decreasing` (its negative), as scikit-bio.
- **Pairwise is deterministic** (scikit-bio and ANCOMBC screen it with the
  global F/χ² test, not a bootstrap), so it takes no `seed`; only Dunnett's
  and the trend test draw.
- **Results:** data frames of class `c('ecodive_ancombc_<test>',
  'data.frame')`: global and trend one row per feature (`feature`, trend's
  `pattern`, `w`, `p_value`, `p_adjusted`, `reject`); pairwise and Dunnett
  one row per feature × comparison (`feature`, `comparison`, `log_fc`, `se`,
  `w`, `p_value`, `p_adjusted`, `reject`). Attributes `method`, `groups`,
  `alpha`, `p_adjust` (+ `n_bootstraps`, `seed`); the trend test also
  `patterns`, `nodes` and `coef` (each feature's constrained estimates of
  its winning pattern). The `pattern` column is new (scikit-bio computes it
  and drops it).
- **scikit-bio's rules kept, quirks included:** the global p-value is
  `2 min(P, 1 - P)` (as ANCOMBC); mdFDR's screening uses BH whatever
  `p_adjust` is; `p_value` in the pairwise/Dunnett tables is the screened
  value (1 for unscreened features, and 0 or NaN become 1, as ANCOMBC);
  each row is adjusted with `p.adjust(n = m D / R)` (fractional `n`, as R's
  `p.adjust()` allows, which is what scikit-bio reproduces); NaN adjusted
  values become 1; non-estimable features get p = 1.
- **The item-22 fit gained two design attributes** for this:
  `attr(design, 'term.labels')` (so `~ .` formulas work) and
  `attr(design, 'xlevels')` (each factor's levels, for the labels). Nothing
  else in `ancombc()` / `ancombc2()` changed.
- **C** (appended to `src/ancombc.c`, reusing its Jacobi SVD):
  `C_ancombc_global` (b' V⁺ b with numpy's 1e-15 pinv cutoff),
  `C_ancombc_dunnett` (exceedance counts), `C_ancombc_trend` (constrained
  estimates, statistics, exceedance counts), `C_ancombc_nonzero` (non-zero
  counts per feature and group, any table type). One task per feature;
  each feature's draws come from PCG stream (seed, feature index), so
  results are identical for any `cpus` (checked for 1–4) and a feature's
  results don't depend on the other features (tested). The mdFDR and
  p-value steps are in R.
- **Samplers** (`src/rng_dist.c`, item 20): N(0, 1) by `rng_norm()`; t(ν)
  as `Z / sqrt(2 G(ν/2) / ν)` with `rng_gamma()` (NumPy's construction);
  failed fits' NaN dof are 999, as scikit-bio and ANCOMBC. Checked by
  distribution: with 200,000 draws the exceedance rate of |t| matches
  `2 pt(-w, ν)` within 5 SE for ν = 3, 12.5, 999 and ∞, w = 0.5, 1.5, 2.5.
- **QP solver** (ecodive's own, MIT): a primal active-set method (Nocedal &
  Wright 2006, Alg. 16.3) for `min 0.5 x'Gx + c'x` s.t. `Cx >= 0`, started
  at x = 0 (always feasible: the constraints are homogeneous), with the
  working-set multipliers from a pseudoinverse (dependent constraints
  allowed). It solves both the observed problems (G = 2V⁻¹ + 1e-10 I, as
  scikit-bio and ANCOMBC) and the bootstrap projections (G = I). Exact up to
  rounding: matches brute force over every active set to 1e-9 on all
  features of two tables, including a custom 5-constraint pattern on 3
  coefficients. V⁻¹ by Cholesky, or (V + 1e-8 I)⁻¹ by Gauss-Jordan when it
  isn't positive definite (`_safe_inverse_spd()`). No port, so no licence
  question (quadprog is GPL; SciPy's SLSQP would have been a large port).
- **`structural_zeros()`:** groups as `ancom()` (factor levels kept, unused
  dropped; else `sort(method = 'radix')`). scikit-bio sorts even a pandas
  Categorical's categories; ecodive keeps a factor's level order. NA in
  `counts` is an error (scikit-bio counts it as 0).

### Reference results

- **scikit-bio 0.7.4**, every fixture case with a grouping (19):
  - Global test: W within 1.8e-8 (relative), p 2.8e-8, decisions
    identical. Pairwise: `log_fc` 8e-9 (absolute), `se` 4.6e-10,
    `p_adjusted` 3.1e-8 (relative), decisions identical. Dunnett's:
    `log_fc`, `se`, `w` likewise. All inherit the fit's 1e-8; on the
    converged cases they agree to 1e-11–1e-13.
  - **Bootstraps on scikit-bio's own null draws** (its generator wrapped to
    record them; cases doc2 (z), doc3 (t) and sparse_bc2 (non-estimable
    features)): Dunnett's `p_adjusted` within 2.5e-12, trend p-values
    identical, decisions identical.
  - **Trend statistic:** within 3.2e-7, except the no-intercept case, 5.4e-5:
    SLSQP (ftol 1e-12) stops short of the exact solution; ecodive's is the
    minimum (feasible, and 20,000 random feasible perturbations never lower
    the objective).
  - **Trend p-values with ecodive's own draws:** 20,000 bootstraps each, on
    doc3, within 5 SE of scikit-bio's.
  - Docstring values reproduced (global, pairwise, trend W; struc_zero's
    example).
  - `struc_zero()`: identical on 6 tables × `neg_lb` (columns by name).
- **ANCOMBC 2.8.1** (`ancombc2(global, pairwise, dunnet, trend = TRUE)`,
  `mdfdr_control$fwer_ctrl_method = p_adjust`; trend with CVXR/ECOS), where
  its EM converges (doc3 with pseudocounts 0.1/0.5/1, sim_bc2_pc; doc2's
  `ancombc(global = TRUE)`): global W within 4e-8, pairwise `log_fc` 3.5e-9
  and `q` 2.7e-8, decisions identical; trend W within 3.9e-5 (ECOS's
  precision); Dunnett's decisions identical except one (sim_bc2_pc), from
  the random screen. ANCOMBC's trend bootstrap is one CVXR solve per
  feature, pattern and draw, so the fixture script runs it with 2
  bootstraps (only the statistic is compared).
- **Timing** (500 samples × 5,000 features, 4 groups, an `ancombc2()` fit):
  global 0.01 s, pairwise 0.2 s, Dunnett 0.3 s (the per-feature
  `p.adjust()` loop in R dominates both), trend 0.7 s on 1 core / 0.25 s on
  4 with 100 bootstraps, 64 s / 22 s with 10,000.

### Open

- `tests/tinytests/expected/ancombc.rds` is now 266 KB (xz), the largest
  fixture: scikit-bio's post-hoc tables and the recorded null draws of
  three cases. It could be trimmed (fewer cases with draws) if CRAN's
  tarball size matters.
- The draft article still lists the ANCOM-BC functions as planned.


## Fixtures slimmed (2026-11-19, the user's request)

The tests should be small and quick, so `tests/tinytests/expected/ancombc.rds`
went from 266 KB to 73 KB (xz), and `data-raw/expected/ancombc/` from 1.8 MB
to 0.6 MB:

- **Smaller tables:** `sim` 24 × 16 (was 42 × 60), `sparse` 21 × 12 (30 ×
  25, its special features kept), `levels` 18 × 10 (36 × 30), `prop` 12 × 8
  (20 × 15); 121 samples over all tables.
- **Fewer cases (22):** `doc3_pc0.1` and `sim_bc_bh` dropped (`sim_bc` now
  uses BH and alpha 0.1).
- **Less stored:** no `w` (tested as `log_fc / se`); `_beta_hat`,
  `_var_hat` and `_vcov_hat` for four cases (`FIT_CASES`); Dunnett's
  tables only where draws are recorded (its `log_fc` and `se` are the
  pairwise table's level-against-reference rows, tested as such); the
  trend test's statistic everywhere, its p-values only with draws; 20
  recorded null draws (was 50); ANCOMBC only the compared columns.
- **ANCOMBC cases** are chosen by the fixture script: those whose primary
  `log_fc` agree with scikit-bio's to 1e-5 (9 cases, 4 with post-hoc
  tests). On the new `sim_bc2_pc`, ANCOMBC's post-hoc step fails (its mdFDR
  calls `p.adjust(n = Inf)` with BY when no feature passes its screen), so
  the script runs it without post-hoc tests there.
- **Lighter Monte Carlo checks:** the t / normal sampler with 50,000 draws
  per setting (was 200,000); the trend p-values in distribution with 4,000
  of ecodive's draws against scikit-bio's 20,000; the many-feature
  threading table 30 × 300 (was 40 × 700).
- **Two tolerances, with reasons in the tests:** scale invariance of
  `log_fc` to 1e-4 (the EM absorbs the rescaling's common shift only to its
  tolerance, 1e-5; on the smaller `prop` table its iteration count changes);
  the trend statistic against scikit-bio to 1e-3 (SLSQP stops up to 1.8e-4
  short; brute force over the active sets now also checks the no-intercept
  case).

Results on the slimmed fixtures: scikit-bio `log_fc` within 2.1e-8, `se`
2.5e-9, p-values 2.5e-8, global W 1.1e-8, pairwise `p_adjusted` 1.2e-8,
decisions identical; ANCOMBC `log_fc` within 2.4e-6, `se` 2.6e-6, p-values
4.4e-5. In the full suite the two test files take 1.0 s and 1.2 s (1,119
results); the suite has 15,185 results at home and 297 in CRAN mode, all OK.


## Sensitivity analysis (2026-11-19, after another session's finding)

Another session reported that ANCOM-BC2 needs ANCOMBC's pseudocount
sensitivity filter. Checked before acting:

- **ANCOMBC** (2.4.0 to 2.14.0): `pseudo_sens = TRUE` is the default in every
  release; 2.14.0 warns when it is off ("essential for reducing false
  positives"). Its procedure (2.10.1 and later, read in 2.14.0's
  `ancombc2.R`): the whole analysis again with pseudocounts 0.1, 0.5 and 1
  added to every value; `passed_ss` when the q-value is significant in all
  four fits or in none; `diff_robust = diff & passed_ss`; the same for the
  global, pairwise and Dunnett's tests, and the trend test takes the global
  test's `passed_ss`. scikit-bio leaves it out of `ancombc2()` but shows the
  same procedure in its docstring ("strongly recommend").
- **Pelto et al. 2025** (Brief. Bioinform. 26, bbaf130): quotes verified
  (ANCOM-BC2's split-data replication 35%, 21% between studies; "highly
  inconsistent", said of ANCOM-BC2, DESeq2 and edgeR). They used default
  settings, so most likely with the filter on: evidence that ANCOM-BC2 is
  inconsistent, not that the filter fixes it.
- **The 2026 preprint** (doi:10.64898/2026.09.01.748571, "A
  configuration-resolved benchmark of differential abundance analysis
  methods for human gut 16S rRNA microbiome data"): exists; search snippets
  say the filter controls the FDR and dropping it inflates it. bioRxiv
  rate-limited every fetch, so its 0.02 / 0.82 figures weren't read
  directly; the docs cite it without numbers. A LOCOM2 preprint (2026)
  reportedly found ANCOM-BC2's FDR inflated even with the filter when
  effects are small or absent.
- The "forest plot and sensitivity plot" the finding mentions don't exist;
  ecodive's plots are a volcano plot, heatmaps and dot charts.

The user's decisions: an argument, **off by default** (scikit-bio's
results stay the default); `reject` kept, with `passed_ss` and
`reject_robust` added; the refits stored so that the post-hoc tests use
them; for `ancombc()` too.

Done: `sensitivity = FALSE` on `ancombc()` and `ancombc2()` (`TRUE` = 0.1,
0.5, 1; or a vector of positive pseudocounts). Each refit is the same call
with that `pseudocount` and `sensitivity = FALSE` (warnings already given by
the main fit are suppressed), kept in `attr(fit, 'sensitivity')` (a list of
`pseudocounts` and `fits`). Significance in each fit is its `reject`. The
post-hoc tests repeat themselves on the stored fits with the same
`groups`, `alpha`, `p_adjust` (and `n_bootstraps`, `seed`); the trend test
takes the global test's check. Plots: robust decisions filled (bold in the
heatmap) and named, significant but not robust hollow (italic). The help
page has a section "False discoveries: the sensitivity analysis", saying
plainly that without it the decisions can include many false discoveries.
Checked against scikit-bio's docstring (doc3: three significant, none
robust) and against separately run fits.

**Default changed (same day, the user's decision):** the user asked why
the reference's default hadn't been followed. ANCOMBC's `ancombc2()` has
`pseudo_sens = TRUE` by default and scikit-bio has no such option; the port
followed scikit-bio, and item 22 listed the analysis as "not carried over"
without asking (it should have been asked, since one reference has it on
by default). Now `ancombc2(sensitivity = TRUE)` by default, as ANCOMBC,
with a warning when `counts` isn't whole numbers and the default
pseudocounts are used (they are on the scale of counts); `ancombc()` stays
`sensitivity = FALSE` (neither reference offers it there). scikit-bio's
columns are unchanged; the fixture comparisons pass `sensitivity = FALSE`.
The default roughly quadruples `ancombc2()`'s time and fit memory.


## Outcome (from the implementing session's final report)

### Item 22 (`ancombc()`, `ancombc2()`), done 2026-11-19

Commits `49a6e59` (work), `c4c9361` (tick), `5f80b89` (merge). Later:
`210a989` / `5856378` (pseudocount sensitivity, on by default in
`ancombc2()`); item 29 made `formula` / `data` default to `NULL`.

### Item 23 (post-hoc tests, `structural_zeros()`), done 2026-11-19

Commits `7d9b859` (work), `a0e4a93` (tick), `3773221` (merge); `c90b9cf`
(fixtures 266 KB → 73 KB).

- Tests (final tree): 15,259 at home, 304 in CRAN mode.
- scikit-bio 0.7.4: fold changes ≤ 2.1e-8, p-values ≤ 2.5e-8, decisions
  identical; global and pairwise ≤ about 3e-8; Dunnett on scikit-bio's own
  draws ≤ 2.5e-12; structural zeros identical on 6 tables. Trend statistic
  ≤ 3e-7 (one case 5.4e-5 where SLSQP stops short; ecodive's active-set QP is
  exact).
- ANCOMBC 2.8.1: fold changes ≤ 2.4e-6 on the 9 cases where it agrees with
  scikit-bio; global and pairwise ≤ about 3e-8 where its EM converges.
- Open: ANCOMBC crashes with BY and no feature passing its screen (worked
  around in the fixture script; ecodive is unaffected).
