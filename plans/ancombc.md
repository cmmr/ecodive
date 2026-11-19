# `ancombc()` and `ancombc2()`: ANCOM with bias correction

**Priority:** 24 (`ancombc()`, `ancombc2()`) and 25 (the post-hoc tests
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
