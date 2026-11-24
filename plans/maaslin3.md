# `maaslin3()`: prevalence and abundance models with a median comparison

**Priority:** 49. **Tier:** port (maaslin3 1.4.0, Bioconductor, **MIT**:
`DESCRIPTION` "MIT + file LICENSE", `LICENSE` "YEAR: 2025 / COPYRIGHT
HOLDER: Harvard School of Public Health", and the full MIT text in
`R/utility_scripts.R`; verified 2026-11-20). lme4, lmerTest, multcomp and
survival, which it calls, are GPL: output-only. Read
[README.md](README.md) first, and [compositional.md](compositional.md).

> **Scope approved by the user (2026-11-20).** MaAsLin3 is a general
> per-feature regression framework (any formula, two models per feature,
> optional random effects), not a compositional method in the narrow
> sense, so porting it widens ecodive from "compositional DA" to
> "per-feature association modelling". The user approved that widening as
> part of making ecodive a one-stop shop for microbiome analysis: its
> family has the highest uptake, its median comparison *is* a
> compositional correction, and the separate prevalence model is new to
> ecodive.


## What it is

For each feature and covariate, MaAsLin3 (Nickols et al. 2026) fits two
models:

- **abundance**: a linear model of \eqn{\log_2} of the feature's
  *non-zero* relative abundance (zeros are missing), and
- **prevalence**: a logistic model of presence/absence, with a data
  augmentation against separation,

then tests each abundance coefficient against **the median coefficient of
that covariate over all features** (the "median comparison", on by default
for abundance), which reads the relative-abundance fit as an
absolute-abundance statement when the typical feature doesn't change. In
IBDMDB, 77% of its associations were with prevalence rather than
abundance. Pelto et al. (2025) found logistic regression on presence /
absence among the three best approaches. Uptake (retrieved 2026-11-19):
MaAsLin3 246 OpenAlex citations (36 in 2025 → 210 in 2026), 155 Europe PMC
mentions; MaAsLin2 2,850 OpenAlex 2022–26, the most cited compositional /
DA paper. MaAsLin2 itself is **not** planned (log-TSS regression with no
compositional correction; compositional.md).


## Names

`maaslin3()`, the package's name in lowercase (rule 1; users know
`maaslin3::maaslin3()`). `DA_TESTS` row: design `formula`, effect
`coef_abundance`, base 2, `zero_free_ok = FALSE` (the prevalence model
needs zeros); aliases `maaslin`, `maaslin3`. "Same names in other
packages": maaslin3 (Bioconductor), the original.


## Design

### Signature

```r
maaslin3(counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
         norm = 'percent', transform = 'log', standardize = TRUE,
         augment = TRUE, median_comparison_abundance = TRUE,
         median_comparison_prevalence = FALSE,
         median_comparison_abundance_threshold = 0,
         median_comparison_prevalence_threshold = 0,
         subtract_median = FALSE, warn_prevalence = TRUE, n_sims = 10000L,
         alpha = 0.1, p_adjust = 'BH', seed = 0L, cpus = n_cpus())
```

- Model convention of [da_test.md](da_test.md) (`groups` = term of
  interest, or alone `~ groups`); the shared `da_design()`. Reference
  levels: a factor's first level, else radix-sorted first (maaslin3
  instead *requires* `reference` for > 2 levels of a non-factor; say so).
- `norm`: maaslin3's `normalization`: `'percent'` = `TSS` (default),
  `'rclr'` = its `CLR` (it centres on the non-zero part only, exactly
  ecodive's rclr), `'none'`. Zeros (values ≤ its `zero_threshold = 0`) are
  missing in the abundance model in every case.
- `transform`: `'log'` (default; \eqn{\log_2}; error on values ≤ 0),
  `'plog'` (zeros → half the smallest observed value, then \eqn{\log_2}),
  `'none'`. Only meaningful with `'percent'` / `'none'` (rclr is already a
  log; `'log'` with rclr is an error, as in maaslin3).
- `standardize`: numeric covariates centred and scaled to SD 1
  (`collapse::fscale`).
- The other arguments keep maaslin3's names, snake case and all, because
  that is what its users type (as zero replacements keep zCompositions'
  codes); `alpha` is its `max_significance` (default **0.1**, the
  reference's), `p_adjust` its `correction` (`'BH'`).
- `n_sims`, `seed`: the median comparison's Monte Carlo (maaslin3 fixes
  10,000 draws from R's RNG); here PCG, so results are reproducible.
- **Not ported:** filters (`min_abundance`, `min_prevalence`,
  `max_prevalence`, `min_variance`: left to the user, as in `ancombc()`),
  plots, files, `feature_specific_covariate`, `unscaled_abundance`
  (spike-ins; could map to `load` later), `group()`, `ordered()`,
  `strata()` terms, contrast tests, random effects (below).

### Algorithm (read from maaslin3 1.4.0 `R/maaslin3.R`, `R/fit.R`, `R/utility_scripts.R`)

1. **Normalize** (`TSSnorm()` / `CLRnorm()`), zeros → `NA`;
   **transform** (`LOG()` = `log2`, `PLOG()`).
2. **Abundance model**, per feature: `lm(expr ~ formula)` on its non-`NA`
   samples; coefficient, SE, t p-value (residual df), from
   `summary(lm)`. A feature whose values are all equal → error "All linear
   values are the same", `NA`s. A factor whose first level has no
   observed sample for this feature → "No data points have the baseline
   factor level", `NA`s (`check_missing_first_factor_level()`).
3. **Prevalence model**, per feature: \eqn{y_s = 1} if present. With
   `augment` and both outcomes present: rows duplicated twice with
   \eqn{y = 1} and \eqn{y = 0}, weight \eqn{p/(2n)} each, original rows
   weight 1, \eqn{p} = `ncol(model.matrix) - 1` (`augment_data()`);
   `glm(binomial, weights)`; Wald z p-values. In C, collapse the three
   copies: sample \eqn{s} has prior weight \eqn{1 + 2a} (\eqn{a = p/2n})
   and response \eqn{(y_s + a)/(1 + 2a)}; the IRLS score and information
   are the same. Use glm's convergence rule (`epsilon = 1e-8` on relative
   deviance change, `maxit = 25`) and its binomial starting values.
4. **Median comparison** (`run_median_comparison_general()`), per
   covariate, over features without errors:
   - `use` = coefficients with p < 0.95; \eqn{m} = their median;
     \eqn{v_m = \frac{\pi}{2}\,\mathrm{var}(\beta_{use}) / n_{use}}
     (code: `0.25 * 2 * pi * var / n`).
   - Monte Carlo: `n_sims` draws of every feature's \eqn{\beta^*_j \sim
     N(\hat\beta_j, se_j^2)}; \eqn{c_j = \mathrm{cov}(\beta^*_j,
     \mathrm{median}(\beta^*_{use}))} (sample covariance, \eqn{n - 1}).
   - Test: \eqn{t_j = (\hat\beta_j - m)/\sqrt{se_j^2 + v_m - 2c_j}}
     (equivalent to maaslin3's `multcomp::glht(rhs = offset)`, with offset
     \eqn{\hat\beta + |m - \hat\beta|\,se/\sqrt{se^2 + v_m - 2c}}). ✱
     Reference distribution: t with residual df for `lm`, normal for `glm`
     (multcomp's documented behaviour as I understand it; **confirm from
     fixtures**).
   - \eqn{|\hat\beta_j - m| <} threshold (divided by `sd()` of the design
     column when it isn't 0/1) → p = 1; \eqn{m} undefined → p `NA` with
     the error "P-value became NA in median comparison".
   - `null` = \eqn{m}; with `subtract_median`, the coefficient becomes
     \eqn{\hat\beta - m}. Without a median comparison, `null` = 0.
5. **q-values** (`add_qvals()`): one `p.adjust()` over the abundance and
   prevalence p-values **together** (errored ones `NA`).
6. **Joint** (`create_combined_pval()`): \eqn{p_{joint} =
   \mathrm{pbeta}(\min(p_{ab}, p_{prev}), 1, 2) = 1 - (1 - \min)^2}; the
   other p if one is `NA`/errored, `NA` if both; \eqn{q_{joint}} =
   `p.adjust(p_joint)`.
7. **`warn_prevalence`** (`flag_abundance_turned_prevalence()`): refit the
   abundance model as TSS + LOG without median comparison; where its q <
   `alpha`, signs agree and \eqn{|\beta_{ab}| > |\beta_{prev}|}, flag the
   prevalence row "Prevalence association possibly induced by stronger
   abundance association".
8. Prevalence sanity flags (`maaslin_fit()`): features present in < 5% of
   samples and fewer than 50 (or absent from < 5% and fewer than 50) with
   |coef| > 15 or p < 1e-10 get maaslin3's warning text.

### C and threads

Abundance: per-feature masked OLS through `src/feature_lm.c`
([da_test.md](da_test.md)); prevalence: per-feature IRLS in
`src/maaslin3.c`; both threaded over feature blocks, deterministic. Median
comparison: features × `n_sims` normals from PCG streams keyed by
(seed, covariate, sim), the covariances accumulated per feature in sim
order; never the full matrix (10,000 × D). `pt()` / `pnorm()` on the main
thread.

### Result

`c('ecodive_maaslin3', 'data.frame')`, one row per feature × non-intercept
covariate (maaslin3 drops the intercept): `feature, covariate,
coef_abundance, se_abundance, null_abundance, p_abundance, q_abundance,
coef_prevalence, se_prevalence, null_prevalence, p_prevalence,
q_prevalence, n_nonzero, flag_abundance, flag_prevalence` (maaslin3's
`error` text, `NA` when none), then `p_value` (= `pval_joint`),
`p_adjusted` (= `qval_joint`), `reject`. maaslin3's "significant results"
rule is `q_joint ≤ α | (q_individual ≤ α & no error)` per model row; here
`reject = q_joint ≤ α | (q_abundance ≤ α & no abundance error) |
(q_prevalence ≤ α & no prevalence error)`, plus `reject_abundance` and
`reject_prevalence`. Abundance coefficients are log2; prevalence ones are
log odds. Every feature has rows (maaslin3 drops features it filters,
including, through `min_variance = 0`, those with fewer than two non-zero
values or constant non-zero values; here they get `NA` and a flag). `plot()`: volcano of either model
(`model = c('abundance', 'prevalence')`) and forest via
[da_plots.md](da_plots.md) helpers.

### Later: random effects

`(1 | subject)` terms via lmerTest / glmer in maaslin3: needs a REML
random-intercept LMM and a Laplace GLMM in C. Share the LMM with
[linda.md](linda.md) and [dirmult_lm.md](dirmult_lm.md); the GLMM is new
work. Separate step, after the fixed-effects port.


## Reference implementations

- **Port from:** maaslin3 1.4.0, functions named above. File header with
  the MIT notice; `inst/COPYRIGHTS` stanza; `Authors@R` `cph` "Harvard
  School of Public Health".
- **Fixtures:** `data-raw/expected/maaslin3.r` runs maaslin3 (all filters
  off, `plot_* = FALSE`) on (a) two groups, (b) three groups (a factor) +
  a numeric covariate, (c) CLR / NONE transform variants; store
  `fit_data_abundance$results`, `fit_data_prevalence$results`. For exact
  median-comparison p-values, ✱ reproduce maaslin3's `rnorm()` stream in
  the generator (`set.seed()` before the call; `replicate(nsims,
  rnorm(n_coefs, coefs, sigmas))` per covariate in `unique(paras$name)`
  order) and pass those draws to an internal `maaslin3_details(sims = )`.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Two models in one table: the `_abundance` and `_prevalence` columns sit
  side by side, and `plot(model = )` picks which to draw, highlighting
  from `reject_abundance` or `reject_prevalence`. Store `log_base` (`2`
  for the abundance coefficients, as MaAsLin 3 reports log2 — check;
  prevalence coefficients are log odds) and `alpha`.
- `plot()`: [da_plots.md](da_plots.md)'s volcano and forest types, one
  panel per model, titled from `method` via `DA_TESTS`, for example
  "Volcano plot of MaAsLin 3", with the model and covariate in the
  subtitle. The input is a count table, so no `display`.

## Tests

- Abundance coefficients, SEs, p at 1e-10; prevalence at 1e-6 (IRLS
  tolerance); median-comparison p at 1e-8 with injected draws, at
  Monte-Carlo tolerance (1e-2 relative) with PCG draws; q, joint p/q,
  flags and `reject` identical.
- `identical()` for `cpus` 1, 2, 4.
- Edge cases: a feature present everywhere ("All logistic values are the
  same"), absent from the baseline level, present in one sample;
  separation (augmentation keeps estimates finite); n ≤ P; NA in `data`.
- Smoke line.


## Docs

Rd page in *Differential abundance*: the two models; how to read the
median comparison (relative → absolute, assuming the typical feature is
unchanged; turn it off for absolute data such as `qmp_scale()` output);
the joint p-value; a "Differences from maaslin3" section (no filters;
reference levels; PCG Monte Carlo; every feature reported). Glossary:
`transform`, `augment`, `subtract_median`, `warn_prevalence`, `n_sims`,
the `median_comparison_*` names. NEWS.


## References

- Nickols, W. A., Kuntz, T., Shen, J., et al. (2026). MaAsLin 3: refining
  and extending generalized multivariable linear models for meta-omic
  association discovery. *Nature Methods*, 23(3), 554–564.
  doi:10.1038/s41592-025-02923-9 (246 OpenAlex, retrieved 2026-11-19).
- Mallick, H., et al. (2021). Multivariable association discovery in
  population-scale meta-omics studies. *PLoS Computational Biology*,
  17(11), e1009442. doi:10.1371/journal.pcbi.1009442 (MaAsLin2).
- Pelto, J., et al. (2025). *Briefings in Bioinformatics*, 26(2),
  bbaf130. doi:10.1093/bib/bbaf130.


## Decisions made while implementing (item 49, 2026-11-24)

Asked of the user at the start:

- **Random intercepts: included now** (`~ g + (1 | id)`), for both models,
  rather than in a later step: abundance through `src/feature_lmm.c`
  (REML at the exact optimum, Satterthwaite df, as lmerTest), prevalence
  through a new Laplace logistic mixed model (as `glmer()`), plus
  maaslin3's `small_random_effects` (the grouping as a fixed effect in the
  prevalence model).
- **`p_adjust = 'holm'`, `alpha = 0.05`** (ecodive's DA convention; maaslin3:
  `'BH'`, `max_significance = 0.1`). The docs give maaslin3's settings.

Settled from maaslin3 1.4.0's source (MIT) and from the outputs of lm(),
glm(), lme4 2.1.0, lmerTest 3.2-1, multcomp 1.4-25 and R's binomial()
family (GPL, run only; their code was not opened, except that printing
`binomial()$initialize` once showed its starting values, which the plan
already named):

- **Version.** maaslin3 1.4.0 (Bioconductor 3.23) installed on R 4.3.3
  from its tarball with `Depends: R (>= 4.4)` relaxed; it needs collapse
  >= 2.1 (2.0.10 fails in `rowbind()`), so collapse 2.1.9 and data.table
  1.18 came from CRAN. Every fixture is from 1.4.0.
- **The median comparison's test** (`multcomp::glht(rhs = offset)`): a t
  test on the fit's residual df for `lm`, the normal for `glm`/`glmer`,
  Satterthwaite's F (= t^2) for `lmer` (`lmerTest::contest`). Here
  `2 pt(-|t|, df)` with t computed as glht does, `(b - offset) / se`;
  multcomp's own numerical integral (mvtnorm) is off by up to ~1e-8
  relative, the largest difference in the fixed-effects fixtures.
- **Draw order.** maaslin3 calls `replicate(10000, rnorm(n, coefs,
  sigmas))` per covariate in design order, over the rows without an error,
  only when the covariate's median is defined. The tests regenerate those
  standard normals with `set.seed(seed); rnorm()` and give them to
  `maaslin3_details(sims = )`, reproducing maaslin3's medians to 1e-14 and
  covariances to 6e-16.
- **glm's IRLS**: convergence `|dev - dev_old| / (|dev| + 0.1) < 1e-8`,
  25 iterations, starts `(w y + 0.5) / (w + 1)`; the augmented copies
  collapse into one weighted row per sample (the first iteration combines
  the three copies' working weights and responses, which differ only
  there); the deviance is summed over the copies, so the iteration count
  is glm's; QR tolerance 1e-11 (glm's `min(1e-7, epsilon / 1000)`); the
  standard errors come from the last iteration's QR (weights before the
  final update). The logit's inverse is `e / (1 + e)` with `e = exp(eta)`,
  saturating to `eps` / `1 / eps` beyond |eta| = 30, and its derivative is
  `eps` there: settled from `binomial()$linkinv(c(29, 30, 31))` and
  `mu.eta`. Warnings become flags exactly as maaslin3 records them:
  "glm.fit: algorithm did not converge", "glm.fit: fitted probabilities
  numerically 0 or 1 occurred" (`mu` within 10 eps of 0 or 1; the later
  warning wins).
- **glmer's criterion with the augmentation**: `glmer(..., devFunOnly =
  TRUE)` equals, to 1e-15, -2 log L of the **original** samples only (the
  copies' fractional weights round to 0 binomial trials) + |u|^2 +
  log det(I + theta^2 Z'WZ), with the conditional modes u and working
  weights W of the augmented (weighted) data. So the augmentation only
  enters through u and the log determinant, and at theta = 0 the criterion
  is the plain, unaugmented likelihood.
- **glmer's standard errors**: 2 Hess^-1 of that criterion in (theta,
  beta), fixed-effects block (matches lme4's `vcov()` to its
  finite-difference Hessian, ~1e-6); but for a **singular fit** (theta <
  1e-4, lme4's `isSingular()`), lme4 stores no derivatives and uses the
  penalized least squares' (RX) covariance, X'WX - sum_j theta^2 / (1 +
  theta^2 S_j) (X'W 1_j)(X'W 1_j)' inverted, with the augmented weights.
  Both are implemented (checked on fits with theta = 0, 1e-7, 4e-5, 9e-5).
- **lm()'s per-feature failures**: a feature present at a single level of
  a factor makes lm() fail ("contrasts can be applied only to factors with
  2 or more levels"), for every coefficient; reproduced (status 7). The
  baseline check comes first, as maaslin3's.
- **The `warn_prevalence` refit** is maaslin3's percent + log abundance
  model without a median comparison; with the default `norm` and
  `transform` it is the main fit before its median comparison, so it isn't
  refitted. maaslin3 forbids `warn_prevalence` with CLR or a transform
  other than LOG: ported as errors naming `warn_prevalence = FALSE`.
- **The flag's tie**: a two-group prevalence coefficient that is 0 in
  theory (equal presence counts) is -1.4e-15 in glm and 0 here; its sign
  decides `warn_prevalence`'s flag, so such ties can flag differently. The
  fixture tables avoid exact ties (`unbalance()` in the generator).

Made while implementing:

- **Signature:** `maaslin3(counts, formula = NULL, data = NULL, groups =
  NULL, margin = 1L, norm = 'percent', transform = 'log', standardize =
  TRUE, augment = TRUE, median_comparison_abundance = TRUE,
  median_comparison_prevalence = FALSE,
  median_comparison_abundance_threshold = 0,
  median_comparison_prevalence_threshold = 0, subtract_median = FALSE,
  warn_prevalence = TRUE, small_random_effects = FALSE, n_sims = 10000L,
  alpha = 0.05, p_adjust = 'holm', seed = 0L, cpus = n_cpus())`. No
  argument shadows `da_test()` or `da_sensitivity()`
  (`test-wrapper_args.r` passes).
- **`transform`: `'log'` and `'none'` only.** maaslin3's PLOG replaces
  zeros, but zeros are already missing (`zero_threshold = 0`), so it equals
  LOG; maaslin3 itself refuses PLOG unless `zero_threshold < 0`.
- **`norm`: `'percent'` (TSS, proportions), `'rclr'` (its CLR: natural
  log, non-zero part), `'none'`.** `transform = 'log'` with `'rclr'` is an
  error, as maaslin3's.
- **Every feature is reported.** A feature with fewer than two non-zero
  values or constant non-zero values (maaslin3's `min_variance = 0` drops
  them) gets "All linear values are the same" for abundance and its
  prevalence model as usual, and that prevalence p-value joins the
  adjustments. With CLR, maaslin3's `min_abundance = 0` filter also drops
  features whose CLR values are all <= 0 (found from outputs); kept here.
  The fixture tables are made so that maaslin3 keeps every feature
  (checked in the generator).
- **Result:** the plan's columns plus `dof_abundance` (after
  `se_abundance`: residual or Satterthwaite df, for the forest plot's t
  intervals), and `reject_abundance` / `reject_prevalence` after `reject`.
  `reject` = q_joint <= alpha, or either model's own q <= alpha without a
  flag (the "<4 average observations" flag doesn't count), which is
  maaslin3's "significant results" over its two rows; each model's own
  reject is its q <= alpha without a flag. Flags keep maaslin3's texts;
  lme4's error texts for too few or too many random-intercept groups among
  a feature's samples ("grouping factors must have > 1 sampled level",
  "number of levels of each grouping factor must be < number of
  observations"), and ecodive's own "The mixed model did not converge" and
  "The model could not be fitted". Attributes: `method`, `formula`,
  `groups`, `design`, `random` (mixed), `norm`, `transform`, `alpha`,
  `p_adjust`, `log_base` (2 with log, e with rclr, NA otherwise), and per
  feature `mean` (mean transformed abundance where present), `prevalence`,
  `resid_sd`, `random_sd` (mixed). `print()` is `print_table_result`.
- **`DA_TESTS` row:** `MaAsLin 3, maaslin3, formula, coef_abundance, 2,
  FALSE, TRUE`, alias `maaslin` (letters only, so `'MaAsLin3'` and
  `'MaAsLin 3'` match too). `match_da_test('MaAsLin2')` (any other digit)
  is an error saying ecodive has MaAsLin 3; without it, letters-only
  matching sent MaAsLin 2 to maaslin3. `zero_free_ok = FALSE`
  (`da_sensitivity()`'s message now says "needs the table's zeros (it
  models counts, or presence)"). `da_sensitivity()` has nothing to vary for
  maaslin3 (no `pseudocount` or `scale_sd`; zeros refused): its three
  errors are tested.
- **Ordering:** rows by feature (table order), then covariate (design
  order), as the other model methods; maaslin3 sorts by q.
- **Reference levels:** `da_design()`'s (a factor's own order with unused
  levels dropped, else radix-sorted). maaslin3 requires `reference` for a
  non-factor with more than two levels. NA in `data` is an error (maaslin3
  drops the samples per fit), as in the other DA methods.
- **The intercept** is required (maaslin3 drops the first coefficient
  whatever it is); without one, every column is reported.
- **C (`src/maaslin3.c`):** `C_maaslin3_fit` reads the table by feature
  with `da_table_init()` (raw values), the per-sample divisors (TSS totals
  in long double, feature order; CLR's `exp(1)^mean(log x)` with R's
  two-pass mean) once; threads over features (interleaved, each feature
  independent): masked OLS (`flm_*`, tol 1e-7), or `lmm_fit()` on the
  feature's samples, its groups renumbered and its full-rank columns
  (`lmm_design_init_mem()`, added to `feature_lmm.c` so worker threads can
  set up a design in their own memory); `m3_glm()` (IRLS on the weighted
  QR); `m3_glmm()` (Newton on (theta, beta) with the analytic gradient
  of the Laplace criterion through the implicit conditional modes, the
  Hessian by central differences of that gradient, a backtracking line
  search, theta reflected to >= 0, start theta = 1 and the glm fit, as
  lme4's start). `C_maaslin3_median`: blocks of 256 draws, each draw's
  normals from stream `(model * 65536 + k) << 32 + draw`, the median by
  quickselect, per-block sums of z and z (m - c) (c the observed median,
  for accuracy), added in block order; covariance = se (S_zm - S_z S_m / N)
  / (N - 1).
- **Speed** (4 cores, -O2): 200 samples x 2,000 features, `~ g + age`,
  10,000 draws: 1.2 s on 1 core, 0.39 s on 4 (the median comparison is
  80%); with `(1 | id)` (50 subjects): 3.9 s and 1.1 s. maaslin3 fits one
  R model per feature.

Reference results (`data-raw/expected/maaslin3.r`, fixture
`tests/tinytests/expected/maaslin3.rds`, 53 KB; 10 cases on 5 tables: two
groups; a three-level factor + numeric + character with a feature absent
at the reference level, one only at the reference level and a rare one
(holm); CLR/NONE; NONE/LOG; TSS/NONE; both median comparisons with
thresholds and subtract_median; no augmentation, no median comparison,
unstandardized; a random intercept (12 subjects x 4); small groups (15 x
3) with and without small_random_effects), maaslin3's draws regenerated:

- Fixed effects: abundance coefficients 1.1e-13, standard errors 1.2e-15
  relative; prevalence coefficients 2.4e-13 (9.4e-10 on a coefficient near
  0, floor 1e-6), standard errors 6.3e-15; prevalence p 7.1e-14; median
  comparison p 3.4e-8 (multcomp's integral), and the q and joint values
  through it; nulls 1.8e-14. Without augmentation, a separated feature's
  25 non-converging iterations: 3.3e-6. Every flag and `reject` identical.
- Random intercept, abundance (lmerTest): coefficients 8.3e-7, standard
  errors 6.3e-8, p 4.4e-6 (lme4's optimizer).
- Random intercept, prevalence, against the best of glmer's nloptwrap,
  bobyqa and Nelder_Mead fits (lowest criterion): coefficients 6.6e-7
  absolute, standard errors 1.8e-6 relative (22 of 40 and 45 of 60 fits
  singular). Against maaslin3's own glmer fits, up to 1.4e-3 (12 x 4) and
  0.30 (15 x 3) absolute: maaslin3 keeps nloptwrap's fit, which stops
  short on flat criteria (ecodive's criterion is lower). `reject` still
  identical.
- With PCG's own draws: -log p within 2% of maaslin3's (Monte Carlo
  error, amplified in the tails: p up to 8% relative at p ~ 1e-6).
- Threads 1, 2, 4 (and mixed 1, 4), dense, `Matrix`, `slam`, `margin = 2`
  and double input: identical; `da_test(counts, 'maaslin3', ...)` and
  `'MaAsLin 3'` identical to the direct call.

Not done: `group()`, `ordered()` and `strata()` terms (clogit), contrast
tests, `feature_specific_covariate`, `unscaled_abundance` (could map to
`load`), several random effects or slopes, lme4's convergence-check
warnings (maaslin3 turns them into flags), `evaluate_only`. The draft
article `vignettes/articles/compositional.Rmd` doesn't mention
`maaslin3()` yet.
