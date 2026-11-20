# `maaslin3()`: prevalence and abundance models with a median comparison

**Priority:** 42. **Tier:** port (maaslin3 1.4.0, Bioconductor, **MIT**:
`DESCRIPTION` "MIT + file LICENSE", `LICENSE` "YEAR: 2025 / COPYRIGHT
HOLDER: Harvard School of Public Health", and the full MIT text in
`R/utility_scripts.R`; verified 2026-11-20). lme4, lmerTest, multcomp and
survival, which it calls, are GPL: output-only. Read
[README.md](README.md) first, and [compositional.md](compositional.md).

> **SCOPE DECISION FOR THE USER.** MaAsLin3 is a general per-feature
> regression framework (any formula, two models per feature, optional
> random effects), not a compositional method in the narrow sense. Porting
> it widens ecodive from "compositional DA" to "per-feature association
> modelling". The review recommends it anyway (reconciliation table): its
> family has the highest uptake, its median comparison *is* a
> compositional correction, and the separate prevalence model is new to
> ecodive. Confirm before building.


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


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- **`da` record** "MaAsLin 3" via `call_step()`: `formula`, transform
  and normalization settings, `p_adjust`, `alpha`, median-comparison
  settings, `seed` (its Monte Carlo covariance).
- Two models in one table: a `model` column (`"abundance"`,
  `"prevalence"`) lets `plot()` facet or colour by model; store `log_base`
  (`2`, as MaAsLin 3 reports log2 coefficients — check) and `alpha`.
- `plot()`: [da_plots.md](da_plots.md)'s volcano and forest types, one
  panel per model, titled "MaAsLin 3: <covariate>".

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
