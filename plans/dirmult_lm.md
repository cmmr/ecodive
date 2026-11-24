# `dirmult_lm()`: linear models on Dirichlet draws, with scale models (ALDEx3)

**Priority:** 42. **Tier:** port (ALDEx3 1.4.0, CRAN, **MIT**:
`DESCRIPTION` "MIT + file LICENSE", `LICENSE` "YEAR: 2026 / COPYRIGHT
HOLDER: ALDEx3 authors"; verified 2026-11-20; already the source of
`dirmult_ttest()`'s scale model and credited in `src/dirmult.c` and
`inst/COPYRIGHTS`). Read [README.md](README.md) first,
[compositional.md](compositional.md), and [dirmult_ttest.md](dirmult_ttest.md).


## What it is

`dirmult_ttest()` compares two groups. ALDEx3 (McGovern & Silverman 2026)
fits **any linear model** to every Monte Carlo draw of the Dirichlet
posterior, after adding a **scale model** to each draw:
\eqn{\log_2 W = \log_2(\text{composition}) + \log_2(\text{scale})}. The
scale model states what is known about the samples' total loads and how
uncertain it is: the clr assumption with uncertainty (`clr.sm`), no
systematic difference (`tss.sm`), a measured or predicted load per sample
(`sample.sm`: flow cytometry, qPCR, spike-ins, or the MLP load predictor of
Nishijima et al. 2025), or a prior on a design coefficient
(`coefficient.sm`). Nixon et al. (2025) and Dos Santos et al. (2026) show
scale models restore FDR control; Nishijima et al. found microbial load
"the major determinant of gut microbiome variation". compositional.md's
2026-11-20 revision replaced the planned scikit-bio `dirmult_lme()` port
with this one.

Uptake is small (ALDEx3 / SR-MEM: 4 OpenAlex citations, 5–7 Europe PMC
mentions; Nixon 2025: 27 OpenAlex; retrieved 2026-11-19), but it extends
code ecodive already has and turns `load` (from `qmp_scale()`) into a
prior. Mixed models (SR-MEM) are a later step.


## Names

`dirmult_lm()`: the regression sibling of `dirmult_ttest()` (rule 2:
family prefix; compositional.md names it so). Not `aldex()` (rule 6, as
for `dirmult_ttest()`). `DA_TESTS` row: design `formula`, effect
`estimate`, base 2, `zero_free_ok = FALSE`; alias `aldex3` (**ask**, as
for `dirmult_ttest`'s `aldex2`).


## Design

### Signature

```r
dirmult_lm(counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
           pseudocount = 0.5, scale_model = c('clr', 'tss'),
           scale_sd = 0.5, load = NULL, load_sd = NULL, n_draws = 128L,
           se_type = c('hc3', 'hc0', 'classical'), alpha = 0.05,
           p_adjust = 'BH', seed = 0L, cpus = n_cpus())
```

- Model convention of [da_test.md](da_test.md); design from `da_design()`.
  ALDEx3 takes any model matrix; ecodive's treatment coding and level order
  apply.
- `pseudocount`: the Dirichlet prior, 0.5, as ALDEx3 (`Y + 0.5`) and
  `dirmult_ttest()`.
- **`scale_model`**, not `scale` (glossary: `scale` standardises
  features). `'clr'` = `clr.sm()`, `'tss'` = `tss.sm()`.
- **`scale_sd`**: their `gamma`, log2 units. Default **0.5**, ALDEx3's
  `clr.sm(gamma = 0.5)`, where `dirmult_ttest()` defaults to 0 (scikit-bio
  parity). **Ask the user** which default to use for both.
- **`load`**, **`load_sd`**: `sample.sm()`. `load` (glossary, as
  `qmp_scale()`): per-sample total load, > 0, or a column of `data`;
  `load_sd`: the log2 SD of its error, one value or one per sample
  (required with `load`). With `load`, `scale_model`/`scale_sd` are
  ignored (an error if both given explicitly). Docs show predicted load
  from the MLP (Nishijima 2025; not ported: its licence is ambiguous, see
  compositional.md) passed this way with `load_sd` from its prediction
  error.
- `coefficient.sm()` (priors on design coefficients, `c.mu`/`c.sd`): a
  later option, if asked for.
- `n_draws`: 128, as `dirmult_ttest()`; ALDEx3 defaults to 2,000 (docs:
  use ≥ 1,000 with a scale model). **Ask the user.**
- `se_type`: ALDEx3's `test` (`"t.HC3"` default, `"t.HC0"`, `"t"`).
- `p_adjust`: ALDEx3's `p.adjust.method = "BH"`.

### Algorithm (read from ALDEx3 1.4.0 `R/aldex.R`, `aldex-internals.R`, `fflm.R`, `scale.R`, `dirichlet.R`)

For each draw \eqn{d}:

1. **Composition**: \eqn{G_{sf} \sim \text{Gamma}(x_{sf} + 0.5, 1)},
   \eqn{\log_2 \text{comp}_{sf} = (\log G_{sf} - \log\sum_g G_{sg})/\log 2}
   (`rLogDirichlet()`).
2. **Scale** \eqn{\log_2 S_{s}}:
   - clr: \eqn{-\text{mean}_f \log_2\text{comp}_{sf} + x_s^\top\Lambda},
     \eqn{\Lambda \sim N(0, \gamma^2 I_P)}, one per design column
     including the intercept;
   - tss: \eqn{x_s^\top\Lambda};
   - load: \eqn{N(\log_2 \text{load}_s, \text{load\_sd}_s^2)},
     independent per sample.
3. \eqn{Y = \log_2 W = \log_2\text{comp} + \log_2 S} (n × D).
4. **OLS** (`fflm()`): \eqn{\Theta = (X^\top X)^{-1}X^\top Y};
   \eqn{e = Y - X\Theta}; with \eqn{A = (X^\top X)^{-1}X^\top},
   \eqn{h = \text{diag}(XA)}:
   - classical: \eqn{se_p = \sqrt{[(X^\top X)^{-1}]_{pp}\,\sum_s e_s^2/(N-P)}};
   - HC0: \eqn{se_p = \sqrt{\sum_s A_{ps}^2 e_s^2}};
   - HC3: \eqn{se_p = \sqrt{\sum_s A_{ps}^2 e_s^2/(1-h_s)^2}}.

   \eqn{t = \Theta/se}; `p_upper = pt(t, N - P, lower.tail = FALSE)`;
   `p_lower = 1 - p_upper` (computed exactly so, for parity).
5. **Aggregate** (`aldex.pvals()`, two-sided): per draw
   \eqn{\ell = \min(1, 2p_{lower})}, \eqn{u = \min(1, 2p_{upper})};
   `p_value` = \eqn{\min(\overline{\ell}, \overline{u})} (means over
   draws); `p_adjusted` = \eqn{\min(\overline{\text{adj}(\ell)},
   \overline{\text{adj}(u)})}, where `adj` is `p.adjust(p_adjust)` across
   features, per covariate and draw. **Note:** ALDEx3 averages *adjusted*
   p-values; `dirmult_ttest()` adjusts *averaged* ones. Follow each
   reference; document the difference.
6. `estimate`, `se`: means over draws (`summary.aldex()`).

### Reuse of `dirmult_ttest()`'s C

- Factor phase A of `src/dirmult.c` into a shared draw routine: the same
  PCG stream per (seed, draw, table sample), `d * (n + 1) + s`, so for the
  same `seed` both functions see the **same Gamma draws** for a sample.
  `dirmult_ttest()` must stay bit-identical (its tests guard this).
- Scale normals from stream `d * (n + 1) + n` (as `dirmult_ttest()`'s one
  shift), drawn in design-column order (clr/tss) or sample order (load).
- Phase B: the OLS of [da_test.md](da_test.md)'s `src/feature_lm.c`
  (precompute \eqn{A}, \eqn{h}, \eqn{(X^\top X)^{-1}} once; blocks of
  features per thread). `pt()` on the main thread; BH per covariate per
  draw needs every feature's p: sort per covariate in threads (one task
  per covariate), then running sums in draw order on the main thread.
- Memory: one n × D draw, plus P × D running sums (×4: estimate, se,
  \eqn{\ell}, u, and their adjusted sums).

### Result

`c('ecodive_dirmult_lm', 'data.frame')`, one row per feature × design
column (intercept included, as `ancombc()`; ALDEx3's summary hides it by
default): `feature, covariate, estimate` (log2), `se, p_value,
p_adjusted, reject`. Attributes: `formula`, `groups`, `design`,
`scale_model`, `scale_sd`, `load_sd`, `n_draws`, `pseudocount`,
`se_type`, `p_adjust`, `alpha`, `seed`. `plot()`: volcano / forest
(estimate ± qt·se, labelled as approximate) / ma via
[da_plots.md](da_plots.md). No `overlap` statistic: ALDEx3's
`aitchison.mean()` is commented "Private function from ALDEx2" (GPL
provenance); don't port it.

### Later: mixed models (SR-MEM)

ALDEx3's `method = "lme4" / "nlme" / "blmm"` (McGovern & Silverman 2026).
`blmm` uses TMB/C++ (`src/blmm.cpp`, part of the MIT package, but TMB and
RcppEigen are dependencies ecodive can't take). Plan a C random-intercept
REML (shared with [linda.md](linda.md) and [maaslin3.md](maaslin3.md));
scikit-bio's BSD `dirmult_lme()` (`_lme.py`) is a permissive reference for
the closed-form random-intercept fit. Separate step.


## Reference implementations

- **Port from:** ALDEx3 1.4.0 as above. Extend the existing ALDEx3 stanza
  in `inst/COPYRIGHTS` and the header of the new C file.
- **Fixtures:** `data-raw/expected/dirmult_lm.r`: ALDEx3 accepts `scale`
  as an N × nsample matrix, and its Gamma draws come from R's `rgamma()` in
  a known order (`aldex.sampler()`: samples, then features, `nsample` at
  a time). So: `set.seed(1)`, run `aldex(Y, X, scale = <logScale matrix>,
  test = ...)`; `set.seed(1)` again and regenerate the same `rgamma()`
  calls to export the draws; feed both to an internal `dirmult_lm_details(
  gamma_draws =, scale_draws =)` (as `dirmult_ttest_details()`). Cover
  `se_type` × {two groups, three groups + numeric}, and `clr.sm`,
  `tss.sm`, `sample.sm` draws captured the same way.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Title: the type label, " of " and the `DA_TESTS` name for `method`
  (`"dirmult_lm"`), for example "Volcano plot of Dirichlet-Multinomial
  Linear Model". The input is a count table, so no `display`.
- Attributes plots read: `alpha`, `log_base` (`2`, as ALDEx3),
  `conf_level` if intervals are returned; and the scale assumption that
  qualifies every estimate (`scale_model`: `"clr"`, `"tss"`, or `"load"`
  when `load` was given; `scale_sd`), which the subtitle shows.
- `plot()`: [da_plots.md](da_plots.md)'s types per covariate.

## Tests

- With injected draws: `estimate`, `se`, `p_value`, `p_adjusted` match
  ALDEx3 at 1e-10 (it solves normal equations; ecodive uses QR).
- With `scale_sd = 0`, `scale_model = 'clr'`, `~ group` (two levels) and
  the same `seed`: `estimate` for the group column equals
  `dirmult_ttest()`'s `log2_fc` (contrast reversed as needed) at 1e-12;
  same draws, same mean difference.
- `identical()` for `cpus` 1, 2, 4; dense = sparse = `margin = 2`.
- Errors: N ≤ P, `load` ≤ 0 or `NA`, `load` without `load_sd`, `load`
  and `scale_sd` both given, invalid `se_type`; rank-deficient design.
- Smoke line.


## Docs

Rd page in *Differential abundance*: the scale model (with a table of
which to use when, following ALDEx3's own guidance), units (log2), the
p-value aggregation, the relation to `dirmult_ttest()`, a γ-sensitivity
pointer to [da_sensitivity.md](da_sensitivity.md), and "Same names in
other packages": ALDEx3's `aldex()` (the original). Glossary:
`scale_model`, `load_sd`, `se_type` (new). NEWS.


## References

- McGovern, K. C., & Silverman, J. D. (2026). Scale reliant mixed effects
  models enhance microbiome data analysis. *Microbiome*, 14, 117.
  doi:10.1186/s40168-026-02377-x (4 OpenAlex, retrieved 2026-11-19).
- Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). *Genome
  Biology*, 26, 139. doi:10.1186/s13059-025-03609-3.
- Dos Santos, S. J., et al. (2026). *PLoS Computational Biology*, 22(9),
  e1014728. doi:10.1371/journal.pcbi.1014728.
- Nishijima, S., et al. (2025). Fecal microbial load is a major determinant
  of gut microbiome variation and a confounder for disease associations.
  *Cell*, 188(1), 222–236. doi:10.1016/j.cell.2024.10.022.
- Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
  standard errors in the linear regression model. *The American
  Statistician*, 54, 217–224.


## Decisions made while implementing (item 42, 2026-11-24)

Asked of the user at the start:

- **`scale_sd` defaults: `dirmult_lm()` 0.5, `dirmult_ttest()` stays 0.**
  ALDEx3's `clr.sm(gamma = 0.5)`. `dirmult_lm()` rejects on averaged
  p-values only (no min/max interval envelope), so its decisions settle as
  `n_draws` grows, which was the reason for `dirmult_ttest()`'s 0
  (`b2aaf8f`).
- **`n_draws = 1000`**: ALDEx3's stated minimum with a scale model (its
  default is 2,000).
- **Alias `aldex3`: added**, as `DA_DIGIT_ALIASES <- c(aldex3 =
  'dirmult_lm')` in `R/da_test.r`, matched by `match_da_test()`'s
  digit-aware exact step (as `'ANCOM-BC2'`). It is not in `alt_ids`, whose
  letters-only matching would make it `aldex`, `dirmult_ttest`'s alias:
  `'ALDEx3'` and `'aldex 3'` give `dirmult_lm`; `'aldex'` and `'ALDEx2'`
  still give `dirmult_ttest`. `'dirmult'` alone now matches both (an
  error, tested).
- **Random intercepts: included** (`~ g + (1 | id)`), through linda's
  `linda_random()` / `linda_formula()` and `src/feature_lmm.c`
  (`lmm_fit()`, REML at the exact optimum, Satterthwaite df), ALDEx3's
  `method = "lme4"`.

Made while implementing:

- **Signature:** `dirmult_lm(counts, formula = NULL, data = NULL, groups =
  NULL, margin = 1L, pseudocount = 0.5, scale_model = c('clr', 'tss'),
  scale_sd = 0.5, load = NULL, load_sd = NULL, n_draws = 1000L, se_type =
  c('hc3', 'hc0', 'classical'), alpha = 0.05, p_adjust = 'holm', seed = 0L,
  cpus = n_cpus())`. `p_adjust = 'holm'` per the DA convention (ALDEx3:
  `'BH'`; the docs say so). No argument shadows `da_test()` or
  `da_sensitivity()` (`test-wrapper_args.r` passes).
- **`DA_TESTS` row:** `Dirichlet-Multinomial Linear Model, dirmult_lm,
  formula, estimate, 2, FALSE, TRUE`, no `alt_ids`. The effect column is
  the plan's `estimate` (not `log2_fc`: a numeric covariate's coefficient
  is per unit, not a fold change).
- **Result:** `feature, covariate, estimate, se, dof, p_value, p_adjusted,
  reject`, one row per feature × design column, intercept included (the
  plan; ALDEx3's `summary()` drops it). `dof` (not in the plan) is `n − P`,
  or with a random intercept Satterthwaite's averaged over the draws; the
  forest plot's intervals use it. Attributes: `method`, `formula` (with
  `(1 | id)` when mixed), `groups`, `design`, `random` (mixed only),
  `scale_model` (`"clr"`, `"tss"`, `"load"`), `scale_sd` (`NA` with
  `load`), `load_sd` (per sample, with `load` only), `se_type`
  (`"classical"` when mixed), `n_draws`, `pseudocount`, `alpha`,
  `p_adjust`, `seed`, `log_base = 2`, and per feature `mean` (mean clr,
  log2, over all samples and draws: the MA plot's x) and `resid_sd` (mean
  over draws of sqrt(RSS / (n − P)), or the mixed model's σ: the effect
  plot's x). `print()` is `print_table_result`.
- **`load` / `load_sd`:** `load` by `qmp_load()` (values > 0, by name or
  position, or a column of `data`); `load_sd` one value, one per sample, or
  a column of `data`, ≥ 0, required with `load`. `scale_model` or
  `scale_sd` given explicitly with `load` is an error (`missing()`
  checks), as is `load_sd` without `load`.
- **`se_type` with a random intercept:** only `'classical'` may be given
  (an explicit `'hc0'`/`'hc3'` is an error; the default is ignored), as
  ALDEx3 ignores `test` for lme4.
- **Rank-deficient designs are an error** (as ALDEx3's `solve()` fails):
  `da_design()`'s rank warning is muffled and replaced by the error. `n ≤
  P` is `da_design()`'s error.
- **HC3 with a leverage-1 sample** (a group of one sample): the HC3
  weights are `NaN` when `1 − h < 1e-10`, so every standard error and
  p-value is `NaN`, with a warning naming `se_type = "classical"`. ALDEx3
  divides ~0 by ~0 and reports meaningless numbers.
- **C (`src/dirmult_lm.c`, `C_dirmult_lm`):** per draw, phase A (threads
  over samples) fills one n × D buffer of Y from the same PCG streams as
  `dirmult_ttest()` (`d * (n + 1) + s`, features in order); phase B
  (threads over blocks of 64 features, or 4 for mixed models) fits every
  feature: A = (X'X)⁻¹X' = R⁻¹Q₁', the leverages and diag (X'X)⁻¹ come once
  from `flm_decompose()` (`src/feature_lm.c`), then coefficients A y,
  residuals, and classical / HC0 / HC3 standard errors in sample order; or
  `lmm_fit()` per feature. The scale model's normals come from stream
  `d * (n + 1) + n`: P of them in design-column order (`clr`, `tss`; none
  when `scale_sd = 0`), or n in sample order (`load`). On the main thread:
  t, `pt()` (R's can warn, so not in workers), the doubled one-sided tails
  capped at 1 (NaN kept), running sums in draw order, and `p.adjust()`
  called through R once per draw, covariate and tail (about 20% of the run
  time for 100 × 1,000; a C version of the adjustments is possible later).
  The `clr` mode skips the log-sum-exp (the closure cancels from the clr);
  `tss`, `load` and given scales use it, as `rLogDirichlet()`.
- **Test hook:** `dirmult_lm_details(..., gamma_draws =, scale_draws =)`:
  Gamma draws (n × D × n_draws, as `dirmult_ttest_details()`) and a log2
  scale matrix (n × n_draws, as `aldex(scale = <matrix>)`). Returns each
  draw's estimate, se, df, p_lower, p_upper (P × D × n_draws), `scale` (the
  log2 scale added to each sample's log2 composition) and `shift` (Λ, P ×
  n_draws, or the load's deviates times `load_sd`, n × n_draws).
- **`dirmult_input()`** takes the method's name for its "at least 2
  features" message; `da_plot()` takes `extra` for the subtitle
  (`"scale model: clr, scale_sd = 0.5"` or `"scale model: measured
  load"`).
- **Plots:** `plot.ecodive_dirmult_lm(x, covariate, type, labels,
  n_features, conf_level, numerator, denominator, ...)`, as linda's:
  volcano, effect (`resid_sd`, 0/1 covariates), MA (`mean`), forest
  (estimate ± qt((1 + conf_level)/2, dof) se, labelled "approximate"),
  rank. Titles "<Type> plot of Dirichlet-Multinomial Linear Model".
- **`da_sensitivity()`:** works unchanged (`scale_sds` and `pseudocounts`
  set `scale_sd` / `pseudocount`; `zeros` refused, `zero_free_ok = FALSE`);
  its docs list `dirmult_lm()`.
- **Not done:** `coefficient.sm()` (priors on each coefficient's shift),
  custom scale functions, correlated load errors (`s.cov`), one-sided
  tests, the `nlme` and `blmm` engines, `aldex.effect()`'s diagnostics
  (`overlap` stays unported, GPL provenance). The draft article
  `vignettes/articles/compositional.Rmd` doesn't mention `dirmult_lm()`
  yet.

How `dirmult_lm(~ group)` agrees with `dirmult_ttest()`: same `seed`, same
Gamma draws per sample (shared streams), so with `scale_model = 'clr'`,
`scale_sd = 0` and two groups, the group column's `estimate` equals
`dirmult_ttest(contrast = c(second, first))`'s `log2_fc`, `resid_sd` its
`pooled_sd` and `mean` its `mean`, all to 1.8e-15 (test bound 1e-12, two
seeds); and each draw's classical standard error is the pooled-variance
t-test's, `pooled_sd · sqrt(1/n1 + 1/n2)` (1e-12). The p-values differ
(Welch and envelope there, OLS/HC3 and ALDEx3's tail averaging here), and
with `scale_sd > 0` the shifts differ (one per design column here, one per
contrast there).

Reference results (`data-raw/expected/dirmult_lm.r`, fixture
`tests/tinytests/expected/dirmult_lm.rds`, 160 KB; ALDEx3 1.4.0, lme4
1.1-35.1, lmerTest 3.2-1, ALDEx2 1.34.0, R 4.3.3):

- **ALDEx3 on the same draws** (15 samples × 20 features, 16 draws; ALDEx3's
  compositions regenerated from `rgamma()` after the same `set.seed()` and
  given as `gamma_draws`; checked against its `logComp` to 3.6e-15): 36
  linear-model cases, two designs (`~ grp2`; `~ grp + age`) × classical /
  HC0 / HC3 × scales (ALDEx3's own `clr.sm(gamma = 0.5)` draws given as
  `scale_draws`, with BH and holm; ecodive's own `clr`, `tss` and `load`
  scales given to ALDEx3 as its scale matrix, with BH). Largest relative
  differences: estimate 2.1e-13, se 7.2e-15, p_value 6.6e-14, p_adjusted
  5.7e-14; per draw, estimate 5.3e-12 (ALDEx3 solves the normal
  equations; on a draw, ecodive is 7e-14 from `qr.coef()` and ALDEx3
  6e-13), se 3.9e-14; `reject` identical. Test bound 1e-10.
- **Random intercept** (`method = "lme4"`, two designs with `(1 | id)`, 5
  subjects): per draw, estimates within 4.2e-5 and standard errors 2.3e-5
  of lmerTest's, Satterthwaite df 2.9e-4 (lmerTest refitted for 4 draws:
  `aldex()` doesn't return df); averaged estimate 3.1e-5, p-values 1.4e-5,
  adjusted 1.7e-5 relative; `reject` identical. The differences are lme4's
  optimizer stopping short of the REML optimum, as linda found.
- **The scale models' draws** (tests, no reference): `scale` equals
  `-colMeans(log2 comp) + X Λ` (clr), `X Λ` (tss), `log2 load + shift`
  (load) to 1e-12; Λ passes a KS test against N(0, 0.7²) (4,000 draws ×
  P), doesn't depend on the table, and is 0 with `scale_sd = 0`; draw d
  doesn't depend on `n_draws`.
- **ALDEx3 with its own 2,000 draws** (deep table, counts ~2e4, `clr.sm`
  0.5) against ecodive's 2,000: estimates within 0.0094 (the Monte Carlo
  error of the mean shift, sd 0.5/sqrt(2000) = 0.011: the differences are
  the same for every feature of a covariate), per-draw SD of the estimates
  0.996–1.039 times ALDEx3's, se within 0.09%, -log10 p correlation
  0.99998. Test bounds 0.05, ±10%, 2%, 0.95.
- **ALDEx2 1.34.0 `aldex.glm()`** (output only; 1,000 clr draws, `denom =
  'all'`) against `scale_sd = 0`, `se_type = 'classical'`, 500 draws:
  estimates correlate 0.999999 (largest difference 0.0027), standard
  errors within 0.15%, -log10 p correlate 0.999997. Test bounds 0.9999,
  0.01, 1%, 0.999.
- **Kernel against `lm.fit()`** with 19 design columns (the heap path):
  estimates and classical / HC0 / HC3 standard errors to 5e-15 (bound
  1e-10).
- **Threads and inputs:** identical for `cpus` 1, 2, 4 (and mixed 1, 4);
  dense, `Matrix`, `slam`, `margin = 2` and double input identical;
  `da_test(counts, 'dirmult_lm', ...)` and `da_test(counts, 'ALDEx3', ...)`
  identical to the direct call; `da_sensitivity()`'s runs identical to the
  direct calls.
- **Speed** (4 cores, -O2, while another build ran): 100 × 1,000, `~ g +
  age`, 1,000 draws: 4.8 s on 1 core, 3.2 s on 4 (`dirmult_ttest()`: 1.8
  s); `p_adjust = 'none'` 2.6 s. 500 × 5,000, 128 draws: 13.4 s on 1 core,
  6.8 s on 4 (`dirmult_ttest()` 5.4 s). Mixed (25 subjects), 100 draws: 2.8
  s on 4 cores against 0.3 s without the random intercept.


## Outcome (from the implementing session's final report)

**Item 42, done 2026-11-24.** Commits `6fd2b9e` (work), merged with item
47 in `699fdc2` (no conflicts); this tick.

- Tests: 28,594 at home, 409 in CRAN mode on the merged tree (28,492 and
  408 before the merge; baseline 28,140 and 406). `R CMD check
  --no-manual --ignore-vignettes`: one NOTE (installed size).
- ALDEx3 1.4.0 on the same draws, 36 linear-model cases: estimate 2.1e-13,
  se 7.2e-15, p 6.6e-14, adjusted p 5.7e-14 (per draw: estimate 5.3e-12,
  ALDEx3's normal equations); `reject` identical. Random intercepts against
  lme4 / lmerTest: 4e-5 (df 3e-4), lme4's optimizer tolerance.
- With their own draws: ALDEx3 (2,000 draws) estimates within 0.0094 (the
  scale shift's Monte Carlo error), per-draw spread 0.996–1.04; ALDEx2
  1.34.0 `aldex.glm()` (outputs only) estimates r = 0.999999, se within
  0.15%.
- `dirmult_lm(~ group)` with `scale_sd = 0` equals `dirmult_ttest()`'s
  `log2_fc`, `pooled_sd` and `mean` to 1.8e-15 (shared draws).
- Speed: 100 × 1,000, 1,000 draws: 3.2 s on 4 cores; 500 × 5,000, 128
  draws: 6.8 s.
