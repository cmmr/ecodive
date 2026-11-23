# Effect, MA and forest plots for differential-abundance results

**Priority:** 28. **Tier:** formula (plots over existing results). The
new `dirmult_ttest()` dispersion column follows ALDEx3 1.4.0's
`aldex.effect()` (`R/aldex-summary.R`, **MIT**, verified); ALDEx2 (GPL ≥ 3)
is output-only; CoDaSeq's `codaSeq.effect()` is a cross-check only (see
below). Read [README.md](README.md) first, and
[compositional.md](compositional.md) for the shared conventions.


## What it is

ecodive draws only a **volcano plot** for `dirmult_ttest()` and
`ancombc()` / `ancombc2()` results. Gloor, Macklaim & Fernandes (2016)
argue that volcano and Bland–Altman plots "do not show what they are
colloquially believed to show", and propose the **effect plot** (ALDEx2's
"MW"): each feature's between-group difference against its within-group
dispersion, with guide lines where the two are equal (|effect| = 1).
ALDEx2's `aldex.plot(type = "MW" / "MA" / "volcano")` made the trio
standard. This plan adds four displays as a `type` argument:

| `type` | x | y | For |
|:-------|:--|:--|:----|
| `'volcano'` (default, as now) | log fold change | −log10 adjusted p | both |
| `'effect'` | within-group dispersion | log fold change; lines y = ±x | both (ANCOM-BC: two-level comparisons only) |
| `'ma'` | mean clr abundance | log fold change | both |
| `'forest'` | log fold change ± interval, one row per feature | features, sorted by effect | both |

Evidence (retrieved 2026-11-19): Gloor et al. 2016, 198 OpenAlex
citations; ALDEx2 1,507 OpenAlex; 27 Europe PMC hits name `aldex.plot` or
an ALDEx2 effect plot (an undercount); the ANCOM-BC2 vignette's LFC ± SE
bar charts flag robust taxa, which is the forest plot here. The review
ranks these plots first among all its recommendations
(`reports/Compositional methods for ecodive.md`, visualization #1).

**Correction to the review:** Crossref lists the 2016 paper's authors as
Gloor, Macklaim and **Fernandes** (not Egozcue).


## Names

- No new functions. compositional.md (2026-11-20): "Alternative displays
  of the same result are a `type` argument, as `plot.ecodive_mrpp(type =)`
  does." Types are lowercase and spelled out: `'effect'` (Gloor's own name;
  docs say it is ALDEx2's `"MW"`), `'ma'`, `'forest'`.
- Rank plots of effects belong to [log_ratio_plots.md](log_ratio_plots.md),
  not here.


## Design

### Step 1: new `dirmult_ttest()` columns

With \eqn{z} the clr values of one draw (natural logs, as `src/dirmult.c`
computes them), side means \eqn{m_1, m_2}, sums of squared deviations
\eqn{q_1, q_2} (already in `welch_worker()`), and \eqn{n_1, n_2} samples:

- **`pooled_sd`**: per draw \eqn{\sqrt{(q_1 + q_2)/(n_1 + n_2 - 2)}}, the
  pooled within-group SD, averaged over the draws (arithmetic mean of the
  per-draw SDs, as ALDEx3's `rowMeans(sqrt(pooled.var))`), divided by
  \eqn{\log 2}. Unaffected by `scale_sd` (the shift is shared within a
  side). `NaN` only when \eqn{n_1 + n_2 = 2}.
- **`mean_clr`**: per draw \eqn{(n_1 m_1 + n_2 m_2)/(n_1 + n_2)}, the mean
  clr over the samples in the contrast, averaged over draws, divided by
  \eqn{\log 2}. With two named groups, samples outside the contrast are
  not drawn and not included (documented). ALDEx2's `rab.all` is a median
  over all samples and draws; this is a mean (simpler, exact).
- **`effect`**: per draw \eqn{\Delta_d / \text{pooled sd}_d} (with the
  scale shift in \eqn{\Delta_d}), averaged over draws: ALDEx3's
  `cohens.d`. Not `log2_fc / pooled_sd` (a mean of ratios, not a ratio of
  means). **Ask the user** whether to add this column; recommended, since
  "effect > 1" filtering is how ALDEx2 users rank features.

Why not ALDEx2's `diff.win`: its documentation defines it as "the
per-feature maximum median difference between Dirichlet instances within
conditions" (random within-group pairs, then medians), and ALDEx2 is GPL,
so only its documentation and output may be used. ALDEx3's pooled SD is MIT,
deterministic, and exact. For normal data the median absolute difference of
two draws is \eqn{0.954\sigma}, so the two are close; say so in the docs.
**Do not port** ALDEx3's `overlap` / `aitchison.mean()`: ALDEx3's own
comment calls the latter a "Private function from ALDEx2" (GPL provenance).

Column order (statistics before p-values, compositional.md): `feature, t,
log2_fc, ci_lower, ci_upper, pooled_sd, mean_clr, [effect,] p_value,
p_adjusted, reject`. C: accumulate per-draw values into running sums in
draw order on the main thread, as `delta_sum` is (no shared accumulators);
`details = TRUE` adds `pooled_sd` and `mean_clr` matrices to `draws`.

### Step 2: new `ancombc()` / `ancombc2()` attributes

Per-feature values (not per covariate), so attributes, named by feature:

- **`resid_sd`**: \eqn{\sqrt{\text{RSS}_f / (n_{obs,f} - \text{rank}_f)}}
  from the fit whose residuals the variance worker already forms
  (`src/ancombc.c`, around the `r2 = (y - X beta - theta)^2` loop):
  ANCOM-BC's θ-centred initial fit; ANCOM-BC2's final unbiased fit.
- **`mean_clr`**: mean over the feature's observed samples of
  \eqn{\log(x + pc)} minus the sample's mean log over its observed
  features: the clr (natural log), or with zeros as missing (ANCOM-BC2,
  `pseudocount = 0`) the rclr. One pass in C.

The sensitivity fits carry them too.

### Plot methods

```r
plot.ecodive_dirmult_ttest(x, type = c('volcano', 'effect', 'ma', 'forest'),
                           labels = 10, n_features = 30, ...)
plot.ecodive_ancombc(x, covariate = NULL,
                     type = c('volcano', 'effect', 'ma', 'forest'),
                     labels = 10, n_features = 30, conf_level = 0.95, ...)
```

- `type` checked with `test_plot_type()` (R/plot_tests.r). `ancombc`
  keeps `covariate` second so `plot(fit, 'age')` still works;
  `dirmult_ttest`'s `type` goes second, as `plot.ecodive_mrpp()`
  (unreleased 2.4.0 API, so moving `labels` is acceptable; say so in NEWS).
- **`'effect'`**: x = `pooled_sd` (dirmult) or `resid_sd` (ANCOM-BC), from
  0; y = `log2_fc` / `log_fc`; dashed `abline(0, 1)` and `abline(0, -1)`;
  dotted `h = 0`. For ANCOM-BC, only covariates whose design column is 0/1
  (a factor level against its reference), else the error "`type = 'effect'`
  needs a two-level comparison; use 'ma' or 'forest'". (For a 0/1 column
  with no other covariates, `resid_sd` is the pooled within-group SD.)
- **`'ma'`**: x = `mean_clr`; y = the log fold change; dotted `h = 0`.
- **`'forest'`**: the `n_features` features with the smallest
  `p_adjusted` (ties in table order), sorted by effect, largest at the
  top; a point and horizontal segment per feature; dotted `v = 0`.
  Intervals: dirmult, `ci_lower`/`ci_upper` (the subtitle says they are
  the envelope of the draws' intervals); ANCOM-BC, `log_fc ± q·se` with
  \eqn{q} = `qnorm((1 + conf_level)/2)` (`ancombc()`) or `qt(...,
  dof[f])` (`ancombc2()`), unadjusted for multiplicity (say so). `NaN`
  coefficients are left out.
- Highlighting everywhere as `plot_volcano()`: `PLOT_HIGHLIGHT` for
  `reject`, `PLOT_MUTED` otherwise; with `reject_robust`, robust filled and
  non-robust hollow; the `labels` most significant rejected features named
  (scatter types); subtitle counts as now.
- Shared helpers in `R/plot_tables.r` beside `plot_volcano()` /
  `plot_ranked()`: `plot_da_scatter(x, y, ...)` (effect, MA) and
  `plot_forest(est, lo, hi, labels, reject, robust, n, args, main, xlab)`.
  The planned `da_ttest()`, `da_lm()`, `linda()` and `dirmult_lm()` plots
  ([da_test.md](da_test.md) etc.) reuse them.


## Reference implementations

- **ALDEx3 1.4.0** (CRAN, MIT, "Copyright (c) 2026 ALDEx3 authors",
  already credited for ecodive's scale model): `aldex.effect()` and
  `.aldex_pooled_var()` in `R/aldex-summary.R`. Add the function names to
  the `src/dirmult.c` header and the existing `inst/COPYRIGHTS` stanza.
- **ALDEx2** (Bioconductor, GPL ≥ 3; `DESCRIPTION` of 1.35.0 verified):
  output-only. Its `aldex.effect` / `aldex.plot` Rd pages were read for the
  column meanings above.
- **CoDaSeq 0.99.7** `codaSeq.effect()`: a different estimator (median of
  `(a - b) / max(|a - a'|, |b - b'|)` over random resamples, R's RNG). Its
  `LICENSE` file is the single line "Licensed under the MIT License" with no
  copyright holder, and the package depends on ALDEx2; **not ported**; use
  only as a sanity check.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Needs item 27's `log_base` attribute on `dirmult_ttest()` and ANCOM-BC
  results. The new columns (`pooled_sd`, `mean_clr`, `effect`) and
  attributes (`resid_sd`, `mean_clr`) are computed when the result is
  made; `plot()` only reads them.
- **Titles**: the type label, " of ", and the method's display name from
  `DA_TESTS` (by `method`), for example "Effect plot of ANCOM-BC2", with
  the covariate or contrast in the subtitle (`attr(x, 'contrast')`, or the
  `covariate` chosen). No `display` needed.
- **Guides and units**: `alpha` → the significance line (volcano) and
  highlight rule; `conf_level` → the forest plot's interval label;
  `log_base` → "log2 fold change" vs "log fold change"; `|effect| = 1`
  guides (effect plot) are fixed.
- **Highlights**: `reject`, plus the `reject_robust` / `passed_ss`
  columns when they exist (solid vs open symbols, as the existing
  ANCOM-BC plots do).
- Tests: each `type`'s title names the method; axis label follows
  `log_base`.

## Tests

- `pooled_sd`, `mean_clr`, `effect`: an R reference computed from
  `dirmult_ttest_details(gamma_draws = )` matches at 1e-13; identical with
  `scale_sd = 0` and `0.5` for the same `seed` (`pooled_sd`, `mean_clr`
  only); identical for `cpus` 1, 2, 4; dense = sparse = `margin = 2`.
- Fixture (`data-raw/expected/da_plots.r`): on a deep table (counts ~1e4,
  so draws barely vary) `pooled_sd` within 1% of ALDEx3's `pooled.SD`
  and `effect` of its `cohens.d`; Spearman correlation > 0.99 of
  `mean_clr` with ALDEx2's `rab.all` (outputs only).
- ANCOM-BC attributes: `resid_sd` equals `sigma(lm())` per feature on the
  same logs and θ for a design without missing values (base R, 1e-10);
  `mean_clr` equals `colMeans(normalize(x, 'clr', pseudocount = pc))`.
- Each `type` draws on `pdf(NULL)` and returns `x` invisibly; `'volcano'`
  unchanged from today; errors for an unknown `type`, `'effect'` on a
  numeric covariate, `n_features < 1`.
- No smoke lines needed (no new exports); update existing tests that
  index `dirmult_ttest()` columns by position.


## Docs

- `plot.ecodive_dirmult_ttest` and `plot.ecodive_ancombc` Rd: each type,
  how to read the effect plot (points beyond the lines differ by more than
  their dispersion), the ALDEx2 mapping (`"MW"` = `'effect'`, `"MA"` =
  `'ma'`), and that ANCOM-BC's forest intervals are Wald intervals.
- `dirmult_ttest()` `@return`: the new columns; `ancombc()` `@return`:
  the new attributes. NEWS under 2.4.0.


## References

- Gloor, G. B., Macklaim, J. M., & Fernandes, A. D. (2016). Displaying
  variation in large datasets: plotting a visual summary of effect sizes.
  *Journal of Computational and Graphical Statistics*, 25(3), 971–979.
  doi:10.1080/10618600.2015.1131161 (198 OpenAlex citations, retrieved
  2026-11-19).
- Fernandes, A. D., et al. (2014). Unifying the analysis of
  high-throughput sequencing datasets. *Microbiome*, 2, 15.
  doi:10.1186/2049-2618-2-15 (1,507 OpenAlex, retrieved 2026-11-19).
- Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions of
  microbiomes with covariate adjustments and repeated measures. *Nature
  Methods*, 21, 83–91. doi:10.1038/s41592-023-02092-7.
- McGovern, K. C., & Silverman, J. D. (2026). Scale reliant mixed effects
  models enhance microbiome data analysis. *Microbiome*, 14, 117.
  doi:10.1186/s40168-026-02377-x (ALDEx3).


## Decisions made while implementing (2026-11-20)

Asked of the user at the start:

- **Cohen's d: added, as `effect`** (the plan's name; the user chose it
  over ALDEx3's `cohens_d`). It is each draw's shifted difference over
  that draw's pooled sd, averaged (ALDEx3's `cohens.d`), not `log2_fc /
  pooled_sd`. The docs say it plays the part of ALDEx2's `effect`, which is
  a different estimator. It shares its name with `DA_TESTS`'s registry
  column `effect` (which names each method's log-fold-change column, here
  still `log2_fc`); the two never meet in code.
- **`mean`, not `mean_clr`, everywhere:** `dirmult_ttest()`'s column (log2
  units), `ancombc()` / `ancombc2()`'s attribute (natural logs; the rclr
  with zeros as missing values), matching item 29's `da_ttest()` column
  and `da_lm()` attribute, so that one MA-plot path reads them all.
- **Item 29's plots were extended and restyled too:** `da_ttest()` (all
  four types), `da_lm()` (as ANCOM-BC, `covariate` second, `type` third),
  `da_wilcoxon()` (volcano and MA). `da_wilcoxon()` results gained a
  `mean` column (after `median_diff`) so that its MA plot has an x axis:
  the C worker copies both sides' values together before ranking and
  takes R's `mean()` of them (`r_mean()`), so it is identical to
  `da_ttest()`'s `mean` (tested), and `C_da_tests` now returns 9 columns.

Made while implementing:

- **Columns:** `feature, t, log2_fc, ci_lower, ci_upper, pooled_sd, mean,
  effect, p_value, p_adjusted, reject`. In C, `welch_worker()` writes each
  draw's pooled sd and mean beside the Welch statistics, and the main
  thread adds them (and the effect) into running sums in draw order; the
  averages are divided by `log 2` after the draws, as `log2_fc`.
  `details = TRUE` adds `pooled_sd` and `mean` (natural logs) to `draws`.
  With one sample a side, `pooled_sd` and `effect` are `NaN` (0/0); with
  one sample on one side only, the other side's variance alone (n - 2 = 1
  degree of freedom), as ALDEx3's formula gives.
- **ANCOM-BC attributes:** `resid_sd` and `mean` follow `bias` / `n_iters`
  among the attributes. `resid_sd` = sqrt(RSS / (n_obs − rank)): RSS from
  `var_worker()`'s squared residuals (the same fit as the variances:
  ANCOM-BC's θ-centred initial fit, ANCOM-BC2's final fit), `NaN` for
  failed fits; the rank is the masked design's with missing values, else
  the global design's (the SVD's non-zero singular values, numpy's
  cutoff). `mean` is computed once in C from the logs before ANCOM-BC2's
  centring. Tested: ANCOM-BC's `resid_sd` equals `sigma(lm())` of the
  logs less θ (θ is orthogonal to the design, so these are the same
  residuals), ANCOM-BC2's equals `sigma(lm())` of the centred logs less
  `sample_fractions`, both to 1e-10; for one two-level factor it is the
  pooled within-group sd; `mean` equals `colMeans(normalize(x, 'clr',
  pseudocount = 1))` and, with zeros as missing, an R rclr.
- **Titles:** "<Volcano / Effect / MA / Forest> plot of <`DA_TESTS`
  name>" for every DA result (`da_plot_title(x, type)`), so the old
  "Dirichlet-multinomial t-test" and "ANCOM-BC2" volcano titles changed to
  "Volcano plot of Dirichlet-Multinomial t-test" and "Volcano plot of
  ANCOM-BC2". The comparison moved from the x label to the subtitle:
  "a (n = 6) vs all others (n = 12): 4 of 30 features differentially
  abundant[, k robust to the pseudocount (filled)]" or "grpb: ...". The
  count is over all features (the volcano used to count only the plotted
  ones; rejected features are always plotted, so only the denominator can
  differ, when coefficients are `NaN`).
- **Axis labels:** dirmult and ANCOM-BC from `lfc_label()`; the `da_*()`
  plots now say "Log fold change (difference in means, clr)" when
  `log_base` is finite, else "Difference in means (percent)" (item 27
  left them norm-based). Dispersion axes: "Pooled within-group standard
  deviation (log2)" (dirmult), "Residual standard deviation (natural
  log)" (ANCOM-BC), "(clr, natural log)" / "(percent)" for the `da_*()`
  ones; MA axes "Mean clr abundance (log2)", "Mean abundance (clr|rclr,
  natural log)", "Mean abundance (<norm units>)".
- **Forest plot:** the interval's description goes in the x label ("Log2
  fold change; bars: envelope of the draws' 95% intervals", "...; bars:
  unadjusted 95% Wald intervals", "...; bars: 95% confidence intervals"),
  because the subtitle was too long to fit. Bars whose ends are `NaN`
  (one sample a side, `NaN` dof) are not drawn; the point is. Light
  dotted row lines help read the names. Like the ANCOM-BC heatmap, it
  widens the left margin to the names and restores it (tested).
- **`plot.ecodive_ancombc(conf_level = 0.95)`** is validated like the
  methods' own `conf_level`; ANCOM-BC2's Wald intervals use `qt()` with
  each feature's `dof` (so `NaN` dof gives no bar).
- **Effect plot of model results:** allowed when the covariate's design
  column holds exactly the values 0 and 1 (`da_plot_binary()`), so a
  numeric 0/1 covariate also qualifies; error "`type = "effect"` needs a
  two-level comparison (a design column of 0s and 1s), and `age` isn't
  one; use "ma" or "forest"."
- **Results made before these columns** (saved objects): the plots that
  need them stop with "This result has no `pooled_sd` column, which
  `type = "effect"` plots: fit it again with this version of ecodive.";
  the volcano and forest plots still work.
- Shared helpers: `da_plot()` (dispatch, in `R/plot_da.r`, with
  `da_plot_title()`, `da_plot_contrast()`, `da_norm_units()`,
  `da_effect_label()`, `da_plot_need()`, `da_plot_binary()`), and in
  `R/plot_tables.r` `plot_volcano()` (now given its subtitle),
  `plot_da_scatter()`, `plot_forest()`, `da_plot_sub()`,
  `da_plot_marks()`, `da_plot_labels()`, `plot_pct()`. `linda()`,
  `dirmult_lm()` and later DA methods can call `da_plot()` with their
  vectors.
- **Licensing:** the `src/dirmult.c` header and `R/dirmult_ttest.r` name
  ALDEx3's `aldex.effect()` and `.aldex_pooled_var()`; the existing
  `inst/COPYRIGHTS` stanza for ALDEx3 says what follows them; the
  `Authors@R` comment reads "the scale model and effect sizes of
  dirmult_ttest()". `overlap` / `aitchison.mean()` were not ported, as the
  plan says.
- **References** (`data-raw/expected/da_plots.r`, fixture
  `tests/tinytests/expected/da_plots.rds`, 120 KB):
  - ALDEx3 1.4.0, on its own Dirichlet draws given to
    `dirmult_ttest_details(gamma_draws = 2^logComp)` (32 draws, 12
    samples × 20 features): `estimate`, `pooled.SD` and `cohens.d` match
    `log2_fc`, `pooled_sd` and `effect` to 7.8e-16 (largest relative
    difference) without a scale model and to 3.9e-15 with one (test
    bound 1e-10). For the
    scale case ALDEx3 was given the scale matrix `-colMeans(logComp) +
    x * shift`, with ecodive's own shifts for those draws (clr.sm()'s form
    with the treatment's shift fixed); `aldex()` takes a matrix `scale`,
    and both calls start from the same R seed, so they draw the same
    compositions (checked with `identical()`).
  - ALDEx3 with its own 2,000 draws on a deep table (counts ~1e4) versus
    ecodive's 256: `pooled_sd` within 0.21%, `effect` within 0.34% (test
    bounds 1% and 2%), `log2_fc` within 0.0019.
  - ALDEx2 1.34.0 (`aldex.clr(mc.samples = 1000, denom = 'all')`,
    `aldex.effect()`), outputs only, on the deep table: Spearman 0.994
    between `mean` and `rab.all` (mean − rab.all from −0.13 to 0.12);
    Pearson 0.994 between `log2_fc` and `diff.btw` (largest difference
    0.20); Spearman 0.94 between `pooled_sd` and `diff.win`, with
    `diff.win` 1.22 to 1.54 times `pooled_sd` (median 1.41): larger than
    the plan's 0.954σ argument suggested, since `diff.win` is the larger
    of the two groups' values and a median over random pairs (the docs
    now say so instead of calling them close); Pearson 0.998 between
    `effect` and ALDEx2's `effect`, ALDEx2's being smaller (median ratio
    0.74). Test bounds: 0.99, 0.99, 0.9, 0.99.
  - CoDaSeq was not run (a cross-check only, and the ALDEx2 comparison
    covers the same ground).
- **Article:** `vignettes/articles/compositional.Rmd` gained an effect plot
  of `res_dm` after the volcano plot, with a sentence on the other types,
  citing Gloor, Macklaim & Fernandes (2016) († : closed access; no free
  version found on arXiv).


## Outcome (from the implementing session's final report)

**Item 28, done 2026-11-20.** Commits `48cf7ad` (work), `a6c7a2c` (tick).

- Tests: 18,780 at home, 342 in CRAN mode.
- ALDEx3 1.4.0 `aldex.effect()` on the same draws: ≤ 3.9e-15 (with a scale
  model). ALDEx2 1.34.0 (outputs): effect r = 0.998, mean ρ = 0.994;
  `diff.win` is 1.2–1.5× `pooled_sd` (docs say so).
- ANCOM-BC `resid_sd` equals `sigma(lm())` to 1e-10.
- Plot titles are now "<Type> plot of <method>".
