# Effect, MA and forest plots for differential-abundance results

**Priority:** 35. **Tier:** formula (plots over existing results). The
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

- Needs item 34's `log_base` attribute on `dirmult_ttest()` and ANCOM-BC
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
