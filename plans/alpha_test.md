# `alpha_test()` wrapper, and `alpha_wilcoxon()`, `alpha_kruskal()`, `alpha_ttest()`, `alpha_lm()`, `alpha_cor()`

**Priority:** 48. **Tier:** formula + derived: classical tests on ecodive's
own per-column kernels (`src/da_tests.c`, `src/feature_lm.c`,
`src/feature_lmm.c`); Dunn's test as scikit-posthocs (MIT) computes it;
q2-diversity's pairwise Kruskal-Wallis approach (BSD-3); base R, lme4,
lmerTest, emmeans, dunn.test and FSA output-only. Read
[README.md](README.md) first, then [da_test.md](da_test.md) (the same
shape for count tables), [pairwise_tests.md](pairwise_tests.md) and
[linda.md](linda.md) (the mixed model).


## What it is

The alpha counterpart of `dist_test()`: does a grouping or a covariate
explain differences in alpha diversity? Nearly every microbiome paper runs
one of these on Shannon or observed richness. Europe PMC microbiome full
texts 2022–26 (retrieved 2026-11-20): "alpha diversity" with
Wilcoxon/Mann-Whitney 15,920; with Kruskal-Wallis 12,841; with a mixed
model 2,852 (baselines: Bray-Curtis 28,581; LEfSe 19,755; MaAsLin2 2,232).

ecodive computes alpha diversity but leaves the test to the user, who
reaches for base R plus FSA or dunn.test (Dunn), lme4/lmerTest (repeated
measures) and emmeans (pairwise contrasts). The sister package rbiom does
exactly this in `adiv_stats()` with emmeans and a mixed model; these
functions should let rbiom drop emmeans for alpha statistics, and give
QIIME 2 users `alpha-group-significance` (Kruskal-Wallis plus pairwise) and
`alpha-correlation`.

Several metrics are tested at once: `x` is a vector, or a matrix with one
column per metric (as [alpha_metrics.md](alpha_metrics.md)'s
`alpha_div(counts, metric = c(...))` returns), and the result has one row
per metric (× pair, term or coefficient).


## Names

- **One function per named test, plus a wrapper with a registry,** as
  README's rule and the `da_test()` precedent: `alpha_wilcoxon()`,
  `alpha_kruskal()`, `alpha_ttest()`, `alpha_lm()` (linear model, ANOVA and
  the random-intercept mixed model), `alpha_cor()`. Five tests justify a
  wrapper (`dist_test()` has six), and rbiom's `adiv_stats(test = )`
  selector maps straight onto `alpha_test(test = )`.
- The prefix `alpha_` sorts the family with `alpha_div()` /
  `alpha_div_rarefied()` (rule 2). `alpha_kruskal`, not `alpha_kw` (rule
  4). The mixed model is `alpha_lm()` with `(1 | id)` in `formula`, as
  `linda()` does, not a separate `alpha_lmm()`: one model, and the formula
  says which.
- **Registry `ALPHA_TESTS`** in `R/alpha_test.r`, with `list_alpha_tests()`
  and `match_alpha_test()`, shaped as `DIST_TESTS` (columns `name, id,
  design, groups, alt_ids`):

  | id | name | design | groups | alt_ids |
  |:---|:-----|:-------|:-------|:--------|
  | `alpha_wilcoxon` | Wilcoxon Rank-Sum Test | groups | categorical | `wilcoxon wilcox mannwhitney mannwhitneyu ranksum` |
  | `alpha_kruskal` | Kruskal-Wallis Test | groups | categorical | `kruskal kruskalwallis kw dunn` |
  | `alpha_ttest` | Welch t-test | groups | categorical | `ttest welch studentttest` |
  | `alpha_lm` | Linear Model | formula | categorical numeric | `lm anova linearregression ols lmm lmer mixedmodel linearmixedmodel emmeans` |
  | `alpha_cor` | Correlation Test | groups | numeric | `cor correlation spearman pearson` |

- **Alias collision (ask the user):** `wilcoxon`, `ttest`, `lm`, ... and
  the display names "Wilcoxon Rank-Sum Test", "Welch t-test", "Linear
  Model" are also `DA_TESTS`'. README asks for aliases unique across
  registries. Recommended: allow the same names in `DA_TESTS` and
  `ALPHA_TESTS` (the same test on different inputs, behind different
  wrappers; each matcher looks in its own registry first), and exempt that
  pair in README's rule and in the alias-uniqueness test. Redirects: a
  `DIST_TESTS` or `COUNT_ORDS` name given to `match_alpha_test()` says "use
  dist_test()" etc.; `match_dist_test('kruskal')` says "a test of alpha
  diversity values: use alpha_test()".
- **Wrapper shadowing** (new README row): `alpha_test()` takes `x`, `test`,
  `groups`, `data`, `formula`, so a method must not take an argument named
  `test` or any prefix of `test`, `groups`, `data`, `formula` other than
  the pass-through ones (`t`, `te`, `g`, `gr`, `d`, `da`, `f`, `fo`, ...).
  `alpha_cor()`'s `method` is allowed (the selector is `test`, as in
  `dist_test()`). Extend `test-wrapper_args.r`.
- **Arguments:** reuse `groups`, `data`, `formula`, `contrast`, `pairwise`,
  `p_adjust`, `conf_level`, `var_equal`, `by` (permanova's meaning),
  `method` (mantel's: the correlation), `cpus`. New: `posthoc` (which
  pairwise procedure follows Kruskal-Wallis). The first argument is `x`
  (glossary: a first argument of more than one kind, here a vector or a
  matrix); `alpha` is taken (the significance level elsewhere).
- Same names in other packages: none known.


## Design

### Signatures

```r
alpha_test(x, test, groups = NULL, data = NULL, formula = NULL, ...)

alpha_wilcoxon(
  x, groups, data = NULL, contrast = NULL, pairwise = FALSE,
  p_adjust = 'holm', conf_level = 0.95, cpus = n_cpus() )

alpha_ttest(
  x, groups, data = NULL, contrast = NULL, var_equal = FALSE,
  pairwise = FALSE, p_adjust = 'holm', conf_level = 0.95, cpus = n_cpus() )

alpha_kruskal(
  x, groups, data = NULL, pairwise = FALSE, posthoc = 'dunn',
  p_adjust = 'holm', cpus = n_cpus() )

alpha_lm(
  x, formula = NULL, data = NULL, groups = NULL, by = 'terms',
  pairwise = FALSE, p_adjust = 'holm', conf_level = 0.95, cpus = n_cpus() )

alpha_cor(x, groups, data = NULL, method = 'spearman', conf_level = 0.95,
          cpus = n_cpus())
```

### Input (shared helper `alpha_input()`)

- `x`: a numeric vector (one value per sample, named by sample), or a
  numeric matrix or data frame with samples in rows and one column per
  metric. Column names are the metric labels (`'1'`, `'2'`, ... when
  absent); a vector is one metric, labelled by `deparse(substitute(x))`
  only when that is a plain name, else `'x'`. Attributes are ignored, so
  `alpha_div_rarefied()` output and [hill_extrapolation.md](hill_extrapolation.md)'s
  `hill_standardize()` matrix work.
- Samples are matched to `data` by name when both are named, else by
  position: `test_group_values()` / `da_design()`, with messages naming
  `x` (their `ref` argument).
- **Non-finite values** (`NA`, `NaN`, `Inf`, which `chao1()` gives when a
  sample has no doubletons) are missing: each metric's test uses the
  samples where it is finite, and one warning per metric says how many were
  left out. Samples with `NA` in `groups` or a model variable are left out of
  every metric (as `lm()`'s `na.omit`), with a warning.
- All results use the internal codes for `NaN`: the kernels already skip
  `NaN` (the rclr path of `da_tests.c`).

### `groups`, `formula`, `contrast`, `pairwise`

- As `da_test()`: group tests take `groups` (+ `data`); `alpha_lm()` takes
  `formula` (+ `data`), with `groups` naming its term of interest (the one
  `pairwise` compares), and `groups` alone meaning `~ groups`. The wrapper
  turns a one-term `formula` into `groups` for group tests, and errors for
  more terms, as `da_test()` does.
- **Two-sample tests** (`alpha_wilcoxon`, `alpha_ttest`): with two groups,
  one row per metric. `contrast = c(a, b)` picks two of several groups (`b
  = NA` pools the rest, as `da_wilcoxon()`); with more than two groups and
  neither `contrast` nor `pairwise = TRUE`, an error says to use one of them
  or `alpha_kruskal()`. `pairwise = TRUE` tests every pair of levels (level
  order, `"a-b"`, as pairwise_tests.md), each on its two groups' samples,
  and adds `p_adjusted` across the pairs **of each metric**.
- **Omnibus tests** (`alpha_kruskal`, `alpha_lm`): the overall test is the
  table; `pairwise = TRUE` adds a `pairwise` attribute (a data frame, as
  the dist tests' `pairwise` field) with `p_adjusted` across each metric's
  pairs.
- `p_adjust` follows `pairwise_method()` (R/pairwise.r): Holm by default;
  giving it without `pairwise = TRUE` is an error. Adjustment is never
  across metrics or covariates; rbiom, which adjusts its whole table, calls
  `p.adjust()` on `p_value` itself. The docs say so.

### The tests (per metric; effects are group1 minus group2)

1. **`alpha_wilcoxon()`** = `wilcox.test(a, b, conf.int = TRUE,
   conf.level)` as of R 4.6.0, through `C_da_tests` (W, ties, medians; the
   exact tie-aware distribution already in `wilcox_exact_ties()`), with
   `norm = 'none'`. Columns: `metric, group1, group2, n1, n2, w, auc,
   median_diff, location_shift, ci_lower, ci_upper, p_value` (+
   `p_adjusted`). `location_shift` is the Hodges-Lehmann estimate: the
   median of the \eqn{n_1 n_2} differences. Its interval in the exact
   no-ties branch is the order-statistic interval from `qwilcox()`
   quantiles (Hollander & Wolfe). ✱ In R's other branches (normal
   approximation; ties under R 4.6's exact conditional test), R documents
   only that the interval and estimate come from inverting the test; settle
   the rule from `wilcox.test()` outputs (as da_test.md did for ties), and
   say in the docs where ecodive's estimate is the plain median of
   differences. Differences: O(\eqn{n_1 n_2}) memory per thread; above
   \eqn{10^7}, use a selection algorithm.
2. **`alpha_ttest()`** = `t.test(a, b, var.equal, conf.level)`, from
   `C_da_tests` (means, variances by R's algorithms) and `da_two_sample()`'s
   arithmetic. Columns: `metric, group1, group2, n1, n2, mean1, mean2, diff,
   ci_lower, ci_upper, t, dof, cohen_d, p_value` (+ `p_adjusted`).
   `cohen_d` = `diff` / the pooled SD (\eqn{n_1 + n_2 - 2} denominator), as
   rbiom's `.effect.size`.
3. **`alpha_kruskal()`** = `kruskal.test()` (= scipy's
   `mstats.kruskalwallis`, which q2 uses): \eqn{H = \frac{12}{N(N+1)}\sum
   R_j^2/n_j - 3(N+1)}, divided by \eqn{1 - \sum(t^3 - t)/(N^3 - N)};
   \eqn{\chi^2_{k-1}} p. Columns: `metric, n, n_groups, h, dof,
   epsilon_squared, p_value`; `epsilon_squared` = \eqn{H (N + 1)/(N^2 -
   1)} (Tomczak & Tomczak 2014). Pairwise, `posthoc`:
   - `'dunn'` (default): Dunn (1964) on the overall ranks,
     \eqn{z = (\bar R_a - \bar R_b) / \sqrt{(N(N+1)/12 - \sum(t^3 -
     t)/(12(N-1)))(1/n_a + 1/n_b)}}, two-sided normal p (scikit-posthocs'
     `posthoc_dunn()`, FSA's `dunnTest()`; dunn.test reports one-sided
     \eqn{P(Z \ge |z|)} unless `altp = TRUE`, its docs say). Columns `metric,
     group1, group2, n1, n2, mean_rank1, mean_rank2, z, p_value, p_adjusted`.
   - `'kruskal'`: Kruskal-Wallis rerun on each pair's samples, q2-diversity's
     `alpha-group-significance` (which adjusts with BH; ecodive's default
     stays Holm). Columns `..., h, p_value, p_adjusted`.
   - **Ask the user** which is the default: Dunn is the textbook follow-up
     and what papers cite ("Kruskal-Wallis with Dunn's post hoc test");
     subset-and-rerun is q2's and ecodive's own rule for the dist tests.
4. **`alpha_lm()`**: OLS per metric through `src/feature_lm.c` (its QR
   handles each metric's own rows, `NaN` = missing). `by`, as `permanova()`:
   - `'terms'` (default): sequential F-tests, `anova(lm())`. From Q'y
     (`flm_fit()`): each term's sum of squares is the sum of its pivoted
     columns' squared effects. Columns `metric, term, df, sum_sq, mean_sq,
     f, den_df, p_value`, plus a `Residuals` row per metric as `anova()`.
   - `'margin'`: `drop1(lm(), test = 'F')` for the terms `drop.scope()`
     allows (refit without each).
   - `'overall'`: `summary(lm())`'s F-statistic.
   - `'onedf'`: each coefficient's t-test, `summary(lm())$coefficients` and
     `confint()`: `metric, covariate, estimate, se, t, dof, ci_lower,
     ci_upper, p_value` (as `da_lm()`).

   Attributes `r_squared`, `adj_r_squared` (per metric). `pairwise = TRUE`:
   every pair of `groups`' levels as a contrast of **one** fitted model (as
   `emmeans::pairs(emmeans(fit, ~ g), adjust = 'none')`: with treatment
   coding, \eqn{L b} with \eqn{L} picking the two levels' coefficients),
   SE \eqn{\sqrt{L V L'}}, residual df. Requires that `groups` is a
   categorical term with no interaction in `formula` (error otherwise; see
   later steps for `at`). **Ask the user:** rbiom refits a model per pair;
   emmeans (and this plan) use the full model, so covariates and other
   groups inform the residual variance.
5. **Mixed model:** `formula` with one `(1 | id)` (parsed by
   `linda_random()`; move it to `R/da_design.r` as a shared
   `formula_random()`), fitted per metric by `C_feature_lmm` (REML, exact
   optimum, Satterthwaite df per coefficient). `by = 'onedf'` comes free.
   `'margin'` and `pairwise` need **contrasts**: add `lmm_contrast()` to
   `src/feature_lmm.c` returning \eqn{L b}, \eqn{L V L'} and the gradient
   of \eqn{L V L'} in \eqn{(\theta, \sigma)} (the same analytic
   derivatives as the per-coefficient df). One row: Satterthwaite's df.
   Several rows (an F-test): Fai & Cornelius (1996) as lmerTest documents
   it (Kuznetsova et al. 2017): eigen-decompose \eqn{L V L'}, a
   Satterthwaite df \eqn{\nu_k} per eigen-direction, \eqn{E = \sum
   \nu_k/(\nu_k - 2)} over \eqn{\nu_k > 2}, denominator df \eqn{2E/(E - q)}
   (✱ the rule when \eqn{E \le q}: settle from lmerTest outputs).
   `by = 'margin'` = lmerTest's `drop1()`; `by = 'terms'` (lmerTest's type I)
   and `'overall'` are an error for mixed models in the first version.
   Attributes `random`, `random_sd`, `resid_sd` (per metric), as `linda()`.
   A singular fit (random variance 0) gives the fixed-effects fit, as
   `linda()`.
6. **`alpha_cor()`** = `cor.test(x, y, method, exact = FALSE, conf.level)`
   per metric, with `method` `'spearman'` (default, q2's) or `'pearson'`:
   `metric, n, estimate, statistic, dof, ci_lower, ci_upper, p_value`
   (Pearson's Fisher-z interval; Spearman's p from the t approximation,
   which is scipy's `spearmanr` and so q2's). R's exact Spearman p (AS 89)
   is a later step, written from Best & Roberts (1975).

### C, threads, determinism

- **Ranks by group:** a small kernel beside `C_da_tests` (or a new mode
  of it): per metric, average ranks by exact equality, each group's rank
  sum and size, and the ties' \eqn{\sum(t^3 - t)}. Kruskal-Wallis and Dunn
  come from these; pairwise Wilcoxon reruns `C_da_tests` on each pair's
  samples (its own ranks).
- **Threads over metrics** (and over pairs for the pairwise reruns), each
  task's sums in a fixed order: results identical for any `cpus`. With few
  metrics this is cheap anyway; the C matters for large tables (e.g.
  10,000 samples, many pairs, the Hodges-Lehmann differences).
- p-values (`pwilcox`, `pt`, `pf`, `pchisq`, `pnorm`) on the main thread, as
  `da_tests.c` does.
- No randomness: no `seed`.

### Result

`c('ecodive_alpha_<test>', 'data.frame')`, rows in metric order (then
pairs in level order, or terms in formula order). Attributes: `method`
(the function's id), `p_adjust` (when `pairwise`), `conf_level` where
intervals are given, `var_equal`, `posthoc`, `by`, `formula`, `groups`,
`pairwise` (the post-hoc table, omnibus tests), `r_squared`,
`adj_r_squared`, and the mixed-model ones above. Printing:
`print_table_result` with a `print_df_label()` header from the registry
name ("Kruskal-Wallis Test"); the `pairwise` attribute prints after the
main table with "Pairwise comparisons (p-values adjusted by Holm's
method):", as `print.ecodive_test()` does.


## Reference implementations

- **Base R** (`stats`): `wilcox.test`, `t.test`, `kruskal.test`, `lm`,
  `anova`, `drop1`, `confint`, `cor.test`, `p.adjust`. Called directly in
  the tests (as da_test.md's tests do). R ≥ 4.6.0 for tied Wilcoxon
  p-values; older R through `da_tests.c`'s own count (as `da_wilcoxon()`'s
  tests do).
- **Dunn:** scikit-posthocs 0.18.0 (PyPI; MIT, "Copyright (c) 2026 Maksim
  Terpilovskii", LICENSE read), `posthoc_dunn()` in
  `src/scikit_posthocs/_posthocs.py` (read 2026-11-20): ties from the
  overall ranks, \eqn{\sum(t^3 - t)/(12(N-1))}, two-sided p. Only the
  formula is used (Dunn 1964), so no notice is needed; cite it.
  Output-only (GPL): dunn.test 1.4.2 (GPL-2; `altp = TRUE`), FSA 0.10.1
  (GPL ≥ 2; `dunnTest(method = 'none')`).
- **q2-diversity** (github.com/qiime2/q2-diversity, dev; BSD-3, "Copyright
  (c) 2016-2026, QIIME 2 development team"; `_alpha/_visualizer.py` read
  2026-11-20): `alpha_group_significance()` runs
  `scipy.stats.mstats.kruskalwallis` over all groups and on each pair,
  then `multipletests(method = 'fdr_bh')` over the pairs;
  `alpha_correlation()` runs `scipy.stats.spearmanr` / `pearsonr`. Approach
  only, nothing ported.
- **Mixed models, contrasts (GPL, output-only):** lme4 2.1-0, lmerTest
  3.2-1 (`summary()`, `drop1()`, `contest()`, `anova(type = 2)`), emmeans
  2.0.4 (`pairs(emmeans(fit, ~ g, lmer.df = 'satterthwaite'), adjust =
  'none')`; its default for `lmer` is Kenward-Roger, so set it). All
  versions from their CRAN `DESCRIPTION` (2026-11-20).
- **rbiom** 3.1.0 (MIT; the user's own): `adiv_stats()`'s documented
  outputs are the target for what rbiom needs; nothing to port.
- **Fixtures:** `data-raw/expected/alpha_test.r` (dunn.test, FSA, lme4,
  lmerTest, emmeans) and `alpha_test.py` (scikit-posthocs; q2's calls
  replayed with scipy and statsmodels), text under
  `data-raw/expected/alpha_test/`, into `tests/tinytests/expected/alpha_test.rds`.
  Data: `ex_counts`-derived Shannon / observed / Pielou for 3–4 groups with
  ties; a repeated-measures design (15 subjects × 2–3 visits, unbalanced,
  one singular-fit metric).


## Plotting

Follows [README.md](README.md) § "Result attributes and plots". The input
`x` has no title, so these results set no `display`.

- `plot(res, metric = NULL)`: for a pairwise or two-sample result, a forest
  plot of each pair's effect (`location_shift`, `diff`, a contrast's
  `estimate`) with its interval, pairs labelled "a-b", significant ones
  (`p_adjusted`, else `p_value`, ≤ 0.05) highlighted; for Dunn, the `z`
  values. One metric per plot (default the first).
- Title: the registry name, " of " and the metric's `METRICS` display name
  when the metric label is a `METRICS` id ("Kruskal-Wallis Test of
  Shannon Diversity Index"), else the label as given.
- Reads: `conf_level` (interval label), `p_adjust` (legend), the effect and
  CI columns. No recomputation.
- Box or strip plots of the values by group need the values, which results
  don't keep; rbiom draws those. A later `plot()` could take `x` and
  `groups` again, as `plot(ord, groups, data)` does, if asked for.


## Tests (`tests/tinytests/test-alpha_test.r`)

- Each test equals base R per metric at 1e-12: Wilcoxon (exact without
  ties, n 5 vs 6; normal, 60 vs 60; ties), the Hodges-Lehmann estimate and
  exact interval; t-test (Welch and Student, `cohen_d` by hand);
  Kruskal-Wallis with ties; `lm` for every `by` with a 3-level factor, a
  numeric covariate and an interaction; `cor.test(exact = FALSE)` both
  methods.
- Dunn equals scikit-posthocs and FSA at 1e-12, dunn.test (`altp = TRUE`)
  too; `posthoc = 'kruskal'` equals q2's replayed pairwise H and p.
- Mixed model: coefficients, SEs, Satterthwaite df, `drop1()` F and df, and
  emmeans pairwise contrasts against the fixtures at about 1e-6 (as
  `linda()`: lme4 stops short of the optimum); a singular fit; unbalanced
  subjects.
- Matrix input = column-by-column vector calls (`identical()`); a data
  frame works; non-finite values are left out per metric, with the
  warning.
- `alpha_test(x, id, ...)` is `identical()` to the direct call for every
  registered test; one-term `formula` → `groups`; multi-term `formula` to a
  group test errors; aliases match; redirects from `dist_test()` and back.
- Errors: more than two groups without `contrast` or `pairwise`; `p_adjust`
  without `pairwise`; `pairwise` with an interacting or numeric `groups`;
  numeric `groups` for a categorical test and the reverse (via
  `test_group_type()`); a group with one sample (t-test `NaN`, warning, as
  `da_ttest()`); unknown `posthoc`, `by`, `method`.
- `identical()` for `cpus` 1, 2 and 4.
- Smoke: one line per export (`alpha_test`, `list_alpha_tests`,
  `match_alpha_test`, and the five tests).


## Docs

- One Rd page for the wrapper (`alpha_test`, with the test table, the
  `groups` / `formula` convention, the alias note) and one shared page
  for the five tests, as `da_wilcoxon()`'s. Explain which test for which
  design: Wilcoxon / Kruskal-Wallis for groups without covariates;
  `alpha_lm()` for covariates; `(1 | subject)` for repeated measures (the
  most common error in alpha statistics is treating repeated samples as
  independent); that alpha diversity depends on depth, so test rarefied
  (`alpha_div_rarefied()`) or coverage-standardized
  ([hill_extrapolation.md](hill_extrapolation.md)'s `hill_standardize()`)
  values, or add log depth as a covariate.
- "Differences from ...": dunn.test's one-sided default; q2's BH default;
  emmeans' Kenward-Roger default for `lmer`; rbiom's whole-table
  adjustment.
- `vignettes/stat_tests.Rmd` (or `adiv.Rmd`): a short "Testing alpha
  diversity" section. pkgdown: Hypothesis Tests. NEWS; NAMESPACE.
- parameter_names.md §1: `posthoc` (new); `x` gains the meaning "alpha
  values: a vector, or a matrix with a column per metric".


## Later steps (not in the first version)

- **Paired and one-sample tests** (rbiom's "mean compared with `mu`"):
  `subjects` (glossary, planned for `ctf()`) for the signed-rank test and
  the paired t-test, and `mu`.
- **`alternative`** (glossary) on the two-sample tests and `alpha_lm()`'s
  coefficients (rbiom's `alt`).
- **`at` and `trend`** on `alpha_lm()`, for what rbiom gets from emmeans:
  group differences at given values of a numeric covariate when it
  interacts with `groups` (`emmeans(fit, ~ g, at = )`), and slopes of a
  numeric term per group and their differences (`emtrends()`); both are
  contrasts \eqn{L b} on the same engine (and `lmm_contrast()` for mixed
  models).
- Exact Spearman p-values (Best & Roberts 1975, AS 89); the stratified (van
  Elteren) Wilcoxon, shared with `da_wilcoxon()`.
- The ranks-by-group kernel could serve a `da_kruskal()`.


## Open questions for the user

1. Shared names and aliases between `DA_TESTS` and `ALPHA_TESTS` (above).
2. Kruskal-Wallis post hoc default: Dunn or q2's rerun per pair.
3. `alpha_lm(pairwise = TRUE)`: contrasts of one model (emmeans) or a model
   per pair (rbiom today).


## References

- Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association*, 47, 583–621. doi:10.1080/01621459.1952.10483441
- Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics*, 6, 241–252. doi:10.1080/00401706.1964.10490181
- Hollander, M., Wolfe, D. A., & Chicken, E. (2014). *Nonparametric Statistical Methods* (3rd ed.). Wiley. (Hodges-Lehmann estimate and interval)
- Tomczak, M., & Tomczak, E. (2014). The need to report effect size estimates revisited. *Trends in Sport Sciences*, 21, 19–25.
- Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 67(1), 1–48. doi:10.18637/jss.v067.i01
- Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017). lmerTest package: tests in linear mixed effects models. *Journal of Statistical Software*, 82(13), 1–26. doi:10.18637/jss.v082.i13
- Fai, A. H.-T., & Cornelius, P. L. (1996). Approximate F-tests of multiple degree of freedom hypotheses in generalized least squares analyses of unbalanced split-plot experiments. *Journal of Statistical Computation and Simulation*, 54, 363–378. doi:10.1080/00949659608811740
- Best, D. J., & Roberts, D. E. (1975). Algorithm AS 89: the upper tail probabilities of Spearman's rho. *Applied Statistics*, 24, 377–379. doi:10.2307/2347111
- Bolyen, E., et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*, 37, 852–857. doi:10.1038/s41587-019-0209-9
