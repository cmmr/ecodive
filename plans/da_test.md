# `da_test()` wrapper, and `da_wilcoxon()`, `da_ttest()`, `da_lm()`

**Priority:** 29. **Tier:** derived (a registry wrapper over ecodive's own
functions) + formula (three classical per-feature tests; base R's
`wilcox.test()`, `t.test()` and `lm()` are the references, run directly in
the tests). Read [README.md](README.md) first, and
[compositional.md](compositional.md) for the shared conventions.


## What it is

1. **`da_test(counts, method, ...)`** runs any differential-abundance
   method by name, as `count_ord()` runs ordinations. compositional.md
   (2026-11-20) records that the condition it set for a DA wrapper is met:
   `ancom()`, `ancombc()`, `ancombc2()` and `dirmult_ttest()` exist, and
   `linda()`, `maaslin3()`, `fastancom()`, `dirmult_lm()`, `rademu()` are
   planned. The wrapper is also what [da_sensitivity.md](da_sensitivity.md)
   loops over.
2. **Three plain per-feature tests** on a normalized table, with p-value
   adjustment. Wirbel et al. (2024): "Only classic statistical methods
   (linear models, the Wilcoxon test, t-test), limma, and fastANCOM
   properly control false discoveries at relatively high sensitivity"; with
   confounders, adjusted linear models mitigate them. Pelto et al. (2025, 69
   gut datasets): the best consistency with sensitivity came from
   "ORM/Wilcoxon, MaAsLin2/t-test, and LogR". Both used TSS (or log-TSS)
   rather than clr; see **Norm** below. Citations (retrieved 2026-11-19):
   Wirbel 58 Crossref / 50 Europe PMC; Pelto 35 Crossref.


## Names

- `da_` is a **prefix** for differential-abundance helpers (decided in
  compositional.md, 2026-11-20); rule 3 concerns names that *end* in
  `da`, so there is no clash. `da_test` matches `dist_test` /
  `count_ord`; `da_wilcoxon`, `da_ttest`, `da_lm` say which test
  (rule 4: spelled out, not `da_wt`). `da_ttest()` (Welch on one
  normalized table) and `dirmult_ttest()` (Welch on Monte Carlo draws)
  each say how they differ from the other.
- Registry **`DA_TESTS`** in `R/da_test.r`, with `list_da_tests()` /
  `match_da_test()`, the same shape as `COUNT_ORDS` / `list_count_ords()` /
  `match_count_ord()` (columns `name, id, alt_ids`, plus the ones below;
  `id` is the function name). Add `'DA_TESTS'` to README's alias-check loop.

| id | name (also matched) | design | effect | base | zero_free_ok | extra aliases (letters only) |
|:---|:-----|:-------|:-------|:-----|:-------------|:-----------------------|
| `ancom` | Analysis of Composition of Microbiomes | groups | — | — | TRUE | |
| `ancombc` | ANCOM-BC | formula | `log_fc` | e | TRUE | |
| `ancombc2` | ANCOM-BC2 | formula | `log_fc` | e | TRUE | |
| `dirmult_ttest` | Dirichlet-Multinomial t-test | groups | `log2_fc` | 2 | FALSE | `aldex aldex2` (**ask**: rule 6 bars the *function* name; an alias may still mislead) |
| `da_wilcoxon` | Wilcoxon Rank-Sum Test | groups | `median_diff` | norm | TRUE | `wilcoxon wilcox mannwhitney ranksum` |
| `da_ttest` | Welch t-test | groups | `diff` | norm | TRUE | `ttest welch` |
| `da_lm` | Linear Model | formula | `estimate` | norm | TRUE | `lm` |
| later: `linda`, `maaslin3`, `fastancom`, `dirmult_lm`, `rademu` | | formula | per plan | | per plan | `maaslin`, `aldex3`, `emu` |

(`base` "norm": natural log for clr/rclr, none otherwise.) Check every
alias against all registries with README's snippet: e.g. `ttest` must not
also prefix-match `dirmult_ttest` (its letters, `dirmultttest`, don't
start with `ttest`, so it doesn't).

`design` drives the wrapper; `effect`/`base` and `zero_free_ok` (may the
input be a `replace_zeros()` table?) drive `da_sensitivity()`.
`match_da_test()` redirects distance tests ("`permanova` is a distance
test: use dist_test()") as `match_count_ord()` does.


## Design

### One convention for `groups` and `formula`

Today `ancom()` / `dirmult_ttest()` take `groups` (+ `data`), and
`ancombc()` / `ancombc2()` take `formula` + `data`, with `groups` naming
the post-hoc grouping (a term of `formula` with ≥ 3 levels). Proposal
(**a decision for the user**, it touches `ancombc()`):

- **Group methods** (`design = 'groups'`): `f(counts, groups, data = NULL,
  ...)`, one categorical variable, as now.
- **Model methods** (`design = 'formula'`): `f(counts, formula = NULL,
  data = NULL, groups = NULL, ...)`. `groups` names **the term of
  interest**: the one plotted by default, compared post hoc, and
  summarised by `da_sensitivity()`. `groups` alone means `formula = ~
  groups` (a column name of `data`, or a per-sample vector, which becomes a
  `data` column `group`), as `permanova()`'s "the same as `formula = ~
  groups`". For `ancombc()`: `formula` gains the default `NULL`, and the
  "at least three levels" check moves from the fit to the post-hoc tests
  that need it (`ancombc_global()` etc. fail with that message instead).
  Other covariates: the first non-intercept term is the default term of
  interest when `groups` is `NULL`.
- **The wrapper** `da_test(counts, method, groups, data = NULL, formula =
  NULL, ...)` passes `groups`, `data` and `formula` through on purpose (as
  `dist_test()` passes `groups`, `data`): to a group method, `formula = ~
  status` (one term) becomes `groups = 'status'`, more terms are an error
  ("`ancom` compares one grouping; use a model method such as `da_lm` or
  `ancombc2` for `formula`"); to a model method, both go through.
  `quote(counts)` as in `count_ord()`. The result is the method's own,
  unchanged.
- **README shadowing table**, new row: `da_test()` | `counts`, `method`,
  `groups`, `data`, `formula` | an argument named `method`, or any prefix
  of `method`, `counts`, `groups`, `data`, `formula` other than the
  pass-through `groups`, `data`, `formula` (`m`, `c`, `co`, `g`, `gr`,
  `d`, `da`, `f`, `fo`, ...). Extend `test-wrapper_args.r` with
  `shadowed(da_test, list_da_tests('func'), passed = c('groups', 'data',
  'formula'))` and the formals check. Existing DA functions pass.

### The three tests

```r
da_wilcoxon(counts, groups, data = NULL, contrast = NULL, margin = 1L,
            norm = 'clr', pseudocount = NULL, alpha = 0.05,
            p_adjust = 'BH', cpus = n_cpus())
da_ttest(counts, groups, data = NULL, contrast = NULL, margin = 1L,
         norm = 'clr', pseudocount = NULL, var_equal = FALSE,
         conf_level = 0.95, alpha = 0.05, p_adjust = 'BH', cpus = n_cpus())
da_lm(counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
      norm = 'clr', pseudocount = NULL, conf_level = 0.95, alpha = 0.05,
      p_adjust = 'BH', cpus = n_cpus())
```

- **`contrast`**: as `dirmult_ttest()` (reuse `test_group_values()` and
  `dirmult_contrast()`): `c(treatment, reference)`, `NA` pools the rest,
  `NULL` = first group against the rest.
- **`norm`**: any `validate_norm()` choice; default `'clr'` with the
  shared `pseudocount` rule (`log_pseudocount()`), computed per sample in
  C without a dense copy (as `normalize()`). `'rclr'`: zeros are
  **missing** (left out of that feature's test, not ranked as 0).
  `'percent'` is TSS: Pelto's Wilcoxon configuration (ranks of TSS) and
  Wirbel's. Wirbel found clr/rclr *lowered* the Wilcoxon test's accuracy
  and log-TSS *raised* the t-test's and LM's. **Ask the user**: keep
  `'clr'` as the default (compositional reading, the plan's brief), and
  whether to add a `'log'` norm (log of proportions with the pseudocount
  rule) to `normalize()` so log-TSS is available (it would belong with
  [normalization_factors.md](normalization_factors.md)).
- **`p_adjust`**: over features, `NaN` p-values left out, as
  `ancombc()`; `da_lm()` within each covariate. Default `'BH'` (the
  benchmarks' setting) though ecodive's other DA functions default to
  `'holm'`: **ask the user**.
- **`da_wilcoxon()`** = `stats::wilcox.test(x_trt, x_ref, exact = NULL,
  correct = TRUE)` per feature: \eqn{W = R_1 - n_1(n_1+1)/2} (rank sum of
  the treatment side, average ranks for ties by exact equality); exact
  `pwilcox()` p-value when \eqn{n_1, n_2 < 50} and no ties, else the normal
  approximation with tie correction and continuity correction
  \eqn{0.5\,\mathrm{sign}(z)}, two-sided. Columns: `feature, w, auc`
  (\eqn{W/(n_1 n_2)}, the probability of superiority), `median_diff`
  (treatment median minus reference median, in `norm` units), `p_value,
  p_adjusted, reject`. **Later step:** `strata` (glossary name) for a
  stratified (van Elteren) test, Wirbel's adjusted Wilcoxon; coin (GPL-2)
  is output-only; write it from van Elteren (1960).
- **`da_ttest()`** = `stats::t.test(x_trt, x_ref, var.equal = var_equal,
  conf.level = conf_level)`. Columns: `feature, diff` (mean difference, a
  natural-log fold change against each sample's geometric mean when `norm`
  is clr), `ci_lower, ci_upper, t, dof, pooled_sd, mean, p_value,
  p_adjusted, reject` (`pooled_sd`, `mean` feed the effect / MA plots of
  [da_plots.md](da_plots.md)).
- **`da_lm()`** = `summary(lm(y ~ ., ...))$coefficients` per feature,
  every design column including `(Intercept)` (as `ancombc()`). Design from
  `ancombc_design()`, renamed `da_design()` and moved to `R/da_design.r`
  (treatment contrasts, radix-sorted levels, the same checks), shared by
  every model method. Columns: `feature, covariate, estimate, se, t,
  ci_lower, ci_upper, p_value, p_adjusted, reject`; attributes `formula`,
  `groups`, `norm`, `pseudocount`, `dof`, `design`, `mean_clr`,
  `resid_sd` (per feature, for the plots). A factor with several levels
  gets one t-test per level, not Wirbel's `anova()` F; say so (for one
  degree of freedom they are the same).

### C

- `src/da_tests.c`: one task per block of features; each normalizes its
  samples' values for the feature from the sparse table (the
  `ecomatrix_t` by feature: transpose once), then the statistic.
  Statistics in threads, `pwilcox()` / `pt()` / `qt()` on the main
  thread (as `ancom.c`: Rmath's `pwilcox()` keeps static state).
- **A shared per-feature OLS kernel**, `src/feature_lm.c`: given \eqn{X}
  (n × P), \eqn{A = X^+}, leverages \eqn{h}, and a block of responses,
  return coefficients, RSS, classical / HC0 / HC3 standard errors, with a
  masked-design path for missing values (from `src/ancombc.c`'s). `da_lm()`
  is its first user; [linda.md](linda.md), [fastancom.md](fastancom.md)
  and [dirmult_lm.md](dirmult_lm.md) reuse it. QR, not normal equations
  (parity with `lm()` at 1e-10).

### Result and plots

Classes `c('ecodive_da_wilcoxon', 'data.frame')` etc., the shared
per-feature structure. `plot()` methods with `type = c('volcano', ...)`
from [da_plots.md](da_plots.md): `da_wilcoxon` volcano / ma; `da_ttest`
volcano / effect / ma / forest; `da_lm` as `ancombc`.


## Reference implementations

Base R only: `stats::wilcox.test`, `t.test`, `lm`, `p.adjust`. No porting,
no license issue. The registry and wrapper copy `R/count_ord.r`'s
structure.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Every per-feature result has a top-level `method` attribute (its
  `DA_TESTS` id). `da_method_name()` reads the registry's `name` column
  for it (as `ord_method_name()` reads the ordination registries); that is
  where plots get the display name. `da_test()` adds nothing: the result
  is the method's own, so `da_test(x, 'ancombc2')` and `ancombc2(x)` give
  identical results, attributes included (test with `identical()`).
- `da_wilcoxon()`, `da_ttest()`, `da_lm()`: attributes `method`,
  `contrast`, `alpha` (behind `reject`), `conf_level` (behind the
  intervals; `da_ttest()`, `da_lm()`) and `log_base` (`exp(1)`:
  differences of natural-log transformed values).
- **Plots**: each new class gets a `plot()` method sharing
  [da_plots.md](da_plots.md)'s types (volcano default). The title is the
  type label, " of " and `da_method_name(attr(x, 'method'))`, for example
  "Volcano plot of Welch t-test"; the inputs are count tables with no
  title, so these results set no `display`.

## Tests

- Each test equals base R per feature at 1e-12 (statistics, p, CI) for
  `norm` in `'clr'`, `'percent'`, `'rclr'`; ties; exact and normal
  Wilcoxon branches (n 5 vs 6 without ties; 60 vs 60; tied data).
- `da_lm()` equals `summary(lm())` at 1e-10, including factors with 3
  levels, interactions, numeric covariates; rank-deficient design warns.
- Invariant to multiplying a sample by a constant (clr, percent).
- `identical()` for `cpus` 1, 2, 4; dense = `Matrix` = `slam` =
  `margin = 2`.
- Wrapper: `da_test(x, id, ...)` is `identical()` to the direct call for
  every registered method; `groups` → `formula` and one-term `formula` →
  `groups` translation; multi-term `formula` to a group method errors;
  aliases match and are unique across all registries; `list_da_tests()`
  shapes as `list_count_ords()`; `ancombc(counts, data = meta, groups =
  'status')` equals `ancombc(counts, ~ status, meta)`.
- Edge cases: a side with one sample (t-test `NaN`, warning), constant
  feature (`NaN` p, not rejected), all-zero feature with rclr, NA in
  `groups`, invalid `norm`/`p_adjust`/`alpha`.
- One smoke line per new export (`da_test`, `list_da_tests`,
  `match_da_test`, `da_wilcoxon`, `da_ttest`, `da_lm`).


## Docs

- Rd pages in "Compositional Data" → *Differential abundance*:
  `da_test` (with the method table, the `groups`/`formula` convention, and
  "Same names in other packages": none known), `da_wilcoxon`, `da_ttest`,
  `da_lm` (one page or three). Say which benchmarks endorse which test and
  with which normalization, and that these tests compare relative (clr:
  geometric-mean-relative) abundances, not absolute ones.
- `parameter_names.md` §1: `var_equal` (new; `t.test`'s `var.equal`).
- NEWS; pkgdown reference group.


## References

- Wirbel, J., Essex, M., Forslund, S. K., & Zeller, G. (2024). A realistic
  benchmark for differential abundance testing and confounder adjustment
  in human microbiome studies. *Genome Biology*, 25, 247.
  doi:10.1186/s13059-024-03390-9 (58 Crossref, retrieved 2026-11-19).
- Pelto, J., Auranen, K., Kujala, J. V., & Lahti, L. (2025). Elementary
  methods provide more replicable results in microbial differential
  abundance analysis. *Briefings in Bioinformatics*, 26(2), bbaf130.
  doi:10.1093/bib/bbaf130 (35 Crossref, retrieved 2026-11-19).
- Nearing, J. T., et al. (2022). Microbiome differential abundance methods
  produce different results across 38 datasets. *Nature Communications*,
  13, 342. doi:10.1038/s41467-022-28034-z.
- van Elteren, P. H. (1960). On the combination of independent two-sample
  tests of Wilcoxon. *Bulletin of the International Statistical
  Institute*, 37, 351–361.


## Decisions made while implementing (2026-11-20)

Asked of the user at the start:

- **`groups` / `formula` convention: adopted**, including for `ancombc()` /
  `ancombc2()`. `formula` and `data` default to `NULL`; `groups` alone is
  `formula = ~ groups` (a column name of `data`, or one value per sample,
  which becomes the variable `group`); with `formula`, `groups` must name
  one of its terms (a column of `data`). The checks that it is categorical
  with at least three levels moved to the post-hoc tests
  (`posthoc_setup()` already made them). `attr(fit, 'groups')` stays what
  was given (`NULL` without it), so `ancombc_global(fit)` without a
  `groups` still asks for one; the "first non-intercept term" default is
  applied where a term is needed (`plot.ecodive_da_lm()`, and
  `da_sensitivity()` should do the same), not stored.
- **`norm`: `'clr'` default, and a new `'log'` norm** in `normalize()` and
  everywhere `norm` is offered: \eqn{\log((x + pc) / \sum (x + pc))}, the
  closure *after* the pseudocount (so it is the clr shifted by one constant
  per sample; the zeros' value goes in `clr_vec` like the clr's, so every
  C consumer of `norm` handles it unchanged). Same pseudocount rule as the
  clr (`log_pseudocount()`, "cannot perform log normalization"). `simper()`
  refuses it as it refuses clr/rclr (negative values).
- **`p_adjust` default `'holm'`**, as the other DA methods (not the
  benchmarks' `'BH'`); the docs point to `'BH'`.
- **`aldex` / `aldex2` aliases for `dirmult_ttest`: added.** The docs say
  it is an ALDEx2-style method, not ALDEx2.

Made while implementing:

- **Registry** `DA_TESTS` in `R/da_test.r` with columns `name, id, design,
  effect, base, zero_free_ok, alt_ids` (`base` is a string: `'e'`, `'2'`,
  `'norm'`, or `NA` for ancom). Aliases: ancom `analysisofcompositionsofmicrobiomes`;
  ancombc `ancombiascorrected`, the long name; dirmult_ttest
  `dirichletmultinomialwelchttest aldex aldex2`; da_wilcoxon `wilcoxon wilcox
  wilcoxontest mannwhitney mannwhitneyu ranksum`; da_ttest `ttest welch
  studentttest`; da_lm `lm linearregression ols`. None is a name or alias in
  any other registry (tested).
- **ANCOM-BC2 vs ANCOM-BC:** letters-only matching makes both `ancombc`, so
  `match_da_test()` first tries an exact match keeping digits (`'ANCOM-BC2'`
  is ancombc2, `'ANCOM-BC'` ancombc); `exact_names('DA_TESTS')` therefore
  has `ancombc` twice (both ids), which the redirect checks don't mind.
- **Redirects both ways:** `match_da_test()` sends distance tests to
  `dist_test()` and ordinations to `count_ord()` / `dist_ord()`;
  `match_count_ord()`, `match_dist_ord()` and `match_dist_test()` send DA
  names to `da_test()` (so `match_dist_test('t-test')` now says "use
  da_test()"; its test was updated).
- **`da_test()` signature** `(counts, method, groups = NULL, data = NULL,
  formula = NULL, ...)`; `seed` is dropped for methods without one (as
  `count_ord()`); a group method given `formula` takes its single variable
  (LHS variables count as terms; `~ log(age)` is an error), and given both
  `groups` and `formula` is an error.
- **`method` attribute** added to `ancom()` and `dirmult_ttest()` results
  (their docs say so); `da_method_name()` reads the registry.
- **`da_design()`** in `R/da_design.r` (renamed `ancombc_design()`, takes
  `what` for its "needs more samples than design columns" message). The
  formula built from `groups` is `as.formula(call('~', as.name(groups)))`,
  so it prints `~status`.
- **Inputs:** any input type; finite; non-negative except with `norm =
  'none'` (an already transformed table). At least one feature.
- **Results:** columns as the plan, with `mean` (not `mean_clr`) for
  `da_ttest()`'s column and `da_lm()`'s attribute, since `norm` needn't be
  clr. Attributes: `method`, `contrast`, `n` (two-sample), `norm` (its
  name), `pseudocount`, `var_equal` + `conf_level` (t-test), `alpha`,
  `p_adjust`, `log_base` (`exp(1)` for clr/rclr/log, `NA` otherwise);
  `da_lm()` also `formula`, `groups`, `design`, and per feature `dof`,
  `mean`, `resid_sd`. Missing values are `NaN` throughout.
- **Base R's errors become `NaN`, never rejected:** a side with no values
  (rclr), Welch with one value on a side (with the warning `dirmult_ttest()`
  gives, when a whole side of the contrast has one sample), t.test()'s
  "essentially constant" rule, `dof <= 0` in `da_lm()`.
- **Wilcoxon:** C returns W, the ties' \eqn{\sum t^3 - t}, and both medians
  (R's `median()` algorithm); the exact / normal branch, continuity
  correction and `pwilcox()` run in R on the main thread. `auc` = W / (n1
  n2) with the feature's own n (rclr).
- **Wilcoxon with ties follows R 4.6.0 (fixed after the first push, when
  the suite failed on the user's R 4.6):** R 4.6.0's `wilcox.test()`
  (NEWS: "can now perform exact (conditional) inference in case of ties",
  Streitberg-Rohmel shift algorithm) gives an exact p-value whenever both
  sides have under 50 values, ties or not; R 4.3.3 used the normal
  approximation with ties. ecodive now gives R 4.6's answer on every R
  version: `wilcox_exact_ties()` in `src/da_tests.c` counts every n1-subset
  of the doubled average ranks by its sum (shifted and divided by their
  gcd; the smaller side; in each worker thread, about 4 MB each), and the
  two-sided p-value is `min(1, 2 * min(P(W <= w), P(W >= w)))`. That rule
  was chosen by comparing three candidates with R 4.6.2's output (built
  from source in the scratchpad; R's own code wasn't read, being GPL): it
  agreed to 4e-15 on 399 random tied cases, the others didn't. All values
  tied gives 1, as R 4.6. The tests compare with `wilcox.test()` directly
  on R >= 4.6.0, and with an independent R count of the subsets
  (`exact_ties_p()`) on older R; a brute-force `combn()` check covers the
  distribution itself. Worst case (49 vs 49, many ties) is about 3 ms per
  feature on one thread.
- **R 4.6.2's `pwilcox()` is wrong when vectorised over different sizes**
  (e.g. `pwilcox(c(30, 40, 25, 50), c(11, 10, 9, 12), c(10, 11, 11, 10))`
  gives 0.0375 for the third, 0.0337 called alone; R 4.3.3 gives 0.0337
  either way). rclr's missing values give features different sizes, so
  `da_pwilcox()` calls it once per (m, n) pair, with scalar sizes (checked
  identical to one-at-a-time calls on R 4.6.2). Worth reporting to R core.
- **t-test:** C returns each side's n, mean and variance by R's own `mean()`
  / `var()` algorithms (long double sums with the refinement step), and the
  pooled mean; the rest is `t.test.default()`'s arithmetic in R.
- **`src/feature_lm.c`:** a Householder QR with dqrdc2's limited pivoting
  (tol 1e-7, columns whose remaining norm falls below tol × their own move
  to the end and are aliased), written from the algorithm's description;
  no BLAS, so worker threads can call it. API: `flm_init()`,
  `flm_decompose()` (rows subset = the masked design), `flm_fit()` (Q'y,
  coefficients in design order with `NaN` where aliased, RSS as the sum of
  the last n − rank squares of Q'y), `flm_se()` (classical, HC0, HC3; the
  HC ones need `q1`). `C_feature_lm` exposes it to the tests, which check it
  against `lm()` and the sandwich formulas. The complete design is
  decomposed once; features with missing values (rclr) get their own QR.
- **`da_lm()` per factor level:** one t-test per design column, as the
  plan; `(Intercept)` rows included.
- **Plots:** volcano only (item 28 adds the other types):
  `plot.ecodive_da_wilcoxon()` / `_da_ttest()` (x: difference in medians or
  means, with the norm in the label) and `plot.ecodive_da_lm(x, covariate =
  NULL)` (default: the first column of `groups`'s term, else the first
  non-intercept column), all titled "Volcano plot of <registry name>"
  (`da_plot_title()`), in `R/plot_da.r`. Printing: `print_table_result`,
  header from `da_method_name()`.
- **Not done (later steps):** the stratified (van Elteren) Wilcoxon; the
  effect / MA / forest plots (item 28); `log_base` on the *existing* DA
  results (item 27; the new ones have it).

What later items can rely on:

- **Item 34 (`da_sensitivity()`):** `list_da_tests()` gives `design`,
  `effect` (column name), `base` and `zero_free_ok` per method;
  `da_test(counts, id, groups, data, formula, ...)` runs any of them with
  one call shape; every result has `method`, `alpha`, `reject` (and
  `reject_robust` / `passed_ss` for ancombc's sensitivity analysis); model
  methods' term of interest is `attr(res, 'groups')`, else the first
  non-intercept term. Effects: `log_fc` (ancombc*), `log2_fc`
  (dirmult_ttest), `median_diff` / `diff` / `estimate` (da_*, in `norm`
  units, `log_base` attribute). ancom has no effect column.
- **Item 27 (`log_base`):** `da_wilcoxon()`, `da_ttest()`, `da_lm()` already
  set it; `ancombc()`, `ancombc2()` and `dirmult_ttest()` still need it.
- **Item 28 (DA plots):** `da_ttest()` has `pooled_sd` and `mean` columns;
  `da_lm()` has `resid_sd`, `mean` and `dof` attributes and `ci_lower` /
  `ci_upper` columns; `da_plot_title()` / `da_method_name()` give titles.
- **linda / fastancom / dirmult_lm:** `da_design()`, `da_table_init()` /
  `da_table_column()` (normalized table by feature, any norm, rclr zeros as
  `NaN`) and the `flm_*` kernel (with HC0 / HC3).
