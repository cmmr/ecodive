# `da_test()` wrapper, and `da_wilcoxon()`, `da_ttest()`, `da_lm()`

**Priority:** 36. **Tier:** derived (a registry wrapper over ecodive's own
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
