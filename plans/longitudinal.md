# Longitudinal: `first_differences()`, `first_distances()`, `volatility()`, `nmit()`

**Priority:** 56 (phases 1 and 2) and 63 (`feature_volatility()`, after `random_forest()`). **Tier:** port (q2-longitudinal, Python, **BSD-3**:
`LICENSE` "BSD 3-Clause License, Copyright (c) 2017-2026, QIIME 2
development team", GitHub `master`, checked 2026-11-20; the latest
release tag is 2026.10.0: read that tag when implementing). The volatility chart's confidence interval
follows Vega's `bootstrapCI()` (vega-statistics, **BSD-3**: "Copyright (c)
2015-2023, University of Washington Interactive Data Lab", checked
2026-11-20). Read [README.md](README.md) first.

Related plans: linear mixed models on alpha diversity (q2-longitudinal's
`linear-mixed-effects`) are planned in [alpha_test.md](alpha_test.md);
ordination over time is [ctf.md](ctf.md) and
[principal_response_curves.md](principal_response_curves.md);
`feature-volatility` and `maturity-index` need random forests
([random_forest.md](random_forest.md)) and are phase 3 here.


## What it is

q2-longitudinal (Bokulich et al. 2018) is QIIME 2's toolkit for repeated
measures: per-subject changes over time, which pooled cross-sectional
tests can't see. Demand: **244** OpenAlex citations 2022–26;
"q2-longitudinal" / "volatility analysis" **94** Europe PMC microbiome
full-text mentions, "first differences" / "first distances" **127**
(retrieved 2026-11-20). R users have no equivalent outside QIIME 2.

| q2-longitudinal action | ecodive | Phase |
|:-----------------------|:--------|:------|
| `first-differences` | `first_differences()` | 1 |
| `first-distances` | `first_distances()` | 1 |
| `volatility` | `volatility()` (+ `plot()`) | 1 |
| `nmit` | `nmit()` | 2 |
| `feature-volatility`, `plot-feature-volatility` | `feature_volatility()` | 3 (needs `random_forest()`) |
| `maturity-index` | not planned (see below) | — |
| `linear-mixed-effects` | `alpha_lm()` with `(1 \| subject)` in its formula ([alpha_test.md](alpha_test.md)) | — |
| `pairwise-differences`, `pairwise-distances` | docs recipes | — |
| `anova` | not planned (base R `aov()` / `lm()`) | — |

**NMIT: include, as phase 2.** NMIT (Zhang et al. 2017) tests whether
groups of subjects differ in their **temporal interdependence** of
features: per subject, a feature × feature correlation over that subject's
time points; the distance between two subjects is a matrix norm of the
difference of their correlation matrices; a PERMANOVA on those distances
tests group differences. Low demand (Zhang et al. 2017: 50 OpenAlex
citations, 18 in 2022–26, retrieved 2026-11-20), but it is low effort in
ecodive: the per-subject correlations reuse [cross_assoc.md](cross_assoc.md)'s
kernels (or `src/mantel.c`'s Kendall code), and the test is the existing
`permanova()`.

**Maturity index: not planned.** It trains a random-forest regressor of
age on control subjects and reports microbiota-for-age Z-scores (MAZ). With
`random_forest()` it is a short docs recipe (`random_forest()` on controls,
`predict()` on everyone, then a Z-score per age bin), so no function is
needed unless the user asks.


## Names

- `first_differences()`, `first_distances()`, `volatility()`,
  `feature_volatility()`: q2-longitudinal's action names in snake case
  (rule 5), which QIIME 2 users know. `nmit()`: the method's acronym
  (rule 1). None collides with an ecodive registry alias (check with
  README's snippet).
- **Same names in other packages:** `volatility()` in TTR (finance;
  unrelated meaning). Say so.
- Arguments: `subjects` and `states` (the glossary's names, from
  [ctf.md](ctf.md)) for q2's `individual_id_column` and `state_column`;
  `baseline`; `replicates` for q2's `replicate_handling`; `groups` for the
  volatility chart's `default_group_column`.


## Shared design

- **`subjects`, `states`, `groups`:** one value per sample or a column
  name of `data`, handled as `groups` is everywhere (`per_sample()`,
  `match_rows()`). `states` must be numeric (q2 requires it: "state_column
  must be numeric").
- **Pairs of samples** (q2's `_get_group_pairs()`): for a pair of states
  (s1, s2), each subject contributes the pair (its sample at s1, its sample
  at s2) if it has exactly one of each. Several samples at a state:
  `replicates = 'error'` (default, q2's), `'drop'` (that subject is left
  out of the pair), or `'random'` (one is picked: q2 uses Python's
  unseeded `random.choice`; ecodive uses PCG with `seed`, one stream per
  (subject, state)). A subject missing a state contributes nothing to that
  pair.
- **States in order** (q2's `_first_differences()`): the sorted unique
  states over **all** samples. Without `baseline`, the pairs of states
  are consecutive in that list: (t1, t2), (t2, t3), ... **A subject
  sampled at t1 and t3 but not t2 contributes nothing** (neither (t1, t2)
  nor (t2, t3)). This is q2's behaviour; replicate it for parity and say
  so in the docs, with the workaround (`states` per subject as visit
  numbers). With `baseline = b`, the pairs are (b, t) for every other
  state t (static differences, q2's "baseline" mode). `baseline` must be
  one of the states.
- **Determinism:** only `replicates = 'random'` and the bootstrap in
  `volatility()` draw random numbers, both on PCG with `seed`. Everything
  else is exact arithmetic on pairs; no threads are needed except in
  `nmit()`.


## `first_differences()`

```r
first_differences(x, subjects, states, data = NULL, baseline = NULL,
                  replicates = 'error', margin = 1L, norm = 'none',
                  pseudocount = NULL, seed = 0L)
```

- **`x`:** a numeric vector with one value per sample (e.g. an
  `alpha_div()` result), a column name of `data` (q2's `metric`), or a
  **count table**, whose features are all differenced at once (q2 takes
  one feature of a table at a time; ecodive returns one column per
  feature), after `norm` (`'percent'` gives q2's relative frequencies;
  `'clr'` the compositional choice).
- **Value:** \eqn{x(s_2) - x(s_1)} for each pair; a pair with a missing
  value is dropped (q2 drops `NaN` differences).
- **Result:** a data frame of class `c('ecodive_first_differences',
  'data.frame')`, one row per pair, in state-pair order then subject
  order (q2's order is set iteration order, which is arbitrary: sort, and
  say so): `sample` (the later sample's name, q2's index), `subject`,
  `state_from`, `state` (the later state), `difference` (or one column
  per feature for a count table, named by feature). Attributes `baseline`
  (when used). Printed with `print_table_result`.
- q2 then runs `linear-mixed-effects` on the differences; the docs point
  to `alpha_lm()` with a random intercept per subject
  ([alpha_test.md](alpha_test.md)) or `lm()` for that.

## `first_distances()`

```r
first_distances(dm, subjects, states, data = NULL, baseline = NULL,
                replicates = 'error', seed = 0L)
```

- `dm`: a distance matrix over samples (any ecodive `dist`, or a `dist`).
  Samples are matched by label to `data` / the per-sample vectors.
- For each pair, the distance between its two samples (from the packed
  triangle; never `as.matrix(dm)`).
- **Result:** class `ecodive_first_distances`, columns `sample`,
  `subject`, `state_from`, `state`, `distance`; `display$title`
  "First distances of <dm's title>" via `set_display()` (its plot title
  names its input).
- q2's typical uses: "how fast does each subject's community change"
  (consecutive) and "how far from baseline" (with `baseline`, e.g.
  recovery after antibiotics).

## `volatility()`

```r
volatility(x, subjects, states, groups = NULL, data = NULL,
           conf_level = 0.95, n_resamples = 1000L, seed = 0L)
```

q2's volatility visualisation is an interactive Vega chart (a control
chart). ecodive computes its numbers and draws it with `plot()`:

- **Per (group, state):** `mean`, `n`, and a bootstrap confidence interval
  of the mean: Vega's `ci0` / `ci1` aggregate ops are
  `bootstrapCI(values, 1000, 0.05)` (Vega's `TupleStore.ci()`), i.e.
  `n_resamples` resamples of the group's values with replacement (index
  `floor(u * n)`), their means sorted, and the type-7 quantiles at
  \eqn{\alpha/2} and \eqn{1 - \alpha/2} (d3's `quantile`, R's default
  `quantile(type = 7)`). Vega draws with `Math.random()`, so its interval
  is not reproducible; ecodive's uses PCG with `seed`, one stream per
  (group, state): a seed-like difference only.
- **Global control limits** (Vega spec `DAT_GLOBAL_VALS`): over all
  samples, the mean and the sample SD (n − 1, Vega's `stdev`), and the
  limits mean ± 2 SD and ± 3 SD.
- **Per subject:** the trajectory (`subject`, `state`, `value`), for the
  "spaghetti" lines.
- `x` as in `first_differences()` but a single variable (a vector or a
  column of `data`); for features, call it per feature or use
  `feature_volatility()`.
- **Result:** a list of class `ecodive_volatility`: `means` (data frame
  `group`, `state`, `mean`, `lower`, `upper`, `n`), `global` (named
  numeric: `mean`, `sd`, `lower_3sd`, `lower_2sd`, `upper_2sd`,
  `upper_3sd`), `trajectories`; attribute `conf_level`.

`plot.ecodive_volatility(x, spaghetti = TRUE, limits = TRUE, ...)`: state on
the x axis; per-group mean lines with points and CI error bars, one colour
per group (`plot_palette()`); each subject's trajectory as a thin, faint
line in its group's colour (q2's spaghetti); the global mean as a dashed
line and the ±2 / ±3 SD limits as dotted lines; legend for groups. Title
"Volatility"; the y label is the variable's name when `x` named a column.

## `nmit()` (phase 2)

```r
nmit(counts, subjects, data = NULL, method = 'kendall',
     matrix_norm = 'frobenius', margin = 1L, norm = 'percent',
     pseudocount = NULL, cpus = n_cpus())
```

q2's `_nmit()` / `_temporal_corr()` / `_temporal_distance()`:

1. Per subject, the correlation matrix of the features over that subject's
   samples (pandas `DataFrame.groupby(id).corr(method)`; Kendall's tau-b
   by default, also Pearson and Spearman); `NaN` (a feature constant within
   the subject) becomes 0 (`fillna(0)`).
   - **pandas' diagonal:** for Kendall (and other non-Pearson methods
     computed pairwise), pandas sets the diagonal to 1 even for a constant
     feature; for Pearson (and pandas' Spearman path), a constant feature's
     diagonal is `NaN`, then 0. So the diagonals can differ between
     subjects and add to their distance. **Confirm from q2's outputs** and
     replicate.
   - q2 expects relative frequencies (its input type is
     `FeatureTable[RelativeFrequency]`), hence `norm = 'percent'`.
2. Distance between subjects i and j: the matrix norm of \eqn{C_i - C_j}
   (`numpy.linalg.norm(ord = 'fro')`, or `'nuc'`, the nuclear norm);
   `matrix_norm = 'frobenius'` or `'nuclear'`. The Frobenius norm over the
   full matrix is \eqn{\sqrt{2\sum_{a<b} d_{ab}^2 + \sum_a d_{aa}^2}}:
   compute it from packed triangles plus diagonals, never D × D per pair.
   The nuclear norm needs the singular values of a D × D symmetric matrix
   (its eigenvalues' absolute values): LAPACK `dsyevd` per pair; document
   the cost.
3. **Result:** an `ecodive_dist` over **subjects** (`dist_result()`,
   `method = 'nmit'`, title "NMIT Distance"), labelled by subject id (q2
   labels each subject by its first sample's id: say so). Users then run
   `permanova(dm, groups)` with one group value per subject; the docs show
   how to collapse `data` to subjects.

Cost and threads: S subjects, D features. Phase A, one task per subject:
its D(D − 1)/2 correlations over its n_s samples (Kendall O(n_s log n_s)
per pair; n_s is small, typically 3–20, so O(n_s²) direct counting is
fine and simpler), stored packed (S × D(D − 1)/2 doubles: 50 subjects ×
500 features is 50 MB; warn above a few hundred MB and suggest filtering
features first, as q2's docs do). Phase B, one task per subject row of
the distance triangle, summing over pairs in order. `identical()` for any
`cpus`.

## `feature_volatility()` (phase 3, after `random_forest()`)

q2's `feature_volatility` pipeline: `regress_samples` (random-forest
regression of the state on the features, with RFE feature selection),
relative frequencies, then the volatility chart of the most important
features with per-feature statistics (`_summarize_feature_stats()`: mean
per state, first differences of those means, cumulative average increase
and decrease, net change, global variance, mean, median, SD, CV). ecodive:
`random_forest(counts, y = states)` ([random_forest.md](random_forest.md),
no RFE: its importances), then a table of the top `n_features` features
with those statistics and a `plot()` with the importance bars beside the
volatility chart of a chosen feature. Write its signature when
`random_forest()` exists.


## Reference implementations

- **q2-longitudinal** (BSD-3) in a QIIME 2 environment (or its functions
  imported directly: `_first_differences()`, `_nmit()` take pandas
  objects): `data-raw/expected/longitudinal.py` writes text outputs under
  `data-raw/expected/longitudinal/`, folded into
  `tests/tinytests/expected/longitudinal.rds` by
  `data-raw/expected/longitudinal.r`.
- Data: q2-longitudinal's tutorial data (ECAM subset) are large; use a
  small simulated design: 12 subjects × 5 states, 2 groups, a subject
  missing a middle state, one with a replicate, a constant feature in one
  subject, `NA` in the metric.
- Match: differences and distances exactly (same pairs; `identical()` on
  values); NMIT distances at 1e-12 for Kendall, Pearson, Spearman,
  Frobenius and nuclear; volatility means and SD at 1e-12; the bootstrap
  interval against a base-R implementation of `bootstrapCI()` using the
  same PCG draws (exported through an internal `volatility_details()`),
  and against Vega's only in distribution (coverage).


## Tests

- Fixtures above; q2's skipped-state behaviour (t1 and t3 only)
  reproduced; `baseline` mode; `replicates` 'error' (q2's message
  wording), 'drop', 'random' (same `seed` same pick).
- `first_distances()` equals `as.matrix(dm)[cbind(a, b)]` on small data.
- `nmit()`: invariant to sample order within subjects; equals a base-R
  implementation (`cor()` per subject, `norm(type = 'F')`, `svd()` for
  nuclear); `identical()` for `cpus` 1, 2, 4.
- `volatility()`: global limits equal `mean() ± k * sd()`; `plot()` leaves
  `.Random.seed` and `par()` unchanged.
- Edge cases: one state (error, q2's "state_column must contain at least
  two unique values"); non-numeric states (error); `subjects` equal to
  `states` (q2's error); a subject with one sample (contributes nothing;
  in `nmit()`, an error naming it, since a correlation needs ≥ 3 samples:
  **check** what q2 does with 1–2 samples, likely all-`NaN` → 0).
- Smoke lines for the three phase-1 functions.


## Docs

- A pkgdown section "Longitudinal" (with `ctf()` and
  `principal_response_curves()` linked).
- The q2 behaviours above (consecutive states over all subjects; ties of
  replicate handling; NMIT's labelling), under "Differences from
  q2-longitudinal"; Vega's bootstrap versus PCG.
- Recipes: paired tests between two states (q2's `pairwise-differences`:
  `first_differences(baseline = s1)` then `wilcox.test()`); MAZ with
  `random_forest()`; linear mixed models via `alpha_lm()` with
  `(1 | subject)`.
- NEWS; `baseline`, `replicates`, `matrix_norm`, `n_resamples` in
  parameter_names.md.


## References

- Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
  Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
  longitudinal and paired-sample analyses of microbiome data. *mSystems*,
  3, e00219-18. doi:10.1128/mSystems.00219-18
- Zhang, Y., Han, S. W., Cox, L. M., & Li, H. (2017). A multivariate
  distance-based analytic framework for microbial interdependence
  association test in longitudinal study. *Genetic Epidemiology*, 41,
  769–778. doi:10.1002/gepi.22065
- Satyanarayan, A., Russell, R., Hoffswell, J., & Heer, J. (2016).
  Reactive Vega: a streaming dataflow architecture for declarative
  interactive visualization. *IEEE TVCG*, 22, 659–668.
  doi:10.1109/TVCG.2015.2467091
- Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
  Bootstrap*. Chapman & Hall.
- Subramanian, S., Huq, S., Yatsunenko, T., et al. (2014). Persistent gut
  microbiota immaturity in malnourished Bangladeshi children. *Nature*,
  510, 417–421. doi:10.1038/nature13421 (the maturity index)
