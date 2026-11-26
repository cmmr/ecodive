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


## Decisions made while implementing (item 56, 2026-11-24)

Asked of the user at the start:

- **`nmit()` returns the distances only** (q2's output and this plan's
  design), not a test: NMIT's test is `permanova()` on them, with one group
  per subject; the docs show how to collapse `data` to subjects. Its
  p-values then come from `permanova()`'s own `perm_matrix()` / `perm_ss()`
  machinery, so they match any PERMANOVA of the same distances given the
  same permutations.
- **Subjects with fewer than 3 samples** (q2 computes ±1 or 0
  correlations for them): computed as q2 does, with a warning naming them
  and their sample counts. Confirmed from q2's outputs: one sample gives
  Kendall's diagonal 1 and everything else 0, Pearson and Spearman all 0;
  two samples give ±1 (0 for a tie).
- **`plot()` for `first_differences()` and `first_distances()`**: the values
  against the later state, a box per state, a faint line per subject, a
  dashed 0 line for differences. Subtitle "Consecutive states" or
  "Differences (Distances) from state b".
- **`replicates = 'random'` re-draws for every pair of states**, as q2
  does: PCG stream `((k - 1) S + (s - 1)) 2 + side` for pair of states k,
  subject s (1-based) and side 0 (earlier) or 1 (later), one draw each
  (`C_rng_dist`), the pick `floor(u m) + 1` of the cell's m samples in input
  order.

Made while implementing:

- **`NA` in `subjects`, `states` or `groups` is an error**, as `groups` is
  everywhere else in ecodive (q2 has no rule; its sorted state list breaks
  on NaN). `NA` in `x` drops the pairs (q2's rule) and is left out of
  volatility's means.
- **Subject order**: a factor's levels (those used), else sorted (radix,
  so locale-free). Rows are ordered by pair of states, then subject.
  `nmit()` labels the subjects by subject id in that order (q2: each
  subject's first sample id, in order of first appearance).
- **States come from the samples of `x` / `dm` only** (q2 reads them from
  every metadata row, even samples missing from the distance matrix).
- **`replicates = 'error'`** reports the first replicated (state,
  subject) cell, with q2's wording adapted to ecodive's argument names.
- **Correlation kernels are new code in `src/longitudinal.c`**, not
  `cross_assoc()`'s: its kernels are bound to its file-level state (an x
  table against y columns grouped by missing-value pattern) and Knight's
  merge sort for Kendall. A subject has few samples, so Kendall counts the
  signs of every sample pair directly (exact integers, so identical to
  SciPy's tau-b), and Pearson and Spearman repeat pandas' own loops
  (Welford's updates; average ranks centred on (n + 1)/2), including
  pandas' diagonal rule (Kendall 1; Pearson and Spearman 0 for a constant
  feature, where R's `cor()` gives 1).
- **The nuclear norm** uses `lz_sym_eigen()` (JAMA's tred2/tql2, already in
  `src/lanczos.c`) rather than LAPACK's `dsyevd`, which ecodive doesn't call
  from threads; it computes eigenvectors too, about twice the work of
  eigenvalues alone (D = 100, 50 subjects: 1.6 s on 1 thread, 0.7 s on 4).
- **`nmit()` also takes `clust = NULL`**, as every other function that
  returns a `dist` does; `print()` calls its units "subjects"
  (`dist_plot_name()`).
- **`volatility()`**: `subjects` may be `NULL` (no trajectories; q2's
  `individual_id_column` is optional too); `groups = NULL` is one group,
  `"all"`. `n` counts the values used (Vega's tooltip `count` also counts
  samples whose value is missing; documented). The SD of a single value is
  0, as Vega's `stdev`. Means are `sum / n` in sample order (Vega keeps a
  running mean; the difference is rounding). The x and y labels and the
  legend title (the column names of `states`, `x` and `groups`) are in
  `display`, read only by `plot()`.
- **Bootstrap streams**: one per (group, state) row of `means`, numbered
  from 0 in row order and seeded as `C_rng_dist` seeds, so that the
  internal `volatility_details()` returns the exact uniforms for the
  test's base-R transcription of `bootstrapCI()` (index `floor(u n)`,
  type-7 quantiles).
- **Vega was compared as code read, not run**: its `bootstrapCI()` was read
  from the `vega.min.js` (Vega 4.2.0) that q2-longitudinal 2026.10.0
  bundles; running that file under Node was refused by the session's
  permission classifier. The tests compare ecodive's intervals with a base-R
  transcription exactly (same uniforms) and in distribution (coverage and
  width over 400 groups, with R's `runif()`).
- **Titles**: "First Distances of <dm's title>" (title case, as ecodive's
  other titles), "NMIT Distance", "Volatility"; `first_differences()` sets
  none (its input has none), and its print header and plot say "First
  Differences".
- **A count table in `first_differences()`** is normalized with
  `normalize()` and made dense: the result has a column per feature
  anyway.
- **pkgdown**: a "Longitudinal" section of the four functions and
  `plot_longitudinal`; `ctf()` and `principal_response_curves()` don't
  exist yet, so they aren't linked.
- **Credits**: `Authors@R` gains the QIIME 2 development team (ctb, cph),
  the University of Washington Interactive Data Lab (cph, Vega) and the
  PyData Development Team (cph, pandas' correlation loops);
  `inst/COPYRIGHTS` has one stanza for `R/longitudinal.r` and
  `src/longitudinal.c` with the four BSD-3 notices (SciPy's too).
- **Fixtures**: `data-raw/expected/longitudinal.py` loads
  `q2_longitudinal/_utilities.py` from the 2026.10.0 tag (fetched through
  jsDelivr) on its own, with an empty package object in place of its
  `__init__.py` (which imports the QIIME 2 framework) and an empty
  `q2templates`, and calls `_first_differences()` and `_nmit()` directly,
  under Python 3.13, pandas 2.3.3, NumPy 2.5.5, SciPy 1.18.1, scikit-bio
  0.7.4 and statsmodels 0.14.6 (0.15 no longer has the `multipletests`
  import q2 uses; pandas is kept below 3, as QIIME 2 releases are).


## Outcome (from the implementing session's final report)

- **Commits:** b701ec5 (the four functions, fixtures, tests, docs), a
  merge of `origin/main` (items 49, 50 and 53 had landed meanwhile:
  conflicts in `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `R/print.r`,
  `inst/COPYRIGHTS`, `pkgdown/_pkgdown.yml` and `test-smoke.r`, each
  resolved by keeping both sides; `man/print_ecodive.Rd` regenerated with
  roxygen2 8.1.1 from the merged source), and this one ticking the plans.
- **Tests:** 217 in `test-longitudinal.r` and 7 smoke lines. Before the
  merge: 29,852 at home and 439 in CRAN mode, all passing (baseline at
  2e2aa88: 29,628 and 432). On the merged tree the suite was re-run before
  pushing. `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed
  size, as before).
- **Against q2-longitudinal 2026.10.0's own functions:**
  `_first_differences()` (shannon: consecutive, baseline 0 and 7; every
  feature of a relative-frequency table; distances, consecutive and
  baseline 0): the same pairs and **identical** values (largest difference
  0; features < 1e-15). `_nmit()` for Kendall, Pearson and Spearman with the
  Frobenius and nuclear norms: largest difference 1.8e-14 absolute,
  1.1e-15 relative. Errors (replicates, a baseline that is not a state, one
  state) reproduced with q2's wording adapted.
- **Vega's `bootstrapCI()`:** exact (1e-14) against a base-R transcription
  on the same PCG uniforms (`volatility_details()`); in distribution, over
  400 groups of 8 normal values, coverage within 0.05 and mean width within
  2% of the same algorithm on R's `runif()`. Vega itself was not run (Node
  on the bundled `vega.min.js` was refused by the permission classifier).
- **Determinism:** `nmit()` identical for `cpus` 1, 2, 3 and 4 (40 subjects
  × 60 features, and the fixture data, both norms), and for dense,
  `dgCMatrix`, `simple_triplet_matrix` and transposed (`margin = 2`) input;
  `first_differences()` identical across those input types.
- **Speed:** 50 subjects × 500 features × 10 samples, Frobenius: Kendall
  0.43 s on 1 thread, 0.14 s on 4; Pearson 0.35 / 0.13 s; Spearman 0.24 /
  0.09 s. Nuclear norm with 100 features: 1.6 / 0.7 s.
- **For item 63 (`feature_volatility()`):** `volatility()`'s per-(group,
  state) means and `C_volatility_ci` (Vega's intervals, any number of
  cells), `long_x()` / `long_meta()` / `long_var()` / `long_levels()` for
  `subjects`, `states` and `groups`, `plot.ecodive_volatility()` for the
  chosen feature's chart, and `first_differences(norm = 'percent')` for the
  per-state mean differences of `_summarize_feature_stats()`.


## Decisions made while implementing (item 63)

`feature_volatility()`, 2026-11-26. Asked of the user at the start (one
round):

- **Workflow: q2's 80/20 hold-out, without RFE.** q2's pipeline runs
  `regress_samples`: `train_test_split(test_size = 0.2)`, RFECV on the
  training samples, a forest on the features RFECV keeps, and the accuracy
  of the held-out predictions. ecodive keeps the hold-out (`test_size =
  0.2`) and drops RFECV (as this plan said): one forest on every feature of
  the training samples ranks them all.
- **Ranking: the forest's impurity importance** (`importance =
  'impurity'`, q2's `feature_importances_`); `'permutation'` is the
  alternative.
- **`n_trees = 100`**, q2's default (random_forest()'s is 500).
- **`norm = 'percent'`** for the forest by default (random_forest()'s; q2
  grows it on the table given, usually counts: `norm = 'none'`); the
  statistics and the chart always use relative abundances, as q2's.

Made while implementing:

- **Signature:** `feature_volatility(counts, subjects, states, groups =
  NULL, data = NULL, n_features = 100L, test_size = 0.2, importance =
  'impurity', margin = 1L, norm = 'percent', pseudocount = NULL, n_trees =
  100L, conf_level = 0.95, n_resamples = 1000L, seed = 0L, cpus = n_cpus(),
  ...)`: `volatility()`'s order (`subjects` may be missing or `NULL`);
  `...` takes only `max_features`, `min_samples_leaf`, `min_samples_split`,
  `max_depth`, `criterion`, `max_samples`, `n_folds` and `keep_forest`
  (others are an error naming the allowed ones). `n_folds` defaults to 0:
  q2's fit has no CV without RFE or tuning; the hold-out measures accuracy.
  q2's `importance_threshold` is not ported (filter `$features`).
- **The split** is scikit-learn's ShuffleSplit: `ceiling(test_size * n)`
  held out, the first of a random permutation; the rest train
  in the permutation's order (q2's order, which the bootstrap indices
  refer to). PCG: `C_rf_shuffle()` on stream 2^62 + 2^48. `test_size = 0`
  trains on every sample, and the held-out metrics are `NA`; fewer than two
  training samples is an error. Internally,
  `feature_volatility_details(numpy_random_state = )` draws numpy's legacy
  `RandomState(rs).permutation(n)` (new `C_rf_rng_test()` kind 5:
  `random_interval()` is `mt_randint()`) and scikit-learn's tree seeds, so
  the whole pipeline reproduces q2's.
- **`rf_run()` gained an internal `train` argument** (training samples in
  a given order; the others are predicted by the full-data forest and
  returned as `holdout`), rather than subsetting the table, so sparse
  tables stay sparse and every input type works.
- **Held-out metrics** use random_forest()'s definitions: `mse`,
  `r_squared` (the coefficient of determination; q2's scatter plot's
  "r-squared" is `r^2`, which the tests compare), `slope`, `intercept`,
  `r` (scipy's `linregress()` of predicted on observed, as q2), plus the
  forest's `oob_r_squared`.
- **Ranking ties** follow the table's feature order (random_forest()'s);
  features absent from every training sample keep importance 0 (q2's
  DictVectorizer drops them).
- **Statistics** (`_summarize_feature_stats()`) are over every sample,
  training and held out, as q2's (it summarizes the whole relative table);
  the mean per state is over all groups; `global_cv` is SD / mean, as q2's
  "Global CV (%)" column, which is not a percentage. Relative abundances
  are computed only for the kept features, from `C_normalize()`'s triplets
  for sparse tables (no dense copy of the table).
- **`NA` in `states`, `subjects` or `groups` is an error** (the
  longitudinal family's rule, and q2's `missing_samples = 'error'`
  default).
- **Chart numbers:** `volatility()`'s per-(group, state) computation moved
  into `long_vol_cells()` (many variables at once; `volatility()`'s results
  are unchanged, bit for bit). `C_volatility_ci()` takes an optional third
  `args` element, a stream offset in units of 2^32; feature f's cell k
  (0-based) draws from PCG stream 2^52 + f · cells + k, a stream R doubles
  can name exactly (2^60 + k cannot), so the tests rebuild the intervals
  with `C_rng_dist()`.
- **Result:** a list of class `ecodive_feature_volatility`: `features`,
  `means` (with a `feature` column), `abundances` (samples × kept features),
  `samples` (`sample`, `subject`, `group`, `state`, `set`), `predictions`
  (held-out samples), `performance`, and `forest` (the random_forest()
  result). No `global` element: the control limits are `global_mean` ±
  k · `global_sd` of `features` (the same numbers as volatility()'s
  `global`), and `plot()` forms them. Attributes `method`, `conf_level`,
  `n_resamples`, `test_size`; `display` has the title "Feature
  Volatility", `xlab` (the `states` column) and `legend` (the `groups`
  column).
- **`plot()`**: `layout()` of two panels (widths 1 : 2): grey importance
  bars of the top `n_bars` (20) features, the chosen one (`feature`, a
  name or rank, default 1) in `PLOT_HIGHLIGHT`; then
  `plot.ecodive_volatility()` of that feature (title its name, y label
  "Relative abundance"), `...` going to the chart. `par()` is restored.
- **`print()`** is its own short method: the forest and hold-out summary
  line, then the `features` table via `print_table_lines()`.
- **Credits:** the "QIIME 2 development team" `Authors@R` entry and the
  longitudinal `inst/COPYRIGHTS` stanza are extended (q2-sample-classifier
  has the same copyright line); `R/feature_volatility.r` carries the BSD-3
  text.
- **Fixtures:** `data-raw/expected/feature_volatility.{r,py}` (78 samples
  of 16 subjects at 5 states × 25 features). The Python script loads
  q2-longitudinal's `_utilities.py` and q2-sample-classifier's
  `utilities.py` / `visuals.py` (2026.10.0 tags, fetched through jsDelivr)
  with stubbed `qiime2.plugin`, `q2templates` and a `MetadataColumn`
  stand-in, and calls `_prepare_training_data()`, `_fit_estimator()`
  (`optimize_feature_selection = False`, and `True` once for information),
  `_linear_regress()`, `_importance_filtering()` and
  `_summarize_feature_stats()`, under Python 3.13.18, scikit-learn 1.9.1,
  NumPy 2.5.5, pandas 2.3.3, SciPy 1.18.1, biom-format 2.1.18 and
  statsmodels 0.14.6.


## Outcome (item 63, from the implementing session's final report)

Implemented in `0d00ddb` on 2026-11-26.

- **Files:** `R/feature_volatility.r` (`feature_volatility()`, internal
  `feature_volatility_details()`, `print()`, `plot()`), the shared
  `long_vol_cells()` in `R/longitudinal.r`, a `train` argument to
  `rf_run()`, `C_rf_rng_test()` kind 5 (numpy's legacy permutation) and a
  stream offset in `C_volatility_ci()`; `test-feature_volatility.r` (258
  tests) and two smoke lines; fixtures from
  `data-raw/expected/feature_volatility.{r,py}` (56 KB of text, a 12 KB
  `.rds`).
- **Against q2-longitudinal / q2-sample-classifier 2026.10.0** (scikit-learn
  1.9.1; random_state 3 and 7, `optimize_feature_selection = False`, the
  forest on counts with `norm = 'none'`): the train/test split identical;
  importances, their order, the held-out predictions and q2's accuracy
  (MSE, r, r², slope, intercept) exact (largest difference 3.3e-16, in r²);
  `_importance_filtering()` with `'all'` and 10 features the same features
  in the same order; `_summarize_feature_stats()` within 5e-15 relative
  (cumulative increase and decrease, which sum differences; the global
  statistics within 5e-16, the median exactly). numpy's
  `RandomState(rs).permutation(n)` identical for n up to 1,000. q2's RFE
  run (random_state 3) kept 7 features, all among the top 10 of the same
  forest without RFE.
- **Intervals:** exact (1e-14) against a base-R transcription of Vega's
  `bootstrapCI()` on the same PCG uniforms; means and counts identical to
  `volatility()`'s per feature.
- **Determinism:** `identical()` for `cpus` 1–4 (impurity and permutation
  importance), and for a transposed table (`margin = 2`); a
  `simple_triplet_matrix` gives the same abundances and statistics
  (importances correlate > 0.99: sparse and dense regression forests can
  break tied splits differently, as noted in random_forest.md).
- **Speed** (400 samples × 2,000 features, 22% zeros, 100 trees, all
  features per split): 8.9 s on 1 thread, 2.6 s on 4, for the whole
  pipeline; scikit-learn's forest fit alone on the same 320 training
  samples took 9.6 s and 2.7 s.
- **Tests:** 33,120 at home and 471 in CRAN mode (32,859 and 469 at
  70dff27). `R CMD check --no-manual --ignore-vignettes`: one NOTE
  (installed size 7.1 Mb), as before.
