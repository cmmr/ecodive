# `rarefaction_curve()`: alpha diversity against sequencing depth

**Priority:** 38. **Tier:** derived: resampling loops over ecodive's own
`rarefy()` and `alpha_div()`, as `alpha_div_rarefied()` does, plus two
published closed forms (expected richness, Hurlbert 1971 and Heck et al.
1975; expected coverage, Chao & Jost 2012). The depth grid and summaries
follow QIIME 2's `q2-diversity alpha-rarefaction` (BSD-3, read); vegan's
`rarefy()` / `rarecurve()` (GPL-2) are output-only. Read
[README.md](README.md) first, and [rarefied_diversity.md](rarefied_diversity.md),
whose machinery this reuses.


## What it is

A rarefaction curve plots a sample's alpha diversity against the depth it
is subsampled to. It shows whether sequencing was deep enough (the curve
levels off) and which depth to rarefy to. It is one of the most common
figures in microbiome papers: "rarefaction curve(s)" appears in 8,700
microbiome full texts 2022–26, "alpha rarefaction" (QIIME 2's name) in
624 (Europe PMC, retrieved 2026-11-20; baseline: Bray-Curtis 28,581).
"Coverage-based rarefaction" / "sample coverage": 375.

Two ways to compute a curve, both offered:

- **Expected richness, exactly** (Hurlbert 1971; Heck et al. 1975): the
  mean number of features in a random subsample of \eqn{d} reads without
  replacement, \eqn{E[S_d] = \sum_i \left[1 - \binom{N - x_i}{d} /
  \binom{N}{d}\right]}, with Heck's variance. This is what vegan's
  `rarefy()` and `rarecurve()` draw. No randomness, no iterations.
- **Resampling** for any alpha metric (QIIME 2's approach): rarefy `times`
  times at each depth and summarise the metric over iterations.

Each curve also carries the **expected sample coverage** of a depth-\eqn{d}
subsample (Chao & Jost 2012), \eqn{C_d = 1 - \sum_i \frac{x_i}{N}
\binom{N - x_i}{d} / \binom{N - 1}{d}}: the sample-completeness curve, and
the x axis of coverage-based curves. It is closed form, so it costs no more
than the expected richness.


## Names

- `rarefaction_curve()`: the field's own term, spelled out (rule 4), in
  snake_case (rule 5); sorts beside `rarefy()`. Rejected: `rarecurve`
  (vegan's name; its result is a plot and a list of vectors, and it would
  mask vegan's), `alpha_rarefaction` (QIIME 2's name; too close to
  `alpha_div_rarefied()`, which averages at one depth), `rarefaction_curves`
  (function names are singular, as `ternary_plot()`).
- Arguments reuse `alpha_div_rarefied()`'s: `metric`, `times`, `average`,
  `cutoff`, `digits`, `tree`, `q` ([alpha_metrics.md](alpha_metrics.md)),
  `margin`, `seed`, `cpus`; and `depth` → `depths` (several). New:
  `n_depths` (an exact count of grid points, `n_<unit>` rule; q2's
  `steps`), `max_depth` (q2's name; the largest grid depth), `analytic`.
- Same names in other packages: none known. Same purpose:
  `vegan::rarecurve()`, `iNEXT::iNEXT()` (see
  [hill_extrapolation.md](hill_extrapolation.md)), q2's
  `alpha-rarefaction`.
- No registry row: it repeats metrics, as `alpha_div_rarefied()` does.


## Design

### Signature

```r
rarefaction_curve(
  counts, metric = 'observed', depths = NULL, n_depths = 20L,
  max_depth = NULL, times = 10L, analytic = NULL, average = 'mean',
  cutoff = 10L, digits = 3L, tree = NULL, q = 2, margin = 1L,
  seed = 0L, cpus = n_cpus() )
```

- `metric`: one or more alpha metric names (`match_metric(div = 'alpha')`
  each, so aliases and prefixes work). Every iteration's rarefied table is
  used for all of them, as q2 does.
- `depths`: whole numbers ≥ 1 (sorted, duplicates dropped). With `NULL`,
  `n_depths` points from 1 to `max_depth`, as q2's `np.linspace(min_depth,
  max_depth, num = steps, dtype = int)`: compute \eqn{1 + i (M - 1)/(n -
  1)}, \eqn{i = 0, ..., n - 1}, with the last point exactly \eqn{M}, and
  truncate toward zero; then drop duplicates. The fixture checks this
  against NumPy. Default 20 points (q2: 10, vegan: every integer; say how
  to get each).
- `max_depth`: the largest library size by default (vegan's and q2's
  limit). Above `max_depth` nothing is computed.
- **Each sample's own library size is added to its grid** (vegan: "always
  including 1 and the total sample size"), so each curve ends at the
  sample's observed value. Depths above a sample's total have no row (q2
  gives `NaN` there).
- `times`: iterations per depth for resampled metrics (q2's `iterations`,
  default 10). `average`: `'mean'` (default) or `'median'` of the
  iterations, as `alpha_div_rarefied()` (q2's visualisation shows the
  median).
- `analytic`: `NULL` (default) means the exact expectation for
  `'observed'` and resampling for every other metric; `FALSE` resamples
  `'observed'` too (to compare with q2); `TRUE` for a metric without a
  closed form is an error naming the ones that have one.
- `seed`: iteration \eqn{t} at depth \eqn{d} is exactly `rarefy(counts, d,
  seed = s_t)`, with \eqn{s_t} the sequence of `rarefy(times =)` (`seed`,
  `seed + 1`, ...; `rarefied_setup()`'s `seeds`), **the same at every
  depth**.

### Why the same seed at every depth (and a test of it)

`rarefy()` keeps reads by Knuth's selection sampling: for each read in
turn, one PCG draw, kept when `draw % not_tried < still_need`. With the
same stream, a smaller target never keeps a read the larger one drops
(its `still_need` is never larger, and the draw and `not_tried` are the
same until it stops). So iteration \eqn{t}'s subsamples are **nested**
across depths: each iteration is one accumulation of reads, its observed
richness never decreases with depth, and curves are smooth rather than
jagged. q2 draws each depth independently; the means agree, the
per-iteration curves differ (say so in the docs). Check the dense, `Matrix`
and `slam` paths of `src/rarefy.c` keep this property, and test it:
`all(rarefy(x, d1, s) <= rarefy(x, d2, s))` for \eqn{d_1 < d_2}.

### Resampling

1. `rarefied_setup()`-style checks (move what both need into a shared
   helper): integer counts, metric(s), `times`, `average`, `seed`, `cpus`,
   `tree` from phyloseq-like objects.
2. For \eqn{t = 1..}`times`, for each depth \eqn{d}: `C_rarefy(counts, d,
   s_t, margin, cpus)`, samples whose total is below \eqn{d} set aside (they
   get no row at \eqn{d}), then each metric's function (as
   `alpha_div_rarefied()` calls it, warnings collected once with
   `rarefied_quietly()` / `rarefied_rewarn()`). Accumulate per (sample,
   metric, depth): the running sum and Welford's M2 for `'mean'` and `sd`,
   or the `n × times` block for `'median'` (`C_row_medians`).
3. Memory: one rarefied table at a time, plus the summaries (or the
   median's blocks).
4. Optional, only if profiling shows the per-depth `C_rarefy` calls
   dominate: a `C_rarefy_depths` that rarefies a sample to every depth in
   one pass over its reads (one draw per read, one `still_need` per depth),
   giving tables identical to separate calls.

Results are identical for any `cpus`: rarefaction seeds by sample
position, metrics are per-sample, and the summaries run in iteration order.

### Closed forms (`src/rarefaction_curve.c`)

Per sample, threads over samples, from the sample's frequency counts
\eqn{f_k} (the number of features with exactly \eqn{k} reads; \eqn{K}
distinct \eqn{k}, at most \eqn{\sqrt{2N}}):

- Absence probabilities \eqn{a_k(d) = \binom{N-k}{d}/\binom{N}{d} =
  \prod_{j=0}^{k-1} \frac{N - d - j}{N - j}} (0 when \eqn{N - k < d}),
  built by the running product over \eqn{k = 1, 2, ...} up to twice the
  largest count (the variance needs \eqn{a_{k+l}}), stopping once it
  underflows to 0. This avoids `lchoose()`
  differences of large numbers.
- \eqn{E[S_d] = S - \sum_k f_k a_k(d)}.
- Heck et al. (1975)'s variance, grouped by abundance class:
  \eqn{\mathrm{Var}(S_d) = \sum_k f_k a_k(1 - a_k) + 2\sum_{k<l} f_k f_l
  (a_{k+l} - a_k a_l) + 2\sum_k \binom{f_k}{2}(a_{2k} - a_k^2)}, where
  \eqn{a_{k+l}} is the probability that two features with \eqn{k} and
  \eqn{l} reads are both absent. O(\eqn{K^2}) per depth. The `sd` column is
  its square root (vegan's `se`).
- Coverage: for \eqn{d < N}, \eqn{C_d = 1 - \sum_k f_k \frac{k}{N}
  \binom{N-k}{d}/\binom{N-1}{d}} (Chao & Jost 2012; Chao et al. 2014, Table
  1), with the same running products. At \eqn{d = N}, the reference
  sample's estimate \eqn{\hat C = 1 - \frac{f_1}{N}\frac{(N-1)f_1}{(N-1)f_1
  + 2f_2}} (Chao et al. 2014, Eq. 12; with \eqn{f_2 = 0}, \eqn{(N-1)(f_1 -
  1)/((N-1)(f_1 - 1) + 2)}, their Appendix G), so the curve joins
  [hill_extrapolation.md](hill_extrapolation.md)'s. Share these helpers
  with that plan's C.
- Sums in a fixed order per sample; results identical for any `cpus`.

### Result

A data frame of class `c('ecodive_rarefaction_curve', 'data.frame')`, one
row per sample × metric × depth (depths up to the sample's total), ordered
by sample, metric, depth. Columns:

- `sample` (name, or position as text), `metric` (the id), `depth`;
- `value`: the expected value (closed form) or the iterations' mean or
  median;
- `sd`: Heck's standard deviation (closed form; 0 at the sample's total),
  or the iterations' SD (n − 1 denominator; `NA` with `times = 1`);
- `coverage`: the expected sample coverage at that depth (always closed
  form; it doesn't depend on the metric).

Attributes: `times` and `average` (they qualify `value` and `sd` for the
resampled metrics), `analytic` (the metric ids computed in closed form),
`totals` (each sample's library size, named). Printed by
`print_table_result` with the header "Rarefaction curves".

### Later step (ask the user): coverage-based rarefaction of a table

Chao & Jost (2012) rarefy samples to **equal coverage** rather than equal
depth. With \eqn{C_d} in hand this is cheap: for each sample, the
smallest depth whose expected coverage reaches a target, then `rarefy()`
with one depth per sample (which `qmp_rarefy()` already uses). Proposed as
`rarefy(counts, coverage = 0.95)` (a new last argument; `depth` and
`coverage` exclusive; samples whose own coverage is below the target are
dropped or kept as `drop` says). Any metric can then be computed at equal
coverage, not just Hill numbers.


## Reference implementations

- **q2-diversity** (github.com/qiime2/q2-diversity, dev branch; BSD-3,
  "Copyright (c) 2016-2026, QIIME 2 development team", LICENSE read;
  `_alpha/_visualizer.py` read 2026-11-20): `alpha_rarefaction(table,
  max_depth, metrics = {observed_features, shannon[, faith_pd]}, min_depth
  = 1, steps = 10, iterations = 10)`; the grid is `np.linspace(min_depth,
  max_depth, num = steps, dtype = int)`; each (depth, iteration) runs
  `feature_table.rarefy` and `diversity.alpha`; samples below a depth are
  `NaN`; summaries are quantiles (2, 9, 25, 50, 75, 91, 98%) over
  iterations, and group medians by metadata column. It refuses `osd`,
  `lladser_ci`, `strong`, `esty_ci`, `kempton_taylor_q`, `chao1_ci`. Only
  the grid rule and the design are taken; no code ported, so no notice.
- **vegan** 2.7-8 (GPL-2, `DESCRIPTION` read 2026-11-20), output-only:
  `rarefy(x, sample, se = TRUE)` (Hurlbert's expectation, Heck's SE; its
  docs) and `rarecurve(x, step, tidy = TRUE)`.
- **iNEXT** 3.0.2 (GPL ≥ 3), output-only: the `SC` column of
  `iNEXT(x)$iNextEst$size_based` for \eqn{m \le n} is \eqn{C_d}; one shared
  fixture with [hill_extrapolation.md](hill_extrapolation.md).
- **Fixtures:** `data-raw/expected/rarefaction_curve.r` (vegan, iNEXT) and
  `rarefaction_curve.py` (NumPy's grid for a few `(min, max, steps)`),
  text under `data-raw/expected/rarefaction_curve/`, into
  `tests/tinytests/expected/rarefaction_curve.rds`. Data: `ex_counts`, and
  one deep sample (\eqn{N \approx 10^5}) to test the products' range.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots". The input
is a count table with no title, so the result has no `display`.

```r
plot(x, metric = NULL, groups = NULL, data = NULL,
     type = 'depth', depth = NULL, ...)
```

- `type = 'depth'` (default): `value` against `depth`, one line per
  sample; `'coverage'`: `value` against `coverage` (coverage-based
  curves); `'completeness'`: `coverage` against `depth` (the sample
  completeness curve, Chao et al. 2014, Fig. 4).
- `metric`: which metric (default the first in the result).
- `groups` (+ `data`): colour samples by group, with a legend, and draw
  each group's median curve thicker (q2's view). Matched to the result's
  `sample` labels as elsewhere (`per_sample()`).
- `depth`: a vertical line at a planned rarefaction depth (vegan's
  `sample`), labelled with how many samples it would keep (from `totals`).
- Title: "Rarefaction Curves of " and the metric's `METRICS` display name,
  e.g. "Rarefaction Curves of Observed Features"; for `'completeness'`,
  "Sample Completeness Curves". Y label: the metric's display name. A
  subtitle names the summary: "Expected values" (closed form) or "Mean of
  10 rarefactions" (from `times`, `average`, `analytic`).
- Reads only columns and attributes; nothing is recomputed. Base graphics
  through `plot_do()` (R/plot.r), so `main`, `col`, `lty` ... override.


## Tests (`tests/tinytests/test-rarefaction_curve.r`)

- **Closed forms:** `value` and `sd` for `'observed'` equal vegan's
  `rarefy(se = TRUE)` fixture at 1e-10 (✱ Heck's pair terms cancel; relax
  to 1e-8 only if a deep sample needs it, and say why); equal a direct R
  computation with `lchoose()` on small samples at 1e-12; `coverage`
  equals iNEXT's `SC` at 1e-10 and the direct formula; at the sample's
  total, `value` = observed, `sd` = 0, `coverage` = Eq. 12.
- **Resampling:** with `times = 1`, each value is `identical()` to
  `alpha_div(rarefy(counts, d, seed), metric)`; with `times = k`, the mean
  is `identical()` to the sequential sum over `rarefy(..., seed = s_t)`
  divided by k (as rarefied_diversity's tests); the median equals
  `apply(, 1, median)`.
- **Nesting:** `rarefy(x, d1, s) <= rarefy(x, d2, s)` elementwise for
  dense, `Matrix` and `slam`; resampled observed richness is non-decreasing
  in depth within every iteration.
- **Statistics:** resampled `'observed'` (`analytic = FALSE`, `times =
  2000`) agrees with the closed form within 4 Monte Carlo SEs.
- **Grid:** `depths = NULL` equals NumPy's grid (fixture); duplicates
  dropped; each sample's total added; no rows above a sample's total.
- Several metrics, phylogenetic (`faith` with `tree`), `q` passed to
  `hill`; `margin = 2`; `identical()` for `cpus` 1, 2 and 4.
- **Errors:** non-integer counts; `depths` with 0, NA or 1.5; `n_depths`
  < 2; `max_depth` < 1; `analytic = TRUE` with `'shannon'`; unknown metric
  or `average`.
- **Plot:** runs for every `type`, with and without `groups`; uses the
  title above; `main =` overrides.
- **Smoke:** one line, small `times`.


## Docs

- Rd page `rarefaction_curve` (+ `plot.ecodive_rarefaction_curve`):
  reading a curve (plateau = enough depth), choosing a rarefaction depth
  (the vertical line, `totals`), closed form versus resampling, the nested
  draws versus q2's independent ones, the coverage curve and how it links
  to `hill_extrapolation()`; how to reproduce q2 (`n_depths = 10`,
  `analytic = FALSE`, `average = 'median'`) and vegan (`depths =
  seq(1, max(totals))`).
- "Same names in other packages": none; `vegan::rarecurve()` does the same
  for richness.
- pkgdown: "Utilities & Helpers" beside `rarefy` and `alpha_div_rarefied`;
  `vignettes/adiv.Rmd`: a rarefaction-curve example. NEWS; NAMESPACE
  (`export(rarefaction_curve)`, `S3method(plot, ecodive_rarefaction_curve)`,
  `S3method(print, ...)`).
- parameter_names.md §1: `depths`, `n_depths`, `max_depth`, `analytic`.


## Dependencies

- `rarefy()`'s C (`C_rarefy`) and `rarefied.r`'s `rarefied_setup()`,
  `rarefied_draw()`, `rarefied_quietly()`, `C_row_medians` (item 30).
- `q` exists only after [alpha_metrics.md](alpha_metrics.md); drop it if
  this lands first.
- The coverage helpers are shared with
  [hill_extrapolation.md](hill_extrapolation.md); whichever lands first
  writes them.


## References

- Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and alternative parameters. *Ecology*, 52, 577–586. doi:10.2307/1934145
- Heck, K. L., van Belle, G., & Simberloff, D. (1975). Explicit calculation of the rarefaction diversity measurement and the determination of sufficient sample size. *Ecology*, 56, 1459–1461. doi:10.2307/1934716
- Chao, A., & Jost, L. (2012). Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. *Ecology*, 93, 2533–2547. doi:10.1890/11-1952.1
- Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. *Ecological Monographs*, 84, 45–67. doi:10.1890/13-0133.1
- Bolyen, E., et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*, 37, 852–857. doi:10.1038/s41587-019-0209-9
- Schloss, P. D. (2024). Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere*, 9, e00354-23. doi:10.1128/msphere.00354-23


## Decisions made while implementing

Recorded 2026-11-24 (item 38).

**Asked the user, at the start:**

1. **Several metrics per call: yes.** `metric` takes one or more alpha
   metrics. The result is always the same long data frame (a `metric`
   column), whatever the number of metrics, so the return type never
   depends on the input's length (the user's rule in
   [alpha_metrics.md](alpha_metrics.md), decision 4). Each rarefied table
   is used for every resampled metric. Duplicates (by ID, so `'obs'` and
   `'observed'`) are dropped; the order given is kept.
2. **`times = 10` by default**, as q2-diversity's `iterations` (not
   `alpha_div_rarefied()`'s 100).
3. **Coverage-based rarefaction done now**, as `rarefy(coverage =)`, a new
   last argument (after `warn`). Each sample's depth is the smallest
   \eqn{d \in [1, N]} whose expected coverage reaches the target (by
   bisection: \eqn{C_d} increases with \eqn{d}, and at \eqn{d = N} the
   Eq. 12 estimate is at least \eqn{C_{N-1} = 1 - f_1/N}). A sample whose
   whole-sample estimate falls short gets depth \eqn{N + 1}, so `drop`
   drops it (or `drop = FALSE` returns it unrarefied) exactly as for an
   unreachable `depth`; `warn` says "k samples have an estimated coverage
   below c and will be dropped / returned unrarefied". Then it rarefies with
   one depth per sample, so the draws are those of `rarefy(depth = <those
   depths>)` with the same seed. Giving both `depth` and `coverage` is an
   error ("Give `depth` or `coverage`, not both."); `validate_coverage()`
   takes a single number in (0, 1]. No attribute is added (count tables
   stay plain); `rowSums()` gives the depths. Empty samples are dropped.

**Decided while implementing:**

4. **`bias_corrected`** is a formal too (after `q`), passed to `chao1`, as
   in `alpha_div_rarefied()` (the plan predates it).
5. **Shared helpers in `R/rarefied.r`:** `rarefied_setup()` is split into
   `rarefied_counts()` (counts, margin, integer, seed and cpus checks,
   negatives, doubles), `rarefied_times()`, `rarefied_average()` and
   `rarefied_seeds()`, which `rarefaction_curve()` reuses; the old
   function's checks, order and messages are unchanged.
6. **Validation order:** metrics (a misspelt metric is the error
   reported), counts/margin/integers/seed/cpus/negatives, `times`,
   `depths`, `n_depths`, `max_depth`, `analytic`, `average`, `cutoff`,
   `digits`, `q`, `bias_corrected`. Messages: "`metric` must be one or
   more alpha diversity metric names.", "`depths` must be NULL or whole
   numbers of at least 1.", "`n_depths` must be a whole number of at
   least 2.", "`max_depth` must be NULL or a whole number of at least 1.",
   "`analytic` must be NULL, TRUE or FALSE.", "`analytic = TRUE` needs a
   closed form, which only "observed" has; not "shannon". Use `analytic =
   NULL` to resample the other metrics."
7. **`max_depth` with `depths`:** defaults to `max(depths)`; depths above
   `max_depth` are dropped. A sample's total is added to its depths only
   when it is at most `max_depth` (with `depths = c(10, 400, 700)`, a
   sample of 886 gets rows 10, 400, 700; one of 615, rows 10, 400, 615).
8. **The grid** is computed as NumPy does it: `(0:(n-1)) * ((M - 1) / (n -
   1)) + 1`, the last set to \eqn{M}, `floor()`, `unique()`. Identical to
   `np.linspace(1, M, num = n, dtype = int)` for 12 cases (fixture). No
   error when `max_depth` exceeds every total (q2 refuses): those depths
   just have no rows.
9. **Empty samples** have no rows and keep their 0 in `totals`.
10. **Unnamed samples** are labelled "1", "2", ... in `sample` and
    `totals`; `plot()` treats labels equal to the positions as unnamed, so
    `groups` is matched in order.
11. **Resampled rows at a sample's own total:** `rarefy(counts, N_s)`
    leaves sample \eqn{s} as it is, so these values come from the metric
    on the unrarefied samples, computed once and added in every
    iteration (so the mean is still the sequential sum over iterations
    divided by `times`, and `identical()` to the rebuild). Only the grid
    depths are rarefied: `n_depths` rarefactions per iteration, not
    `n_depths` + the number of samples.
12. **Each rarefaction is `C_rarefy` of the whole table** then the samples
    below the depth dropped (as `rarefied_draw()`), never a subset first,
    since streams are seeded by sample position.
13. **Summaries:** mean = sequential sum / `times`, sd by Welford (as the
    beta mean); median by `C_row_medians` over the `R × times` matrix, sd
    two-pass from it. sd is `NA` with `times = 1`. The median equals
    pandas' `describe()` 50% (q2's summary) exactly (fixture).
14. **Closed forms in C, no `lchoose()`:**
    - absence probabilities by one running product over the sorted
      distinct abundances (`rc_absence()`), stopping at 0;
    - **Heck's variance without cancellation.** As a difference,
      \eqn{a_{k+l} - a_k a_l} lost up to 10 digits for deep samples at
      small depths (and vegan's lchoose() form loses as many). With
      \eqn{g_m = a_{k+m}/a_k}, \eqn{\delta_m = g_m - a_m} follows
      \eqn{\delta_{m+1} = \delta_m s_m + g_m (r_m - s_m)}, \eqn{r_m - s_m =
      -dk/((N-k-m)(N-m))} computed exactly; all terms have one sign. Cost
      per depth: at most the sum of the distinct abundances (≤ N) steps,
      O(K) memory (no K² pair table).
    - **Exact values where the formulas are exact:** `value` = 1 and `sd`
      = 0 at depth 1, `sd` = 0 at the total (where Eq. 12 gives the
      coverage).
15. **Threads over rows**, not samples (`row r` to thread `r % n`), so a
    few samples with a dense grid still use every thread (20 samples × 3,000
    depths: 21 s → 5.3 s on 4 cores). Identical for any `cpus`.
16. **The coverage helpers** (`rc_freq()`, `rc_absence()`,
    `rc_coverage()`, `rc_coverage_full()`) are declared in `ecodive.h` for
    `hill_extrapolation()`; [hill_extrapolation.md](hill_extrapolation.md)
    says so.
17. **Plot:** `legend = 'bottomright'` and `col` are formals (as
    `plot.ecodive_ord()`); without groups every line is grey
    (`PLOT_MUTED`), with groups each group's colour and its median curve 3×
    as thick: at each depth shared by at least two of its samples (all its
    depths for a group of one), the median of x and of y. Numeric `groups`
    become a factor (one colour per value). `depth` with `type =
    "coverage"` is an error (the x axis is coverage). Subtitles: "Expected
    values", "Mean of 10 rarefactions", "Median of 1 rarefaction";
    `type = "completeness"` always "Expected values". The labels come from
    an internal `rc_plot_labels()`, which the tests check.
18. **No `display`:** the input is a count table with no title; the
    print header is "Rarefaction curves" (`print_df_label()`).
19. **References for the closed forms:** besides vegan 2.7.8 and iNEXT
    3.0.2 (output only), the fixture has a 50-digit mpmath computation of
    the published formulas (`data-raw/expected/rarefaction_curve.py`),
    which shows ecodive within 5.4e-13 of the exact values everywhere,
    vegan within 1.4e-12 (richness) and 7e-6 (sd, absolute; 1e-8 relative
    where sd > 0.001), iNEXT within 1.5e-10 (coverage). The tests compare
    with vegan at 1e-12 (richness) and 1e-10 (sd > 0.001), relaxed to
    1e-11 and 1e-7 for the deep samples because of vegan's own error.
20. **q2-diversity:** nothing ported (only the grid rule and design), so
    no notice in `inst/COPYRIGHTS` and no `Authors@R` entry; the R file's
    header names it.
21. **pkgdown:** `rarefaction_curve` under "Utilities & Helpers" beside
    `rarefy` and `alpha_div_rarefied`; its `plot()` method in "Plotting and
    Printing" with the others.
