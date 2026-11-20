# `sourcetracker()` and `feast()`: microbial source tracking

**Priority:** 47 (`sourcetracker()`) and 58 (`feast()`). **Tier:**
- `sourcetracker()` is a port of SourceTracker2 (Python, **BSD-3**,
  github.com/biota/sourcetracker2, master 2.0.1-dev).
- `feast()` is **clean room**: FEAST's code is CC BY-NC-SA 4.0, which is
  non-commercial and can't be ported.

Read [README.md](README.md) first.


## What it is

Source tracking estimates what share of each **sink** sample's community
came from each of several **source** environments, plus an **unknown**
source. Examples are how much of an infant's gut came from the mother's
skin, vagina or stool, or how much of a hospital surface came from
patients' skin or the water.

| Method | Model | Citations (OpenAlex) | Notes |
|:-------|:------|:---------------------|:------|
| SourceTracker (Knights et al. 2011) | Bayesian; a Gibbs sampler reassigns each sink read to a source | about 2,000 in total; 190–260 a year, steady | The standard. v1 is an R script (GPL); v2 is Python (BSD-3) with a QIIME 2 plugin. |
| FEAST (Shenhav et al. 2019) | maximum likelihood by EM, re-estimating the source distributions too | about 720 in total; 140 a year and rising | Now the usual benchmark; papers often run both. 30–300× faster than SourceTracker. |

Neither has replaced the other, so ecodive plans both. SourceTracker comes first
because it can be ported. FEAST is clean room and comes later. STENSL (a sparse
FEAST, An et al. 2022) is cited about 5 times a year and isn't planned.

**Why ecodive:** SourceTracker2 loops over every read in Python, for every
pass, restart and sink. That is slow, and SourceTracker2's parallelism
(`multiprocessing`) isn't even reproducible (see below). In C the same
loop is a tight integer kernel, and (sink, restart) pairs are independent
work units: a natural fit for `run_parallel()`.


## `sourcetracker()`

### Signature

```r
sourcetracker(counts, sources, data = NULL, margin = 1L,
              depth = 1000L, sink_depth = depth,
              alpha1 = 1e-3, alpha2 = 0.1, beta = 10,
              n_restarts = 10L, n_draws = 1L, burnin = 100L, delay = 1L,
              loo = FALSE, features = FALSE, seed = 0L, cpus = n_cpus())
```

- **`sources`:** the samples' environment labels, one per sample, or a
  column name of `data` (resolved with `match_rows()` / `per_sample()`).
  **`NA` marks a sink.** SourceTracker2 needs two metadata columns
  (`SourceSink`, `Env`); one vector with `NA` for sinks says the same
  thing more simply. To leave a sample out, subset `counts`. **Ask the
  user** whether a separate `sinks =` selector is wanted.
- **`depth`, `sink_depth`:** rarefaction depths for the collapsed sources and
  for the sinks. They reuse `rarefy(depth)`'s name; SourceTracker2's flags are
  `--source_rarefaction_depth` and `--sink_rarefaction_depth`, which default to
  1000, with 0 meaning none. In ecodive, `NULL` means no rarefaction. A source or sink
  shallower than its depth is an **error** naming the samples, as in
  SourceTracker2. Rarefy with ecodive's `rarefy()` C code and `seed`.
- **`alpha1`, `alpha2`, `beta`:** the priors, under SourceTracker's own
  names, which users know from its documentation and tutorials:
  - `alpha1` is the pseudocount on known sources;
  - `alpha2` is the unknown-source prior, per read of the sink;
  - `beta` is the pseudocount on environment counts.

  `alpha1` and `alpha2` don't clash with the glossary's `alpha` (none of
  them is a significance level), but **add them to parameter_names.md
  §4** as kept on purpose.
- **`n_restarts`, `n_draws`, `burnin`, `delay`:** the sampler schedule.
  They follow SourceTracker2's defaults (10, 1, 100, 1), and the
  `n_<unit>` rule for the two exact counts. `burnin` and `delay` keep
  SourceTracker's names; **ask the user** whether they should be
  `n_burnin` and `thin`. v1's `delay` default is 10, and the docs say
  so.
- **`loo`:** leave-one-out over **environments**, as SourceTracker2 does.
  Each collapsed, rarefied environment in turn is treated as a sink
  against the others, and the result is an environment × environment
  matrix with 0 on the diagonal. This differs from v1, which leaves out
  individual source samples; the docs say so.
- **`features`:** also return the per-feature contributions. SourceTracker2's
  CLI has these off and its QIIME 2 plugin has them on. Default `FALSE`, since
  they are a sinks × environments × features block.

### Algorithm (SourceTracker2, read from source)

The code is in `sourcetracker/_sourcetracker.py`: `gibbs`, `gibbs_sampler`,
`ConditionalProbability`, `collapse_source_data` and `subsample_dataframe`,
with the pipeline in `_gibbs.py` and defaults in `_gibbs_defaults.py`.

**Preprocessing:**

1. Collapse the source samples by environment using the **mean**, then
   **truncate** to integer (`astype(int32)`, not rounding, although the
   README says "rounded"). Sort the environments by name (C locale).
2. Rarefy each collapsed source to `depth`, without replacement.
3. Rarefy each sink to `sink_depth`. Under `loo`, nothing is rarefied
   again: the sinks are the already rarefied environments.

**Per sink.** Let V be the number of environments plus 1 (the unknown), n
the sink's depth, \eqn{\tau} the number of features, and \eqn{m_{vt}} the
counts of source v.

- \eqn{\tau} counts **every feature column in the table**, including
  features that are zero everywhere, as SourceTracker2 does. Filtering
  empty features therefore changes the result. Match SourceTracker2, and
  document it.
- **Fixed part:** \eqn{q_{vt} = (m_{vt} + \alpha_1)/(m_{v\cdot} +
  \tau\alpha_1)}.
- **State:** each sink read s has an environment \eqn{e_s}. Keep:
  - environment counts \eqn{c_v};
  - unknown-source feature counts \eqn{u_t};
  - \eqn{U = \sum_t u_t}.
- **Initialise each restart:** give every read a uniform environment in 0..V−1.
  Build \eqn{u}, \eqn{U} and \eqn{c} from that assignment. Nothing carries over between restarts.
- **Each pass:**
  1. Shuffle the read order.
  2. Take each read (feature t) out of its environment.
  3. Compute the weights:
     - known v: \eqn{p_v = q_{vt}(c_v + \beta)};
     - unknown: \eqn{p_U = (u_t + \alpha_2 n)/(U + \alpha_2 n \tau)\,(c_U
       + \beta)}.

     The common denominator \eqn{n - 1 + \beta V} is dropped.
  4. Draw \eqn{w \sim U(0, \sum p)} and pick the first v whose cumulative
     sum is ≥ w (NumPy `searchsorted`, side left).
  5. Add the read back to v.
- **Schedule:** `burnin + (n_draws − 1) * delay + 1` passes per restart.
  Record a draw (copy \eqn{c}, and \eqn{e} if `features`) after pass
  `burnin + 1`, then every `delay` passes.

**Outputs:**
- **Proportions:** sum \eqn{c} over all draws of all restarts, then
  divide by the total.
- **`sd`:** the standard deviation (**ddof = 0**, as `np.std`) of the
  per-draw proportions. With the defaults it rests on only 10 draws; the docs say so.
- **Features:** read counts by environment × feature, **summed** over all
  draws. SourceTracker2's docstring says they're averaged, but the code sums them. Match the code,
  and say so in the docs. **Ask the user** whether ecodive should
  return means instead and note the difference.

### Random numbers and reproducibility

- **Exact parity with SourceTracker2 isn't feasible:**
  - its rarefaction goes through scikit-bio and biom's PCG64 subsampler;
  - the sampler draws from NumPy's legacy MT19937 (`choice`, `shuffle`,
    `uniform`);
  - its `multiprocessing` workers either inherit identical RNG states
    (fork) or start unseeded (spawn), so its output depends on `jobs`
    and the OS.

  SourceTracker2's own tests compare with v1 only to `atol = 0.01`.
  Match the **deterministic parts exactly** and the rest **by
  distribution** (see Tests). This is the "seed-like difference" case in
  README.
- **ecodive:** one PCG stream per (sink, restart),
  `pcg32_srandom_r(&rng, seed, sink * n_restarts + restart)`, consumed in a fixed
  order: the initial assignment, then per pass the Fisher–Yates shuffle,
  then one uniform per read. Rarefaction uses `rarefy()`'s own seeding.
  The result is identical for any `cpus`, and independent of which other sinks are in the
  table.
- Use `pcg32_boundedrand_r` for the shuffle and initial assignment, and
  a 53-bit uniform double. If the planned `src/rng_dist.c` (see
  [dirmult_ttest.md](dirmult_ttest.md)) exists by then, use it; otherwise
  create it here.

### Performance

- **Cost per sink:** restarts × passes × n × V, which with the defaults is
  10 × 101 × 1000 × V operations. The work unit is (sink, restart): with
  20 sinks that is 200 units, spread evenly over threads. Reduce draws
  per sink in restart order.
- **Buffers per thread:**
  - the read environments \eqn{e} (n ints) and the read order (n ints);
  - \eqn{c} (V);
  - \eqn{u}, which is sparse: the reads touch at most n features, so
    reset only the touched entries, never all \eqn{\tau}.
- \eqn{q} is environments × features, computed once and shared read-only.
  For a sink, only the columns of features it contains are ever read.
- With `features`, each work unit adds V × (the sink's features) counts.
  Return them sparse: a list with one `dgCMatrix`-free triplet per sink,
  or a dense matrix per sink if that's small. Decide when implementing.

### Result

- A numeric **matrix** of sinks × (environments..., `"Unknown"`). Each row
  holds the mean proportions and sums to 1.
- Attributes:
  - `sd` (same shape);
  - `n_draws` (draws in total);
  - `depth`, `sink_depth`;
  - `alpha1`, `alpha2`, `beta`;
  - `seed`;
  - `features` (only when requested): a list with one environment ×
    feature count matrix per sink.
- Under `loo`, it is an environments × environments matrix, without
  `"Unknown"` in the rows.
- Plain matrix, no new class. That's consistent with
  `replace_zeros_*()`, and it plots directly as a stacked bar chart; the
  docs show how.

### Not ported

- v1's `tune.st()`, which tunes `alpha1` and `alpha2` by leave-one-out. It is GPL and
  not in SourceTracker2. It could be written later as a loop over
  `loo = TRUE` fits.
- SourceTracker2's diagnostics and plots.


## `feast()` (clean room, later)

### Licence

FEAST (github.com/cozygene/FEAST, version 0.1.0, not on CRAN) and its STENSL
branch are **CC BY-NC-SA 4.0** on every branch. That is non-commercial,
so FEAST is **output-only, and even that needs the user's OK**: storing
its outputs as fixtures raises the same question as PHATE's (see
[phate.md](phate.md)). Write `feast()` from the paper (Shenhav et al.
2019, open access via PMC/eScholarship) and FEAST's README and Rd files
only. **Don't read its R or C++ source.** FEAST2 (a fork) has the same licence.

### Model (from the paper)

The sink \eqn{x \sim \mathrm{Mult}(C, \sum_i \alpha_i\gamma_i)} over K
known sources plus one unknown. Each known source sample \eqn{y_i \sim
\mathrm{Mult}(C_i, \gamma_i)}, so the source distributions are
re-estimated together with the mixing proportions.

- **E-step:** \eqn{p(i \mid j) = \alpha_i\gamma_{ij} / \sum_k
  \alpha_k\gamma_{kj}}.
- **M-step:**
  - \eqn{\gamma_{ij} \leftarrow (x_j p(i \mid j) + y_{ij}) / \sum_j(x_j
    p(i \mid j) + y_{ij})}; the unknown source has no \eqn{y} term;
  - \eqn{\alpha_i \leftarrow \sum_j x_j p(i \mid j) / C}.
- **Stopping:** at most 1000 iterations (`EM_iterations`), and a 1e-6
  threshold. The paper doesn't say which quantity the threshold is
  tested on. Use the change in log-likelihood, and document that.
- **Sources are individual samples,** not collapsed environments (FEAST's
  default). An `environments` option could pool them by label.

### Decisions the paper leaves open (ask the user)

1. **Initialisation** of \eqn{\alpha} and of the unknown source isn't in
   the main paper (it may be in the supplement, which the research agent
   didn't read). Read the supplement first. If it's absent, use uniform
   \eqn{\alpha}, and initialise the unknown source to the normalised
   positive residual \eqn{\max(0, x - \sum_i \bar\alpha_i y_i/C_i)}.
   STENSL's paper initialises it in the same spirit.
2. **Rarefaction.**
   - FEAST rarefies each sink and its sources to the minimum depth in
     the group (`COVERAGE`), and doesn't document whether it samples
     with or without replacement.
   - ecodive's selling point is determinism, and FEAST's EM is
     deterministic apart from this step.
   - So the proposed default is `depth = NULL` (no rarefaction), with
     `depth` available to rarefy with `rarefy()`.
   - Say in the docs that FEAST's defaults rarefy.
3. **Per-sink sources.** FEAST can pair each sink with its own sources
   (`different_sources_flag`, which works through an `id` column). Supporting it
   needs a `pairs`-style argument; leave it out at first.

### Implementation notes

- Each sink is an independent work unit. The cost is O(iterations × (K+1) ×
  features present) per sink. Restrict every update to the sink's and
  the sources' non-zero features (the union), since \eqn{\gamma_{ij}} for a feature
  absent from both stays 0.
- The function is deterministic, with no `seed` unless `depth` is set.
- Return the same structure as `sourcetracker()` (a sinks × (sources...,
  `"Unknown"`) matrix of proportions), with attributes `n_iters`
  (per sink) and `log_lik`. FEAST gives no standard deviation.

### Tests

- **Algebra:**
  - EM on a synthetic sink that is an exact mixture of
    deep sources recovers the mixing proportions to about 1e-6;
  - the log-likelihood never decreases.
- **Against FEAST (only if the user approves its outputs as
  fixtures):** proportions within about 0.02 on simulated mixtures, given
  the same rarefied input tables (pass FEAST a pre-rarefied table with
  `COVERAGE` set to its depth, so rarefaction drops out).


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Attribute `method` (`'sourcetracker'` or `'feast'`). The inputs are
  count tables, which have no title, so the plot title is built from the
  registry's display name for `method`: "Source contributions
  (SourceTracker)", "Source contributions (FEAST)"; no `display`
  attribute.
- Store what plots need: per-sink source proportions (and their SDs, for
  SourceTracker's draws) in a fixed source order with an `unknown` source
  last, so a stacked-bar `plot()` of source contributions agrees across
  calls.

## Tests (`sourcetracker()`)

- **Deterministic kernels at 1e-10:** expose internal entry points and
  compare them with SourceTracker2's functions on fixed inputs:
  - `collapse_source_data` (the mean and truncation);
  - the fixed part \eqn{q};
  - the conditional probability vector for a fixed state
    (`ConditionalProbability`);
  - the draw schedule (which passes are recorded).
- **Whole function by distribution:**
  - simulate sinks as known mixtures of
    deep sources, plus an unknown component;
  - run SourceTracker2 with many seeds and store the means and between-seed
    SDs of its proportions (`data-raw/expected/sourcetracker.py`);
  - check that ecodive's proportions over several seeds fall within a few
    SDs, and within 0.02 of the true mixture where the sampler converges.
- **Output-only cross-check:** SourceTracker v1 (R, GPL) on the same tables,
  run in `data-raw/expected/sourcetracker.r`. Check that the
  proportions agree loosely.
- **Determinism:** identical results for `cpus` 1, 2 and 4 and for the same `seed`,
  and different results for another `seed`. A sink's result doesn't change when
  other sinks are added or removed.
- **Edge cases:**
  - one environment;
  - a sink with all its reads in features no source has (mostly unknown);
  - a source shallower than `depth` (error);
  - `depth = NULL`;
  - no sinks (error, unless `loo`);
  - `NA`-only `sources` (error);
  - non-integer counts (error, as `rarefy()`).
- One smoke-test line (small depth, 1 restart, burnin 5).


## Docs

- A new pkgdown group, **"Source Tracking"**, after "Compositional Data",
  with `sourcetracker` and later `feast`.
- **Explain:**
  - sinks, sources and the unknown source;
  - the assumption that a sink is a convex mixture of its sources (the
    FEAST paper notes where this fails);
  - why results depend on rarefaction depth;
  - that FEAST's benchmark found SourceTracker underestimates the
    unknown;
  - how to draw the usual stacked bar chart.
- **"Differences from SourceTracker2":**
  - per-sink seeding, which makes results reproducible across `cpus`;
  - `NA`-for-sink input;
  - v1 vs v2 differences in `delay` and in leave-one-out.
- **"Same names in other packages":** `sourcetracker` (v1's function, the
  same method) and `FEAST` (cozygene).
- NEWS.


## New argument names (add to [parameter_names.md](parameter_names.md))

| Name | Meaning |
|:-----|:--------|
| `sources` | per-sample source-environment labels; `NA` marks a sink |
| `sink_depth` | the rarefaction depth for sinks (`depth` is the sources') |
| `alpha1`, `alpha2`, `beta` | SourceTracker's priors (§4: kept as the reference names) |
| `n_restarts`, `burnin`, `delay` | the Gibbs schedule (`n_draws` already exists) |
| `loo` | leave-one-out over source environments |
| `features` | also return per-feature contributions |


## References

- Knights, D., Kuczynski, J., Charlson, E. S., et al. (2011). Bayesian
  community-wide culture-independent microbial source tracking. *Nature
  Methods*, 8, 761–763.
- Shenhav, L., Thompson, M., Joseph, T. A., et al. (2019). FEAST: fast
  expectation-maximization for microbial source tracking. *Nature
  Methods*, 16, 627–632.
- An, U., Shenhav, L., Olson, C. A., et al. (2022). STENSL: microbial
  source tracking with ENvironment SeLection. *mSystems*, 7, e00995-21.


## Decisions made while implementing (sourcetracker)

Item 47, 2026-11-20. Asked of the user at the start:

1. **No `sinks =` selector** (user's choice): `NA` in `sources` marks a
   sink; to leave a sample out, subset `counts`.
2. **`n_burnin` and `thin`** (user's choice), not SourceTracker's `burnin`
   and `delay`. The docs give the SourceTracker names.
3. **`features` returns means over the draws** (user's choice), so each
   sink's environment × feature matrix sums to its depth. The docs say
   SourceTracker2 sums, and that `* attr(x, "n_draws")` gives its tables.
4. **A class, `ecodive_source`** (user's choice), on the plain matrix:
   `c('ecodive_source', 'matrix', 'array')`, with `print()` and a
   stacked-bar `plot()` (unknown last, in grey; legend outside on the
   right; `horiz`, `legend`, `...` to `barplot()`). It is still a matrix.
   `feast()` should return the same class; both methods are shared.

Made while implementing:

5. **Random streams.** Rarefying sink row r of `counts` uses stream r, the
   stream `rarefy()` gives sample r, so a sink is rarefied exactly as
   `rarefy(counts, sink_depth, seed)` rarefies it (tested). A pooled
   environment e uses stream `1 << 56 | e`; the chain of sink row r and
   restart k `2 << 56 | r << 24 | k`; under `loo`, environment e's chain
   `3 << 56 | e << 24 | k`. The tags keep the streams apart (the plan's
   `sink * n_restarts + restart` would collide with the rarefaction
   streams and change every chain when `n_restarts` changes). `n_restarts`
   is capped at 2^24 - 1. A chain draws: each read's start
   (`pcg32_bounded(V)`), then per pass a Fisher–Yates shuffle
   (`pcg32_bounded`), then one `rng_unif()` per read. Because the stream is
   keyed by the sink's row in `counts`, a sink's result is unchanged when
   sinks *after* it are added or removed (tested), not when rows before it
   are; keying by sample name was considered and rejected (duplicates,
   unnamed tables, and `rarefy()`'s position-based streams).
6. **Rarefaction is done in `src/sourcetracker.c`** with the same Knuth
   selection sampling as `src/rarefy.c` (so identical draws), rather than
   by calling `rarefy()`, to choose the environments' streams.
7. **`sd` is the standard deviation (ddof 0) of the per-draw proportions,**
   as the plan says. SourceTracker2's `mps` divides each draw's counts by
   the total over *all* draws first, so its values are ours divided by the
   number of draws: a scaling bug. The docs say so, and the tests compare
   our `sd` with SourceTracker2's `mps × n_draws` (ratio within 0.75–1.33).
8. **Attributes:** `sd`, `n_draws`, `features` (when asked) and `method`
   only. The plan's `depth`, `sink_depth`, `alpha1`, `alpha2`, `beta` and
   `seed` attributes are left out under README's rule 1 (set 2026-11-20,
   after the plan): they are the caller's arguments, not results.
9. **Environment order:** a factor's levels (unused ones dropped), else
   sorted by name in the C locale (`sort(method = 'radix')`), as
   SourceTracker2's `groupby()`. `"Unknown"` as an environment name is an
   error.
10. **Argument ranges:** `alpha1 >= 0`; `alpha2 > 0` and `beta > 0`
    (SourceTracker2 allows 0, which gives 0/0 weights); `n_burnin >= 0`
    (SourceTracker2 wants > 0; 0 is harmless). A pooled environment with no
    reads is an error when `alpha1 = 0` (0/0) or under `loo` (an empty
    sink); with `alpha1 > 0` it is allowed (a uniform source), as in
    SourceTracker2. A sink with no reads (`sink_depth = NULL`) is an error
    naming it. Too-shallow sources and sinks are errors naming them.
11. **`loo`:** samples with `NA` are ignored; the result is environments ×
    (environments, `"Unknown"`) with 0 in each environment's own column
    (and its `sd`), and with `features`, a zero row for the held-out
    environment, as SourceTracker2's `collate_gibbs_results()`.
12. **Per-feature output** is a list of dense (environments + Unknown) ×
    features matrices, named by sink, with every feature column (so they
    line up with `counts`). Inside C, each (sink, restart) unit accumulates
    only the sink's own features.
13. **Threads:** units are (sink, restart) pairs, striped over the threads
    in blocks (so the user can interrupt between blocks); each unit writes
    only its own outputs, so the result is identical for any `cpus`
    (tested for 1, 2, 4, 7, also with `loo`). `threads_for()` keeps small
    jobs on one thread.
14. **The weights drop SourceTracker2's common denominator**
    `n - 1 + beta V`; `C_st_cond_prob()` (for the tests) divides by it to
    compare with `ConditionalProbability` at 1e-10.
15. **Unnamed tables:** sinks are named by their row number in `counts`.
16. **Reference testing.** SourceTracker2's `_sourcetracker.py` was run
    (numpy 2.5.5, pandas 2.3.3, scikit-bio 0.7.4) by
    `data-raw/expected/sourcetracker.py`: collapse (exact), conditional
    probabilities (max relative difference 4.6e-13), the schedule (exact,
    found by running `gibbs_sampler()` with a one-hot probability per pass),
    and `gibbs()` with 40 seeds per setting (default; 4 restarts × 3 draws
    with burnin 20, delay 5; leave-one-out) on 6 simulated sinks of 400
    reads from 3 sources of 2000 reads and an unknown source.
    **SourceTracker2's own leave-one-out fails** (`gibbs(sources, None)`):
    `sources[lambda x: x != source]` masks values instead of dropping the
    held-out row, so its sampler gets one source too many and
    `collate_gibbs_results()` raises. The loo fixture runs what it means to
    (`gibbs_sampler()` on `sources.drop(index=source)`, then
    `collate_gibbs_results(loo = TRUE)`). The docs' "Differences from
    SourceTracker2" say ecodive's `loo` works.
    **SourceTracker 1.0** (GPL, github.com/danknights/sourcetracker,
    `src/SourceTracker.r`) was run for its outputs only, by
    `data-raw/expected/sourcetracker.r` with `SOURCETRACKER_V1` set: 20 seeds
    of R's generator on the same simulated tables, `delay = 1`, other
    settings as SourceTracker2's defaults. Only its README, `example.r`'s
    calls and the two functions' argument lists were read. Its means are
    within 0.0038 of ecodive's (max z 2.1) and 0.0045 of SourceTracker2's;
    the test checks both within 0.02. (A first attempt to run it was refused
    by the session's permission checks; the user then asked for it, and it
    ran.)
