# `dirmult_ttest()`: Dirichlet-multinomial Welch t-tests on clr draws

**Priority:** 20. **Tier:** port (scikit-bio 0.7.4,
`skbio/stats/composition/_dirmult.py`, BSD-3), with an optional
scale-uncertainty addition ported from ALDEx3 (CRAN 1.4.0, **MIT**;
ALDEx2 itself is GPL and output-only). Read [README.md](README.md) first, and
[compositional.md](compositional.md) for the conventions the
compositional plans share.


## What it is

The ALDEx2 idea (Fernandes et al. 2013, 2014), in scikit-bio's form:

1. Treat each sample's counts as one draw from a multinomial, and sample
   plausible underlying proportions from the Dirichlet posterior
   (counts + a pseudocount). This carries the uncertainty of low counts
   and zeros forward, instead of fixing one pseudocount.
2. clr-transform each draw.
3. Run a Welch t-test per feature between two groups on each draw.
4. Average across draws.

The result is a per-feature effect size (log2 fold change of the clr, with
an interval) and an averaged p-value, which ANCOM doesn't give. It is one
of the most used compositional differential-abundance approaches
(Nearing et al. 2022 found ALDEx2 among the most consistent), and
scikit-bio's permissive implementation makes it portable.

**Name:** `dirmult_ttest()`, scikit-bio's name. Not `aldex()`: ALDEx2
differs in details (it reports expected Welch and Wilcoxon p-values, uses
its own Monte Carlo defaults and its own effect size), and its name
shouldn't be borrowed for a different procedure (README rule 6). The docs
say "ALDEx2-style".

scikit-bio's sibling `dirmult_lme()` (a linear mixed model per feature and
draw) isn't planned: it needs statsmodels' REML fitter (see
compositional.md).


## Design

### Signature

```r
dirmult_ttest(counts, groups, data = NULL, contrast = NULL, margin = 1L,
              pseudocount = 0.5, n_draws = 128L, p_adjust = 'holm',
              seed = 0L, cpus = n_cpus())
```

- `groups`, `data`: categorical, via `test_groups()`.
- **`contrast`:** which samples are compared, as `c(treatment, reference)`.
  - `NULL` (default): scikit-bio's default — treatment is the **first
    group in sorted order** (C-locale, as `np.unique`), reference is
    **all other samples pooled**. With two groups that is first vs second
    (note: the opposite of R's usual "second level vs first"; say so).
  - `c('a', 'b')`: group a vs group b; other samples are ignored.
  - `c('a', NA)`: group a vs all other samples pooled.
  - Named `contrast`, not `treatment` / `reference`: `reference` would
    clash with other meanings (see compositional.md), and `contrast` is the
    R idiom (DESeq2). The two-element order (treatment first) matches how
    the result reads: positive fold change = more in the treatment.
- `pseudocount = 0.5`: added to **every** cell before drawing, as
  scikit-bio (and ALDEx2's prior). Unlike the other compositional
  functions there's no `NULL` rule: the pseudocount is the Dirichlet prior
  here, not a zero patch. Must be > 0.
- `n_draws = 128L` (scikit-bio's `draws`; `n_<unit>` rule).
- `p_adjust`: as `ancom()`; applied once across features to the mean
  p-values.
- `seed`: PCG, as everywhere (README).

### Algorithm (scikit-bio, read from source)

For each draw d = 1..`n_draws`:

1. \eqn{G_{ij} \sim \mathrm{Gamma}(x_{ij} + \mathrm{pseudocount}, 1)},
   independently (a Dirichlet draw without the closure, which the clr
   makes unnecessary). scikit-bio uses `rng.standard_gamma(matrix)`, with
   the pseudocount already added to `matrix`.
2. clr: \eqn{\log G_{ij} - \overline{\log G_{i\cdot}}}.
3. Welch's t-test, treatment minus reference, per feature:
   \eqn{v_1 = s_t^2 / n_1}, \eqn{v_2 = s_r^2 / n_2} (variances with
   `ddof = 1`), \eqn{se = \sqrt{v_1 + v_2}},
   \eqn{\nu = (v_1 + v_2)^2 / (v_1^2/(n_1 - 1) + v_2^2/(n_2 - 1))},
   \eqn{t = \Delta / se}, \eqn{p = 2\,\mathrm{pt}(-|t|, \nu)}, and the
   95% interval \eqn{\Delta \pm \mathrm{qt}(0.975, \nu)\,se}.

Then, across draws:

- `t` = mean t; `log2_fc` = mean \eqn{\Delta} / ln 2;
- `ci_lower` = **minimum** lower bound / ln 2, `ci_upper` = **maximum**
  upper bound / ln 2 (scikit-bio's conservative envelope, not a quantile);
- `p_value` = mean p; `p_adjusted` = `p.adjust(p_value, p_adjust)` across
  features;
- `reject` = `p_adjusted <= 0.05` **and** both interval bounds strictly on
  the same side of 0. scikit-bio hard-codes 0.05; ecodive should too for
  parity, and document it (or add `alpha = 0.05`, as `ancom()` has: ask the
  user).

Groups of one sample give NaN (variance undefined), as scikit-bio.

### Random numbers and reproducibility

- **Exact parity with scikit-bio is impossible**: it draws from NumPy's
  PCG64 Gamma sampler. Match the deterministic parts exactly instead (see
  Tests) and the random parts by distribution.
- ecodive needs a **Gamma sampler on PCG** (none exists yet). Use
  Marsaglia & Tsang (2000) for shape ≥ 1, with the standard
  \eqn{U^{1/a}} boost for shape < 1 (counts of 0 give shape 0.5), and a
  normal sampler (e.g. Marsaglia polar or ziggurat) on PCG. Put them in a
  small shared file (`src/rng_dist.c`) for reuse.
- Determinism for any `cpus`: seed one PCG stream per draw
  (`pcg32_srandom_r(&rng, seed, draw)`), split work by draws (or blocks of
  draws), and reduce the per-draw results in draw order. Within a draw the
  samples are filled in a fixed order from that one stream.

### Performance

n × D Gamma draws per draw; 128 draws × 500 samples × 5,000 features is
3.2e8 draws, a few seconds threaded. Keep one n × D buffer per thread,
reused across its draws; never store all draws. Accumulate per-feature
running sums (t, Δ, p) and running min/max of the interval bounds, per
draw block, combined in block order.

### Optional: scale uncertainty (clean-room addition; ask the user first)

Nixon et al. (2024) show that ALDEx2-style tests implicitly assume the
clr's denominator (the geometric mean) is the same across conditions, and
that adding uncertainty about the **scale** (total load) to each draw
controls false positives. ALDEx2 exposes it as `aldex.clr(gamma = )`.

- Plan it as an argument (it's an option of the same method in ALDEx2):
  `scale_sd = 0` (no scale model, scikit-bio's behaviour). Not `gamma`,
  which `phate()` uses for something else, and not `scale`, which means
  standardisation (parameter_names.md).
- The model (Nixon, Gloor & Silverman 2025, Methods; ALDEx2 vignette):
  per Monte Carlo draw, \eqn{\log_2 W^\perp_n = -\overline{\log_2
  W^\parallel_{\cdot n}} + \Lambda x_n} with \eqn{\Lambda \sim N(0,
  \gamma^2)} drawn **once per draw and shared within a condition**
  (\eqn{x_n} the 0/1 condition code); \eqn{\gamma = 0} gives back the
  clr. Recommended \eqn{\gamma = 0.5}. Base-2 logs: convert to natural
  logs carefully.
- **Port from ALDEx3** (MIT; `clr.sm(X, logComp, gamma = 0.5)` and
  related scale-model functions): read its code, credit it in
  `inst/COPYRIGHTS`. Check whether it draws the noise per condition or per
  sample (ALDEx2's numeric-`gamma` path is undocumented on this) and
  follow ALDEx3. ALDEx2's own code stays unread (GPL); fixtures from
  `ALDEx2::aldex.clr(gamma =)` / ALDEx3 compare by distribution.
- Ship the port first; add this in a second commit.

### Result

The shared per-feature structure (compositional.md):

- `data.frame` of class `c('ecodive_dirmult_ttest', 'data.frame')`, one
  row per feature in table order.
- Columns: `feature`, `t`, `log2_fc`, `ci_lower`, `ci_upper`, `p_value`,
  `p_adjusted`, `reject`. (scikit-bio: `T-statistic`, `Log2(FC)`,
  `CI(2.5)`, `CI(97.5)`, `pvalue`, `qvalue`, `Signif`.)
- Attributes: `contrast` (the two labels; `NA` for pooled), `n` (samples
  per side), `n_draws`, `pseudocount`, `p_adjust`, `seed`.


## Reference implementations

- **Port from:** scikit-bio 0.7.4, `dirmult_ttest()`, `_dirmult_draw()`,
  `_welch_draw_stats()` (`_dirmult.py`), `_check_trt_ref_groups()` and
  `_adjust_pvalues()` (`_utils.py`), BSD-3.
- **Fixtures** (`data-raw/expected/dirmult_ttest.py`):
  - **The deterministic kernel:** give the Welch step fixed clr matrices
    (no draws) and store its per-draw t, Δ, p and bounds; then the
    aggregation over a fixed stack of such matrices. Expose an internal
    entry point (`dirmult_ttest_kernel()`) for the test. Match at 1e-10.
  - **The whole function by distribution:** for a table with large counts
    the Dirichlet noise is small, so `log2_fc` agrees with scikit-bio to
    a loose tolerance (e.g. 1e-2); `reject` agrees on clear-cut features.
  - **The Gamma sampler:** moments and a Kolmogorov–Smirnov check against
    `qgamma` for shapes 0.5, 1, 7.5, 1000.
- **Output-only cross-check:** `ALDEx2::aldex()` with `test = "t"` on the
  same table: effects should correlate closely, not match.


## Tests

- As above; plus identical results for `cpus` 1, 2 and 4 and for the same
  `seed`; different for another `seed`.
- `contrast` forms: default, explicit pair, pooled (`NA`); errors for
  unknown labels, identical labels, a side with no samples.
- A group of one sample gives NaN columns, not an error (as scikit-bio);
  say so in a warning.
- A smoke-test line (small `n_draws`).


## Docs

- Rd page in "Compositional Data", subtitle "Differential abundance".
- Explain the Dirichlet draws and why the interval is a min/max envelope;
  note the default contrast direction; point to `ancom()`.
- NEWS.


## References

- Fernandes, A. D., Macklaim, J. M., Linn, T. G., Reid, G., & Gloor, G. B.
  (2013). ANOVA-like differential expression (ALDEx) analysis for mixed
  population RNA-Seq. *PLoS ONE*, 8, e67019.
- Fernandes, A. D., Reid, J. N. S., Macklaim, J. M., McMurrough, T. A.,
  Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of
  high-throughput sequencing datasets. *Microbiome*, 2, 15.
- Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating
  scale uncertainty into differential expression analyses using ALDEx2.
  *Genome Biology*, 26, 139. doi:10.1186/s13059-025-03609-3
- Nixon, M. P., McGovern, K. C., Letourneau, J., David, L. A., Lazar, N.
  A., Mukherjee, S., & Silverman, J. D. (2024). Scale reliant inference.
  *arXiv*:2201.03616.
- Nearing, J. T., Douglas, G. M., Hayes, M. G., et al. (2022). Microbiome
  differential abundance methods produce different results across 38
  datasets. *Nature Communications*, 13, 342.
- Marsaglia, G., & Tsang, W. W. (2000). A simple method for generating
  gamma variables. *ACM Transactions on Mathematical Software*, 26,
  363–372.


## Decisions made while implementing (2026-11-19)

Asked of the user at the start:

- **Scale model:** yes, as `scale_sd`, in a second commit, with default
  **0.5** (ALDEx3's and Nixon et al.'s recommendation), not 0: the default
  result therefore differs from scikit-bio's by design; `scale_sd = 0`
  gives scikit-bio's procedure. (Details below.)
- **Default contrast:** scikit-bio's rule. Treatment = the first group
  (`sort(method = 'radix')`, as `numpy.unique`; a factor's first level),
  reference = all other samples pooled. With two groups that is first vs
  second, the reverse of R's usual direction; the docs say so.
  `c(a, b)` leaves other groups out; `c(a, NA)` pools.
- **Result:** the planned data frame (`feature`, `t`, `log2_fc`,
  `ci_lower`, `ci_upper`, `p_value`, `p_adjusted`, `reject`) plus **both**
  `alpha = 0.05` (reject's cutoff on `p_adjusted`) and `conf_level = 0.95`
  (each draw's interval, quantile `(1 + conf_level) / 2`, which is exactly
  `0.975` for the default). Attributes `contrast` (named `treatment`,
  `reference`; `NA` reference when pooled), `n` (named, per side),
  `n_draws`, `pseudocount`, `alpha`, `conf_level`, `p_adjust`, `seed`. No
  print method.
- **Samplers:** Marsaglia–Tsang Gamma on the log scale, Marsaglia polar
  normals, 52-bit uniforms, and **one PCG stream per (draw, sample)**
  rather than per draw (see below).

Decided along the way:

- **`src/rng_dist.c`** (shared, MIT, ecodive's own): `rng_t` (a PCG state
  plus the polar method's spare normal), `rng_seed(r, seed, stream)`,
  `rng_unif()` (two PCG outputs, 26 bits each, `(k + 0.5) * 2^-52`: never 0
  or 1, exact), `rng_norm()` (polar; the pair's second value is kept for
  the next call on the same stream), `rng_log_gamma()` (Marsaglia–Tsang for
  shape >= 1, returning `log(d * v)`; shape < 1 as `log G(a + 1) + log(U) /
  a`, so tiny pseudocounts such as 1e-4 stay finite where NumPy's
  `log(standard_gamma())` would give `-Inf`), `rng_gamma()`. The internal
  `C_rng_dist` exposes them for the tests. Reusable as is by later items.
- **Streams:** sample s of draw d uses stream `d * (n + 1) + s` (n = samples
  in the table, s its index in the table); features are drawn in table
  order. So a draw's samples fill in parallel into one shared
  (samples in the contrast) × D buffer, draw d doesn't depend on
  `n_draws` (more draws extend the sequence), samples outside the contrast
  aren't drawn and change nothing, and sparse and dense input agree
  exactly. Stream `d * (n + 1) + n` is left free for the scale model.
- **Per draw:** phase A (threads over samples): log-Gamma row, then the clr
  (sum in feature order / D). Phase B (threads over blocks of 64
  features): means, then variances as sums of squared deviations
  (`ddof = 1`), in sample order, and Welch's difference, se and dof in
  `_welch_draw_stats()`'s order of operations. Then, on the main thread,
  `t`, `pt(-|t|, dof)`, `qt(prob, dof)` (R's are not thread-safe when they
  warn; skipped when dof isn't positive) and the running sums / NaN-
  propagating min / max in draw order, then the averages in scikit-bio's
  order (`/ n_draws`, `* 2 / n_draws` for p, then `/ log(2)`). The main
  thread's `pt()`/`qt()` cost about 0.6 s of the 4.8 s below.
- **Single-sample side:** a warning, NaN in every column but `log2_fc`
  (which is defined), `p_adjusted` NaN (left out of the adjustment, as
  `p.adjust()` and scikit-bio do), `reject` FALSE.
- **Input:** any table type (`validate_counts()`), non-negative, finite,
  at least 2 features; samples of all zeros are allowed (the prior makes
  them uniform; scikit-bio refuses them). `pseudocount` must be > 0 (no
  `NULL` rule; `log_ratio_input()` isn't used because its zero rule doesn't
  apply). `p_adjust` is any of `stats::p.adjust.methods` (by prefix, via
  `pairwise_method()`), applied with `stats::p.adjust()`.
- **The test hook** is the internal `dirmult_ttest_details(...,
  gamma_draws =)`: the whole function on given Gamma draws (an n × D ×
  n_draws array), also returning each draw's diff, se, dof, t, one-sided
  tail and bounds. The fixture script replaces scikit-bio's generator with
  one that hands out the same NumPy draws.
- **Reference results** (scikit-bio 0.7.4, NumPy 2.5.5, SciPy 1.18.1), on
  the same draws, 8 cases × 5 `p_adjust`: largest relative difference
  3.4e-14 over `t`, `log2_fc`, the bounds, `p_value` and `p_adjusted`;
  `reject` identical everywhere. Per draw: diff, se, dof and t to 2.2e-14,
  the tail probability to 7.5e-15, the bounds to 3.3e-13 relative to the
  bound (they are differences near 0). R's `pt()`/`qt()` against SciPy's
  `stdtr()`/`stdtrit()` on a grid: 1.3e-14 / 8.6e-15. With their own draws
  (deep table, 20 samples × 12 features, 200,000 reads each), `log2_fc`
  differs from scikit-bio's by at most 1.4e-3 and `reject` agrees.
- **Samplers by distribution:** moments and Kolmogorov–Smirnov tests against
  `punif`, `pnorm` and `pgamma` (shapes 0.5, 1, 7.5, 1000, and the log draws'
  mean against `digamma()`, including shape 1e-3) all pass; checked further
  by hand with 1e7 draws over 5 seeds.
- **Timing** (4 cores, -O2): 500 samples × 5,000 features (22% zeros), 128
  draws: 16 s on 1 core, 4.8 s on 4. About 50 ns per Gamma draw.

### Scale model (`scale_sd`, second commit)

- **Port:** ALDEx3 1.4.0's `clr.sm()` (CRAN, MIT; its code read, credited in
  `inst/COPYRIGHTS`, the `src/dirmult.c` and `R/dirmult_ttest.r` headers
  and `Authors@R` as "ALDEx3 authors", ctb/cph). ALDEx2's code stayed
  unread. `clr.sm()` draws **one N(0, gamma²) coefficient per row of the
  design per Monte Carlo draw** (`logScale = -colMeans(logComp) + t(X) %*%
  Lambda`), so samples with the same design share the shift: per
  condition, not per sample. Followed: for a contrast the design is
  intercept + treatment indicator; the intercept's shift cancels from
  every difference and isn't drawn, and the treatment's Λ (log2 units) is
  added, times log 2, to every feature's difference in that draw. Within-
  side variances don't change. A pooled reference is one condition (one
  shift), not one per group.
- **Stream:** Λ for draw d is the first polar normal of stream
  `d * (n + 1) + n` (the stream the port left free), so it doesn't depend
  on the Gamma draws (the same with `gamma_draws` given), the Gamma draws
  are the same with and without the scale model, and it is the same for
  any `cpus`.
- **Default 0.5** (the user's choice; ALDEx3's documented example value and
  Nixon et al.'s recommendation). So `dirmult_ttest()`'s default result
  differs from scikit-bio's; the tests compare with scikit-bio at
  `scale_sd = 0`, and the docs and the example say so.
- **Not mirrored:** reversing the contrast with the same seed doesn't
  negate Λ (same stream), so with `scale_sd > 0` the two results are
  equivalent in distribution, not exact mirror images (they are at
  `scale_sd = 0`, which the tests check).
- **The envelope widens with `n_draws` under the scale model** (found while
  testing, kept as scikit-bio's rule and documented): the interval is the
  min/max over the draws' intervals, so it reaches the most extreme Λ
  drawn, and `reject` gets stricter as draws are added. On the deep table
  (true log2 changes about ±1.1 to 1.5 and -0.22), rejections at
  `scale_sd` 0 / 0.25 / 0.5 / 1 were 12/12/12/12, 4/4/4, 4/4/0 and 0/0/0
  for 32/128/512 draws. The averaged p-values converge, slowly (rare large
  shifts drive them). The docs suggest filtering on `p_adjusted` alone (as
  ALDEx2 does) for decisions that don't depend on `n_draws`. A future
  option could be a quantile interval, or an interval from the averaged
  draws, but neither is scikit-bio's or ALDEx3's.
- **Cross-check with ALDEx3** (output only, by distribution;
  `aldex(t(x), ~cond, nsample = 2000, scale = clr.sm, gamma =, test = 't')`
  against 2000 ecodive draws): the per-draw spread of each feature's
  difference matches (ratio of SDs, median 0.995–1.013, range 0.947–1.045,
  over three tables at `gamma` 0 and 0.5; Monte Carlo error about ±3%), the
  mean log2 fold changes differ by at most 0.064 (sparse table; 0.002 on the
  deep one), and -log p-values correlate 0.95–0.99 (ALDEx3's are OLS t, not
  Welch).
- **Cross-check with ALDEx2 1.34.0** (Bioconductor, GPL: output only;
  `aldex(test = 't', effect = TRUE, mc.samples = 512)` against
  `scale_sd = 0`, 512 draws, on four of the fixture tables): ALDEx2's
  `diff.btw` (a median of differences, not a difference of means)
  correlates 0.989–1.000 with `log2_fc` (largest difference 0.56, on the
  sparse three-group table; 0.002 on the deep one), the signs agree on 55
  of 56 features, and -log of ALDEx2's expected Welch p (`we.ep`)
  correlates 0.993–0.999 with -log `p_value`.
- **Not changed:** the draft article `vignettes/articles/compositional.Rmd`
  still lists `dirmult_ttest()` (and `ancom()`) as planned and says the
  scale model is "planned for `dirmult_ttest()`"; its differential-abundance
  chunk is `eval = FALSE`. Left for the article's own update, as item 19 did.
