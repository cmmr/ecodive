# `dirmult_ttest()`: Dirichlet-multinomial Welch t-tests on clr draws

**Priority:** see [TODO.md](TODO.md). **Tier:** port (scikit-bio 0.7.4,
`skbio/stats/composition/_dirmult.py`, BSD-3), with an optional
clean-room addition (ALDEx2-style scale uncertainty; ALDEx2 is
output-only). Read [README.md](README.md) first, and
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
- Write it from Nixon et al. (2024) and ALDEx2's **documentation**
  (ALDEx2 is GPL: no source). Fixtures from `ALDEx2::aldex.clr(gamma =)`
  can only be compared by distribution.
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
- Nixon, M. P., McGovern, K. C., Letourneau, J., David, L. A., Lazar, N.
  A., Mukherjee, S., & Silverman, J. D. (2024). Scale reliant inference.
  *arXiv*:2201.03616 (check for the journal version).
- Nearing, J. T., Douglas, G. M., Hayes, M. G., et al. (2022). Microbiome
  differential abundance methods produce different results across 38
  datasets. *Nature Communications*, 13, 342.
- Marsaglia, G., & Tsang, W. W. (2000). A simple method for generating
  gamma variables. *ACM Transactions on Mathematical Software*, 26,
  363–372.
