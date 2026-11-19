# `ancom()`: analysis of composition of microbiomes

**Priority:** 19. **Tier:** port (scikit-bio 0.7.4,
`skbio/stats/composition/_ancom.py`, BSD-3). Read [README.md](README.md)
first, and [compositional.md](compositional.md) for the conventions the
compositional plans share.


## What it is

ANCOM (Mandal et al. 2015) finds features whose abundance differs between
groups **without assuming the total is constant**. For each feature it
tests every log ratio with every other feature for a group difference; a
feature that really changed makes most of its log ratios differ. Its
statistic **W** is the number of significant log ratios a feature is in,
and features with a large W are declared differentially abundant.

It was the first widely used compositionally aware differential-abundance
method and is still in QIIME 2 (`composition ancom`). ANCOM-BC / ANCOM-BC2
(Lin & Peddada 2020, 2024) are its successors, planned separately as
`ancombc()` / `ancombc2()` (see [ancombc.md](ancombc.md)).

The work is \eqn{D(D-1)/2} independent tests, which suits ecodive's
threaded C well: scikit-bio's version is slow on thousands of features.


## Design

### Signature

```r
ancom(counts, groups, data = NULL, margin = 1L, pseudocount = NULL,
      alpha = 0.05, tau = 0.02, theta = 0.1, p_adjust = 'holm',
      sig_test = 'anova', cpus = n_cpus())
```

- `groups`, `data`: categorical grouping, handled by `test_groups()` as in
  the tests (at least two groups; error if every sample is its own group).
  **Group order:** scikit-bio sorts groups with `np.unique`; use C-locale
  sorting (`sort(method = 'radix')`) to match.
- `pseudocount`: scikit-bio requires strictly positive values and adds
  nothing. ecodive applies the shared `pseudocount` rule (see
  compositional.md), so a table with zeros works, with the usual warning.
  Fixtures use a zero-free table, or scikit-bio on `x + pc`.
- `alpha`, `tau`, `theta`: as scikit-bio; each strictly between 0 and 1.
- `p_adjust`: `'holm'` (default), `'bonferroni'`, `'BH'`, `'BY'`, or
  `'none'`, matched with the pairwise tests' helper (`pairwise_method()`).
- `sig_test`: the per-log-ratio test.
  - `'anova'` (default; scikit-bio's `"f_oneway"`): one-way ANOVA, which
    is Student's t-test for two groups.
  - `'kruskal'`: Kruskal-Wallis (scikit-bio accepts any `scipy.stats`
    name; this is the other one users pick, and the one the original ANCOM
    paper used).
  - These are options of one method in every reference, so they stay an
    argument (README rule), not separate functions.
- Not planned: scikit-bio's `percentiles` table (per-group percentiles of
  the raw table, a reporting aid; users can call `quantile()`), and ANCOM
  II's covariate adjustment / structural-zero detection (Kaul et al. 2017),
  which scikit-bio doesn't have.

### Algorithm (scikit-bio)

1. \eqn{L = \log(x)} (after the pseudocount).
2. For every pair i < j, the log ratio \eqn{L_i - L_j} per sample, and the
   test across groups. Symmetric p-value matrix, 0 on the diagonal.
3. **Adjust each feature's row with the diagonal 0 still in it** (a quirk
   of scikit-bio). For Holm this equals Holm over the other D − 1 values;
   for Bonferroni, BH and BY it doesn't. The R equivalent for all methods
   is `p.adjust(c(0, row_without_diagonal), method)[-1]`. NaN p-values
   leave the family and count as not significant.
4. Diagonal set to 1. \eqn{W_i = \#\{j : p_{ij} < \alpha\}} (strict).
5. The cutoff rule, exactly:

   ```text
   c_start = max(W) / D
   if c_start < theta: reject nothing
   else:
     cutoff   = c_start - c(0.05, 0.10, 0.15, 0.20, 0.25)
     prop_cut = for each cutoff, mean(W > D * cutoff)            # strict
     dels     = abs(prop_cut - c(prop_cut[-1], prop_cut[1])); dels[5] = 0
     if      dels[1] <  tau & dels[2] < tau & dels[3] < tau: nu = cutoff[2]
     else if dels[1] >= tau & dels[2] < tau & dels[3] < tau: nu = cutoff[3]
     else if dels[2] >= tau & dels[3] < tau & dels[4] < tau: nu = cutoff[4]
     else nu = cutoff[5]
     reject = W >= nu * D
   ```

   (0-based indices in scikit-bio translated to R's 1-based.)

### C

- For each pair, compute the test statistic from the differences
  \eqn{L_i - L_j} directly (per-group means, then two-pass sums of
  squares): O(n D²) in all, parallel over blocks of pairs. Don't take the
  shortcut through a grouped Gram matrix (\eqn{s_{ii} + s_{jj} - 2s_{ij}}):
  it cancels badly for nearly proportional pairs and won't match
  scikit-bio at 1e-10 (the same reasoning as in
  [proportionality.md](proportionality.md)).
- R's C API isn't thread-safe, so worker threads store only the statistic
  (F or H); the p-values (`pf()` / `pchisq()` with the known degrees of
  freedom) are computed afterwards on the main thread, vectorised.
- Kruskal-Wallis needs ranks of each difference: O(n log n) per pair,
  directly.
- Memory: the D × D p-value matrix (8 D² bytes). Only row counts (W) are
  needed after adjustment, but Holm/BH need the whole row: process one
  row block at a time from the packed triangle; don't keep a dense D × D
  unless D is modest.

### Result

A per-feature table, the shared structure in compositional.md:

- `data.frame` of class `c('ecodive_ancom', 'data.frame')`, one row per
  feature in table order.
- Columns: `feature`, `W` (integer), `reject` (logical; scikit-bio's
  `Signif`).
- Attributes: `groups` (the levels, in the order used), `sig_test`,
  `p_adjust`, `alpha`, `tau`, `theta`, `cutoff` (the chosen \eqn{\nu}, or
  `NA` if nothing was rejected because `c_start < theta`), `pseudocount`.
- No print method (as `envfit()` and `simper()`).


## Reference implementations

- **Port from:** scikit-bio 0.7.4, `ancom()` and `_log_compare()` in
  `skbio/stats/composition/_ancom.py`; `_adjust_pvalues()` and
  `_check_grouping()` in `_utils.py` (BSD-3). The behaviour above was read
  from that source.
- **Fixtures:** `data-raw/expected/ancom.py` runs scikit-bio; commit the
  text output under `data-raw/expected/ancom/`. Include scikit-bio's own
  test cases: `table1` with `p_adjust=None` gives
  `W = [5,5,2,2,2,2,2]`, `Signif = [T,T,F,F,F,F,F]`; the docstring example
  with Holm gives `W = [0,4,0,1,1,0,1]`, only b2 significant.
- Cover: two groups and three; `'anova'` and `'kruskal'`; each
  `p_adjust`; a case with `c_start < theta`; constant log ratios (zero
  within-group variance; check what scipy returns and match it).


## Tests

- Fixtures exact for `W` and `reject`; p-values (exposed only in tests via
  an internal function) at 1e-10.
- Invariant to multiplying any sample by a constant.
- Identical for `cpus` 1, 2 and 4.
- Errors: one group; all-singleton groups; `alpha`/`tau`/`theta` outside
  (0, 1); unknown `sig_test` or `p_adjust`; NA in `groups`.
- A smoke-test line.


## Docs

- Rd page in "Compositional Data", subtitle "Differential abundance".
- Say that ANCOM gives no effect size or p-value per feature, only W and a
  decision; point to `dirmult_ttest()` for effect sizes and to `ancombc()` /
  `ancombc2()` for bias-corrected estimates.
- Document the row-adjustment quirk and that ecodive keeps it to match.
- NEWS.


## References

- Mandal, S., Van Treuren, W., White, R. A., Eggesbø, M., Knight, R., &
  Peddada, S. D. (2015). Analysis of composition of microbiomes: a novel
  method for studying microbial composition. *Microbial Ecology in Health
  and Disease*, 26, 27663.
- Kaul, A., Mandal, S., Davidov, O., & Peddada, S. D. (2017). Analysis of
  microbiome data in the presence of excess zeros. *Frontiers in
  Microbiology*, 8, 2114.
