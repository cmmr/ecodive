# Rank plots, `plot()` of log ratios, and `log_ratio_features()` (Qurro-style)

**Priority:** 44. **Tier:** port (Qurro 0.9.0, BSD-3; its display rules,
not its code: Qurro is Python/JavaScript/Vega-Lite and nothing in it is
algorithmically hard). Read [README.md](README.md) first, and
[compositional.md](compositional.md) and [log_ratio.md](log_ratio.md).


## What it is

Morton et al. (2019) argue that relative data cannot say which taxa changed
in absolute terms, only how **ratios** between taxa changed, so the
interpretable output of a differential-abundance (DA) analysis is a
**ranking** of features (a reference frame), and the useful follow-up is the
**log ratio of top-ranked over bottom-ranked features**, looked at sample by
sample. Qurro (Fedarko et al. 2020) shows the two together: a **rank plot**
(features sorted by their "differential" or "feature loading") and a
**sample plot** (each sample's log ratio by a metadata field).

Three pieces:

| Piece | What |
|:------|:-----|
| (a) `type = 'rank'` | a rank plot on `dirmult_ttest()`, `ancombc()` / `ancombc2()` results and on ordinations with `loadings` (`pca()`, `rpca()`, `rpca_phylo()`, `rda()`, ...) |
| (b) `plot.ecodive_log_ratio()` | the sample plot: `log_ratio()` / `balance()` values by group (box + points) or against a numeric covariate |
| (c) `log_ratio_features()` | Qurro's "autoselection": the `n` top- and bottom-ranked features, as `numerator` / `denominator` for `log_ratio()` |

**Uptake:** Qurro 176 OpenAlex citations, 73 Europe PMC mentions 2022–26,
declining (0.33× normalised, tracking the QIIME 2 plugins); reference
frames 668 citations. The idea outlives the tool; a tool-agnostic R version
on ecodive's DA results and RPCA loadings is cheap. Moderate demand.


## Names

- **(a)** an argument value, `type = 'rank'`, as `plot.ecodive_test()` uses
  `type`. [da_plots.md](da_plots.md) (item 35) adds `type = c('volcano',
  'effect', 'ma', 'forest')` to the DA plot methods; this item appends
  `'rank'` to that list. If item 35 hasn't landed, add `type = c('volcano',
  'rank')` here, in the position item 35 specifies (second for
  `dirmult_ttest`; after `covariate` for `ancombc`, because
  `plot(fit, 'age')` passes `covariate` by position and is tested).
- **(b)** an S3 method, no new export besides the method.
- **(c)** `log_ratio_features()`: a companion of `log_ratio()` (README rule
  2), named for what it returns. Not `balance_select()` (reserved in
  [coda_lasso.md](coda_lasso.md) for a possible clean-room selbal) and not
  `autoselect()` (Qurro's UI term, meaningless out of context). No
  collisions in `NAMESPACE`.
- **`ancom()` is left out of (a) and (c), on purpose:** its `W` is a count
  of significant log ratios with no direction, so it defines no reference
  frame; its existing plot is already a ranked dot chart of `W`. The docs say
  so.


## Design

### (a) Rank plot

```r
plot.ecodive_dirmult_ttest(x, type = c('volcano', 'effect', 'ma', 'forest', 'rank'),
                           labels = 10, n_features = 30,
                           numerator = NULL, denominator = NULL, ...)
plot.ecodive_ancombc(x, covariate = NULL,
                     type = c('volcano', 'effect', 'ma', 'forest', 'rank'),
                     labels = 10, n_features = 30, conf_level = 0.95,
                     numerator = NULL, denominator = NULL, ...)
plot.ecodive_ord(x, groups = NULL, data = NULL, axes = c(1, 2), loadings = 10,
                 biplot = TRUE, centroids = TRUE, target = TRUE,
                 legend = 'topright', col = NULL, pch = NULL,
                 scaling = 'form', links = NULL, type = 'samples',
                 numerator = NULL, denominator = NULL, ...)
```

(`'effect'`, `'ma'`, `'forest'`, `n_features` and `conf_level` come from
[da_plots.md](da_plots.md), item 35; `scaling` and `links` from
[biplot_scaling.md](biplot_scaling.md), item 43. `numerator` /
`denominator` are used only by `type = 'rank'`.)

- **Ranking:** `dirmult_ttest`: `log2_fc`; `ancombc` / `ancombc2`:
  `log_fc` of `covariate` (default as today: the first non-intercept);
  ordinations: `loadings[, axes[[1]]]` (error if `x` has no `loadings`).
  Non-finite values (ancombc's `NaN`) are left out and counted in the
  subtitle.
- **Drawing** (internal `plot_rank_bars(effect, feature, num, den,
  reject, ...)` in `R/plot_tables.r`, distinct from the existing
  `plot_ranked()` dot chart): sort ascending (Qurro sorts ascending, x = "Feature
  Rankings"); bar \eqn{i} is `rect(i - 0.5, 0, i + 0.5, e[i], border = NA)`,
  one vectorised call; `abline(h = 0)`. x label "Feature rank (of m)"; y
  label as Qurro's "Differential: <column>" / "Feature Loading: <axis>":
  here "log2 fold change, <treatment> / <reference>", "Log fold change,
  <covariate>", "Loading on <axis name>". Title: the method's plot title
  (`ancombc_plot_name()`, "Dirichlet-multinomial t-test", the ordination's)
  + ": ranked features".
- **Colours:** with `numerator` / `denominator` (any selector form that
  `feature_set()` accepts, matched to the result's `feature` column or
  `rownames(loadings)`; overlap → `feature_set()`'s error): numerator bars
  `PLOT_HIGHLIGHT`, denominator `'royalblue3'`, others `'grey85'`, with a
  legend "Numerator (k)" / "Denominator (k)". Without them: features with
  `reject = TRUE` `'grey30'`, others `'grey85'` (ordinations: all
  `'grey50'`). Qurro's own colours (`#f00`, `#00f`, `#e0e0e0`) are not
  copied; only the idea.
- **Labels:** the `labels` selected features (else rejected features) with
  the largest |effect| are named at their bar ends with `plot_label()`;
  `0` for none.
- **`type` on ordinations:** `type = 'samples'` (default) is today's plot;
  `'rank'` the rank plot. **Any other value is passed on to
  `plot.default()`** as now, so `plot(ord, type = 'n')` keeps working (the
  docs say so). DA methods validate `type` strictly, as
  `plot.ecodive_test()` does (`test_plot_type()`).

### (b) `log_ratio()` / `balance()` results get a class

Both return `class = 'ecodive_log_ratio'` (one class; `attr(, 'method')`
says which). What that changes (checked in R 4.x):

- `print()` would show `attr(,"class")` and every attribute: add
  `print.ecodive_log_ratio()`, which prints a header ("Log ratio of 2 over
  3 features; pseudocount 0" / "Balance of ...") and then
  `print(unclass-and-strip(x))`, returning `x` invisibly.
- `x[i]`, `sort()`, `rev()`, `as.numeric()` return plain numerics
  (attributes and class dropped; names kept by `[`). Arithmetic (`2 * x`)
  keeps class and attributes. A data.frame column keeps the class and prints
  through `format.default()`. `is.numeric()` stays `TRUE`, so `groups`,
  `per_sample()`, `envfit()` and `lm()` are unaffected.
- Existing tests compare with `attributes(v) <- NULL` or `is.numeric()`,
  so only tests asserting `attributes()` in full need updating (check
  `test-log_ratio.r`, `test-smoke.r`). Add a smoke line for the class.

```r
plot.ecodive_log_ratio(x, groups = NULL, data = NULL, legend = 'topright',
                       col = NULL, pch = NULL, ...)
```

- `groups` / `data` resolved by `plot_groups(length(x), names(x), groups,
  data)` (the `plot.ecodive_ord()` rule; numeric → covariate).
- **Factor groups:** `boxplot(split(v, f), outline = FALSE, border =
  'grey60')`, then the samples as points, coloured by `ord_plot_look()`,
  spread horizontally **without R's RNG** (README: ecodive never touches
  `.Random.seed`): within a group of \eqn{m}, in value order, offset
  \eqn{0.6 \times ((i \cdot 0.618034) \bmod 1 - 0.5)}. Unused levels keep
  their colour and are left out, as `ord_plot_look()` does.
- **Numeric groups:** scatter, x = covariate (axis title = its column name),
  y = value. No fitted line or correlation: plot methods compute no new
  statistics (`R/plot.r`).
- **No groups:** one box with all samples.
- **Undefined values** (`NA` from zero sums): left out, counted in the
  subtitle "2 samples with an undefined log ratio not shown" (Qurro also
  drops them, `computeBalance()` returning null for sums ≤ 0).
- y label: `log(N / D)` with N, D the feature names joined by " + " when ≤ 3,
  else "k features"; balances: "Balance of N vs D". Title "Log ratio" /
  "Balance". Subtitle: "Numerator: a, b; denominator: c, d, e" (up to 5
  each via `name_list()`).
- Returns `x` invisibly.

### (c) `log_ratio_features()`

```r
log_ratio_features(x, n = 5L, covariate = NULL, axes = 1L)
```

- `x`: a `dirmult_ttest()`, `ancombc()` or `ancombc2()` result, an
  ordination with `loadings`, or a named numeric vector of effects.
  `covariate` as in `plot.ecodive_ancombc()`; `axes` one axis number
  (the glossary's `axes`, length 1 here).
- Drop non-finite effects; `o <- order(e)` (ties keep feature order, which
  is what Qurro's stable ascending sort then `slice(featureCt - n)` /
  `slice(0, n)` gives). `numerator` = the last `n` of `o`, largest effect
  first; `denominator` = the first `n`, smallest first.
- `n`: a whole number ≥ 1 with `2n` ≤ the number of features, else an
  error. (Qurro also takes a percentage; not offered: one meaning per
  argument.)
- `rpca_phylo()`: `loadings` has internal nodes (`"n12"`), which aren't
  columns of `counts`. Use tips only (rows in `attr(x, 'tree')$tip.label`),
  with a message; the docs show node log ratios via
  `log_ratio(attr(ord, 'node_counts'), ...)` from `node_counts = TRUE`.
- Returns `list(numerator =, denominator =)` (character), so
  `do.call(log_ratio, c(list(counts), sel))` works.
- **Docs warning:** choosing features by their DA ranking and then testing
  the same log ratio between the same groups is circular; use it to look,
  or select on one subset and test on another.


## Reference implementations

- **Qurro 0.9.0** (PyPI sdist; `LICENSE.txt`: BSD-3, "Copyright (c) 2018--,
  Qurro development team"; verified 2026-11-20). Read for behaviour only:
  `generate.py` `gen_rank_plot()` (ascending sort, row-number x,
  classification colours None / Numerator / Denominator / Both);
  `support_files/js/feature_computation.js` (`filterFeatures()` autoselect:
  `floor()` counts, negative input swaps top and bottom;
  `extremeFilterFeatures()`; `computeBalance()`: `null` when either sum ≤ 0);
  `display.js` (sample plot: metadata field on x, nominal or quantitative,
  optional boxplot and jitter). No code is ported, so no notice is needed;
  cite the paper (as [log_ratio.md](log_ratio.md) decided for `log_ratio()`).
- Songbird (BSD-3) produces the differentials Qurro ranks; not needed.


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- `log_ratio()` / `balance()` results get their `log_ratio` record in item
  34; this item adds the S3 class and reads the record in `plot()` (title
  "Log Ratio of 5 over 5 features", or the feature names when few).
- `log_ratio_features(x, n)` records where its features came from: the
  `log_ratio` record built from them gets `from = attr(x, 'da')` (or
  `attr(x, 'ord')`), so the strip plot's subtitle can say "top and bottom 5
  by ANCOM-BC2 log fold change".
- The rank plot (`type = 'rank'`) titles from the result's `da` or `ord`
  record and reads `log_base` for its axis label.

## Tests

- `log_ratio_features()` exactly: a named vector with ties (feature order
  wins), `NA`s dropped, `n` too large / 0 / 2.5 → errors; on
  `dirmult_ttest()`, `ancombc2()` (`covariate`), `pca()` (`axes = 2`),
  `rpca_phylo(ex_counts, ex_tree)` (tips only, message).
- `log_ratio()` / `balance()` have class `ecodive_log_ratio`; `print()`
  output (`expect_stdout`) shows no `attr(`; `x[1:2]` is plain numeric.
- `pdf(NULL)`: `invisibly()` for `type = 'rank'` on each DA class and on
  `pca()`, `rpca()`; with `numerator` / `denominator`; `plot(ord, type =
  'n')` still draws nothing and returns invisibly; `plot(lr)`, `plot(lr,
  'Species', iris)`, `plot(lr, 'Sepal.Length', iris)`, with `NA` values
  (subtitle counts them: test the internal subtitle helper's string).
- Errors: `type = 'rank'` on an ordination without loadings; overlapping
  `numerator` / `denominator`; unknown features; `type = 'x'` on DA results.
- `.Random.seed` unchanged after every new plot (the deterministic spread).


## Docs

- `plot.ecodive_dirmult_ttest`, `plot.ecodive_ancombc`, `plot.ecodive_ord`:
  `@param type`, `numerator`, `denominator`; a "Rank plots" paragraph citing
  Morton (2019) and Fedarko (2020).
- New Rd `plot.ecodive_log_ratio` (Plotting group) and `log_ratio_features`
  ("Compositional Data" → "Log ratios", beside `log_ratio`, `balance`).
- `log_ratio()` / `balance()` `@return`: the class and what it changes.
- A worked example: `fit <- ancombc2(...)`; `sel <- log_ratio_features(fit,
  n = 5)`; `plot(fit, type = 'rank', numerator = sel$numerator, denominator
  = sel$denominator)`; `plot(do.call(log_ratio, c(list(counts), sel)),
  'status', meta)`. NEWS.


## References

- Morton, J. T., Marotz, C., Washburne, A., et al. (2019). Establishing
  microbial composition measurement standards with reference frames.
  *Nature Communications*, 10, 2719. \doi{10.1038/s41467-019-10656-5}
- Fedarko, M. W., Martino, C., Morton, J. T., et al. (2020). Visualizing
  'omic feature rankings and log-ratios using Qurro. *NAR Genomics and
  Bioinformatics*, 2, lqaa023. \doi{10.1093/nargab/lqaa023}
- Martino, C., Morton, J. T., Marotz, C. A., et al. (2019). A novel sparse
  compositional technique reveals microbial perturbations. *mSystems*, 4,
  e00016-19. \doi{10.1128/mSystems.00016-19}
