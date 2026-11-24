# Rank plots, `plot()` of log ratios, and `log_ratio_features()` (Qurro-style)

**Priority:** 46. **Tier:** port (Qurro 0.9.0, BSD-3; its display rules,
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
  `type`. [da_plots.md](da_plots.md) (item 28) adds `type = c('volcano',
  'effect', 'ma', 'forest')` to the DA plot methods; this item appends
  `'rank'` to that list. If item 28 hasn't landed, add `type = c('volcano',
  'rank')` here, in the position item 28 specifies (second for
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
[da_plots.md](da_plots.md), item 28; `scaling` and `links` from
[biplot_scaling.md](biplot_scaling.md), item 44. `numerator` /
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


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- `plot.ecodive_log_ratio()` reads the attributes `log_ratio()` /
  `balance()` already set: `method` for the title ("Log ratio" /
  "Balance"), `numerator` and `denominator` for the y label and subtitle,
  as in (b). The input is a count table, so no `display`.
- `log_ratio_features()` returns a plain list; the log ratio built from it
  doesn't say where its features came from.
- The rank plot (`type = 'rank'`) takes the method's name from the
  registry by `method` (`DA_TESTS` for DA results) or, for ordinations,
  from `display$title`, and reads `log_base` for its axis label.

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


## Decisions made while implementing (item 46, 2026-11-24)

Asked of the user at the start (all three answered with the recommended
option):

1. **Every DA result with a signed effect** gets `type = 'rank'` and is
   accepted by `log_ratio_features()`, not only the plan's
   `dirmult_ttest()` / `ancombc()` / `ancombc2()`: also `da_ttest()`,
   `da_wilcoxon()` (`type = c('volcano', 'ma', 'rank')`), `da_lm()` and
   `linda()`. The effect column is read from `DA_TESTS$effect` by the
   result's `method` (`rank_source()`), so a new method needs no special
   case in `log_ratio_features()` if its registry row names a signed
   effect; its own `plot()` method still has to add `'rank'` and pass
   `num` / `den` to `da_plot()`. `ancom()` (no effect column) is refused:
   "ancom() results have no signed effect, so they rank no features into a
   reference frame." **`lefse()`** (item 41, merged while this was being
   built) registered `lda_score`, which is a magnitude (the class is in
   `enriched`), so `DA_TESTS` gained a **`signed`** column (`FALSE` for
   `ancom` and `lefse`, listed by `list_da_tests()`), and `rank_source()`
   refuses unsigned effects with the same message. `plot.ecodive_lefse()`
   gets no `'rank'` type.
2. **Arithmetic drops the class** (the plan said it keeps it):
   `Ops.ecodive_log_ratio` and `Math.ecodive_log_ratio` return a plain
   named numeric vector, as `[`, `sort()` and `rev()` already do, so `-lr`
   or `lr - mean(lr)` never plots or prints as `log(N / D)`. This follows
   README's "arithmetic drops the label" rule for distance matrices.
3. **No Vega from npm.** Qurro's own `feature_computation.js` runs under
   the container's Node.js (autoselection and `computeBalance()`);
   rank-plot positions are tested against a direct stable ascending sort.

Made while implementing:

4. **Files.** The rank-plot code is in a new `R/plot_rank.r` (not
   `R/plot_tables.r`, which other sessions were editing):
   `rank_source()` (the ranking a result defines), `da_rank_covariate()`,
   `rank_layout()` (the bars' coordinates, classes, colours and labels, as
   a data frame that the tests read), `plot_rank_bars()`,
   `rank_selection()`, `da_plot_rank()`, `ord_plot_rank()`, `rank_only()`.
   The class methods and the sample plot are in a new `R/plot_log_ratio.r`;
   `log_ratio_features()` is in `R/log_ratio.r`.
5. **Titles** follow item 28's pattern, "Rank plot of ANCOM-BC2" /
   "Rank plot of Dirichlet-Multinomial t-test" (via `da_plot_title(x,
   'rank')`) and "Rank plot of <the ordination's title>", rather than the
   plan's "<method>: ranked features".
6. **Axis labels and subtitle.** The y label is the effect label the
   method's other DA plots use ("Log2 fold change", "Log fold change
   (difference in means, clr)", ...); the comparison or covariate goes in
   the subtitle with the rejection count, as in the other DA plots
   (`da_plot_sub()`), plus "2 features with no estimate not shown".
   Ordinations: "Loading on PC1"; x: "Feature rank (of m)".
7. **One covariate rule.** `da_rank_covariate()` picks the default
   covariate for the plot methods and `log_ratio_features()` alike:
   ANCOM-BC's first non-intercept covariate, else `da_term_covariate()`
   (the term of interest's first column). It replaces the inline copies in
   `plot.ecodive_ancombc()`, `plot.ecodive_da_lm()` and
   `plot.ecodive_linda()` (same results; the suite is unchanged).
8. **Highlighting.** Without a selection, rejected features `'grey30'`,
   others `'grey85'`; ordinations all `'grey50'` (`PLOT_MUTED`). With
   `numerator` / `denominator`: `PLOT_HIGHLIGHT` and `'royalblue3'`, the
   rest `'grey85'`, legend "Numerator (k)" / "Denominator (k)" in the
   top-left corner (empty in an ascending rank plot); on ordinations,
   `legend` moves it only when given. Sensitivity analyses'
   `reject_robust` isn't shown (only `reject`).
9. **Labels.** DA plots: the `labels` highlighted features (selected, else
   rejected) with the largest |effect|. Ordinations have no `labels`
   argument (the plan's signature has none), so `loadings` (default 10)
   names the features with the largest |loading| (among the selection,
   when given). The y range is padded 5% so the names fit.
10. **`numerator` / `denominator` with another `type`** is an error
    ("are for `type = \"rank\"`"), on DA results and ordinations.
11. **Messages.** `feature_set()` gained `src` (default `` `counts` ``),
    so a rank plot's selection errors name "the result" instead.
12. **`log_ratio_features()` details.** `covariate` on a result without a
    `covariate` column, and `axes` (given explicitly) on anything but an
    ordination, are errors rather than ignored. An unnamed vector gives
    integer indices. The vector form must be a plain numeric vector (an
    `ecodive_log_ratio` or other classed object is refused). Numerator:
    largest first; denominator: smallest first (Qurro returns both in
    ascending order; the sets and the order match after reversing the
    top). Ties keep the result's order (R's `order()` is stable, as is
    V8's `Array.prototype.sort`).
13. **rpca_phylo.** `log_ratio_features()` ranks tips only, with the
    message "Ranking the 6 tips only: the 5 internal nodes of rpca_phylo()
    aren't features of the count table."; the rank plot shows every row of
    `loadings`, nodes included (they are features of `node_counts`).
14. **print()**: the header "Log ratio of 2 over 3 features; pseudocount
    0" (or "Balance of ..."), "Numerator:" and "Denominator:" lines (up to
    5 names), "2 samples are NA (a zero sum)" when any, then the values as
    a plain named vector (a per-sample variable prints in full, like any
    numeric vector).
15. **Sample plot.** Boxes with `col = NA`, `border = 'grey60'`, `boxwex
    = 0.6`, drawn with `add = TRUE` over a `plot.default()` frame, so `...`
    reaches `plot.default()` as everywhere else. Factor groups: the
    plan's legend (`legend = 'topright'`) and the group names on the x
    axis. Numeric groups: plain points (no gradient: x already shows the
    variable). No groups: one "All samples" box with grey points. The
    subtitle is one line ("...; 2 samples with an undefined log ratio not
    shown"); the y label parenthesises a side of several features,
    `log((a + b) / c)`.
16. **No `@family`** tags: no family existed for these pages, and one
    would add a lone "Other" section.
17. **Qurro drops features with no counts** in the table before its rank
    plot; ecodive ranks every feature in the result (the fixture tables
    have none all-zero, so the comparison is like for like; Qurro's
    `process_input()` keeps every ranked feature, in the ranking's order,
    in all 26 cases).
18. **References** (`data-raw/expected/log_ratio_plots.r` (two modes),
    `.py`, `.js`; 26 cases on two tables: 13 synthetic rankings with and
    without ties, and the effects / loadings of `dirmult_ttest()`,
    `ancombc2()` (two covariates), `da_ttest()`, `da_wilcoxon()`,
    `da_lm()`, `linda()`, `pca()` (two axes) and `rpca()`):
    - Qurro 0.9.0's autoselection (`filterFeatures()`, "autoLiteralTop" /
      "autoLiteralBot"): **identical** numerator and denominator in all 26
      cases, ties included.
    - `qarcoal()`: max abs difference **3.5e-18** (one value, `log(45 /
      46)`, where numpy's `log()` and C's differ by an ulp; all others
      bit-identical); `NA` exactly where Qurro drops a sample.
    - The sample plot's `computeBalance()` (`log(N) - log(D)` on
      display.js's sums): max 6.7e-16.
    - `cpus` 1, 2, 4 give identical selections from `dirmult_ttest()`;
      sparse and dense input give identical selections (`pca()`) and log
      ratios.
19. **Nothing is ported from Qurro** (its display rules only; the Node
    driver in `data-raw/` loads Qurro's file at fixture time and isn't
    built into the package), so Qurro still gets no `Authors@R` or
    `inst/COPYRIGHTS` entry; the docs cite the paper.


## Outcome (from the implementing session's final report)

**Item 46, done 2026-11-24.** Commits `32f2e17` (work), `9ead0e5` (merge
with items 41, 44 and 45: `DA_TESTS$signed`, `plot.ecodive_ord()`'s
arguments `scaling`, `links`, `type`, `numerator`, `denominator`), and the
tick.

- Tests on the merged tree: 28,140 at home, 406 in CRAN mode (before this
  item: 26,600 and 396; the new file has 319). `R CMD check --no-manual
  --ignore-vignettes`: 1 NOTE (installed size 5.5 MB), as before.
- Qurro 0.9.0, 26 cases: autoselection (`feature_computation.js` under
  Node.js) identical in every case, ties included; `process_input()` keeps
  every ranked feature in the ranking's order; `qarcoal()` within 3.5e-18
  (one value differs by an ulp of `log()`, numpy vs C); the sample plot's
  `computeBalance()` within 6.7e-16.
- `cpus` 1, 2, 4 and sparse vs dense: identical selections and log ratios.
- New: `log_ratio_features()`; `type = 'rank'` on the plots of
  `dirmult_ttest`, `ancombc` / `ancombc2`, `da_ttest`, `da_wilcoxon`,
  `da_lm`, `linda` and ordinations with loadings; class
  `ecodive_log_ratio` with `print()`, `plot()`, `Ops`, `Math`;
  `DA_TESTS$signed` (`ancom` and `lefse` refused).
- Nothing ported from Qurro; not credited (the docs cite the paper).
