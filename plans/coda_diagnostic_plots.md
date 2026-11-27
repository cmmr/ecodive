# `variation_array()` and `zero_pattern()`: two compositional diagnostics

**Priority:** 60. **Tier:** clean room (easyCODA, propr, zCompositions,
robCompositions: GPL, no code read). The variation array is formulas over
ecodive's own `vlr()`; the zero patterns are bookkeeping. Read
[README.md](README.md) first, and [compositional.md](compositional.md).


## What it is

| Function | Shows | Classic source |
|:---------|:------|:---------------|
| `variation_array()` | for every pair of features, the variance of their log ratio \eqn{\tau_{ij} = \mathrm{var}(\log x_i/x_j)} (above the diagonal) and its mean \eqn{\xi_{ij} = \overline{\log(x_i/x_j)}} (below), plus each clr component's share of the total variance | Aitchison (1986); shown with "% clr var" and "Tot var" in Pawlowsky-Glahn & Egozcue (2011, Fig. 2) |
| `zero_pattern()` | which combinations of zero features occur, how often, and how often each feature is zero | zCompositions' `zPatterns()` (Palarea-Albaladejo & Martín-Fernández 2015); robCompositions' `zeroPatterns()` |

The variation array is the standard first summary of a composition: small
\eqn{\tau_{ij}} means \eqn{i} and \eqn{j} are nearly proportional; the clr
variances decompose the total variance \eqn{\mathrm{totvar} = \frac{1}{2D}
\sum_{i,j}\tau_{ij}}, which the CoDa-dendrogram
([coda_dendrogram.md](coda_dendrogram.md)) splits by balance. ecodive's
`vlr()` already returns \eqn{\tau} (as a `dist`, plotted as a heatmap);
the array adds the means and the variance decomposition in one display.

The zero-pattern display tells users, before choosing a pseudocount or a
`replace_zeros()` method, whether zeros are scattered (count zeros) or
structured (whole groups of features absent together; see
`structural_zeros()`).

**Uptake:** modest. zCompositions has 1,105 citations, mostly for its
imputation, not the pattern plot; the variation array is cited through
Aitchison (1986) and propr / Lovell et al. (2015, 332). Low effort,
diagnostic value; nobody counted papers that show either display.


## Names

- `variation_array()`: Aitchison's own term; spelled out (README rule 4)
  because `variation()` (compositions' name for the \eqn{\tau} matrix alone)
  is too generic. `zero_pattern()`: as listed in [TODO.md](TODO.md); it
  returns each sample's pattern. No collisions in `NAMESPACE`.
- **Compute, then plot:** both return an object with a `plot()` method, as
  `vlr()` and the tests do, because the numbers are useful on their own (the
  array is a table people print; the patterns label samples).
- Arguments: `counts`, `margin`, `pseudocount`, `clust`, `cpus` (as
  `vlr()` / `beta_div()`); plot methods: `digits`, `key` (as
  `plot.ecodive_dist()`), `n_features`, `n_patterns` (counts, `n_<unit>`).


## Design

### `variation_array()`

```r
variation_array(counts, margin = 1L, pseudocount = NULL, clust = NULL,
                cpus = n_cpus())
```

- \eqn{\tau}: **call `vlr()`'s C worker** (same pseudocount rule, divisor
  \eqn{n-1}, sparse input stays sparse). \eqn{\xi_{ij} = \bar\ell_i -
  \bar\ell_j}, \eqn{\ell = \log(x + c)}, from one pass of column means of
  logs (computed in C beside the variances, or in R from `colMeans()` of the
  log table for dense input).
- **Result:** a dense \eqn{D \times D} numeric matrix, class
  `c('ecodive_variation_array', 'matrix', 'array')`, dimnames = features.
  Element `[i, j]`, \eqn{i < j}: \eqn{\tau_{ij}}; \eqn{i > j}:
  \eqn{\xi_{ij} = } mean of \eqn{\log(x_i / x_j)}, **row over column**;
  diagonal `NA`. Attributes: `clr_var` (named), `total_var`,
  `pseudocount`, `n` (samples), `method = 'variation_array'`, and
  `display = list(title = 'Variation Array', clust = ...)` so titles
  work.
  - \eqn{\mathrm{clr\_var}_i = (\sum_j \tau_{ij} - \mathrm{totvar})/D},
    \eqn{\mathrm{totvar} = \frac{1}{2D}\sum_{i,j}\tau_{ij} = \sum_i
    \mathrm{clr\_var}_i} (since \eqn{\sum_j \tau_{ij} = D\,v_i + \sum_k
    v_k}, the clr components summing to 0). Exact identities: tested.
- **`clust`:** as in `beta_div()` (`validate_clust()`, `DIST_CLUST`):
  reorders the features by the clustering of \eqn{\sqrt{\tau}} (a Euclidean
  distance between the features' centred log vectors, so `ward.D2` is
  valid) or by `'pcoa'`, recorded as `order` (top level) and
  `display$clust`.
- **Size:** dense \eqn{D^2}; above 2,000 features warn once, as large-dense
  results do elsewhere (`log_ratio_input(dense_out =)`), and point to
  filtering.
- `print.ecodive_variation_array()`: the matrix rounded to 3 significant
  digits (class and attributes stripped), then "Total variance: …".

**Sign convention unverified:** Pawlowsky-Glahn & Egozcue (2011, Fig. 2)
put means in the lower triangle but the text doesn't say whether
\eqn{\log(x_i/x_j)} or \eqn{\log(x_j/x_i)}; Aitchison (1986) not checked
(book). ecodive states its own convention (row over column) in the Rd.

`plot.ecodive_variation_array(x, digits = 2, key = TRUE, ...)`:

- Feature order: `x`'s own if `clust` was set, else average-linkage
  clustering of \eqn{\tau} (as `plot.ecodive_dist()` does by default, via
  `dist_clust_order()`).
- One `image()` call per triangle on a shared grid (first feature at top):
  upper cells on `hcl.colors(64, 'YlOrRd', rev = TRUE)` over
  \eqn{[0, \max\tau]}; lower cells on `hcl.colors(65, 'Blue-Red 3')`
  symmetric about 0; diagonal blank. Cell values written (`digits`) when
  \eqn{D \le 15}.
- Right margin: a column of bars, `clr_var / total_var` in % ("% clr var",
  as the 2011 figure), then two keys (`dist_plot_key()`), upper half for
  \eqn{\tau}, lower half for \eqn{\xi}. Title `display$title` ("Variation
  Array"); subtitle
  "D features, n samples, total variance 2.55; upper: var(log ratio),
  lower: mean log(row / column)".
- Margins restored with `on.exit(par(opar))`. Returns `x` invisibly.

### `zero_pattern()`

```r
zero_pattern(counts, margin = 1L)
```

- Each sample's **pattern** is the set of its zero features. Key per sample:
  the indices of its **non-zero** features pasted together (sparse input:
  straight from the column indices of a `dgCMatrix`, never densified; dense:
  `which(x[i, ] != 0)`). Unique keys → patterns.
- Pattern ids: 1 = the most frequent; ties by first occurrence in sample
  order. (zCompositions' numbering is not documented; no parity claimed.)
- **Result:** an integer vector, one id per sample, named by sample, class
  `ecodive_zero_pattern`. Attributes: `patterns` (logical matrix, patterns
  × features, `TRUE` = zero, rownames = ids), `n_samples` (per pattern),
  `zeros` (proportion of samples zero, per feature), `n`, `D`. A
  `print()` method: "12 zero patterns among 40 samples; 3 features never
  zero; most frequent pattern: 9 samples, zeros in 'f3', 'f8'".

`plot.ecodive_zero_pattern(x, n_patterns = 30, n_features = 50, ...)`:

- `graphics::layout()` of a main panel, a bar panel above it and one to its
  right (restored on exit).
- **Features shown:** those zero in some but not all samples (always-present
  and always-absent ones are counted in the subtitle, not drawn). If more
  than `n_features`, the `n_features` whose zero proportion is closest to
  0.5 (the most informative). Columns sorted by zero proportion, fewest
  zeros left.
- **Patterns shown:** the `n_patterns` most frequent, top to bottom.
  Patterns are defined on **all** features, so two rows may look identical
  when they differ only in hidden features; the subtitle says how many
  features are hidden.
- Main panel: `image()` with zeros dark (`'grey20'`), non-zeros
  `'grey95'`; top bars: % zeros per shown feature; right bars: % of samples
  with each shown pattern (`barplot(horiz = TRUE)`). Title "Zero patterns";
  subtitle "k of m patterns shown (p% of samples); f features hidden".
- Returns `x` invisibly.
- Microbiome tables often give one pattern per sample; the docs say so and
  that the feature bars are then the useful part.


## Reference implementations

- **None ported.** \eqn{\tau} is `vlr()` (scikit-bio's formula, already
  ported); the rest are identities above.
- **Output-only (GPL):** `compositions::variation()` for \eqn{\tau};
  `zCompositions::zPatterns()` for the number of patterns and their sample
  counts on a small table (ids may differ). Fixtures in
  `data-raw/expected/coda_diagnostics.r` if run; not required.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- `variation_array()`: sets `method = 'variation_array'` and `display =
  list(title = 'Variation Array', clust = ...)` with `set_display()`; its
  plot reads `display$title`, `order`, `clr_var`, `total_var` and `n`.
- `zero_pattern()`: no `display`; its input is a count table, which has no
  title, so the plot titles itself "Zero patterns". It reads `patterns`,
  `n_samples`, `zeros`, `n` and `D`.

## Tests (`test-coda_diagnostics.r`)

- `variation_array()`: upper triangle equals `as.matrix(vlr(x))` exactly
  (same worker); lower equals `outer(m, m, '-')` of column mean logs at
  1e-14; `clr_var` equals `apply(clr, 2, var)` at 1e-12 and sums to
  `total_var`; dense, sparse, `margin = 2` identical; `cpus` 1, 2, 4
  identical; `clust` reorders and sets `order` and `display$clust`;
  `display$title` is "Variation Array"; pseudocount rule and warning as
  `vlr()`.
- `zero_pattern()`: a hand-built 6 × 4 table with known patterns, ids and
  counts; sparse = dense; no zeros → one pattern; `margin = 2`.
- `pdf(NULL)`: both plot methods return invisibly; `n_patterns = 1`,
  `n_features = 1`; D ≤ 15 writes values; `.Random.seed` untouched; `par()`
  restored (`identical(par('mfrow'), c(1L, 1L))` after the zero plot).
- Errors: `n_patterns = 0`, `n_features = -1`, `clust = 'x'`, one feature.
- Smoke lines for both.


## Docs

- Both Rd pages in **"Compositional Data"**, new subtitle **"Diagnostics"**;
  the plot methods in "Plotting".
- `variation_array`: the formulas, the sign convention, the identities with
  `vlr()` and the clr variances, link to `coda_dendrogram()` and the
  covariance biplot (both decompose the same total variance).
- `zero_pattern`: how to read it before `replace_zeros()`, and
  `structural_zeros()`. NEWS.


## References

- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
  Chapman & Hall.
- Pawlowsky-Glahn, V., & Egozcue, J. J. (2011). Exploring compositional
  data with the CoDa-dendrogram. *Austrian Journal of Statistics*, 40,
  103–113. \doi{10.17713/ajs.v40i1&2.202}
- Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., & Bähler,
  J. (2015). Proportionality: a valid alternative to correlation for
  relative data. *PLoS Computational Biology*, 11, e1004075.
  \doi{10.1371/journal.pcbi.1004075}
- Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015). zCompositions —
  R package for multivariate imputation of left-censored data under a
  compositional approach. *Chemometrics and Intelligent Laboratory Systems*,
  143, 85–96. \doi{10.1016/j.chemolab.2015.02.019}


## Decisions made while implementing

- **Clustering distance (the user's choice, 2026-11-27):** \eqn{\sqrt{\tau}}
  for both `variation_array(clust =)` and the plot's default order (average
  linkage), so that `clust = 'average'` and the default plot agree. The plan
  had \eqn{\tau} itself for the plot default; that order can differ from
  `plot(vlr(x))`'s.
- **Printing (the user's choice, 2026-11-27):** `print()` shows a header
  (features, samples, total variance), the top-left block of at most
  10 × 10 features at 3 significant digits, a footer naming how many
  features are not shown, then the shown features' % clr variance. Not the
  whole matrix (README: results never print every value).
- **The means from C:** `C_proportionality()` takes a fourth method code
  (vlr plus attribute `mean_log`, each feature's mean log from the same
  pass that centres the logs), so \eqn{\tau} is bit-identical to `vlr()`
  and sparse input stays sparse. No new entry point for the variation
  array; `init.c` gains only `C_zero_pattern`.
- **`zero_pattern()` in C** (`src/coda_diagnostics.c`): samples are sorted
  by their lists of non-zero features (length, then features, then sample:
  a total order, so deterministic), and runs of equal lists are the
  patterns. Single-threaded (`cpus` isn't an argument: the plan's
  signature). Explicitly stored zeros in sparse input count as zeros (the
  ecomatrix reader drops them).
- **Labels:** a variation array always has dimnames: feature positions
  when `counts` has no feature names (as `dist_reorder()` labels unnamed
  samples).
- **`pseudocount` attribute:** `0` when none was added.
- **Arithmetic** (`Ops`, `Math`) on a variation array returns a plain
  matrix, as for `ecodive_log_ratio`; `[` already does.
- **Zero-pattern plot:** title and subtitle in the outer margins (the
  layout has three panels), wrapped to the device's width; `main` and
  `sub` in `...` replace them. The top bars' axis runs to the largest
  shown proportion of zeros (rounded up), not always 100%. White grid
  lines between cells up to 80 × 80. The subtitle breaks the hidden
  features down into never zero, always zero and beyond `n_features`.
  With no feature that varies, the main panel says so.
- **Variation-array plot:** `dist_plot_key()` (R/plot_dist.r) gained
  optional `x` and `y` placements so the two keys share the right margin;
  its default placement is unchanged. Bars' values are written up to 40
  features; cell values up to 15 (as planned).
- **Fixtures:** `data-raw/expected/coda_diagnostics.r` stores
  compositions 2.0.9's `variation()` and `mvar()` (two positive tables)
  and zCompositions 1.6.2's `zPatterns()` ids (three tables with zeros).


## Outcome (from the implementing session's final report)

- **Commit:** `212a876` (code, docs, tests, fixtures).
- **Files:** `R/coda_diagnostics.r`, `src/coda_diagnostics.c`
  (`C_zero_pattern`), a variation mode in `src/proportionality.c`,
  `test-coda_diagnostics.r` (141 tests), 4 smoke lines,
  `data-raw/expected/coda_diagnostics.r` and its `.rds`.
- **References:** compositions 2.0.9 `variation()`: largest relative
  difference 8e-16 (tested at 1e-12); `mvar()` against `total_var`, 7e-15
  absolute. Upper triangle `identical()` to `vlr()`; lower triangle within
  1.3e-15 of the column means of logs; `clr_var` within 6e-16 (relative)
  of `apply(clr, 2, var)`. zCompositions 1.6.2 `zPatterns()`: the same
  number of patterns, sample counts and partition on all 4 tables tried
  (6 × 4 hand-built, 40 × 8, `ex_counts`, 100 × 15 with 89 patterns).
- **Determinism:** dense, `dgCMatrix`, slam and `margin = 2` results
  `identical()`; `cpus` 1, 2 and 4 `identical()`.
- **Tests:** 34,151 at home (baseline 34,006) and 481 in CRAN mode
  (baseline 477), all passing, against an `-O2` install, R 4.3.3. `R CMD
  check --no-manual --ignore-vignettes`: only the existing installed-size
  NOTE.
- **Plots** checked by eye as PNG: 6, 25 and 120 features (variation
  array); 4 to 80 patterns, `n_patterns = 1, n_features = 1`, and a small
  device where the subtitle wraps (zero patterns).
