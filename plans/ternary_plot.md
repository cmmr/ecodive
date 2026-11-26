# `ternary_plot()`: three-part subcompositions in the simplex

**Priority:** 54. **Tier:** clean room (ggtern GPL-2, Ternary GPL ≥ 2,
robCompositions and compositions GPL ≥ 2: none of their code is read; work
from the papers and their documentation). Read [README.md](README.md)
first, and [compositional.md](compositional.md).


## What it is

A ternary diagram places a three-part composition \eqn{(a, b, c)}, closed to
sum 1, inside an equilateral triangle: each vertex is a pure part, and
gridlines are lines of constant proportion. Two uses dominate:

| Mode | Points | Typical microbiome use |
|:-----|:-------|:-----------------------|
| **samples** | one per sample: a 3-part subcomposition (three features, or three groups of features summed, e.g. Bacteroidota / Bacillota / everything else) | phylum triplets; CoDa exploration |
| **features** | one per feature: its mean relative abundance in each of **three sample groups**, closed across the groups; sized by overall abundance | "where is each ASV enriched" among three sites or treatments (the common microbiome figure) |

Two compositional additions matter in the samples mode:

- **Centring** (von Eynatten, Pawlowsky-Glahn & Egozcue 2002):
  perturbing every sample by the inverse of the sample's centre,
  \eqn{x \ominus g = C(x_1/g_1, x_2/g_2, x_3/g_3)} with \eqn{g} the closed
  geometric mean, moves the centre to the barycentre. Microbiome triplets
  crowd against one vertex; centring spreads them without changing any log
  ratio *difference* between samples. Gridlines are perturbed the same way,
  so they become curves that still read in the original proportions.
- **Confidence ellipses in ilr space**, mapped back: an ellipse around the
  ilr mean becomes a closed curve in the simplex (straight lines in ilr
  space are curves in the triangle).

**Uptake:** "ternary plot" appears in 331 microbiome-filtered Europe PMC
articles 2022–26 (410 under a looser 2019–26 query), flat (1.02×
normalised), far more than "compositional biplot" (15); ggtern has 641
citations. Some hits are soil-texture ternaries inside microbiome papers,
and the samples/features split is unquantified. Good demand, low–medium
effort.


## Names

- `ternary_plot()`. It **draws immediately** and returns the drawn
  coordinates invisibly, like `graphics::barplot()` / `boxplot()`; the
  `_plot` suffix says so. (ecodive's other results are computed first and
  plotted by a `plot()` method; that pattern doesn't fit here, because the
  coordinates are not useful without the picture.) Plain `ternary()` would
  read as a computation, and as the ternary operator.
- Same names elsewhere: `vcd::ternaryplot()` (GPL; different spelling, so
  no masking), `Ternary::TernaryPlot()`, `robCompositions::ternaryDiag()`,
  `compositions::plot.acomp()`. Say so in the Rd.
- Arguments from the glossary: `counts`, `groups`, `data`, `margin`,
  `pseudocount`, `conf_level`, `legend`, `col`, `pch`, `labels` (a count, as
  in the volcano plots). New: `parts` (the three parts; the CoDa term),
  `center` (R's spelling, as `scale(center =)`), `ellipses`, `grid`. Add
  them to [parameter_names.md](parameter_names.md) §1.


## Design

### Signature

```r
ternary_plot(counts, parts = NULL, groups = NULL, data = NULL, margin = 1L,
             pseudocount = NULL, center = FALSE, ellipses = 'none',
             conf_level = 0.95, grid = TRUE, labels = 10, legend = 'topright',
             col = NULL, pch = NULL, ...)
```

- **Mode:** `parts` given → samples mode; `parts = NULL` → features mode,
  which needs `groups` with exactly three levels (error otherwise: "Without
  `parts`, `groups` must have exactly 3 levels (it has 4): each feature is
  placed by its mean abundance in the 3 groups.").
- **`parts`:** a vector selecting exactly 3 features (any `feature_set()`
  form), or a **list of 3** selectors whose features are summed
  (amalgamation; the docs say it is not a log-ratio operation but is what
  "phylum" triplets mean). One list element may be `NA`, meaning every
  feature not in the other two. Names of the list (or the features) label
  the vertices. Overlap → `feature_set()`'s error.
- **`groups` / `data`:** `per_sample()` / `match_rows()` as in
  `plot.ecodive_ord()` (`plot_groups()`); in samples mode they colour the
  points and give one ellipse per group; numeric `groups` → gradient, no
  ellipses.
- **`center`, `ellipses`:** samples mode only (error in features mode,
  where group means are often 0). `ellipses = 'none'`, `'confidence'` (for
  the mean) or `'prediction'` (for a new sample), at `conf_level`.
- **Zeros:** the log is needed only for `center` and `ellipses`; then
  `pseudocount` follows the shared rule (`log_pseudocount()`, as
  `balance()`), applied to the three amalgamated parts, and used for the
  points too so that they match the ellipses. Without logs, zeros are drawn
  on the edges and `pseudocount` is ignored (the subtitle says which).

### Computation (samples mode)

1. \eqn{y_{i\cdot}} = the three part sums (+ pseudocount), closed:
   \eqn{p_{ik} = y_{ik} / \sum_k y_{ik}}. Samples whose three parts are all
   0 are dropped with one warning (counted in the subtitle).
2. **Centre:** \eqn{g_k = \exp(\overline{\log p_{\cdot k}})}, closed;
   \eqn{p'_{ik} = C(p_{ik}/g_k)}.
3. **ilr** with the fixed orthonormal basis
   \eqn{V = \begin{pmatrix} 1/\sqrt2 & -1/\sqrt2 & 0 \\ 1/\sqrt6 & 1/\sqrt6 & -2/\sqrt6\end{pmatrix}}:
   \eqn{z_i = V \log p_i} (rows of V sum to 0, so centring the log is not
   needed). Any orthonormal basis gives the same curve.
4. **Ellipse** for a group of \eqn{m \ge 3} samples with mean \eqn{\bar z},
   covariance \eqn{S} (divisor \eqn{m-1}), \eqn{F = } `qf(conf_level, 2, m - 2)`:
   - confidence: \eqn{r^2 = \frac{2(m-1)}{m(m-2)} F} (Hotelling's region for
     the mean);
   - prediction: \eqn{r^2 = \frac{2(m+1)(m-1)}{m(m-2)} F}.
   Points \eqn{\bar z + r\,L^\top (\cos\theta, \sin\theta)^\top}, \eqn{L}
   = `chol(S)`, 200 values of θ; back-transformed by
   \eqn{p = C(\exp(V^\top z))}, then centred as the points. A group with
   \eqn{m < 3} or singular \eqn{S} gets none, with one warning naming them.
5. **Barycentric → Cartesian:** vertices (0, 0), (1, 0), (1/2, √3/2) for
   parts 1, 2, 3: \eqn{x = p_2 + p_3/2}, \eqn{y = \frac{\sqrt3}{2} p_3}.

### Computation (features mode)

Close every sample over **all** features (relative abundances); for feature
\eqn{j} and group \eqn{g}, \eqn{m_{jg}} = mean over the group's samples;
point \eqn{= C(m_{j1}, m_{j2}, m_{j3})}. Features absent from all three
groups are dropped (counted). Point size: `cex` from 0.4 to 3 in proportion
to \eqn{\sqrt{\bar m_j}} (overall mean relative abundance). Default colour:
the group's palette colour when the feature's largest share is ≥ 0.5
("mostly in that group"), else `PLOT_MUTED`; the docs say this is a visual
convention, not a test, and that `col` (one per feature, by name) can carry
DA results instead. The `labels` most abundant features are named.

### Drawing

- `plot.new()`, `plot.window(c(0, 1), c(0, √3/2), asp = 1)`; the triangle;
  vertex labels outside the corners (`plot_label()`).
- **Grid** (`grid = TRUE`): lines \eqn{p_k = 0.2, 0.4, 0.6, 0.8} for each
  part in `'grey85'`, tick labels on the edge they meet. When centred, each
  line is 101 points of constant *original* proportion, perturbed by
  \eqn{g^{-1}}, closed and drawn with `lines()`: curves whose end points stay
  on the edges, labelled as before. Subtitle: "Centred: perturbed by the
  inverse of the geometric mean (a, b, c)".
- Points (`ord_plot_look()` colours/symbols), ellipses (`lines()` in the
  group colour), legend via `plot_legend()`. Title: "Ternary plot of a, b,
  c" / "Mean relative abundance in A, B, C". All through `plot_do()`, so
  `main`, `sub`, `cex` ... in `...` override.

### Value

Invisibly, a `data.frame` with one row per point: the three closed (and
centred) proportions named after the parts, `x`, `y`, and `group` (samples
mode) or `abundance` (features mode); attributes `center` (\eqn{g}, or
`NULL`), `pseudocount`, `mode`, and `ellipses` (a list of x/y matrices by
group).


## Reference implementations

- **None read.** Formulas: barycentric coordinates (textbook); centring by
  perturbation (von Eynatten et al. 2002; Aitchison 1986); ilr (Egozcue et
  al. 2003); Hotelling regions (any multivariate text).
- **Output-only cross-checks (GPL; run, don't read):**
  `compositions::plot(acomp(x), center = TRUE)` and
  `compositions::ellipses()` for the curve shapes; ggtern's
  `geom_confidence_tern()` only if convenient. Not fixtures: the identities
  below are exact. **Unverified:** whether ggtern's confidence curves use
  Hotelling's F or a χ² radius; ecodive documents its own (F).


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- A plotting function: it draws, then returns its coordinates (see
  Value) invisibly, with no `display` attribute; nothing plots that table
  again.
- Title from the parts and mode, as in Drawing ("Ternary plot of a, b,
  c" / "Mean relative abundance in A, B, C"). Counts have no title, so
  nothing is added from the input, and an earlier zero replacement is not
  named.

## Tests (`test-ternary_plot.r`)

- Barycentric map: the three unit compositions go to the vertices;
  (1/3, 1/3, 1/3) to (1/2, √3/6) exactly.
- Centring: the centred points' closed geometric mean is (1/3, 1/3, 1/3)
  at 1e-14; their ilr mean is 0; ilr differences between samples are
  unchanged.
- Ellipses: every curve point's ilr \eqn{z} satisfies
  \eqn{(z-\bar z)^\top S^{-1}(z-\bar z) = r^2} at 1e-10, for both kinds;
  centred and uncentred curves differ by the same perturbation.
- Amalgamation: a list of selectors sums the right columns; `NA` takes the
  rest; dense, sparse and `margin = 2` agree.
- Features mode: means equal a hand computation on a small table.
- `pdf(NULL)`: returns invisibly (not `x`: the data frame, `withVisible`
  FALSE); `.Random.seed` untouched; numeric `groups`; `grid = FALSE`.
- Errors: `parts` of 2 or 4; a list with two `NA`s; features mode with 2 or
  4 levels; `center` / `ellipses` in features mode; `conf_level` 0 or 1;
  `ellipses = 'x'`; warning for groups with < 3 samples.
- One smoke line.


## Docs

- Rd in the reference group **"Plotting"**, new subtitle **"Compositional
  displays"** (with `composition_barplot`, `coda_dendrogram`, `tree_plot`).
- Explain the two modes, why to centre (and that centred gridlines are
  curves), confidence vs prediction regions, amalgamation, and the colour
  convention of features mode. "Same names in other packages" as above. NEWS.


## References

- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
  Chapman & Hall.
- von Eynatten, H., Pawlowsky-Glahn, V., & Egozcue, J. J. (2002).
  Understanding perturbation on the simplex: a simple method to better
  visualize and interpret compositional data in ternary diagrams.
  *Mathematical Geology*, 34, 249–257. \doi{10.1023/A:1014826205533}
- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal,
  C. (2003). Isometric logratio transformations for compositional data
  analysis. *Mathematical Geology*, 35, 279–300. \doi{10.1023/A:1023818214614}
- Hamilton, N. E., & Ferry, M. (2018). ggtern: ternary diagrams using
  ggplot2. *Journal of Statistical Software*, 87(3). \doi{10.18637/jss.v087.c03}


## Decisions made while implementing

- **Default pseudocount: the whole table** (the user's choice, 2026-11-26).
  With `pseudocount = NULL`, the shared rule (`log_pseudocount()`, as
  `balance()`) looks at the whole of `counts`: none without zeros, else
  half the smallest non-zero value, with its usual warning. The value is
  added once to each of the three part sums. It is used only with
  `center = TRUE` or `ellipses`. A given `pseudocount` is checked (it must
  be positive) even when unused. The `pseudocount` attribute is 0 when no
  logs were taken.
- **No `cpus`** (the user's choice). The sums run in a single-threaded C
  pass over the sparse table (`src/ternary_plot.c`: `C_ternary_sums()`
  for samples mode, `C_ternary_means()` for features mode), in a fixed
  order.
- **F vs χ² (the plan's open question), settled from outputs.**
  ggtern 4.0.0's `geom_confidence_tern()` (run, not read) draws the curve
  \eqn{(z-\bar z)^\top S^{-1}(z-\bar z) =} `qchisq(level, 2)`, with \eqn{S}
  of divisor m − 1, exactly, at m = 6, 20 and 100. That is a normal
  contour of the samples, the large-m limit of the prediction region, and
  not a region for the mean. ecodive keeps Hotelling's F regions, as
  planned, and the Rd says how ggtern differs. `compositions::ellipses()`
  draws the same quadratic form with a radius the caller gives. Its
  vertices are exactly at r², and the points it interpolates between them
  fall up to 0.2% inside. Notably, it draws them **uncentred** even on a
  `plot(acomp, center = TRUE)` plot.
- `compositions::plot(acomp(x), center = TRUE)` (compositions 2.0-9) puts its vertices where
  ecodive does: (0, 0), (1, 0), (1/2, √3/2). Its centred points match
  ecodive's to 2.2e-16.
- **Ellipses** are drawn with `polygon()` through 200 values of θ (no
  repeated end point). The `ellipses` attribute holds 200 × 5 matrices
  (the three proportions, `x`, `y`) rather than x/y only, so the identities
  can be checked from the result. Without `groups`, one region for all the
  samples, named `"all"`. A covariance is treated as singular when its
  smaller eigenvalue is ≤ 1e-10 × its larger. Numeric `groups` with
  `ellipses` warn and draw none.
- **Samples with all three parts 0** are left out on their raw sums
  (before any pseudocount), with one warning and a subtitle count.
- **The result's `group` column** is present only when `groups` is given.
- **Part labels:** the list's names, else a single feature's name,
  `"Other"` for the `NA` part, `"Part k"` for an unnamed set of several
  features. A named vector's names win too. `make.unique()` guards
  duplicates. Parts keep the order given; vertices are parts 1, 2, 3 at
  bottom left, bottom right and top.
- **`labels`** is used in features mode only.
- **Features mode:**
  - Samples whose group is `NA` are in no group. Numeric `groups` become a
    factor. A group with no samples (an unused factor level, or all-zero
    samples) is an error, and all-zero samples are left out with a warning.
  - `abundance` (which sizes the points) is the mean relative abundance
    over all samples of the three groups, not the mean of the three group
    means.
  - Colours come from `plot_palette(3)`; ties for the largest share go to
    the first group. The legend is titled "Mostly in" and adds "No
    majority".
  - `col` is one colour, one per feature of `counts`, or named by feature
    (partial). Giving `col` drops the legend, since the convention no
    longer holds.
  - Points are drawn largest first, so that small ones stay visible.
- **Gridlines:** tick labels run counter-clockwise: part 2 along the
  bottom, part 3 up the right edge, part 1 down the left. When centred,
  gridlines crowd together, so a tick label that would overlap another, or
  a vertex label, is skipped; its line and tick are still drawn.
- **Subtitle:** plot.default's `sub` clipped long subtitles, so the pieces
  ("Centred: ... (g)", "pseudocount c", "zeros drawn on the edges",
  dropped counts, "95% confidence regions") are joined with "; " and
  wrapped into the bottom margin with `mtext()`. A user's `sub` replaces
  them.
- **`...`:** `cex`, `lwd` and `bg` go to `points()`; everything else
  (`main`, `xlim`, ...) goes to `plot.default()`, which sets up the frame.
- `plot_groups()` (R/plot_ord.r) gained a `ref` argument, so its messages
  name `counts` rather than `x`.
- `?ternary_plot` inherits the "Input Types" section.
  `composition_barplot()`'s docs now link to it.
