# Biplot scalings and links: `plot(ord, scaling =, links =)`

**Priority:** 44. **Tier:** formula (no reference code is read or ported;
compositions, robCompositions, easyCODA and coda.base are GPL and
output-only). Read [README.md](README.md) first, and
[compositional.md](compositional.md) for the conventions the compositional
plans share.


## What it is

`plot.ecodive_ord()` draws a **form biplot** today: samples at their
scores, loadings as arrows rescaled to fit, so only arrow directions mean
anything. The compositional biplot that Gloor et al. (2017) recommend to
microbiome users is Aitchison & Greenacre's (2002) **relative variation
biplot**: a PCA of clr data drawn in one of two scalings.

| Scaling | Rows (samples) | Columns (features) | What reads off the plot |
|:--------|:---------------|:-------------------|:------------------------|
| `'form'` (default, today's plot) | principal: \eqn{U\Gamma} (the scores) | standard: \eqn{V}, rescaled to fit | distances between samples approximate Aitchison distances |
| `'covariance'` | standard: \eqn{U\sqrt{n-1}} (unit variance per axis) | principal: \eqn{V\Gamma/\sqrt{n-1} = V\sqrt{\lambda}} | ray lengths ≈ clr standard deviations; **links** between ray tips ≈ SD of \eqn{\log(x_j/x_k)}; cosines of angles between links ≈ correlations of log ratios |
| `'contribution'` | principal (the scores) | contribution coordinates \eqn{V} (equal feature weights), one stated multiplier | features far from the origin contribute most to the plane (Greenacre 2013) |

**Interpretation rules, verified** against the primary source: the
open-access working-paper version of Aitchison & Greenacre (UPF Economics
Working Paper 557, 2001, econ-papers.upf.edu/papers/557.pdf), section 4,
Properties 1–8. The research note had them from a secondary source
(MetricGate docs); they agree, and the paper adds:

- **Property 3:** in the covariance biplot, \eqn{|JJ'|} (the link between
  two column points) approximates the SD of that log ratio, **from below**:
  in the full space the link lengths are exactly the SDs, and projection
  shortens them.
- **Property 4:** cosines of angles between links estimate correlations
  between log ratios (perpendicular links ⇒ uncorrelated log ratios).
- **Property 5:** three column points on a line, at distances λ (A–B) and
  μ (B–C), diagnose the constant log contrast
  \eqn{\mu\log A + \lambda\log C - (\lambda+\mu)\log B}.
- **Property 6:** four points forming a parallelogram A, B, C, D diagnose
  \eqn{\log A - \log B + \log C - \log D} = constant.
- Section 2.2 / 3: because clr data are **row-centred**, links are optimal
  least-squares representations of the pairwise log ratios, "explaining the
  same percentage of variance as the rays do". This holds only for clr
  input with `scale = FALSE`.
- The paper rescales the column points by \eqn{1/\sqrt{n-1}} so they are in
  log-ratio SD units and draws rows as \eqn{U}. ecodive draws rows as
  \eqn{U\sqrt{n-1}} (unit variance), a uniform rescale of the rows that
  changes no reading, so that both sets fill the plot.

The Pawlowsky-Glahn & Egozcue (2011, open access) summary agrees: "Length of
links between the end point of rays are approximately proportional to the
corresponding variance of the pairwise log-ratio" (they say variance; the
2002 paper says SD for the covariance scaling; the SD is what the algebra
below gives). **Not verified:** Greenacre's (2013) contribution-biplot text
(paywalled; abstract only) and easyCODA's code (GPL, not read). easyCODA's
`PLOT.LRA` docs say only that `map = "contribution"` puts "columns ... in
contribution coordinates" and that `rescale` multiplies the column
coordinates. Greenacre's LRA weights parts by their masses; ecodive's PCA is
unweighted, so with equal weights \eqn{1/p}, standard coordinates are
\eqn{\sqrt{p}V} and contribution coordinates (standard × √mass) are exactly
\eqn{V}. So `'contribution'` and `'form'` have the same column geometry up to
a constant: what differs is the stated, fixed multiplier and the marking of
above-average contributors. Say so in the docs.

**Uptake:** Aitchison & Greenacre 2002, 830 OpenAlex citations; contribution
biplots 91. 53 Europe PMC microbiome papers 2020–26 mention "clr / compositional
/ Aitchison biplot", but only 15 use the exact phrase "compositional biplot"
in 2022–26, and that is declining (0.23× normalised). Many clr-PCA figures
are never named, so these undercount; nobody counted papers that read links.
Cheap because ecodive already stores everything needed.


## Names

- `scaling`: the glossary's reserved name for score scaling
  ([parameter_names.md](parameter_names.md) §1). Values `'form'`,
  `'covariance'`, `'contribution'` (Aitchison & Greenacre's and Greenacre's
  terms). Also accept vegan's numbers as aliases: `1` = `'form'` (vegan's
  distance biplot), `2` = `'covariance'` (its correlation biplot). Update the
  glossary row from "reserved" to "used by `plot.ecodive_ord()`".
- `links`: the paper's term. No collision with any existing argument.
- No new exports: only new arguments to an existing S3 method.


## Design

### Signature (after this item; [log_ratio_plots.md](log_ratio_plots.md) later adds `type`)

```r
plot.ecodive_ord(x, groups = NULL, data = NULL, axes = c(1, 2), loadings = 10,
                 biplot = TRUE, centroids = TRUE, target = TRUE,
                 legend = 'topright', col = NULL, pch = NULL,
                 scaling = 'form', links = NULL, ...)
```

New arguments go after `pch` so every existing positional call
(`plot(ord, 'Species', iris)`) still works.

### What is stored, and what is missing

For \eqn{n} samples (rows of `x`), axes \eqn{a, b} = `axes`:

| Method | scores | `loadings` | `eig` | notes |
|:-------|:-------|:-----------|:------|:------|
| `pca()` | \eqn{US} (prcomp's `x`) | \eqn{V} (unit) | \eqn{s^2/(n-1)} | `attr(x, 'norm')`, `attr(x, 'scale')`: added by this plan (see "Plotting") |
| `rpca()`, `rpca_phylo()` | \eqn{u\,\mathrm{diag}(s)} | \eqn{v} (unit) | \eqn{s^2/(n-1)} | `prop` is relative to the rank-`n_axes` completed matrix |

So \eqn{\lambda_k} = `eig[k]`, \eqn{\sqrt{\lambda_k} = s_k/\sqrt{n-1}},
and **no matrix needs adding**: every scaling is a column rescale of stored
matrices. Let \eqn{S} = `x[, axes]`, \eqn{V} = `loadings[, axes]`,
\eqn{\ell = \sqrt{\lambda_{axes}}}.

- `'form'`: rows \eqn{S}; columns \eqn{V}, times today's
  `plot_arrow_mul()` factor. Unchanged.
- `'covariance'`: rows `sweep(S, 2, ell, '/')`; columns
  `sweep(V, 2, ell, '*')`. Drawn **at their true coordinates** (no fitting
  multiplier: link lengths must be readable in log-ratio SD units), as rays
  from the origin with labels at the tips. `xlim`/`ylim` cover both sets.
- `'contribution'`: rows \eqn{S}; columns \eqn{V \cdot c}, with \eqn{c} the
  `plot_arrow_mul()` factor **rounded down to 1, 2 or 5 × 10^k**, written in
  the subtitle ("feature coordinates × 20"). Features whose squared
  coordinate on either plotted axis exceeds \eqn{1/p} (the average
  contribution; \eqn{p} = `nrow(loadings)`) are drawn in `PLOT_HIGHLIGHT`;
  the others in `PLOT_MUTED`, unlabelled.

In every scaling, \eqn{\mathrm{rows} \times \mathrm{cols}^\top} on the
plotted axes equals \eqn{S V^\top} (for `'form'` and `'contribution'`,
before the multiplier): the test checks this.

**Which features are drawn** (`loadings`, a count): by length in the drawn
coordinates (today's `ord_plot_loadings()` rule) for `'form'` and
`'covariance'` (so covariance shows the most variable features); for
`'contribution'`, by contribution to the plane,
\eqn{(\lambda_a v_{ja}^2 + \lambda_b v_{jb}^2)/(\lambda_a + \lambda_b)}.
Features named in `links` are always drawn.

### Links

- `links`: `NULL`, a vector of at least two features (names or indices into
  `rownames(loadings)`): every pair among them is linked; or a two-column
  matrix of feature pairs. Unknown names → error listing up to 5
  (`name_list()`).
- Only with `scaling = 'covariance'`; otherwise error: "`links` needs
  `scaling = 'covariance'`: only there do link lengths approximate the
  standard deviations of log ratios."
- Drawn as dashed `segments()` in `PLOT_HIGHLIGHT`, with the planar length
  (`plot_num(len, 2)`) written at the midpoint with `plot_label()`. The
  docs say the full-space value is larger ("from below").
- **When links mean log ratios:** `pca()` with `attr(x, 'norm') == 'clr'`
  and `attr(x, 'scale')` FALSE; `rpca()` and `rpca_phylo()` (their
  completed rclr matrix is row-centred; for `rpca_phylo()` the features are
  branch-length-weighted node counts, which the docs say). Any other `pca()`
  input: one warning, "links of `pca(norm = 'hellinger')` show standard
  deviations of differences of hellinger values, not of log ratios." Methods
  other than these three: `scaling` other than `'form'` is an error naming
  the three.
- **Low-variance plane:** with `scaling = 'covariance'` the subtitle ends
  "plane explains 62% of the variance". When `links` are drawn and
  `sum(prop[axes]) < 0.5`, one warning: "The two axes explain 38% of the
  variance, so link lengths may greatly underestimate the log-ratio standard
  deviations." For `rpca()` the docs note that `prop` is relative to the
  completed rank-`n_axes` matrix and so overstates the fit to the data.

### Code layout

A helper `ord_plot_scaling(x, axes, scaling, n_load)` in `R/plot_ord.r`
returns `list(rows, cols, mult, sub)`; `plot.ecodive_ord()` uses `rows` in
place of `pts` and draws `cols` by scaling. `ord_plot_links(cols, links)`
validates and draws. One-axis results ignore `scaling` (strip chart).


## Reference implementations

- **None ported.** The formulas are the SVD algebra above and in Aitchison &
  Greenacre (2002, §2–4) and Greenacre (2013).
- **Output-only cross-check (optional, GPL):** easyCODA `PCA()` +
  `PLOT.PCA(map = 'asymmetric')`, or `compositions::biplot` of a clr
  `princomp`. Not needed as fixtures: the identities below are exact.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Nothing new on the result for the scaling: it is a display choice, an
  argument of `plot()`. The title is the ordination's `display$title`,
  set by `ord_result()` (e.g. "Principal Component Analysis"); the
  subtitle names the scaling ("covariance biplot").
- `plot()` decides whether links are log-ratio SDs from top-level
  attributes: `method` in the rpca family, or `method == "pca"` with
  `norm == "clr"` and `scale == FALSE`; anything else gets a warning
  instead of link annotations.
- So `pca()`'s result gains top-level `norm` and `scale` attributes (the
  transform that qualifies its loadings), passed through `ord_result()`.
- Tests: the link check follows `method`, `norm` and `scale`, not the
  result's column names.

## Tests (`tests/tinytests/test-plot.r`, section "Biplot scalings")

Numeric helpers, exactly (at 1e-12):

- `pca(x, norm = 'clr', n_axes = D - 1)` on a zero-free table: the
  distance between covariance column points over **all** axes equals
  `sqrt(as.matrix(vlr(x)))`; their lengths equal `apply(clr, 2, sd)`.
  On two axes every link is ≤ the full value.
- For each scaling: `rows %*% t(cols) / mult` equals `S %*% t(V)`;
  covariance rows have `apply(rows, 2, var) == 1`.
- Contribution: `colSums(loadings^2) == 1`; the drawn features match the
  planar-contribution order computed by hand; `mult` is in {1, 2, 5} × 10^k.
- Same three on `rpca()` and `rpca_phylo(ex_counts, ex_tree)`.

Drawing (into `pdf(NULL)`): `invisibly(ord, scaling = s)` for each scaling
and the aliases `1`, `2`; with `links = c('a', 'b', 'c')` and a pair matrix;
existing calls unchanged; `.Random.seed` untouched.

Errors and warnings: `scaling = 'x'`, `scaling = 3`; a non-`'form'`
scaling on `pcoa()`; `links` without covariance; one link feature; an
unknown feature; the `norm = 'hellinger'` warning; the low-plane warning
(construct a table where PC1+PC2 < 50%).


## Docs

- `plot.ecodive_ord()` Rd: `@param scaling`, `@param links`, and a new
  section **"Compositional biplots"**: the three scalings; Properties 3–6 in
  plain words; when links are log ratios; "from below"; the low-variance
  caveat; that `'contribution'` differs from `'form'` only in its fixed
  multiplier and marking for an unweighted PCA. Example:
  `plot(pca(counts, norm = 'clr'), scaling = 'covariance', links = c(...))`.
- `pca()` / `rpca()` Rd "Biplots" paragraphs point to it.
- `parameter_names.md` §1: the `scaling` row; NEWS.


## References

- Aitchison, J., & Greenacre, M. (2002). Biplots of compositional data.
  *JRSS C*, 51, 375–392. \doi{10.1111/1467-9876.00275}. Working paper:
  UPF Economics Working Paper 557 (2001).
- Greenacre, M. (2013). Contribution biplots. *Journal of Computational and
  Graphical Statistics*, 22, 107–122. \doi{10.1080/10618600.2012.702494}
- Pawlowsky-Glahn, V., & Egozcue, J. J. (2011). Exploring compositional
  data with the CoDa-dendrogram. *Austrian Journal of Statistics*, 40,
  103–113. \doi{10.17713/ajs.v40i1&2.202}
- Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
  (2017). Microbiome datasets are compositional: and this is not optional.
  *Frontiers in Microbiology*, 8, 2224. \doi{10.3389/fmicb.2017.02224}


## Decisions made while implementing

Item 44, 2026-11-24. No reference code was read or ported, and no GPL
source text reached the session (vegan, easyCODA and compositions were only
run, for fixtures).

Asked of the user at the start, and answered:

1. **`scaling` also covers `rda()`, `cca()` and `ca()`,** with `vegan`'s
   scalings 1, 2 and 3 (the user chose this over the plan's PCA-only
   scope). They store scaling 2, so `scaling` defaults to **`NULL`**, "as
   stored": `'form'` for `pca()`, `rpca()` and `rpca_phylo()`, and scaling 2
   for `rda()`, `cca()` and `ca()`. Every existing plot is unchanged.

Made while implementing:

2. **Values.** `'form'` = `1` = `'sites'`; `'covariance'` = `2` =
   `'species'`; `'symmetric'` = `3`; `'contribution'`. Names are matched
   case-insensitively and by unique prefix (`'cov'`; `'co'` and `'s'` are
   ambiguous). `vegan`'s names work for both families, and `'form'` /
   `'covariance'` for `rda()` / `cca()` / `ca()`. `'symmetric'` on the PCA
   family and `'contribution'` on `rda()` / `cca()` / `ca()` are errors
   naming what each takes (the plan's `scaling = 3` error is the first).
   `vegan`'s negative scalings and `hill` / `correlation` options are not
   offered.
3. **Other ordinations:** any non-`NULL` `scaling` is an error naming the six
   methods that take one (the plan allowed `'form'`; with a `NULL` default,
   accepting it on `pcoa()` would only suggest it did something).
   `dbrda()` and `capscale()` were not extended, though the same per-axis
   rule would apply to their sites, constraints and centroids.
4. **`vegan`'s conversions,** from the stored scaling 2, with
   \eqn{f = \sqrt{\lambda/\mathrm{trace}}} for `rda()` (`trace` is
   `tot.chi`, conditional inertia included) and \eqn{\sqrt{\lambda}} for
   `cca()` and `ca()`: scaling 1 multiplies the site scores by \eqn{f} and
   divides the species scores by it, scaling 3 uses \eqn{\sqrt{f}}.
   **Constraint arrows and centroids scale as the sites**: established from
   `vegan` 2.7.8's outputs (and 2.6-4's), not from its code. Constraint
   arrows are still rescaled to fit, so only their per-axis proportions
   change.
5. **Drawing.** `rda()`'s species stay arrows scaled to fit, and `cca()` /
   `ca()`'s species points at their coordinates, in every scaling. The
   PCA family's covariance and contribution biplots draw rays from the
   origin at their own coordinates; the plot's limits cover the tips (times
   1.1, for the labels). Links are dashed `PLOT_HIGHLIGHT` segments labelled
   with `plot_num(length, 2)` at their midpoints.
6. **Contribution multiplier:** `plot_arrow_mul()` of the drawn features,
   rounded down to 1, 2 or 5 × 10^k (`ord_round_mult()`). Highlighting and
   labels mark features with \eqn{v_{ja}^2 > 1/p} or \eqn{v_{jb}^2 > 1/p};
   the others are `PLOT_MUTED` and unlabelled.
7. **Subtitles** only when `scaling` is given: "form biplot", "covariance
   biplot, plane explains 62% of the variance" (the share omitted if `prop`
   is `NA`), "contribution biplot, feature coordinates × 20" (`x` where the
   locale is not UTF-8, from `print_sym()`), and "scaling 1" for `rda()` /
   `cca()` / `ca()`. They follow the existing subtitle (R²) after "; ";
   `sub =` still replaces the whole.
8. **`links`:** a vector (names, positions or a factor) links every pair, in
   `utils::combn()` order; a two-column matrix or data frame gives pairs.
   Duplicated features in a vector are dropped; a pair of a feature with
   itself, `NA`s, logicals and other shapes are errors. Linked features are
   drawn even with `loadings = 0`. On a result outside the PCA family,
   `links` is an error ("drawn only for the results of pca(), rpca() and
   rpca_phylo()") before the `scaling` check.
9. **Warnings,** once per call: `pca()` not of unscaled clr (wording per
   `norm`: "hellinger values", "the untransformed values", "relative
   abundances", "rclr values, whose zeros stay 0"; with `scale = TRUE`,
   "standardized clr values"); and the plane under 50%. A `pca()` result
   with no `norm` attribute gets a generic warning.
10. **`pca()` gains top-level `norm` (the transform's lower-case name, from
    `validate_norm()`'s code) and `scale`,** after `loadings`. Its attribute
    test was updated, and the tests that compared `pca(x, norm = )` with
    `pca(normalize(x))` (three in `test-normalize.r`, one in `test-pca.r`)
    now compare them without `norm`.
11. **Code.** `ord_plot_scaling(x, axes, scaling, n_load, links)` returns
    more than the plan's four fields: `rows`, `cols` (drawn coordinates,
    times `mult`), `mult`, `style` (`arrows` / `rays` / `points`), `col`,
    `label`, `row_mul` and `col_mul` (per-axis multipliers, which
    `plot()` applies to centroids and constraints), `sub`, `scaling` and
    `links` (a data frame of `from`, `to` and planar `length`). Helpers:
    `ord_scaling_arg()`, `ord_round_mult()`, `ord_links_norm_msg()`,
    `ord_link_pairs()`, `ord_plot_rays()`, `ord_plot_links()`, with
    `ORD_PCA_FAMILY` and `ORD_VEGAN_FAMILY`. `ord_plot_loadings()` was
    removed: its rule (longest first) is now inside `ord_plot_scaling()`.
    `plot()` still returns `x` invisibly; the tests check the helper's
    coordinates.
12. **One-axis results** (strip charts) validate the `scaling` value but
    otherwise ignore `scaling` and `links`.
13. **Covariance rows are \eqn{S/\sqrt{\lambda}}** (= \eqn{U\sqrt{n-1}}).
    `rpca()`'s and `rpca_phylo()`'s scores are column-centred, so their
    rows have unit variance too (tested).
14. **Fixtures** (`data-raw/expected/biplot_scaling.r` →
    `tests/tinytests/expected/biplot_scaling.rds`, 95 kB), all output-only:
    `vegan` 2.7.8 `scores(scaling = 1, 2, 3)` of an RDA of clr values and of
    `rda()` (partial, and unconstrained), `cca()` and `ca()` of `dune`;
    easyCODA 0.50's unweighted `PCA(CLR(x))`; compositions 2.0-9's
    `princomp.acomp()`. The plan called the CoDa cross-check optional; it
    was done. The relations, for \eqn{n} samples, \eqn{p} parts and
    `vegan`'s constant \eqn{c = ((n-1)\,\mathrm{tot})^{1/4}}:
    form rows = `vegan` sites (1) × \eqn{\sqrt{(n-1)\mathrm{tot}}/c} =
    easyCODA `rowpcoord` × \eqn{\sqrt{p}}; form columns = species (1) /
    \eqn{c}; covariance rows = sites (2) × \eqn{\sqrt{n-1}/c} = `rowcoord` ×
    \eqn{\sqrt{(n-1)/n}}; covariance columns = species (2) ×
    \eqn{\sqrt{\mathrm{tot}}/c} = `colpcoord` × \eqn{\sqrt{n/(n-1)}};
    contribution columns = `colcoord` × \eqn{\sqrt{\mathrm{colmass}}} (the
    loadings); compositions' `sdev^2` = `eig`. robCompositions and coda.base
    were installed but not used: they add no biplot coordinates beyond these.
15. **Tests** (`test-plot.r`, "Biplot scalings", about 980 new
    expectations) and two `test-smoke.r` lines (a covariance biplot with
    links, and `rda()` in scaling 1).
16. **Docs:** `plot.ecodive_ord()`'s `@param scaling`, `@param links`, the
    "Compositional biplots" and "Scalings of rda(), cca() and ca()"
    sections, references and examples; "Biplots" paragraphs in `pca()` and
    `rpca()`; notes in `rpca_phylo()`, `rda()`, `cca()` and `ca()`; NEWS;
    `parameter_names.md` (`scaling` in §1 and the 2026-11-20 table, and
    `links`).

Reference results (largest relative differences, over every axis pair):
`vegan` 2.7.8 `scores()` of a clr RDA against the form and covariance
biplots, 6e-14; of `rda()` (partial), unconstrained `rda()`, `cca()` and
`ca()` in scalings 1, 2 and 3 (sites, species, constraints, centroids),
2e-14, 5e-14, 5e-14 and 3e-14; easyCODA 0.50, 2e-14; compositions 2.0-9,
8e-14. Link lengths over all axes against `sqrt(vlr(x))`, 5e-16; ray lengths
against the clr values' standard deviations, 6e-16.
