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
