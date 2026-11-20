# `coda_dendrogram()`: balances of a sequential binary partition, as a tree

**Priority:** 62. **Tier:** clean room (`compositions::CoDaDendrogram()` is
GPL ≥ 2: its documentation was read, its code was not). Read
[README.md](README.md) first, [compositional.md](compositional.md), and
[log_ratio_transforms.md](log_ratio_transforms.md) for `ilr(basis =)`.


## What it is

A sequential binary partition (SBP) splits the parts in two, then splits
each group again, until every part is alone: \eqn{D - 1} splits, each
giving one **balance** (an ilr coordinate). The CoDa-dendrogram
(Egozcue & Pawlowsky-Glahn 2005; Pawlowsky-Glahn & Egozcue 2011) draws the
partition as a dendrogram and puts on it, for each balance:

1. **the partition** itself (dendrogram links between parts);
2. **the balance mean**, where the balance's axis meets the link to its
   parent;
3. **the balance variance**, as the length of that link; since the
   coordinates of an orthonormal basis have variances summing to the total
   variance, **the links decompose the total variance**;
4. optionally, **quantile box plots** (5%, 25%, 50%, 75%, 95%) on each
   balance's axis, and **several groups overlaid**.

Each balance axis spans the interval \eqn{(-u, u)} whatever its drawn
length, so boxes can be compared between balances ("as long as required to
connect the parts in a group. Each one represents an interval (−u, u),
where u is user defined", 2011, §5).

Verified from the open-access 2011 paper (§5, Figure 1) and compositions'
`CoDaDendrogram` Rd text ("Each coordinate is represented in a horizontal
axis, which limits correspond to the values given in the parameter range.
The vertical bar going up from each one of these coordinate axes represent
the variance of that specific coordinate, and the contact point the
coordinate mean"). The two orientations differ: the paper draws the parts
as rows on the left (balance axes vertical, variance bars horizontal);
compositions draws parts along the bottom (axes horizontal, bars vertical).
**ecodive uses the paper's orientation**, because taxon names then read
horizontally. **Reconstructed, not verified:** the exact placement rule
below (child axis offset from its parent by its variance) is inferred from
the two texts and Figure 1; check against the 2011 figure when building.

**Uptake:** Egozcue & Pawlowsky-Glahn (2005) 740 citations;
Pawlowsky-Glahn & Egozcue (2011) 110; but only 3 Europe PMC microbiome
papers use a CoDa-dendrogram (2019–21), none since. Its microbiome
descendants are the balance trees of gneiss (274) and PhILR (359), both
declining. Medium effort; distinctive; low demand. Honest ranking: last
among the compositional displays except `tree_plot()`.


## Names

- `coda_dendrogram()`: the method's name in snake_case (README rule 5),
  with "coda" as in `coda_lasso()`; compositions' `CoDaDendrogram()` has
  a different spelling, so no masking. No collisions in `NAMESPACE`.
- Draws immediately and returns the node summaries invisibly (the
  `_plot` / display rule in [ternary_plot.md](ternary_plot.md); "dendrogram"
  already names a picture).
- Arguments: `counts`, `basis` (the glossary's "ILR basis or SBP", as
  `ilr()`), `tree` (as `philr()`), `groups`, `data`, `margin`,
  `pseudocount`, `legend`, `col`, `labels`, `cpus`. New: `limits` (the
  \eqn{(-u, u)} of every axis; compositions calls it `range`, which shadows
  `base::range()` in reading), `type`.


## Design

### Signature

```r
coda_dendrogram(counts, basis = NULL, tree = NULL, groups = NULL, data = NULL,
                margin = 1L, pseudocount = NULL, limits = NULL,
                type = c('boxplot', 'mean'), labels = TRUE,
                legend = 'topright', col = NULL, cpus = n_cpus(), ...)
```

### The partition

- `basis`: an SBP, \eqn{(D-1) \times D} of 1 / −1 / 0, validated by
  `ilr_basis()` exactly as `ilr()` does (sequential, both signs per row).
  An orthonormal non-SBP basis is an error here: it has no tree.
- `tree`: a `phylo` whose tips include the features; `philr_tree()`
  prunes and resolves it as `philr()` does, and each internal node's
  first child is the numerator (philr's convention). Pass
  `attr(philr(x, tree), 'tree')` to reuse a PhILR tree. (Only
  `part_weights = 'uniform'` PhILR coordinates are plain balances; the docs
  say weighted ones are not drawn.)
- Neither: features clustered by `hclust(sqrt(vlr(counts)), 'ward.D2')`
  (\eqn{\sqrt{\tau_{ij}}} is the Euclidean distance between the features'
  sample-centred log vectors, scaled, so Ward is valid), and the merge tree
  turned into an SBP (each merge = one balance, left branch numerator).
  ecodive's choice; Pawlowsky-Glahn & Egozcue (2011, §6) suggest the
  variation array and the biplot to guide partitions, and compositions'
  example clusters `dist(t(clr(x)))` with Ward. Say so.
- Both given → error.

Internal helpers (`R/coda_dendrogram.r`): `sbp_from_hclust(merge)`,
`sbp_from_phylo(phy)`, and `sbp_nodes(sbp)` returning, per row, its parent
row, its two child groups, and the leaf order (preorder, numerator first).

### The numbers

- Balances: `ilr(counts, basis = sbp, margin, pseudocount)` (the shared
  pseudocount rule; bit-identical to `balance()` for each row).
- Per balance and per group (all samples when `groups = NULL`): `n`,
  `mean`, `var` (divisor \eqn{n-1}), `var_prop` = var / Σ var over the
  group's balances (= its total variance), and `quantile(type = 7)` at
  0.05, 0.25, 0.5, 0.75, 0.95 (the paper's \eqn{p_{0.05}, Q_1, Q_2, Q_3,
  p_{0.95}}).
- `limits = NULL`: \eqn{u} = the largest |5% or 95% quantile| over all
  balances and groups, rounded up with `pretty()`; or `c(-u, u)` given.
  Means or quantiles outside are clamped to the axis end and drawn with an
  open marker.

### Layout (paper's orientation)

- **Parts** are rows: leaf \eqn{k} (preorder, numerator side first) at
  \eqn{y = D - k + 1}, its label at the left (`labels = TRUE`; a count or
  `FALSE` as elsewhere).
- Every balance \eqn{b} has a **vertical axis** spanning its two children's
  attachment heights, \eqn{y_{\mathrm{den}}} to \eqn{y_{\mathrm{num}}}
  (a leaf's attachment height is its row; a balance's is the height of its
  mean, below), mapped linearly to \eqn{(-u, u)}: \eqn{-u} at the
  denominator child, \eqn{+u} at the numerator child.
- **Mean:** \eqn{y_b = y_{\mathrm{den}} + (\bar b + u)/(2u)\,(y_{\mathrm{num}} -
  y_{\mathrm{den}})}: the point where the balance's link leaves its axis.
- **Variance:** the link from the axis at \eqn{y_b} toward the parent's axis
  is horizontal, thick (`lwd = 3`), of length \eqn{\mathrm{var}(b)}. So the
  axis of the root sits at \eqn{x = W - \mathrm{var}(\mathrm{root})} and
  each child's axis at \eqn{x_{\mathrm{child}} = x_{\mathrm{parent}} -
  \mathrm{var}(\mathrm{child})}; \eqn{W} is the largest sum of variances
  along any root-to-node path, so every axis has \eqn{x \ge 0}. Leaves join
  their parent's axis by thin `'grey60'` lines from \eqn{x = 0} (no
  quantitative meaning).
- The x axis is labelled "Variance (total = 2.55)", ticks by `pretty()`
  measured **from the right** (0 at \eqn{W}).
- **Box plots** (`type = 'boxplot'`): on each axis, at the axis' x, a
  narrow box (width 2% of \eqn{W}) from Q1 to Q3, whiskers to 5% and 95%,
  a tick at the median, all in the axis' \eqn{(-u, u)} scale.
- **Groups:** the layout (axes' x and the means' y) is computed from **all
  samples pooled**; each group then adds, in its palette colour
  (`ord_plot_look()`), its mean (a point on the axis), its variance (a
  thick segment from that point toward the parent, of length its variance,
  offset vertically by a few pixels per group so overlays don't hide each
  other) and its box (side by side, left to right in level order). The
  paper overlays separate dendrograms instead ("the most variable sample has
  to be plotted first"); the pooled layout is ecodive's simplification, so
  that axes line up. Say so in the docs.
- Title "CoDa-dendrogram"; subtitle "Partition: Ward clustering of
  log-ratio variances" / "from the tree" / "given"; legend by
  `plot_legend()`. `plot_do()` passes `...` (`main`, `xlab`, `cex`, ...).

### Value

Invisibly, a `data.frame` with one row per balance × group: `balance` (the
SBP row name, else `b1`, `b2`, ... in preorder), `numerator`,
`denominator` (features joined by ", "), `group` (`NA` without groups),
`n`, `mean`, `var`, `var_prop`, `q05`, `q25`, `q50`, `q75`, `q95`, `x`,
`y` (pooled layout); attribute `sbp` (the partition used) and
`pseudocount`.


## Reference implementations

- **None ported.** compositions' Rd (read for behaviour) and the 2011 paper
  (open access) define the display; the numbers are `ilr()` and base
  `quantile()`.
- **Output-only (GPL):** `compositions::CoDaDendrogram(acomp(x), signary =
  sbp, type = 'lines')` for a visual check of a 5-part example; nothing to
  compare numerically beyond means and variances, which the tests check
  against identities.


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- A plotting function: returns its layout invisibly with a `plot`
  attribute `ecodive_step("CoDa-Dendrogram", partition =, ...)` (partition
  from `"basis"`, `"tree"` or `"ward"`).
- Title from the input's `transform` record; when the partition comes from
  a `philr()` result, use its `transform` record's weights in the
  subtitle.

## Tests (`test-coda_dendrogram.r`)

- Variances: Σ var(b) over balances = `sum(apply(clr, 2, var))` (total
  variance) at 1e-12; equals `attr(variation_array(x), 'total_var')` if
  item 53 has landed. Means equal `colMeans(ilr(x, sbp))` exactly.
- `sbp_from_hclust()` on a hand-built merge matrix gives the expected SBP;
  every SBP produced passes `ilr_basis()`'s checks. `sbp_from_phylo()` on
  `ex_tree` gives balances equal to `philr(ex_counts, ex_tree)`'s uniform
  coordinates at 1e-12 (column order mapped by node label).
- Layout: \eqn{x_{\mathrm{child}} = x_{\mathrm{parent}} - \mathrm{var}} exactly;
  all \eqn{x \ge 0}; leaves at rows 1..D; each mean's y within its axis.
- Groups: per-group means and variances equal hand computations; a
  singleton group has `var = NA` and no variance segment (one warning).
- `pdf(NULL)`: returns invisibly (`withVisible`), `type = 'mean'`, `labels
  = FALSE`, `limits = c(-2, 2)` with clamping; `.Random.seed` untouched.
- Errors: both `basis` and `tree`; a non-SBP basis; a tree missing features;
  `limits` not symmetric or not increasing; `type = 'x'`; fewer than 3
  features.
- Smoke line.


## Docs

- Rd in "Plotting" → "Compositional displays". Explain the four elements,
  the \eqn{(-u, u)} convention, that link lengths add up to the total
  variance, how to read group differences (means apart, variances similar →
  a discriminating balance; 2011, §5), the pooled-layout simplification and
  the default Ward partition. Link `variation_array()`, `ilr()`, `philr()`,
  `balance()`, `tree_plot()`. NEWS.


## References

- Egozcue, J. J., & Pawlowsky-Glahn, V. (2005). Groups of parts and their
  balances in compositional data analysis. *Mathematical Geology*, 37,
  795–828. \doi{10.1007/s11004-005-7381-9}
- Pawlowsky-Glahn, V., & Egozcue, J. J. (2011). Exploring compositional
  data with the CoDa-dendrogram. *Austrian Journal of Statistics*, 40,
  103–113. \doi{10.17713/ajs.v40i1&2.202}
- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal,
  C. (2003). Isometric logratio transformations for compositional data
  analysis. *Mathematical Geology*, 35, 279–300. \doi{10.1023/A:1023818214614}
- Morton, J. T., Sanders, J., Quinn, R. A., et al. (2017). Balance trees
  reveal microbial niche differentiation. *mSystems*, 2, e00162-16.
  \doi{10.1128/mSystems.00162-16}
