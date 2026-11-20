# `tree_plot()`: a base-R phylogenetic tree with values on nodes and tips

**Priority:** 63. **Tier:** first principles (ape, ggtree, phytools,
philr and metacoder are GPL or LGPL: none of their code is read; the
rectangular layout is textbook). Read [README.md](README.md) first, and
[philr.md](philr.md) and [rpca_phylo.md](rpca_phylo.md) for the trees they
return.


## What it is

ecodive reads trees (`read_tree()`) and computes per-node quantities, but
cannot draw a tree:

- `rpca_phylo()` returns **loadings for every node** (tips and clades,
  named `"n<id>"`) and the tree it used (`attr(x, 'tree')`, internal nodes
  labelled with those names);
- `philr()` returns **one balance per internal node** (columns named by
  node label) and its tree (`attr(x, 'tree')`), the PhILR / gneiss
  "balance tree" (Silverman et al. 2017; Morton et al. 2017).

`tree_plot()` draws a rooted tree as a rectangular **phylogram** (x =
distance from the root) or **cladogram** (tips aligned), and puts values on
nodes and tips as coloured, sized points. That is the figure PhILR and gneiss
papers show (balances, or their group differences, on the tree), and the
natural display of `rpca_phylo()`'s clade loadings. It is also the base
that a later heat-tree or tree + heatmap display would need.

**Uptake:** PhILR 359 citations (83 Europe PMC mentions 2022–26), gneiss 274
(27); both declining (0.45× and 0.43× normalised). Phylogenetic RPCA is
young. The demand is real but small; the cost (a tree drawer from scratch)
is the largest of the compositional displays, hence last.


## Names

- `tree_plot()`: draws immediately and returns the layout invisibly (the
  `_plot` rule in [ternary_plot.md](ternary_plot.md)). Not `plot.phylo()`:
  that is ape's S3 method, and registering ecodive's would mask ape's for
  every user who has both. Not `plot_tree()` (phyloseq's, GPL;
  `ggtree`-style) nor `treePlot()` (phylobase). No collisions in
  `NAMESPACE`.
- Arguments: `tree` (as `philr()`, `rpca_phylo()`, the UniFrac metrics),
  `groups`, `data`, `axes` (one axis of an ordination, as
  [log_ratio_plots.md](log_ratio_plots.md)'s `log_ratio_features()`),
  `labels` (a count), `legend`, `col`. New: `values`, `type`
  (`'phylogram'` / `'cladogram'`), `node_labels` (a count). Add to
  [parameter_names.md](parameter_names.md) §1.


## Design

### Signature

```r
tree_plot(tree, values = NULL, type = c('phylogram', 'cladogram'),
          groups = NULL, data = NULL, axes = 1L, labels = NULL,
          node_labels = 5L, legend = TRUE, col = NULL, ...)
```

### Inputs

- `tree`, one of:
  - a `phylo` (`read_tree()`, ape, or the `tree` attribute of a result);
  - an `rpca_phylo()` ordination: tree = `attr(x, 'tree')`, values =
    `attr(x, 'loadings')[, axes]` (one axis; error if `length(axes) != 1`);
  - a `philr()` matrix (`attr(x, 'method') == 'philr'`): tree =
    `attr(x, 'tree')`; values = the column means of the coordinates, or,
    with `groups` of exactly two levels (`per_sample()` / `match_rows()`
    on the matrix's sample names), mean of the second level minus mean of
    the first (a balance "effect"; the subtitle says which). Respect the
    result's orientation (`margin = 2` results are transposed: coordinates
    are rows; detect from `attr(x, 'tree')$node.label` matching
    `rownames(x)`).
- `values`: a numeric vector named by tip labels and/or internal node labels
  (`tree$node.label`); `NA` for none. Names not in the tree → error listing
  up to 5. Unnamed of length `Ntip + Nnode` → by node number (ape order:
  tips 1..Ntip, then internal nodes).

### Validation (`tree_check()`, internal)

`inherits(tree, 'phylo')`; `edge` an integer two-column matrix; `Ntip =
length(tip.label)` ≥ 2; nodes numbered 1..(Ntip + Nnode); exactly one
**root** (a parent never a child; ape numbers it Ntip + 1 but don't assume
it); every other node a child exactly once; no cycles (the preorder walk
below visits every node exactly once, else "`tree` is not a valid rooted
tree"). Polytomies are allowed (unlike `philr()`).

### Layout (first principles; `tree_layout()`, internal, `R/plot_tree.r`)

1. **Children:** `split(edge[, 2], factor(edge[, 1], levels = 1:N))`, in
   edge-matrix order (that order sets the tips' vertical order, as ape's
   "cladewise" order does).
2. **Preorder** by an explicit stack (no recursion, so 100,000-tip trees
   don't overflow R's stack): push the root; pop a node, append it, push its
   children in reverse. **Postorder** = reverse of preorder with children
   reversed, or simply process preorder backwards (every child appears after
   its parent).
3. **y (vertical):** tips in preorder get \eqn{y = N_{tip}, N_{tip}-1,
   \dots, 1} (first tip at the top). Internal nodes, in reverse preorder:
   \eqn{y = (y_{\text{first child}} + y_{\text{last child}})/2}.
4. **x (horizontal):**
   - phylogram: \eqn{x_{\text{root}} = } `tree$root.edge` (or 0), then in
     preorder \eqn{x_c = x_p + \ell_{pc}} with \eqn{\ell} the branch length.
     Needs `edge.length` with no `NA` and none negative; if `edge.length` is
     `NULL` (as `read_tree()` leaves all-zero trees), fall back to the
     cladogram with the subtitle "no branch lengths: cladogram".
   - cladogram: in reverse preorder, \eqn{h = 0} for tips and \eqn{1 +
     \max h_{\text{children}}} for internal nodes; \eqn{x = h_{\text{root}} - h}
     (all tips at the right).
5. **Segments:** one horizontal per edge, \eqn{(x_p, y_c) \to (x_c, y_c)};
   one vertical per internal node, \eqn{(x_p, \min y_c) \to (x_p, \max y_c)};
   a root edge, if any, from \eqn{(0, y_{\text{root}})}. Each drawn with a
   **single vectorised `segments()` call**.

### Drawing

- `plot.new()`, `plot.window(xlim = c(0, max(x) * (1 + label room)), ylim =
  c(0.5, Ntip + 0.5))`, label room from `strwidth()` of the tip labels;
  branches in `'grey40'`, `lwd = 1`.
- **Values:** points at the nodes (and tips) that have one, `pch = 21`,
  fill on a diverging scale symmetric about 0 (`hcl.colors(65, 'Blue-Red
  3')`, or `col` as the colour ramp), size `cex` from 0.6 to 2.5 in
  proportion to |value| (sqrt-scaled). With `legend = TRUE`, a colour key
  in the right margin (`dist_plot_key()`), titled from the source
  ("Loading on RPCA1", "Mean balance", "Balance: B − A").
- **Tip labels:** `labels = NULL` → all when `Ntip` ≤ 60, none above; a
  number → that many, those with the largest |value| (or the first, without
  values); `TRUE` / `FALSE`. Drawn with `text(pos = 4, cex = 0.7)`.
- **Node labels:** the `node_labels` internal nodes with the largest
  |value| are named with `plot_label()` (their `node.label`).
- Title from the source: "Tree", "Phylogenetic RPCA loadings, RPCA1",
  "PhILR balances". `plot_do()` for `main`, `sub`, `cex` ... in `...`.
- Values of 20,000-tip trees draw in well under a second: no per-node R
  loops beyond the preorder walk (vectors indexed by node number).

### Value

Invisibly, a `data.frame` with one row per node, in node-number order:
`node`, `label` (tip label or node label), `parent` (`NA` at the root),
`is_tip`, `x`, `y`, `value`; so users can add to the plot (`text(lay$x,
lay$y, ...)`).

### Reuse

`coda_dendrogram()` (item 62) needs a preorder walk and child lists for
SBPs built from trees; whichever item is built first writes `tree_check()`
and the walk in `R/plot_tree.r`, and the other reuses them.
`coda_dendrogram()`'s variance-driven x positions are its own.


## Reference implementations

- **None read.** The rectangular layout (tips evenly spaced in traversal
  order, internal nodes at the midpoint of their children, x by cumulative
  branch length or by height) is the standard construction (see e.g.
  Felsenstein 2004 on drawing trees; chapter not checked). ape's `plot.phylo()`
  (GPL) may be **run** to compare pictures by eye; its node coordinates
  (`get("last_plot.phylo", envir = .PlotPhyloEnv)`) can be an output-only
  check that tip y and node x agree on `ex_tree` — optional, not a fixture.


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- A plotting function: no record of its own. Title from what is drawn:
  `step_name(attr(x, 'ord'))` for an `rpca_phylo()` result ("Phylogenetic
  Robust Aitchison PCA: loadings on axis 1"), the `transform` record's
  name for a `philr()` result.
- Return the layout (node coordinates) invisibly, so values can be
  overlaid by a later call.

## Tests (`test-tree_plot.r`)

- `tree_layout()` exactly, on `ex_tree` and on hand-written Newick strings
  via `read_tree()`: tip y in preorder; each internal y the midpoint of its
  first and last child; phylogram x equals the summed branch lengths from
  the root (`root.edge` added); cladogram tips all at `max(x)` and each
  internal x = its parent's + 1; a polytomy; a caterpillar tree of 5,000
  tips (no recursion error).
- Values: by tip and node name, by position; unknown names error;
  `rpca_phylo(ex_counts, ex_tree, n_axes = 2)` with `axes = 2` uses the
  second column; a `philr()` result with and without two-level `groups`
  (difference of means equals a hand computation; three levels → error).
- `pdf(NULL)`: returns the layout invisibly; `type = 'cladogram'`;
  `labels = FALSE`; no branch lengths (fallback note); `.Random.seed`
  untouched.
- Errors: not a `phylo`; two roots; a cycle; negative branch lengths in a
  phylogram; one tip.
- Smoke line.


## Docs

- Rd in "Plotting" → "Compositional displays" (with a pointer from
  `rpca_phylo()` and `philr()` Rd "See also"). Explain both layouts, how
  values are matched, the PhILR group-difference reading (which balances
  differ between groups; positive = numerator clade relatively more
  abundant in the second level), and that it is not ape's `plot.phylo()`.
  NEWS.


## References

- Silverman, J. D., Washburne, A. D., Mukherjee, S., & David, L. A.
  (2017). A phylogenetic transform enhances analysis of compositional
  microbiota data. *eLife*, 6, e21887. \doi{10.7554/eLife.21887}
- Morton, J. T., Sanders, J., Quinn, R. A., et al. (2017). Balance trees
  reveal microbial niche differentiation. *mSystems*, 2, e00162-16.
  \doi{10.1128/mSystems.00162-16}
- Martino, C., McDonald, D., Cantrell, K., et al. (2022). Compositionally
  aware phylogenetic beta-diversity measures better resolve microbiomes
  associated with phenotype. *mSystems*, 7, e00050-22.
  \doi{10.1128/msystems.00050-22}
- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer.
