# `philr()` and `philr_distance()`: phylogenetic isometric log-ratio transform

**Priority:** 32. **Tier:** formula, clean room.
philr (Bioconductor, GPL-3) is output-only; the formulas are in the
open-access paper (Silverman et al. 2017, eLife, CC-BY). scikit-bio's
`tree_basis()` (BSD-3) is a permissive reference for the unweighted basis,
but its orientation and node order differ from philr's (see below). Read
[README.md](README.md) first, and [compositional.md](compositional.md) for
the conventions the compositional plans share.


## What it is

PhILR (Silverman et al. 2017) builds an ILR basis from a rooted binary
phylogeny: each **internal node** gives one coordinate, the balance between
the taxa of its two child clades. The coordinates are then real-valued,
orthogonal (before weighting), and each one has an evolutionary meaning
("clade A vs its sister clade B"). Two weightings make it more than a
rotated clr:

- **Part (taxon) weights** \eqn{p_j} down-weight rare, noisy taxa through a
  weighted clr (Egozcue & Pawlowsky-Glahn 2016).
- **Balance (ILR) weights** scale each balance by the branch lengths below
  its node, so splits between distant clades count more.

Euclidean distance on PhILR coordinates is a compositional, phylogenetic
beta diversity, the CoDa counterpart of UniFrac. **With uniform weights it
equals the Aitchison distance** (any orthonormal ILR does), so the
weighting is what makes it phylogenetic.

Two functions, because they return different things and the distance
joins `beta_div()`:

| Function | Returns |
|:---------|:--------|
| `philr()` | the PhILR coordinates (samples × internal nodes) |
| `philr_distance()` | Euclidean distances between samples' PhILR coordinates, with the paper's weights; a phylogenetic beta metric |


## Names

- `philr()`: the package's name, same meaning (no conflict under rule 6;
  masks `philr::philr` when both are attached; say so).
- `philr_distance()`: not a second `philr`. In `METRICS`:
  `PhILR Distance, philr_distance, phylo YES, weighted YES, int_only NO,
  true_metric YES, beta, alt_ids: philr phylogeneticilr`. The alias
  `philr` lets `beta_div(x, 'philr', tree = t)` work. Check the aliases
  with the README snippet.
- Argument names are philr's, in snake_case: `part_weights`,
  `ilr_weights`.


## Design

### Signatures

```r
philr(counts, tree = NULL, margin = 1L, pseudocount = NULL,
      part_weights = 'uniform', ilr_weights = 'uniform', cpus = n_cpus())

philr_distance(counts, tree = NULL, margin = 1L, pseudocount = NULL,
               pairs = NULL, cpus = n_cpus())
```

- `tree`: a rooted, **strictly binary** `phylo`, found as the UniFrac
  functions find it (`validate_counts()` pulls it from phyloseq, rbiom or
  TreeSummarizedExperiment objects; `validate_tree()` matches tips to
  features and pads absent tips with zero columns).
  - Polytomies are an error that suggests `ape::multi2di()` (philr's
    documented advice; it adds zero-length branches). Don't resolve them
    silently.
  - Unrooted trees: error (an unrooted tree's root is a trifurcation).
- `pseudocount`: the shared rule (compositional.md). philr's own default
  is 0 (zeros are an error there); the paper used 1. Fixtures pass an
  explicit value.
- `part_weights` (philr's `part.weights`): `'uniform'` (default),
  `'gm.counts'`, `'anorm'`, `'anorm.x.gm.counts'`, `'enorm'`,
  `'enorm.x.gm.counts'`, or a numeric vector named by feature.
- `ilr_weights` (philr's `ilr.weights`): `'uniform'` (default), `'blw'`,
  `'blw.sqrt'`, `'mean.descendants'`, or a numeric vector named by
  internal node.
  - These are options of one method in philr and the paper, so they stay
    arguments (README rule).
- **`philr_distance()` has no weight arguments.** `beta_div()` passes
  arguments with `mget(metric$params, environment())`, so every formal of a
  metric must be a formal of `beta_div()`. Rather than widen `beta_div()`,
  `philr_distance()` is *the* PhILR distance of the paper: part weights
  `'enorm.x.gm.counts'`, ILR weights `'blw.sqrt'` (Silverman et al. 2017,
  Methods; the package's docs call `enorm.x.gm.counts` "often gives good
  results"). Any other weighting is `euclidean(philr(x, tree, ...))`.
  Document the identity
  `philr_distance(x, t) == euclidean(philr(x, t, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt'))`.
  - **Ask the user** whether they'd rather widen `beta_div()` (adding
    `part_weights` / `ilr_weights` formals that other metrics ignore) and
    default `philr_distance()` to philr's uniform weights. Recommended: the
    fixed-weights design above, because with uniform weights the "PhILR
    distance" is just `aitchison()`.

### Formulas (Silverman et al. 2017, Methods; verified from the article)

For one sample with closed composition \eqn{x} and part weights \eqn{p}:

1. **Shift:** \eqn{y_j = x_j / p_j}.
2. For internal node i with "+" clade \eqn{C^+} and "−" clade \eqn{C^-}:
   \eqn{n^\pm = \sum_{j \in C^\pm} p_j}, and the weighted geometric mean
   \eqn{\log g_p(y^\pm) = \sum_{j \in C^\pm} p_j \log y_j / n^\pm}.
3. **Balance:** \eqn{y^*_i = \sqrt{n^+ n^- / (n^+ + n^-)}\;\bigl(\log g_p(y^+) - \log g_p(y^-)\bigr)}.
4. **ILR weight:** \eqn{y^*_i \leftarrow y^*_i \, f(d^+, d^-)}, where
   \eqn{d^\pm} are the branch lengths of node i's two child edges, and
   \eqn{f} is 1 (`uniform`), \eqn{d^+ + d^-} (`blw`),
   \eqn{\sqrt{d^+ + d^-}} (`blw.sqrt`), or \eqn{d^\pm} plus each child's
   mean distance to its tips (`mean.descendants`). philr's
   `calculate.blw` first gives zero-length **tip** branches the smallest
   non-zero branch length in the tree.
5. **Part weights** (computed from the whole table, so one sample's
   coordinates depend on the others; say so):
   - `gm.counts`: the geometric mean over samples of each taxon's counts
     (the paper writes \eqn{\prod_k (c_{jk} + 1)^{1/N}}; **confirm** whether
     philr adds 1);
   - `enorm`: the Euclidean norm over samples of the taxon's closed
     relative abundances; `anorm`: its Aitchison norm (**confirm the
     definition**);
   - `*.x.gm.counts`: the product of the two.

**Computation is O(D) per sample**, no basis matrix needed: with
\eqn{w_j = p_j \log y_j}, one postorder pass accumulates
\eqn{\sum p_j} and \eqn{\sum w_j} per clade; each node's balance follows
from its two children's sums. Parallel over samples. This is also how
`philr_distance()` works: coordinates per sample (dense, D − 1 values),
then the squared-Euclidean kernel of `src/beta_div.c` over sample pairs.
Respect `pairs` as the other metrics do.

### Orientation and order (match philr's outputs)

- philr builds its sign matrix with `phylo2sbp()`; its docs say the "+"
  clade is "given by the default of `phangorn::Children`", i.e. (inferred)
  the child listed **first** in `tree$edge` order. scikit-bio's
  `tree_basis()` instead puts the **second** child on top and orders rows
  level-order. Match **philr**, and confirm on fixtures.
- Coordinates (columns for `margin = 1`) are the internal nodes in ape
  node-number order (root first: node `Ntip + 1`). Name them by
  `tree$node.label` when the tree has internal labels, else `n1`, `n2`, ...
  (philr's vignette labels nodes this way with
  `ape::makeNodeLabel(method = "number", prefix = "n")`). Check how philr
  names unlabelled trees, and ecodive's `read_tree()` node numbering
  (must be ape-compatible preorder) on the fixtures.

### Results

- `philr()`: a dense base `matrix`, input orientation; attributes `method =
  'philr'`, `part_weights` (the numeric vector used, named by feature),
  `ilr_weights` (named by node), `pseudocount`. Not the basis (D × (D − 1)
  can be huge); document how philr's `return.all` maps to these.
- `philr_distance()`: a `dist`, as every metric.


## Reference implementations

- **Output-only:** philr 1.38.0 (GPL-3). `data-raw/expected/philr.r` runs
  `philr(x, tree, part.weights = , ilr.weights = , pseudocount = )` for
  every weighting, on a tree with labelled and one with unlabelled
  internal nodes, plus `dist()` of the paper's setting. Compare at 1e-10.
  Read only philr's manual and vignette, never its R code.
- **Paper:** Silverman, Washburne, Mukherjee & David (2017), eLife,
  Methods (CC-BY): every formula above.
- **Permissive cross-check:** scikit-bio `tree_basis()` (BSD-3): with
  uniform weights, `philr()` coordinates equal `clr(x) %*% t(tree_basis)`
  up to row order and sign.
- **To confirm from outputs** (the docs don't say): `gm.counts`' "+1",
  `anorm`'s definition, the "+" child, names on unlabelled trees.


## Tests

- Fixtures as above.
- Uniform weights: Euclidean distances on `philr()` equal `aitchison()`
  (same pseudocount); coordinates equal `ilr()` with the matching SBP basis.
- `philr_distance()` equals `euclidean(philr(..., paper weights))`, and
  `beta_div(x, 'philr', tree = t)` gives it.
- Scale invariance; identical for `cpus` 1, 2 and 4; sparse and dense
  agree; `pairs` respected.
- Errors: no tree; polytomy (message suggests `ape::multi2di()`); tips
  missing from the tree; negative branch lengths; unknown weight names;
  weight vectors with missing names.
- Smoke-test lines.


## Docs

- `philr()` in "Compositional Data" → "Log-ratio transforms";
  `philr_distance()` under Beta Diversity → Phylogenetic metrics, next to
  UniFrac, and in `list_metrics(phylo = TRUE)`.
- Explain the weights, the uniform-weights identity with Aitchison, and
  `ape::multi2di()` for polytomies.
- NEWS.


## References

- Silverman, J. D., Washburne, A. D., Mukherjee, S., & David, L. A.
  (2017). A phylogenetic transform enhances analysis of compositional
  microbiota data. *eLife*, 6, e21887.
- Egozcue, J. J., & Pawlowsky-Glahn, V. (2016). Changing the reference
  measure in the simplex and its weighting effects. *Austrian Journal of
  Statistics*, 45, 25–44.
- Washburne, A. D., Silverman, J. D., Leff, J. W., et al. (2017).
  Phylogenetic factorization of compositional data yields lineage-level
  associations in microbiome datasets. *PeerJ*, 5, e2969.
