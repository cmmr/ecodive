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


## Decisions made while implementing (2026-11-19)

Asked of the user at the start:

1. **`philr_distance()`'s weights:** the paper's, fixed
   (`'enorm.x.gm.counts'`, `'blw.sqrt'`), with no weight arguments;
   `beta_div()` is unchanged (recommended option).
2. **Non-binary trees: resolved automatically** (the user chose this over
   the recommended error). Each node with children c_1 .. c_k (k > 2) keeps
   c_1 and gets a new child node, joined by a zero-length branch, holding
   c_2 .. c_k, and so on down a ladder, exactly as
   `ape::multi2di(random = FALSE)` does (checked: same coordinates, same
   column order). An ape-unrooted tree (three children at the root) is
   resolved the same way. A warning gives the number of nodes resolved.
3. **Tips that aren't features: pruned** (recommended). Internal nodes
   without kept tips go; nodes left with one child are collapsed into it,
   their lengths added (a single-child root is dropped with its length), as
   `ape::keep.tip()` / `drop.tip()` do (checked against philr on
   `keep.tip()`'s tree). No warning. Features missing from the tree are an
   error.
4. **Zeros:** the shared rule (recommended): `pseudocount = NULL` is none
   without zeros, else half the smallest non-zero value with a warning.
   `pseudocount = 0` is refused, as in `ilr()`.

Made while implementing:

5. **Weight options:** all of philr's (six part weightings, four ILR
   weightings) plus numeric vectors. Names may use `_` for `.`
   (`'enorm_x_gm_counts'`); matching is exact, not by prefix. Numeric part
   weights must be positive (one per feature, by name, extra names
   ignored, or by position); numeric ILR weights non-negative (by
   coordinate name or position). Computed `anorm` weights of 0 (a feature
   constant across samples, or one sample) are an error; philr would give
   NaN.
6. **Column order:** the input's internal node numbers (ape order, root
   first) when the tree is used unchanged, which is philr's order even when
   the numbers aren't in preorder (fixture `renumbered`); when steps 2-3
   changed the tree, preorder of the new tree, which is the numbering that
   `multi2di()` and `keep.tip()` produce.
7. **Coordinate names:** the node labels when every node used has a
   non-empty one and the names are distinct; nodes added at polytomy `X`
   are `X.1`, `X.2`, ... down the ladder. Otherwise `n1`, `n2`, ... in
   column order (philr itself returns no column names for an unlabelled
   tree; `n<k>` is the vignette's `makeNodeLabel(prefix = 'n')` naming).
   Bootstrap values as labels are usually duplicated, so they fall back to
   `n<k>`.
8. **Attributes:** `method = 'philr'`, `part_weights` (named by feature, in
   the table's order), `ilr_weights` (named by coordinate), `tree` and
   `pseudocount`. `tree` is the binary tree used, as a `phylo`: tip f is
   feature f, internal node D + k is column k (labelled with the column
   names), edges in cladewise preorder with each node's first child the
   "+" clade, and its branch lengths as pruned (ladder branches 0, tip
   zeros not replaced). `philr(x, attr(res, 'tree'))` gives identical
   coordinates (tested). philr's `V` (D × (D − 1)) is not kept.
9. **Branch lengths** are needed only for `'blw'`, `'blw.sqrt'` and
   `'mean.descendants'` (missing or NA lengths, or none positive, is an
   error for those only). Negative or non-finite lengths are always an
   error. A root edge is ignored. ecodive gives no message when it replaces
   zero-length tip branches (philr prints a note); the docs say so.
10. **Computation:** `src/philr.c`, written for this item (MIT). The tree is
    prepared once in C (`C_philr_tree()`: checks, pruning, collapsing,
    resolution, preorder numbering, `blw` and `mean.descendants` weights;
    linear in the tree's size, no recursion). Coordinates
    (`C_philr_coords()`) take one sample per task: logs into a per-thread
    scratch vector of 2D doubles, `w_j = p_j (log x_j − log p_j)`, then
    clade sums from the tips up (reverse preorder), so O(D) per sample, no
    basis, no BLAS, and identical for any `cpus` (tested 1–4). Per-feature
    sums for the part weights come from `C_philr_stats()` (one pass over
    the stored values plus the unstored zeros' share; a second for the
    Aitchison norm). Sparse and dense input agree to 1e-14 (the order of
    the part-weight sums differs). scikit-bio's `tree_basis()` was used as a
    cross-check only, not ported, so `inst/COPYRIGHTS` and `Authors@R` are
    unchanged.
11. **`philr_distance()`** computes the coordinates, then `C_beta_div()`'s
    Euclidean kernel on them, so `pairs` and `clust` work as for every
    metric, and the `dist` record carries the pseudocount used. It is
    identical to `euclidean(philr(x, tree, 'enorm.x.gm.counts',
    'blw.sqrt'))` (tested). In `METRICS` with the plan's row; the aliases
    `philr` and `phylogeneticilr` collide with nothing in any registry
    (README snippet run). `philr()` itself joins no wrapper.
12. **Speed** (4-core container, -O2): 2,000 samples × 5,000 taxa, sparse
    (5% non-zero), paper weights: `philr()` 0.23 s on 1 thread, 0.15 s on
    4; philr 1.28.0 161 s on the dense table (same coordinates to 7.6e-15).
    `philr_distance()` is dominated by the O(n² D) Euclidean kernel on the
    dense coordinates: 12.9 s on 1 thread, 2.8 s on 4 (as `aitchison()`).
13. **`log_ratio_input()`** gained a `check` hook, called with the feature
    names before the pseudocount is resolved, so that tree errors come
    before the zeros warning.

How the plan's "confirm" items were settled (all from philr 1.28.0's
outputs, never its code):

- **The "+" clade** is the first child in `tree$edge` order (philr's
  documented `phangorn::Children` default): on every fixture, including a
  tree whose internal nodes are numbered out of preorder. scikit-bio's
  `tree_basis()` puts the second child on top: its coordinates are ours
  negated, rows matched by node name.
- **`gm.counts`** is the geometric mean over samples of `counts +
  pseudocount`, with **no extra +1** (candidates with and without +1
  compared on a table with pseudocounts 0.5 and 1: exact without).
- **`enorm`** is the Euclidean norm of each feature's column after closing
  the **column** (over samples) to 1, not each sample: `sqrt(sum(x_j^2)) /
  sum(x_j)` with `x = counts + pseudocount`.
- **`anorm`** is the Aitchison norm of each feature's column over samples,
  `sqrt(sum((log x_kj − mean_k log x_kj)^2))` (closure cancels).
- **Part weights are not closed** (philr's docs: "should not be closed");
  uniform weights are 1, so the uniform transform is the ordinary ilr.
- **The balance formula** (paper, Methods, with the shift `y = x / p`) holds
  for every weighting to ~1e-15.
- **`blw`** is the sum of the two child branch lengths, after zero-length
  **tip** branches get the smallest positive branch length of the whole
  tree, internal branches included (fixture `unlabelled`, whose smallest
  positive branch is internal). **`mean.descendants`** adds each child's
  mean distance to its tips, computed with those replaced lengths.
- **Names on unlabelled trees:** philr returns none (see 7).
- **`read_tree()`'s numbering** is ape's: the scikit-bio fixtures are read
  with `read_tree()` in the tests and match philr's results from
  `ape::read.tree()` column for column.
- philr 1.28.0 itself accepts non-binary and unrooted trees without
  complaint (its sign matrix is then not an orthonormal basis) and returns
  NaN when the tree has tips that aren't columns; it errors on a numeric
  `part.weights` vector (`if (part.weights == "uniform")` on a vector), so
  numeric weights are tested by giving back the computed ones.

**Clean-room note.** That last error message, printed by `try()` while
exploring, included the deparsed condition of philr's part-weights `if`
chain. It appeared after the part-weight formulas above had been settled
from outputs (gm.counts, enorm and anorm candidates compared numerically in
the earlier scripts); nothing was taken from it, and later scripts printed
only `conditionMessage()`.

Reference results (fixtures `data-raw/expected/philr.r` and `philr.py`):

- **philr 1.28.0** (Bioconductor 3.18, the newest for R 4.3.3; needed
  tidytree 0.4.6 from the CRAN archive for treeio to load), every
  weighting (24 combinations) on labelled and unlabelled trees, plus a
  renumbered tree, a polytomy / unrooted / pruned tree (against philr on
  `multi2di(keep.tip(...), random = FALSE)`) and two Newick trees:
  coordinates to 6.9e-15 relative (2.1e-14 absolute), part and ILR weights
  to 1e-12, `dist()` of the paper's coordinates to 4.4e-16 relative.
- **scikit-bio 0.7.4 `tree_basis()`**: `clr %*% t(basis)` equals the
  negated uniform coordinates to 3.1e-16 relative (rows matched by name).
- **Identities:** uniform-weight distances equal `aitchison()` (1.8e-15),
  coordinates equal `ilr()` with the tree's SBP (1e-12 tested).

For `dpcoa()` (item 33): `C_philr_tree()` is a general "binary tree for
these features" preparation (prune to features, collapse single-child
nodes summing lengths, resolve polytomies, internal nodes numbered in
preorder so every child's id is above its parent's, tips mapped to feature
indices, cladewise edge list), and `philr_tree()` wraps it with the name and tip checks;
patristic distances or clade sums for DPCoA can start from its `c1`, `c2`
and `edge_len`.
