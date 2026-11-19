# `rpca_phylo()`: phylogenetic robust Aitchison PCA

**Priority:** 19. **Tier:** port (gemelli 0.0.13,
`phylogenetic_rpca()`, BSD-3; its tree code comes from scikit-bio 0.7.4,
BSD-3, and iow 1.0.8 `bp`, BSD-3), built on ecodive's existing `rpca()`.
Read [README.md](README.md) first, and [compositional.md](compositional.md)
for the conventions the compositional plans share.


## What it is

Phylogenetic RPCA (Martino et al. 2022, gemelli) is `rpca()` on an
expanded table whose features are **every node of the tree**: each node
holds the summed counts of the tips below it, weighted by its branch
length. Clades that shift together then show up as one feature, and the
result is the compositional counterpart of a UniFrac PCoA, with feature
loadings for every clade. It's QIIME 2's `gemelli phylogenetic-rpca`.

**Name:** `rpca_phylo()`, so it sorts with `rpca()` (rule 2), not
gemelli's `phylogenetic_rpca`. Aliases in `COUNT_ORDS`:
`phylogeneticrpca`, `phylorpca`, `phylogeneticdeicode`. It joins
`count_ord()`: add the row
`Phylogenetic Robust Aitchison PCA, rpca_phylo, phylogeneticrpca phylorpca`.
`tree` passes through `count_ord()`'s `...`; it is no prefix of `counts` or
`method`.


## Design

### Signature

```r
rpca_phylo(counts, tree = NULL, margin = 1L, n_axes = 3L, max_iters = 5L,
           tol = 1e-5, min_sample_count = 0, min_feature_count = 0,
           min_feature_frequency = 0, min_depth = 0L, distances = FALSE,
           node_counts = FALSE, cpus = n_cpus())
```

- Every argument shared with `rpca()` means exactly the same (same
  defaults, same validation via `ord_arg()`), as gemelli shares them.
  - Note: gemelli **ignores** `max_iterations` in both `rpca` and
    `phylogenetic_rpca` (OptSpace always runs 5). ecodive's `rpca()`
    already documents how it handles this; do the same.
- `tree`: as the UniFrac functions (`validate_counts()` /
  `validate_tree()`).
- `min_depth` (gemelli's): drop internal nodes with at most this many
  descendant tips. Its net effect (from gemelli's `tree_topology_filter`
  then `prune`) is to **delete the internal-node features with ≤
  `min_depth` tips**; the other nodes' counts are unchanged. Default 0
  removes nothing.
- `node_counts = FALSE`: if `TRUE`, also return gemelli's counts-by-node
  table (raw subtree sums, not multiplied by branch length) as attribute
  `node_counts` (sparse `dgCMatrix` if Matrix is available, else a base
  matrix).
- `distances`: as `rpca()` (attribute `dm`).

### Algorithm (gemelli 0.0.13, read from source)

1. **Filter the tip table** exactly as `rpca()` (its filters, in its
   order).
2. **Shear and collapse the tree:** keep the root, the tips present in
   the filtered table, and their ancestors; remove single-child nodes,
   adding their length to the child. "Present" means a non-zero total
   after filtering (`validate_tree()` pads absent tips with zero columns;
   don't count those). The root is always kept, even with one child.
   - Edge case: if the root has one present child, scikit-bio 0.7.4's
     `prune` merges them and the Newick root length enters a real
     feature. ecodive's `read_tree()` discards root lengths, so this can't
     match; document it.
3. **`min_depth`:** count descendant tips per node (one pass) and drop
   internal nodes with ≤ `min_depth` tips.
4. **Count expansion:** each node (tips, internal nodes, root) gets the
   sum of its tips' counts. Its branch length is its own edge's; missing
   lengths (the root's) are 0. Node features are named by tip label, or
   `n<id>` for internal nodes, where gemelli numbers them by scikit-bio's
   `assign_ids` (postorder, root last). **Match gemelli's naming.** (Using
   the tree's own internal labels instead would be a documented
   difference; only if the user asks.)
   - gemelli casts counts to **int64** (proportions become all zeros).
     ecodive should **not** truncate: document it as a difference, and
     use integer inputs in fixtures (recommended by the research; ask the
     user if parity on non-integers matters).
   - Nodes with zero total count are dropped. **Zero-length nodes are
     not** (a scalar-vs-array quirk in scikit-bio's code that gemelli's
     tests confirm: the root, length 0, is kept).
5. **Robust clr with branch lengths:** value \eqn{\log(c \cdot b)} minus
   the sample's mean of \eqn{\log(c \cdot b)} over observed entries. An
   entry is missing if its count is 0 **or its length is 0**, so the root
   and zero-length nodes become **columns with no observed values** that
   still count in the column total m (they affect OptSpace's scaling,
   the final centering, and get non-zero loadings).
6. **OptSpace, then PCA:** `rpca()`'s code path unchanged.

### C

- **Refactor `C_rpca`** into (a) a filtering stage and (b) a core that
  takes a prebuilt sparse observed-entry table plus the **total** column
  count. Today `C_rpca` drops features with no non-zero count; the core
  must allow all-missing columns to match gemelli (step 5). Plain
  `rpca()` keeps its current behaviour through (a).
- **Tree work** reuses `src/ecotree.c` (`new_ecotree()` gives each non-root
  node its parent and edge length) and the leaf-to-root accumulation
  already in `src/unifrac.c`:
  - mark present nodes, count present children, collapse single-child
    chains (summing lengths), drop `min_depth` nodes;
  - a threaded sparse expansion, samples × nodes, walking up from each
    sample's non-zero tips, storing \eqn{\log(c \cdot b)} for observed
    entries (and raw sums when `node_counts = TRUE`).
- Size estimate (from the research): 350–500 new C lines, 100–150 lines
  of `rpca.c` refactoring, ~250 R lines.

### Result

`ecodive_ord` via `ord_result()`, exactly `rpca()`'s structure
(`method = 'rpca_phylo'`, prefix `"RPCA"`, `eig`, `prop`, `trace`,
`loadings`, `cost`, `n_iters`, optional `dm`), with these differences:

- `loadings` has one row per kept **node** (tips and internal nodes).
- New attribute `tree`: the sheared, collapsed (and `min_depth`-filtered)
  `phylo`, internal nodes labelled with the names used in `loadings`, so
  loadings can be mapped onto it.
- `node_counts` when requested.

Out of scope: gemelli's optional **taxonomy** output (tax2tree consensus
names per node; tax2tree's license notice is defective, so it would need a
clean-room rewrite) and **phylogenetic CTF** (see [ctf.md](ctf.md)).


## Reference implementations

- gemelli 0.0.13: `phylogenetic_rpca()` (`gemelli/rpca.py`),
  `bp_read_phylogeny()`, `tree_topology_filter()`, `fast_unifrac()`
  (`gemelli/preprocessing.py`, `gemelli/utils.py`); scikit-bio 0.7.4
  `vectorize_counts_and_tree()` and `TreeNode.assign_ids/prune`; iow `bp`
  `shear` / `collapse`. All BSD-3. Add scikit-bio's and iow's notices
  beside gemelli's in `inst/COPYRIGHTS` if code is ported closely.
  (gemelli's `tree_topology_filter` came from biocore/wol; check its
  license, or rewrite it from the behaviour above, which is trivial.)
- **Fixtures:** `data-raw/expected/rpca_phylo.py` runs gemelli on gemelli's
  test tree (`tree.nwk`) and a random table; **generate fresh outputs**
  (gemelli's own QIIME 2 phylo fixtures are overwritten by several tests,
  and `labeled-phylogeny.nwk` is missing upstream). Compare scores,
  loadings, eigenvalues and distances at about 1e-10 (sign-fixed with
  `ord_flip()`), as `rpca()`'s tests do.


## Tests

- Fixtures; `min_depth` 0 and > 0; `node_counts` matches gemelli's
  counts-by-node table.
- With a star tree (every tip on the root) of unit lengths, the node
  table is the tips plus one all-missing root column; check against
  gemelli.
- Same results through `count_ord(x, 'rpca_phylo', tree = t)`.
- Identical for `cpus` 1, 2 and 4.
- Plain `rpca()` results unchanged by the refactor (compare to stored
  fixtures before and after).
- Errors as `rpca()`, plus tree errors.
- A smoke-test line.


## Docs

- Rd page under Ordination (counts), cross-linked with `rpca()`,
  `philr()` and the UniFrac family; "Differences from gemelli" (no int
  truncation; root-length edge case; ignored `max_iterations`).
- NEWS.


## References

- Martino, C., Morton, J. T., Marotz, C. A., et al. (2019). A novel sparse
  compositional technique reveals microbial perturbations. *mSystems*, 4,
  e00016-19.
- Martino, C., McDonald, D., Cantrell, K., et al. (2022). Compositionally
  aware phylogenetic beta-diversity measures better resolve microbiomes
  associated with phenotype. *mSystems*, 7, e00050-22.
