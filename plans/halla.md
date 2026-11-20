# `halla()`: hierarchical all-against-all association testing (HAllA)

**Priority:** 84. **Tier:** port (HAllA 0.8.40, Python, **MIT**: the
repository's `LICENSE` "The HAllA software is licensed under the MIT
license. Copyright (c) 2019 Harvard School of Public Health", and PyPI's
metadata for 0.8.40, `license = MIT`; checked 2026-11-20. The GitHub `master`
`setup.py` still says 0.8.30: read the 0.8.40 sdist from PyPI when
implementing, and record any difference). SciPy's average-linkage code
(`scipy/cluster/_hierarchy.pyx`, **BSD-3**) is ported for the trees. HAllA's
optional helpers through rpy2, the R packages `eva` (GPD tails) and `XICOR`,
are not ported (see "Scope"). Read [README.md](README.md) first, then
[cross_assoc.md](cross_assoc.md), whose pairwise code this item reuses.


## What it is

HAllA (Ghazi et al. 2022) finds **blocks of associated features** between
two paired tables (taxa × metabolites, taxa × host genes, taxa × clinical
variables). Testing all p × q pairs and reporting each significant one
(what `cross_assoc()` does; HAllA calls it "AllA") gives long lists of
redundant hits, because features within each table are correlated. HAllA:

1. computes all pairwise associations and their BH q-values (AllA);
2. clusters each table's features hierarchically by their own similarity;
3. descends the two trees together from the roots, reporting a block
   (a cluster of X × a cluster of Y) as soon as at least `1 − fnr` of its
   cells are FDR-significant, otherwise splitting the side whose split
   gives purer halves (weighted Gini impurity of the significant/not
   indicator);
4. ranks the blocks by their best q-value.

Its signature figure is the **hallagram**: the association heatmap with
features in tree order, significant cells dotted and the top blocks boxed
and numbered. Uptake: Ghazi et al. (2022) has **88** OpenAlex citations
(5 in 2022 rising to 28 in 2026), and "HAllA" has **112** Europe PMC
microbiome full-text mentions 2022–26 (retrieved 2026-11-20). Modest but
rising; it is the bioBakery's multi-omics association tool, used alongside
MaAsLin.


## Names

- `halla()`: the method's own name, lowercase (rule 1), as users know
  `halla` (the command) and `HAllA` (the Python class). No registry alias
  starts with `halla` (check with README's snippet).
- AllA is `cross_assoc()`; the docs say so, and `halla()`'s result
  contains the same pairwise table.
- **Same names in other packages:** none in R.


## Design

### Signature

```r
halla(x, y, data = NULL, method = 'spearman', margin = 1L,
      norm = 'percent', norm_y = 'none', pseudocount = NULL,
      alpha = 0.05, p_adjust = 'BH', fnr = 0.2, rank_by = 'best',
      gini_tolerance = 0.02, linkage = 'average', max_freq = 1,
      cpus = n_cpus())
```

Inputs, `method`, `norm`, `norm_y`, `pseudocount`, `alpha`, `p_adjust` mean
exactly what they mean in `cross_assoc()` (same validation code). HAllA's
names and defaults (`config.yaml` of 0.8.40, as read from GitHub master;
re-check against the sdist):

| ecodive | HAllA | Default | Meaning |
|:--------|:------|:--------|:--------|
| `method` | `pdist_metric` | `'spearman'` | association; `'pearson'` too (Scope) |
| `alpha` | `fdr_alpha` | 0.05 | FDR level for a significant cell |
| `p_adjust` | `fdr_method` | `'BH'` (`fdr_bh`) | p-value adjustment over all pairs |
| `fnr` | `fnr_thresh` | 0.2 | a block is reported when ≥ 1 − `fnr` of its cells are significant |
| `rank_by` | `rank_cluster` | `'best'` | blocks ranked by their smallest (`'best'`) or mean (`'average'`) q-value |
| `gini_tolerance` | `gini_uncertainty_level` | 0.02 | Gini scores closer than this count as tied (see step 4) |
| `linkage` | `linkage_method` | `'average'` | hierarchical linkage (SciPy's) |
| `max_freq` | `max_freq_thresh` | 1 | drop a feature whose most frequent value covers ≥ this fraction of samples (1 = drop only constant features) |

`fnr`, `rank_by`, `gini_tolerance`, `linkage`, `max_freq` are new names
(add them to parameter_names.md; **ask the user**). `max_freq` uses `>=`
as HAllA's `keep_feature()` keeps `count.max()/count.sum() < thresh`.
`alpha` here is the cell-level FDR, not a block-level error rate; the docs
say what HAllA does and does not control.

No `seed`: with the closed-form p-values (Spearman, Pearson) HAllA is
deterministic.

### Algorithm (HAllA 0.8.40, `halla/main.py`, `utils/tree.py`, `utils/stats.py`, `utils/similarity.py`, `utils/data.py`, `hierarchy.py`)

1. **Preprocess** (`preprocess()`): drop features with `max_freq` (above),
   per table; no transform (HAllA's `transform_funcs` default is none);
   discretisation is bypassed for continuous data with Spearman/Pearson.
   HAllA joins the two tables on their shared sample columns; ecodive
   matches by name as `cross_assoc()` does.
2. **Pairwise table** (AllA, `_compute_pairwise_similarities()`): for each
   X feature and Y feature, `scipy.stats.spearmanr` (or `pearsonr`) on the
   pairwise non-missing samples (`remove_missing_values()`); a pair where
   either side is constant after removing missing values gets similarity 0
   and p = 1 (HAllA's convention; `cross_assoc()` gives `NA` there, so
   `halla()` maps `NA` to 0 and 1 and documents it). p-values are SciPy's
   t-approximation (`does_return_pval('spearman')` is `True`, so no
   permutations). q-values: `statsmodels.multipletests(method='fdr_bh')`
   over all p × q pairs, and its `reject` flags. ecodive uses
   `p.adjust('BH')` and `reject = p_adjusted <= alpha`; statsmodels' BH
   reject compares the sorted p-values with `alpha * k / m`, which should
   agree, including a q exactly equal to `alpha`: **confirm from the
   fixtures**.
   **Reuse `cross_assoc()`'s C kernels**, called internally.
3. **Trees** (`HierarchicalTree`): within each table, the similarity of
   every pair of features with the same function (`spd.pdist` with the
   Spearman function), turned into a distance \eqn{1 - |r|}
   (`similarity2distance(set_abs = True)`), clipped at 0, then
   `scipy.cluster.hierarchy.linkage(method = 'average')` and `to_tree()`.
   - The within-table similarities reuse `cross_assoc()`'s kernel with
     `y = x` (packed triangle; ranks once).
   - **Port SciPy's linkage** (average linkage through the nearest-neighbour
     chain, `nn_chain` in `_hierarchy.pyx`, BSD-3, ~150 lines) rather than
     calling `stats::hclust()`: with tied distances (common: Spearman on
     zero-inflated taxa) R's and SciPy's merge order can differ, and the
     blocks depend on the tree. Return the trees in `hclust` form
     (`merge`, `height`, `order`, `labels`) so `plot(as.dendrogram())`
     works; `order` is SciPy's `to_tree().pre_order()` (left before right,
     left being the first cluster id of the linkage row).
   - Threaded: the distance triangle as in `cross_assoc()`; the linkage
     itself is O(m²) sequential per table, the two tables in two tasks.
4. **Blocks** (`compare_and_find_dense_block()`), recursion from the two
   roots (X node, Y node):
   - Block = leaves(X node) × leaves(Y node), in pre-order. If the
     fraction of rejected cells ≥ 1 − `fnr` (`is_densely_associated()`,
     `>=`), report it and stop.
   - Else bifurcate one side (`bifurcate_one()`): if both are leaves, stop;
     if one is a leaf, split the other; else compute the weighted Gini
     impurity of splitting X (its two children × all of Y) and of splitting
     Y, \eqn{\sum_b (|b|/N)(1 - p_b^2 - (1-p_b)^2)}. If `gini_tolerance >
     0`, the two differ by less than it, **and** \eqn{|X| \cdot |Y| > 25},
     split X if it has strictly more leaves than Y, else Y (HAllA's
     `x_tree.get_count() > y_tree.get_count()`); otherwise split X if its
     Gini is strictly smaller, else Y (so Y on equality).
   - Recurse over `itertools.product(x_branches, y_branches)` in order
     (left before right, X branches outer). Blocks are appended in
     discovery order.
   - Then sort by the best (or mean) q-value of the block's cells.
     Python's sort is stable: ties keep discovery order. Use R's
     `order(..., method = 'radix')` (stable) on the same key.
   - Implement in C or R? The recursion touches each cell O(depth) times;
     with prefix sums of the reject indicator in tree order (a 2-D
     cumulative sum over the pre-ordered p × q 0/1 matrix), any block's
     count is O(1). Write it in R first with the prefix sums; move to C
     only if profiling says so.
5. **Not ported:** `trim_block()` (unused by `HAllA.run()`), the
   diagnostic lattice plots, file output, `splitting_diagnostic_mode`.

### Scope (what to leave out, and why)

- **Pearson** and **Spearman** only, at first. HAllA's other metrics:
  - `nmi` / `mi` (with quantile discretisation through scikit-learn's
    `KBinsDiscretizer`, and chi-squared p-values from
    `scipy.stats.chi2_contingency`, which applies Yates' correction on
    2 × 2 tables): HAllA switches to MI automatically when a table has
    categorical features. A later phase if the user wants mixed
    categorical metadata; both references are BSD-3, so portable.
  - `dcor` (distance correlation) has no closed-form p-value: HAllA runs
    permutations with early stopping and a generalised-Pareto tail fitted
    by R's `eva` package through rpy2 (eva 0.2.7 is **GPL (>= 2)**, CRAN
    `DESCRIPTION` checked 2026-11-20: output-only). Its early
    stopping and `np.random.seed(seed)` per pair are not reproducible
    across threads. If wanted later: plain permutations on PCG, no GPD.
  - `xicor` calls the R package XICOR (0.4.1, **Apache (>= 2)**, CRAN
    `DESCRIPTION` checked 2026-11-20: portable, keeping the Apache notice
    as `src/pcg_basic.c` does). A later phase if asked; it has a
    closed-form asymptotic p-value (Chatterjee 2021).
- **The `HAllA(... seed=0)` argument** matters only for the permutation
  metrics.

### Result

A list of class `ecodive_halla`:

- `associations`: the `cross_assoc()` table (all pairs, with `q_value` as
  `p_adjusted` and `reject`), so `as.matrix()` and the AllA heatmap work on
  it;
- `blocks`: a data frame, one row per reported block in rank order:
  `block` (rank, 1..), `x_features`, `y_features` (list columns of names,
  in tree order), `n_x`, `n_y`, `best_q` (or `mean_q`, per `rank_by`),
  `fraction_significant`, `mean_association` (mean signed r over the
  block's cells, for reading the sign);
- `x_tree`, `y_tree`: `hclust` objects of the two tables' features;
- attributes `method`, `alpha`, `p_adjust`, `fnr`, `rank_by`,
  `n_samples` (rule 1: they qualify the blocks).

`print()`: counts (features kept per table, significant pairs, blocks) and
the first blocks with their members abbreviated, as HAllA's
`sig_clusters.txt` (`cluster_rank`, `cluster_X`, `cluster_Y`,
`best_adjusted_pvalue`).


## Reference implementations

- **HAllA 0.8.40** (PyPI, MIT): `halla.HAllA(...).load(X, Y).run()`;
  outputs `all_associations.txt` (similarity, p, q per pair),
  `sig_clusters.txt` (blocks), `X_linkage.npy` / `Y_linkage.npy`. Its
  `stats.py` imports `rpy2` and the R package `eva` at module load even
  for Spearman: the fixture environment needs R + `eva` installed in the
  sandbox (a scratch library), not in the user's library.
- **SciPy** `cluster.hierarchy.linkage(method='average')`: tree fixtures,
  including tied distances.
- Fixtures: `data-raw/expected/halla.py` writes text outputs under
  `data-raw/expected/halla/`, folded into `tests/tinytests/expected/halla.rds`
  by `data-raw/expected/halla.r`. Data: HAllA's own synthetic generator
  (`scripts/synthetic_data.py`, MIT: blocks of associated features with
  known truth) at small sizes (e.g. 50 samples, 30 × 30 features, 3
  blocks), plus a real-data-like case with zero-inflated taxa (ties) and a
  metadata table with `NA`.
- Match: similarities and p-values at 1e-12 / 1e-10; q-values and rejects
  exactly; linkage matrices (merge order and heights at 1e-12); blocks
  (member sets and rank order) exactly.


## Plotting

`plot.ecodive_halla(x, n_blocks = 50, mask = FALSE, dots = TRUE,
col = NULL, key = TRUE, ...)`: the **hallagram** (HAllA's
`generate_hallagram()`), in base graphics:

- Rows: the X features that appear in the top `n_blocks` blocks (HAllA's
  `block_num = 50`, `trim = True`), in X-tree pre-order; columns likewise
  for Y.
- Cells coloured by the association on the diverging `"Blue-Red 3"` scale
  symmetric about 0 (HAllA: `RdBu_r` for signed metrics, `vmin = -vmax`);
  with `mask = TRUE`, cells outside the shown blocks are blank.
- `dots`: a white dot with a black ring on each FDR-significant cell
  (HAllA's `signif_dots`).
- Each top block outlined in black and numbered by rank at its centre
  (white text with a dark outline, size scaled to the block); blocks
  ranked below `n_blocks` whose features are all shown are outlined in
  grey (HAllA's `show_lower`). The block outlines need the shown features
  of each block to be contiguous, which tree order guarantees for whole
  blocks; a block only partly shown gets a dotted edge, as HAllA draws.
- Title "HAllA: Spearman Correlation" (from `method`'s display name in
  `cross_assoc()`'s internal table); the colour key from `dist_plot_key()`
  on the left or right. X labels on the right, Y labels at the bottom, as
  in HAllA's figure.
- Reads only the result (`associations`, `blocks`, trees). Returns `x`
  invisibly.

A "clustermap" variant with the two dendrograms in the margins (HAllA's
`generate_clustermap()`) is optional: `plot(x, dendrograms = TRUE)` if
cheap.


## Tests

- Fixtures above (similarities, p, q, trees, blocks, ranks).
- The linkage port equals `stats::hclust(method = 'average')` merge
  heights (as multisets) on tie-free distances, and SciPy's merge order on
  tied ones.
- Block search: hand-made reject matrices with known answers (a full
  block found at the root; a checkerboard that splits to single cells; the
  `gini_tolerance` size rule with |X||Y| ≤ 25 and > 25); prefix-sum counts
  equal brute-force counts.
- `identical()` for `cpus` 1, 2, 4.
- Edge cases: no significant pair (no blocks; plot says so, no error); a
  single feature in one table; a constant feature (dropped by `max_freq`);
  `fnr = 0`; mismatched sample names.
- `plot()` keeps `.Random.seed` and `par()` unchanged; `n_blocks = 1`.
- One smoke-test line (20 samples, 5 × 5 features).


## Docs

- Rd page in "Associations between tables", next to `cross_assoc()`.
- Explain AllA versus HAllA, what `fnr` means (blocks may contain up to
  `fnr` non-significant cells), that the block-level error rate is not
  formally controlled (Ghazi et al. evaluate it by simulation), and that
  features with mostly zeros make Spearman ties.
- "Differences from HAllA": no MI/NMI, dcor or XICOR yet; constant pairs
  are 0 / p = 1 as in HAllA, but reported `NA` in `cross_assoc()`; sample
  matching by name; no output files; the hallagram in base graphics.
- NEWS; the new argument names in parameter_names.md; `inst/COPYRIGHTS`
  stanzas for HAllA (MIT) and SciPy (BSD-3), and `Authors@R` credit for
  HAllA's authors if not already there.


## References

- Ghazi, A. R., Sucipto, K., Rahnavard, A., Franzosa, E. A., McIver, L. J.,
  Lloyd-Price, J., Schwager, E., Weingart, G., Moon, Y. S., Morgan, X. C.,
  Waldron, L., & Huttenhower, C. (2022). High-sensitivity pattern
  discovery in large, paired multiomic datasets. *Bioinformatics*,
  38(Suppl. 1), i378–i385. doi:10.1093/bioinformatics/btac232
- Müllner, D. (2011). Modern hierarchical, agglomerative clustering
  algorithms. arXiv:1109.2378.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
  rate. *JRSS B*, 57, 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Breiman, L., Friedman, J., Olshen, R., & Stone, C. (1984).
  *Classification and Regression Trees*. Wadsworth (Gini impurity).
