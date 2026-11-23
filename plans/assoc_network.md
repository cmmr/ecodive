# `assoc_network()`, `assoc_network_compare()`, `assoc_network_diff()`: association networks

**Priority:** 33. **Tier:** clean room, from textbook formulas and the NetCoMi
paper (Peschel et al. 2021, *Brief Bioinform*). NetCoMi is GPL-3 and
output-only: its paper, its Rd documentation (quoted below as "NetCoMi docs")
and its outputs may be used, its R source may not. igraph (GPL ≥ 2,
checked on CRAN on 2026-11-20) is likewise output-only. Read [README.md](README.md), then
[compositional.md](compositional.md), then [sparcc.md](sparcc.md), which
defines the input structure.


## What it is

NetCoMi made the analysis layer *on top of* an association matrix
standard. That layer covers: sparsify into edges; centralities and hubs;
clusters; global properties; comparison of two groups' networks (Jaccard
index of central nodes, adjusted Rand index of clusterings, permutation
tests); and differential associations (Fisher's z). NetCoMi has 672
citations, rising every year from 76 in 2022 to 170 in 2026 (OpenAlex,
2026-11-19). This plan gives ecodive that layer for every association result
it returns: `rho_p()` now; `sparcc()`, `spiec_easi_*()` and `secom_*()`
later.


## Names

- The `assoc_network*` family (rule 2), not `network()`, which would mask
  the network package (compositional.md). `_compare` and `_diff` are
  companions that take networks or associations. They are plain functions:
  they return lists or tables, so no wrapper (README).
- No registry has an alias starting with `assoc` or `net` (checked
  2026-11-20).

### An association wrapper: now justified (recommendation)

compositional.md left a wrapper for `vlr()` / `rho_p()` / `phi_s()` to
"revisit if more association measures arrive". Five are now planned
(`sparcc`, `spiec_easi_mb`, `spiec_easi_glasso`, `secom_pearson`,
`secom_distance`), making eight that share the first argument `counts` and
one output structure. The wrapper is also *needed*: `assoc_network_compare()`'s
group-permutation test must re-run "whatever association made these networks"
on permuted samples, and a registry is the clean way to call a method by
name with the caller's settings. **Recommend** adding it, as part of this item
or just before it:

- `feature_assoc(counts, method, ...)` with registry `FEATURE_ASSOCS` and
  `list_feature_assocs()` / `match_feature_assoc()`, shaped like `count_ord()`.
  **`FEATURE_ASSOCS` itself landed with item 32** (`R/feature_assoc.r`, rows
  for `vlr`, `rho_p`, `phi_s`, `sparcc`; the range as numeric `lower` and
  `upper` columns); see sparcc.md's "What item 33 can rely on".
- Registry columns: `id`, `name`, `alt_ids`, `fun`, plus `kind`
  (`'dissimilarity'` for vlr and phi_s; `'similarity'` for the rest), `range`
  (`'[-1,1]'`, `'[0,1]'`, `'[0,Inf)'`), `fisher_z` (TRUE only for
  Pearson-type correlations: `sparcc`, `secom_pearson`), `sparse` (TRUE
  for `spiec_easi_*`: non-zero entries are the selected graph) and
  `strict` (TRUE for `secom_pearson`: its threshold keeps |r| > t).
- A row in the registry means "over features": `dist_clust.r` and
  `plot_dist.r` look up `method` here instead of testing names. These
  facts never sit on the result as flags (README.md § "Result attributes
  and plots").
- **Not `assoc()`:** that masks `vcd::assoc()` (mosaic association plots),
  which is widely used and means something else. `feature_assoc` follows the
  `<input>_<kind>` pattern of `dist_ord` / `count_ord` and says what is
  associated.
- Add its row to README.md's shadowing table: members must not take
  `method` or any prefix of `method` / `counts`. Every planned signature
  complies; `test-wrapper_args.r` enforces it.
- Ask the user about the name before implementing.


## Design

### Input: association results (from sparcc.md)

An `ecodive_dist` over features whose `method` has a `FEATURE_ASSOCS` row
of kind `similarity`, with a top-level `n_samples` and optionally a
parallel `p_value` attribute. `assoc_network()` accepts results of kind
`similarity` with values in [−1, 1] or [0, 1]. For `vlr()` / `phi_s()`
(unbounded dissimilarities) it errors: "use rho_p(), or convert with
1/(1 + x)". A plain symmetric matrix with unit diagonal is also accepted
(converted with `as.dist()`), so users can bring their own correlations;
it then has no `method` or `n_samples`, so `fisher_z` and `n_samples` must
be supplied.

### `assoc_network()`

```r
assoc_network(x, threshold = NULL, alpha = NULL, p_adjust = 'BH',
              sign = 'signed', clusters = 'greedy',
              hubs = 'eigenvector', hub_quantile = 0.95,
              weighted_degree = FALSE)
```

1. **Edges.** \eqn{r^*_{ij} = r_{ij}} if it passes every filter given, else 0.
   **`assoc_pairs()` (landed 2026-11-20, `R/assoc_pairs.r`) already does
   this step**, with the same `threshold` (inclusive), `alpha` and
   `p_adjust` (adjusted over all pairs, then filtered), for `sparcc()`,
   `rho_p()` and correlation matrices; build `edges` from its `from`, `to`,
   `association`, `p_value` and `p_adjusted` rather than filtering again.
   It returns the pairs strongest first; sort them as `edges` needs. Its
   only gap is `threshold = NULL` (a result's own sparsity, for
   SPIEC-EASI and SECOM), which this item adds:
   - `threshold`: \eqn{|r_{ij}| \ge} threshold (NetCoMi docs: "greater than
     or equal"). `NULL` = the result's own sparsity, where it has one
     (SPIEC-EASI's selected graph: non-zero entries; SECOM: its
     cross-validated threshold), else 0.3 (NetCoMi's default `thresh`, and
     the SparCC paper's network cutoff).
   - `alpha`: needs a `p_value` attribute. The edge is kept if the p-value
     adjusted by `p_adjust` over all pairs is ≤ alpha. This is NetCoMi's
     testing route (its default adjustment is "adaptBH"/"lfdr", which the docs
     contradict). ecodive uses `p.adjust()` methods only; say so.
2. **Dissimilarity and similarity** (Peschel et al. 2021, §2.2, after van
   Dongen & Enright): `sign = 'signed'` gives \eqn{d = \sqrt{0.5(1-r^*)}};
   `'unsigned'` gives \eqn{d = \sqrt{1-r^{*2}}}; `'positive'` (NetCoMi's
   "signedPos") is signed with negative edges dropped. Then \eqn{s = 1 - d}.
   Edges keep their sign for plotting. These are NetCoMi's names except
   `positive`; it is an option in every reference, so it stays an argument.
3. **Centralities.** Use NetCoMi's definitions (docs; defaults `centrLCC =
   TRUE`, all normalised), on the **largest connected component**; nodes
   outside it get 0. A tie for the largest component goes to the one with the
   lowest-index node.
   - Degree: the number of neighbours (or \eqn{\sum s} with `weighted_degree`),
     divided by n − 1.
   - Betweenness: Brandes (2001), with Dijkstra on edge lengths \eqn{d},
     times \eqn{2/(n^2-3n+2)}. Path lengths are compared with a relative
     tolerance of 1e-12 when counting equal shortest paths; document it.
   - Closeness: **harmonic**, \eqn{\sum_{u} 1/\mathrm{dist}(v,u)} over the
     component, divided by n − 1 (NetCoMi docs; the paper's table reads
     differently).
   - Eigenvector: the leading eigenvector of the similarity matrix \eqn{S}
     (non-negative, so Perron), scaled to max 1. Use power iteration with a
     fixed start of all ones, or `src/lanczos.c`. Test against `eigen()`.
   - n is the size of the component (verify against NetCoMi's outputs; the
     docs don't say).
4. **Hubs.** `hubs` names one or more of `'degree'`, `'betweenness'`,
   `'closeness'`, `'eigenvector'`. A node is a hub when it is above the
   empirical `hub_quantile` quantile (R type 7; **unverified** that NetCoMi
   uses type 7) of **every** named centrality (NetCoMi's `hubPar`,
   `hubQuant = 0.95`). A log-normal-fit quantile (`lnormFit`) is not planned
   now.
5. **Clusters.** `'greedy'` is Clauset–Newman–Moore fast greedy modularity
   (NetCoMi's default for association networks) on weights \eqn{s}, written
   from the 2004 paper. Ties in the merge gain go to the pair with the
   smallest (i, j). igraph may break ties differently: compare
   **modularity**, and partitions only where there are no ties. `'none'`
   puts every node in cluster 0. A singleton node gets its own cluster.
6. **Global properties** (NetCoMi Table 4): `n_nodes`, `n_edges`, `density`,
   `positive_fraction`, `n_components`, `lcc_size`, `average_path_length`
   (LCC; normalised by the mean edge dissimilarity, as NetCoMi's
   `sPathNorm = TRUE`), `clustering_coefficient` (global transitivity,
   unweighted, which differs from NetCoMi's weighted Barrat version: say
   so), and `modularity` (Newman 2006, weighted, of the clusters).

**Result**: a list of class `ecodive_network`:
- `edges`: a data.frame `from`, `to` (feature names, from before to in
  feature order), `association` (\eqn{r}), `similarity`, `distance`, `sign`
  (±1), `p_value`/`p_adjusted` when available;
- `nodes`: a data.frame `feature`, `degree`, `betweenness`, `closeness`,
  `eigenvector`, `cluster` (integer), `hub` (logical), `in_lcc`;
- `global`: a named numeric vector;
- attributes `method` (the input's association method, e.g. `"sparcc"`),
  `n_samples`, `threshold` (the cutoff used, including one the result
  chose itself) and `display`; see "Plotting" below.

Plus a `print` method (counts, hubs) and `as.data.frame()` (the `nodes`
table).

Cost: all-pairs Dijkstra is O(V(E + V log V)) on the LCC. Run it in C,
threaded over source nodes. Each source's dependency vector is summed into
the betweenness in source order after all threads finish (no shared
accumulator), so results are identical for any `cpus`.

### `assoc_network_compare()`

```r
assoc_network_compare(x, y, quantile = 0.75, permutations = 1000L,
                      counts = NULL, groups = NULL, data = NULL,
                      assoc_args = list(), network_args = list(),
                      p_adjust = 'BH', seed = 0L, cpus = n_cpus())
```

`x`, `y` are `ecodive_network`s over the same features (matched by name;
error if the sets differ).
- **Differences:** each node's centralities and each global property, x − y.
- **Jaccard index** of the "most central" sets (above the empirical
  `quantile` of each centrality; NetCoMi's `jaccQuant = 0.75`) and of the
  hub sets. P(J ≤ j) and P(J ≥ j) are exact hypergeometric probabilities
  conditional on the two set sizes (k = |A ∩ B| ~ Hyper(N, |A|, |B|), and J
  is monotone in k). NetCoMi follows Real & Vargas (1996), whose model may
  not be conditional: **compare with NetCoMi's outputs** and document any
  difference.
- **Adjusted Rand index** (Hubert & Arabie 1985) of the two clusterings
  (singletons and nodes outside the LCC kept as their own clusters; check
  against NetCoMi's outputs). Test H0: ARI = 0 by permuting one network's
  cluster labels over nodes, `permutations` times (PCG): `p = (1 + #{ARI_b ≥
  ARI})/(B + 1)`, one-sided as NetCoMi docs say.
- **Group-permutation tests** (only when `counts` and `groups` are given):
  - `groups` must have two levels; x was built from the first level, y from
    the second.
  - Each permutation shuffles the group labels (keeping group sizes; reuse
    `perm_matrix()` / `src/permute.c`), re-estimates **both** associations
    through `feature_assoc()` with the networks' `method` (which must
    match) and the settings in `assoc_args`, and rebuilds both networks
    with `assoc_network()` and the settings in `network_args`. Results
    don't record the settings they were made with, so the caller passes
    the same ones again.
  - p-values: `perm_p(|diff_b|, |diff|)`, adjusted over nodes by `p_adjust`.
    NetCoMi defaults `permTest = FALSE`, `nPerm = 1000`.
  - Each permutation is a task with its own seed stream. A SparCC stream must
    not depend on which thread runs it, so results are identical for any
    `cpus`.
  - Document the cost: 2 × B association runs.
- **Result:** a list of class `ecodive_network_compare` with
  - `nodes` (feature, each centrality for x and y, difference, and
    `p_value` / `p_adjusted` when tested);
  - `global` (property, x, y, difference, p-value);
  - `jaccard` (set, jaccard, p_less, p_greater);
  - `ari` (value, p_value);
  - `edges_x`, `edges_y` (the two networks' `edges`, which its plot
    draws).

### `assoc_network_diff()`: differential associations (Fisher's z)

```r
assoc_network_diff(x, y, alpha = 0.05, p_adjust = 'BH')
```

`x`, `y` are two **association results** (not networks) from the same
method, with `fisher_z = TRUE` in its registry row.
- \eqn{z = (\mathrm{atanh}\,r_x - \mathrm{atanh}\,r_y)/\sqrt{1/(n_x-3) +
  1/(n_y-3)}}, with n from each result's `n_samples`, or from a parallel
  `n_pairs` attribute where a method uses pairwise-complete samples (SECOM).
  The p-value is two-sided normal; values with \eqn{|r| \ge 1} give `NA` and
  a warning.
- For SparCC this is an approximation, since its correlations aren't sample
  Pearson correlations: say so. NetCoMi likewise restricts Fisher's z to
  correlations.
- Refuse `rho_p()` (concordance, not Pearson) and `spiec_easi_*()` (partial
  correlations need n − 3 − (p − 2) and are a different test).
- **Result:** a data.frame of class `ecodive_network_diff`, one row per pair:
  `feature1`, `feature2`, `association_x`, `association_y`, `z`,
  `p_value`, `p_adjusted`, `reject`. This matches the per-feature result
  convention, with `feature1`/`feature2` in place of `feature`.
- NetCoMi's non-parametric resampling test and Discordant are not planned.

### `plot()` methods (base graphics)

- **Layout: circular** (chosen because it is deterministic, needs no
  dependency, and is easy to share between two networks). Nodes go
  clockwise from 12 o'clock, grouped by cluster (larger clusters first), and
  within a cluster by decreasing degree. Isolated nodes come last, or are
  left out with `isolated = FALSE`. Edges are straight segments; an option
  draws them as quadratic Bézier curves pulled toward the centre, via
  `graphics::xspline()`.
- `layout` also accepts an n × 2 matrix with row names (features), so
  users can supply any layout. A Fruchterman–Reingold layout would need
  iterations and a seed; not now.
- Styling, as NetCoMi:
  - node colour by cluster (`plot_palette()`);
  - node size by the first `hubs` centrality;
  - hubs drawn with a border and bold labels;
  - edge colour by sign, using the two ends of `hcl.colors(2, "Blue-Red 3")`
    as `plot.ecodive_dist()` does (colour-blind-safer than NetCoMi's
    green/red);
  - edge width ∝ |r|.
- `plot.ecodive_network_compare()`: the two networks side by side on **one
  shared layout** (built from the union of nodes, ordered by x's clusters),
  as NetCoMi does.
- `plot.ecodive_network_diff()`: the differential network, edges with
  `reject`, coloured by which group's association is larger.
- All return their input invisibly.


## Reference implementations (output-only)

- **NetCoMi** (GPL-3, `stefpeschel/NetCoMi`). Fixtures come from
  `netConstruct(..., measure = "pearson", sparsMethod = "threshold",
  thresh = 0.3, dissFunc = "signed")` on a fixed association matrix passed as
  `dataType = "correlation"`, then `netAnalyze()` and `netCompare()`. Use a
  supplied matrix, so the fixtures test only this layer. Save edges,
  centralities, hubs, clusters, global properties, Jaccard probabilities and
  ARI.
- Deterministic parts (edges, degree, betweenness, closeness, eigenvector,
  modularity, ARI, Fisher z): match at 1e-10.
- Clusters: match exactly when the fixture has no ties in the greedy merge;
  otherwise compare modularity.
- Permutation p-values: match in distribution only (NetCoMi uses R's RNG).
  Supply a fixed permutation matrix through an internal hook for exact
  checks.
- igraph outputs (`betweenness()`, `closeness()`, `cluster_fast_greedy()`)
  serve as a second check of the textbook formulas.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Title: `assoc_network()` sets `display$title` with `set_display()` and
  `title_of("Association Network", <the input's display$title>)`, e.g.
  "Association Network of SparCC Correlation". A plain matrix input has
  no title, so the title is "Association Network". The subtitle is the
  top-level `threshold` and `global['n_edges']`.
- `assoc_network_diff()` titles its result the same way from x's title
  ("Differential Associations of SparCC Correlation");
  `assoc_network_compare()` from x's network title ("Comparison of
  Association Network of SparCC Correlation"), with panels labelled by
  the levels of `groups`, or "x" and "y".
- What the plots read: the edge `sign` (colour) and `association` (width)
  columns, node `cluster`, `degree` and `hub` (order, size, labels), and
  `threshold` (legend). The circular layout is built from those columns,
  so no layout is stored; the compare plot builds one from x's clusters
  and draws `edges_x` and `edges_y` on it.

## Tests

- Hand-built graphs with known answers: a star, a path, a ring, two
  cliques joined by one edge, a disconnected graph with isolated nodes.
- All four centralities against base-R reference implementations (Floyd–Warshall
  for paths) at 1e-12. Fisher z against a direct formula. Hypergeometric
  Jaccard against `phyper()`.
- NetCoMi fixtures as above.
- Edge cases: no edges (every centrality 0, modularity 0, plot draws nodes
  only); every pair an edge; |r| = 1; a single feature; `alpha` given with
  no `p_value` (error); features that differ between x and y (error).
- `identical()` for `cpus` 1, 2 and 4, for `assoc_network()` and the
  permutation tests.
- Plot methods run on the smoke data without error (`pdf(NULL)`).
- Smoke-test lines for all three functions.


## Docs

- One Rd page for `assoc_network()` (with `plot.ecodive_network`), one for
  `_compare` and one for `_diff`.
- Each says what NetCoMi calls the same thing, and lists the "Differences
  from NetCoMi": the BH family only, no log-normal hub fit, unweighted
  transitivity, the hypergeometric Jaccard test, the circular layout, and
  the colours.
- `feature_assoc()` gets its own page if accepted.
- pkgdown "Compositional Data" → *Associations and networks*. NEWS.
- Add `threshold`, `quantile`, `hub_quantile`, `sign`, `clusters`, `hubs`,
  `weighted_degree`, `assoc_args` and `network_args` to
  parameter_names.md.


## References

- Peschel, S., Müller, C. L., von Mutius, E., Boulesteix, A.-L., & Depner, M.
  (2021). NetCoMi: network construction and comparison for microbiome data
  in R. *Briefings in Bioinformatics*, 22, bbaa290. doi:10.1093/bib/bbaa290
- Brandes, U. (2001). A faster algorithm for betweenness centrality. *J Math
  Sociol*, 25, 163–177. doi:10.1080/0022250X.2001.9990249
- Clauset, A., Newman, M. E. J., & Moore, C. (2004). Finding community
  structure in very large networks. *Phys Rev E*, 70, 066111.
  doi:10.1103/PhysRevE.70.066111
- Newman, M. E. J. (2006). Modularity and community structure in networks.
  *PNAS*, 103, 8577–8582. doi:10.1073/pnas.0601602103
- Hubert, L., & Arabie, P. (1985). Comparing partitions. *J Classif*, 2,
  193–218. doi:10.1007/BF01908075
- Real, R., & Vargas, J. M. (1996). The probabilistic basis of Jaccard's
  index of similarity. *Syst Biol*, 45, 380–385. doi:10.1093/sysbio/45.3.380
- van Dongen, S., & Enright, A. J. (2012). Metric distances derived from
  cosine similarity and Pearson and Spearman correlations. arXiv:1208.3145
  (cited by NetCoMi for the signed/unsigned distances; not checked here).


## Decisions made while implementing (item 33, 2026-11-20)

Asked of the user at the start:

1. **`feature_assoc(counts, method, ...)`** is added now, with
   `list_feature_assocs()` / `match_feature_assoc()`, shaped like
   `count_ord()` (`R/feature_assoc.r`). No `fun` column: each `id` is its
   function's name. Aliases: `logratiovariance` (vlr), `rho rhop propr`
   (rho_p), `phi phis phs` (phi_s), `fastspar` (sparcc); `'proportionality'`
   matches both rho_p and phi_s and so is an error. Checked: no other
   registry has these as exact names, nor any name starting `assoc` or
   `net` (tested). Its row is in README.md's shadowing table, and
   `test-wrapper_args.r` checks it.
2. **Clustering coefficient: NetCoMi's** (the mean of Barrat's weighted
   local coefficient over nodes with two or more neighbours), not the
   plan's unweighted transitivity.
3. **Jaccard p-values: NetCoMi's model** (see 9), not the plan's
   conditional hypergeometric.
4. **One `permutations`** (default 1000) for the ARI test and the
   group-permutation test; the docs give the cost.

Made while implementing. Every NetCoMi convention below was read from its
**outputs** (and its Rd); no NetCoMi or igraph source was opened:

5. **Reference: NetCoMi 1.2.0**, from the jsDelivr mirror of
   `stefpeschel/NetCoMi@1.2.0`. 1.3.0 needs R ≥ 4.6 (this container has
   4.3.3) and mia; its NEWS lists no change to these computations. Its
   imports that the correlation-input path never calls (mia, mixedCCA,
   phyloseq, SPRING, SpiecEasi, SummarizedExperiment, Biobase, qgraph)
   were replaced by stand-in packages (one erroring function per imported
   name, `data-raw/expected/assoc_network.r` says how); igraph 1.6.0,
   WGCNA, fdrtool, orca, pulsar and the rest were real (Ubuntu and CRAN).
   **Exposure record:** read NetCoMi's `DESCRIPTION`, `NEWS.md`, the
   import lines of `NAMESPACE`, and the Rd of `netConstruct`,
   `netAnalyze`, `netCompare` and `diffnet`; listed (not printed) the
   names of its namespace and called two internal functions,
   `.calcJaccard()` and `.sigTestRand()`, for their outputs only. No R
   source was opened or printed; no error message showed code.
6. **Agreement.** Fixtures: `data-raw/expected/assoc_network.r` →
   `tests/tinytests/expected/assoc_network.rds` (75 KB): 12 pairs of
   correlation matrices (8–30 features) over all three signs, both degree
   kinds, five hub settings and thresholds 0.2–0.4; edge cases; the
   Jaccard grid; `diffnet()`; igraph on five hand-built graphs. Plus a
   sweep (not stored) of 160 networks (2 seeds × 40 settings × 2 groups,
   8–40 features). Largest differences, all networks:
   degree 8.9e-16, betweenness 0, closeness 7.8e-16, eigenvector 3.2e-14,
   global properties 8.9e-16; hubs identical and clusterings identical
   (same partition, same 0s) in all 160. Jaccard indices and both
   probabilities: 2e-16. ARI value: identical. Fisher's z p-values:
   2e-6 relative at p ≈ 1e-11, because NetCoMi computes `1 - pnorm()`
   (tested absolutely at 1e-14). igraph: betweenness, harmonic sums,
   eigenvector, Barrat and fast greedy modularity at 1e-10 or better;
   partitions identical except the star, whose two partitions differ in
   modularity only by rounding (0 vs −7e-17). Base-R implementations
   (Floyd–Warshall, path counting, `eigen()`, Barrat's and Newman's
   formulas): 1e-12.
7. **NetCoMi's definitions, as its outputs show them** (all now in the
   Rd):
   - degree ÷ (N − 1), N all features (igraph's normalised degree of the
     whole network); **weighted degree is the raw strength**;
   - betweenness ÷ (k − 1)(k − 2)/2, k the LCC size;
   - **closeness = avDiss_LCC · Σ 1/dist ÷ (k − 1)**: shortest paths are
     divided by the LCC's *average dissimilarity*, the mean of d over all
     the LCC's pairs counting **1 for non-edges** (its `avDiss` with
     `avDissIgnoreInf = FALSE`);
   - **average path length** = mean finite shortest path over connected
     pairs ÷ the whole network's average dissimilarity (1 without edges);
     `average_dissimilarity` was added to `global` for it;
   - **modularity is unweighted** (Newman) on the weighted fast greedy
     clusters, with cluster 0 as one group;
   - **hubs: the type-7 quantile over all N features** (zeros outside the
     LCC included), strictly above it; the same for the "most central" sets
     of `_compare()` (so `quantile` isn't over the LCC, as the plan had it);
   - **clusters: a one-node cluster is 0** (not only isolated features),
     and **igraph never returns one cluster of all nodes**: when the best
     modularity is at the final merge of a connected network it returns the
     partition one merge earlier (shown with K3, K4, K6, a star, and
     disconnected unions in igraph 1.6.0). `C_network_greedy()` does the
     same. Ties in the merge gain go to the first pair (i, j); no fixture
     differed from igraph through a tie;
   - `positive_fraction` with `sign = 'positive'` counts the dropped
     negative edges in its denominator (NetCoMi's `pep` with "signedPos");
   - an edge with **d = 0 (r = 1) is not used by shortest paths** (NetCoMi's
     graph of dissimilarities reads 0 as "no edge"), and one with **s = 0
     (r = −1, signed) is no edge at all**; the plan's 1e-12 length was
     dropped;
   - clustering coefficient 0 when no node has two neighbours.
8. **Kept different from NetCoMi** (documented): no switch to weighted
   degree for a complete network (NetCoMi switches silently; the fixture
   checks ecodive's `weighted_degree = TRUE` against it); an empty network
   gives zeros (NetCoMi errors); hub and most-central comparisons ignore
   differences within a relative 1e-12 of the quantile (NetCoMi's hub for
   two exactly tied eigenvector values of 1 depends on igraph's rounding;
   the fixture `edge_r1` happens to agree); only these global properties
   (no natural/vertex/edge connectivity, graphlets, or LCC-only copies);
   p-values adjusted with `p.adjust()` only.
9. **Jaccard probabilities** (`.calcJaccard()` outputs for every set size
   up to 6): the shared count is **Binomial(|A ∪ B|, 1/3)**: each feature
   in either set is in both with probability 1/3 (Real & Vargas's
   equiprobable model); they don't depend on N. `p_less = P(K ≤ k)`,
   `p_greater = P(K ≥ k)`; two empty sets give J = 0 and probabilities 1.
   Matched at 2e-16.
10. **ARI** uses the `cluster` columns as they are (cluster 0 is one
    cluster; matched NetCoMi's `randInd` exactly). **Its test is two-sided**
    in NetCoMi's outputs (|ARI_b| ≥ |ARI|; compared over 20,000
    permutations on three cases: 0.2940 vs 0.2944, 0.6564 vs 0.6624, 1 vs
    1, where the one-sided values would be 0.19, 0.75, 0.48), though its
    Rd describes a one-sided test. ecodive is two-sided, permutes y's labels
    over the features with `perm_source` streams (new `C_perm_indices()`
    in `src/permute.c`), and adds 1 to both counts (NetCoMi doesn't, and
    gives p = 0).
11. **Group-permutation test**: as NetCoMi's Rd describes (relabel samples
    keeping group sizes, re-estimate, |difference| ≥ observed, pseudo-count),
    through `feature_assoc()` and `assoc_network()` with `assoc_args` /
    `network_args`. Permutation b's association runs get the seed
    `pair_seed(seed, b)` (`C_pair_seed()`), replacing any `seed` in
    `assoc_args`. It is **not compared with NetCoMi** (NetCoMi re-estimates
    with its own measures; ecodive has no Pearson association, and NetCoMi's
    SparCC is SpiecEasi's), so it is tested exactly against networks rebuilt
    by hand for the same permutations (the internal hook `group_perms` of
    `network_compare_core()`; `ari_perms` likewise). Identical for `cpus`
    1, 3, 4. Results: `nodes` is long (feature × centrality), which keeps
    `p_value`/`p_adjusted` one column each; a `clusters` table (each
    network's clusters and hubs) was added for the compare plot.
12. **C** (`src/assoc_network.c`): CSR adjacency; Dijkstra with a lazy
    binary heap and Brandes' accumulation, equal paths at a relative 1e-12;
    sources in fixed blocks of 32, each block's betweenness in its own row,
    rows summed in block order (identical for `cpus` 1, 2, 3, 4, 8,
    tested); paths run on the whole network once (LCC values are the same,
    and the average path length needs every component). Fast greedy with
    sorted link lists, ΔQ recomputed from e and a. Eigenvector by
    `lanczos()` (k = 1), sign fixed and scaled to max 1. Timing (4 cores):
    3,000 features, 168,000 edges: 3.8 s (1 thread 8.6 s); paths 1.8 s,
    `assoc_select()` 0.8 s, greedy 0.8 s.
13. **`assoc_pairs()` refactored**: its selection is now `assoc_select()`
    (shared with `assoc_network()`, with a `strict` option for methods
    whose cutoff is `>`); its results are unchanged (its 63 tests pass).
    `threshold = NULL` takes the method's own sparsity (`sparse` registry
    column: non-zero entries) or a numeric `threshold` attribute (with
    `strict`), else 0.3; tested by flipping `sparse` for sparcc.
14. **Inputs**: `vlr()`/`phi_s()` refused with the plan's message; a matrix
    must be within [−1, 1]; values beyond ±1 (SparCC) are clamped for d and
    s but `association` keeps the raw value. `_diff()` takes `n_samples`
    (one or two numbers) for matrices, uses a parallel `n_pairs` when a
    result has one (SECOM), and refuses methods with `fisher_z = FALSE`.
15. **Plots**: circular layout as planned; labels placed radially outside
    the nodes; edge colours from positions 2 and 10 of
    `hcl.colors(11, "Blue-Red 3")` (the very ends read as black at thin
    widths); the compare plot titles its panels by group under one outer
    title; ASCII `>=`/`<=` in subtitles (the Unicode signs failed on
    `png()` in a C locale). `hubs` and `hub_quantile` are stored as
    attributes (they qualify the `hub` column, rule 1), and the node size
    reads the first of `hubs`.
16. **Names in results**: `global` gained `average_dissimilarity`; the
    compare result is `nodes`, `global`, `jaccard`, `ari` (named vector
    `value`, `p_value`), `clusters`, `edges_x`, `edges_y`.

### What `spiec_easi_*` (item 65) and `secom_*` (item 75) need to add

- A `FEATURE_ASSOCS` row each (the function name as `id`): spiec_easi_*
  with `kind = similarity`, `lower = -1`, `upper = 1`, `fisher_z = FALSE`,
  **`sparse = TRUE`** (then `assoc_network(threshold = NULL)` makes edges of
  the non-zero entries); secom_pearson with `fisher_z = TRUE`, **`strict =
  TRUE`**, and its chosen cutoff as a numeric top-level **`threshold`**
  attribute (then `threshold = NULL` uses it, with `>`), plus a parallel
  **`n_pairs`** attribute (pairwise-complete sample counts), which
  `assoc_network_diff()` already uses and `dist_reorder()` already keeps
  aligned; secom_distance's kind and range as its plan decides.
- Their functions must take `counts` first and no argument named
  `method` or a prefix of `method`/`counts` (`test-wrapper_args.r`
  checks); `assoc_network_compare()` re-runs them through
  `feature_assoc()`, passing `seed` only if they take one.
- Test that `assoc_network()` of their results gives the expected edges,
  and add their aliases to the cross-registry check in
  `test-assoc_network.r`.


## Outcome (from the implementing session's final report)

**Item 33, done 2026-11-20.** Commits `a4d81fd` (work, with
`feature_assoc()`), `e28f5b8` (tick), `834f9b2` (merge).

- Tests: 20,351 at home, 362 in CRAN mode (merged tree).
- NetCoMi 1.2.0 + igraph 1.6.0, 160 networks: eigenvector 3.2e-14, degree
  8.9e-16, betweenness exact; hubs and clusterings identical; Jaccard and its
  p-values 2e-16; ARI test matches in distribution.
- Clean room: two NetCoMi internals were called for outputs only (recorded
  above).
