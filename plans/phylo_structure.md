# `mpd()`, `mntd()`, `mpd_ses()`, `mntd_ses()`: phylogenetic community structure (NRI / NTI)

**Priority:** 72. **Tier:** clean room against picante (GPL-2,
output-only), with **phylocom 4.2 (BSD-3) as a permissive reference** for
the definitions and the null-model loop. picante 1.8.4 (CRAN, published
2026-08-23) is `License: GPL-2` (CRAN page, verified 2026-11-20): read only
its Rd pages, vignette and outputs, never its `R/` or `src/`. phylocom 4.2
(github.com/phylocom/phylocom, by Webb, Ackerly & Kembel, the same authors)
has no LICENSE file, but its manual (`doc/phylocom_manual.tex`, section
"Legal") says "distributed under a BSD ('3-clause') license", (c) 2009
Campbell O. Webb, David D. Ackerly, Steven W. Kembel (verified
2026-11-20). Its `src/combase.c` (`MeanDistance()`,
`MeanMinimumDistance()`, the comstruct null loop, `IndependentSwap()`,
`TrialSwap()`) may be read and ported; it was read while writing this plan.
Read [README.md](README.md) first, and [dpcoa.md](dpcoa.md) and
[philr.md](philr.md) for the tree code this reuses.


## What it is

Webb (2000) and Webb et al. (2002) measure whether the taxa in a sample are
more closely related (phylogenetic clustering) or more distantly related
(overdispersion) than expected by chance:

- **MPD**, the mean pairwise phylogenetic (patristic) distance among the
  sample's taxa: sensitive to tree-wide (deep) structure.
- **MNTD**, the mean distance from each taxon to its nearest relative in the
  same sample: sensitive to structure near the tips.
- Both optionally **abundance-weighted** (picante `abundance.weighted`,
  phylocom `-a`).
- **Standardized effect sizes** against a null model that randomizes the
  sample or the tree: SES = (obs − mean(null)) / sd(null).
  **NRI = −SES(MPD)**, **NTI = −SES(MNTD)**, so positive values
  mean clustering (Webb et al. 2002). A rank-based p-value comes with them.

They are the standard "phylogenetic community structure" summary in
microbial ecology (habitat filtering versus competition arguments), usually
via picante's `ses.mpd()` / `ses.mntd()` (Kembel et al. 2010).

**Demand** (retrieved 2026-11-20): picante (Kembel et al. 2010) 3,272
OpenAlex citations 2022–26; Webb et al. (2002) 1,608; Europe PMC microbiome
full-text mentions of "NRI" AND "NTI": 35 (a floor: papers also write
"net relatedness index", "SES.MPD", ...). Baseline: Bray-Curtis 28,581
mentions.

**Why ecodive helps.** The null randomizations are the cost: picante loops
in R over 999 draws × samples, each indexing the dense D × D cophenetic
matrix (`cophenetic()` of a 20,000-tip tree is 3.2 GB). ecodive computes
MPD and MNTD from the tree directly in O(S log S) per sample and draw (S =
the sample's richness), threads the draws, and gives identical results for
any `cpus`.


## Names

- `mpd()`, `mntd()`: the standard acronyms (rule 1), as picante, phylocom and
  PhyloMeasures name them. No clash with base R or vegan.
- `mpd_ses()`, `mntd_ses()`: rule 2 (`<base>_<modifier>`), so the family
  sorts together; picante's `ses.mpd()` / `ses.mntd()` with the dot dropped
  and the order fixed. The result has an `nri` / `nti` column, and the docs
  name NRI / NTI prominently (aliases in "Same names in other packages").
- **Ask the user:** whether to also export `nri()` and `nti()` as thin
  wrappers (returning only the index). Recommended: no; one function per
  computation, with the familiar names in the docs and the result.
- **Abundance weighting is an argument, `weighted = FALSE`** (picante
  `abundance.weighted`, phylocom `-a`): every reference treats it as an
  option (README rule). It means "weight by abundance", not UniFrac's
  `weighted_*` family; the docs say so.
- **METRICS:** `mpd()` / `mntd()` don't join `alpha_div()` in this item,
  because `alpha_div()` passes only its own formals (`mget(metric$params)`)
  and has no `weighted`. **Ask the user** whether they want METRICS rows
  instead, which needs UniFrac-style fixed-weighting ids (`mpd`,
  `weighted_mpd`, `mntd`, `weighted_mntd`); [assembly_processes.md](assembly_processes.md)
  asks the same for `beta_mpd()` / `beta_mntd()`. Decide both together.
- New glossary rows ([parameter_names.md](parameter_names.md) §1):
  `weighted` (abundance weighting of a phylogenetic mean), `null_model`,
  `n_swaps` (an exact count of swaps per null matrix, `n_<unit>` rule).
  `n_draws` (existing: Monte Carlo draws) is the number of null
  communities, picante's `runs`. `taxa_dm` is `dpcoa()`'s.


## Design

### Signatures

```r
mpd(counts, tree = NULL, taxa_dm = NULL, weighted = FALSE, margin = 1L,
    cpus = n_cpus())

mntd(counts, tree = NULL, taxa_dm = NULL, weighted = FALSE, margin = 1L,
     cpus = n_cpus())

mpd_ses(counts, tree = NULL, taxa_dm = NULL, weighted = FALSE,
        null_model = 'taxa_labels', n_draws = 999L, n_swaps = 1000L,
        margin = 1L, seed = 0L, cpus = n_cpus())

mntd_ses(...)   # the same arguments as mpd_ses()
```

- `tree`: a `phylo` with branch lengths, found as `dpcoa()` finds it
  (embedded in phyloseq / rbiom / TSE objects). `taxa_dm`: any taxon
  dissimilarity instead (picante's `dis` can be any matrix, e.g. functional
  trait distances); exactly one of the two, as in `dpcoa()`.
- `null_model`, picante's names in snake_case: `'taxa_labels'`,
  `'richness'`, `'frequency'`, `'sample_pool'`, `'phylogeny_pool'`,
  `'independent_swap'`, `'trial_swap'`. Exact match, with the dotted
  spellings accepted (`'taxa.labels'`), as `philr()` accepts its weight
  names.
- `n_swaps`: picante's `iterations` (swaps per null matrix), used only by
  the two swap models.
- Counts may be proportions (unweighted MPD / MNTD use only presence).

### Definitions (Webb 2002; picante's Rd; phylocom `combase.c`)

For a sample with taxa i = 1..S, abundances x_i, and patristic distances
d_ij:

- **MPD**, unweighted: Σ_{i<j} d_ij / (S(S − 1)/2).
- **MPD**, weighted: Σ_i Σ_j x_i x_j d_ij / Σ_i Σ_j x_i x_j, over **all
  ordered pairs including i = j** (phylocom: "conspecifics are counted as
  part of the average distance (we count individuals here)"), i.e.
  Σ_{i,j} x_i x_j d_ij / (Σ x)². Scale-free, so counts and proportions give
  the same value. ✱ Confirm on picante's outputs that it uses the same
  denominator (the alternative is Σ_{i≠j} x_i x_j); this is the first item
  of the protocol below.
- **MNTD**: m_i = min_{j≠i} d_ij; unweighted mean of m_i, or weighted Σ x_i
  m_i / Σ x_i.
- **Fewer than 2 taxa:** ✱ picante's value (NA expected; phylocom gives
  MNTD 0). Match picante.

On a tree, with W = Σ x and w_e the abundance below edge e (length L_e):
Σ_{i,j} x_i x_j d_ij = 2 Σ_e L_e w_e (W − w_e), and for presence (x = 1)
Σ_{i<j} d_ij = Σ_e L_e s_e (S − s_e). Edges above the taxa's most recent
common ancestor contribute 0, so only the subtree spanning the sample's taxa
matters.

### Null models (picante's `ses.mpd` Rd and `randomizeMatrix` Rd; phylocom's manual)

| `null_model` | What a draw randomizes | Keeps |
|:-------------|:-----------------------|:------|
| `taxa_labels` | one random permutation of the tips' labels, applied to every sample | the community matrix |
| `richness` | each sample's values moved to S random columns of the table (all features, including ones absent everywhere) | each sample's richness and abundance multiset |
| `frequency` | each feature's values permuted across samples | each feature's occurrence frequency |
| `sample_pool` | each sample's taxa redrawn uniformly without replacement from the taxa present in ≥ 1 sample | richness |
| `phylogeny_pool` | as `sample_pool`, from all tips of the tree (or all rows of `taxa_dm`) | richness |
| `independent_swap` | `n_swaps` checkerboard swaps (Gotelli 2000) | row and column totals of presence |
| `trial_swap` | `n_swaps` attempted swaps (Miklós & Podani 2004) | as above |

Details to confirm on picante's outputs (statistically: means, variances,
which taxa can appear) before coding, since picante draws with R's RNG and
can't be matched draw for draw:

- ✱ `taxa_labels` and `phylogeny_pool` draw from **every tip of the tree**,
  not only the table's features (picante shuffles `dis`'s labels, and
  `dis = cophenetic(tree)` has every tip). So ecodive must keep tips that
  aren't features for these two models (unlike `philr()`, which prunes),
  and the docs repeat picante's advice to prune the tree to the study's
  taxa first.
- ✱ Under `sample_pool` / `phylogeny_pool` / `richness`, the sample's own
  abundances follow the shuffle (matters only for `weighted = TRUE`).
- ✱ Swaps: each draw starts from the observed matrix (independent chains,
  needed for threading; phylocom instead continues one chain after a
  burn-in). Abundance handling under swaps with `weighted = TRUE`: settle
  from picante's docs and outputs; if it can't be pinned down, refuse
  `weighted = TRUE` with the swap models in the first version and say so.
  Port the swap moves from phylocom's `IndependentSwap()` / `TrialSwap()`
  (BSD-3; credit in `inst/COPYRIGHTS` and `Authors@R` if code is ported).

### Statistics (picante's `ses.mpd` Rd, Value section)

Per sample: `n_taxa`, `obs`, `null_mean`, `null_sd` (n − 1 divisor, as R's
`sd()`), `z` = (obs − null_mean)/null_sd, `nri` (or `nti`) = −z, `rank`
(rank of obs among obs and the draws; ✱ ties as `rank()`'s default
"average", confirm), `p_value` = rank / (n_draws + 1) (picante's
`mpd.obs.p`, the lower tail: small p = clustering; document that), and
`n_draws`. A draw where MPD / MNTD is NA (fewer than 2 taxa, possible under
`frequency`) is skipped in mean and sd (✱ check picante: `na.rm`?), and
`n_draws` counts the draws used.

### Algorithm and C

New `src/phylo_structure.c` (MIT; with a phylocom notice only for any ported
swap code).

1. **Tree preparation** (R + C, once): reuse `C_philr_tree()` (prune,
   collapse single-child nodes, resolve polytomies with zero-length
   branches, preorder ids). Collapsing and resolving don't change
   patristic distances, as `dpcoa()` relies on. For `taxa_labels` /
   `phylogeny_pool`, call it with **every tip** as a "feature" (the table's
   features first, then the other tips), so no tip is pruned. Branch lengths
   are required; negative or missing ones are an error.
2. **LCA structure:** node depths from the root (sums of branch lengths,
   computed root-down in preorder), an Euler tour and a sparse table for
   O(1) lowest-common-ancestor queries, and each tip's DFS (preorder)
   position. Memory O(T log T) for T tips: 50,000 tips ≈ 15 MB. Then
   d_ij = depth_i + depth_j − 2 depth(lca(i, j)).
3. **Per sample (observed or one draw):** sort the S present tips by DFS
   position (S log S), build the **virtual tree** (the S tips and the LCAs of
   DFS-adjacent pairs; at most 2S − 1 nodes) with a stack, edge lengths =
   depth differences. On it:
   - MPD: one postorder pass for s_e (or w_e), then the edge sums above:
     O(S).
   - MNTD: two passes (bottom-up nearest present tip below each node, with
     the best two from different children; then top-down nearest outside),
     giving every tip's nearest other present tip: O(S).
   Cost O(S log S) per sample and draw, never O(T) and never a T × T
   matrix. A tip pair at distance d that is deep in the tree loses about
   eps × depth / d relative precision to the depth subtraction; fixtures on
   deep trees with short tip branches check that MNTD stays within 1e-10.
4. **`taxa_dm`:** O(S²) per sample over the packed triangle (MPD: weighted
   pair sums; MNTD: row minima). Dense `taxa_dm` is the user's choice.
5. **Draws:** draw r uses its own PCG stream, `rng_seed(seed, r)`
   (`src/rng_dist.c`), and `pcg32_bounded()` for shuffles (Fisher-Yates,
   or Floyd's algorithm for S of D without replacement). Models that act on
   one sample (`richness`, `sample_pool`, `phylogeny_pool`) use the stream
   `(seed, r * n + sample)`, so each (draw, sample) is independent work;
   `taxa_labels`, `frequency` and the swaps build one null matrix per draw
   and score every sample on it.
6. **Threading:** tasks are fixed blocks of draws (e.g. 8 draws per task,
   independent of `cpus`); each task writes its draws' values into a
   `n_draws × n` buffer (999 × 2,000 samples = 16 MB; process samples in
   blocks beyond that), and the mean, sd and rank are summed per sample in
   draw order after `run_parallel()`. Results are identical for any `cpus`.
   Per-thread scratch: the virtual-tree arrays (O(max S)), and for the
   matrix models one null matrix (CSR, the table's nnz; swaps need a dense
   presence bitset, n × D bits, 1.25 MB for 500 × 20,000).
7. **Cost** (one thread, estimated; measure when implemented): 500 samples,
   S ≈ 400, 999 draws: ≈ 500 × 999 × 400 × log₂400 ≈ 1.7e9 elementary steps,
   a few seconds; picante needs one `cophenetic()` (O(T²) memory) and
   ~500,000 R-level `mpd()` sample evaluations.

`mpd()` / `mntd()` are step 3 for the observed table only, parallel over
samples.

### Results

- `mpd()`, `mntd()`: a numeric vector named by sample, as `faith()`
  returns. No attributes besides names.
- `mpd_ses()`, `mntd_ses()`: `c('ecodive_mpd_ses', 'data.frame')` (resp.
  `ecodive_mntd_ses`), one row per sample: `sample, n_taxa, obs, null_mean,
  null_sd, z, nri` (`nti`), `rank, p_value`. Attributes (results only):
  `method` (`'mpd_ses'`), `null_model`, `n_draws`, `n_swaps` (swap models
  only), `weighted`. Printing: `print.ecodive_mpd_ses <- print_table_result`
  and a header label in `print_df_label()`.


## Reference implementations

- **Output-only:** picante 1.8.4 `mpd()`, `mntd()`, `ses.mpd()`,
  `ses.mntd()` (`dis = cophenetic(tree)`, `abundance.weighted` both ways,
  each `null.model`, `runs = 999`), via `data-raw/expected/phylo_structure.r`.
  Cases: `ex_counts` / `ex_tree`; picante's own `phylocom` example data
  (`data(phylocom)`: its `phylo` and `sample`); a 300-tip `rtree()` with a
  sparse 40-sample table, extra tips not in the table, and one
  single-taxon sample; a deep tree with short tip branches (precision); a
  trait-distance `taxa_dm`. Store obs values (deterministic, 1e-10) and the
  null summaries (for statistical checks only).
- **Permissive:** phylocom 4.2 (BSD-3) `comstruct -m 0..4 [-a]` outputs as a
  second check of the observed MPD / MNTD (phylocom prints 4 decimals and
  computes in `float`, so ~1e-5 only).
- **Black-box protocol (do first):** (a) weighted MPD denominator with and
  without the diagonal; (b) MPD / MNTD for S = 0, 1; (c) `mpd.obs.rank`
  under ties (a 2-taxon sample under `taxa_labels` has many tied null
  values); (d) whether `taxa_labels` / `phylogeny_pool` ever place
  non-table tips in a sample (a tree with 10 extra tips: count draws whose
  null mean moves); (e) NA draws under `frequency`; (f) abundances under the
  pool models and the swaps. Record the answers in this plan.
- **Later, optional:** exact null moments without randomization for the
  pool models (Tsirogiannis, Sandel & Cheliotis 2012; Tsirogiannis & Sandel
  2016, PhyloMeasures, GPL-3: clean room from the papers): E[MPD] and
  Var[MPD] for S taxa drawn uniformly from the tree have closed forms
  computable in O(T). Would be `null_model = 'exact'`; ask the user.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Title: "Net Relatedness Index (NRI)" / "Nearest Taxon Index (NTI)"; the
  input is a count table, so nothing is appended.
- `plot()` of an SES result: a dot chart of `nri` (`nti`) per sample, ordered,
  with reference lines at ±1.96 and points filled where `p_value` < 0.05 or
  > 0.95; `groups =` (optional, as `plot_ranked()` takes groups elsewhere)
  colours and splits samples. Subtitle from attributes: null model and
  `n_draws`. Reads only columns and attributes.
- `mpd()` / `mntd()` return plain vectors: no `plot()` method.


## Tests

- Observed MPD / MNTD (both weightings) against picante at 1e-10, tree and
  `taxa_dm` routes; tree route equals the `taxa_dm` route on
  `cophenetic()` (computed in the fixture script, stored).
- An independent R computation in the test from `as.matrix()` of a
  patristic matrix built by summing edges (plain R, no ape): MPD / MNTD
  equal at 1e-12.
- Invariances: `weighted = TRUE` is unchanged by scaling a sample; polytomy
  resolution and single-child nodes don't change results; extra tips don't
  change observed values.
- Nulls: `identical()` for `cpus` 1, 2, 4 (all seven models); same `seed`
  gives identical results, a different one different results; null means
  and sds within Monte Carlo error of picante's (e.g. |Δ mean| < 4 SE over
  9,999 draws, done at home only); `richness` keeps each sample's richness,
  `frequency` each feature's frequency, swaps both (checked on the null
  matrices through an internal test hook, or by construction).
- Edge cases: S = 0, 1, 2; all samples identical; a zero-length tree;
  `n_draws = 1` (sd NA); a feature missing from the tree (error); both or
  neither of `tree` / `taxa_dm`; unknown `null_model` (error lists them).
- Smoke-test lines for all four functions (CRAN mode: `n_draws = 9`).


## Docs

Rd pages in a new "Phylogenetic community structure" group of the
reference index (with Faith's PD linked): the sign convention (NRI = −SES,
positive = clustering), which null model answers which question (Kembel
2009; Hardy 2008), the tree-pruning advice, the lower-tail p-value, and
"Differences from picante" (ecodive's RNG, so null values differ draw by
draw; independent swap chains; trees used without a D × D matrix; any ✱
item that couldn't be matched). "Same names in other packages":
`picante::mpd`, `picante::mntd`, `PhyloMeasures::mpd.query`. NEWS.


## References

- Webb, C. O. (2000). Exploring the phylogenetic structure of ecological
  communities: an example for rain forest trees. *American Naturalist*,
  156, 145–155. doi:10.1086/303378
- Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
  Phylogenies and community ecology. *Annual Review of Ecology and
  Systematics*, 33, 475–505. doi:10.1146/annurev.ecolsys.33.010802.150448
- Webb, C. O., Ackerly, D. D., & Kembel, S. W. (2008). Phylocom: software
  for the analysis of phylogenetic community structure and trait
  evolution. *Bioinformatics*, 24, 2098–2100.
  doi:10.1093/bioinformatics/btn358
- Kembel, S. W., Cowan, P. D., Helmus, M. R., et al. (2010). Picante: R
  tools for integrating phylogenies and ecology. *Bioinformatics*, 26,
  1463–1464. doi:10.1093/bioinformatics/btq166
- Kembel, S. W. (2009). Disentangling niche and neutral influences on
  community assembly: assessing the performance of community phylogenetic
  structure tests. *Ecology Letters*, 12, 949–960.
  doi:10.1111/j.1461-0248.2009.01354.x
- Gotelli, N. J. (2000). Null model analysis of species co-occurrence
  patterns. *Ecology*, 81, 2606–2621.
  doi:10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2
- Miklós, I., & Podani, J. (2004). Randomization of presence-absence
  matrices: comments and new algorithms. *Ecology*, 85, 86–92.
  doi:10.1890/03-0101
- Tsirogiannis, C., Sandel, B., & Cheliotis, D. (2012). Efficient
  computation of popular phylogenetic tree measures. *Algorithms in
  Bioinformatics* (WABI 2012), LNCS 7534, 30–43.
  doi:10.1007/978-3-642-33122-0_3
- Bender, M. A., & Farach-Colton, M. (2000). The LCA problem revisited.
  *LATIN 2000*, LNCS 1776, 88–94. doi:10.1007/10719839_9


## Decisions made while implementing

Answers from the user (one round of questions, 2026-11-27):

- **No `nri()` / `nti()` wrappers.** The SES results carry `nri` / `nti`
  columns, and the docs name NRI / NTI prominently.
- **METRICS rows: both join.** `mpd()`, `weighted_mpd()`, `mntd()` and
  `weighted_mntd()` are four exported alpha metrics with fixed weighting,
  UniFrac-style (ids `mpd`, `weighted_mpd`, `mntd`, `weighted_mntd`;
  aliases `wmpd`, `wmntd`, `mnnd`). The SES functions keep `weighted =`.
  Item 73 should give `beta_mpd()` / `beta_mntd()` METRICS rows the same
  way (`beta_mpd`, `weighted_beta_mpd`, ...), as the user decided both
  together.
- **`null_model = 'exact'` included now**, for MPD (see below).

My own calls:

- **`alpha_div()` gains `taxa_dm = NULL`, as its last argument** (after
  `cpus`, so positional calls are unaffected): `alpha_div()` passes a
  metric's formals with `mget()`, so every formal of a metric must exist
  in `alpha_div()`.
- **Black-box protocol, from picante 1.8.4's outputs** (source never read):
  - (a) Weighted MPD's denominator includes the diagonal:
    sum over i, j of x_i x_j d_ij / W² (as phylocom's `MeanDistance()`).
  - (b) S = 0 or 1: picante's `mpd()` / `mntd()` give `NA` (phylocom's
    MNTD gives 0). `ses.*` gives `NaN` null means; ecodive gives `NA`.
  - (c) `mpd.obs.rank` averages ties (half-integer ranks seen).
  - (d) `taxa.labels` and `phylogeny.pool` draw from every tip of `dis`
    (null means match the mean over all 40 tips of a tree with 20 table
    taxa, not over the table's).
  - (e) NA draws (frequency) are left out of mean and sd (`na.rm`), and
    sort last in the rank; `runs` (the `p` denominator) still counts them.
    ecodive divides by the draws used plus one instead (documented).
  - (f) Abundances follow the shuffled taxa under the pool models and
    `richness`; under the swaps, each taxon keeps its values (column
    multisets kept, row sums not), as in phylocom's swap code. So
    `weighted = TRUE` is supported for every model.
  - Also found: `ses.mpd(null.model = 'sample.pool')` gives exactly the
    same communities as `'richness'` (identical results with the same
    seed), i.e. it draws from every column of the table, empty ones
    included. ecodive's `sample_pool` draws only from taxa present in some
    sample, as documented; they agree when no column is empty.
  - Also found: picante's `richness` and `frequency` shuffles are not
    uniform (each position's chance of a value varies 0.465–0.536 where
    0.5 is expected, over 20,000 draws), so its null means for these two
    models are off by 5–11 standard errors from a uniform shuffle's.
    ecodive's are uniform; the tests check them against exact moments
    (richness, sample_pool) and a base-R uniform shuffle (frequency)
    instead of picante.
  - The draws of a seeded `ses.mpd()` call are exactly those of
    `taxaShuffle(dis)` / `randomizeMatrix(samp, model)` called `runs`
    times after the same `set.seed()`, so the fixture script rebuilds
    them and the tests feed them to ecodive's internal hook.
- **Internal hook:** `phylo_ses(..., null_tables = )`, a list of null
  tables (samples × tips, matched by column name) used instead of draws,
  as `assoc_network_compare()`'s `ari_perms` and `dmm()`'s `init`. Given
  picante's tables, ecodive's `null_mean`, `null_sd`, `z`, `rank` and
  `p_value` equal picante's to ≤ 1.8e-14, for all 6 rebuildable models ×
  both weightings × MPD / MNTD × tree and `taxa_dm` routes.
- **Ties** in the rank: a null value within 1e-12 (relative) of `obs`
  counts as a tie, so values equal in exact arithmetic tie whatever the
  rounding (documented).
- **Swaps:** each draw is a chain of its own from the observed table (as
  picante; phylocom continues one chain). Independent swap counts
  successful swaps, trial swap attempts. When every pair of samples is
  nested no swap is possible: picante's independent swap would loop
  forever; ecodive detects it in O(nnz) (samples ranked by richness: each
  taxon must be in exactly the first k samples) and warns, using the
  observed table. The swap state is a dense samples × taxa `int` map per
  thread (4 bytes per cell), cleared entry by entry after each draw.
  `pcg32_bounded()` replaces phylocom's biased `random() % n`.
- **Streams:** `taxa_labels`, `frequency` and the swaps use PCG stream
  `r` for draw r; `richness`, `sample_pool` and `phylogeny_pool` stream
  `r * n + i` for sample i. Draws are summed per sample in draw order
  (Welford mean and sum of squares) after each batch of draws, so results
  are identical for any `cpus`, and batch size can depend on `n` and
  `cpus`.
- **Virtual tree:** instead of sorting the 2S − 1 nodes, the tree is built
  with a stack holding the root path to the latest tip (the classic
  construction), using only the LCAs of tips adjacent in preorder (a
  range minimum over the LCAs of globally adjacent tips, from a sparse
  table, with a log table). Nodes leave the stack in postorder, which is
  reversed for the numbering. Tips are radix-sorted (8-bit digits) above
  48. This made draws about 4× faster than a heapsort-and-sort first
  version. `ptree_t`, `ptree_init()`, `ptree_ws_init()`, `ptree_vtree()`
  and `ptree_pair_sum()` are declared in `src/ecodive.h` for item 73.
- **Exact moments** (MPD only, both weightings): derived from first
  principles, because the Tsirogiannis et al. papers were not reachable
  (Wiley returned 403, Springer redirects to a login) and PhyloMeasures
  (GPL-3, archived from CRAN) could not be installed here (the install was
  refused by the session's permission classifier). With distances centred
  on their mean μ over distinct pairs: E[Z] = μ Σa,
  Var[Z] = A2 P2/(T)_2 + A1 (R2 − P2)/(T)_3 + A0 (2 P2 − 4 R2)/(T)_4
  (see `src/phylo_structure.c`). Centring keeps the cancellation mild. On
  a tree, Σ_{a<b} d, each tip's row sum, and Σ_{a<b} d² take
  two passes. Checked against full enumeration of every ordered draw on a
  7-tip and a 9-tip tree (≤ 4.4e-16), and against `phylogeny_pool` /
  `taxa_labels` draws. With `'exact'`, `rank`, `p_value` and `n_draws` are
  `NA`. **No exact MNTD**: I found no tractable exact variance (the pairwise
  terms need the joint distribution of two tips' nearest neighbours), so
  `mntd_ses(null_model = 'exact')` is an error saying so.
- **Result:** `n_draws` column = draws used (non-NA); `p_value = rank /
  (n_draws + 1)`. A sample with fewer than 2 taxa: `obs`, `z`, `nri`,
  `rank`, `p_value` `NA` (null mean/sd still given where defined).
  Attributes `method`, `null_model`, `n_draws`, `n_swaps` (swaps only),
  `weighted`; `display$title` "Net Relatedness Index (NRI)" / "Nearest
  Taxon Index (NTI)".
- **`null_model` matching:** exact after lower-casing and dropping `.`,
  `_` and spaces, so picante's `'taxa.labels'`, `'independentswap'` work.
- **phylocom `comstruct`:** its source (BSD-3) was fetched through
  jsDelivr and compiled, but running the binary was refused by the
  session's permission classifier, so no phylocom outputs were compared;
  picante's are the reference.
- **Reference index:** a "Phylogenetic community structure" subtitle under
  Community Assembly, pointing to `faith()`.
- **README:** `mpd()`, `mntd()` in the Alpha diversity row; `mpd_ses()`,
  `mntd_ses()` in the Community assembly row (no new row).
