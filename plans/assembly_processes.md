# `beta_mntd()`, `beta_mpd()`, `beta_nti()`, `raup_crick_bray()`, `assembly_processes()`: the Stegen community-assembly framework

**Priority:** 73. **Tier:** **clean room**. Every reference is GPL or
unlicensed, all verified 2026-11-20:

| Reference | Licence (how verified) | Use |
|:----------|:-----------------------|:----|
| Stegen's scripts, github.com/stegen/Stegen_etal_ISME_2013 (`bNTI_Local_Machine.r`, `Raup_Crick_Abundance.r`; last push 2018-02-23) | **none** (GitHub API: `license: null`; no LICENSE file): all rights reserved | don't open the `.r` files; README and outputs only |
| iCAMP 1.9.1 (CRAN, 2026-09-24) | GPL-2 (CRAN page) | output-only |
| NST 3.1.10 (CRAN, 2022-06-05) | GPL-2 (CRAN page) | output-only |
| picante 1.8.4 (`comdistnt()`, `comdist()`) | GPL-2 (CRAN page) | output-only |
| vegan 2.7-8 (`raupcrick()`, Chase's incidence version) | GPL-2 (CRAN page) | output-only |
| phylocom 4.2 (`comdist`, `comdistnt`) | BSD-3 (manual, "Legal") | permissive, but only the unweighted, one-directional βMNTD (see [phylo_structure.md](phylo_structure.md)) |

Work from the papers (Stegen et al. 2012, 2013, 2015; Chase et al. 2011;
Ning et al. 2020; the methods sections are detailed), the packages' Rd
pages and vignettes, and their outputs. Read [README.md](README.md) and
[phylo_structure.md](phylo_structure.md) first: this plan reuses its tree
preparation, LCA structure and null-draw threading.


## What it is

Stegen et al. (2012, 2013, 2015) estimate which ecological processes
shaped the difference between each pair of communities:

1. **βMNTD** (Fine & Kembel 2011; Stegen 2012): the abundance-weighted mean
   phylogenetic distance from each taxon of one sample to its nearest
   relative in the other.
2. **βNTI** = (βMNTD_obs − mean(βMNTD_null)) / sd(βMNTD_null), the null
   shuffling taxa across the tips of the tree (999 draws). |βNTI| > 2 means
   phylogenetic turnover is more (> +2, **heterogeneous selection**) or less
   (< −2, **homogeneous selection**) than chance.
3. For the other pairs, **RC_bray** (Stegen 2013, extending Chase et al.
   2011's Raup-Crick to abundances): where the observed Bray-Curtis falls in
   a null distribution of Bray-Curtis between communities assembled at
   random from the regional pool, rescaled to −1..1. RC_bray > +0.95:
   **dispersal limitation**; < −0.95: **homogenizing dispersal**; otherwise
   **undominated** ("drift").
4. The proportions of pairs in each class summarize a study (or each
   group's within-group pairs).

iCAMP (Ning et al. 2020) refines this by applying it within phylogenetic
**bins** (βNRI from βMPD instead of βNTI) and weighting bins by abundance;
it is planned as a later phase.

**Demand** (retrieved 2026-11-20): Stegen et al. (2013) 2,024 OpenAlex
citations 2022–26, rising from 317 (2022) to 479 (2026); Stegen et al.
(2012) 893; iCAMP (Ning et al. 2020) 1,478, rising from 136 to 420; Chase
et al. (2011) 412. Europe PMC microbiome full-text mentions 2022–26:
"βNTI" / "betaNTI" / "βMNTD" 620; "RCbray" / "Raup-Crick" 479; "iCAMP"
313. Baseline: Bray-Curtis 28,581.

**Why ecodive helps.** βNTI and RC_bray need 999 null values for every
pair of samples: 999 × n(n − 1)/2 βMNTD or Bray-Curtis evaluations
(4.9 million for 100 samples, 499 million for 1,000). Stegen's scripts run
them in R loops and take hours to days; iCAMP and NST parallelize with
`parallel` clusters but still index dense cophenetic matrices and aren't
reproducible across worker counts. This is exactly the threaded,
deterministic C ecodive is built for.


## Names

- `beta_mntd()`, `beta_mpd()`: the literature's βMNTD / βMPD; picante
  calls them `comdistnt()` / `comdist()`. Rule 2 isn't violated: they are
  the beta (between-sample) forms of `mntd()` / `mpd()`, named as the
  literature does.
- `beta_nti()`: βNTI (rule 1, the standard acronym, with the Greek letter
  spelled).
- `raup_crick_bray()`: RC_bray. Rule 2 with base `raup_crick`, the
  analytic presence/absence metric planned in
  [vegdist_metrics.md](vegdist_metrics.md); spelled out because `rc` would
  be cryptic (rule 4). Aliases in the docs: RCbray, RC_bray, "abundance-based
  Raup-Crick".
- `assembly_processes()`: the classification. Ning's iCAMP calls the
  whole-community version **QPEN** (`iCAMP::qpen()`: quantitative process
  estimation from a whole-community null model analysis); **ask the
  user** whether to name it `qpen()` instead (the name iCAMP users know) or
  keep the descriptive name. Recommended: `assembly_processes()`, with QPEN
  in the title and docs, because the Stegen framework is mostly cited
  without that acronym.
- **Weighting:** `weighted = TRUE` for `beta_mntd()` / `beta_mpd()`
  (Stegen's abundance-weighted βMNTD; picante's default is unweighted, so
  ecodive's default differs and the docs say so). As in
  [phylo_structure.md](phylo_structure.md), **ask the user** whether these
  two should instead join `METRICS` / `beta_div()` with fixed-weighting ids
  (`beta_mntd`, `weighted_beta_mntd`, ...); decide with that plan's
  question. Joining is attractive here: βMNTD is ordinated and tested like
  any distance.
- Argument names: `tree`, `taxa_dm`, `weighted`, `n_draws`, `pairs`,
  `seed`, `cpus` as elsewhere; new: `nti_threshold`, `rc_threshold`.


## Design

### Signatures

```r
beta_mntd(counts, tree = NULL, taxa_dm = NULL, weighted = TRUE,
          margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())

beta_mpd(counts, tree = NULL, taxa_dm = NULL, weighted = TRUE,
         margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())

beta_nti(counts, tree = NULL, taxa_dm = NULL, weighted = TRUE,
         n_draws = 999L, margin = 1L, pairs = NULL, seed = 0L,
         cpus = n_cpus())

raup_crick_bray(counts, n_draws = 999L, margin = 1L, pairs = NULL,
                seed = 0L, cpus = n_cpus())

assembly_processes(counts, tree = NULL, taxa_dm = NULL, groups = NULL,
                   data = NULL, n_draws = 999L, nti_threshold = 2,
                   rc_threshold = 0.95, beta_nti = NULL, rc_bray = NULL,
                   margin = 1L, seed = 0L, cpus = n_cpus())
```

- `pairs` as in every beta metric: compute only some pairs (e.g. within
  treatments), the rest `NA`. This is the main lever on cost for large n.
- `assembly_processes()` accepts precomputed `beta_nti` / `rc_bray`
  results (so users can rerun the classification with other thresholds),
  else computes them. `groups` restricts to within-group pairs (built into
  `pairs`) and reports proportions per group, as most papers do.
- `raup_crick_bray()` takes **integer counts** (individuals are drawn):
  `assert_integer_counts()`. Stegen used rarefied counts; the docs say so
  and point to `rarefy()`.

### Definitions

**βMNTD** (Stegen et al. 2012, eq. 1; picante `comdistnt(abundance.weighted
= TRUE, exclude.conspecifics = FALSE)`), with f_ik the relative abundance
of taxon i in sample k and d_ij patristic distances:

βMNTD_kl = ½ [ Σ_{i∈k} f_ik min_{j∈l} d_ij + Σ_{j∈l} f_jl min_{i∈k} d_ij ].

Unweighted: f = 1/S. A taxon present in both samples has distance 0.
✱ Confirm on picante's outputs that the weights are relative abundances
(closure per sample) and the ½ symmetrization. `exclude.conspecifics` is
not offered (Stegen doesn't use it).

**βMPD** (picante `comdist`): Σ_{i∈k} Σ_{j∈l} f_ik f_jl d_ij (unweighted:
the mean over the S_k × S_l pairs). Needed for iCAMP's βNRI (later phase)
and cheap to offer.

**βNTI** (Stegen et al. 2013, Methods): null draw r applies one random
permutation of the tree's tip labels to the whole table (picante's
`taxaShuffle()`, phylo_structure.md's `taxa_labels`) and recomputes βMNTD
for every pair; βNTI_kl = (obs − mean_r) / sd_r (n − 1 divisor). Stegen
pruned the tree to the table's taxa, so the shuffle is among the table's
features: ecodive prunes extra tips here (philr_tree() behaviour), and
says so. ✱ Confirm against iCAMP's `bNTIn.p()` outputs (Rd: "taxa
shuffle") that the shuffle is over the tree's tips and is shared by all
pairs within a draw (it doesn't change the per-pair distribution, only
cross-pair correlation, so this is seed-like either way).

**RC_bray** (Stegen et al. 2013, Methods and Supplementary; Chase et al.
2011): with the regional pool = all features of the table,
occupancy o_i = number of samples containing taxon i, and regional
abundance a_i = Σ over samples of counts:

1. For sample k (richness S_k, total N_k), a null community: draw S_k
   distinct taxa without replacement, probability ∝ o_i; give each one
   individual; distribute the remaining N_k − S_k individuals among the
   chosen taxa with replacement, probability ∝ a_i.
2. Null Bray-Curtis between the null communities of k and l, n_draws times.
3. RC = (#{null > obs} + ½ #{null = obs}) / n_draws, rescaled
   RC_bray = 2 (RC − 0.5), in −1..1.

✱ Every step above is from the papers' descriptions and must be confirmed
statistically on outputs before coding: iCAMP's `RC.pc()` and NST's
`tNST()` RC output (both document Stegen's RC_bray), and, if the user
agrees, Stegen's script run as a black box on fixtures (running it is
output-only use). Points to settle: (a) weights ∝ occupancy vs ∝ regional
abundance for the taxa draw; (b) whether a_i is total counts or mean
relative abundance; (c) tie handling (½) and the equality tolerance; (d)
whether nulls are fresh per pair or per sample (see C below); (e) the
direction of the rescaling.

**Classification** (Stegen et al. 2013, 2015): βNTI > +t: heterogeneous
selection; βNTI < −t: homogeneous selection; else RC_bray > +r: dispersal
limitation; RC_bray < −r: homogenizing dispersal; else undominated. Default
t = 2 (Stegen; iCAMP uses 1.96, a documented option), r = 0.95. Ties at the
thresholds: strict inequalities, as written in the papers (✱ check
iCAMP's `qpen()` on a crafted case).

### C, memory, threading and cost

New `src/assembly.c` (MIT). Shared with [phylo_structure.md](phylo_structure.md):
tree preparation through `C_philr_tree()`, node depths, LCA sparse table.

**Observed βMNTD for all pairs (and per null draw):**

- Per sample l, the **nearest-l-tip table** NN_l(t) for every tip t: one
  bottom-up pass (nearest present tip below each node) and one top-down
  pass (nearest outside), O(T) per sample, T = tips (after pruning, the
  features). Then βMNTD_kl = ½ [Σ_{i∈k} f_ik NN_l(i) + Σ_{j∈l} f_jl
  NN_k(j)]: O(S_k + S_l) lookups per pair.
- Memory: n × T doubles (1,000 × 20,000 = 160 MB), so compute tables for
  a block of samples at a time (e.g. 64) and do all pairs touching the
  block; block size fixed, not tied to `cpus`.
- `taxa_dm` route: O(S_k S_l) per pair over the packed triangle.
- Under a taxa shuffle π_r, sample k's taxa sit at tips π_r(i): the same
  passes with permuted marks.

**βMPD** for all pairs: per sample, the edge-abundance vector a_e (fraction
of the sample below edge e); βMPD_kl = Σ_e L_e [a_ek(1 − a_el) + a_el(1 −
a_ek)] = A_k + A_l − 2 Σ_e L_e a_ek a_el. The pair term is a weighted inner
product over edges, the same shape as weighted UniFrac: reuse
`src/unifrac.c`'s per-sample edge weight matrix and pair loop.

**Null draws, deterministic for any `cpus`:**

- Draw r uses PCG stream `rng_seed(seed, r)` for the shuffle (βNTI), or
  `(seed, r * n + k)` for sample k's null community (RC_bray).
- **Per-sample nulls for RC_bray (✱ settle d above):** null community r of
  sample k is drawn once and reused for every pair (k, ·) in draw r. Each
  pair's null distribution is exactly Stegen's (two independent null
  communities with the right richness and totals); only the correlation
  between pairs differs, which is seed-like. This turns 999 × pairs × 2
  community draws into 999 × n. If the user wants fresh draws per pair,
  use stream `(seed, r, pair)` instead, at ~n/2 times the drawing cost.
- **Rounds of draws.** Process draws in rounds of B (e.g. 32, fixed). In a
  round, tasks are (draw, sample block) units writing null βMNTD (or null
  Bray-Curtis) values into a B × n_pairs buffer (32 × 499,500 × 8 B =
  128 MB at n = 1,000; 5 MB at n = 200; B shrinks for huge n); after the
  round, a parallel pass over pair blocks folds the B values into each
  pair's running mean and sum of squares (Welford) **in draw order**, and
  for RC counts the null > obs and null = obs. Identical for any `cpus`.
- **Null communities** (RC_bray): the taxa draw is weighted sampling
  without replacement (Efraimidis-Spirakis keys, O(T) per community, or a
  Fenwick tree, O(S log T)); the N − S remaining individuals are a
  multinomial over the S chosen taxa, drawn by conditional binomials
  (O(S) per community, not O(N)). That needs a PCG binomial sampler in
  `src/rng_dist.c`: BTPE (Kachitvichyanukul & Schmeiser 1988) written from
  the paper, or ported from NumPy's `random/src/distributions/distributions.c`
  (NumPy is BSD-3; its `numpy/random/LICENSE.md` says components carry
  their own licences, so check that file's header first). Null communities
  are stored sorted by taxon (CSR) so that Bray-Curtis is a merge, O(S_k +
  S_l).

**Cost** (estimates; measure when implemented). S ≈ 500 taxa per sample,
T = 5,000 features, 999 draws, one thread at ~1e9 simple steps/s:

| n samples | pairs | βNTI: NN tables + lookups | RC_bray: draws + Bray-Curtis |
|----------:|------:|:-------------------------|:----------------------------|
| 100 | 4,950 | 0.5e9 + 5e9 ≈ 6 s | 0.1e9 + 5e9 ≈ 5 s |
| 500 | 124,750 | 2.5e9 + 1.2e11 ≈ 2 min | 1.2e11 ≈ 2 min |
| 1,000 | 499,500 | 5e9 + 5e11 ≈ 8 min | 5e11 ≈ 8 min |

Divide by the thread count (the work is embarrassingly parallel within a
round). For n ≥ 500 the docs recommend `pairs` (within-group pairs) or
fewer draws for exploration. By comparison, Stegen's R script evaluates
each null βMNTD with an R loop over cophenetic-matrix rows.

### Results

- `beta_mntd()`, `beta_mpd()`: `ecodive_dist` (via a `dist_result()`-like
  call; title "Beta Mean Nearest Taxon Distance" / "Beta Mean Pairwise
  Distance", from a registry when they join METRICS), with `clust` and
  `pairs` support.
- `beta_nti()`: `ecodive_dist` of βNTI values (negative values are normal;
  `dist` allows them, and `plot.ecodive_dist()`'s heatmap needs a diverging
  palette for this `method`: a registry column, not a flag). `method =
  'beta_nti'`. Attributes parallel to the triangle: `beta_mntd` (observed),
  `null_mean`, `null_sd`: add them to `DIST_PARALLEL` in
  `R/dist_clust.r`. Plus `n_draws`, `weighted`.
- `raup_crick_bray()`: `ecodive_dist` of RC_bray in −1..1, `method =
  'raup_crick_bray'`; attributes `bray` (observed Bray-Curtis, parallel),
  `n_draws`.
- `assembly_processes()`: `c('ecodive_assembly', 'data.frame')`, one row
  per computed pair: `sample_1, sample_2, group` (if `groups`), `beta_nti,
  rc_bray, process` (factor with the five levels, in the order homogeneous
  selection, heterogeneous selection, homogenizing dispersal, dispersal
  limitation, undominated). Attribute `proportions`: a data frame of the
  five proportions, overall and per group. Attributes `nti_threshold`,
  `rc_threshold`, `n_draws`. `print` via `print_table_result`, with a
  header label.


## Later phase: iCAMP (Ning et al. 2020)

Worth doing once the above is in, because its uptake is rising fastest
(136 → 420 citations a year) and it reuses every kernel here, but it is a
**separate item** of **high** effort:

1. Bin taxa by the phylogeny (iCAMP's `taxa.binphy.big()`: split the tree
   into clades with a minimum bin size and a phylogenetic distance limit,
   ✱ defaults believed to be 24 taxa and `ds = 0.2`; small bins merged with
   their nearest). The binning rule and defaults must be pinned down from
   iCAMP's Rd and outputs.
2. Per bin, βNRI (null: shuffle taxa within the bin's tips, βMPD) and
   RC_bray (null within the bin's taxa); classify each bin × pair with
   thresholds 1.96 and 0.95.
3. Per pair, weight each bin's process by the bin's relative abundance in
   the pair; sum to process importances.
4. Phylogenetic signal check within bins (iCAMP's `ps.bin()`, a Mantel
   correlogram): depends on `mantel_correlogram()`, planned.

Output-only fixtures from iCAMP 1.9.1 (`icamp.big()`), and the per-bin
statistics are the same kernels with masks. NST (normalized stochasticity
ratio, Ning et al. 2019) is another candidate on the same null-community
code; not planned here.


## Reference implementations

- **Output-only:** picante `comdistnt(comm, cophenetic(tree),
  abundance.weighted = TRUE)` and `comdist()` for observed βMNTD / βMPD
  (deterministic, 1e-10); iCAMP `bNTIn.p()` (βNTI), `RC.pc()` (RC_bray),
  `qpen()` (classification) for statistical comparison of the nulls;
  vegan `raupcrick()` for Chase's incidence-based version, if the user
  wants it (below). `data-raw/expected/assembly_processes.r`, on
  `ex_counts` / `ex_tree`, a 30-sample × 400-taxon rarefied table with a
  random tree (two groups with different composition so all five classes
  occur), and a tiny case for exact null enumeration (3 taxa: the full
  null distribution of RC_bray can be computed exactly in R, in the test).
- No permissive implementation of βNTI or RC_bray was found
  (2026-11-20).
- **Statistical checks:** with 9,999 draws, ecodive's null mean and sd per
  pair within 4 Monte Carlo SE of iCAMP's (at home only), and the
  proportions of the five classes close to `qpen()`'s.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- `beta_nti()`, `raup_crick_bray()`: the `ecodive_dist` heatmap, titled
  "Beta Nearest Taxon Index" / "Raup-Crick (Bray-Curtis)", with a diverging
  palette centred at 0 and the thresholds marked on the legend (from the
  registry and `n_draws`).
- `assembly_processes()`: `plot()` draws a stacked bar per group (or one
  bar) of the five proportions, from the `proportions` attribute; title
  "Community Assembly Processes"; `type = 'scatter'`: βNTI against RC_bray
  per pair with the threshold lines, coloured by `process`.


## Tests

- Observed βMNTD / βMPD equal picante at 1e-10 (both weightings, tree and
  `taxa_dm`), and an independent plain-R computation from a patristic
  matrix in the test.
- βMNTD of a sample with itself is 0; symmetric; `pairs` respected;
  `beta_mpd(weighted = FALSE)` equals the mean of the cross block of the
  patristic matrix.
- `identical()` results for `cpus` 1, 2, 4 for `beta_nti()` and
  `raup_crick_bray()`; same seed identical, different seed different.
- RC_bray: the 3-taxon case against its exact null distribution
  (Monte Carlo tolerance); two identical samples give RC_bray near −1;
  disjoint samples near +1.
- Classification: a crafted `beta_nti` / `rc_bray` pair of inputs hitting
  each class and each threshold; `groups` proportions sum to 1.
- Edge cases: n = 2; a sample with one taxon; a taxon in every sample;
  non-integer counts to `raup_crick_bray()` (error); missing tree tips
  (error); `n_draws = 1` (sd NA, βNTI NA).
- Smoke-test lines with `n_draws = 9`.


## Docs

Rd pages in a new "Community assembly" group of the reference index, next to
the phylogenetic-structure functions. Cover: the framework's logic and its
assumptions (phylogenetic signal in niche preferences, which the user
should check; Stegen et al. 2015's caveats), the thresholds, rarefied
integer counts for RC_bray, the cost and `pairs`, and "Differences from
Stegen's scripts and iCAMP" (ecodive's RNG and per-sample nulls; ✱ items
not matched exactly). Also mention Chase et al.'s incidence-based Raup-Crick
(vegan `raupcrick()`); **ask the user** whether to add it as
`raup_crick_null()` (same null code with presence only; low extra cost).
NEWS.


## References

- Stegen, J. C., Lin, X., Konopka, A. E., & Fredrickson, J. K. (2012).
  Stochastic and deterministic assembly processes in subsurface microbial
  communities. *ISME Journal*, 6, 1653–1664. doi:10.1038/ismej.2012.22
- Stegen, J. C., Lin, X., Fredrickson, J. K., et al. (2013). Quantifying
  community assembly processes and identifying features that impose them.
  *ISME Journal*, 7, 2069–2079. doi:10.1038/ismej.2013.93
- Stegen, J. C., Lin, X., Fredrickson, J. K., & Konopka, A. E. (2015).
  Estimating and mapping ecological processes influencing microbial
  community assembly. *Frontiers in Microbiology*, 6, 370.
  doi:10.3389/fmicb.2015.00370
- Chase, J. M., Kraft, N. J. B., Smith, K. G., Vellend, M., & Inouye, B. D.
  (2011). Using null models to disentangle variation in community
  dissimilarity from variation in α-diversity. *Ecosphere*, 2, art24.
  doi:10.1890/ES10-00117.1
- Fine, P. V. A., & Kembel, S. W. (2011). Phylogenetic community structure
  and phylogenetic turnover across space and edaphic gradients in western
  Amazonian tree communities. *Ecography*, 34, 552–565.
  doi:10.1111/j.1600-0587.2010.06548.x
- Ning, D., Yuan, M., Wu, L., et al. (2020). A quantitative framework
  reveals ecological drivers of grassland microbial community assembly in
  response to warming. *Nature Communications*, 11, 4717.
  doi:10.1038/s41467-020-18560-z
- Ning, D., Deng, Y., Tiedje, J. M., & Zhou, J. (2019). A general framework
  for quantitatively assessing ecological stochasticity. *PNAS*, 116,
  16892–16898. doi:10.1073/pnas.1904623116
- Zhou, J., & Ning, D. (2017). Stochastic community assembly: does it
  matter in microbial ecology? *Microbiology and Molecular Biology
  Reviews*, 81, e00002-17. doi:10.1128/MMBR.00002-17
- Kachitvichyanukul, V., & Schmeiser, B. W. (1988). Binomial random variate
  generation. *Communications of the ACM*, 31, 216–222.
  doi:10.1145/42372.42381
- Efraimidis, P. S., & Spirakis, P. G. (2006). Weighted random sampling
  with a reservoir. *Information Processing Letters*, 97, 181–185.
  doi:10.1016/j.ipl.2005.11.003
