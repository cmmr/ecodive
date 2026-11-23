# `dpcoa()`: double principal coordinate analysis

**Priority:** 26, second of the clean-room methods (after `philr()`). **Tier:** clean
room. No permissive implementation exists, so write it from Pavoine,
Dufour & Chessel (2004). The reference, ade4 `dpcoa()` (GPL), is
output-only: run it for fixtures, but don't read its source. Read
[README.md](README.md) first.


## What it is

DPCoA ordinates samples and taxa in one space, using a phylogeny (or any
taxon dissimilarity). Taxa are placed by a PCoA of their patristic
distances (in fact √ of them, so that the space is Euclidean). Each
sample sits at the abundance-weighted centroid of its taxa. Distances
between samples are related to Rao's quadratic entropy, a phylogenetic
diversity. It's in phyloseq (`ordinate(method = "DPCoA")`) and mia
(`runDPCoA`).


## Design

- **Signature:** `dpcoa(counts, tree = NULL, margin = 1L, n_axes = 3L, cpus = n_cpus())`.
  - `tree` is a `phylo` object, pulled from phyloseq, rbiom or TSE
    objects as the UniFrac functions do (see `validate_counts()` /
    `validate_tree()`).
  - Optionally, a `taxa_dm` argument takes a taxon dissimilarity instead
    of a tree. Decide this with the user; a tree is the common case.
  - Counts become relative abundances per sample, as ade4 weights them.
- **Result:** an `ecodive_ord` matrix of sample coordinates, columns
  `DPCoA1..`.
  - `eig` and `prop` hold the eigenvalues and inertia proportions of the
    sample cloud.
  - Attribute `taxa` holds the taxa coordinates (the matrix that phyloseq
    plots as "species"), plus `rao` (each sample's Rao diversity) and
    `trace`.
- **Wrapper:** joins `COUNT_ORDS` (it takes a count table, plus a tree via
  `...`), with aliases `dpcoa` and `doubleprincipalcoordinateanalysis`.
  - Check that the `dpco` prefix doesn't collide with `pco` (pcoa's alias)
    — it doesn't, since prefix matching runs from the start of the name.
- **Masking:** `ade4::dpcoa`, and phyloseq's `DPCoA` (a different case).


## Algorithm (from the paper; outline)

1. Patristic distances between taxa from the tree. ecodive's tree code
   (`src/ecotree.c`, used by UniFrac) can compute branch-path sums. Take
   Δ = √(patristic distance), or apply ade4's Euclidean check and
   correction; read Pavoine et al. §2.
2. A PCoA of Δ, **weighted by the taxa's overall relative abundances**
   (not equal weights), gives the taxa coordinates. This is a weighted
   PCoA: centre with the weights, then scale rows and columns by √w.
3. Sample coordinates are the abundance-weighted means of their taxa's
   coordinates. Sample weights are their share of the total, as ade4 uses.
4. A weighted PCA of the sample coordinates gives the DPCoA axes. Its
   total inertia is the among-sample part of Rao's apportionment.

For large numbers of taxa, avoid the dense taxa × taxa matrix where
possible: use block Lanczos (`src/lanczos.c`) with a matrix-free product
over the tree, as `pcoa()` avoids an n × n matrix.


## Fixtures and tests

- **Output-only (GPL / AGPL):** ade4 1.7-24 `dpcoa(df, dis)`, with
  `$li` (sample coordinates), `$dls` (taxa coordinates), `$eig` and
  `$RaoDiv`; and phyloseq `ordinate(physeq, "DPCoA")` as a cross-check.
  Write `data-raw/expected/dpcoa.r` using ex_counts and ex_tree, a random
  tree, and a non-ultrametric tree.
- Compare the coordinates (sign-matched), eigenvalues and Rao diversities
  at 1e-10.
- Rao's quadratic entropy per sample equals an independent R calculation
  in the test, computed from Δ.
- A taxon absent from every sample; a sample with one taxon; a tree with
  extra tips (pruned, as the UniFrac functions do).
- Results are identical for any `cpus`; a smoke-test line.


## Docs

- An Rd page; add it to the Ordination section of the reference index,
  next to the UniFrac-aware methods.
- A vignette note that DPCoA is the phylogeny-aware ordination of counts.
- NEWS.


## References

- Pavoine, S., Dufour, A.-B., & Chessel, D. (2004). From dissimilarities
  among species to dissimilarities among communities: a double principal
  coordinate analysis. *Journal of Theoretical Biology*, 228, 523–537.
- Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified
  approach. *Theoretical Population Biology*, 21, 24–43.
- Purdom, E. (2011). Analysis of a data matrix and a graph: metagenomic
  data and the phylogenetic tree. *Annals of Applied Statistics*, 5,
  2326–2358. (DPCoA for microbiome data.)


## Decisions made while implementing (2026-11-19)

Asked of the user at the start:

1. **`taxa_dm`: added** (recommended). `dpcoa(counts, tree = NULL,
   taxa_dm = NULL, margin = 1L, n_axes = 3L, cpus = n_cpus())`: exactly one
   of `tree` and `taxa_dm` (both is an error; an embedded tree is ignored
   when `taxa_dm` is given). `taxa_dm` is used as given, like ade4's `dis`,
   as a `dist` or symmetric matrix; labels are matched to the taxa (extra
   taxa dropped, missing ones an error), and without labels it must have
   one row per taxon, in the table's order. It must be Euclidean: its
   smallest PCoA eigenvalue may not fall below `-sqrt(eps) * trace` (the
   threshold `pcoa()` uses for rounding), found by Lanczos on the negated
   Gower product, which stops as soon as a Ritz value proves it; otherwise
   an error, as ade4 gives ("an Euclidean matrix is needed").
2. **Trees: accepted quietly** (recommended). Tree preparation reuses
   `philr_tree()` / `C_philr_tree()` (item 25). Polytomies, an unrooted
   tree and single-child nodes don't change path lengths, so DPCoA's
   results don't depend on how they are resolved: no warning (unlike
   `philr()`). Tips that aren't taxa are pruned quietly; a taxon missing
   from the tree is an error (philr_tree()'s message); missing (`NULL` or
   `NA` on a kept branch) branch lengths are an error. A root edge is
   ignored (it can't change a path length).
3. **Attributes: `taxa` and `rao`** (both recommended; the user didn't
   choose `rao_decomp` or the weights). After `method`, `eig`, `prop`:
   `taxa` (taxa × axes, ade4's `dls`, on the same `n_axes` axes as the
   samples; named by taxon), `rao` (ade4's `RaoDiv`, named by sample), and
   `trace` (the total inertia, ade4's between-samples diversity, which
   `prop` divides by; kept as `pca()` keeps `trace`). The taxa PCoA is
   never truncated (exact, as ade4 keeps every axis) and never formed.
4. **Weights as ade4; taxa absent from every sample kept** (recommended).
   Samples are weighted by their share of the grand total, taxa likewise;
   empty samples are an error (ade4: "Remove empty samples"; ecodive uses
   `cca_empty_msg()` with 'DPCoA'). A taxon absent from every sample, which
   ade4 refuses (`dudi.pco`: "Non convenient row.w (p=0)"), is kept with
   weight 0: it doesn't affect the axes, and its coordinates are its exact
   projection (tested against plain R; the others match ade4 run without
   it).

Made while implementing:

5. **Conventions (from ade4's outputs):** tree dissimilarities are
   √(patristic distance), as Pavoine et al. and phyloseq do (phyloseq's
   docs: square root of the cophenetic distance, Cailliez-corrected only if
   `is.euclid()` fails, which for √patristic can only be rounding). Rao is
   ½ Σ p_k p_l δ²_kl. `eig` is ade4's (the samples' w-weighted variance per
   axis, not divided by n − 1). Sample coordinates are ade4's `li`, taxa
   coordinates its `dls`. Axis signs follow `lz_sign_index()` on the sample
   coordinates, taxa flipped with them. ade4's `c1` (taxa scores on its own
   PCoA axes, basis-dependent), `RaoDis` and `RaoDecodiv` are not returned
   (the within-sample diversity is `weighted.mean(rao, w)`).
6. **The algorithm (derived for this item, MIT, `src/dpcoa.c`):** with
   Q = P − p̄1' (taxa × samples, columns summing to 0), the samples' Gram
   matrix is G = W^½ Q' K Q W^½ for any kernel K whose induced squared
   distances are δ²; the parts of K that aren't inner products cancel. For
   a tree, K_kl = the length the root's paths to k and l share, so K q is
   a clade sum up the tree times branch lengths, then a path sum down it:
   O(D), no D × D matrix and no taxa coordinates (equivalently, the taxa
   are embedded exactly in "branch space", x_ke = √L_e (1[k below e] −
   ā_e)). For `taxa_dm`, K is the Gower matrix −½ J δ² J, through
   `gower.c`'s `gower_mult()`. Block Lanczos (`lanczos()`, threaded) finds
   G's top k eigenpairs on the sample side, basis kept orthogonal to √w;
   each product is P'·(K·(P·u)) with a taxon-major copy of P for the first
   step and the sample-major one for the last, in fixed 64-row tasks, so
   results are identical for any `cpus` (tested 1-4). Rao per sample: over
   branches, L a(1 − a) (O(D) per sample, per-thread scratch); with
   `taxa_dm`, over pairs of the sample's non-zero taxa. Trace: tree,
   Σ_i w_i Σ_e L_e (a_ie − ā_e)² directly; `taxa_dm`, the pooled Rao of p̄
   less the weighted mean of `rao` (one pass over the packed distances).
7. **No variation:** an among-sample diversity at or below 1e-12 times the
   total (pooled) diversity is taken as 0, so samples with identical
   compositions give no axes (with `ord_keep()`'s warning) on every
   platform, rather than rounding-noise axes.
8. **Limits:** `n_axes` is cut to `min(n − 1, D − 1)`; one sample gives no
   axes but still `rao`; at least 2 taxa; non-negative finite values
   (proportions fine); taxa need names. Axes with eigenvalues below
   `sqrt(eps)` times the first are dropped by `ord_keep()` with its warning.
9. **Registry:** `COUNT_ORDS` row "Double Principal Coordinate Analysis",
   id `dpcoa`, aliases `doubleprincipalcoordinateanalysis` and
   `doubleprincipalcoordinatesanalysis`. The README snippet finds no
   collision in any registry; `'double'` prefix-matches only `dpcoa`.
   (A step-record exclusion for `taxa_dm` went with the step records,
   which were removed on 2026-11-20.)
10. **Masking:** documented in "Same names in other packages" of `?dpcoa`
    and `?count_ord` (`ade4::dpcoa()`).

11. **Speed** (4-core container, -O2): 300 samples × 1,000 taxa (5% non-zero,
    random tree): ecodive 0.02 s, ade4 252 s (same results to 1.4e-14).
    Full GlobalPatterns (26 samples × 19,216 taxa): 0.03 s. 5,000 samples ×
    20,000 taxa, sparse (2% non-zero): 0.50 s on 1 thread, 0.29 s on 4,
    identical.

Reference results (fixtures `data-raw/expected/dpcoa.r`, 37 KB):

- **ade4 1.7-24 `dpcoa()`**, nf = 4, on six cases: ex_counts / ex_tree; a
  60-tip `rtree()` (non-ultrametric) with 15 extra tips, shuffled taxa and
  an absent taxon; an unrooted `rcoal()` tree with a polytomy, proportions
  and a single-taxon sample; ade4's `humDNAm` (√ of its distances, through
  `taxa_dm`) and `ecomor` (`dist.taxo()`); and GlobalPatterns' 150 most
  abundant taxa. Sample and taxa coordinates (sign-matched), eigenvalues,
  Rao and the between-samples diversity agree to ≤ 1.3e-14 relative
  (tested at 1e-12); the tree cases agree through `taxa_dm` too, and with
  a plain-R `eigen()` implementation in the test.
- **phyloseq 1.46.0 `ordinate(method = "DPCoA")`** on the GlobalPatterns
  subset: the same agreement (it calls ade4; no correction was applied).
- Sparse (slam, dgCMatrix) and dense input agree to ~1e-15; `margin = 2`
  is identical.

**Clean-room note.** No ade4 or phyloseq source was read; only help pages
(`?ade4::dpcoa`, phyloseq's DPCoA documentation) and outputs. No GPL code
reached this session through error messages: the ones seen ("Remove empty
samples", "Non convenient row.w (p=0)", "an Euclidean matrix is needed")
print no code.


## Outcome (from the implementing session's final report)

**Item 26, done 2026-11-19.** Commits `6ed9799` (work), `e7ccd62` (tick).

- Tests: 16,495 at home, 312 in CRAN mode.
- ade4 1.7-24, six cases (trees and `taxa_dm`): coordinates, eigenvalues,
  Rao diversities ≤ 1.4e-14 relative. phyloseq 1.46.0 DPCoA: the same on
  GlobalPatterns.
- Speed: 300 × 1,000 in 0.02 s against ade4's 252 s; all 19,216
  GlobalPatterns taxa in 0.03 s. No taxa × taxa matrix is formed.
