# `dpcoa()`: double principal coordinate analysis

**Priority:** 20, the first of the clean-room methods. **Tier:** clean
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

- **Signature:** `dpcoa(counts, tree = NULL, margin = 1L, k = 3L, cpus = n_cpus())`.
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
