# `phate()`: Potential of Heat-diffusion for Affinity-based Trajectory Embedding

**Priority:** 43. **Tier:** clean room, with a **license warning on the
reference itself**. Read [README.md](README.md) first.


## Licensing (read first)

- PHATE's repository LICENSE.md is now the **Yale "Non-Commercial
  License"** (checked 2026-11-18 at
  [LICENSE.md](https://raw.githubusercontent.com/KrishnaswamyLab/PHATE/master/LICENSE.md)).
  It says "Any commercial use of the Software requires a separate
  commercial license from the copyright holder at Yale University", and
  requires derivatives to be shared on GitHub.
  - PyPI `phate` 2.0.0 still says GPL-2, and the old source headers say
    GPLv2.
  - Its dependency graphtools is **GPL-3**, and phateR is GPL-2.
- **No PHATE code may be read for porting.** Implement strictly from Moon
  et al. (2019) and the method's public documentation.
- **Fixtures:** before running phate to generate any, read the
  non-commercial terms and confirm with the user that generating fixtures
  for an MIT package is acceptable. If not, validate only on quality
  measures and the paper's own worked examples.


## What it is

PHATE (Moon et al. 2019) has five steps:

1. Build an adaptive α-decay kernel over the k nearest neighbours.
2. Row-normalise it into a diffusion operator P.
3. Raise P to a time t, chosen automatically by von Neumann entropy (the
   "knee").
4. Take potential distances: the distances between −log(P^t) rows.
5. Embed them with metric MDS (SMACOF).

It's used in single-cell work to show trajectories, with some microbiome
use for time series and gradients.


## Design

- **Signature:** `phate(dm, n_axes = 2L, n_neighbors = 5L, decay = 40, t = 'auto', gamma = 1, n_landmarks = NULL, seed = 0L, cpus = n_cpus())`.
  - The defaults are PHATE's documented ones (knn = 5, decay = 40,
    t = "auto", gamma = 1). Verify them against its documentation, not its
    code.
  - `gamma` is the potential transform: 1 is log, 0 is √.
  - `n_landmarks` (PHATE's speed-up for large n) is out of scope for v1.
- **Result:** an `ecodive_ord` matrix with columns `PHATE1..`; `eig` and
  `prop` are `NA`. Attributes `t` (as used), `decay`, `n_neighbors`,
  `stress` (from the MDS step).
- **Wrapper:** joins `ORDINATIONS`, with alias `phate`.
- **Masking:** `phateR::phate`.


## Build from ecodive's own code

- Neighbours: `knn_from_dist()`.
- Diffusion operator powers: sparse or dense, depending on n. For n up to
  a few thousand, dense P^t by repeated squaring is fine; parallelise the
  matrix products with a deterministic split.
- Von Neumann entropy for automatic t: needs P's eigenvalues, from
  `src/lanczos.c` or a dense symmetric solver.
- The final embedding: ecodive's `mds()` engine (the SMACOF port in
  `src/nmds.c`), started from classical MDS as PHATE does.


## Tests

- Quality measures, not coordinates:
  - neighbours kept;
  - correlation of the embedding's distances with the potential
    distances;
  - trajectory recovery on a synthetic branching dataset (as in the
    paper's examples).
- If fixtures are approved (see Licensing): compare the potential
  distances at 1e-8 for a fixed t, and the automatic t exactly.
- Results are identical for any `cpus`, and the same `seed` gives the same
  result; a smoke-test line.


## Docs

- An Rd page that notes the method's origin, and that ecodive's
  implementation is independent.
- The method table in `vignettes/ordination.Rmd`, and NEWS.


## References

- Moon, K. R., van Dijk, D., Wang, Z., Gigante, S., Burkhardt, D. B.,
  Chen, W. S., Yim, K., van den Elzen, A., Hirn, M. J., Coifman, R. R.,
  Ivanova, N. B., Wolf, G., & Krishnaswamy, S. (2019). Visualizing
  structure and transitions in high-dimensional biological data. *Nature
  Biotechnology*, 37, 1482–1492.
