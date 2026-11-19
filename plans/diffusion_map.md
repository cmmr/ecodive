# `diffusion_map()`: diffusion maps

**Priority:** 30. **Tier:** port from pydiffmap (MIT) or scanpy (BSD-3).
Read [README.md](README.md) first.


## What it is

Diffusion maps (Coifman & Lafon 2006) build a Markov chain on a kernel
over the samples, and take the leading non-trivial eigenvectors of its
transition matrix, scaled by powers of their eigenvalues. Distances on the
map approximate "diffusion distances", which suit continuous gradients and
trajectories. Common in single-cell work (scanpy, destiny), and growing in
microbiome time series.


## Design

- **Signature:** `diffusion_map(dm, n_axes = 3L, n_neighbors = 15L, epsilon = 'bgh', alpha = 0.5, t = 1, cpus = n_cpus())`.
  - **Kernel:** Gaussian on `dm`, restricted to each sample's
    `n_neighbors` nearest neighbours.
  - **`epsilon`:** a number, or the Berry–Harlim–Giannakis automatic
    choice (`'bgh'`), as pydiffmap offers.
  - **`alpha`:** the density normalisation (Coifman–Lafon α; 0.5 is the
    Fokker–Planck form, 1 the Laplace–Beltrami form). `t` is the diffusion
    time.
  - **Pick one convention and document it.** scanpy's default kernel is
    UMAP's fuzzy graph; its `method = 'gauss'` uses Haghverdi et al.'s
    (2016) adaptive Gaussian; destiny is different again. Stating the
    convention matters more than matching all of them.
- **Result:** an `ecodive_ord` matrix with columns `DC1..`.
  - `eig` holds the eigenvalues of the transition matrix (excluding the
    trivial 1); `prop` is `NA`.
  - Attributes `epsilon`, `alpha`, `t`, `n_neighbors`.
- **Wrapper:** joins `ORDINATIONS`, with aliases `diffusionmap`, `diffmap`,
  `dmap`. Check them against existing names.
- **Masking:** none (destiny's is `DiffusionMap`).


## Reference implementations

- **Port from:** pydiffmap 0.2.0.1 (MIT), the textbook Coifman–Lafon form
  with α normalisation and the Berry–Harlim–Giannakis bandwidth choice.
  Unmaintained since 2021, but the maths is stable.
  - Alternatives: scanpy 1.13 `tl.diffmap` (BSD-3; the maths is in its
    `DPT` class); linfa-reduction `DiffusionMap` (MIT/Apache).
- **Not this:** scikit-learn's `SpectralEmbedding` is Laplacian
  eigenmaps, not a diffusion map.
- **Output-only fixture:** destiny `DiffusionMap` (GPL-3) is optional, as a
  quality comparison only, since its kernel differs.
- **Python fixture:** pydiffmap on a precomputed kernel
  (`data-raw/expected/diffusion_map.py`): eigenvalues and eigenvectors at
  1e-10, sign-matched.
- **ecodive code to reuse:**
  - `knn_from_dist()` (`src/knn.c`) for neighbours;
  - block Lanczos (`src/lanczos.c`) for the top eigenvectors of the
    symmetrised kernel D^{-1/2} K D^{-1/2}, which has the same
    eigenvalues as the transition matrix.


## Tests

- Matches pydiffmap at 1e-10 (sign-matched) on a swiss roll and a Bray
  table.
- Eigenvalues lie in (0, 1], with the trivial one dropped.
- A disconnected kNN graph gives an error or warning, as `isomap()` does.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

An Rd page, an entry in the method table of `vignettes/ordination.Rmd`, the
reference index, and NEWS.


## References

- Coifman, R. R., & Lafon, S. (2006). Diffusion maps. *Applied and
  Computational Harmonic Analysis*, 21, 5–30.
- Berry, T., Harlim, J., & Giannakis, D. (2015). Nonparametric forecasting
  of low-dimensional dynamical systems. *Physical Review E*, 91, 032915.
- Haghverdi, L., Büttner, M., Wolf, F. A., Buettner, F., & Theis, F. J.
  (2016). Diffusion pseudotime robustly reconstructs lineage branching.
  *Nature Methods*, 13, 845–848.
