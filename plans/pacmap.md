# `pacmap()`: Pairwise Controlled Manifold Approximation

**Priority:** 29. **Tier:** port from PaCMAP (Apache-2.0). Read
[README.md](README.md) first.


## What it is

PaCMAP (Wang, Huang, Rudin & Shaposhnik 2021) is a neighbour-embedding
method like t-SNE and UMAP. It optimises three kinds of pairs (near
neighbours, mid-near pairs and further pairs), with weights that change
over the course of the fit. That aims to keep both local and global
structure. Growing use in single-cell work and, increasingly, microbiome
work.


## Design

- **Signature:** `pacmap(dm, n_axes = 2L, n_neighbors = 10L, mn_ratio = 0.5, fp_ratio = 2, n_iters = 450L, init = 'pcoa', learning_rate = 1, seed = 0L, cpus = n_cpus())`.
  - The argument names follow PaCMAP's (`n_neighbors`, `MN_ratio`,
    `FP_ratio`, `num_iters`), mapped to ecodive's conventions:
    `n_iters` is an exact count, as for tsne.
  - `init`: `'pcoa'` (ecodive's equivalent of PaCMAP's PCA start), or a
    matrix.
- **Result:** an `ecodive_ord` matrix with columns `PaCMAP1..`; `eig` and
  `prop` are `NA`. Attributes `n_neighbors`, `mn_ratio`, `fp_ratio`,
  `n_iters`.
- **Wrapper:** joins `ORDINATIONS`, with aliases `pacmap` and
  `pairwisecontrolledmanifoldapproximation`.
- **Licensing:** PaCMAP's LICENSE file is **Apache-2.0** (© 2020 Wang,
  Huang, Rudin, Shaposhnik; checked). PyPI's classifiers also say MIT, but
  the LICENSE file governs.
  - Ported files keep the Apache header and notices, as
    `src/pcg_basic.c` does.
  - Add a `Files:` stanza with the full Apache text to `inst/COPYRIGHTS`,
    and the authors as `cph` in `DESCRIPTION`.
- **Masking:** none in R.


## Reference implementations

- **Port from:** PaCMAP 0.9.1 (Python with numba; maintained),
  [LICENSE](https://raw.githubusercontent.com/YingfanWang/PaCMAP/master/LICENSE).
  - It uses Annoy for approximate neighbours. ecodive takes `dm`, so use
    exact neighbours from `knn_from_dist()` and **document that
    difference**.
  - It samples mid-near and further pairs at random. Use PCG with `seed`,
    and generate the pairs per sample from (seed, sample index), so that
    results don't depend on the thread count.
- **Python fixture:** PaCMAP on `dm` (it accepts precomputed distances via
  `distance = 'precomputed'` — verify this), written by
  `data-raw/expected/pacmap.py`.
  - Compare the pair-selection rules and one optimisation step from a
    fixed start and fixed pairs, exactly.
  - Compare final layouts by quality only: neighbours kept and global
    structure (seed-equivalent differences are fine).
- **ecodive code to reuse:** `src/umap.c`'s SGD layout loop and the
  thread pool. PaCMAP updates all points each iteration (full-batch Adam),
  which parallelises deterministically.


## Tests

- Pair selection and one step match the Python fixture exactly.
- Quality on the GEMS-like and swiss-roll data matches PaCMAP's within
  tolerance.
- Results are identical for any `cpus`, and the same `seed` gives the same
  result; a smoke-test line.


## Docs

- An Rd page, the method table in `vignettes/ordination.Rmd`, the
  reference index, and NEWS.
- The `inst/COPYRIGHTS` stanza (Apache).


## References

- Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). Understanding
  how dimension reduction tools work: an empirical approach to deciphering
  t-SNE, UMAP, TriMap, and PaCMAP for data visualization. *Journal of
  Machine Learning Research*, 22(201), 1–73.
