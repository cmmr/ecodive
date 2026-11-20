# `largevis()`: LargeVis

**Priority:** 70 (last). **Tier:** port from uwot's C++ headers (BSD-2).
**Superseded** by UMAP. Read [README.md](README.md) first.


## Status

LargeVis (Tang et al. 2016) is a neighbour-embedding method that came
between t-SNE and UMAP. It builds a kNN graph with t-SNE-style
perplexity-calibrated weights, then optimises a layout by negative
sampling with a LargeVis-specific gradient. UMAP has largely replaced it,
and the original C++ code has been abandoned since 2016. It's last on the
list, but cheap: ecodive's `umap()` already follows the uwot headers that
contain LargeVis's gradient.


## Licensing

- **Port from:** uwot's C++ headers in `inst/include/uwot/`:
  `gradient.h` (`class largevis_gradient`), `perplexity.h` and `epoch.h`.
  Each is **BSD-2** ("Copyright 2020 James Melville"; checked 2026-11-18),
  **even though the uwot R package is GPL-3**.
  - Port only these headers. Don't read uwot's `R/` or `src/` wrappers,
    which are GPL-3; for the R-level defaults, use `lvish()`'s
    *documentation*.
  - `src/umap.c` already carries the uwot BSD-2 notice (see
    `inst/COPYRIGHTS`). Extend that stanza to cover the new files.
- The original lferry007/LargeVis is Apache-2.0, with its copyright line
  left blank. It uses GSL (GPL) for random numbers, so don't port it.


## Design

- **Signature:** `largevis(dm, n_axes = 2L, perplexity = 50, n_neighbors = NULL, n_epochs = NULL, learning_rate = 1, negative_sample_rate = 5, repulsion_strength = 7, init = 'pcoa', seed = 0L, cpus = n_cpus())`.
  - The defaults are those of uwot's `lvish()` documentation: perplexity
    50, repulsion_strength 7 (LargeVis's γ), negative_sample_rate 5, and
    n_neighbors 3 × perplexity. Verify them against the documentation.
  - Share `umap()`'s argument names where they mean the same thing
    (`n_epochs`, `learning_rate`, `negative_sample_rate`,
    `repulsion_strength`, `init`, `seed`). `repulsion_strength` is also
    `lvish()`'s name; don't call it `gamma`, which means something else in
    `phate()`.
- **Result:** an `ecodive_ord` matrix with columns `LargeVis1..`; `eig`
  and `prop` are `NA`. Attributes `perplexity`, `n_epochs`.
- **Wrapper:** joins `ORDINATIONS`, with aliases `largevis` and `lvish`.
- **Masking:** none (uwot's is `lvish`).


## Build

- Perplexity calibration: reuse `tsne()`'s binary search (`src/tsne.c`)
  or uwot's `perplexity.h`; use whichever ecodive's tsne already matches.
- Layout: `src/umap.c`'s batch SGD, with the LargeVis gradient swapped
  in. Keep the deterministic, thread-count-independent batching that
  `umap()` uses.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- The title comes from `ord_result()` automatically: pass
  `from = dv$title` (the input distance matrix's title, from
  `dist_vector()`), e.g. "LargeVis of Bray-Curtis Dissimilarity". The
  existing `plot.ecodive_ord()` uses it.
- The existing plot types need no extra attributes.
- Test: `display$title` for one call.

## Tests

- Affinities equal uwot's (an `lvish()` fixture, GPL, output-only) at
  1e-12.
- The SGD from the same start and seeds is bit-identical on Windows
  x86_64, as for `umap()`.
- Final quality (neighbours kept) is comparable to `lvish()`'s.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

- An Rd page that says UMAP is usually preferable.
- The method table in `vignettes/ordination.Rmd`, and NEWS.


## References

- Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016). Visualizing large-scale
  and high-dimensional data. *Proceedings of the 25th International
  Conference on World Wide Web*, 287–297.
