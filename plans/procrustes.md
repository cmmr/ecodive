# `procrustes()` and `procrustes_test()`: Procrustes rotation and its permutation test

**Priority:** 6. **Tier:** port from SciPy (BSD-3). Read
[README.md](README.md) first.

These are **two functions**: the rotation, and the test. vegan and QIIME 2
both keep them separate.


## What they are

- **Procrustes analysis** rotates, reflects and scales one ordination to
  fit another as closely as possible. Its residual sum of squares, m², says
  how similar two configurations of the same samples are: for example, the
  same samples under two beta diversity metrics, or 16S against
  metagenomics.
- **PROTEST** (Jackson 1995; Peres-Neto & Jackson 2001) tests whether the
  fit is better than chance, by permuting one configuration's rows.

QIIME 2 has `procrustes-analysis` and `partial-procrustes`; vegan has
`procrustes()` and `protest()`.


## Design

- **`procrustes(x, y, k = NULL, symmetric = FALSE, scale = TRUE)`.**
  - `x` is the target and `y` is rotated to fit it.
  - Each may be an ordination (`ecodive_ord`), a plain matrix, or a `dist`
    (then `pcoa()` is applied first).
  - Rows are matched by name when both are named, as elsewhere in
    ecodive. `k` axes are used; the default is
    `min(ncol(x), ncol(y))`.
  - Zero-pad the narrower configuration, as vegan and SciPy do.
  - `symmetric` and `scale` follow vegan's meanings, which are options of
    the one method, so they stay arguments.
- **Its result** (decide the class: a list of class `ecodive_procrustes`
  is likely):
  - the rotated `y` (an `ecodive_ord` with `method = 'procrustes'`);
  - `ss` (m²), `rotation`, `translation`, `scale`;
  - per-sample `residuals`;
  - `symmetric`.
- **`procrustes_test(x, y, k = NULL, permutations = 999L, strata = NULL, seed = 0L, cpus = n_cpus())`.**
  - Returns an `ecodive_test` of class `c('ecodive_procrustes_test', 'ecodive_test')`.
  - `statistic = c(r = sqrt(1 - ss))`, where `ss` is the **symmetric**
    m², as vegan computes it.
  - `p_value = (k + 1) / (n_perm + 1)`, where k is the number of permuted
    r ≥ observed, as vegan does.
  - The permuted values go in attribute `r_perm`.
  - Also return `ss` and the `procrustes()` result. It's a plain function,
    not part of `dist_test()`, since it takes two configurations.
- **Name:** the test is `procrustes_test()`, not vegan's `protest`, so that the `procrustes_*` family sorts together and "protest" isn't misread. The docs name it PROTEST (Jackson 1995).
- **Masking:** `vegan::procrustes`; vegan's `protest` has a different name, so it isn't masked.


## Reference implementations

- **Port from:** SciPy 1.18.1 (BSD-3,
  [LICENSE](https://github.com/scipy/scipy/blob/main/LICENSE.txt)).
  - `scipy.linalg.orthogonal_procrustes`: the SVD of XᵀY. Use it for the
    rotation in both modes.
  - `scipy.spatial.procrustes`: the **symmetric form only**. It centres
    and scales both configurations to unit norm, allows reflections, and
    gives disparity = m². vegan defaults to `symmetric = FALSE`, so
    implement the non-symmetric scaling on top of `orthogonal_procrustes`,
    following Gower & Dijksterhuis (2004) and vegan's documentation.
- **ecopy `procrustes_test` (MIT):** useful only as a cross-check. Its
  statistic is `1 − W²`, the same as vegan's symmetric ss. But its p-value
  is "proportion smaller" with no +1, and it takes no seed. **Use vegan's
  definitions instead.**
- **Output-only fixtures (GPL):** vegan `procrustes()` (both `symmetric`
  settings) and `protest()` with a supplied permutation matrix, so that
  p-values match exactly.
- **C:** the rotation is a small k × k SVD (k ≤ about 10). The permutation
  loop is the cost: parallelise it over permutations, each one a k × k SVD
  of a permuted cross-product, with a deterministic split.


## Tests

- Fixtures from vegan, written by `data-raw/expected/procrustes.r`.
  - Data: varespec PCoA against an NMDS of the same samples; two
    different-width configurations; a reflected case; an exact match
    (ss = 0).
  - Compare the rotated coordinates, rotation, translation, scale and ss
    at 1e-10; and `procrustes_test()`'s r, every permuted r, and the p-value
    (identical, given the shared permutation matrix).
- SciPy agreement: the symmetric disparity equals
  `scipy.spatial.procrustes`'s. Store it with a small Python generator if
  that helps, but the vegan fixtures are enough.
- `dist` input equals running `pcoa()` first.
- Mismatched samples give an error, or are matched by name.
- Results are identical for any `cpus`.
- A smoke-test line for each function.


## Docs

- Rd pages for both, cross-linked.
- Add them to the reference index, in a new "Comparing ordinations"
  subsection or beside `envfit()`.
- A short vignette example comparing Bray-Curtis and UniFrac PCoAs.
- NEWS.


## References

- Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
  40, 33–51.
- Gower, J. C., & Dijksterhuis, G. B. (2004). *Procrustes Problems*.
  Oxford University Press.
- Jackson, D. A. (1995). PROTEST: a PROcrustean randomization TEST of
  community environment concordance. *Écoscience*, 2, 297–303.
- Peres-Neto, P. R., & Jackson, D. A. (2001). How well do multivariate
  data sets match? *Oecologia*, 129, 169–178.
