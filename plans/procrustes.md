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

- **`procrustes(ord, target, axes = NULL, symmetric = FALSE, dilate = TRUE)`.**
  - `ord` is rotated (and reflected, and optionally dilated) to fit
    `target`. vegan's `procrustes(X, Y)` rotates `Y` onto `X`, so
    ecodive's `procrustes(ord, target)` corresponds to vegan's
    `procrustes(target, ord)`; say so in the docs.
  - Each may be an ordination (`ecodive_ord`), a plain matrix, or a `dist`
    (then `pcoa()` is applied first).
  - Rows are matched by name when both are named, as elsewhere in
    ecodive.
  - `axes` selects **which** axes (column indices) of each input are used,
    as in `envfit()`. The default `NULL` uses the first
    `min(ncol(ord), ncol(target))`.
  - If one input has fewer of the requested axes, zero-pad it, as vegan
    and SciPy do.
  - `symmetric` follows vegan's meaning. `dilate` is vegan's `scale`:
    whether to allow a uniform dilation of `ord` (Gower's term). It's
    renamed because `scale` means "standardise features" everywhere else
    in ecodive. Both are options of the one method, so they stay
    arguments.
- **Its result** (decide the class: a list of class `ecodive_procrustes`
  is likely):
  - the rotated `ord` (an `ecodive_ord` with `method = 'procrustes'`);
  - `ss` (m²), `rotation`, `translation`, `dilation` (the dilation factor, vegan's `scale`;
    1 when `dilate = FALSE`);
  - per-sample `residuals`;
  - `symmetric`.
- **`procrustes_test(ord, target, axes = NULL, permutations = 999L, strata = NULL, seed = 0L, cpus = n_cpus())`.**
  - `ord`, `target` and `axes` as in `procrustes()`.
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


## Decisions made while implementing

Asked of the user at the start (2026-11-18):

- **Result class (user):** `procrustes()` returns the rotated `ord` itself,
  an `ecodive_ord` matrix with `method = 'procrustes'` (eig and prop `NA`),
  not a list. The rest are attributes, in this order: `ss`, `rotation`,
  `translation`, `dilation`, `residuals`, `symmetric`, and `target` (the
  target as fitted: centred, scaled if symmetric, zero-padded; vegan's `X`).
  `print.ecodive_ord()` then works unchanged.
- **dist inputs (user):** a new `n_axes` argument (default `3L`, as
  `pcoa()`'s) is passed to `pcoa()` when `ord` or `target` is a `dist`, so
  a dist input gives exactly `procrustes(pcoa(d, n_axes), ...)`. Ignored for
  matrix inputs. `procrustes_test()` takes it too.
- **Mismatched samples (user):** when both inputs have row names, `target`
  is reordered to `ord`'s samples; different sample sets are an error naming
  the missing samples (as `mantel()` does). Without names on both, the rows
  must be in the same order and number.
- **Default `axes` (user):** as the plan says, the first
  `min(ncol(ord), ncol(target))` axes. vegan pads the narrower input with
  zeros; documented as a difference, and `axes = 1:k` reproduces vegan.

Made while implementing:

- `axes` takes whole numbers only (not names, since the two inputs'
  axis names usually differ), up to the wider input's column count; axes
  past an input's width are zero columns. No warning (vegan warns).
- The rotated result's column names come from `target` (it is in target's
  frame), or from `ord` for axes where `target` is zero-padded; unnamed
  columns are `Dim1`, `Dim2`, ...
- `rotation` is always k × k after padding. When `ord` is the narrower
  input, vegan returns a (ncol(ord) × k) rotation; its rows equal the first
  rows of ecodive's (the padded rows are arbitrary). Tests compare those
  rows.
- `translation` follows vegan: `xmean - dilation * ymean %*% rotation`, so
  `dilation * ord %*% rotation + translation` gives the fitted coordinates.
  With `symmetric = TRUE` it is exactly 0, since both inputs are centred
  before scaling (vegan's is ~1e-18 there).
- Inputs with no variation (all samples at one point) are an error, as are
  fewer than 2 samples, NA/infinite coordinates, and non-numeric inputs.
- The observed fit is computed in R with `svd()` (SciPy's
  `orthogonal_procrustes()`: X'Y = U S V', rotation V U', dilation
  tr(S) / tr(Y'Y)). The permutation loop is C (`src/procrustes.c`,
  `C_protest`): per permutation the k × k cross-product X' Y[p, ] and its
  nuclear norm by one-sided Jacobi, one permutation per thread, so results
  are identical for any `cpus`. Permutation p takes `ord`'s row i to be
  `ord[p[i], ]`, as vegan's `Y[take, ]` (vegan permutes its `Y`, which is
  ecodive's `ord`).
- `procrustes_test()` result fields: `statistic` (`r = sqrt(1 - ss)`),
  `p_value`, `z`, `n_perm`, `ss`, `n`, `n_axes` (axes compared, after
  padding), `procrustes` (the symmetric fit), `permutation`; attribute
  `r_perm`. Title "Procrustes randomization test (PROTEST)". No
  `alternative` (vegan's test is one-sided); ties within
  `sqrt(.Machine$double.eps)` count, which reproduces vegan's p-values.
- `strata` is a vector with one value per sample (matched by name), since
  there is no `data` argument.
- Fixtures: `data-raw/expected/procrustes.r` (vegan 2.7.8, permute 0.10.0)
  and `procrustes.py` (SciPy 1.18.1, symmetric disparity). Cases: varespec
  PCoA (3 axes) vs NMDS (2 axes), padded and same-width; a reflected,
  rotated, dilated, shifted copy (ss = 0); an exact match; unrelated noise
  (p = 0.66), each with symmetric and dilate both ways.
- Docs: the two functions sit beside `envfit()` in the pkgdown Ordination
  section; a "Comparing two ordinations" section in
  `vignette('ordination')` compares Bray-Curtis and unweighted UniFrac PCoAs
  of 30 simulated samples of `ex_tree`'s taxa (`ex_counts` has only 4
  samples); a "Comparing ordinations" row in the README table.
- Credit: "The SciPy developers" added to `Authors@R` (ctb, cph), with a
  `Files: src/procrustes.c R/procrustes.r` stanza in `inst/COPYRIGHTS`.


## Outcome (from the implementing session's final report)

**Item 6, done 2026-11-18.** Commits `a1c2a66` (work), `0756364` (tick).

- Tests: 7,186 at home, 194 in CRAN mode.
- vegan 2.7.8 `procrustes()`, 20 cases: rotated coordinates, rotation and
  residuals about 5e-16 relative, translation 1e-17, dilation exact, m²
  about 2e-15.
- vegan `protest()` with a shared 199-row matrix: r and every permuted r
  ≤ 2.8e-16; p-values identical.
- SciPy 1.18.1 `scipy.spatial.procrustes`: symmetric disparity ≤ 2.8e-17.
- Speed: 2,000 samples, 999 permutations: 0.03 s against vegan's 0.25 s.
- Default `axes` is the first `min(ncol)` (vegan zero-pads; `axes = 1:k`
  reproduces it). SciPy credited (`Authors@R`, `inst/COPYRIGHTS`).
