# `principal_response_curves()`: principal response curves

**Priority:** 57. **Tier:** *derived* (a partial RDA). **Requires
[rda.md](rda.md) first.** Read [README.md](README.md) first.

**Name:** spelled out, not `prc`. `prc` is easily mistaken for
`prcomp()`, and it would prefix-match `pca`'s `prcomp` alias in
`count_ord()`.


## What it is

Principal response curves (van den Brink & ter Braak 1999) summarise how
treated communities diverge from a control over time. The method is a
partial RDA,

    Y ~ treatment:time + Condition(time)

with treatment contrasts against the control. The result is drawn as one
curve per treatment: its canonical coefficient at each time point. Each
taxon's weight shows how much it follows the curve. Used for experiments
sampled over time. vegan `prc()`.


## Design

- **Signature:** `principal_response_curves(counts, treatment, time, data = NULL, control = NULL, margin = 1L, norm = 'none', pseudocount = NULL, n_axes = 1L, cpus = n_cpus())`.
  - `n_axes` defaults to 1, the method's convention (one curve per
    treatment). Document it as the exception to the rule that ordinations
    with eigenvalues default to 3 axes
    ([parameter_names.md](parameter_names.md) §3.1).
  - `treatment` and `time` are factors, or column names of `data`.
  - `control` names the control level; the default is the first level,
    as in vegan.
- **Result:** a plain function, since the result isn't a sample × axis
  matrix. A list of class `ecodive_prc` with:
  - `coefficients`: a data frame with one row per treatment × time
    (control 0), and columns `treatment`, `time` and `PRC1..` (one per axis);
  - `weights`: the taxa × `n_axes` species weights;
  - `eig`, `prop`;
  - the underlying `rda()` result as `rda`.

  Plotting is out of scope; `coefficients` is ready for `ggplot`.
- **Masking:** none (vegan's is `prc`).


## Reference implementations

- **No permissive code:** build it on ecodive's partial `rda()`.
- **Output-only fixtures (GPL):** vegan 2.7-x `prc()` on its `pyrifos`
  example (`pyrifos`, with `week` and `dose`): `summary(prc)` coefficients
  and species weights.
- vegan's `prc` documentation explains the scaling: the coefficients are
  scaled relative to the species weights, so match `summary()`'s
  `scaling = 2` output. Check this numerically.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- The result is a list, not made by `ord_result()`, and plotting is out
  of scope for v1 (see Design). The `rda` element keeps the title its
  constructor gives it.
- If a `plot()` method is added later: its title is the method's name,
  "Principal Response Curves" (count tables have no title, so no
  `display` is needed), and the response-curve plot reads `coefficients`
  (one row per treatment × time, in time order).

## Tests

- The `pyrifos` fixture: coefficients and species weights at 1e-10, with
  each axis's sign matched.
- The result's `rda` equals
  `rda(counts, data = d, formula = ~ treatment:time + Condition(time))`.
- Edge cases: a single non-control treatment; unbalanced time points.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

An Rd page, an entry in the reference index under Ordination, and NEWS.


## References

- van den Brink, P. J., & ter Braak, C. J. F. (1999). Principal response
  curves: analysis of time-dependent multivariate responses of biological
  community to stress. *Environmental Toxicology and Chemistry*, 18,
  138–148.


## Decisions made while implementing

Item 57, 2026-11-25. No reference code was read or ported: vegan 2.7.8's
`prc()`, `summary()`, `anova()` and `permutest()` were only run, for
fixtures and to establish what they compute.

Asked of the user at the start, and answered:

1. **`plot()` and `print()` are included** (the plan had left plotting
   out). `plot.ecodive_prc()` draws the PRC diagram: each treatment's
   coefficients over time (the control's flat line at 0), with every
   taxon's weight as a tick on a scale at the right, the `species` (default
   10) largest labelled. Title "Principal Response Curves"; subtitle the
   axis's share of the variance within time points and, for axis 1 with
   tests, the first-axis p-value. `print()` gives the variance split, the
   axes, the tests, and the coefficients and largest weights of the first
   two axes.
2. **`scaling` argument, default `'symmetric'`** (vegan's scaling 3), with
   1 / 2 / 3, `'sites'` / `'species'` / `'symmetric'`, and `'form'` /
   `'covariance'` for 1 / 2 (item 44's `ord_scaling_arg()`; `'contribution'`
   is an error). Note: vegan 2.7.8's `plot.prc()` defaults to
   `"symmetric"`, but its `summary.prc()` to `"sites"` (scaling 1); the
   question to the user said both default to symmetric. All three are
   tested against `summary(prc, scaling =)`.
3. **Permutation tests of the first axis and of all axes**, as
   `anova(prc, first = TRUE)` and `anova(prc)`.
4. **Signs:** each axis is flipped so that its coefficient with the largest
   absolute value is positive (the first such, if tied); the weights and
   the stored `rda` element (scores, loadings, biplot scores, centroids)
   flip with it.

Made while implementing:

5. **Signature:** `principal_response_curves(counts, treatment, time, data
   = NULL, control = NULL, margin = 1L, norm = 'none', pseudocount = NULL,
   n_axes = 1L, scaling = 'symmetric', permutations = 999L, strata = NULL,
   plots = NULL, seed = 0L, cpus = n_cpus())`. No `scale` argument (the
   plan has none; vegan's `correlation` and `const` are not offered).
   `treatment`, `time`, `strata` and `plots` follow `per_sample()`'s rule.
   `time` becomes a factor (a factor's levels keep their order; numbers sort
   numerically). `control` defaults to the first level, as in vegan; any
   other level may be named, and the coefficients' `treatment` keeps its
   own level order.
6. **What vegan computes, found from its outputs:** the coefficients are
   each treatment × time cell's mean site score less the control's at that
   time (the "lc" and "wa" scores give the same cell means, to 2e-14), in
   the requested scaling; the weights are `scores(display = "species",
   scaling =)`. The model is `rda(Y ~ treatment:time + Condition(time))`:
   ecodive calls `rda()` with that formula (or `~ treatment` with one time
   point) and quietly drops its aliasing message; the scalings use item
   44's per-axis factors \eqn{\sqrt{\lambda/\mathrm{trace}}}.
7. **The first-axis test** (from `permutest(first = TRUE)`'s `num` and
   `den`): F = λ₁ / (RSS / ν), where for each permutation λ₁ is the first
   constrained eigenvalue and RSS the residual of the full permuted model
   (not total − λ₁). Permutations follow `permanova()`'s scheme (scatter
   the design's rows, `v[p[a]] = x[a]`, i.e. vegan's `Y[p, ]`; make them
   orthogonal to the conditions; re-orthonormalize), reproduced by hand in
   R for whole-plot, free and within-week permutations before coding.
8. **C kernel** (`src/principal_response_curves.c`, `C_prc_perm()`): for
   each permutation, the orthonormal permuted design Q (n × m), then L = A'Q
   by `pca.c`'s `pca_cross()` (sparse, threaded, fixed summation order),
   then C = L'L (m × m), its trace and its largest eigenvalue
   (`lz_sym_eigen()`). Permutations go in batches of about 256 design
   columns (capped at 2²⁴ doubles for L); designs and eigenproblems are
   one permutation per thread. Results are identical for any `cpus`
   (tested at 1, 2 and 4). Both tests come from one pass. The conditions'
   sum of squares is ‖A'Qz‖².
9. **`plots`** (whole-plot permutation, as permute's `Plots(type =
   "free")` with `Within(type = "none")`): every plot must have exactly one
   sample at each time point; each permutation of the plots
   (`C_perm_indices()`, PCG, seeded) moves each plot's sample at time t to
   another plot's sample at time t. Only one of `strata`, `plots` and a
   `permutations` matrix may be given. `strata` permutes freely within
   blocks (`perm_source_*`), as in `permanova()`.
10. **The overall test equals `permanova()` of the Euclidean distances
    only with vegan's coding of the constraints:** the indicators of the
    non-control cells, `model.matrix(~ time + time:treatment)` after time's
    columns (tested: same F, p and permuted F). `formula = ~
    treatment:time + Condition(time)` spans the same model, so has the same
    F, but its aliased model matrix keeps other columns, and permuting
    those gives other permuted F values. Documented.
11. **Times with no control samples:** their coefficients are `NA`, with a
    warning. vegan drops some other treatment's column there and reports
    coefficients relative to that treatment (an artefact of aliasing), so
    those coefficients, and possibly the permuted statistics (which depend
    on which aliased column is kept), differ from vegan's. Cells with no
    samples give `NA` (as in vegan).
12. **Result:** a list of class `ecodive_prc`: `coefficients` (long data
    frame: `treatment`, `time`, `PRC1..`; control first, then the other
    treatments in level order, time in level order within each; control
    rows 0), `weights` (taxa × axes), `eig`, `prop` (as `rda()`'s: of the
    variance not explained by time), `test` (rows `PRC1` and `all`;
    columns `df`, `variance`, `F`, `z`, `p_value`), `df_residual`,
    `residual`, `F_perm` (n_perm × 2), `n_perm`, `permutation` (`'free'`,
    `'within strata'`, `'whole plots'`, `'supplied matrix'`), `scaling`,
    `control`, `rda`. No `display`: count tables have no title, and the
    `rda` element keeps its constructor's title.
13. **Edge cases:** `n_axes` above the number of constrained axes is
    reduced with a warning; no constrained variation is an error; a
    saturated model (one sample per cell) gives `NA` F and p-values; fewer
    than 3 samples, a single treatment level, an unknown `control`, NAs and
    conflicting permutation options are errors.
14. **Fixtures:** `data-raw/expected/principal_response_curves.r` →
    `tests/tinytests/expected/principal_response_curves.rds` (108 kB, xz):
    vegan's `pyrifos` with 9 cases (whole-plot, free and within-week
    `shuffleSet()` permutations; Hellinger; CLR of back-transformed counts
    with pseudocount 1; a control that is not the first level; two doses;
    an empty cell; 20 samples dropped), each in scalings 1, 2 and 3 for
    axes 1 and 2. vegan 2.7.8 needs `shuffleSet()` objects; the fixture
    stores plain matrices.
15. **Docs:** `man/principal_response_curves.Rd`, `man/plot.ecodive_prc.Rd`;
    `rda()`'s note on the methods built on it; NEWS (2.4.0); pkgdown
    (Ordination, after `forward_select`; `plot.ecodive_prc` under Plotting);
    `parameter_names.md` (`treatment`, `time`, `control`, `plots`, and
    `scaling`). README and vignette tables were left alone (the README's
    table does not list the newer methods either).

Reference results:

- **vegan 2.7.8** `prc()` / `summary()`, 9 cases on `pyrifos` (see 14),
  scalings 1, 2 and 3, axes 1 and 2: coefficients and species weights
  within 7e-15 (relative to the largest), eigenvalues within 5e-15; `NA`
  cells in the same places.
- **`anova(prc, first = TRUE)` and `anova(prc)`**, with `shuffleSet()`
  permutations (whole ditches, free, within weeks; plain, Hellinger, CLR,
  another control): F within 4e-14, every permuted F within 1e-13
  (relative), p-values identical.
- **ecodive's own methods:** the `rda` element equals `rda(counts, data =
  d, formula = ~ treatment:time + Condition(time))` up to each axis's sign
  (1e-12); the overall test equals `permanova()` of the Euclidean
  distances with vegan's coding of the constraints (F, every permuted F,
  p); with one time point, `eig` equals `rda(counts, treatment)`'s and the
  test `permanova()`'s; strata and whole-plot permutations equal the same
  permutations given as a matrix.
- **Speed:** 480 samples (40 plots × 12 times, 5 doses) × 2,000 taxa
  (14% non-zero), both tests with 999 whole-plot permutations: 3.1 s on 4
  threads; `vegan::prc()` with `anova(first = TRUE)` and `anova()`: 697 s.
  Eigenvalues equal to 2e-16.
- **Tests:** 31,598 at home and 456 in CRAN mode (from 30,895 and 450).
  `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed size, as
  before).


## Outcome (from the implementing session's final report)

**Item 57, done 2026-11-26.** Commits `e4ce79a` (work) and the tick.

- `principal_response_curves()`: a partial `rda()` with coefficients and
  weights in vegan's scalings 1, 2, 3 (default `'symmetric'`), permutation
  tests of the first axis and all axes (`plots` for whole-plot,
  `strata` for within-block permutations), `plot()` (the PRC diagram) and
  `print()`.
- Tests: 31,598 at home and 456 in CRAN mode (from 30,895 and 450);
  `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed size).
- vegan 2.7.8: coefficients and weights within 7e-15, eigenvalues 5e-15,
  F 4e-14, every permuted F 1e-13; p-values identical.
- Speed: 480 × 2,000, both tests with 999 whole-plot permutations, 3.1 s
  on 4 threads against vegan's 697 s.
- No code read or ported; no GPL source text seen.
- Open: the overall test equals `permanova()` only with vegan's coding of
  the constraints (decision 10); at a time point with no control samples,
  coefficients are `NA` where vegan reports aliasing artefacts.
