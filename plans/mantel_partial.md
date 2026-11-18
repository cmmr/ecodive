# `mantel_partial()`: split from `mantel(given =)`

**Priority:** 1 (housekeeping). **Tier:** refactor of existing, tested
code, so no new reference is needed. Read [README.md](README.md) first.


## Why

ecodive gives each named method variant its own function, as the UniFrac
family does. The partial Mantel test (Smouse, Long & Sokal 1986) is a
distinct, named test: vegan has `mantel.partial()` separate from
`mantel()`. Today ecodive runs it through `mantel(dm, y, given = z)`.
`mantel()` is new in the unreleased 2.4.0, so splitting it costs no
released API.


## Design

- **`mantel(dm, y, data = NULL, method = 'pearson', permutations = 999L, strata = NULL, alternative = 'greater', seed = 0L, cpus = n_cpus())`**:
  the same as now, without `given`.
- **`mantel_partial(dm, y, given, data = NULL, method = 'pearson', permutations = 999L, strata = NULL, alternative = 'greater', seed = 0L, cpus = n_cpus())`**:
  `given` is required, and accepts what `given` accepts now (a dist, a
  matrix, a vector, or a column name of `data`).
- **Shared code:** move the body into an internal `mantel_core(..., given)`
  that both functions call. The C code (`src/mantel.c`) is unchanged.
- **Result:** keep the class `c('ecodive_mantel', 'ecodive_test')` and the
  `partial` field, so printing and `as.data.frame()` work for both. Make
  the title read "Partial Mantel test" when partial.
- **`dist_test()`:** add a `DIST_TESTS` row:
  `Partial Mantel Test, mantel_partial, numeric, partialmantel mantelpartial`.
  `dist_test()` must map `groups` to `y` for it, as it does for `mantel`.
  Edit the `y <- if (test$id == 'mantel') ...` line in `R/stat_tests.r` to
  cover both, and `given` passes through `...`.
- **Old calls:** `mantel(..., given = z)` is an error pointing to
  `mantel_partial()`. Because 2.4.0 is unreleased, no deprecation period is
  needed; check with the user if 2.4.0 has been released by then.
- **Masking:** `vegan::mantel.partial` doesn't clash with
  `mantel_partial`.


## Tests

- Move the partial cases in `tests/tinytests/test-mantel.r` (those using
  `ex$z` / `given =`) to `mantel_partial()`.
- Keep the vegan `mantel.partial()` fixture comparisons, which cover
  Pearson, Spearman, Kendall, ties and `y | z`, at their current
  tolerances (statistic and `r_perm` 1e-12, p-values identical).
- New checks:
  - `mantel(..., given = )` gives the pointer error;
  - `mantel_partial()` without `given` is an error;
  - `dist_test(dm, 'mantel_partial', 'ph', meta, given = 'lat')` equals
    `mantel_partial()`, and `match_dist_test('mantel.partial')` finds it;
  - results are identical for `cpus` 1 and 4;
  - a smoke-test line.


## Docs

- Give `mantel_partial()` its own Rd file, and remove `given` from
  `mantel()`'s docs.
- Update `vignettes/stat_tests.Rmd` if it mentions `given`.
- Add `mantel_partial` to the "Hypothesis Tests" section of
  `pkgdown/_pkgdown.yml`.
- Reword the existing `mantel()` NEWS bullet in place, since 2.4.0 is
  unreleased.


## Reference

Smouse, P. E., Long, J. C., & Sokal, R. R. (1986). Multiple regression and
correlation extensions of the Mantel test of matrix correspondence.
*Systematic Zoology*, 35, 627–632.


## Decisions made while implementing (2026-11-18)

- **2.4.0 is unreleased** (CRAN has 2.3.1), so there is no deprecation period.
- **`mantel(given =)` gives R's own error** ("unused argument (given = ...)"),
  not a custom one pointing to `mantel_partial()`. The user chose this over a
  trailing `...` or a hidden formal argument, to keep `mantel()`'s signature
  clean (as the `k` -> `n_axes` rename did). The test checks for "unused
  argument".
- **`as.data.frame()` names the row `'mantel_partial'`** for a partial test
  (a small `as.data.frame.ecodive_mantel()` method), so that rows bound
  together from `dist_test()` carry each test's `dist_test()` id. The class
  stays `c('ecodive_mantel', 'ecodive_test')`.
- **`mantel_partial()` without `given`** (missing or `NULL`) is an error that
  says to use `mantel()` for a simple test; `dist_test(dm,
  'mantel_partial', ...)` without `given` gives the same error.
- **Tests:** the partial cases moved to a new
  `tests/tinytests/test-mantel_partial.r`, which shares the
  `expected/mantel.rds` fixture with `test-mantel.r`.
- **Reference check:** besides the stored fixtures, 54 fresh cases (n = 12
  to 60; Pearson, Spearman, Kendall; tied and untied distances) against
  live `vegan::mantel.partial()` 2.6-4 with the same permutation matrices:
  statistics within 7e-14 (relative), permuted statistics within 4e-16, and
  every p-value identical. Re-running `data-raw/expected/mantel.r` with
  vegan 2.6-4 reproduces the stored fixture (made with 2.7.6) to 1e-12.
- **Docs:** `vignettes/stat_tests.Rmd` gains a `mantel_partial()` row and a
  short example; `vignettes/ecodive.Rmd`, `README.md` and the draft
  benchmark article list it too.
