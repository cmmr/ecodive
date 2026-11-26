# `mantel_regression()`: multiple regression on distance matrices (MRM)

**Priority:** 69. **Tier:** clean room. Write it from Legendre, Lapointe &
Casgrain (1994) and Lichstein (2007). ecodist `MRM()` (GPL) is
output-only. Read [README.md](README.md) first.

**Name:** `mantel_regression()` puts it in the `mantel_*` family (Smouse
et al. 1986 introduced it as a "multiple regression extension of the
Mantel test"). It also avoids `mrm`, which is easily mistaken for
`mrpp()`. Give it the aliases `mrm` and `multipleregressionondistancematrices`
in any registry it joins (it is a plain function at first).


## What it is

MRM regresses the community distances, unfolded into a vector of pairs, on
one or more explanatory distance matrices (geography, environment, time).
Significance comes from permuting the samples of the response matrix,
since the pairs aren't independent. It's used to model distance decay and
to separate the effects of environment and space.


## Design

- **Signature:** `mantel_regression(dm, ..., data = NULL, formula = NULL, method = 'pearson', permutations = 999L, seed = 0L, cpus = n_cpus())`.
  - Explanatory distances come through `...` (named dist objects or
    matrices, matched by sample name), or through `formula` with `data`.
    In a formula, numeric columns become Euclidean differences, as
    `mantel()`'s `y` does.
  - `method` uses the same vocabulary as the rest of the Mantel family.
    `'pearson'` uses ordinary least squares on the pair vectors (ecodist's
    `mrank = FALSE`). `'spearman'` ranks every vector first, the rank
    regression that is the Spearman analogue (ecodist's `mrank = TRUE`).
    Document this mapping.
- **Result:** a list of class `ecodive_mantel_regression`:
  - `coefficients`: a data frame with `term`, `estimate`, `p_value`
    (permutation);
  - `r_squared`, `F`, and the overall model `p_value`;
  - `n_perm`;
  - attribute `f_perm`.

  It is a plain function, not part of `dist_test()` (it can take several
  matrices).
- **Masking:** none (ecodist's is `MRM`).


## Algorithm (from the papers)

- Unfold each lower triangle into vectors of length n(n − 1)/2. Use
  `dist_vector()` in `R/ord_utils.r`, or work directly on packed triangles
  in C.
- Fit OLS once: X is small (number of pairs × p), so form XᵀX and Xᵀy.
- For each permutation, permute the rows and columns of the **response**
  matrix only, refit the slopes and R², and count how many permuted
  pseudo-F or |t| are ≥ the observed ones.
  - Each permutation costs O(pairs × p). Parallelise over permutations,
    with a deterministic split, using `perm_source_*` from
    `src/permute.c`.
  - Keep sums in a fixed order.
- Check the p-value definition used for coefficients (two-sided on t, or
  on |b|) against ecodist's documented output.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- The result is a list of class `ecodive_mantel_regression`, not an
  `ecodive_test`, so no constructor titles it. Set `display$title` with
  `set_display()` and `title_of()` from the response's title: "Multiple
  Regression on Distance Matrices of Bray-Curtis Dissimilarity"; without
  one, just the method's name.
- Coefficient tables and plots name the predictors by the `term` column:
  the names given in `...` or the formula's terms.
- `plot()` draws the coefficients with their permutation p-values, read
  from `coefficients`.

## Fixtures and tests

- **Output-only (GPL):** ecodist 2.1.3
  `MRM(dist(y) ~ dist(x1) + dist(x2), nperm = ...)`, with `$coef`, `$r.squared`
  and `$F.test`. ecodist uses R's RNG, so compare statistics exactly, and
  p-values either via a supplied permutation order (if achievable) or
  statistically. Write `data-raw/expected/mantel_regression.r`.
- With one explanatory matrix, `r_squared` equals `mantel()`'s r², and
  the slope equals an OLS fit on the vectors.
- Coefficients and R² equal `lm()` on the unfolded vectors at 1e-12,
  computed live in the test.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

An Rd page next to `mantel()` in the reference index, and NEWS.


## References

- Smouse, P. E., Long, J. C., & Sokal, R. R. (1986). Multiple regression
  and correlation extensions of the Mantel test of matrix correspondence.
  *Systematic Zoology*, 35, 627–632.
- Legendre, P., Lapointe, F.-J., & Casgrain, P. (1994). Modeling brain
  evolution from behavior: a permutational regression approach.
  *Evolution*, 48, 1487–1499.
- Lichstein, J. W. (2007). Multiple regression on distance matrices: a
  multivariate spatial analysis tool. *Plant Ecology*, 188, 117–131.


## Decisions made while implementing (2026-11-26)

ecodist 2.1.3's `MRM()` was characterised from its documentation (its Rd
page) and outputs only; its source was not read. What it does, as
`mantel_regression()` now does:

- **Statistics:** OLS with an intercept on the unfolded lower triangles;
  `mrank = TRUE` ranks each vector (average ranks) first. Coefficients
  include the intercept (`Int`), with a p-value.
- **P-values:** each coefficient (the intercept too) is tested two-sided on
  its pseudo-t, `b / se(b)` (|t_perm| >= |t|), as the Rd page says; R^2 and
  F upper-tailed (always the same p-value). The observed samples count as
  one of `nperm`, so `nperm = 1000` is `permutations = 999`.
- **Permutations:** `MRM()` draws them in C from R's RNG, so there is no R
  function to override with `assignInNamespace()` (its namespace has only
  `mrmperm`, a native symbol). Instead, from its outputs: it consumes exactly
  `(n - 1) * nperm` uniforms, and of eight candidate shuffles (with and
  without composing, inverting, and counting the observed one), only one
  reproduced its p-values: nperm - 1 *cumulative* Sattolo shuffles
  (`for i = n..2: j = floor(U * (i - 1)) + 1; swap(i, j)`), each applied to
  the already-permuted response. `data-raw/expected/mantel_regression.r`
  rebuilds them as a permutations matrix (rows = absolute permutations of the
  response's samples), which gives `MRM()`'s p-values exactly.

The user's choices (asked at the start):

- **Factor, character and logical terms** become 0 (same value) / 1
  (different) distances.
- **`strata`** is an argument, as in the rest of the Mantel family (an error
  together with a permutations matrix).
- **`plot()`** draws the standardized coefficients (`std_estimate`, a new
  column: `b * sd(x) / sd(y)` over the pairs, ranks for Spearman), filled
  where `p_value <= alpha`, dashed line at 0, subtitle "R2 = ..., p = ...".
  `alpha` is a `plot()` argument, as in `mantel_correlogram()`.

My own decisions:

- **Signature:** the plan's, plus `strata`, in this order: `dm, ..., data,
  formula, method, permutations, strata, seed, cpus`. `dm` may also be a
  two-sided formula (`dm ~ ph + temp`), as in `permanova()` (`formula_dm()`).
  A formula passed positionally (into `...`) or as `data` is an error that
  says to name it.
- **Term forms:** in `formula`, numeric variables give `|x_i - x_j|`, numeric
  matrices (`cbind(lat, lon)`) Euclidean distances, factors 0/1; terms may be
  expressions (`log(depth)`), but interactions and `Condition()` are errors.
  `model_design()` aligns `data` and checks the variables (missing, NA,
  infinite), and then each term label is evaluated in the aligned data. In
  `...`: dist objects or distance matrices (matched by sample name, as
  `mantel()`'s `y`), coordinates (a data frame, or a matrix that isn't
  square-symmetric with a zero diagonal, as `mantel_correlogram()`), one
  value per sample (vectors, through `per_sample()`, so they follow the same
  rules as `strata`), or a column name of `data`. Formula terms come first,
  then `...`'s. Unnamed `...` terms are named by their deparsed expression
  (or the column name). Duplicate names are an error.
- **Collinear terms** (by `qr()`'s rank, tolerance 1e-7, after centring:
  so a constant term counts) are an error naming them, not `lm()`'s `NA`
  coefficients. A constant `dm` is an error. Fewer than `p + 2` pairs is an
  error (no residual degrees of freedom).
- **Algorithm:** centre the terms and orthonormalize them once (QR), so
  every statistic of a permutation follows from `Q'y_perm` (`m` dot products,
  one pass over the pairs in `C_mantel_regression()` in `src/mantel.c`, the
  explanatory pair `(i, j)` meeting response pair `(p[i], p[j])`): slopes
  `R^-1 Q'y`, `SSR = |Q'y|^2`, `RSS = TSS - SSR`, the intercept from the
  means. The observed RSS is taken from the residuals, exactly. Each
  permutation is one thread's, summed in a fixed order: identical for any
  `cpus`. It reuses `perm_source_*`, `perm_matrix()`, `perm_p()` (with its
  `sqrt(.Machine$double.eps)` tie tolerance), `mantel_dist()`,
  `mantel_align()`, `mantel_is_dist()`, `per_sample()` and `model_design()`.
  It doesn't use `mantel_inputs()`, which handles exactly one `y` and one
  `given`.
- **Result:** the plan's list, with `coefficients` (`term`, `estimate`,
  `std_estimate`, `t`, `p_value`; the intercept first, as `(Intercept)`),
  `r_squared`, `F`, `p_value` (F's), `df` (`c(p, N - p - 1)`), `n_perm`,
  `method`, `n`, `permutation`; attributes `f_perm` and `t_perm` (permutations
  x coefficients). Title "Multiple Regression on Distance Matrices of
  <dm's title>", or the name alone. Its own `print()` (in
  `R/mantel_regression.r`, not `R/print.r`), and `as.data.frame()` gives
  `coefficients`.
- **No `method = 'logistic'`** (ecodist's logistic MRM); documented under
  "Differences".
- **Not in any registry** (a plain function, as the plan says), so no
  aliases yet.
- **Fixture:** the oribatid mites (vegan's `mite`, `mite.env`, `mite.xy`),
  Hellinger distances on geographic distance, substrate density, water
  content, and 0/1 microtopography and shrub matrices; 5 cases, Pearson and
  Spearman, one to three terms, 199 shared permutations.
- **Docs:** `vignettes/stat_tests.Rmd` gains a short example (iris);
  `mantel()`, `mantel_partial()` and `mantel_correlogram()` link to it.


## Outcome (from the implementing session's final report)

**Item 69, done 2026-11-26.** Commit `d496ac7` (function), plus the tick.

- Tests: 33,088 at home and 473 in CRAN mode before merging the other
  sessions' work (32,859 / 469 at 70dff27). `test-mantel_regression.r` has
  225 tests.
- ecodist 2.1.3 `MRM()`, with its permutations rebuilt as a matrix: the 5
  stored mite cases and 300 fresh random cases (n 5–40, 1–4 terms, `mrank`
  FALSE and TRUE, tied and untied, nperm 20–200): every coefficient, R^2 and
  F p-value identical; coefficients within 1.6e-14 of the largest
  coefficient (7.8e-12 relative, on near-zero ones); R^2 within 3e-14
  absolute. ecodist's R^2 has absolute rounding errors (relative 3e-10 when
  R^2 ~ 1e-5); ecodive's R^2 is within 5.4e-14 (relative) of `lm()`'s.
- cpus 1, 2, 3, 4, 8 identical.
- Speed (n = 200, 3 terms, 999 permutations): ecodist 0.10 s, ecodive
  0.06 s on one thread.
- `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed size, the
  container's debug build flags).


## Follow-up: backward elimination (2026-11-26)

The user asked for backward elimination of terms, which Lichstein (2007)
uses with MRM and `ecodist::MRM()` lacks. Decisions:

- **Arguments:** `backward = FALSE` and `alpha = 0.05` on
  `mantel_regression()` (after `permutations`), not a separate function:
  the selected model is an ordinary `mantel_regression()` result. `alpha`
  is the significance level, as in `forward_select()`.
- **Procedure:** while some term's p-value is above `alpha`, remove the term
  with the largest p-value and refit. Permutation p-values tie often, so
  ties go to the smallest |t|, then the last term. The intercept is never
  removed. Every refit uses the same permutations (the same supplied matrix,
  or the same `seed` and `strata`, which `perm_source_*` reproduces), so the
  steps are deterministic and identical for any `cpus`.
- **All terms removed:** the intercept-only model is returned: estimate the
  mean distance, its pseudo-t (identical for every permutation, so p = 1),
  `r_squared = 0`, `F` and `p_value` `NA`, `df = c(0, N - 1)`. `plot()` of it
  is an error ("no terms to plot").
- **Result:** the final model's fields, plus `alpha` and `elimination` (one
  row per removed term: `step`, `term`, `t`, `p_value`, `r_squared_before`,
  `r_squared_after`; zero rows if none was removed). Both are present only
  with `backward = TRUE`. `print()` adds the elimination table.
- **No p-value correction for the selection:** documented as descriptive,
  as after any stepwise selection.
- **`backward = TRUE` with `permutations = 0`** is an error.
- **Code:** the fit moved, unchanged, into `mrm_fit()`, which also handles
  zero terms; `mantel_regression()` checks the inputs, then loops.
- **Tests:** each step is replayed by refitting the remaining terms with
  the same permutations (removed term = largest p, smallest |t| among ties;
  p, t and R² identical); the final model is `identical()` to a plain fit of
  the kept terms; supplied matrix, strata, cpus 1–8; none removed, all
  removed (against `lm(y ~ 1)`'s t); printing and errors.
