# `mantel_regression()`: multiple regression on distance matrices (MRM)

**Priority:** 59. **Tier:** clean room. Write it from Legendre, Lapointe &
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
