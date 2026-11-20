# `mirkat()`: microbiome regression-based kernel association test

**Priority:** 60. **Tier:** clean room. Write it from Zhao et al. (2015),
with p-values from Davies (1980) or Liu, Tang & Zhang (2009). MiRKAT (GPL)
and CompQuadForm (GPL) are output-only. Read [README.md](README.md) first.


## What it is

MiRKAT tests whether overall microbiome composition is associated with an
outcome (continuous or binary), adjusting for covariates.

- It turns each distance matrix into a kernel: K = −½ H D² H, made
  positive semi-definite.
- It computes a variance-component score statistic, Q = rᵀ K r, where r
  are the residuals of the covariate-only model (a linear model, or a
  logistic one for binary outcomes).
- Q is a mixture of χ² variables, so its p-value comes from Davies' method
  or moment matching. Permutation p-values are also possible.
- **Omnibus MiRKAT** combines several distance kernels (e.g. UniFrac
  variants and Bray-Curtis) by taking the minimum p-value, calibrated by
  permutation.

It's in demand in microbiome methods work.

**Scope:** only the original MiRKAT, for continuous and binary outcomes.
The rest of the family (MiRKAT-S for survival, MMiRKAT for multivariate
outcomes, GLMM-MiRKAT, MiRKAT-iQ and others) would be separate functions
named `mirkat_<variant>()`, later.


## Design

- **Signature:** `mirkat(dm, groups, data = NULL, formula = NULL, p_method = 'davies', permutations = 999L, seed = 0L, cpus = n_cpus())`.
  - `dm` is one dist, or a named list of dists for the omnibus test.
  - `groups` is the outcome: a vector, or a column name of `data`.
  - **The outcome type is inferred from `groups`,** not given by an
    argument (no other ecodive test asks for the variable's type, and
    `dist_test()` already checks it from the values):
    - numeric → continuous (a linear covariate model);
    - categorical with two categories → binary (a logistic model);
    - categorical with more than two categories → an error, saying that
      MiRKAT handles continuous and binary outcomes only.

    Document that a numeric 0/1 vector is treated as continuous, so a
    binary outcome should be a factor, character or logical vector. Record
    the inferred type in the result (`outcome`).
  - `formula` holds one-sided covariates (`~ age + sex`) over `data`.
  - `p_method` is `'davies'` (exact, as MiRKAT's default), `'liu'`
    (moment matching) or `'permutation'`.
- **Result:** an `ecodive_test` of class
  `c('ecodive_mirkat', 'ecodive_test')`: `statistic = c(Q = ...)`,
  `p_value`, `z` (`NA` unless permuted), `n_perm`, and `outcome`
  (`'continuous'` or `'binary'`, as inferred).
  - For several kernels: a `kernels` data frame of each kernel's Q and p,
    plus the omnibus p as the core `p_value`.
- **Wrapper:** joins `dist_test()` (it tests `dm` against one variable),
  with `groups` taking `numeric` or `categorical` (binary).
  - Add a `DIST_TESTS` row:
    `Microbiome Regression-based Kernel Association Test, mirkat, categorical numeric, kernelassociation`.
  - `dist_test()`'s group-type check must allow only two categories for
    `mirkat`. Add the same check in `mirkat()` itself, as part of
    inferring the outcome type.
- **Masking:** none (MiRKAT's is `MiRKAT`).


## Reference implementations: licensing traps

- **No permissive implementation of MiRKAT.**
- **Don't port `chi2comb` (limix, MIT) or `chiscore` (MIT).** Their READMEs
  say chi2comb repackages CompQuadForm's `davies`, which is GPL ≥ 2. The
  MIT label doesn't remove that.
- Davies' original `qfc.c` (AS 155) has no stated license, so it's
  unverified and should be treated as unavailable.
- **Write the p-value methods from the papers:**
  - Davies (1980): numerical inversion of the characteristic function.
    Write it from the paper's algorithm, not from any code.
  - Liu, Tang & Zhang (2009): four-moment matching to a scaled noncentral
    χ². Simple and robust, so a good first implementation.
  - Imhof (1961) as an alternative exact method.
- **Output-only fixtures:** MiRKAT 1.2.3 (GPL) `MiRKAT(y, X, Ks, out_type,
  method = 'davies' / 'liu' / 'permutation')`; CompQuadForm `davies()` and
  `liu()` for the mixture p-values themselves.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- The title is automatic: `new_ecodive_test(from = dv$title)` gives
  "Microbiome Regression-based Kernel Association Test of Bray-Curtis
  Dissimilarity" for one kernel.
- Store what the plot needs on the result (see Result).
- Omnibus over several kernels: no single input title, so the title is
  the method's name alone. The `kernels` table gets a `metric` column,
  each kernel's name in the `dm` list (else its `method`), for
  per-kernel labels.

## Tests

- Q equals the formula rᵀKr, computed in R in the test.
- The eigenvalues of the mixture are correct.
- Davies and Liu p-values match CompQuadForm's on fixed eigenvalue sets
  (stored as a fixture) at 1e-8.
- End-to-end p-values match MiRKAT for continuous and binary outcomes,
  with and without covariates, for one and several kernels. The outcome
  type comes from `groups` alone: a numeric `groups` matches MiRKAT's
  `out_type = 'C'`, and a two-category one its `out_type = 'D'`.
- Outcome inference: numeric → continuous; a two-level factor, character
  or logical → binary; three or more categories → a clear error; a
  numeric 0/1 vector → continuous (as documented).
- The permutation p-value agrees with Davies' statistically.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

- An Rd page; a "Kernel tests" note in `vignettes/stat_tests.Rmd`; an
  entry in the reference index under Hypothesis Tests.
- NEWS.


## References

- Zhao, N., Chen, J., Carroll, I. M., Ringel-Kulka, T., Epstein, M. P.,
  Zhou, H., Zhou, J. J., Ringel, Y., Li, H., & Wu, M. C. (2015). Testing
  in microbiome-profiling studies with MiRKAT, the microbiome
  regression-based kernel association test. *American Journal of Human
  Genetics*, 96, 797–807.
- Davies, R. B. (1980). Algorithm AS 155: the distribution of a linear
  combination of χ² random variables. *Applied Statistics*, 29, 323–333.
- Liu, H., Tang, Y., & Zhang, H. H. (2009). A new chi-square approximation
  to the distribution of non-negative definite quadratic forms in
  non-central normal variables. *Computational Statistics & Data
  Analysis*, 53, 853–856.
- Imhof, J. P. (1961). Computing the distribution of quadratic forms in
  normal variables. *Biometrika*, 48, 419–426.
