# `splsda()`: sparse PLS discriminant analysis

**Priority:** 37. **Tier:** clean room. Write it from Lê Cao, Boitard &
Besse (2011). mixOmics `splsda()` (GPL) is output-only. **Build after
`plsda()`** ([plsda_oplsda.md](plsda_oplsda.md)), whose NIPALS core this
extends. Read [README.md](README.md) first.

Following the one-function-per-method rule, it is its own function, not a
`sparse = TRUE` argument of `plsda()`.


## What it is

Sparse PLS-DA adds an L1 (soft-threshold) penalty to the PLS-DA loading
vectors. Each component then uses only a chosen number of features
(`keepX`), which selects discriminating taxa at the same time as it
ordinates. It's popular with mixOmics users in microbiome work.


## Design

- **Signature:** `splsda(counts, groups, data = NULL, n_features = 10L, margin = 1L, norm = 'none', pseudocount = NULL, scale = TRUE, n_axes = 2L, max_iters = 500L, tol = 1e-6, cpus = n_cpus())`.
  - `n_features` is the exact number of features kept per component: one
    integer, or one per component (mixOmics' `keepX`).
  - Tuning `n_features` by cross-validation is out of scope for v1.
- **Result:** as `plsda()`, with columns `sPLSDA1..` and attribute
  `selected` (each component's chosen features).
- **Wrapper:** joins `COUNT_ORDS`, with aliases `splsda` and
  `sparseplsda`. Check `splsda` against `plsda` — no prefix clash, since
  `s` comes first.
- **Masking:** `mixOmics::splsda`.


## Reference implementations

- **No suitable permissive implementation.** `sparse-pls` 0.1.3 (Python,
  MIT) is L1 *regression*, has no discriminant form, and cites no paper.
- **Write from Lê Cao et al. (2011),** and Lê Cao et al. (2008) for sPLS.
  - Within each NIPALS iteration, soft-threshold the X weight vector so
    that `n_features` features stay non-zero.
  - Deflate as PLS-DA does: mixOmics uses "regression" deflation of X for
    DA. Confirm this from the paper and mixOmics' documented behaviour.
- **Output-only fixtures (GPL):** mixOmics 6.36 `splsda(X, Y, ncomp,
  keepX)`, with `$variates$X`, `$loadings$X`, `selectVar()` and
  `$prop_expl_var`. Also sgPLS (GPL), optionally.


## Tests

- With `n_features` equal to all features, `splsda()` equals `plsda()`.
- Selected features and loadings match mixOmics at 1e-8 (sign-matched),
  using the same scaling.
- Ties at the soft-threshold boundary are broken deterministically.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

An Rd page next to `plsda()`, and NEWS.


## References

- Lê Cao, K.-A., Boitard, S., & Besse, P. (2011). Sparse PLS discriminant
  analysis: biologically relevant feature selection and graphical displays
  for multiclass problems. *BMC Bioinformatics*, 12, 253.
- Lê Cao, K.-A., Rossouw, D., Robert-Granié, C., & Besse, P. (2008). A
  sparse PLS for variable selection when integrating omics data.
  *Statistical Applications in Genetics and Molecular Biology*, 7(1), 35.
