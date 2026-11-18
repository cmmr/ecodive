# `principal_response_curves()`: principal response curves

**Priority:** 15. **Tier:** *derived* (a partial RDA). **Requires
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
