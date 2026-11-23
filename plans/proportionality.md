# `vlr()`, `rho_p()` and `phi_s()`: proportionality between features

**Priority:** 17. **Tier:** port (`vlr()` from
scikit-bio's `vlr` / `pairwise_vlr`, BSD-3) + formula (`rho_p()`, `phi_s()`;
propr is GPL and output-only). Read [README.md](README.md) first, and
[compositional.md](compositional.md) for the conventions the compositional
plans share.


## What it is

Correlating relative abundances is the best-known compositional mistake:
closure alone induces negative correlations, and the correlation between
two taxa changes when a third is dropped (Pearson 1897; Lovell et al.
2015). **Proportionality** asks instead whether two taxa keep a constant
ratio across samples, which doesn't depend on the total and is coherent
under subsetting. All three measures are built on the variance of the
log ratio:

| Function | Measure | Formula | Range | Kind |
|:---------|:--------|:--------|:------|:-----|
| `vlr()` | variance of the log ratio (Aitchison 1986; Lovell et al. 2015) | \eqn{\mathrm{vlr}_{ij} = \mathrm{var}(\log x_i - \log x_j)} | [0, ∞) | dissimilarity; 0 = perfectly proportional |
| `rho_p()` | proportionality correlation ρp (Erb & Notredame 2016) | \eqn{1 - \mathrm{vlr}_{ij} / (\mathrm{var}(z_i) + \mathrm{var}(z_j))} | [−1, 1] | similarity; 1 = perfectly proportional |
| `phi_s()` | symmetric φs (Quinn et al. 2017; propr's `"phs"`) | \eqn{\mathrm{vlr}_{ij} / \mathrm{var}(z_i + z_j)} | [0, ∞) | dissimilarity; 0 = perfectly proportional |

\eqn{z} is the log-ratio-transformed abundance: the clr (default) or an
alr against a chosen feature. \eqn{\mathrm{vlr}} itself doesn't depend on
the transform, because the denominator cancels in \eqn{z_i - z_j}.

`vlr()` is unbounded and scale-dependent; ρp rescales it into a
correlation-like coefficient that propr made popular (its default
metric); φs is the symmetric version of Lovell's φ. They are separate
named measures in the literature, so they are separate functions.

Lovell's original asymmetric φ (\eqn{\mathrm{vlr}_{ij} / \mathrm{var}(z_i)})
is **not** planned: it isn't symmetric, so it can't be a `dist`, and
Quinn et al. (2017) replace it with φs.


## Names

- `vlr()`: the established acronym (scikit-bio's function name).
- `rho_p()`: ρp, not `rho()` (Spearman's ρ).
- `phi_s()`: φs, not `phi()` (the φ coefficient) or `phs()`.
- Docs aliases: "proportionality", "propr".


## Design

### Signatures

```r
vlr  (counts, margin = 1L, pseudocount = NULL, cpus = n_cpus())
rho_p(counts, margin = 1L, denominator = NULL, pseudocount = NULL, cpus = n_cpus())
phi_s(counts, margin = 1L, denominator = NULL, pseudocount = NULL, cpus = n_cpus())
```

- `denominator`: `NULL` (default) uses the clr (the geometric mean of all
  features), as propr's default `ivar = "clr"`. A feature (name or index)
  uses the alr against it, as propr's `ivar = <index>`; that feature is
  then left out of the result. The same argument name as `alr()` and
  `log_ratio()`.
  - propr's `ivar = "iqlr"` (a clr over the features with mid-range
    variance) isn't planned; mention it.
- `pseudocount`: the shared rule (see compositional.md). The recommended
  input is a `replace_zeros_*()` result. Quinn et al. (2017) say propr
  "automatically replaces all zero values with 1" (check that propr 5.x
  still does, from its documentation or outputs only: it is GPL); the docs
  should say how to reproduce it. Replacing only the zeros by 1 isn't
  `pseudocount = 1` (which adds 1 everywhere), so fixtures use zero-free
  tables.
- Pairs are taken **between features** (columns when `margin = 1`); the
  variance is over samples.

### Computation

Everything comes from \eqn{\mathrm{vlr}} and each feature's variance:

- With \eqn{L} the log of the (pseudocounted) table and each feature
  centered over samples, \eqn{\mathrm{vlr}_{ij} = \lVert \tilde L_i - \tilde
  L_j \rVert^2 / (n - 1)}: the **squared Euclidean distance between
  centered log columns**, divided by \eqn{n - 1}. Reuse the
  squared-Euclidean kernel of `src/beta_div.c` over features (the
  transposed problem), which is already parallel and pair-blocked.
- With \eqn{v_i = \mathrm{var}(z_i)} (clr or alr, one pass):
  - \eqn{\rho_{p,ij} = 1 - \mathrm{vlr}_{ij} / (v_i + v_j)};
  - \eqn{\phi_{s,ij} = \mathrm{vlr}_{ij} / (2v_i + 2v_j - \mathrm{vlr}_{ij})},
    since \eqn{\mathrm{var}(a + b) = 2\mathrm{var}(a) + 2\mathrm{var}(b) -
    \mathrm{var}(a - b)}.
- **Don't** form the covariance matrix and use
  \eqn{v_i + v_j - 2c_{ij}} for \eqn{\mathrm{vlr}}: it cancels badly for
  nearly proportional pairs, which are exactly the interesting ones. The
  direct difference keeps 1e-10 agreement.
- The divisor cancels in ρp and φs. For `vlr()` it is \eqn{n - 1}:
  scikit-bio's `vlr()` and `pairwise_vlr()` default to `ddof=1`. Document
  it.
- scikit-bio closes each sample first; that doesn't change any log ratio,
  so skip it.
- Dense log values are needed (every zero becomes a value after the
  pseudocount), so the input is densified per block of features, never
  as a whole n × D copy when that is large. The output is the dominant
  cost: D(D − 1)/2 doubles (20,000 features: 1.6 GB). Warn above 1e8
  pairs, as `normalize()` does, and suggest filtering rare features.

### Result

- A `dist` object **over features** (packed lower triangle, as
  `beta_div()` returns over samples), labelled with the feature names, with
  attributes `method` (`'vlr'`, `'rho_p'`, `'phi_s'`), `denominator`
  (`NULL` or the feature), `pseudocount`, and the usual `Size`, `Labels`,
  `Diag = FALSE`, `Upper = FALSE`.
- `vlr()` and `phi_s()` are dissimilarities, so they work directly as `dm`
  in `pcoa()`, `dist_ord()`, `hclust()` ... to ordinate or cluster
  **features**.
- `rho_p()` is a **similarity** in a `dist` container. This keeps one
  output structure and its memory savings, but `as.matrix()` puts 0 on the
  diagonal (where ρp is 1), and passing it as `dm` would be wrong. The docs
  must say: threshold it (`which(x > 0.7)`) for networks, and use
  `1 - rho_p(x)` (equal to \eqn{\mathrm{vlr}_{ij} / (v_i + v_j)}) if a
  dissimilarity is wanted.
  - **Open question for the user** (ask before implementing): keep the
    `dist` container for `rho_p()`, or return a dense symmetric matrix with
    1 on the diagonal (propr's choice, and safer, but 8 D² bytes)? The
    recommendation is the `dist`, for consistency and memory.


## Wrapper

None for now (see compositional.md): three functions with one signature
could share one, as propr's `metric` does. Revisit if more feature
association measures are added.


## Later (not in this item)

- **Permutation FDR**, as propr's `updateCutoffs()`: permute each feature
  across samples independently, recompute, and report the FDR at a grid of
  cutoffs. A separate function (`rho_p_fdr()`, following the
  `<base>_<modifier>` rule) once `rho_p()` exists; write it from propr's
  documentation and Quinn et al. (2017), as propr is GPL.
- **Differential proportionality** (propd; Erb et al. 2017). Out of scope.


## Reference implementations

- **`vlr()`:** scikit-bio 0.7.4, `vlr(x, y, ddof=1, robust=False)` and
  `pairwise_vlr(mat, ids=None, ddof=1, robust=False, validate=True)` in
  `skbio/stats/composition/_base.py` (BSD-3). Natural log; zeros give NaN
  (no pseudocount). `pairwise_vlr()` computes
  \eqn{c_{ii} + c_{jj} - 2c_{ij}} from `np.cov` of the logs (the
  cancelling route; ecodive computes differences directly, so agreement is
  to rounding, ~1e-12 on well-conditioned fixtures) and its `robust=True`
  raises `NotImplementedError`. `vlr(robust=True)` masks zeros in either
  feature and takes the variance of the rest: **not planned** for
  `vlr()` (it would make each pair use different samples; mention it).
  Fixtures from scikit-bio's tests: `vlr([1,1,2],[3,5,8])` =
  `0.0655828061998637`; `pairwise_vlr([[1,1,2],[3,5,8],[13,21,55]])` has
  condensed sum `0.2857382286903922` (its features are the **columns**).
  Generated fixtures committed as text under
  `data-raw/expected/proportionality/`.
- **`rho_p()`, `phi_s()`:** formulas from Erb & Notredame (2016) and Quinn
  et al. (2017, Eqs. 4 and 5; note \eqn{\phi_s = (1 - \rho_p)/(1 + \rho_p)}).
  propr (GPL-2; archived from CRAN on 2022-11-10, now GitHub tpq/propr
  5.1.8) is **output-only**: install it from GitHub; fixtures from
  `propr::propr(counts, metric = "rho")` and `metric = "phs"`, with
  `ivar = "clr"` and with an alr index, matched at 1e-10 after the same
  zero handling.


## Tests

- Fixtures as above, on a zero-free table and on a sparse table with a
  pseudocount.
- Identities, at 1e-12: `vlr()` equals `var(log(x[, i] / x[, j]))` (with
  the matched `ddof`); `rho_p()` equals the explicit formula; `phi_s()`
  equals `var(z_i - z_j) / var(z_i + z_j)`; `1 - rho_p()` equals
  `vlr / (v_i + v_j)`; `vlr()` is the same with any `denominator`.
- Exactly proportional features give `vlr` 0, `rho_p` 1, `phi_s` 0.
- Scale invariance: multiplying a sample by a constant changes nothing.
- Identical for `cpus` 1, 2 and 4; sparse and dense agree.
- Errors: fewer than 2 samples or 2 features; a `denominator` not in the
  table.
- Smoke-test lines.


## Docs

- One Rd page per function, cross-linked; "Compositional Data" group,
  subtitle "Proportionality".
- A short explanation of why correlation of relative abundances is
  misleading, with Lovell et al. (2015).
- NEWS.


## References

- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
  Chapman & Hall.
- Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., &
  Bähler, J. (2015). Proportionality: a valid alternative to correlation
  for relative data. *PLoS Computational Biology*, 11, e1004075.
- Erb, I., & Notredame, C. (2016). How should we measure proportionality
  on relative gene expression data? *Theory in Biosciences*, 135, 21–36.
- Quinn, T. P., Richardson, M. F., Lovell, D., & Crowley, T. M. (2017).
  propr: an R-package for identifying proportionally abundant features
  using compositional data analysis. *Scientific Reports*, 7, 16252.


## Decisions made while implementing (item 17, 2026-11-19)

Asked of the user at the start:

1. **`rho_p()` returns a `dist`** (as `vlr()` and `phi_s()`), not a dense
   matrix: one output structure, half the memory. The docs say that
   `as.matrix()` puts 0 on the diagonal, that it isn't a `dm`, and to
   threshold it or use `1 - rho_p(x)`.
2. **Zeros:** `pseudocount = NULL` is the shared rule (`log_pseudocount()`:
   none without zeros; else half the smallest non-zero value, added to
   every value, with a warning). An explicit `pseudocount = 0` is accepted
   only for a zero-free table (e.g. `replace_zeros()` output, which goes in
   silently with `pseudocount` attribute 0); with zeros it is an error:
   "`counts` has zeros, and `pseudocount = 0` leaves their logs undefined.
   Use `pseudocount = NULL` or a positive pseudocount, or replace_zeros()
   first." (checked in C, `nnz < n * D`). The Zeros section shows how to
   reproduce propr's handling (zeros alone replaced by 1: `x[x == 0] <- 1`).
3. **Large D: no size check or warning** (the user's choice, against the
   plan's 1e8-pair warning). The docs' Size section gives the cost
   (20,000 features, 1.6 GB) and suggests filtering rare features.
4. **`denominator` is one feature or a set** (`feature_set()`: names,
   indices or a logical vector, any order, stored increasing). One feature
   is the alr, and that feature is left out of the result (its alr is 0);
   a set of two or more uses their mean log (geometric mean) per sample and
   keeps every feature. A set of every feature equals the clr (to
   rounding). The `denominator` attribute is the set's names (indices when
   unnamed), absent for the clr. With one feature, at least 2 others are
   needed.

Made while implementing:

5. **C, sparse, no dense copy** (`src/proportionality.c`,
   `C_proportionality()`). The table is transposed into per-feature lists of
   centred logs of the non-zero values; every zero of feature f has the one
   centred value `log(pc) - mean_f`. A pair walks both lists in sample
   order (`pair_ss()`), then adds `n_both_zero * (z_i - z_j)^2`. So the
   memory is O(nnz + n + D) besides the output, and no n x D dense copy is
   ever made (rather than the plan's "densify per block"). `beta_div()`'s
   squared-Euclidean kernel wasn't reused: it runs over samples with a
   per-sample zero value, and here the zero value is per feature (the
   transposed problem). Explicit zeros are dropped by `new_ecomatrix()`, so
   dense, Matrix and slam input give **identical** results.
6. **Each variance from differences**, never `c_ii + c_jj - 2c_ij`:
   SS_ij = Σ(a_i − a_j)²; vlr = SS / (n − 1); V_f = Σ(a_f − g)² with g the
   centred reference; ρp = 1 − SS / (V_i + V_j); φs = SS / (2V_i + 2V_j −
   SS). Divisor n − 1 (scikit-bio's `ddof = 1`; no `ddof` argument).
   - φs's denominator (var(z_i + z_j)) is ≤ 0 only by rounding, for exactly
     reciprocal pairs (z_i = −z_j, e.g. the two features of a two-feature
     `denominator`): those give `Inf` (`NaN` if SS is 0 too); near-reciprocal
     pairs give a very large, rounding-dominated value. Documented.
   - ρp is `NaN` for two features both constant relative to the reference
     (0/0). Documented.
7. **Threads:** one task per row of the packed triangle (pairs i, j > i),
   rows dealt round-robin; each pair is summed in one fixed order and
   written to its own slot, so results are **identical for any `cpus`**
   (tested for 1, 2, 3, 4, 8 on 150 × 300, sparse and dense, both margins).
   `threads_for(cpus, pairs × n, 1e6)` keeps small jobs serial. Timing:
   500 × 3,000, 10% non-zero, 4.5 M pairs: 3.0 s on 1 thread, 0.8–1.0 s on
   4; dense 500 × 3,000: 1.0 s on 4 threads (`cor()` takes 2.6 s).
8. **Labels:** the kept features' names; for unnamed features, no labels,
   except after an alr drop, when the kept features' positions
   (`"1", "3", ...`) label the result, so it stays unambiguous.
9. **Errors:** fewer than 2 samples or 2 features; fewer than 2 features
   besides a one-feature `denominator`; the `feature_set()` errors (unknown
   names, out-of-range indices, duplicates, NA, empty, names without
   feature names); and the shared count checks.
10. **References:**
    - **scikit-bio 0.7.4** (`data-raw/expected/proportionality.py` / `.r`,
      hex floats): `pairwise_vlr()` and per-pair `vlr()` on 6 tables
      (scikit-bio's own 3 × 3 example; counts 12 × 6; counts with zeros
      30 × 20 plus the default pseudocount; wide 10 × 60; Dirichlet
      proportions 15 × 25; nearly proportional features, 20 × 8, log-ratio
      sd 1e-4). Max abs difference **3.6e-15** (largest values ~10); on
      the nearly proportional table, 1.1e-15 absolute, 3.7e-12 relative,
      from scikit-bio's covariance route (tested at 1e-14 absolute
      everywhere, 1e-12 relative except there). scikit-bio's test values
      (`0.0655828061998637`, condensed sum `0.2857382286903922`) are
      matched at 1e-14.
    - **propr: not run.** Its GitHub repository isn't accessible from this
      session, and installing the CRAN-archived 4.2.6 was not permitted in
      the session either. `rho_p()` and `phi_s()` are tested instead against
      base-R implementations of the published formulas (`1 - var(z_i -
      z_j) / (var(z_i) + var(z_j))`, `var(z_i - z_j) / var(z_i + z_j)`)
      for the clr, an alr at either end, a three-feature set and the full
      set, on five tables (1e-12, or 1e-10 on the nearly proportional one,
      where base R's `var()` rounds), plus the identities `φs = (1 − ρp) /
      (1 + ρp)` and `1 − ρp = vlr / (v_i + v_j)`. The docs say these are
      the coefficients propr computes with `metric = "rho"` / `"phs"`, as
      its paper defines them, without claiming a numeric comparison. A
      later session with propr available should add fixtures from
      `propr(counts, metric = "rho")` / `"phs"`, with `ivar = "clr"` and an
      index, on a zero-free table.
11. **Licensing:** `src/proportionality.c` carries scikit-bio's BSD-3
    text; `R/proportionality.r` points to it; the scikit-bio stanza in
    `inst/COPYRIGHTS` lists both files; the `Authors@R` comment for
    scikit-bio mentions the variance of log ratios.
12. **Docs:** `vlr.Rd` (why not correlation, Zeros, Size; alias
    "proportionality"), `rho_p.Rd` (alias "propr"; similarity caveats;
    `iqlr` and `updateCutoffs()` named as not offered), `phi_s.Rd`
    (Lovell's asymmetric φ not offered; reciprocal pairs). `log_ratio()`
    and `balance()` now link to all three. pkgdown "Compositional Data" →
    "Proportionality". NEWS. `parameter_names.md` updated. Regenerating
    `balance.Rd` also picked up a "Same names" sentence that had changed
    in `R/log_ratio.r` without its Rd being rebuilt.

Still for later (unchanged): the permutation FDR (`rho_p_fdr()`), an
`iqlr` option, and a `metric` wrapper.


## Outcome (from the implementing session's final report)

**Item 17, done 2026-11-19.** Commits `60d304e` (work), `4600cb8` (tick).

- Tests: 10,206 at home, 239 in CRAN mode.
- scikit-bio 0.7.4 `vlr` / `pairwise_vlr`: ≤ 3.6e-15 absolute (3.7e-12
  relative on a nearly proportional table, where scikit-bio loses
  precision).
- `rho_p()` / `phi_s()`: base-R formulas only (≤ 1e-12). **propr not
  compared:** 4.2.6 (CRAN archive) won't load on R 4.3.3 and GitHub was
  blocked; deferred by the user (2026-11-19).
- Speed: 500 × 3,000 (4.5 M pairs) in 3.0 s (1 thread), about 1 s (4).
