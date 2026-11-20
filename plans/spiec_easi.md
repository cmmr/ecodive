# `spiec_easi_mb()`, `spiec_easi_glasso()`: sparse inverse-covariance networks (SPIEC-EASI)

**Priority:** 56. **Tier:** clean room for the SPIEC-EASI pipeline and StARS.
SpiecEasi is GPL (≥ 3; also on Bioconductor 2.0.0) and output-only, and so
are its back ends huge (GPL-2) and pulsar (GPL ≥ 3), checked on CRAN on
2026-11-20. The solvers are a **port** from scikit-learn 1.9.1 (BSD-3, verified
2026-11-20 from its `COPYING`): `_graphical_lasso()` in
`sklearn/covariance/_graph_lasso.py`, and `enet_coordinate_descent_gram()`
in `sklearn/linear_model/_cd_fast.pyx`. Read [README.md](README.md), then
[sparcc.md](sparcc.md) (output structure) and
[assoc_network.md](assoc_network.md).


## What it is

SPIEC-EASI (Kurtz et al. 2015) clr-transforms the counts, then estimates
a **sparse conditional-dependence graph**: two taxa are linked if they are
associated *given all the others*. That removes the indirect correlations a
correlation network shows. The penalty is chosen by StARS stability
selection. It has 2,019 citations, 1,335 of them in 2022–26, and 627
microbiome full-text mentions (OpenAlex/Europe PMC, 2026-11-19). In its own
simulations S-E(MB) "beat SparCC and CCREPE ... under almost all tested
benchmark scenarios". A 2025 preprint finds SparCC more robust on realistic
data. Note both in the docs.


## Names

The paper defines two named variants, "S-E(MB)" and "S-E(glasso)", which are
different estimators (per-node lasso regressions vs a penalised likelihood).
So there is **one function each** (README rule), not SpiecEasi's
`method = "mb" / "glasso"`:

- `spiec_easi_mb()`: Meinshausen–Bühlmann neighbourhood selection;
- `spiec_easi_glasso()`: graphical lasso.

`spiec_easi` is SpiecEasi's `spiec.easi()` in snake_case (rule 5). The base
name is the method's (rule 1), and the variant is a suffix (rule 2), so the
family sorts together. SpiecEasi's other modes (`"slr"`, Ising/Poisson
models) are not planned. No registry alias starts with `spi` (checked
2026-11-20). Both join `feature_assoc()` (assoc_network.md) if it is
accepted. Their arguments avoid `method` and any prefix of `method` /
`counts`.


## Design

### Signatures

```r
spiec_easi_mb(counts, margin = 1L, pseudocount = 1, rule = 'or',
              n_lambdas = 20L, lambda_min_ratio = 1e-3, lambda = NULL,
              n_subsamples = 20L, subsample_size = NULL,
              max_instability = 0.05, max_iters = 100L, tol = 1e-4,
              clust = NULL, seed = 0L, cpus = n_cpus())

spiec_easi_glasso(counts, margin = 1L, pseudocount = 1,
                  n_lambdas = 20L, lambda_min_ratio = 1e-3, lambda = NULL,
                  n_subsamples = 20L, subsample_size = NULL,
                  max_instability = 0.05, max_iters = 100L, tol = 1e-4,
                  clust = NULL, seed = 0L, cpus = n_cpus())
```

- `pseudocount = 1`: the paper adds "a unit pseudo count to the original
  count data". It is part of the method, so it is the default. This is an
  exception to the shared `NULL` rule, as `dirmult_ttest()`'s prior is;
  document it.
- `n_lambdas`, `lambda_min_ratio`, `lambda`: the glossary names (as in
  `coda_lasso()`). The defaults are SpiecEasi's `nlambda = 20`,
  `lambda.min.ratio = 1e-3`. **Unverified:** read `?SpiecEasi::spiec.easi`
  and take the path itself from `se$lambda` in the fixtures.
- `n_subsamples` (StARS `rep.num`, said to default to 20) and
  `subsample_size` (`NULL`: ⌊10√n⌋ when n > 144, else ⌊0.8 n⌋, pulsar's
  documented rule; the paper says "80% in the standard setting").
  **Unverified**; confirm from SpiecEasi/pulsar docs.
- `max_instability`: StARS's β = 0.05 (a cap, hence the `max_` rule).
- `rule` (MB only): `'or'` / `'and'` combine the two neighbourhoods. Every
  reference treats this as an option (huge `sym`, SpiecEasi), so it stays an
  argument. The default follows SpiecEasi (said to be `'or'`; **unverified**).
- `max_iters`, `tol`: the solver's (sklearn's `max_iter = 100`,
  `tol = 1e-4`).

### Pipeline (Kurtz et al. 2015, Methods)

1. \eqn{Z = \mathrm{clr}(x + \text{pseudocount})}, per sample (eq. 1). Then
   each feature is standardised over samples. The paper's eq. 5 uses the
   covariance \eqn{\hat\Gamma}. huge, which SpiecEasi calls, is believed to
   work on the **correlation** matrix. **Unverified**: decide from the
   fixtures (do `se$est$path` estimates match a fit on `cor(Z)` or on
   `cov(Z)`?) and document.
2. **λ path:** \eqn{\lambda_\max = \max_{i\ne j}|S_{ij}|} (sklearn's
   `alpha_max()`, the point above which all edges vanish), then
   `n_lambdas` values log-spaced down to `lambda_min_ratio · λmax`.
   **Check the spacing against `se$lambda`.**
3. **Estimator per λ, warm-started down the path:**
   - **glasso:** port sklearn's `_graphical_lasso(mode = "cd")`: block
     coordinate descent over columns; each column is a lasso on the Gram
     matrix (`enet_coordinate_descent_gram`, `beta = 0`); the dual-gap stop
     (`_dual_gap()`); the 0.95 shrink of the starting off-diagonals. sklearn
     **doesn't penalise the diagonal**. Friedman et al. (2008) do (start with
     \eqn{W = S + \lambda I}), and so does R's glasso by default. Keep an
     internal `penalize_diagonal` switch, set it to whatever matches
     SpiecEasi's fixtures, and say which in the docs. The weights are partial
     correlations \eqn{-\Theta_{ij}/\sqrt{\Theta_{ii}\Theta_{jj}}}.
   - **MB:** for each node, a lasso of \eqn{Z_i} on \eqn{Z_{\neg i}}, using the
     same coordinate-descent kernel on the Gram matrix \eqn{S}. Eq. 3 writes
     \eqn{\tfrac1n\|\cdot\|^2 + \lambda\|\beta\|_1}; huge/glmnet use
     \eqn{\tfrac1{2n}}. The λ scale must match the fixtures' path, so take
     the factor from them. Edges come from `rule`. Weights follow the paper:
     an AND edge gets the average of \eqn{\beta_{ij}} and \eqn{\beta_{ji}};
     an OR-only edge gets its single non-zero β. Clip to [−1, 1] and warn
     if any value falls outside.
4. **StARS** (Liu, Roeder & Wasserman 2010), from the paper:
   - `n_subsamples` subsamples of `subsample_size` samples, without
     replacement. Draw them with PCG, one stream per subsample.
   - Fit the whole path on each, recomputing the clr's standardisation on
     the subsample (check against the fixtures).
   - \eqn{\hat\theta_{ij}(\lambda)} is the edge's frequency;
     \eqn{D(\lambda) = \binom{p}{2}^{-1}\sum_{i<j} 2\hat\theta(1-\hat\theta)}.
   - Monotonise from the largest λ: \eqn{\bar D(\lambda_k) = \max_{l\le k}
     D(\lambda_l)}. Select the **smallest λ with \eqn{\bar D \le}
     `max_instability`**.
   - The result is the full-data fit at that λ. The subsample frequencies at
     that λ are reported as edge stability. The paper uses them to rank edges.
   - If even the largest λ exceeds the cap, take the largest and warn. If
     the selection is the smallest λ on the path, warn that the path may be
     too short (lower `lambda_min_ratio`), as SpiecEasi does.

### Threads, determinism, memory

- One task per subsample (plus one for the full data). Each runs its whole
  path sequentially with warm starts, so it is independent of `cpus`. Sum
  edge counts per λ in subsample order after the join.
- The glasso keeps D × D dense matrices, about 3 D² doubles per thread
  (D = 1,000: 24 MB per thread). Cap the thread count by memory, and give
  the cost in the docs. MB needs only the Gram matrix (shared, read-only)
  plus O(D) per node.
- Every sum runs in fixed order. sklearn's coordinate descent runs in
  cyclic order (`random = False`); keep it.

### Result (feature-correlation structure from sparcc.md)

An `ecodive_dist` over features, a **similarity**: partial correlations
(glasso) or symmetrised β (MB), 0 for non-edges. Its `dist` record holds:

- `metric = "spiec_easi_mb"` / `"spiec_easi_glasso"`;
- `similarity = TRUE`, `features = TRUE`;
- `n_samples`;
- `pseudocount`, `rule`, `lambda_selected`, `n_subsamples`, `subsample_size`,
  `max_instability`, `seed`;
- `sparse = TRUE`, telling `assoc_network(threshold = NULL)` to take
  non-zero entries as edges.

Attributes:

- `stability`: parallel packed edge frequencies at the selected λ, listed
  in the record's `parallel` (see sparcc.md);
- `lambda`: the path;
- `instability`: \eqn{\bar D} per λ.

There are no p-values: SPIEC-EASI gives none. A `plot()` of the λ path
against instability would be a small extra (`plot.ecodive_dist` stays the
heatmap); leave it for later.


## Reference implementations

- **Port:** scikit-learn 1.9.1, as above. It needs the BSD-3 text in the C
  file's header and a new `Files:` stanza in `inst/COPYRIGHTS`. scikit-learn
  is already in `Authors@R`; extend its comment ("graphical lasso,
  coordinate descent"). Fixture generator
  `data-raw/expected/spiec_easi.py`: `sklearn.covariance.graphical_lasso(S,
  alpha)` on fixed S at several α. Match precision and covariance at 1e-10
  with `penalize_diagonal = FALSE`. Use `sklearn.linear_model.Lasso`
  (`fit_intercept = False`, on the standardised Z) for the MB node
  regressions.
- **Output-only:** SpiecEasi (`spiec.easi(X, method = "mb" / "glasso",
  nlambda, lambda.min.ratio, pulsar.params = list(rep.num, seed))`), in
  `data-raw/expected/spiec_easi.r`. Save:
  - the λ path (match at 1e-12);
  - the full-data estimates at each λ (`se$est$path`, `getOptBeta`,
    `getOptiCov`): edges identical, values at 1e-6, since the solvers'
    tolerances differ;
  - StARS's `summary` and selected index.

  StARS subsamples come from R's RNG, so for exact selection checks an
  internal `spiec_easi_details()` takes a subsample index matrix. The
  generator records pulsar's subsample indices if its outputs expose them
  (**unverified**); otherwise it recomputes the StARS rule in base R on
  indices it draws itself, using the huge path on each subsample.


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- **`dist` record** per function ("SPIEC-EASI (neighbourhood selection)",
  "SPIEC-EASI (graphical lasso)"): `similarity = TRUE`, `features = TRUE`,
  the chosen `lambda`, `stars_threshold`, `n_subsamples`, `seed`; the
  input's `transform` record beside it.
- Store the StARS path (lambda vs instability) as an attribute, so a
  `plot(type = 'stability')` can show the selection; the default plot is
  the `plot.ecodive_dist()` heatmap of the selected network's weights.

## Tests

- glasso and lasso kernels against the sklearn fixtures (1e-10).
- The KKT conditions hold at every λ (subgradient check, 1e-8).
- At λ ≥ λmax the graph is empty; on a tiny problem with λ = 0 and n > D,
  glasso gives `solve(S)`.
- StARS rule on hand-made \eqn{D(\lambda)} sequences (monotonisation, cap
  never reached, selection at the end of the path).
- SpiecEasi fixtures as above, on two tables (D = 20, n = 100; D = 60,
  n = 40, so p > n).
- `identical()` for `cpus` 1, 2 and 4; the same `seed` gives the same
  result.
- Errors and edge cases: D < 3; n < 4; a constant feature (zero variance
  after clr: drop it with a warning, or error; decide and test);
  `subsample_size` ≥ n.
- Smoke-test lines.


## Docs

- One Rd page for both functions (`@rdname`), explaining conditional
  dependence vs correlation, MB vs glasso, StARS, and the edge stability.
- Note the cost (n_subsamples × path).
- "Differences from SpiecEasi": PCG subsamples; sklearn's solver; the
  diagonal-penalty choice; the defaults.
- "Same names": none (SpiecEasi's is `spiec.easi()`).
- NEWS. Add `n_subsamples`, `subsample_size`, `max_instability` and `rule`
  to parameter_names.md.


## References

- Kurtz, Z. D., Müller, C. L., Miraldi, E. R., Littman, D. R., Blaser, M. J.,
  & Bonneau, R. A. (2015). Sparse and compositionally robust inference of
  microbial ecological networks. *PLoS Computational Biology*, 11,
  e1004226. doi:10.1371/journal.pcbi.1004226
- Meinshausen, N., & Bühlmann, P. (2006). High-dimensional graphs and
  variable selection with the Lasso. *Annals of Statistics*, 34,
  1436–1462. doi:10.1214/009053606000000281
- Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
  covariance estimation with the graphical lasso. *Biostatistics*, 9,
  432–441. doi:10.1093/biostatistics/kxm045
- Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to
  regularization selection (StARS) for high dimensional graphical models.
  *NeurIPS 23*. arXiv:1006.3316 (no DOI)
- Pedregosa, F., et al. (2011). Scikit-learn: machine learning in Python.
  *JMLR*, 12, 2825–2830 (no DOI).
