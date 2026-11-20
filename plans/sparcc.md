# `sparcc()`: sparse correlations for compositional data (SparCC)

**Priority:** 39. **Tier:** clean room. Write it from Friedman & Alm
(2012, *PLoS Comput Biol*, CC-BY). The original Python 2 code
(bitbucket.org/yonatanf/sparcc) is gone and its license is unknown. FastSpar
and SpiecEasi's `sparcc()` are GPL-3. The MIT/BSD Python ports (dlegor/SparCC,
fast_sparCC) say they follow the original closely, so where their code comes
from is unclear. **Don't read any of their source.** FastSpar and SpiecEasi
are output-only: their README/Rd documentation and their outputs may be used.
Read [README.md](README.md) first, then [compositional.md](compositional.md).
[proportionality.md](proportionality.md) covers the `vlr()` machinery this
plan reuses.


## What it is

SparCC estimates the **Pearson correlations between the log *absolute*
abundances** of taxa (the "basis") from relative abundances alone. That is
what users mean when they report a "co-occurrence network". The basis
correlations can't be identified from compositions in general. SparCC makes
them approximately identifiable by assuming that **most pairs are uncorrelated**
(a sparse network) and that there are many taxa. Uptake: 3,045 OpenAlex
citations, 1,848 of them in 2022–26, and 1,363 microbiome full-text mentions;
FastSpar adds 434 citations (retrieved 2026-11-19; `research_notes/.../networks_typing.md`).

### Why it is not `rho_p()` (corrects compositional.md of 2026-11-19)

compositional.md said `rho_p()` "answers the same question without SparCC's
sparsity assumption". That was wrong. Both start from the variation matrix
\eqn{t_{ij} = \mathrm{var}(\log x_i/x_j)} (`vlr()`), but they normalise it
differently and estimate different things:

- \eqn{\rho_p = 1 - t_{ij}/(v_i + v_j) = 2\,\mathrm{cov}(z_i,z_j)/(v_i+v_j)},
  with \eqn{v_i} the variance of the clr (or alr) coordinate. This is Lin's
  concordance coefficient of the clr values. It is 1 only when the log
  ratio is **constant** (slope 1 on the log scale), and it depends on the
  clr reference. It needs no assumption, but it is not a correlation of
  abundances. Pairs related by a power law with an exponent other than 1, or
  with unequal variances, score low.
- SparCC solves for the **basis variances** \eqn{\omega_i^2} instead, and
  returns \eqn{\rho_{ij} = (\omega_i^2+\omega_j^2-t_{ij})/(2\omega_i\omega_j)}:
  a Pearson-form correlation of the log basis abundances. This is valid only
  under the sparsity assumption.
- When \eqn{\omega^2 \approx v}, SparCC ≈ the Pearson correlation \eqn{r} of the
  clr values, while \eqn{\rho_p = r \cdot 2\sigma_i\sigma_j/(\sigma_i^2+\sigma_j^2)
  \le |r|}.

So `rho_p()` answers "do these two keep a constant ratio?" and `sparcc()`
answers "are their absolute abundances correlated?". The docs of both
functions should say this in one paragraph and link to each other.


## Names

- `sparcc()`: the method's own name (rule 1), as SpiecEasi and NetCoMi
  (`measure = "sparcc"`) call it. No registry has an alias starting `spar`
  (checked against METRICS, ORDINATIONS, COUNT_ORDS, DIST_TESTS and
  ZERO_REPLACEMENTS on 2026-11-20; re-check DA_TESTS once it exists).
- **Same names in other packages:** `SpiecEasi::sparcc()`, with the same
  meaning.
- It joins the association wrapper proposed in
  [assoc_network.md](assoc_network.md) (`feature_assoc()`), if accepted. Its
  arguments must avoid `method` and any prefix of `method` or `counts`; the
  ones below do.


## Design

### Signature

```r
sparcc(counts, margin = 1L, n_draws = 100L, max_exclusions = 20L,
       exclusion_threshold = 0.1, permutations = 0L, clust = NULL,
       seed = 0L, cpus = n_cpus())
```

- `n_draws`: Dirichlet posterior draws; the result is the per-pair
  **median** over draws. The paper used 100. FastSpar's `--iterations` and
  SpiecEasi's `iter` are this count (their defaults are said to be 50 and 20,
  **unverified**: check `fastspar --help` and `?SpiecEasi::sparcc`).
- `max_exclusions`: a cap on the exclusion iterations (the `max_<unit>`
  rule). The paper used 20. FastSpar's `--exclude_iterations` and SpiecEasi's
  `inner_iter` are said to default to 10 (**unverified**).
- `exclusion_threshold`: a pair with \eqn{|\rho| >} this is "strongly
  correlated" (the paper's 0.1, also FastSpar's `--threshold`). Don't name it
  `threshold`: that is `assoc_network()`'s edge cutoff (compositional.md
  glossary).
- `permutations`: the number of null datasets for p-values (0 = none); see
  below. A count only: a supplied permutation matrix makes no sense here,
  because every feature is permuted independently.
- `clust`, `cpus`: as `vlr()`. `seed`: PCG, as `dirmult_ttest()`.
- **No `pseudocount`.** SparCC handles zeros with its Dirichlet prior (uniform,
  so every count gets +1). Raw **integer counts** are required: warn when
  values aren't whole numbers, as `ancombc()` does, because the posterior's
  spread depends on the counts' size.
- **Defaults: ask the user** whether to use the paper's (100 / 20 / 0.1,
  proposed) or FastSpar's. The docs give the settings that reproduce FastSpar
  and SpiecEasi.

### Algorithm (Friedman & Alm 2012, Methods; equations read from the paper's figures)

For each draw \eqn{d = 1..n_\mathrm{draws}}:

1. **Fractions.** The posterior given counts \eqn{N_s} is \eqn{\mathrm{Dir}(N_s+1)}.
   Only log ratios are used, and these don't depend on the normalisation, so
   draw \eqn{g_{sf} = \log G(N_{sf}+1)} with `rng_log_gamma()` (`src/rng_dist.c`).
   Each (draw, sample) pair gets its own stream,
   `rng_seed(seed, d * (n + 1) + s)`, as `src/dirmult.c` does, and features are
   drawn in order so that sparse and dense input give the same numbers.
2. **Variation matrix.** \eqn{t_{ij} = \mathrm{var}_s(g_{si}-g_{sj})}. Use the
   difference kernel of `src/proportionality.c` (`pair_ss()`), here on a dense
   n × D block, because every value is a draw. Never use
   \eqn{c_{ii}+c_{jj}-2c_{ij}} (see proportionality.md). The divisor cancels
   in \eqn{\rho}. Use n − 1 as `vlr()` does; it matters only for `variances`.
3. **Basis variances.** \eqn{t_i = \sum_{j\ne i} t_{ij} \approx d\,\omega_i^2 +
   \sum_{j\ne i}\omega_j^2} with \eqn{d = D-1} (eq. 13). The matrix is
   \eqn{M = (D-2)I + \mathbf{1}\mathbf{1}^\top}, so
   \eqn{\omega_i^2 = (t_i - \sum_k t_k/(2D-2))/(D-2)}. This needs D ≥ 4 (the
   paper's minimum).
4. **Correlations.** \eqn{\rho_{ij} = (\omega_i^2+\omega_j^2-t_{ij})/(2\omega_i\omega_j)}
   (eq. 10).
5. **Exclusion loop**, at most `max_exclusions` times:
   - Take the pair with the largest \eqn{|\rho|} that isn't already excluded.
     Stop if its \eqn{|\rho| \le} `exclusion_threshold`; otherwise add it to
     the excluded set \eqn{E}.
   - A component whose pairs are all excluded is dropped from this draw
     (its correlations are `NA`). Stop if only 3 components remain.
     "Re-estimate fractions relative to the subset" changes no log ratio, so
     \eqn{t} is just restricted.
   - Recompute \eqn{t_i^{(n+1)} = \sum_{j \notin c_i, j\ne i} t_{ij}} (eq. 15).
   - Solve the reduced system. The paper prints only eq. 15; its consequence
     for the matrix is **derived here, not printed**: \eqn{M = (D-2)I +
     \mathbf{1}\mathbf{1}^\top - \sum_{(i,j)\in E} u u^\top} with
     \eqn{u = e_i + e_j}, i.e. \eqn{M_{ii} = D-1-|c_i|} and \eqn{M_{ij} = 0} for
     excluded pairs. Each exclusion is a rank-1 update. Solve with Woodbury
     (\eqn{A^{-1}} is closed form, at most 21 columns): O(D·K + K³), never a
     D × D factorisation. Check against a dense `solve()` in the tests.
   - Recompute every \eqn{\rho} with the new \eqn{\omega^2}.
   - Ties for the largest \eqn{|\rho|}: take the pair (i < j) first in
     packed order. **Check against FastSpar's outputs** on data with ties.
6. **Across draws:** the per-pair median of the non-`NA` values (R's median:
   the mean of the two middle values when the count is even). Also take the
   per-feature median \eqn{\omega^2} (`variances`).

Edge cases to settle in the tests and document:
- \eqn{\omega_i^2 \le 0} from the approximation gives `NaN` for that feature's
  pairs in that draw. Leave them out of the median, and warn once with the
  count.
- \eqn{|\rho| > 1} can happen when too many pairs are excluded (the paper
  says so). Keep the value and warn.

### p-values (FastSpar's approach, from its README and paper only)

The paper's "bootstrap" builds null datasets in which "all components are
uncorrelated", runs SparCC on each with the same settings, and reports the
fraction at least as extreme. FastSpar (Watts et al. 2019) replaces the
biased estimator with Phipson & Smyth (2010).

- **Null dataset:** permute each feature's counts across samples
  independently (Fisher–Yates on PCG, stream per (permutation, feature)).
  **Unverified:** that this is what FastSpar's `fastspar_bootstrap` does.
  Check its outputs: is each column of a bootstrap table a permutation of the
  original column (same multiset), or a resample with replacement? Follow
  whatever it does, and record the finding in this plan.
- Each null dataset runs the full SparCC: `n_draws` draws, its own streams,
  offset past the observed run's.
- `p = (1 + #{|r_b| ≥ |r_obs|})/(B + 1)`, two-sided, via `perm_p(alternative =
  'two.sided')`, which already has Phipson & Smyth's +1. No adjustment here:
  `assoc_network()` adjusts (see [assoc_network.md](assoc_network.md)).
- Cost: (B + 1) × n_draws full runs. Document it, and suggest
  `n_draws = 20L` for the null runs. (Ask the user whether to add a separate
  `n_null_draws`; FastSpar uses fewer draws for its bootstraps.)

### Threads and memory

- Parallelise over (draw) tasks. Each task holds an n × D block of log
  draws, the packed D(D−1)/2 \eqn{t}, and its \eqn{\rho}.
- Medians need every draw's value per pair. Store the draws as a
  `n_draws × pairs` float64 array when it fits; otherwise do a second pass
  per block of pairs that regenerates the draws (deterministic streams make
  this exact).
- Sum in fixed order; one task per draw. Results are `identical()` for any
  `cpus`.
- No sparse input path: every value is drawn. Densify per sample.

### Result (the shared "feature correlations" structure)

As `rho_p()`: a `dist` **over features** of class `c("ecodive_dist", "dist")`,
a **similarity** in a dist container. Its `dist` record (`ecodive_step`)
holds:

- `"SparCC Correlation"`, `metric = "sparcc"`;
- **`similarity = TRUE`** (new; see below);
- `n_samples`;
- `n_draws`, `max_exclusions`, `exclusion_threshold`, `permutations`, `seed`.

Attributes:

- `p_value`: a numeric vector **parallel to the packed triangle** (same
  length and order), present only when `permutations > 0`;
- `variances`: the named per-feature median \eqn{\omega^2}.

This is the convention that `spiec_easi_*()` and `secom_*()` follow and that
`assoc_network()` reads. It needs two small changes to existing code:

- **`similarity` flag.** `dist_clust_order()` (R/dist_clust.r) and
  `plot.ecodive_dist()` (R/plot_dist.r, the `features <- m %in% c('vlr',
  'rho_p', 'phi_s')` line) test for the metric `rho_p` by name. Make them
  read the record's `similarity` flag (and a `features = TRUE` flag), and set
  both flags in `proportionality()` for `rho_p`. `vlr` and `phi_s` get
  `features = TRUE` only.
- **`dist_reorder()`** keeps attributes but doesn't reorder parallel packed
  vectors. Make it reorder every attribute named in the record's
  `parallel = c("p_value")`, so `clust = "average"` keeps p-values aligned.


## Reference implementations (output-only)

- **FastSpar** (C++, GPL-3; Watts et al. 2019): `fastspar --iterations
  --exclude_iterations --threshold`, `fastspar_bootstrap`, `fastspar_pvalues`.
  Fixtures come from its median correlation (and covariance) outputs.
- **SpiecEasi** `sparcc(data, iter, inner_iter, th)` (GPL-3).
- Results are stochastic, so parity is distributional. FastSpar, SpiecEasi
  and ecodive draw different random numbers. Compare medians on data where
  the posterior is tight (counts ≥ 1,000 per cell; the draw noise is then
  ~1e-3), at 0.02 absolute. **Deterministic parity:** `sparcc_details()`, an
  internal function, takes a user-supplied list of log-fraction matrices
  (one per draw) in place of the Dirichlet draws. The fixture generator
  computes the same thing with a base-R implementation of steps 2–6. Match
  it at 1e-12 (exclusion order and \eqn{\rho}).


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- The `dist` record described above is the record (`similarity = TRUE`,
  `features = TRUE`, `n_samples`, `n_draws`, `max_exclusions`,
  `exclusion_threshold`, `permutations`, `seed`), plus the input's
  `transform` record copied beside it.
- `plot()`: the existing `plot.ecodive_dist()` heatmap, which must read
  `similarity` (diverging colour scale centred on 0, diagonal 1) and
  `features` (axis labels say "features"); title from the record's name.
- `p_value` is reordered with the triangle by `dist_reorder()` (`clust`).
- Downstream: `assoc_network()` keeps this record as its input record.

## Tests

- `sparcc_details()` against the base-R reference at 1e-12, including:
  exclusions that drop a component; the 3-component stop; ties.
- The closed-form \eqn{\omega^2} equals `solve(M, t)`, and Woodbury equals a
  dense solve after 1..20 exclusions (1e-10).
- Distributional comparison with FastSpar fixtures (0.02, as above), and
  the null check: on per-feature permuted data, edges with |r| > 0.3 are
  rare (the paper's "no significant correlations in shuffled data").
- Invariance: scaling every count by 10 tightens the result but leaves the
  `sparcc_details()` path unchanged; feature order doesn't matter (a
  permuted input gives a permuted output).
- `identical()` for `cpus` 1, 2 and 4; same `seed` gives the same result;
  sparse equals dense.
- p-values: in (0, 1], multiples of 1/(B+1), and aligned after `clust`.
- Errors: D < 4; n < 3; negative or non-finite counts; a warning for
  non-integer counts.
- A smoke-test line (D = 10, n_draws = 5).


## Docs

- Rd page in "Compositional Data" → *Associations and networks*.
- Explain what is estimated, the sparsity and many-taxa assumptions (the
  paper's rule of thumb: an effective number of components ≥ 50 for standard
  correlation; SparCC is more robust), the "Why not rho_p()" paragraph, the
  similarity-in-a-dist caveats (copy `rho_p()`'s), and the cost of p-values.
- "Differences from FastSpar/SpiecEasi": RNG, the defaults, how the median
  is taken, and that negative variances are left out.
- Link to `assoc_network()`. NEWS. Add `n_draws` to parameter_names.md (it
  is already in compositional.md's table) along with `max_exclusions` and
  `exclusion_threshold`.


## References

- Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from
  genomic survey data. *PLoS Computational Biology*, 8, e1002687.
  doi:10.1371/journal.pcbi.1002687 (CC-BY; 3,045 citations, 1,848 in
  2022–26, OpenAlex 2026-11-19)
- Watts, S. C., Ritchie, S. C., Inouye, M., & Holt, K. E. (2019). FastSpar:
  rapid and scalable correlation estimation for compositional data.
  *Bioinformatics*, 35, 1064–1066. doi:10.1093/bioinformatics/bty734
- Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
  zero. *Statistical Applications in Genetics and Molecular Biology*, 9, 39.
  doi:10.2202/1544-6115.1585
- Erb, I., & Notredame, C. (2016), for \eqn{\rho_p}; see proportionality.md.
