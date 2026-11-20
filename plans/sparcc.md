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
a **similarity** in a dist container. Attributes:

- `method = "sparcc"`, and `display$title` "SparCC Correlation" (from
  `dist_result()`);
- `n_samples`;
- `p_value`: a numeric vector **parallel to the packed triangle** (same
  length and order), present only when `permutations > 0`;
- `variances`: the named per-feature median \eqn{\omega^2}.

That it is a similarity over features is a fact about the method, kept in
the `FEATURE_ASSOCS` registry (planned in
[assoc_network.md](assoc_network.md); add it with this item if it hasn't
landed), not a flag on the result.

This is the convention that `spiec_easi_*()` and `secom_*()` follow and that
`assoc_network()` reads. It needs two small changes to existing code:

- **Registry lookups.** `dist_clust_order()` (R/dist_clust.r) and
  `plot.ecodive_dist()` (R/plot_dist.r, the `features <- m %in% c('vlr',
  'rho_p', 'phi_s')` line) test `method` against `vlr`, `rho_p` and
  `phi_s` by name. Replace those tests with lookups of `method` in the
  registry: over features when it has a row, and a similarity when its
  `kind` is `'similarity'` (`rho_p`, `sparcc`, ...; not `vlr` or `phi_s`).
- **`dist_reorder()`** keeps attributes but doesn't reorder parallel packed
  vectors. Make it reorder every attribute named in an internal constant
  (`DIST_PARALLEL = c("p_value", "stability", "n_pairs")`) with the
  triangle, so `clust = "average"` keeps p-values aligned.


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


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Title: `dist_result()` sets `display$title` "SparCC Correlation"
  automatically.
- `plot()`: the existing `plot.ecodive_dist()` heatmap, which looks up
  `method` in the registry: a similarity gets a diverging colour scale
  centred on 0 (diagonal 1), and over features the axis labels say
  "features".
- `p_value` is reordered with the triangle by `dist_reorder()` (`clust`),
  since it is in `DIST_PARALLEL`.
- Downstream: `assoc_network()` reads `method`, `n_samples` and
  `p_value`, and titles its network "Association Network of SparCC
  Correlation".

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


## Decisions made while implementing (item 39, 2026-11-20)

Asked of the user at the start:

1. **Defaults: the paper's** (`n_draws = 100L`, `max_exclusions = 20L`,
   `exclusion_threshold = 0.1`). Checked from the references' help:
   FastSpar's are 50 / 10 / 0.1 (`fastspar --help`) and SpiecEasi's
   `sparcc(iter = 20, inner_iter = 10, th = 0.1)` (its Rd). The docs give
   both.
2. **No `n_null_draws`.** The null runs use `n_draws`, so null and observed
   medians carry the same Monte Carlo noise. The P-values section suggests
   a smaller `n_draws` for the whole call when p-values are too slow.
3. **Registry only, no wrapper.** `FEATURE_ASSOCS` (internal,
   `R/feature_assoc.r`) has rows for `vlr`, `rho_p`, `phi_s` and `sparcc`,
   with columns `name` (the title), `id`, `alt_ids` (all empty), `kind`,
   `lower`, `upper`, `fisher_z`, `sparse` and `strict`. The range is two
   numeric columns rather than assoc_network.md's `range` string: `'[-1,1]'`
   contains the table's separator, and plots need numbers. `feature_assoc()`,
   `list_feature_assocs()` and `match_feature_assoc()` are left to item 40.

Made while implementing:

4. **References: run, output only.** FastSpar 1.0.0 and SpiecEasi 1.1.1
   came from the jsDelivr mirror of their GitHub repositories. The first
   attempts to build FastSpar and `R CMD INSTALL` SpiecEasi were refused by
   the session's permission classifier; the user then allowed them.
   FastSpar was built with `autoreconf -fi && ./configure && make` (the
   mirror loses `autogen.sh`'s executable bit); SpiecEasi's CRAN
   dependencies `huge` and `pulsar` went into the scratch library.
   **Exposure record:** no source file of either package was opened. Read:
   SpiecEasi's `DESCRIPTION` (its `Imports`), its Rd for `sparcc` and
   `sparccboot` (via `tools::Rd2txt`), FastSpar's `README.md` (grep for
   usage and p-values) and `--help` output, the programs' messages, and the
   top-level directory listings. Fixtures: `data-raw/expected/sparcc_refs.r`
   (input tables; SpiecEasi's outputs), `data-raw/expected/sparcc_fastspar.sh`
   (FastSpar's), folded into `sparcc.rds` by `data-raw/expected/sparcc.r`
   (with the base-R cases). Results, same settings, largest absolute
   difference over all pairs:

   | Table | Reference | Settings | ecodive draws | Max diff | Mean diff |
   |:------|:----------|:---------|--------------:|---------:|----------:|
   | 40 × 20, counts ~ e^9 | FastSpar | 100 draws, 10 excl. | 100 | 0.0026 | 0.0007 |
   | same | FastSpar | 17 (its 20; see 21) | 100 | 0.0014 | 0.0004 |
   | same | SpiecEasi | 100, 10 | 100 | 0.0017 | 0.0005 |
   | 40 × 8, counts ~ e^8 | FastSpar | 200, 5 (its 10) | 200 | 0.0012 | 0.0005 |
   | same | SpiecEasi | 200, 10 | 200 | 0.0013 | 0.0005 |
   | 60 × 30, counts ~ 20 | FastSpar | 1000, 20 | 1000 | 0.0188 | 0.0025 |

   Variances against FastSpar's covariance diagonal: within 0.4% (tight)
   and 1.8% (moderate). Tests: 0.02 absolute (tight), mean < 0.005 and max
   < 0.05 (moderate). SpiecEasi and FastSpar differ from each other by
   0.0021 on the tight table. The moderate table's 0.019 is Monte Carlo
   noise: single draws spread identically in both (median per-draw sd of
   ρ 0.0582 FastSpar, 0.0578 ecodive; of the variances 0.212 and 0.209),
   and two ecodive seeds differ by 0.018 (mean 0.0025, the 0.0026 expected
   for independent 1,000-draw medians). FastSpar's own seeds differ less
   (1 vs 2: 8e-4; 1 vs 100: 0.0038), as they would if FastSpar's seeds
   shared most of their iterations' streams (not checked further).
   Earlier checks kept as well:
   - the deterministic steps against a base-R implementation
     (`data-raw/expected/sparcc.r`, dense `solve()` of the variance
     system, `var()`, `median()`): 6 cases, per-draw ρ, variances,
     exclusion order, stop reason and medians. Largest difference
     **7.2e-15** (tested at 1e-12; exclusion order identical);
   - a tight posterior (counts around e^9): the median over 40 Dirichlet
     draws is within 0.02 of SparCC of the observed log fractions;
   - recovery of planted correlations (ρ > 0.8 and < −0.8; 90% of the
     other pairs within ±0.3), and the null check (feature-shuffled data,
     n = 200: fewer than 0.5% of pairs with |r| > 0.3).

   Ties can't be compared with the references, which always draw. Other
   notes: FastSpar's long option is `--exclude_iterations` (its `--help`
   says `--exclusion_iterations`, which it rejects), and SpiecEasi's
   `inner_iter = 0` doesn't give the basic SparCC (0.17 from both FastSpar
   and ecodive with no exclusions, which agree to 0.0017).
5. **The variance step is least squares.** With B the (unsigned)
   incidence matrix of the pairs kept, M = B'B and t_i = B't, so the
   variances are the least-squares fit of t_ij ≈ ω_i² + ω_j², and
   ρ_ij is proportional to that fit's residual. So a pair whose exclusion
   would drop a component (its last pair) or make M singular (the kept
   pairs become bipartite) has leverage 1 and a residual, so ρ, of
   exactly 0. **With any threshold above rounding error, neither can
   happen.** They happen only at `exclusion_threshold = 0`, where rounding
   (|ρ| around 1e-16) decides. There they are handled, and tested against
   dense solves using the C code's own exclusion list (30 tables, every
   stop reason reached):
   - a component left with no pair is dropped from the draw (its ρ and
     variance are `NA`);
   - fewer than 4 components left: stop, keeping the last estimates
     (the plan's "stop at 3");
   - M singular (a Cholesky pivot of Woodbury's K × K matrix ≤ 1e-10):
     undo that exclusion and stop.

   The fixtures use a threshold of 1e-6 so that rounding can't decide;
   their first version used 0, and the base-R and C versions then
   disagreed on which zero-correlation pair came next. The docs'
   "Other implementations" section says this.
6. **Solving.** No exclusions: ω_i² = (t_i − Σt_k/(2m − 2))/(m − 2) over
   the m active components. With K exclusions: Woodbury on M = A − UU',
   A = (m − 2)I + 11' (A⁻¹ = (I − 11'/(2m − 2))/(m − 2)), so
   S = I − U'A⁻¹U (K × K, entries from how many components two excluded
   pairs share), Cholesky, ω² = y + A⁻¹U S⁻¹ U'y with y = A⁻¹t. O(D² + K³)
   per iteration, never an m × m factorisation. Matched against dense
   `solve()` after 0 to 20 exclusions (the fixtures, and the threshold-0
   tests at 1e-10 relative).
7. **The exclusion step.** The largest |ρ| strictly above the threshold,
   among pairs not yet excluded and with both components active; `NaN`
   (a non-positive variance) is never chosen; ties go to the first pair
   in packed order (tested with three identical features). Excluded pairs
   keep their ρ, recomputed from each new set of variances (eq. 10).
8. **Draws.** `rng_log_gamma(N + 1)` per value (shape ≥ 1, so always
   Marsaglia–Tsang). Stream `(b · n_draws + d)(n + 1) + s` for sample s
   of draw d of dataset b (0 = observed, as the plan's `d(n + 1) + s`;
   1..B the null datasets). Draw d doesn't depend on `n_draws` (tested:
   the first 20 of 30 draws equal a 20-draw run's). Each draw's logs are
   centred per feature and t_ij = Σ(a_si − a_sj)²/(n − 1), summed in
   sample order, on a dense feature-major D × n block per thread.
9. **Null datasets** (`permutations`): **FastSpar's resampling.** From
   `fastspar_bootstrap`'s outputs (99 tables of the moderate table): no
   feature's row is a permutation of its original (0 of 2,970), every value
   comes from the same feature's own counts (2,970 of 2,970), rows have
   fewer distinct values (resampling with replacement), and each feature
   uses its own sample indices. So each feature's counts are resampled
   across the samples with replacement, independently. ecodive does the
   same: stream `2^62 + (b − 1)D + f`, sample j = ⌊u n⌋ with `rng_unif()`
   for each position; then the full SparCC with the dataset's own draw
   streams. (The first version shuffled, Fisher–Yates, as the plan
   proposed; changed once FastSpar's outputs were seen.) SpiecEasi's
   documented `sparccboot()` null permutes each feature instead
   (`apply(data[indices, ], 2, sample)`); the docs say so.
10. **p-values** are computed in C exactly as `perm_p(alternative =
    'two.sided')` would from the null medians (hit when |r_b| ≥ |r| −
    √ε; a missing null median is not a hit and B + 1 stays the
    denominator; a missing r gives `NA`), so the null medians never need
    a B × pairs matrix. Tested equal (`identical()`) to `perm_p()` on the
    null medians that `sparcc_details()` returns. **This is FastSpar's
    estimate exactly:** recomputing (1 + #{|r_b| ≥ |r|})/(B + 1) from the
    99 bootstrap correlation tables matches `fastspar_pvalues`' output on
    all 435 pairs (difference 0; with > instead of ≥, or k/B, it is 0.01).
    The README's "accounts for the possibility of repetitious data
    permutations" didn't show in these outputs. Against FastSpar's own
    p-values (B = 99, 20 draws, 10 exclusions, its own nulls): Spearman
    0.91, the planted pair at 0.01 in both, 92% agreement on p ≤ 0.05
    (21 pairs in both, 8 and 7 in one only), as random nulls allow.
11. **Medians.** Per pair over the draws' non-missing values (`NA` for a
    dropped component, `NaN` for a non-positive variance), by quickselect;
    the mean of the middle two for an even count. All missing gives `NA`,
    with a warning. Each feature's `variances` is the median of its
    variances where active, non-positive values included. Tested
    `identical()` to `apply(..., median, na.rm = TRUE)` on the details.
12. **Memory.** The draws' ρ are stored pair-major (`n_draws` × block)
    up to 256 MB. Beyond that the pairs are taken in blocks, and each
    later block draws again (same streams), computing t only for its pairs
    and ρ from each draw's variances, kept from the first pass. Tested
    `identical()` for blocks of 1, 7, 100 and 299 pairs, with and without
    permutations. Per thread: n × D + 2 × pairs doubles plus pairs bytes.
13. **Threads:** one task per draw (`threads_for()`, at most `n_draws`
    threads), then one per pair for the medians. Identical for `cpus` 1, 2,
    3, 4 and 8; Matrix, slam, `margin = 2`, double and data-frame input
    give identical results. Timing on this container (4 cores): 200
    samples × 500 features, 100 draws: 4.7 s on 1 thread, 1.5 s on 4; 100
    × 2,000, 20 draws: 3.7 s on 4 (4.4 s with forced 200,000-pair blocks).
14. **Input.** Own checks (`sparcc_input()`): numeric, finite,
    non-negative; a warning when not whole numbers; errors for fewer than
    4 features or 3 samples. The counts are copied densely (D × n) in C;
    there is no sparse path, as planned. `max_exclusions` above the number
    of pairs is treated as that number (its K × K scratch would otherwise
    be allocated at the requested size).
15. **Result.** As planned: `method = "sparcc"`, `n_samples`,
    `variances` and, with permutations, `p_value`; title "SparCC
    Correlation" from the registry's `name`. `variances` follows the
    labels: reordered by `clust` and named by `Labels` (positions for
    unnamed features). `proportionality()` now takes its titles from the
    registry too (same strings as before).
16. **Warnings:** non-positive variances (how many draws, how many
    values); pairs `NA` in every draw; medians beyond ±1. `ex_counts`
    (4 samples) warns about non-positive variances, so the examples
    simulate 60 samples instead.
17. **Shared code changed** (the plan's two small changes):
    - `dist_clust_order()` clusters any `kind = 'similarity'` registry
      method as `1 - x` (was: `method == 'rho_p'`), via `is_similarity()`;
    - `dist_plot_name()` treats any registry method as over features (was
      the `vlr`/`rho_p`/`phi_s` test). `plot.ecodive_dist()` now draws a
      similarity on the diverging scale, at least from −1 to 1, with the
      diagonal at 1 (it was blank, and `rho_p()`'s range came from the
      data);
    - `DIST_PARALLEL <- c('p_value', 'stability', 'n_pairs')` in
      `R/dist_clust.r`; `dist_reorder()` reorders those attributes with
      the triangle when their length matches it (tested: p-values stay with
      their pairs after `clust = 'average'`).
18. **`sparcc_details()`** (internal) returns the result plus per-draw ρ
    (`n_draws` × pairs), variances (`n_draws` × D), exclusions (packed,
    1-based, in order), `stop` per draw (`'max_exclusions'`,
    `'threshold'`, `'three'`, `'singular'`) and the null medians. It
    takes `log_fractions` (a list of n × D matrices, one per draw) and
    `block_pairs`.
19. **Docs:** `sparcc.Rd` (what is estimated, the assumptions and the
    paper's rule of thumb, zeros, "Why not rho_p()", P-values and their
    cost, Other implementations, Size, Same names); `rho_p.Rd` gained the
    matching paragraph and links to `sparcc()`; `plot.ecodive_dist.Rd`
    mentions `sparcc()`. pkgdown: "Compositional Data" → "Associations and
    networks". NEWS. `parameter_names.md` §6 gains `max_exclusions` and
    `exclusion_threshold`. No registry has an alias starting `spar`
    (tested; DA_TESTS doesn't exist yet).
20. **roxygen:** run with roxygen2 7.3.1 (the container's) when the
    docs were final, and once more after the reference comparisons
    changed `sparcc()`'s docs. It reformatted about 70 Rd files that were
    last written by 8.1.1, and added `RoxygenNote: 7.3.1` to
    `DESCRIPTION`; committed as written, per the README rule.
21. **FastSpar stops at D − 3 exclusions; ecodive and SpiecEasi don't.**
    FastSpar's output is identical for every `--exclude_iterations` from
    D − 3 up (D = 8, 10, 12, 20: from 5, 7, 9, 17), and it prints
    "Exclusion iterations ending early; we need at least 4 non-excluded
    pairs for SparCC". SpiecEasi keeps going and agrees with ecodive (D =
    8, 10 exclusions: 0.0013), while FastSpar then differs from both by up
    to 0.43. ecodive keeps the plan's rule (stop when fewer than 4
    components are left, which by 5 can't happen above rounding), so it
    matches SpiecEasi. The cap matters only for fewer than
    `max_exclusions + 3` features (22 at the defaults); the docs say to
    give `max_exclusions = D - 3` to reproduce FastSpar there. Tested both
    ways against the fixtures.

### What item 40 (`assoc_network()`) can rely on

- A `sparcc()` result is an `ecodive_dist` over features with `method =
  "sparcc"`, top-level `n_samples`, `variances` (named, in label order),
  and `p_value` (only with `permutations > 0`) parallel to the triangle,
  in `[1/(B+1), 1]` and on a grid of `1/(B+1)`; `NA` where the
  correlation is `NA`.
- Values are correlations, normally in [−1, 1]; medians beyond ±1 are
  possible (warned) and should be clamped or kept as the user prefers,
  e.g. before Fisher's z.
- `FEATURE_ASSOCS` (`R/feature_assoc.r`): `feature_assoc_info(id)` gives a
  row as a list; `is_similarity(d)` tests a result. `fisher_z` is `TRUE`
  for `sparcc` only. Add rows (and `fun`, if the wrapper needs it) for
  `spiec_easi_*` and `secom_*`; `alt_ids` is there for
  `match_feature_assoc()`.
- `dist_reorder()` keeps `p_value`, `stability` and `n_pairs` aligned
  through `clust`, but **not** per-feature attributes: `sparcc()` reorders
  its own `variances`.
- Re-running "the same association" on permuted samples (for
  `assoc_network_compare()`): `sparcc()` needs only `counts` and its
  arguments; its `seed` streams are by (dataset, draw, sample), so a
  caller wanting independent runs should vary `seed`.
