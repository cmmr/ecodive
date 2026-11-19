# `ctf()`: compositional tensor factorization

**Priority:** 31. **Tier:** port (gemelli 0.0.13,
`ctf()`, BSD-3). The largest of the gemelli ports, with several apparent
bugs in the reference to decide on first. Read [README.md](README.md)
first, and [compositional.md](compositional.md) for the conventions the
compositional plans share.


## What it is

CTF (Martino et al. 2021) is RPCA for **repeated measures**: samples are
arranged in a tensor **subjects × features × states** (time points or
conditions), robust-clr transformed, and factorised by a CP decomposition
that treats zeros as missing. It separates differences **between
subjects** (a subject ordination with feature loadings) from **changes
over states** (state loadings), and gives per-subject trajectories, which
an ordinary ordination of all samples confounds. It's QIIME 2's
`gemelli ctf`.

**Name:** `ctf()`, the method's established acronym. Not mistakable for
another ecodive function (`ca`, `cca` differ). Docs aliases: "compositional
tensor factorization", "gemelli CTF". **No wrapper:** it needs subjects and
states besides the table, and returns subject-level results.


## Decide before writing code (ask the user)

Reading gemelli's source (not running it) found behaviour that looks
wrong. **Confirm each by running gemelli 0.0.13 in a sandbox** on a small
table where the right answer is known, then ask the user, per item,
whether ecodive should replicate it (parity, documented) or fix it
(documented difference; fixtures then compare only what isn't affected):

1. **Label misalignment:** the subject and feature factors are row-sorted
   by their first column and then labelled in the **original** order, so
   biplot coordinates may be attached to the wrong subjects/features.
   (gemelli's tests check only types and the absence of NaN.)
2. **Eigenvalue pairing:** eigenvalues sorted descending but factor
   columns ascending before \eqn{X = \text{subjects}\, S\, \text{features}^T}.
3. **Proportion explained** is \eqn{s / \sum s}, not \eqn{s^2 / \sum s^2}.
4. **State loadings** labelled PC1..PCk in ascending eigenvalue order.
5. **Power-method "deflation":** each start is deflated by
   (Σ eigenvalues) × its **own** outer product (numpy broadcasting), not
   by the components found so far.
6. **Duplicate (subject, state) samples:** the first sample's column is
   replaced by the mean **truncated to an integer**; the other duplicates
   stay in the table (non-in-place `drop`), only the first enters the
   tensor.
7. **Filters vs metadata:** if a filter drops a sample still in the
   metadata, tensor construction likely crashes.
8. **Phylogenetic CTF** (not planned here): branch lengths keep their
   original order after features are reordered.

Recommendation: **fix** 1–4 and 6 (they change what users read off the
plot), document each under "Differences from gemelli", and match gemelli
exactly on everything upstream of them (the tensor, the rclr, the
factorisation's eigenvalues and factor columns as sets). Replicate 5 (it
defines the algorithm's starting point and the results depend on it).
Handle 7 by filtering metadata with the table.


## Design

### Signature

```r
ctf(counts, subjects, states, data = NULL, margin = 1L, n_axes = 3L,
    max_iters = 25L, max_power_iters = 25L, n_starts = 25L, tol = 1e-7,
    min_sample_count = 0, min_feature_count = 0, min_feature_frequency = 0,
    cpus = n_cpus())
```

- `subjects`, `states`: per-sample vectors or column names of `data`
  (gemelli's `individual_id_column`, `state_column`), handled as `groups`
  is (`per_sample()`). One state variable only (gemelli builds a 3-way
  tensor).
- `max_iters` (gemelli's `max_iterations_als`), `max_power_iters`
  (`max_iterations_rptm`), `n_starts` (`n_initializations`; the glossary's
  name for random starts, as in `mds()` / `nmds()`), `tol` (`tol_als`).
  gemelli's fixed `tol_rtpm = 1e-5` stays internal.
- Filters as `rpca()`.
- gemelli errors if the table has no zeros at all; match it, with a clear
  message.

### Algorithm (gemelli 0.0.13, as read)

1. Filter (as `rpca()`); restrict metadata to the kept samples.
2. Features ordered by descending total (ties by ascending id). rclr per
   sample over the whole table; non-finite → 0 (0 means missing in the
   tensor).
3. Subjects and states ordered by **ascending summed sequencing depth**
   (not time order; gemelli's sort isn't stable: use a stable sort and
   document ties).
4. **Initialisation (modified robust tensor power method):** factors
   start as ones, eigenvalues as ones. For each component: re-seed the
   generator with 42; draw `n_starts` random starts (normalised); power
   iterations update each mode from the other modes' previous vectors
   until the summed change < 1e-5 or `max_power_iters`; keep the start
   with the largest |eigenvalue|.
5. **Masked ALS:** for each component c, set its eigenvalue to 0, form
   the masked residual without c, update each mode as (residual
   contracted with the other modes) / (mask contracted with their
   squares) (zero denominators → 1); normalise only the last mode, whose
   norm becomes the eigenvalue. Stop when
   ‖T − mask ⊙ reconstruction‖ / ‖T‖ < `tol` or `max_iters`.
6. Subject × feature matrix from the factors, double-centred and SVD'd:
   subject scores, feature loadings, eigenvalues.
7. **Trajectories:** per component, subject loading × state loading, one
   row per (subject, state), centred; sample distances are Euclidean
   between those points, cut to the real samples. Feature trajectories
   likewise.

### Random numbers: an exception to README's PCG rule

gemelli's starts come from numpy's `RandomState(42)` (MT19937, doubles as
`((a >> 5) * 2^26 + (b >> 6)) / 2^53`), re-seeded per component, drawn
per start in the order subjects, features, states. To match it, port
**MT19937** (public-domain reference code by Matsumoto & Nishimura) for
this function only, fixed at seed 42 as gemelli is (no `seed` argument;
or `seed = 42L` if the user prefers one, defaulting to gemelli's). The
research verified a pure-Python MT19937 against numpy over 1,000 draws.
Note this exception in README.md once agreed.

### C

- **Never form the dense tensor:** every contraction runs over observed
  entries, because the subtracted terms are rank-1 and their contractions
  factorise. Parallel over observed entries in fixed blocks, reduced in
  block order.
- Reuse the QR/SVD code from `rpca()` for step 6.
- Size estimate: 900–1,300 C lines, 350–450 R lines.

### Result

`ecodive_ord` via `ord_result()` with **subjects as rows** (prefix
`"PC"`, as gemelli's biplot; `method = 'ctf'`), `eig`, `prop` (decided
under item 3), plus:

- `loadings`: features × axes.
- `state_loadings`: states × axes.
- `trajectories`: a data frame, one row per sample (subject, state), with
  `subject`, `state` and one column per axis.
- `feature_trajectories`: as gemelli's feature-trajectory table.
- `dm`: a `dist` between samples (step 7).
- `n_iters`.

Out of scope: `phylogenetic_ctf` (later, as `ctf_phylo()`, sharing
[rpca_phylo.md](rpca_phylo.md)'s tree expansion), and TEMPTED (gemelli's
temporal tensor decomposition with smooth time loadings; a large separate
port).


## Reference implementations

- gemelli 0.0.13: `ctf()` (`gemelli/ctf.py`), `build_sparse` /
  `tensor_rclr` (`gemelli/preprocessing.py`), `TensorFactorization`
  (`gemelli/factorization.py`), BSD-3. `khatri_rao` there is "modified from
  TensorLy" (check TensorLy's license, BSD-3 at the time of writing,
  before porting it; it's a one-liner to rewrite).
- **Fixtures** (`data-raw/expected/ctf.py`): a small longitudinal table
  (e.g. 8 subjects × 4 states × 50 features, some missing samples), plus
  the sandbox runs that settle the "Decide" list. Compare eigenvalues and
  factor columns (as sets) at about 1e-8 (iterative), and the fixed outputs
  per the decisions.


## Tests

- Fixtures; a missing (subject, state) cell; duplicate samples (per the
  decision); `n_axes = 2` (gemelli pads a third axis of zeros in some
  outputs: match or document).
- The MT19937 stream matches numpy's `RandomState(42).random_sample()`
  (store 1,000 values).
- Identical for `cpus` 1, 2 and 4.
- Errors: unknown columns; a single subject or state; a table with no
  zeros (gemelli's error).
- A smoke-test line.


## Docs

- Rd page under Ordination → "Several tables or time series", with
  `rpca_joint()`; a "Differences from gemelli" section listing every
  decision above.
- NEWS.


## References

- Martino, C., Shenhav, L., Marotz, C. A., et al. (2021). Context-aware
  dimensionality reduction deconvolutes gut microbial community dynamics.
  *Nature Biotechnology*, 39, 165–168.
- Anandkumar, A., Ge, R., Hsu, D., Kakade, S. M., & Telgarsky, M. (2014).
  Tensor decompositions for learning latent variable models. *Journal of
  Machine Learning Research*, 15, 2773–2832.
- Matsumoto, M., & Nishimura, T. (1998). Mersenne twister. *ACM
  Transactions on Modeling and Computer Simulation*, 8, 3–30.
