# documentation

documentation

## Pairwise comparisons

With `pairwise = TRUE`, each pair of groups is also tested on its own,
as q2-diversity's `beta_group_significance(pairwise=True)` does: the
distances are cut down to the two groups' samples, and the same test is
run on them, with the same number of `permutations`, the same `strata`
(cut down too), and the test's other arguments. The overall test is
unchanged.

- **Permutations:** Each pair has random permutations of its own, drawn
  from `seed` and the pair's position in the table, never from the
  number of `cpus`. So a pair's p-value is reproducible, but differs (by
  sampling error) from that of the same test run by hand on the pair's
  samples with the same `seed`. A supplied `permutations` matrix is cut
  down to each pair's samples by following each permutation's cycles:
  each of the pair's samples goes to the next of the pair's samples on
  its cycle. That gives a permutation of the pair's samples, uniformly
  distributed when the full one is, which keeps any restriction the
  matrix has, such as permuting only within strata.

- **Multiple testing:** The pairs' permutation p-values are adjusted
  together with `stats::p.adjust(method = p_adjust)`. The default,
  Holm's method, controls the family-wise error rate whatever the
  dependence between the pairs' tests. `'BH'` controls the false
  discovery rate instead, as q2-diversity's q-values do. The raw
  p-values are kept as well.

- **Pairs that can't be tested:** A pair the test can't be run on, such
  as two groups of one sample each, gets `NA`s, with a warning, rather
  than an error. A test's other warnings for a pair, such as all of its
  distances being 0, are passed on, naming the pair.

The `pairwise` field of the result is a data frame with one row per pair
of groups, named like `"a-b"`, in the order of the groups' factor levels
(sorted values, for other vectors). Its columns are the test's
statistics for the pair, `z`, `n` (the pair's number of samples),
`n_perm`, `p_value` (unadjusted) and `p_adjusted`. The `p_adjust` field
names the method, and [`print()`](https://rdrr.io/r/base/print.html)
shows both.
