// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Sample permutations for the permutation tests (permanova(),
 * mrpp(), ...): either supplied, as a matrix, or random,
 * shuffled freely within strata.
 *
 * Permutation k (k >= 1) is drawn from its own PCG stream,
 * seeded with (seed, k), so permutations can be generated in
 * any order, by any thread, with the same results.
 * Permutation 0 is the identity: the observed statistics.
 */

#include "ecodive.h"


//======================================================
// Sets up `ps` for n samples. All memory is R_alloc()'d.
//
// sexp_perms:  n x n_perm integer matrix of 0-based
//              permutations (one per column), or NULL for
//              random ones.
// sexp_strata: each sample's stratum, 0 .. n_strata - 1
//              (integer), or NULL for free permutation.
// seed:        the seed for random permutations.
//======================================================
void perm_source_init(
    perm_source_t *ps, int n, SEXP sexp_perms, SEXP sexp_strata, int seed ) {

  ps->n     = n;
  ps->perms = isNull(sexp_perms) ? NULL : INTEGER(sexp_perms);
  ps->seed  = (uint64_t) seed;

  int *order = (int *) R_alloc(n, sizeof(int));
  int  ns    = 1;
  int *sbeg;

  if (isNull(sexp_strata)) {
    for (int a = 0; a < n; a++) order[a] = a;
    sbeg    = (int *) R_alloc(2, sizeof(int));
    sbeg[0] = 0;
    sbeg[1] = n;
  }
  else {
    const int *strata = INTEGER(sexp_strata);
    ns = 0;
    for (int a = 0; a < n; a++)
      if (strata[a] < 0 || strata[a] >= n) error("Invalid strata codes."); // # nocov
      else if (strata[a] + 1 > ns) ns = strata[a] + 1;
    sbeg = (int *) R_alloc(ns + 1, sizeof(int));
    for (int s = 0; s <= ns; s++) sbeg[s] = 0;
    for (int a = 0; a < n; a++) sbeg[strata[a] + 1]++;
    for (int s = 0; s < ns; s++) sbeg[s + 1] += sbeg[s];
    int *fill = (int *) R_alloc(ns, sizeof(int));
    memcpy(fill, sbeg, (size_t) ns * sizeof(int));
    for (int a = 0; a < n; a++) order[fill[strata[a]]++] = a;
  }

  ps->order = order;
  ps->sbeg  = sbeg;
  ps->ns    = ns;
}


//======================================================
// Permutation k into p (length n). k = 0 is the identity.
// Safe to call from worker threads: no R API, and `ps` is
// only read.
//
// The meaning of p is the caller's: permanova() scatters
// (permuted row a is original row p[a], as vegan's G[p, p]),
// mrpp() gathers (sample a takes group g[p[a]], as vegan's
// grouping[p]). A supplied matrix is used as given.
//======================================================
void perm_source_get(const perm_source_t *ps, int k, int *p) {

  int n = ps->n;

  if (k == 0) {
    for (int a = 0; a < n; a++) p[a] = a;
    return;
  }

  if (ps->perms) {
    memcpy(p, ps->perms + (size_t) (k - 1) * n, (size_t) n * sizeof(int));
    return;
  }

  // Fisher-Yates shuffle within each stratum.
  pcg32_random_t rng;
  pcg32_srandom_r(&rng, ps->seed, (uint64_t) k);

  for (int a = 0; a < n; a++) p[a] = a;

  for (int s = 0; s < ps->ns; s++) {
    const int *idx = ps->order + ps->sbeg[s];
    int        len = ps->sbeg[s + 1] - ps->sbeg[s];
    for (int t = len - 1; t > 0; t--) {
      int u   = (int) pcg32_bounded(&rng, (uint32_t) t + 1);
      int tmp = p[idx[t]];
      p[idx[t]] = p[idx[u]];
      p[idx[u]] = tmp;
    }
  }
}


//======================================================
// The seed for pair `pair` (1, 2, ...) of a pairwise test
// (pairwise = TRUE), from the test's own seed: a 31-bit
// draw from PCG, seeded with (seed ^ a constant, pair).
// The constant keeps these draws apart from the streams of
// any test's own permutations, whose states are seed
// values sign-extended to 64 bits.
//======================================================
SEXP C_pair_seed(SEXP sexp_seed, SEXP sexp_pair) {

  pcg32_random_t rng;
  uint64_t       state = (uint64_t) asInteger(sexp_seed) ^ 0x9E3779B97F4A7C15ULL;

  pcg32_srandom_r(&rng, state, (uint64_t) asInteger(sexp_pair));

  return ScalarInteger((int) (pcg32_random_r(&rng) >> 1));
}


//======================================================
// Permutations 1 .. n_perm of n samples, as an n x n_perm
// integer matrix (1-based, so R can index with them):
// column k is perm_source_get()'s permutation k, the same
// one the C tests use for (seed, k), with `strata`
// (0-based codes, or NULL). For tests whose statistic is
// computed in R.
//======================================================
SEXP C_perm_indices(SEXP sexp_n, SEXP sexp_n_perm, SEXP sexp_strata, SEXP sexp_seed) {

  int n      = asInteger(sexp_n);
  int n_perm = asInteger(sexp_n_perm);
  if (n < 1 || n_perm < 0) error("Invalid arguments to C_perm_indices()."); // # nocov

  perm_source_t ps;
  perm_source_init(&ps, n, R_NilValue, sexp_strata, asInteger(sexp_seed));

  SEXP res = PROTECT(allocMatrix(INTSXP, n, n_perm));
  int *p   = INTEGER(res);
  for (int k = 1; k <= n_perm; k++) {
    int *col = p + (size_t) (k - 1) * n;
    perm_source_get(&ps, k, col);
    for (int a = 0; a < n; a++) col[a]++;
  }

  UNPROTECT(1);
  return res;
}
