// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * capda()'s permutation test: the trace statistic of canonical analysis
 * of principal coordinates (Anderson & Willis 2003), the sum of the
 * squared canonical correlations between the first m PCoA axes and the
 * groups, for the observed groups and each permutation of them.
 *
 * With Q the n x m PCoA axes scaled to unit length (centered, orthonormal
 * columns) and H the projection onto the group means, the trace is
 * tr(Q' H Q) = sum over groups k of |S_k|^2 / n_k, where S_k is the sum
 * of the rows of Q in group k: O(n m) per permutation, with no n x n
 * matrix.
 *
 * Each permutation's trace is computed on its own, in a fixed order, so
 * results are the same for any number of threads.
 */

#include "ecodive.h"

#define CAP_BLOCK 16  // Permutations per task.

static const double *cp_Qt;     // m x n: Q transposed (rows of Q contiguous).
static const int    *cp_grp;    // n: groups, 0 .. g - 1.
static const double *cp_inv;    // g: 1 / n_k.
static perm_source_t cp_ps;
static int           cp_n, cp_m, cp_g, cp_blk0, cp_blk1, cp_total;
static double       *cp_work;   // Per thread: g x m sums.
static int          *cp_iwork;  // Per thread: n.
static double       *cp_out;    // n_perm + 1.


static void *trace_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = cp_n, m = cp_m, g = cp_g;

  double *S = cp_work  + (size_t) thread_i * g * m;
  int    *p = cp_iwork + (size_t) thread_i * n;

  for (int blk = cp_blk0 + thread_i; blk < cp_blk1; blk += n_threads) {

    int k1 = (blk + 1) * CAP_BLOCK;
    if (k1 > cp_total) k1 = cp_total;

    for (int k = blk * CAP_BLOCK; k < k1; k++) {

      // Sample a takes the group of sample p[a].
      perm_source_get(&cp_ps, k, p);

      memset(S, 0, (size_t) g * m * sizeof(double));
      for (int a = 0; a < n; a++) {
        double       *s = S + (size_t) cp_grp[p[a]] * m;
        const double *q = cp_Qt + (size_t) a * m;
        for (int j = 0; j < m; j++) s[j] += q[j];
      }

      double tr = 0;
      for (int c = 0; c < g; c++) {
        const double *s  = S + (size_t) c * m;
        double        ss = 0;
        for (int j = 0; j < m; j++) ss += s[j] * s[j];
        tr += ss * cp_inv[c];
      }
      cp_out[k] = tr;
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// Q:      n x m double matrix, orthonormal centered
//         columns (the PCoA axes scaled to unit length).
// grp:    each sample's group, 0 .. g - 1 (integer).
// g:      number of groups.
// perms:  n x n_perm integer matrix of 0-based
//         permutations, or NULL for random ones.
// strata: each sample's stratum, 0 .. n_strata - 1, or
//         NULL.
// opts:   c(n_perm, seed, n_threads).
//
// Returns the trace for the observed groups, then for
// each permutation (length n_perm + 1).
//======================================================
SEXP C_capda_trace(
    SEXP sexp_Q, SEXP sexp_grp, SEXP sexp_g, SEXP sexp_perms,
    SEXP sexp_strata, SEXP sexp_opts ) {

  int n         = nrows(sexp_Q);
  int m         = ncols(sexp_Q);
  int g         = asInteger(sexp_g);
  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];

  if (!isReal(sexp_Q) || !isInteger(sexp_grp) || XLENGTH(sexp_grp) != n ||
      n < 2 || m < 1 || g < 1 || n_perm < 0 || n_threads < 1)
    error("Invalid arguments to C_capda_trace()."); // # nocov

  const double *Q   = REAL(sexp_Q);
  const int    *grp = INTEGER(sexp_grp);

  int *cnt = (int *) R_alloc(g, sizeof(int));
  for (int c = 0; c < g; c++) cnt[c] = 0;
  for (int a = 0; a < n; a++)
    if (grp[a] < 0 || grp[a] >= g) error("Invalid arguments to C_capda_trace()."); // # nocov
    else cnt[grp[a]]++;

  double *inv = (double *) R_alloc(g, sizeof(double));
  for (int c = 0; c < g; c++) inv[c] = cnt[c] ? 1.0 / cnt[c] : 0;

  double *Qt = (double *) R_alloc((size_t) n * m, sizeof(double));
  for (int a = 0; a < n; a++)
    for (int j = 0; j < m; j++)
      Qt[(size_t) a * m + j] = Q[a + (size_t) j * n];

  perm_source_init(&cp_ps, n, sexp_perms, sexp_strata, seed);

  cp_Qt = Qt; cp_grp = grp; cp_inv = inv;
  cp_n = n; cp_m = m; cp_g = g;
  cp_total = n_perm + 1;

  int n_blocks = (int) (((size_t) cp_total + CAP_BLOCK - 1) / CAP_BLOCK);
  if (n_threads > n_blocks) n_threads = n_blocks;

  cp_work  = (double *) R_alloc((size_t) g * m * n_threads, sizeof(double));
  cp_iwork = (int *)    R_alloc((size_t) n * n_threads,     sizeof(int));

  SEXP sexp_out = PROTECT(allocVector(REALSXP, cp_total));
  cp_out = REAL(sexp_out);

  // Rounds of blocks, checking for interrupts between
  // them: each round has at least ~1e8 operations per
  // thread, or every remaining block.
  double blk_work  = ((double) n * (m + 2) + (double) g * m) * CAP_BLOCK + 1;
  double rounds    = 1e8 * n_threads / blk_work;
  int    per_round = (rounds < n_blocks) ? (int) rounds : n_blocks; // no overflow
  if (per_round < n_threads) per_round = n_threads;

  for (cp_blk0 = 0; cp_blk0 < n_blocks; cp_blk0 = cp_blk1) {

    cp_blk1 = (n_blocks - cp_blk0 > per_round) ? cp_blk0 + per_round : n_blocks;

    int nb = cp_blk1 - cp_blk0;
    int nt = threads_for(nb < n_threads ? nb : n_threads, nb * blk_work, 1e6);

    // run_parallel() starts threads only when told of at
    // least 100 tasks.
    run_parallel(trace_worker, nt, nt > 1 ? 100 : 1);

    if (cp_blk1 < n_blocks) R_CheckUserInterrupt();
  }

  UNPROTECT(1);
  return sexp_out;
}
