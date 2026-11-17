// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Multi-response permutation procedure (MRPP): the weighted
 * mean within-group distance, delta, for the observed groups
 * and for each permutation of them.
 *
 * Group g, with n_g samples, has P_g = n_g (n_g - 1) / 2
 * pairs, and mean within-group distance
 *
 *     delta_g = S_g / P_g,   S_g = sum of d_ij over its pairs.
 *
 * With weights w_g (summing to W over the groups with at
 * least two samples),
 *
 *     delta = sum_g (w_g / W) delta_g = sum_g c_g S_g,
 *
 * where c_g = w_g / (W P_g) is computed in R. So each
 * permutation costs one comparison per pair: add d_ij to
 * pair (i, j)'s sum when i and j are in the same group, and
 * scale each column's sum by its group's c_g. That is the
 * kernel of permanova()'s one-factor path, with d_ij in place
 * of d_ij^2.
 *
 * Permutations: as vegan::mrpp(), which takes grouping[p],
 * permuted sample a is in group g[p[a]].
 *
 * Determinism: permutation k (k >= 1) comes from its own PCG
 * stream, seeded with (seed, k), or from a supplied matrix.
 * Permutation 0 is the identity: the observed delta. The
 * permutations are split into fixed blocks of MRPP_FW; each
 * block is computed by one thread, and every sum runs in a
 * fixed order. So results are bit-identical for any number
 * of threads.
 *
 * The permutations (permute.c) and the group kernel
 * (group_kernel.h) are shared with permanova().
 */

#include "ecodive.h"
#include "group_kernel.h"

#define MRPP_W  GK_W           // Lanes per chunk.
#define MRPP_FW (4 * MRPP_W)   // Permutations per block.


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the parts of the output that each block owns.
//------------------------------------------------------
static int           mr_n;       // Samples.
static const double *mr_d;       // Packed lower triangle of D.
static const int    *mr_grp;     // Each sample's group, 0 .. g - 1.
static const double *mr_coef;    // c_g for each group.
static perm_source_t mr_ps;      // Where permutations come from.
static int           mr_blk0;    // Blocks for this round:
static int           mr_blk1;    // ... mr_blk0 to mr_blk1 - 1.
static int           mr_total;   // Permutations, including the identity.
static double       *mr_out;     // mr_total deltas.
static double       *mr_work;    // Per-thread workspace ...
static size_t        mr_wlen;    // ... of this many doubles.
static int          *mr_iwork;   // Per-thread integer workspace (n each).


//------------------------------------------------------
// Each task is one block of MRPP_FW permutations. L holds
// each sample's group (as a double, so comparisons
// vectorize alongside the sums) under each permutation:
// chunk h is L + h * n * MRPP_W, row-major (n x MRPP_W).
// Unused lanes are group 0, and ignored.
//------------------------------------------------------
static void *mrpp_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = mr_n;
  int nc        = MRPP_FW / MRPP_W;

  double *L = mr_work + (size_t) thread_i * mr_wlen;     // nc x n x W
  double *s = L + (size_t) n * MRPP_FW;                  // MRPP_FW
  int    *p = mr_iwork + (size_t) thread_i * n;          // n

  for (int blk = mr_blk0 + thread_i; blk < mr_blk1; blk += n_threads) {

    int k0 = blk * MRPP_FW;
    int nb = mr_total - k0;
    if (nb > MRPP_FW) nb = MRPP_FW;

    memset(L, 0, (size_t) n * MRPP_FW * sizeof(double));

    for (int b = 0; b < nb; b++) {
      double *lc = L + (size_t) (b / MRPP_W) * n * MRPP_W + b % MRPP_W;
      perm_source_get(&mr_ps, k0 + b, p);
      for (int a = 0; a < n; a++)
        lc[(size_t) a * MRPP_W] = mr_grp[p[a]];
    }

    for (int b = 0; b < MRPP_FW; b++) s[b] = 0;

    for (int j = 0; j < n - 1; ) {

      int nd = n - j - 1;
      const double *d0 = mr_d + dist_col_start(j, n);

      if (nd >= 2) {
        const double *d1 = mr_d + dist_col_start(j + 1, n);
        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * MRPP_W;
          group_chunk2(
            d0, d1, vc + (size_t) (j + 1) * MRPP_W,
            vc + (size_t) j * MRPP_W, vc + (size_t) (j + 1) * MRPP_W,
            nd, s + h * MRPP_W, mr_coef );
        }
        j += 2;
      }
      else {
        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * MRPP_W;
          group_chunk(
            d0, vc + (size_t) (j + 1) * MRPP_W, vc + (size_t) j * MRPP_W,
            nd, s + h * MRPP_W, mr_coef );
        }
        j += 1;
      }
    }

    for (int b = 0; b < nb; b++)
      mr_out[k0 + b] = s[b];
  }

  return NULL;
}


//======================================================
// R interface.
//
// dm:     packed lower triangle of D (double), checked.
// n:      number of samples.
// grp:    each sample's group, 0 .. g - 1 (integer).
// coef:   c_g = w_g / (W P_g) for each group (double); 0
//         for groups left out of delta.
// perms:  n x n_perm integer matrix of 0-based permutations,
//         or NULL for random ones.
// strata: each sample's stratum, 0 .. n_strata - 1, or NULL.
// opts:   c(n_perm, seed, n_threads).
//
// Returns list(delta, class_sum): delta has n_perm + 1
// values, the observed one first; class_sum is S_g for the
// observed groups.
//======================================================
SEXP C_mrpp(
    SEXP sexp_dm, SEXP sexp_n, SEXP sexp_grp, SEXP sexp_coef,
    SEXP sexp_perms, SEXP sexp_strata, SEXP sexp_opts ) {

  int n         = asInteger(sexp_n);
  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];
  int g         = LENGTH(sexp_coef);

  if (n < 2 || n_perm < 0 || n_threads < 1 || g < 1 ||
      XLENGTH(sexp_dm) != (R_xlen_t) n * (n - 1) / 2 ||
      LENGTH(sexp_grp) != n)
    error("Invalid arguments to C_mrpp()."); // # nocov

  const int *grp = INTEGER(sexp_grp);
  for (int a = 0; a < n; a++)
    if (grp[a] < 0 || grp[a] >= g)
      error("Invalid arguments to C_mrpp()."); // # nocov

  mr_n     = n;
  mr_d     = REAL(sexp_dm);
  mr_grp   = grp;
  mr_coef  = REAL(sexp_coef);
  mr_total = n_perm + 1;


  //--------------------------------------------------
  // Permutations: sample a takes group grp[p[a]], as
  // vegan's grouping[p].
  //--------------------------------------------------
  perm_source_init(&mr_ps, n, sexp_perms, sexp_strata, seed);


  //--------------------------------------------------
  // Observed within-group sums, S_g.
  //--------------------------------------------------
  SEXP sexp_cs = PROTECT(allocVector(REALSXP, g));
  double *cs = REAL(sexp_cs);
  for (int c = 0; c < g; c++) cs[c] = 0;
  for (int j = 0; j < n - 1; j++) {
    const double *d  = mr_d + dist_col_start(j, n);
    int           gj = grp[j];
    for (int i = j + 1; i < n; i++)
      if (grp[i] == gj) cs[gj] += d[i - j - 1];
  }


  //--------------------------------------------------
  // Blocks, and workspace for each thread.
  //--------------------------------------------------
  int n_blocks = (int) (((size_t) mr_total + MRPP_FW - 1) / MRPP_FW);
  if (n_threads > n_blocks) n_threads = n_blocks;

  mr_wlen  = (size_t) n * MRPP_FW + MRPP_FW;
  mr_work  = (double *) R_alloc(mr_wlen * n_threads, sizeof(double));
  mr_iwork = (int *)    R_alloc((size_t) n * n_threads, sizeof(int));

  SEXP sexp_delta = PROTECT(allocVector(REALSXP, mr_total));
  mr_out = REAL(sexp_delta);


  //--------------------------------------------------
  // Rounds of blocks, checking for interrupts between
  // them. Each round has at least ~1e8 comparisons per
  // thread, or every remaining block.
  //--------------------------------------------------
  double pairs     = (double) n * (n - 1) / 2;
  double blk_work  = pairs * MRPP_FW + 1;
  double rounds    = 1e8 * n_threads / blk_work;
  int    per_round = (rounds < n_blocks) ? (int) rounds : n_blocks; // no overflow
  if (per_round < n_threads) per_round = n_threads;

  // run_parallel() starts threads only when told of at
  // least 100 tasks.
  int tasks_hint = (blk_work > 1e5) ? 100 : 1;

  for (mr_blk0 = 0; mr_blk0 < n_blocks; mr_blk0 = mr_blk1) {

    mr_blk1 = (n_blocks - mr_blk0 > per_round) ? mr_blk0 + per_round : n_blocks;

    int nt = (mr_blk1 - mr_blk0 < n_threads) ? mr_blk1 - mr_blk0 : n_threads;
    run_parallel(mrpp_worker, nt, (mr_blk1 - mr_blk0 > 1) ? tasks_hint : 1);

    if (mr_blk1 < n_blocks) R_CheckUserInterrupt();
  }

  const char *names[] = { "delta", "class_sum", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_delta);
  SET_VECTOR_ELT(sexp_result, 1, sexp_cs);

  UNPROTECT(3);
  return sexp_result;
}
