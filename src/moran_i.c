// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// moran_i(): the spatial weights and the permutation test of Moran's I.
//
// The weight transformations (row-standardised, binary, variance-
// stabilising) and the sums s0, s1 and s2 follow libpysal 4.15.1's
// W.set_transform(), W.s0, W.s1 and W.s2 (libpysal/weights/weights.py),
// and the statistic follows esda 2.10.0's Moran (esda/moran.py),
// restructured for ecodive: the weights are stored as the packed lower
// triangle of their symmetric part, which is all that I needs. Those
// sources are under this license:
//
//   libpysal: Copyright (c) 2007-, PySAL Developers.
//   esda: Copyright 2018 PySAL Developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//
//   3. Neither the name of the copyright holder nor the names of its
//      contributors may be used to endorse or promote products derived
//      from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*
 * Moran's I of centred values z under weights W is
 *
 *   I = n / s0 * z' W z / z' z,
 *
 * and z' W z = z' S z with S = (W + W') / 2, so only S's packed lower
 * triangle is kept (the diagonal is ignored: there are no self-weights).
 * Each permutation's z' S z is one pass over the triangle, summed in a
 * fixed order, so results are the same for any number of threads.
 */

#include "ecodive.h"

#define MORAN_BLOCK 8  // Permutations per task.


//======================================================
// The weights.
//
// w:      the raw weights: an n x n matrix (any diagonal,
//         which is ignored), or, with `inverse`, packed
//         distances (a dist's lower triangle, all > 0)
//         whose weights are 1 / d.
// style:  0 row-standardised, 1 binary, 2 variance-
//         stabilising, 3 none.
//
// Returns list(S, s = c(s0, s1, s2), islands): S the
// packed lower triangle of (W + W') / 2 after the
// transformation; s0 = sum(W), s1 = sum((W + W')^2) / 2,
// s2 = sum((rowSums(W) + colSums(W))^2), as libpysal's;
// and the number of samples whose raw weights are all 0.
//======================================================

static R_INLINE double raw_w(const double *w, int inverse, int n, size_t pk, int i, int j) {
  return inverse ? 1.0 / w[pk] : w[i + (size_t) j * n];
}

SEXP C_moran_weights(SEXP sexp_w, SEXP sexp_n, SEXP sexp_inverse, SEXP sexp_style) {

  int n       = asInteger(sexp_n);
  int inverse = asLogical(sexp_inverse);
  int style   = asInteger(sexp_style);

  if (!isReal(sexp_w) || n < 2 || style < 0 || style > 3 ||
      XLENGTH(sexp_w) != (inverse ? (R_xlen_t) n * (n - 1) / 2 : (R_xlen_t) n * n))
    error("Invalid arguments to C_moran_weights()."); // # nocov

  const double *w = REAL(sexp_w);

  // Per sample: raw row sum (of the transformed-to-binary
  // weights for style 1), raw sum of squares, then the
  // row's scale a, and the transformed row and column sums.
  double *rs = (double *) R_alloc(n, sizeof(double));
  double *sq = (double *) R_alloc(n, sizeof(double));
  double *a  = (double *) R_alloc(n, sizeof(double));
  double *rc = (double *) R_alloc(n, sizeof(double));
  int    *nz = (int *)    R_alloc(n, sizeof(int));
  for (int i = 0; i < n; i++) { rs[i] = 0; sq[i] = 0; rc[i] = 0; nz[i] = 0; }


  //--------------------------------------------------------
  // Pass 1: row sums, in pair order (column j, rows i > j).
  //--------------------------------------------------------
  size_t pk = 0;
  for (int j = 0; j < n - 1; j++)
    for (int i = j + 1; i < n; i++, pk++) {
      double u = raw_w(w, inverse, n, pk, i, j);   // w_ij
      double v = raw_w(w, inverse, n, pk, j, i);   // w_ji
      if (!R_FINITE(u) || !R_FINITE(v) || u < 0 || v < 0)
        error("The spatial weights must be finite and not negative.");
      if (u != 0) nz[i] = 1;
      if (v != 0) nz[j] = 1;
      if (style == 1) { u = (u != 0); v = (v != 0); }
      rs[i] += u; sq[i] += u * u;
      rs[j] += v; sq[j] += v * v;
    }

  int islands = 0;
  for (int i = 0; i < n; i++) islands += !nz[i];

  double g = 1;
  for (int i = 0; i < n; i++)
    switch (style) {
      case 0:  a[i] = rs[i] > 0 ? 1 / rs[i] : 0;       break;
      case 2:  a[i] = sq[i] > 0 ? 1 / sqrt(sq[i]) : 0; break;
      default: a[i] = 1;
    }

  // Variance-stabilising: rows scaled to unit length, then
  // all scaled so that the weights sum to n.
  if (style == 2) {
    double big_q = 0;
    for (int i = 0; i < n; i++) big_q += rs[i] * a[i];
    g = big_q > 0 ? n / big_q : 0;
  }


  //--------------------------------------------------------
  // Pass 2: the transformed weights' symmetric part and
  // their sums.
  //--------------------------------------------------------
  R_xlen_t n_pairs = (R_xlen_t) n * (n - 1) / 2;
  SEXP sexp_S = PROTECT(allocVector(REALSXP, n_pairs));
  double *S = REAL(sexp_S);

  double s0 = 0, s1 = 0;
  pk = 0;
  for (int j = 0; j < n - 1; j++)
    for (int i = j + 1; i < n; i++, pk++) {
      double u = raw_w(w, inverse, n, pk, i, j);
      double v = raw_w(w, inverse, n, pk, j, i);
      if (style == 1) { u = (u != 0); v = (v != 0); }
      u *= a[i] * g;
      v *= a[j] * g;
      double t = u + v;
      S[pk] = t / 2;
      s0 += t;
      s1 += t * t;
      rc[i] += t;   // row i's sum plus column i's sum
      rc[j] += t;
    }

  double s2 = 0;
  for (int i = 0; i < n; i++) s2 += rc[i] * rc[i];

  SEXP sexp_s = PROTECT(allocVector(REALSXP, 3));
  REAL(sexp_s)[0] = s0;
  REAL(sexp_s)[1] = s1;
  REAL(sexp_s)[2] = s2;

  SEXP res = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, sexp_S);
  SET_VECTOR_ELT(res, 1, sexp_s);
  SET_VECTOR_ELT(res, 2, ScalarInteger(islands));

  SEXP nms = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(nms, 0, mkChar("S"));
  SET_STRING_ELT(nms, 1, mkChar("s"));
  SET_STRING_ELT(nms, 2, mkChar("islands"));
  setAttrib(res, R_NamesSymbol, nms);

  UNPROTECT(4);
  return res;
}



//======================================================
// The permutation test: z' S z (S's diagonal 0), with z
// permuted, for the identity and each permutation.
//======================================================

static const double *mi_S;      // Packed lower triangle.
static const double *mi_z;      // n: centred values.
static perm_source_t mi_ps;
static int           mi_n, mi_blk0, mi_blk1, mi_total;
static double       *mi_work;   // Per thread: n.
static int          *mi_iwork;  // Per thread: n.
static double       *mi_out;    // n_perm + 1.


// sum over i > j of S_ij z_i z_j, column by column.
static double moran_quad(const double *S, const double *z, int n) {
  double total = 0;
  for (int j = 0; j < n - 1; j++) {
    double acc = 0;
    for (int i = j + 1; i < n; i++) acc += S[i - j - 1] * z[i];
    total += z[j] * acc;
    S += n - j - 1;
  }
  return 2 * total;
}


static void *moran_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = mi_n;

  double *zp = mi_work  + (size_t) thread_i * n;
  int    *p  = mi_iwork + (size_t) thread_i * n;

  for (int blk = mi_blk0 + thread_i; blk < mi_blk1; blk += n_threads) {

    int k1 = (blk + 1) * MORAN_BLOCK;
    if (k1 > mi_total) k1 = mi_total;

    for (int k = blk * MORAN_BLOCK; k < k1; k++) {
      // Sample a takes the value of sample p[a].
      perm_source_get(&mi_ps, k, p);
      for (int a = 0; a < n; a++) zp[a] = mi_z[p[a]];
      mi_out[k] = moran_quad(mi_S, zp, n);
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// S:      packed lower triangle (C_moran_weights()'s S).
// z:      n centred values.
// perms:  n x n_perm integer matrix of 0-based
//         permutations, or NULL for random ones.
// opts:   c(n_perm, seed, n_threads).
//
// Returns z' S z for the observed values, then for each
// permutation (length n_perm + 1).
//======================================================
SEXP C_moran_perm(SEXP sexp_S, SEXP sexp_z, SEXP sexp_perms, SEXP sexp_opts) {

  int n         = length(sexp_z);
  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];

  if (!isReal(sexp_S) || !isReal(sexp_z) || n < 2 || n_perm < 0 || n_threads < 1 ||
      XLENGTH(sexp_S) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_moran_perm()."); // # nocov

  perm_source_init(&mi_ps, n, sexp_perms, R_NilValue, seed);

  mi_S = REAL(sexp_S); mi_z = REAL(sexp_z);
  mi_n = n;
  mi_total = n_perm + 1;

  int n_blocks = (int) (((size_t) mi_total + MORAN_BLOCK - 1) / MORAN_BLOCK);
  if (n_threads > n_blocks) n_threads = n_blocks;

  mi_work  = (double *) R_alloc((size_t) n * n_threads, sizeof(double));
  mi_iwork = (int *)    R_alloc((size_t) n * n_threads, sizeof(int));

  SEXP sexp_out = PROTECT(allocVector(REALSXP, mi_total));
  mi_out = REAL(sexp_out);

  // Rounds of blocks, checking for interrupts between
  // them: each round has at least ~1e8 operations per
  // thread, or every remaining block.
  double blk_work  = ((double) n * (n + 3) / 2) * MORAN_BLOCK + 1;
  double rounds    = 1e8 * n_threads / blk_work;
  int    per_round = (rounds < n_blocks) ? (int) rounds : n_blocks; // no overflow
  if (per_round < n_threads) per_round = n_threads;

  for (mi_blk0 = 0; mi_blk0 < n_blocks; mi_blk0 = mi_blk1) {

    mi_blk1 = (n_blocks - mi_blk0 > per_round) ? mi_blk0 + per_round : n_blocks;

    int nb = mi_blk1 - mi_blk0;
    int nt = threads_for(nb < n_threads ? nb : n_threads, nb * blk_work, 1e6);

    // run_parallel() starts threads only when told of at
    // least 100 tasks.
    run_parallel(moran_worker, nt, nt > 1 ? 100 : 1);

    if (mi_blk1 < n_blocks) R_CheckUserInterrupt();
  }

  UNPROTECT(1);
  return sexp_out;
}
