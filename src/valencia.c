// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// VALENCIA's nearest-centroid classification (valencia()): a C port of the
// similarity, assignment and score steps of Valencia.py (VALENCIA, commit
// 8559d45 of github.com/ravel-lab/VALENCIA), restructured for ecodive: the
// Yue-Clayton similarity is computed sparsely, over each sample's non-zero
// taxa, and samples are split over threads (see plans/community_typing.md).
// Valencia.py is under this licence:
//
// MIT License
//
// Copyright (c) 2023 ravel-lab
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "ecodive.h"


// Yue & Clayton's (2005) similarity of a sample's proportions p to a
// centroid q, over the union of their taxa (missing values are 0):
//
//   theta = sum(p q) / (sum((p - q)^2) + sum(p q))
//         = sum(p q) / (sum(p^2) + sum(q^2) - sum(p q))
//
// The second form needs only the sample's non-zero taxa, given each
// centroid's sum(q^2). Its denominator is at least (sum(p^2) + sum(q^2)) / 2,
// since sum(p q) <= (sum(p^2) + sum(q^2)) / 2, so it suffers no cancellation.

static int     vl_n_samples;
static int     vl_n_cents;
static int    *vl_pos_vec;
static int    *vl_otu_vec;
static double *vl_val_vec;
static int    *vl_col;     // each feature's centroid column, or -1
static double *vl_cents;   // n_cents x n_taxa, column-major
static double *vl_q2;      // each centroid's sum(q^2)
static double *vl_totals;  // each sample's total, or NULL for its sum
static double *vl_buf;     // n_threads x n_cents
static double *vl_scores;  // n_samples x n_cents, column-major
static int    *vl_best;    // 1-based, or NA
static double *vl_matched; // proportion of the sample in centroid taxa


static void *vl_worker (void *arg) {

  int     thread_i  = ((worker_t *) arg)->i;
  int     n_threads = ((worker_t *) arg)->n;
  int     n         = vl_n_samples;
  int     K         = vl_n_cents;
  double *pq        = vl_buf + (size_t) thread_i * K;

  for (int s = thread_i; s < n; s += n_threads) {

    int beg = vl_pos_vec[s], end = vl_pos_vec[s + 1];

    double total = 0;
    if (vl_totals) total = vl_totals[s];
    else           for (int i = beg; i < end; i++) total += vl_val_vec[i];

    if (!(total > 0)) {
      for (int k = 0; k < K; k++) vl_scores[s + (size_t) n * k] = R_NaN;
      vl_best[s]    = NA_INTEGER;
      vl_matched[s] = NA_REAL;
      continue;
    }

    for (int k = 0; k < K; k++) pq[k] = 0;
    double pp = 0, matched = 0;

    for (int i = beg; i < end; i++) {
      double p = vl_val_vec[i] / total;
      pp += p * p;
      int col = vl_col[vl_otu_vec[i]];
      if (col < 0) continue;
      matched += p;
      const double *q = vl_cents + (size_t) col * K;
      for (int k = 0; k < K; k++) pq[k] += p * q[k];
    }

    // The largest similarity, the first on ties, skipping NaN (a centroid
    // of zeros against a sample of zeros), as pandas' idxmax() does.
    int best = -1;
    for (int k = 0; k < K; k++) {
      double theta = pq[k] / (pp + vl_q2[k] - pq[k]);
      vl_scores[s + (size_t) n * k] = theta;
      if (!ISNAN(theta) && (best < 0 || theta > vl_scores[s + (size_t) n * best]))
        best = k;
    }

    vl_best[s]    = best < 0 ? NA_INTEGER : best + 1;
    vl_matched[s] = matched;
  }

  return NULL;
}



//======================================================
// R interface.
//
// sexp_col:    for each feature, its 1-based centroid column, or 0.
// sexp_cents:  the centroids (n_cents x n_taxa, doubles).
// sexp_totals: each sample's total, or NULL to sum its counts.
//
// Returns list(scores (n_samples x n_cents), best (1-based or NA),
// matched).
//======================================================

SEXP C_valencia (SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_col,
                 SEXP sexp_cents, SEXP sexp_totals, SEXP sexp_n_threads) {

  init_n_ptrs(20);

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  int n      = em->n_samples;
  int n_otus = em->n_otus;
  int K      = nrows(sexp_cents);
  int T      = ncols(sexp_cents);

  vl_n_samples = n;
  vl_n_cents   = K;
  vl_pos_vec   = em->pos_vec;
  vl_otu_vec   = em->otu_vec;
  vl_val_vec   = em->val_vec;
  vl_cents     = REAL(sexp_cents);
  vl_totals    = isNull(sexp_totals) ? NULL : REAL(sexp_totals);

  int *col_in = INTEGER(sexp_col);
  vl_col = (int *) safe_malloc((size_t) (n_otus + 1) * sizeof(int));
  for (int j = 0; j < n_otus; j++) vl_col[j] = col_in[j] - 1;

  vl_q2 = (double *) safe_malloc((size_t) (K + 1) * sizeof(double));
  for (int k = 0; k < K; k++) {
    double q2 = 0;
    for (int t = 0; t < T; t++) {
      double q = vl_cents[k + (size_t) K * t];
      q2 += q * q;
    }
    vl_q2[k] = q2;
  }

  vl_buf = (double *) safe_malloc((size_t) n_threads * (K + 1) * sizeof(double));

  SEXP sexp_scores  = PROTECT(allocMatrix(REALSXP, n, K));
  SEXP sexp_best    = PROTECT(allocVector(INTSXP, n));
  SEXP sexp_matched = PROTECT(allocVector(REALSXP, n));
  vl_scores  = REAL(sexp_scores);
  vl_best    = INTEGER(sexp_best);
  vl_matched = REAL(sexp_matched);

  run_parallel(vl_worker, n_threads, n);

  free_all();

  SEXP sexp_result = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(sexp_result, 0, sexp_scores);
  SET_VECTOR_ELT(sexp_result, 1, sexp_best);
  SET_VECTOR_ELT(sexp_result, 2, sexp_matched);

  UNPROTECT(4);
  return sexp_result;
}
