// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The biplot scores follow pcoa_biplot() from scikit-bio
// (skbio/stats/ordination/_principal_coordinate_analysis.py),
// which carries this notice:
//
//   Copyright (c) 2013--, scikit-bio development team.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or
//   without modification, are permitted provided that the
//   following conditions are met:
//
//   * Redistributions of source code must retain the above
//     copyright notice, this list of conditions and the following
//     disclaimer.
//
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
//   * Neither the names scikit-bio, skbio, or biocore nor the names
//     of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written
//     permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
//   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
//   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
//   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
//   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//   POSSIBILITY OF SUCH DAMAGE.

/*
 * Biplot scores for a PCoA: the projection of each feature
 * (column of a table Y, n samples x p features) onto the
 * ordination's axes. As in scikit-bio's pcoa_biplot()
 * (Legendre & Legendre, Numerical Ecology, eqs. 9.44 and 9.55):
 *
 *     Z      = the coordinates, each axis centered and divided
 *              by its standard deviation (n - 1 denominator;
 *              a zero one taken as 1)
 *     S      = (1 / (n - 1)) * Y' * Z
 *     U_proj = sqrt(n - 1) * S * diag(eig^(-1/2))
 *
 * R passes P = Z * diag(eig^(-1/2)) / sqrt(n - 1), so here
 * U_proj = Yc' * P, with Yc the table with each feature
 * centered (scikit-bio multiplies Y itself; Z's columns sum
 * to 0, so the two agree in exact arithmetic, but centering
 * keeps rounding error proportional to each feature's spread
 * rather than to its mean).
 *
 * The table is read and transformed with ecodive's own
 * ecomatrix and normalize(), then copied by feature (CSC),
 * keeping only the ordination's samples. It stays sparse: a
 * feature with zeros costs O(nnz * k), its zeros summed in
 * closed form,
 *
 *     sum_{zeros} (0 - mean) P_a = -mean * (sum_all P_a - sum_{non-zeros} P_a).
 *
 * A feature without zeros, and every feature under CLR (where
 * a zero stands for its sample's CLR value), is centered
 * densely, O(n * k). A feature with the same value in every
 * sample scores exactly 0.
 *
 * Determinism: each feature's scores are summed by one
 * thread, over its samples in a fixed order, so results are
 * identical for any number of threads, and for dense and
 * sparse versions of the same table.
 */

#include "ecodive.h"


static int           bp_n;     // Samples (ordination rows).
static int           bp_k;     // Axes.
static int           bp_m;     // Features.
static const int    *bp_cp;    // CSC column starts (m + 1).
static const int    *bp_ci;    // ... each value's ordination row.
static const double *bp_cv;    // ... values (after `norm`).
static const double *bp_z0;    // CLR: each row's value for a zero, or NULL.
static const double *bp_P;     // n x k projection, column-major.
static const double *bp_Ps;    // k column sums of P.
static double       *bp_buf;   // n doubles per thread.
static double       *bp_out;   // m x k scores, column-major.


//------------------------------------------------------
// Scores of features w->i, w->i + w->n, ...
//------------------------------------------------------
static void *biplot_worker(void *arg) {

  worker_t *w = (worker_t *) arg;
  int       n = bp_n, k = bp_k, m = bp_m;
  double   *d = bp_buf + (size_t) n * w->i;

  for (int t = w->i; t < m; t += w->n) {

    int           beg = bp_cp[t];
    int           len = bp_cp[t + 1] - beg;
    const int    *ci  = bp_ci + beg;
    const double *cv  = bp_cv + beg;

    if (bp_z0 || len == n) {

      // Dense: zeros stand for z0 (CLR) or 0.
      for (int a = 0; a < n; a++) d[a] = bp_z0 ? bp_z0[a] : 0;
      for (int e = 0; e < len; e++) d[ci[e]] = cv[e];

      double mean = 0, lo = d[0], hi = d[0];
      for (int a = 0; a < n; a++) {
        mean += d[a];
        if (d[a] < lo) lo = d[a];
        if (d[a] > hi) hi = d[a];
      }
      mean /= n;

      for (int j = 0; j < k; j++) {
        double s = 0;
        if (lo != hi) {
          const double *p = bp_P + (size_t) n * j;
          for (int a = 0; a < n; a++) s += (d[a] - mean) * p[a];
        }
        bp_out[(size_t) m * j + t] = s;
      }
    }

    else {

      // Sparse: the n - len zeros in closed form.
      double mean = 0, lo = 0, hi = 0;
      for (int e = 0; e < len; e++) {
        mean += cv[e];
        if (cv[e] < lo) lo = cv[e];
        if (cv[e] > hi) hi = cv[e];
      }
      mean /= n;

      for (int j = 0; j < k; j++) {
        double s = 0;
        if (lo != hi) {
          const double *p  = bp_P + (size_t) n * j;
          double        pn = 0;
          for (int e = 0; e < len; e++) {
            s  += (cv[e] - mean) * p[ci[e]];
            pn += p[ci[e]];
          }
          s -= mean * (bp_Ps[j] - pn);
        }
        bp_out[(size_t) m * j + t] = s;
      }
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// counts:  A table in any form ecodive's counts take
//          (dense, Matrix, slam).
// margin, norm, pseudocount: as for counts elsewhere.
// P:       n x k projection matrix (see above).
// map:     For each sample of counts, its row in the
//          ordination (0-based), or -1 to leave it out.
//          Every row must be some sample's.
// n_threads
//
// Returns the m x k scores.
//======================================================
SEXP C_pcoa_biplot(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_pseudocount,
    SEXP sexp_P, SEXP sexp_map, SEXP sexp_n_threads ) {

  int    n           = nrows(sexp_P);
  int    k           = ncols(sexp_P);
  int    norm        = asInteger(sexp_norm);
  double pseudocount = asReal(sexp_pseudocount);
  int    n_threads   = asInteger(sexp_n_threads);

  if (n < 1 || k < 1 || n_threads < 1)
    error("Invalid arguments to C_pcoa_biplot()."); // # nocov

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int        ns  = em->n_samples;
  int        m   = em->n_otus;
  const int *map = INTEGER(sexp_map);

  if (XLENGTH(sexp_map) != ns) {
    free_all();                                      // # nocov
    error("Invalid arguments to C_pcoa_biplot()."); // # nocov
  }

  // Before normalizing, which could hide them (as 'binary'
  // does); and again after, for what `norm` makes.
  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1 && norm) normalize(em, norm, n_threads, pseudocount);
    for (int j = 0; j < ns; j++) {
      if (map[j] < 0) continue;
      for (int e = em->pos_vec[j]; e < em->pos_vec[j + 1]; e++)
        if (!R_FINITE(em->val_vec[e])) {
          free_all();
          error(pass == 0
            ? "`counts` must not contain NA, NaN, or infinite values."
            : "`counts` must not contain NA, NaN, or infinite values (after `norm`).");
        }
    }
  }


  //--------------------------------------------------
  // CSC copy, by feature, of the ordination's samples.
  //--------------------------------------------------
  int *cp = (int *) R_alloc((size_t) m + 1, sizeof(int));
  memset(cp, 0, ((size_t) m + 1) * sizeof(int));

  for (int j = 0; j < ns; j++) {
    if (map[j] < 0) continue;
    for (int e = em->pos_vec[j]; e < em->pos_vec[j + 1]; e++)
      cp[em->otu_vec[e] + 1]++;
  }
  for (int t = 0; t < m; t++) cp[t + 1] += cp[t];

  size_t  nnz  = (size_t) cp[m];
  int    *ci   = (int *)    R_alloc(nnz + 1, sizeof(int));
  double *cv   = (double *) R_alloc(nnz + 1, sizeof(double));
  int    *fill = (int *)    R_alloc((size_t) m + 1, sizeof(int));
  memcpy(fill, cp, (size_t) m * sizeof(int));

  // Samples in increasing order of ordination row, so that
  // each feature's values are summed in that order whatever
  // the order of samples in `counts`.
  int *order = (int *) R_alloc((size_t) n + 1, sizeof(int));
  for (int a = 0; a < n; a++) order[a] = -1;
  for (int j = 0; j < ns; j++)
    if (map[j] >= 0) {
      if (map[j] >= n || order[map[j]] >= 0) {
        free_all();                                      // # nocov
        error("Invalid arguments to C_pcoa_biplot()."); // # nocov
      }
      order[map[j]] = j;
    }

  for (int a = 0; a < n; a++) {
    int j = order[a];
    if (j < 0) { free_all(); error("Invalid arguments to C_pcoa_biplot()."); } // # nocov
    for (int e = em->pos_vec[j]; e < em->pos_vec[j + 1]; e++) {
      int t = em->otu_vec[e];
      ci[fill[t]] = a;
      cv[fill[t]] = em->val_vec[e];
      fill[t]++;
    }
  }

  double *z0 = NULL;
  if (em->clr_vec) {
    z0 = (double *) R_alloc((size_t) n, sizeof(double));
    for (int a = 0; a < n; a++) z0[a] = em->clr_vec[order[a]];
  }

  free_all();


  //--------------------------------------------------
  // Scores, features in parallel.
  //--------------------------------------------------
  const double *P  = REAL(sexp_P);
  double       *Ps = (double *) R_alloc((size_t) k, sizeof(double));
  for (int j = 0; j < k; j++) {
    Ps[j] = 0;
    for (int a = 0; a < n; a++) Ps[j] += P[(size_t) n * j + a];
  }

  SEXP sexp_out = PROTECT(allocMatrix(REALSXP, m, k));

  double work = ((double) nnz + (double) m) * k + (z0 ? (double) m * n * k : 0);
  int    nt   = threads_for(n_threads, work, 1e6);
  if (nt > m) nt = (m > 0) ? m : 1;

  bp_n   = n;
  bp_k   = k;
  bp_m   = m;
  bp_cp  = cp;
  bp_ci  = ci;
  bp_cv  = cv;
  bp_z0  = z0;
  bp_P   = P;
  bp_Ps  = Ps;
  bp_buf = (double *) R_alloc((size_t) n * nt, sizeof(double));
  bp_out = REAL(sexp_out);

  run_parallel(biplot_worker, nt, m);

  UNPROTECT(1);
  return sexp_out;
}
