// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// cross_assoc(): the correlation of every feature of x with every column of
// y, measured on the same samples. Kendall's tie-corrected variance (its z
// statistic) follows scipy.stats.kendalltau(method = 'asymptotic') from
// SciPy 1.18 (scipy/stats/_stats_py.py), under this license:
//
//   Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
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
//   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The rest of this file is Copyright 2026 ecodive authors, under the MIT
// license above.
//
// The columns of y are grouped by which samples they have (their pattern of
// missing values); R prepares each pattern's y (see cross_assoc_patterns()):
//
//   Pearson, Spearman: each column's values (ranks), centred, residualised
//     on the pattern's covariates when there are any, and scaled to unit
//     length: Y (m x q_P, stored row-major, one sample's values together).
//   Kendall: each column's values as integer codes 0 .. K - 1, and their
//     tie sums.
//
// Each task is one feature of x, in every pattern. The feature is read into
// a dense column (from a sparse or dense table), taken on the pattern's
// samples, ranked (Spearman), centred and residualised (Pearson, Spearman),
// giving xr and its norm. Since every y column sums to 0 (and is orthogonal
// to the covariates), for any constant c
//
//   r_j = sum_s (xr_s - c) / norm * Y_sj
//
// so with c the value of the feature's zeros, only its non-zero samples
// enter the sums: O(nnz x q) per feature. Kendall's tau-b is Knight's
// algorithm on y's codes in x's order (as mantel.c's).
//
// Permutations relabel y's rows (y_perm[s] = y[p[s]], R's y[p, ]) with one
// permutation of a pattern's samples for all its pairs. With covariates
// (ter Braak), the residualised y is permuted, and then is no longer
// orthogonal to the covariates, so each sum loses (Q'd) . (Q' P y_j), with
// Q the covariates' orthonormal basis and d = (xr - c) / norm.
//
// Every sum is over the samples in their order, and each feature's results
// are computed by one task, so they are identical for any number of
// threads; permutation hits are integers, kept per pair.

#include "ecodive.h"

#include <float.h> // DBL_EPSILON

#define CA_PEARSON  0
#define CA_SPEARMAN 1
#define CA_KENDALL  2

// Permutations per round: their indices (and, with covariates, the
// products Q' P y) are made first, then every feature runs through them.
#define CA_ROUND 256

typedef struct {
  int           m;        // Samples.
  const int    *rows;     // Their 0-based samples in x.
  int           qp;       // Columns of y.
  const int    *cols;     // Their 0-based columns of y.
  const double *Y;        // m x qp, row-major (Pearson, Spearman).
  const int    *codes;    // m x qp, row-major (Kendall).
  const double *ystat;    // 3 x qp: ytie, y0, y1 (Kendall).
  int           kq;       // Columns of Q (0: no covariates).
  const double *Q;        // m x kq, column-major.
  perm_source_t ps;
  int          *perm;     // This round: m x CA_ROUND permutations.
  double       *W;        // This round: kq x qp x CA_ROUND, Q' P Y.
} ca_pat_t;

static int           ca_method;
static int           ca_n, ca_p, ca_q;
static int           ca_npat;
static ca_pat_t     *ca_pat;
static const double *ca_xd;          // Dense x: n x p, or NULL.
static const int    *ca_xp;          // Sparse x: column starts (p + 1),
static const int    *ca_xi;          //   0-based samples,
static const double *ca_xv;          //   values.
static double       *ca_r;           // p x q (feature-major) correlations.
static double       *ca_z;           // p x q Kendall z, or NULL.
static int          *ca_hits;        // p x q hit counts, or NULL.
static int           ca_k0, ca_k1;   // This round's permutations, k0 .. k1 - 1.
static size_t        ca_wlen;        // Per-share workspace, in doubles.
static double       *ca_work;
static int           ca_mmax, ca_qmax, ca_kmax;


//------------------------------------------------------
// Stable merge sort of idx[0 .. m) by v[idx].
//------------------------------------------------------
static void ca_sort_index (const double *v, int *idx, int *tmp, int m) {

  for (int w = 1; w < m; w *= 2) {
    for (int lo = 0; lo < m; lo += 2 * w) {
      int mid = (lo + w < m)     ? lo + w     : m;
      int hi  = (lo + 2 * w < m) ? lo + 2 * w : m;
      int l = lo, r = mid, o = lo;
      while (l < mid && r < hi)
        tmp[o++] = (v[idx[r]] < v[idx[l]]) ? idx[r++] : idx[l++];
      while (l < mid) tmp[o++] = idx[l++];
      while (r < hi)  tmp[o++] = idx[r++];
    }
    memcpy(idx, tmp, (size_t) m * sizeof(int));
  }
}


//------------------------------------------------------
// One share's scratch space.
//------------------------------------------------------
typedef struct {
  double *xcol;   // n: the feature, dense.
  double *xs;     // mmax: on the pattern's samples.
  double *xr;     // mmax: ranked / residualised.
  double *d;      // mmax: the sums' weights ...
  double *acc;    // qmax
  double *qd;     // kmax: Q' d.
  int    *dpos;   // mmax: ... and their samples.
  int    *idx;    // mmax
  int    *tmp;    // mmax
  int    *seq;    // mmax
  int    *gbeg;   // mmax + 1
} ca_ws_t;

static ca_ws_t ca_ws (int share) {
  double *w = ca_work + (size_t) share * ca_wlen;
  size_t  m = (size_t) ca_mmax;
  ca_ws_t s;
  s.xcol = w;  w += ca_n;
  s.xs   = w;  w += m;
  s.xr   = w;  w += m;
  s.d    = w;  w += m;
  s.acc  = w;  w += ca_qmax;
  s.qd   = w;  w += ca_kmax;
  int *iw = (int *) w;
  s.dpos = iw; iw += m;
  s.idx  = iw; iw += m;
  s.tmp  = iw; iw += m;
  s.seq  = iw; iw += m;
  s.gbeg = iw;
  return s;
}

static size_t ca_ws_len (void) {
  size_t m = (size_t) ca_mmax;
  size_t ints = 5 * m + 1;
  return (size_t) ca_n + 4 * m + ca_qmax + ca_kmax + (ints + 1) / 2 + 1;
}


//------------------------------------------------------
// Feature f, dense, into ws->xcol.
//------------------------------------------------------
static void ca_read_feature (int f, ca_ws_t *ws) {
  if (ca_xd) {
    memcpy(ws->xcol, ca_xd + (size_t) f * ca_n, (size_t) ca_n * sizeof(double));
    return;
  }
  memset(ws->xcol, 0, (size_t) ca_n * sizeof(double));
  for (int t = ca_xp[f]; t < ca_xp[f + 1]; t++) ws->xcol[ca_xi[t]] = ca_xv[t];
}


//------------------------------------------------------
// Pearson, Spearman: the weights d and their samples, for
// the feature in ws->xcol and pattern P. Returns their
// count, or -1 when the feature is constant (no
// correlation).
//------------------------------------------------------
static int ca_prep (const ca_pat_t *P, ca_ws_t *ws, int want_qd) {

  int     m  = P->m;
  double *xs = ws->xs, *xr = ws->xr;

  int has_zero = 0, constant = 1;
  for (int s = 0; s < m; s++) {
    xs[s] = ws->xcol[P->rows[s]];
    if (xs[s] == 0) has_zero = 1;
    if (xs[s] != xs[0]) constant = 0;
  }
  if (constant) return -1;

  if (ca_method == CA_SPEARMAN) {
    // Average ranks of ties, as rank() and scipy.stats.rankdata().
    int *idx = ws->idx;
    for (int s = 0; s < m; s++) idx[s] = s;
    ca_sort_index(xs, idx, ws->tmp, m);
    for (int lo = 0; lo < m; ) {
      int hi = lo + 1;
      while (hi < m && xs[idx[hi]] == xs[idx[lo]]) hi++;
      double rk = (double) (lo + 1 + hi) / 2;
      for (int t = lo; t < hi; t++) xr[idx[t]] = rk;
      lo = hi;
    }
  }
  else {
    memcpy(xr, xs, (size_t) m * sizeof(double));
  }

  // The value of the zeros (all the same, ranked or not), or the mean.
  double mean = 0;
  for (int s = 0; s < m; s++) mean += xr[s];
  mean /= m;
  double c = mean;
  if (has_zero)
    for (int s = 0; s < m; s++) if (xs[s] == 0) { c = xr[s]; break; }

  // The norm of the centred, residualised feature (in d, for now).
  double *res = ws->d;
  double  ss0 = 0;
  for (int s = 0; s < m; s++) { res[s] = xr[s] - mean; ss0 += res[s] * res[s]; }

  double ss = ss0;
  if (P->kq > 0) {
    for (int l = 0; l < P->kq; l++) {
      const double *q = P->Q + (size_t) l * m;
      double        a = 0;
      for (int s = 0; s < m; s++) a += q[s] * res[s];
      for (int s = 0; s < m; s++) res[s] -= a * q[s];
    }
    ss = 0;
    for (int s = 0; s < m; s++) ss += res[s] * res[s];
    // Explained by the covariates (within rounding): no partial correlation.
    if (!(ss > 1e-20 * ss0)) return -1;
  }
  if (!(ss > 0)) return -1;

  double inv = 1 / sqrt(ss);

  if (want_qd)
    for (int l = 0; l < P->kq; l++) {
      const double *q = P->Q + (size_t) l * m;
      double        a = 0;
      for (int s = 0; s < m; s++) a += q[s] * (xr[s] - c);
      ws->qd[l] = a * inv;
    }

  int L = 0;
  for (int s = 0; s < m; s++) {
    if (xr[s] == c) continue;
    ws->dpos[L] = s;
    ws->d[L]    = (xr[s] - c) * inv;
    L++;
  }
  return L;
}


//------------------------------------------------------
// The sums r_j = sum_t d_t Y[p[s_t], j] into ws->acc.
// p NULL: the identity.
//------------------------------------------------------
static void ca_sums (const ca_pat_t *P, ca_ws_t *ws, int L, const int *p) {

  int     qp  = P->qp;
  double *acc = ws->acc;
  for (int j = 0; j < qp; j++) acc[j] = 0;

  for (int t = 0; t < L; t++) {
    int           s   = p ? p[ws->dpos[t]] : ws->dpos[t];
    const double *row = P->Y + (size_t) s * qp;
    double        dt  = ws->d[t];
    for (int j = 0; j < qp; j++) acc[j] += dt * row[j];
  }
}


static inline double ca_clamp (double r) {
  return (r > 1) ? 1 : (r < -1) ? -1 : r;
}


//------------------------------------------------------
// Kendall: x's order on pattern P's samples (ws->idx), its
// tie groups (ws->gbeg) and their sums. Returns the
// number of groups.
//------------------------------------------------------
static int ca_kendall_prep (const ca_pat_t *P, ca_ws_t *ws, double *xs3) {

  int     m  = P->m;
  double *xs = ws->xs;
  int    *idx = ws->idx;

  for (int s = 0; s < m; s++) { xs[s] = ws->xcol[P->rows[s]]; idx[s] = s; }
  ca_sort_index(xs, idx, ws->tmp, m);

  int    ng = 0;
  double t0 = 0, t1 = 0, t2 = 0;
  for (int lo = 0; lo < m; ) {
    int hi = lo + 1;
    while (hi < m && xs[idx[hi]] == xs[idx[lo]]) hi++;
    double c = hi - lo;
    t0 += c * (c - 1) / 2;
    t1 += c * (c - 1) * (c - 2);
    t2 += c * (c - 1) * (2 * c + 5);
    ws->gbeg[ng++] = lo;
    lo = hi;
  }
  ws->gbeg[ng] = m;
  xs3[0] = t0; xs3[1] = t1; xs3[2] = t2;
  return ng;
}


//------------------------------------------------------
// Kendall's tau-b of column j (in x's order, read through
// permutation p) and its z, as scipy.stats.kendalltau().
//------------------------------------------------------
static void ca_kendall (
    const ca_pat_t *P, ca_ws_t *ws, int ng, const double *xs3,
    int j, const int *p, double *tau, double *z ) {

  int     m   = P->m;
  int     qp  = P->qp;
  int    *seq = ws->seq;
  double  tot = (double) m * (m - 1) / 2;
  double  xtie = xs3[0], ytie = P->ystat[3 * j];

  if (xtie == tot || ytie == tot) { *tau = NAN; *z = NAN; return; }

  for (int t = 0; t < m; t++) {
    int s = ws->idx[t];
    if (p) s = p[s];
    seq[t] = P->codes[(size_t) s * qp + j];
  }

  double ntie = 0;
  for (int g = 0; g < ng; g++) {
    size_t lo = (size_t) ws->gbeg[g], len = (size_t) ws->gbeg[g + 1] - lo;
    if (len < 2) continue;
    kendall_sort_count(seq + lo, ws->tmp, len);
    ntie += kendall_tied_sorted(seq + lo, len);
  }
  double dis = (double) kendall_sort_count(seq, ws->tmp, (size_t) m);

  double cmd = tot - xtie - ytie + ntie - 2 * dis;
  *tau = ca_clamp(cmd / sqrt(tot - xtie) / sqrt(tot - ytie));

  double mm  = (double) m * (m - 1);
  double y0  = P->ystat[3 * j + 1], y1 = P->ystat[3 * j + 2];
  double var = (mm * (2.0 * m + 5) - xs3[2] - y1) / 18 +
               (2 * xtie * ytie) / mm + xs3[1] * y0 / (9 * mm * (m - 2));
  *z = cmd / sqrt(var);
}


//------------------------------------------------------
// Observed correlations: features i, i + n, ...
//------------------------------------------------------
static void *ca_worker (void *arg) {

  worker_t *w  = (worker_t *) arg;
  ca_ws_t   ws = ca_ws(w->i);

  for (int f = w->i; f < ca_p; f += w->n) {

    ca_read_feature(f, &ws);
    double *rf = ca_r + (size_t) f * ca_q;

    for (int g = 0; g < ca_npat; g++) {
      const ca_pat_t *P = ca_pat + g;

      if (ca_method == CA_KENDALL) {
        double xs3[3];
        int    ng = ca_kendall_prep(P, &ws, xs3);
        for (int j = 0; j < P->qp; j++) {
          double tau, z;
          ca_kendall(P, &ws, ng, xs3, j, NULL, &tau, &z);
          rf[P->cols[j]] = tau;
          ca_z[(size_t) f * ca_q + P->cols[j]] = z;
        }
        continue;
      }

      int L = ca_prep(P, &ws, 0);
      if (L < 0) {
        for (int j = 0; j < P->qp; j++) rf[P->cols[j]] = NAN;
        continue;
      }
      ca_sums(P, &ws, L, NULL);
      for (int j = 0; j < P->qp; j++) {
        // A constant y column is all zeros: NA, set in R.
        rf[P->cols[j]] = ca_clamp(ws.acc[j]);
      }
    }
  }
  return NULL;
}


//------------------------------------------------------
// This round's permutations (and Q' P Y): permutations
// k0 + i, k0 + i + n, ...
//------------------------------------------------------
static void *ca_perm_worker (void *arg) {

  worker_t *w = (worker_t *) arg;

  for (int k = ca_k0 + w->i; k < ca_k1; k += w->n) {
    int b = k - ca_k0;
    for (int g = 0; g < ca_npat; g++) {
      ca_pat_t *P = ca_pat + g;
      int      *p = P->perm + (size_t) b * P->m;
      perm_source_get(&P->ps, k, p);
      if (P->kq == 0) continue;
      double *W = P->W + (size_t) b * P->kq * P->qp;
      for (int j = 0; j < P->qp; j++)
        for (int l = 0; l < P->kq; l++) {
          const double *q = P->Q + (size_t) l * P->m;
          double        a = 0;
          for (int s = 0; s < P->m; s++) a += q[s] * P->Y[(size_t) p[s] * P->qp + j];
          W[(size_t) j * P->kq + l] = a;
        }
    }
  }
  return NULL;
}


//------------------------------------------------------
// This round's permuted correlations: features i, i + n,
// ...; hits where |r_b| >= |r| - sqrt(eps), as perm_p().
//------------------------------------------------------
static void *ca_hit_worker (void *arg) {

  worker_t *w   = (worker_t *) arg;
  ca_ws_t   ws  = ca_ws(w->i);
  double    tol = sqrt(DBL_EPSILON);

  for (int f = w->i; f < ca_p; f += w->n) {

    ca_read_feature(f, &ws);
    const double *rf = ca_r    + (size_t) f * ca_q;
    int          *hf = ca_hits + (size_t) f * ca_q;

    for (int g = 0; g < ca_npat; g++) {
      const ca_pat_t *P = ca_pat + g;

      if (ca_method == CA_KENDALL) {
        double xs3[3];
        int    ng = ca_kendall_prep(P, &ws, xs3);
        for (int b = 0; b < ca_k1 - ca_k0; b++) {
          const int *p = P->perm + (size_t) b * P->m;
          for (int j = 0; j < P->qp; j++) {
            double r = rf[P->cols[j]];
            if (isnan(r)) continue;
            double tau, z;
            ca_kendall(P, &ws, ng, xs3, j, p, &tau, &z);
            if (fabs(tau) >= fabs(r) - tol) hf[P->cols[j]]++;
          }
        }
        continue;
      }

      int L = ca_prep(P, &ws, P->kq > 0);
      if (L < 0) continue;

      for (int b = 0; b < ca_k1 - ca_k0; b++) {
        const int *p = P->perm + (size_t) b * P->m;
        ca_sums(P, &ws, L, p);
        const double *W = P->kq ? P->W + (size_t) b * P->kq * P->qp : NULL;
        for (int j = 0; j < P->qp; j++) {
          double r = rf[P->cols[j]];
          if (isnan(r)) continue;
          double rb = ws.acc[j];
          if (W)
            for (int l = 0; l < P->kq; l++) rb -= ws.qd[l] * W[(size_t) j * P->kq + l];
          if (fabs(ca_clamp(rb)) >= fabs(r) - tol) hf[P->cols[j]]++;
        }
      }
    }
  }
  return NULL;
}


//======================================================
// R interface.
//
// sexp_x:     list(n, p, dense (n x p double, or NULL),
//             starts, samples, values (sparse; when dense
//             is NULL)).
// sexp_pats:  list of patterns: list(rows, cols, Y, codes,
//             ystat, Q, perms, strata, seed), see ca_pat_t;
//             rows and cols 0-based; perms an m x n_perm
//             integer matrix of 0-based permutations or
//             NULL; strata 0-based codes or NULL.
// sexp_q:     columns of y.
// sexp_nperm: permutations.
//
// Returns list(r, z, hits): p x q (feature-major) vectors;
// z only for Kendall, hits only with permutations.
//======================================================
SEXP C_cross_assoc(
    SEXP sexp_x, SEXP sexp_pats, SEXP sexp_method, SEXP sexp_q,
    SEXP sexp_nperm, SEXP sexp_cpus ) {

  ca_method = asInteger(sexp_method);
  ca_n      = asInteger(VECTOR_ELT(sexp_x, 0));
  ca_p      = asInteger(VECTOR_ELT(sexp_x, 1));
  ca_q      = asInteger(sexp_q);
  int n_perm    = asInteger(sexp_nperm);
  int n_threads = asInteger(sexp_cpus);

  SEXP dense = VECTOR_ELT(sexp_x, 2);
  ca_xd = isNull(dense) ? NULL : REAL(dense);
  if (!ca_xd) {
    ca_xp = INTEGER(VECTOR_ELT(sexp_x, 3));
    ca_xi = INTEGER(VECTOR_ELT(sexp_x, 4));
    ca_xv = REAL(VECTOR_ELT(sexp_x, 5));
  }

  ca_npat = LENGTH(sexp_pats);
  ca_pat  = (ca_pat_t *) R_alloc(ca_npat > 0 ? ca_npat : 1, sizeof(ca_pat_t));
  ca_mmax = 1; ca_qmax = 1; ca_kmax = 1;

  double work = 0;
  for (int g = 0; g < ca_npat; g++) {
    SEXP      sp = VECTOR_ELT(sexp_pats, g);
    ca_pat_t *P  = ca_pat + g;
    SEXP rows = VECTOR_ELT(sp, 0), cols = VECTOR_ELT(sp, 1);
    SEXP Y    = VECTOR_ELT(sp, 2), codes = VECTOR_ELT(sp, 3);
    SEXP ys   = VECTOR_ELT(sp, 4), Q = VECTOR_ELT(sp, 5);
    P->m     = LENGTH(rows);
    P->rows  = INTEGER(rows);
    P->qp    = LENGTH(cols);
    P->cols  = INTEGER(cols);
    P->Y     = isNull(Y)     ? NULL : REAL(Y);
    P->codes = isNull(codes) ? NULL : INTEGER(codes);
    P->ystat = isNull(ys)    ? NULL : REAL(ys);
    P->kq    = isNull(Q)     ? 0    : ncols(Q);
    P->Q     = isNull(Q)     ? NULL : REAL(Q);
    P->perm  = NULL;
    P->W     = NULL;
    if (P->m  > ca_mmax) ca_mmax = P->m;
    if (P->qp > ca_qmax) ca_qmax = P->qp;
    if (P->kq > ca_kmax) ca_kmax = P->kq;
    double lg = (ca_method == CA_KENDALL) ? log2((double) P->m + 2) : 1;
    work += (double) P->m * P->qp * lg;

    if (n_perm > 0) {
      perm_source_init(&P->ps, P->m, VECTOR_ELT(sp, 6), VECTOR_ELT(sp, 7),
                       asInteger(VECTOR_ELT(sp, 8)));
      P->perm = (int *) R_alloc((size_t) P->m * CA_ROUND, sizeof(int));
      if (P->kq > 0)
        P->W = (double *) R_alloc((size_t) P->kq * P->qp * CA_ROUND, sizeof(double));
    }
  }

  size_t pq = (size_t) ca_p * ca_q;
  SEXP sexp_r = PROTECT(allocVector(REALSXP, pq));
  SEXP sexp_z = PROTECT(ca_method == CA_KENDALL ? allocVector(REALSXP, pq) : R_NilValue);
  SEXP sexp_h = PROTECT(n_perm > 0 ? allocVector(INTSXP, pq) : R_NilValue);
  ca_r    = REAL(sexp_r);
  ca_z    = isNull(sexp_z) ? NULL : REAL(sexp_z);
  ca_hits = isNull(sexp_h) ? NULL : INTEGER(sexp_h);
  for (size_t t = 0; t < pq; t++) ca_r[t] = NA_REAL;
  if (ca_z)    for (size_t t = 0; t < pq; t++) ca_z[t] = NA_REAL;
  if (ca_hits) memset(ca_hits, 0, pq * sizeof(int));

  // Each share's scratch, by share index (at most n_threads shares).
  if (n_threads < 1) n_threads = 1;
  ca_wlen = ca_ws_len();
  ca_work = (double *) R_alloc(ca_wlen * n_threads, sizeof(double));

  // Observed.
  int nt = threads_for(n_threads, work * ca_p, 1e6);
  if (nt > ca_p) nt = ca_p > 0 ? ca_p : 1;
  run_parallel(ca_worker, nt, (nt > 1) ? 100 : 1);

  // Permutations, in rounds.
  for (int k0 = 1; k0 <= n_perm; k0 += CA_ROUND) {

    if (user_interrupted()) { UNPROTECT(3); error("Interrupted."); } // # nocov

    ca_k0 = k0;
    ca_k1 = (k0 + CA_ROUND <= n_perm + 1) ? k0 + CA_ROUND : n_perm + 1;

    int ntp = threads_for(n_threads, (double) (ca_k1 - ca_k0) * ca_n * (1 + ca_kmax * ca_qmax), 1e6);
    if (ntp > ca_k1 - ca_k0) ntp = ca_k1 - ca_k0;
    run_parallel(ca_perm_worker, ntp, (ntp > 1) ? 100 : 1);

    int nth = threads_for(n_threads, work * ca_p * (ca_k1 - ca_k0), 1e6);
    if (nth > ca_p) nth = ca_p > 0 ? ca_p : 1;
    run_parallel(ca_hit_worker, nth, (nth > 1) ? 100 : 1);
  }

  SEXP res = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, sexp_r);
  SET_VECTOR_ELT(res, 1, sexp_z);
  SET_VECTOR_ELT(res, 2, sexp_h);
  UNPROTECT(4);
  return res;
}
