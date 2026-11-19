// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The search follows bioenv() from scikit-bio
// (skbio/stats/distance/_bioenv.py), which carries this notice:
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

#include "ecodive.h"

#include <limits.h> // INT_MAX


//======================================================
// BIOENV (Clarke & Ainsworth 1993): for each subset size
// k = 1 .. max_vars, the subset of the n x p standardized
// variables z whose distances among the samples have the
// highest correlation with the community distances (the
// first target), as in scikit-bio's bioenv().
//
// Subsets of each size are visited in lexicographic order
// (as itertools.combinations() and utils::combn()), and
// the first with the highest correlation wins, as in
// scikit-bio and vegan.
//
// Each subset is computed by one thread, start to finish,
// and threads take contiguous runs of subsets whose bests
// are then combined in order, so results are identical
// for any thread count.
//
// The targets are fixed: standardized values (Pearson),
// standardized ranks (Spearman), or integer codes of the
// distinct values (Kendall), as mantel.c takes them. With
// a second target (`given`), the correlation is the
// partial one, given rxz, the first target's correlation
// with the second.
//======================================================


#define BE_PEARSON  0
#define BE_SPEARMAN 1
#define BE_KENDALL  2

#define BE_EUCLIDEAN   0
#define BE_MANHATTAN   1
#define BE_GOWER       2
#define BE_MAHALANOBIS 3

// Merge sort: insertion-sorted runs of this length.
#define BE_RUN 32

static int           be_n;       // Samples.
static size_t        be_np;      // Pairs, n(n - 1) / 2.
static int           be_p;       // Variables.
static int           be_k;       // This round's subset size.
static int           be_m;       // Targets, 1 or 2.
static int           be_method;
static int           be_metric;
static const double *be_z;       // n x p variables, row-major.
static const double *be_t[2];    // Pearson / Spearman targets.
static const int    *be_tc[2];   // Kendall targets' codes.
static double        be_n2[2];   // Kendall: pairs of pairs tied in each target.
static double        be_rxz;     // Correlation of the targets.
static const double *be_binom;   // (p + 1) x (p + 1) binomial coefficients.

static double        be_r0;      // This round's subsets: ranks r0 .. r1 - 1.
static double        be_r1;
static double       *be_best_r;  // Per thread: best correlation ...
static double       *be_best_i;  // ... and its rank.

static double       *be_dwork;   // Per-thread workspaces.
static size_t        be_dwlen;
static int          *be_iwork;
static size_t        be_iwlen;


#define BINOM(a, b) be_binom[(size_t) (a) * (size_t) (be_p + 1) + (size_t) (b)]


//------------------------------------------------------
// The subset of rank r (0-based) among all subsets of
// size k of 0 .. p - 1, in lexicographic order.
//------------------------------------------------------
static void unrank(double r, int k, int *comb) {
  int c = 0;
  for (int i = 0; i < k; i++) {
    for (;;) {
      double cnt = BINOM(be_p - c - 1, k - i - 1);
      if (r < cnt) { comb[i] = c++; break; }
      r -= cnt;
      c++;
    }
  }
}


// The next subset in lexicographic order.
static void next_comb(int *comb, int k, int p) {
  int i = k - 1;
  while (i >= 0 && comb[i] == p - k + i) i--;
  if (i < 0) return;
  comb[i]++;
  for (int j = i + 1; j < k; j++) comb[j] = comb[j - 1] + 1;
}


//------------------------------------------------------
// Stable sort of the indices 0 .. len - 1 by d, with
// buf (len ints) as scratch.
//------------------------------------------------------
static void order_by(const double *d, int *a, int *buf, size_t len) {

  for (size_t i = 0; i < len; i++) a[i] = (int) i;

  for (size_t lo = 0; lo < len; lo += BE_RUN) {
    size_t hi = (lo + BE_RUN < len) ? lo + BE_RUN : len;
    for (size_t i = lo + 1; i < hi; i++) {
      int    v  = a[i];
      double dv = d[v];
      size_t k  = i;
      while (k > lo && d[a[k - 1]] > dv) { a[k] = a[k - 1]; k--; }
      a[k] = v;
    }
  }

  int *src = a, *dst = buf;
  for (size_t w = BE_RUN; w < len; w *= 2) {
    for (size_t lo = 0; lo < len; lo += 2 * w) {
      size_t mid = (lo + w < len)     ? lo + w     : len;
      size_t hi  = (lo + 2 * w < len) ? lo + 2 * w : len;
      size_t l = lo, r = mid, o = lo;
      while (l < mid && r < hi)
        dst[o++] = (d[src[r]] < d[src[l]]) ? src[r++] : src[l++];
      while (l < mid) dst[o++] = src[l++];
      while (r < hi)  dst[o++] = src[r++];
    }
    int *tmp = src; src = dst; dst = tmp;
  }
  if (src != a) memcpy(a, src, len * sizeof(int));
}


//------------------------------------------------------
// The distances among the samples on the variables in
// comb (k of them), in dist order: d (N doubles). w:
// n * k + k * k + k doubles of scratch, for Mahalanobis.
// Returns 0, or 1 if the variables' covariance is
// singular (Mahalanobis).
//------------------------------------------------------
static int subset_dist(const int *comb, int k, double *d, double *w) {

  int           n = be_n, p = be_p;
  const double *z = be_z;
  size_t        s = 0;

  if (be_metric == BE_MAHALANOBIS) {

    // Covariance of the (centred) variables and its
    // Cholesky factor L; then each pair's distance,
    // |L^-1 (z_i - z_j)|, from the difference of its
    // samples: pairs with the same differences have the
    // same distance, and for one variable the distances
    // are in the order of Euclidean's, ties and all.
    double *x = w;           // n x k, row-major.
    double *L = w + (size_t) n * k;
    double *u = L + (size_t) k * k;

    for (int i = 0; i < n; i++)
      for (int a = 0; a < k; a++)
        x[(size_t) i * k + a] = z[(size_t) i * p + comb[a]];

    for (int a = 0; a < k; a++)
      for (int b = 0; b <= a; b++) {
        double c = 0;
        for (int i = 0; i < n; i++) c += x[(size_t) i * k + a] * x[(size_t) i * k + b];
        L[a * k + b] = c / (n - 1);
      }

    for (int a = 0; a < k; a++) {
      double diag = L[a * k + a];
      for (int b = 0; b <= a; b++) {
        double v = L[a * k + b];
        for (int c = 0; c < b; c++) v -= L[a * k + c] * L[b * k + c];
        if (b < a) {
          L[a * k + b] = v / L[b * k + b];
        }
        else {
          if (!(v > 1e-10 * diag)) return 1;
          L[a * k + a] = sqrt(v);
        }
      }
    }

    for (int j = 0; j < n - 1; j++) {
      const double *xj = x + (size_t) j * k;
      for (int i = j + 1; i < n; i++) {
        const double *xi = x + (size_t) i * k;
        double acc = 0;
        for (int a = 0; a < k; a++) {
          double v = xi[a] - xj[a];
          for (int c = 0; c < a; c++) v -= L[a * k + c] * u[c];
          u[a] = v / L[a * k + a];
          acc += u[a] * u[a];
        }
        d[s++] = sqrt(acc);
      }
    }

    return 0;
  }

  for (int j = 0; j < n - 1; j++) {
    const double *zj = z + (size_t) j * p;
    for (int i = j + 1; i < n; i++) {
      const double *zi = z + (size_t) i * p;
      double acc = 0;
      if (be_metric == BE_EUCLIDEAN) {
        for (int a = 0; a < k; a++) { double dev = zi[comb[a]] - zj[comb[a]]; acc += dev * dev; }
        d[s++] = sqrt(acc);
      }
      else {
        for (int a = 0; a < k; a++) acc += fabs(zi[comb[a]] - zj[comb[a]]);
        d[s++] = (be_metric == BE_GOWER) ? acc / k : acc;
      }
    }
  }

  return 0;
}


//------------------------------------------------------
// Centre d and scale it to unit length. Returns 0 if it
// is constant (no correlation), else 1.
//------------------------------------------------------
static int standardize(double *d, size_t N) {

  double mean = 0;
  for (size_t s = 0; s < N; s++) mean += d[s];
  mean /= (double) N;

  double ss = 0, mx = 0;
  for (size_t s = 0; s < N; s++) {
    d[s] -= mean;
    ss += d[s] * d[s];
    if (fabs(d[s]) > mx) mx = fabs(d[s]);
  }

  if (!(mx > 0)) return 0;

  // Very small or large values: scale first, so that the
  // sum of squares neither underflows nor overflows.
  if (!(isfinite(ss) && ss >= 1e-290)) {
    ss = 0;
    for (size_t s = 0; s < N; s++) { d[s] /= mx; ss += d[s] * d[s]; }
  }

  double f = 1 / sqrt(ss);
  for (size_t s = 0; s < N; s++) d[s] *= f;
  return 1;
}


//------------------------------------------------------
// The correlation of the subset with the community
// distances: partial, given the second target, if any.
//------------------------------------------------------
static double subset_cor(const int *comb, int k, double *d, double *w, int *iw) {

  size_t N = be_np;
  double r[2];

  if (subset_dist(comb, k, d, w)) return NAN;

  if (be_method == BE_KENDALL) {

    int *ord = iw, *buf = iw + N, *seq = iw + 2 * N;

    order_by(d, ord, buf, N);

    // Pairs of pairs tied in d.
    double n1 = 0;
    size_t run = 1;
    for (size_t s = 1; s <= N; s++) {
      if (s < N && d[ord[s]] == d[ord[s - 1]]) { run++; continue; }
      n1 += (double) run * (double) (run - 1) / 2;
      run = 1;
    }

    double n0 = (double) N * (double) (N - 1) / 2;
    if (!(n0 - n1 > 0)) return NAN;

    for (int c = 0; c < be_m; c++) {

      const int *tc = be_tc[c];
      for (size_t s = 0; s < N; s++) seq[s] = tc[ord[s]];

      // Sort within d's tie groups, so that only pairs
      // between groups can be inversions.
      double n3 = 0;
      size_t lo = 0;
      for (size_t s = 1; s <= N; s++) {
        if (s < N && d[ord[s]] == d[ord[s - 1]]) continue;
        if (s - lo > 1) {
          kendall_sort_count(seq + lo, buf, s - lo);
          n3 += kendall_tied_sorted(seq + lo, s - lo);
        }
        lo = s;
      }

      double disc = (double) kendall_sort_count(seq, buf, N);
      double n2   = be_n2[c];
      double sg   = 2 * (n0 - n1 - n2 + n3 - 2 * disc);
      double rc   = sg / (sqrt(2 * (n0 - n1)) * sqrt(2 * (n0 - n2)));
      r[c] = (rc > 1) ? 1 : rc;
    }
  }
  else {

    if (be_method == BE_SPEARMAN) {

      // Average ranks of tied values, written over d.
      int *ord = iw, *buf = iw + N;
      order_by(d, ord, buf, N);
      size_t lo = 0;
      for (size_t s = 1; s <= N; s++) {
        if (s < N && d[ord[s]] == d[ord[s - 1]]) continue;
        double rk = ((double) lo + (double) s + 1) / 2;
        for (size_t t = lo; t < s; t++) d[ord[t]] = rk;
        lo = s;
      }
    }

    if (!standardize(d, N)) return NAN;

    for (int c = 0; c < be_m; c++) {
      const double *t = be_t[c];
      double acc = 0;
      for (size_t s = 0; s < N; s++) acc += d[s] * t[s];
      r[c] = (acc > 1) ? 1 : (acc < -1) ? -1 : acc;
    }
  }

  if (be_m == 1) return r[0];

  // Partial: x the community distances, y the subset's,
  // z the given ones.
  double rxy = r[0], ryz = r[1], rxz = be_rxz;
  if (fabs(ryz) >= 1 - 1e-12) return NAN;
  return (rxy - rxz * ryz) / sqrt(1 - rxz * rxz) / sqrt(1 - ryz * ryz);
}


//------------------------------------------------------
// Worker: a contiguous run of this round's subsets.
//------------------------------------------------------
static void *bioenv_worker(void *arg) {

  worker_t *w   = (worker_t *) arg;
  double   *d   = be_dwork + (size_t) w->i * be_dwlen;
  double   *mw  = d + be_np;
  int      *iw  = be_iwork + (size_t) w->i * be_iwlen;
  int      *cmb = iw + 3 * be_np;
  int       k   = be_k;

  double tot = be_r1 - be_r0;
  double lo  = be_r0 + floor(tot * w->i / w->n);
  double hi  = be_r0 + floor(tot * (w->i + 1) / w->n);

  double best_r = -INFINITY, best_i = -1;

  if (lo < hi) unrank(lo, k, cmb);

  for (double r = lo; r < hi; r++) {
    double c = subset_cor(cmb, k, d, mw, iw);
    if (c > best_r) { best_r = c; best_i = r; }
    next_comb(cmb, k, be_p);
  }

  be_best_r[w->i] = best_r;
  be_best_i[w->i] = best_i;

  return NULL;
}


//======================================================
// R interface.
//
// z:        n x p matrix of the standardized variables.
// targets:  List of 1 or 2 targets, in dist order.
//           Pearson / Spearman: standardized values or
//           ranks (double). Kendall: codes 0 .. K - 1 of
//           the distinct values (integer).
// ties:     Kendall: tied pairs in each target; else NULL.
// rxz:      With 2 targets, their correlation.
// opts:     c(method, metric, max_vars, n_threads).
//
// Returns list(r, best): the best correlation of each
// size, and a max_vars x max_vars integer matrix whose
// column k holds that size's subset (1-based variable
// indices, NA padded). Sizes with no correlation (every
// subset undefined) have r NA and an NA subset.
//======================================================
SEXP C_bioenv(SEXP sexp_z, SEXP sexp_targets, SEXP sexp_ties, SEXP sexp_rxz, SEXP sexp_opts) {

  int n         = nrows(sexp_z);
  int p         = ncols(sexp_z);
  int method    = INTEGER(sexp_opts)[0];
  int metric    = INTEGER(sexp_opts)[1];
  int max_vars  = INTEGER(sexp_opts)[2];
  int n_threads = INTEGER(sexp_opts)[3];
  int m         = LENGTH(sexp_targets);

  size_t N = (size_t) n * (size_t) (n - 1) / 2;

  if (n < 3 || p < 1 || max_vars < 1 || max_vars > p || n_threads < 1 ||
      m < 1 || m > 2 || N > INT_MAX / 3)
    error("Invalid arguments to C_bioenv()."); // # nocov

  be_n      = n;
  be_np     = N;
  be_p      = p;
  be_m      = m;
  be_method = method;
  be_metric = metric;
  be_rxz    = (m == 2) ? asReal(sexp_rxz) : 0;

  // Row-major copy of z, so that each sample's values are
  // together.
  const double *zc = REAL(sexp_z);
  double       *zr = (double *) R_alloc((size_t) n * p, sizeof(double));
  for (int i = 0; i < n; i++)
    for (int v = 0; v < p; v++)
      zr[(size_t) i * p + v] = zc[(size_t) v * n + i];
  be_z = zr;

  for (int c = 0; c < m; c++) {
    SEXP t = VECTOR_ELT(sexp_targets, c);
    if ((size_t) XLENGTH(t) != N) error("Invalid arguments to C_bioenv()."); // # nocov
    if (method == BE_KENDALL) {
      be_tc[c] = INTEGER(t);
      be_n2[c] = REAL(sexp_ties)[c];
    }
    else {
      be_t[c] = REAL(t);
    }
  }

  // Binomial coefficients, exact while below 2^53.
  double *binom = (double *) R_alloc((size_t) (p + 1) * (p + 1), sizeof(double));
  for (int a = 0; a <= p; a++)
    for (int b = 0; b <= p; b++)
      binom[(size_t) a * (p + 1) + b] =
        (b == 0) ? 1 : (a == 0 || b > a) ? 0 :
        binom[(size_t) (a - 1) * (p + 1) + b - 1] + binom[(size_t) (a - 1) * (p + 1) + b];
  be_binom = binom;

  for (int k = 1; k <= max_vars; k++)
    if (BINOM(p, k) > 4503599627370496.0) // 2^52
      error("Too many subsets of %d variables to search.", k);

  // Workspaces.
  be_dwlen = N + (size_t) n * max_vars + (size_t) max_vars * max_vars + max_vars;
  be_iwlen = 3 * N + (size_t) max_vars;
  be_dwork = (double *) R_alloc(be_dwlen * n_threads, sizeof(double));
  be_iwork = (int *)    R_alloc(be_iwlen * n_threads, sizeof(int));
  be_best_r = (double *) R_alloc(n_threads, sizeof(double));
  be_best_i = (double *) R_alloc(n_threads, sizeof(double));

  SEXP sexp_r    = PROTECT(allocVector(REALSXP, max_vars));
  SEXP sexp_best = PROTECT(allocMatrix(INTSXP, max_vars, max_vars));
  double *out_r  = REAL(sexp_r);
  int    *best   = INTEGER(sexp_best);
  for (size_t u = 0; u < (size_t) max_vars * max_vars; u++) best[u] = NA_INTEGER;


  //--------------------------------------------------
  // Each size, in rounds of subsets, checking for
  // interrupts between them. Each round has at least
  // ~1e8 operations per thread, or every remaining one.
  //--------------------------------------------------
  for (int k = 1; k <= max_vars; k++) {

    be_k = k;
    double total    = BINOM(p, k);
    double per_subs = (double) N * (k + 1) + 1;
    if (method != BE_PEARSON)    per_subs += 2 * (double) N * log2((double) N + 2) * m;
    if (metric == BE_MAHALANOBIS) per_subs += ((double) n + N / 2) * k * k;

    double per_round = floor(1e8 * n_threads / per_subs);
    if (per_round < n_threads) per_round = n_threads;

    double size_r = -INFINITY, size_i = -1;

    for (be_r0 = 0; be_r0 < total; be_r0 = be_r1) {

      be_r1 = (total - be_r0 > per_round) ? be_r0 + per_round : total;

      // Threads only for enough work: at least ~1e6
      // operations each. run_parallel() starts them only
      // when told of at least 100 tasks.
      int nt = threads_for(n_threads, per_subs * (be_r1 - be_r0), 1e6);
      if (nt > be_r1 - be_r0) nt = (int) (be_r1 - be_r0);
      run_parallel(bioenv_worker, nt, (nt > 1) ? 100 : 1);

      // Threads' runs are in order: the first strictly
      // greater keeps the first-found rule.
      for (int t = 0; t < nt; t++)
        if (be_best_r[t] > size_r) { size_r = be_best_r[t]; size_i = be_best_i[t]; }

      if (be_r1 < total) R_CheckUserInterrupt();
    }

    if (size_i < 0) {
      out_r[k - 1] = NA_REAL;
    }
    else {
      int *cmb = be_iwork + 3 * be_np;
      unrank(size_i, k, cmb);
      out_r[k - 1] = size_r;
      for (int a = 0; a < k; a++) best[(size_t) (k - 1) * max_vars + a] = cmb[a] + 1;
    }
  }

  SEXP sexp_out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(sexp_out, 0, sexp_r);
  SET_VECTOR_ELT(sexp_out, 1, sexp_best);

  UNPROTECT(3);
  return sexp_out;
}
