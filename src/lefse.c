// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// LEfSe, lefse(), is a port of LEfSe 1.1.2 (lefse/lefse.py,
// lefse_run.py and lefse_format_input.py), under this license:
//
//   Copyright (c) 2012, Nicola Segata and Curtis Huttenhower
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//
//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// LEfSe calls R through rpy2 for three steps: kruskal.test(), coin's
// wilcox_test() and MASS's lda(). Those packages are GPL, so their code
// was not read: the statistics here are reproduced from their documented
// behaviour and outputs (rank_tests.c's Kruskal-Wallis H; the rank-sum
// statistic, whose normal approximation without continuity correction is
// coin's asymptotic p-value, computed in R; and the discriminant below).
//
// The entry points read the table by feature (da_tests.c's
// da_table_t) and, unless `total` is NA, scaling each sample to sum to
// `total` as lefse_format_input.py's numerical_values() does (value *
// (total / sum), the sum in feature order; a zero-total sample stays 0;
// and a feature whose relative standard deviation is below 1e-10 is
// rounded to 6 decimals):
//
//   C_lefse_kw         Each feature's class means and Kruskal-Wallis H.
//   C_lefse_columns    The rescaled values of the features to score.
//   C_lefse_wilcoxon   For the given features and pairs of subgroups: the
//                      rank sum W, the ties' sum of t^3 - t, both medians
//                      (numpy's), and whether both subgroups hold one and
//                      the same value.
//   C_lefse_lda        The LDA effect sizes of those features: the
//                      jitter of test_lda_r(), then n_boots bootstrap
//                      subsamples, each fitted by the discriminant below.
//
// The discriminant is MASS::lda(tol = 1e-10)'s first axis as its outputs
// show it (and as lda.c computes it for full-rank tables, following
// scikit-learn): each feature scaled by its pooled within-group standard
// deviation (divisor m - 1); a feature whose deviation is below tol is
// left out (MASS stops: "appears to be constant within groups"); the SVD of
// the scaled residuals Z / sqrt(m - g) keeps the directions whose singular
// value is above tol (MASS warns "variables are collinear" otherwise), which
// handles more features than samples; and the first right singular vector
// of the sphered, prior-weighted group means gives the coefficients w.
// With u = w / |w|, each pair of classes (a, b) adds, for each feature j,
// (|m_aj - m_bj| + |u_j| |(m_a - m_b)' u|) / 2, with the subsample's class
// means (0 for a class that is absent, whose effect counts as 0).
//
// Threads: features (C_lefse_kw, C_lefse_wilcoxon) or bootstraps
// (C_lefse_lda). Bootstrap i draws its samples from PCG stream (seed, i),
// and the jitter of feature f (its index in the table) from stream (seed,
// 2^32 + f); each bootstrap's sums are kept apart and averaged in order, so
// results are identical for any number of threads.

#include "ecodive.h"

#define LF_TOL      1e-10   // LEfSe's tol for MASS::lda().
#define LF_TRIES    1000    // Redraws of a bootstrap subsample.
#define LF_JITTER   ((uint64_t) 1 << 32)



//======================================================
// The scaled table, by feature.
//======================================================
typedef struct {
  da_table_t tab;
  double    *mul;     // n multipliers, or NULL: values as given.
} lf_table_t;

static void lf_init (SEXP sexp_counts, SEXP sexp_margin, double total, int n_threads, lf_table_t *t) {

  da_table_init(sexp_counts, sexp_margin, 0, 0, n_threads, &t->tab);
  t->mul = NULL;
  if (ISNAN(total)) return;

  int n = t->tab.n, D = t->tab.D;
  double *sum = (double *) R_alloc((size_t) n + 1, sizeof(double));
  memset(sum, 0, ((size_t) n + 1) * sizeof(double));
  for (int f = 0; f < D; f++)
    for (int e = t->tab.fpos[f]; e < t->tab.fpos[f + 1]; e++)
      sum[t->tab.fsam[e]] += t->tab.fval[e];
  for (int s = 0; s < n; s++) sum[s] = (sum[s] == 0) ? 0.0 : total / sum[s];
  t->mul = sum;
}

static void lf_column (const lf_table_t *t, int f, double *y) {

  int n = t->tab.n;
  da_table_column(&t->tab, f, y);
  if (!t->mul) return;

  for (int s = 0; s < n; s++) y[s] *= t->mul[s];

  // numerical_values(): a (nearly) constant feature is rounded.
  long double m = 0, v = 0;
  for (int s = 0; s < n; s++) m += y[s];
  m /= n;
  for (int s = 0; s < n; s++) v += (y[s] - m) * (y[s] - m);
  double mean = (double) m, sd = sqrt((double) (v / n));
  if (mean != 0 && sd / mean < 1e-10)
    for (int s = 0; s < n; s++) y[s] = nearbyint(y[s] * 1e6) / 1e6;
}



//======================================================
// Class means and Kruskal-Wallis H of every feature.
//======================================================
static const lf_table_t *kw_tab;
static const int        *kw_cls;
static int               kw_k;
static const int        *kw_nc;
static double           *kw_buf;     // per thread: n + k doubles
static vg_t             *kw_vg;      // per thread: n pairs
static long double      *kw_ld;      // per thread: k sums
static double           *kw_stat, *kw_means;

static void *kw_worker (void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = kw_tab->tab.n, D = kw_tab->tab.D, k = kw_k;

  double *y    = kw_buf + (size_t) thread_i * (n + k);
  double *rsum = y + n;
  vg_t   *vg   = kw_vg + (size_t) thread_i * n;

  for (int f = thread_i; f < D; f += n_threads) {
    lf_column(kw_tab, f, y);

    long double *sm = kw_ld + (size_t) thread_i * k;
    for (int c = 0; c < k; c++) sm[c] = 0;
    for (int s = 0; s < n; s++) sm[kw_cls[s]] += y[s];
    for (int c = 0; c < k; c++) kw_means[c + (size_t) f * k] = (double) (sm[c] / kw_nc[c]);

    for (int s = 0; s < n; s++) { vg[s].v = y[s]; vg[s].g = kw_cls[s]; }
    kw_stat[f] = rank_kruskal_h(vg, n, kw_nc, k, rsum);
  }
  return NULL;
}


// cls: 0-based class of each sample; k classes, none empty.
// Returns list(stat = D values of H, means = k x D).
SEXP C_lefse_kw (
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_total, SEXP sexp_cls,
    SEXP sexp_k, SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  lf_table_t t;
  lf_init(sexp_counts, sexp_margin, asReal(sexp_total), n_thr, &t);
  int n = t.tab.n, D = t.tab.D, k = asInteger(sexp_k);
  if (LENGTH(sexp_cls) != n) error("`cls` must have one value per sample.");  // # nocov

  int *nc = (int *) R_alloc((size_t) k, sizeof(int));
  memset(nc, 0, (size_t) k * sizeof(int));
  for (int s = 0; s < n; s++) nc[INTEGER(sexp_cls)[s]]++;

  SEXP stat  = PROTECT(allocVector(REALSXP, D));
  SEXP means = PROTECT(allocMatrix(REALSXP, k, D));

  int nt = threads_for(n_thr, (double) D * n * 30, 1e6);

  kw_tab   = &t;
  kw_cls   = INTEGER(sexp_cls);
  kw_k     = k;
  kw_nc    = nc;
  kw_stat  = REAL(stat);
  kw_means = REAL(means);
  kw_buf   = (double *) R_alloc((size_t) nt * (n + k), sizeof(double));
  kw_vg    = (vg_t *)   R_alloc((size_t) nt * n, sizeof(vg_t));
  kw_ld    = (long double *) R_alloc((size_t) nt * k, sizeof(long double));

  run_parallel(kw_worker, nt, D < 100 ? 100 : D);

  SEXP res = PROTECT(allocVector(VECSXP, 2));
  SEXP nms = PROTECT(allocVector(STRSXP, 2));
  SET_VECTOR_ELT(res, 0, stat);  SET_STRING_ELT(nms, 0, mkChar("stat"));
  SET_VECTOR_ELT(res, 1, means); SET_STRING_ELT(nms, 1, mkChar("means"));
  setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}



//======================================================
// The rescaled values of some features, dense (n x p),
// for the LDA.
//======================================================
SEXP C_lefse_columns (
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_total, SEXP sexp_feat,
    SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  lf_table_t t;
  lf_init(sexp_counts, sexp_margin, asReal(sexp_total), n_thr, &t);
  int n = t.tab.n, p = LENGTH(sexp_feat);

  SEXP x = PROTECT(allocMatrix(REALSXP, n, p));
  for (int j = 0; j < p; j++) lf_column(&t, INTEGER(sexp_feat)[j], REAL(x) + (size_t) j * n);
  UNPROTECT(1);
  return x;
}



//======================================================
// Rank sums between pairs of subgroups.
//======================================================
static int cmp_dbl (const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

// numpy.median() of x[0 .. n - 1] (sorted here).
static double np_median (double *x, int n) {
  qsort(x, (size_t) n, sizeof(double), cmp_dbl);
  if (n % 2) return x[n / 2];
  return (x[n / 2 - 1] + x[n / 2]) / 2;
}

static const lf_table_t *wx_tab;
static const int        *wx_feat, *wx_pos, *wx_sam, *wx_pa, *wx_pb;
static int               wx_nf, wx_np;
static double           *wx_buf;     // per thread: 3n doubles
static vg_t             *wx_vg;      // per thread: n pairs
static double           *wx_w, *wx_ties, *wx_m1, *wx_m2;
static int              *wx_same;

static void *wx_worker (void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = wx_tab->tab.n;

  double *y  = wx_buf + (size_t) thread_i * 3 * n;
  double *a  = y + n;
  double *b  = a + n;
  vg_t   *vg = wx_vg + (size_t) thread_i * n;

  for (int i = thread_i; i < wx_nf; i += n_threads) {
    lf_column(wx_tab, wx_feat[i], y);

    for (int p = 0; p < wx_np; p++) {
      size_t o = i + (size_t) p * wx_nf;
      int sa = wx_pa[p], sb = wx_pb[p];
      int n1 = wx_pos[sa + 1] - wx_pos[sa], n2 = wx_pos[sb + 1] - wx_pos[sb];
      for (int e = 0; e < n1; e++) a[e] = y[wx_sam[wx_pos[sa] + e]];
      for (int e = 0; e < n2; e++) b[e] = y[wx_sam[wx_pos[sb] + e]];

      // Both subgroups hold one and the same value.
      int same = (n1 > 0 && n2 > 0 && a[0] == b[0]);
      for (int e = 1; same && e < n1; e++) if (a[e] != a[0]) same = 0;
      for (int e = 1; same && e < n2; e++) if (b[e] != b[0]) same = 0;
      wx_same[o] = same;

      int N = 0;
      for (int e = 0; e < n1; e++) { vg[N].v = a[e]; vg[N++].g = 1; }
      for (int e = 0; e < n2; e++) { vg[N].v = b[e]; vg[N++].g = 0; }
      double ties;
      double r1 = rank_sum_ties(vg, N, &ties);
      wx_w[o]    = r1 - (double) n1 * (n1 + 1) / 2;
      wx_ties[o] = ties;
      wx_m1[o]   = n1 ? np_median(a, n1) : R_NaN;
      wx_m2[o]   = n2 ? np_median(b, n2) : R_NaN;
    }
  }
  return NULL;
}


// feat: 0-based features; sub: 0-based subgroup of each sample (-1: none),
// nsub subgroups; pa, pb: the pairs of subgroups to compare. Returns
// list(w, ties, med1, med2, same), each n_feat x n_pairs.
SEXP C_lefse_wilcoxon (
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_total, SEXP sexp_feat,
    SEXP sexp_sub, SEXP sexp_nsub, SEXP sexp_pa, SEXP sexp_pb, SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  lf_table_t t;
  lf_init(sexp_counts, sexp_margin, asReal(sexp_total), n_thr, &t);
  int n = t.tab.n, nsub = asInteger(sexp_nsub);
  int nf = LENGTH(sexp_feat), np = LENGTH(sexp_pa);
  const int *sub = INTEGER(sexp_sub);

  // Each subgroup's samples, in order.
  int *pos = (int *) R_alloc((size_t) nsub + 2, sizeof(int));
  int *sam = (int *) R_alloc((size_t) n + 1, sizeof(int));
  memset(pos, 0, ((size_t) nsub + 2) * sizeof(int));
  for (int s = 0; s < n; s++) if (sub[s] >= 0) pos[sub[s] + 1]++;
  for (int g = 0; g < nsub; g++) pos[g + 1] += pos[g];
  int *fill = (int *) R_alloc((size_t) nsub + 1, sizeof(int));
  memcpy(fill, pos, (size_t) nsub * sizeof(int));
  for (int s = 0; s < n; s++) if (sub[s] >= 0) sam[fill[sub[s]]++] = s;

  const char *names[5] = { "w", "ties", "med1", "med2", "same" };
  SEXP res = PROTECT(allocVector(VECSXP, 5));
  SEXP nms = PROTECT(allocVector(STRSXP, 5));
  for (int j = 0; j < 5; j++) {
    SET_VECTOR_ELT(res, j, allocMatrix(j == 4 ? LGLSXP : REALSXP, nf, np));
    SET_STRING_ELT(nms, j, mkChar(names[j]));
  }
  setAttrib(res, R_NamesSymbol, nms);

  int nt = threads_for(n_thr, (double) nf * np * n * 30, 1e6);

  wx_tab  = &t;
  wx_feat = INTEGER(sexp_feat);
  wx_pos  = pos;
  wx_sam  = sam;
  wx_pa   = INTEGER(sexp_pa);
  wx_pb   = INTEGER(sexp_pb);
  wx_nf   = nf;
  wx_np   = np;
  wx_w    = REAL(VECTOR_ELT(res, 0));
  wx_ties = REAL(VECTOR_ELT(res, 1));
  wx_m1   = REAL(VECTOR_ELT(res, 2));
  wx_m2   = REAL(VECTOR_ELT(res, 3));
  wx_same = LOGICAL(VECTOR_ELT(res, 4));
  wx_buf  = (double *) R_alloc((size_t) nt * 3 * n + 1, sizeof(double));
  wx_vg   = (vg_t *)   R_alloc((size_t) nt * n + 1, sizeof(vg_t));

  if (nf > 0 && np > 0) run_parallel(wx_worker, nt, nf < 100 ? 100 : nf);

  UNPROTECT(2);
  return res;
}



//======================================================
// The LDA effect size.
//======================================================

// Householder QR of the n x p A (n >= p), in place, serially; R (p x p,
// upper triangle) is written to R. v holds n doubles.
static void lf_qr (double *A, int n, int p, double *R, double *v) {

  memset(R, 0, (size_t) p * p * sizeof(double));

  for (int a = 0; a < p; a++) {
    double *col = A + (size_t) a * n + a;
    int     len = n - a;

    double norm = 0;
    for (int i = 0; i < len; i++) norm += col[i] * col[i];
    norm = sqrt(norm);

    double r = (col[0] > 0) ? -norm : norm;
    double vn = 0;
    if (norm > 0) {
      v[0] = col[0] - r;
      for (int i = 1; i < len; i++) v[i] = col[i];
      for (int i = 0; i < len; i++) vn += v[i] * v[i];
    }

    if (vn > 0) {
      double f = sqrt(2 / vn);
      for (int i = 0; i < len; i++) v[i] *= f;
      for (int j = a + 1; j < p; j++) {
        double *cj = A + (size_t) j * n + a;
        double  s  = 0;
        for (int i = 0; i < len; i++) s += v[i] * cj[i];
        for (int i = 0; i < len; i++) cj[i] -= s * v[i];
      }
    }
    else {
      r = 0;
    }

    R[a + (size_t) a * p] = r;
    for (int j = a + 1; j < p; j++) R[a + (size_t) j * p] = A[a + (size_t) j * n];
  }
}


static const double *ld_X;          // n x p, the (jittered) values
static const int    *ld_cls;        // n, 0-based classes
static const int    *ld_pa, *ld_pb; // the pairs of classes
static const int    *ld_boots;      // m x n_boots given samples, or NULL
static int           ld_n, ld_p, ld_k, ld_m, ld_np, ld_min_cl, ld_n_boots;
static uint64_t      ld_seed;
static int          *ld_drawn;      // m x n_boots: the samples used
static double       *ld_out;        // p x n_pairs x n_boots
static double       *ld_dbuf;       // per thread: LD_DOUBLES doubles
static int          *ld_ibuf;       // per thread: LD_INTS ints
static long double  *ld_acc;        // per thread: k sums

#define LD_C        ((size_t) (ld_m < ld_p ? ld_m : ld_p))
#define LD_DOUBLES  (3 * (size_t) ld_m * ld_p + 2 * LD_C * LD_C + 4 * (size_t) ld_p + \
                     (size_t) ld_k * ld_p + (size_t) ld_k * ld_k + 2 * (size_t) ld_m + \
                     LD_C * ld_k + LD_C + 3 * (size_t) ld_k + 10)
#define LD_INTS     (2 * (size_t) ld_m + 4 * (size_t) ld_k + (size_t) ld_p + 4)


// LEfSe's contast_within_classes_or_few_per_class(): whether the drawn
// samples miss a class, have fewer than min_cl of one, or some feature has
// at most min_cl distinct values (at most 1 when min_cl is 1) within a
// class. `ord` gets the draws ordered by class, `vals` m doubles.
static int ld_redraw (const int *idx, int *cnt, int *start, int *ord, double *vals) {

  int m = ld_m, k = ld_k;
  memset(cnt, 0, (size_t) k * sizeof(int));
  for (int j = 0; j < m; j++) cnt[ld_cls[idx[j]]]++;
  for (int c = 0; c < k; c++) if (cnt[c] == 0 || cnt[c] < ld_min_cl) return 1;

  start[0] = 0;
  for (int c = 0; c < k; c++) start[c + 1] = start[c] + cnt[c];
  int *fl = start + k + 1;
  for (int c = 0; c < k; c++) fl[c] = start[c];
  for (int j = 0; j < m; j++) ord[fl[ld_cls[idx[j]]]++] = idx[j];

  for (int f = 0; f < ld_p; f++) {
    const double *x = ld_X + (size_t) f * ld_n;
    for (int c = 0; c < k; c++) {
      int len = start[c + 1] - start[c];
      for (int e = 0; e < len; e++) vals[e] = x[ord[start[c] + e]];
      qsort(vals, (size_t) len, sizeof(double), cmp_dbl);
      int distinct = len ? 1 : 0;
      for (int e = 1; e < len; e++) if (vals[e] != vals[e - 1]) distinct++;
      if (distinct <= ld_min_cl) return 1;
    }
  }
  return 0;
}


static void *ld_worker (void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = ld_n, p = ld_p, k = ld_k, m = ld_m;
  size_t c_max = LD_C;

  double *S    = ld_dbuf + (size_t) thread_i * LD_DOUBLES;  // m x p subsample
  double *Z    = S  + (size_t) m * p;                       // m x q, or q x m
  double *sc   = Z  + (size_t) m * p;                       // q x r scaling
  double *R    = sc + (size_t) m * p;                       // c x c
  double *V    = R  + c_max * c_max;                        // c x c
  double *sd   = V  + c_max * c_max;                        // p
  double *u    = sd + p;                                    // p
  double *tmp  = u  + p;                                    // m + p + 2
  double *M    = tmp + m + p + 2;                           // k x p means
  double *Bt   = M  + (size_t) k * p;                       // r x g
  double *Vb   = Bt + c_max * k;                            // g x g
  double *d    = Vb + (size_t) k * k;                       // c
  double *db   = d  + c_max + 1;                            // k
  double *wcl  = db + k + 1;                                // k (sqrt(m prior / (g - 1)))
  double *xbar = wcl + k + 1;                               // p
  double *vals = xbar + p;                                  // m + 1
  long double *acc = ld_acc + (size_t) thread_i * k;        // k sums

  int *idx   = ld_ibuf + (size_t) thread_i * LD_INTS;       // m
  int *cnt   = idx + m;                                     // k
  int *start = cnt + k;                                     // k + 1 (+ k)
  int *gc    = start + 2 * k + 1;                           // k: compact codes
  int *incl  = gc + k;                                      // p
  int *ord   = incl + p;                                    // m

  for (int b = thread_i; b < ld_n_boots; b += n_threads) {

    // The subsample.
    if (ld_boots) {
      memcpy(idx, ld_boots + (size_t) b * m, (size_t) m * sizeof(int));
    }
    else {
      rng_t rng;
      rng_seed(&rng, ld_seed, (uint64_t) b);
      for (int t = 0; t < LF_TRIES; t++) {
        for (int j = 0; j < m; j++) idx[j] = (int) pcg32_bounded(&rng.pcg, (uint32_t) n);
        if (!ld_redraw(idx, cnt, start, ord, vals)) break;
      }
    }
    memcpy(ld_drawn + (size_t) b * m, idx, (size_t) m * sizeof(int));

    // Classes present, and their compact codes.
    memset(cnt, 0, (size_t) k * sizeof(int));
    for (int j = 0; j < m; j++) cnt[ld_cls[idx[j]]]++;
    int g = 0;
    for (int c = 0; c < k; c++) gc[c] = cnt[c] ? g++ : -1;

    // Subsample, class means (rows by class code, 0 when absent),
    // residuals and pooled sds (divisor m - 1).
    memset(M, 0, (size_t) k * p * sizeof(double));
    for (int f = 0; f < p; f++) {
      const double *x = ld_X + (size_t) f * n;
      double       *s = S + (size_t) f * m;
      for (int c = 0; c < k; c++) acc[c] = 0;
      for (int j = 0; j < m; j++) {
        s[j] = x[idx[j]];
        acc[ld_cls[idx[j]]] += s[j];
      }
      for (int c = 0; c < k; c++)
        if (cnt[c]) M[c + (size_t) f * k] = (double) (acc[c] / cnt[c]);
      double ss = 0;
      for (int j = 0; j < m; j++) {
        double r = s[j] - M[ld_cls[idx[j]] + (size_t) f * k];
        s[j] = r;
        ss  += r * r;
      }
      sd[f] = (m > 1) ? sqrt(ss / (m - 1)) : 0;
    }

    // The first discriminant direction u (0 where it can't be had).
    memset(u, 0, (size_t) p * sizeof(double));
    int q = 0;
    for (int f = 0; f < p; f++) if (sd[f] >= LF_TOL) incl[q++] = f;

    if (g >= 2 && m - g >= 1 && q > 0) {

      double fz = 1 / sqrt((double) (m - g));
      int r = 0;

      if (q <= m) {
        // Z (m x q), its QR, and the SVD of R.
        for (int a = 0; a < q; a++) {
          const double *s = S + (size_t) incl[a] * m;
          double       *z = Z + (size_t) a * m;
          double        h = fz / sd[incl[a]];
          for (int j = 0; j < m; j++) z[j] = s[j] * h;
        }
        lf_qr(Z, m, q, R, tmp);
        lda_jacobi_svd(R, q, q, V, d, tmp);
        for (int a = 0; a < q; a++) if (d[a] > LF_TOL) r++;
        for (int a = 0; a < r; a++)
          for (int i = 0; i < q; i++)
            sc[i + (size_t) a * q] = V[i + (size_t) a * q] / (sd[incl[i]] * d[a]);
      }
      else {
        // Z' (q x m): its columns become U D, Z's right singular
        // vectors times their singular values.
        for (int a = 0; a < q; a++) {
          const double *s = S + (size_t) incl[a] * m;
          double        h = fz / sd[incl[a]];
          for (int j = 0; j < m; j++) Z[a + (size_t) j * q] = s[j] * h;
        }
        lda_jacobi_svd(Z, q, m, V, d, tmp);
        for (int a = 0; a < m; a++) if (d[a] > LF_TOL) r++;
        for (int a = 0; a < r; a++)
          for (int i = 0; i < q; i++)
            sc[i + (size_t) a * q] = Z[i + (size_t) a * q] / (sd[incl[i]] * d[a] * d[a]);
      }

      if (r > 0) {
        // B' (r x g): the sphered group means about the weighted centre
        // xbar = sum_c prior_c M_c, prior_c = n_c / m.
        for (int i = 0; i < q; i++) {
          int f = incl[i];
          double xb = 0;
          for (int c = 0; c < k; c++)
            if (cnt[c]) xb += ((double) cnt[c] / m) * M[c + (size_t) f * k];
          xbar[f] = xb;
        }
        int gi = 0;
        for (int c = 0; c < k; c++) {
          if (gc[c] < 0) continue;
          wcl[gi++] = sqrt(m * ((double) cnt[c] / m) / (g - 1));
        }
        for (int a = 0; a < r; a++) {
          double *bt = Bt + a;
          gi = 0;
          for (int c = 0; c < k; c++) {
            if (gc[c] < 0) continue;
            double s = 0;
            for (int i = 0; i < q; i++)
              s += (M[c + (size_t) incl[i] * k] - xbar[incl[i]]) * sc[i + (size_t) a * q];
            bt[(size_t) gi * r] = wcl[gi] * s;
            gi++;
          }
        }
        lda_jacobi_svd(Bt, r, g, Vb, db, tmp);

        if (db[0] > 0) {
          // w = scaling v_1, with v_1 = B''s first column / db[0].
          double nw = 0;
          for (int i = 0; i < q; i++) {
            double s = 0;
            for (int a = 0; a < r; a++) s += sc[i + (size_t) a * q] * Bt[a];
            u[incl[i]] = s / db[0];
            nw += u[incl[i]] * u[incl[i]];
          }
          nw = sqrt(nw);
          if (nw > 0) for (int f = 0; f < p; f++) u[f] /= nw;
          else        memset(u, 0, (size_t) p * sizeof(double));
        }
      }
    }

    // Each pair's effect sizes.
    for (int pr = 0; pr < ld_np; pr++) {
      int ca = ld_pa[pr], cb = ld_pb[pr];
      double *out = ld_out + ((size_t) b * ld_np + pr) * p;
      double e = 0;
      int both = cnt[ca] && cnt[cb];
      if (both) {
        for (int f = 0; f < p; f++) e += (M[ca + (size_t) f * k] - M[cb + (size_t) f * k]) * u[f];
        e = fabs(e);
      }
      for (int f = 0; f < p; f++) {
        double gm = fabs(M[ca + (size_t) f * k] - M[cb + (size_t) f * k]);
        out[f] = (gm + (both ? fabs(u[f]) * e : 0)) * 0.5;
      }
    }
  }
  return NULL;
}


// x: n x p dense values of the features to score (already jittered when
// `jittered` is TRUE); feat: their 0-based indices in the table (for the
// jitter's streams); cls: 0-based classes (k); pa, pb: the pairs of
// classes; m samples per bootstrap; boots: NULL or m x n_boots 0-based
// samples. Returns list(means = p x n_pairs (each pair's effect, averaged
// over the bootstraps), x = the jittered values, boots = the samples
// drawn, m x n_boots, 0-based).
SEXP C_lefse_lda (
    SEXP sexp_x, SEXP sexp_feat, SEXP sexp_jittered, SEXP sexp_cls, SEXP sexp_k,
    SEXP sexp_pa, SEXP sexp_pb, SEXP sexp_n_boots, SEXP sexp_m, SEXP sexp_min_cl,
    SEXP sexp_boots, SEXP sexp_seed, SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  int n = nrows(sexp_x), p = ncols(sexp_x), k = asInteger(sexp_k);
  int n_boots = asInteger(sexp_n_boots), m = asInteger(sexp_m);
  int np = LENGTH(sexp_pa);
  const int *cls = INTEGER(sexp_cls);

  SEXP x = PROTECT(duplicate(sexp_x));
  double *X = REAL(x);

  // test_lda_r()'s jitter: in each class where a feature has at most
  // max(n_c / 2, 4) distinct values, each value v becomes
  // |v + N(0, max(0.05 v, 0.01))|.
  if (!asLogical(sexp_jittered)) {
    int *nc = (int *) R_alloc((size_t) k, sizeof(int));
    memset(nc, 0, (size_t) k * sizeof(int));
    for (int s = 0; s < n; s++) nc[cls[s]]++;
    double *vals = (double *) R_alloc((size_t) n + 1, sizeof(double));
    for (int f = 0; f < p; f++) {
      double *xf = X + (size_t) f * n;
      rng_t rng;
      rng_seed(&rng, (uint64_t) (uint32_t) asInteger(sexp_seed), LF_JITTER + (uint64_t) INTEGER(sexp_feat)[f]);
      for (int c = 0; c < k; c++) {
        int len = 0;
        for (int s = 0; s < n; s++) if (cls[s] == c) vals[len++] = xf[s];
        qsort(vals, (size_t) len, sizeof(double), cmp_dbl);
        int distinct = len ? 1 : 0;
        for (int e = 1; e < len; e++) if (vals[e] != vals[e - 1]) distinct++;
        double lim = nc[c] * 0.5;
        if (lim < 4) lim = 4;
        if (distinct > lim) continue;
        for (int s = 0; s < n; s++) {
          if (cls[s] != c) continue;
          double v  = xf[s];
          double sd = v * 0.05;
          if (sd < 0.01) sd = 0.01;
          xf[s] = fabs(v + rng_norm(&rng) * sd);
        }
      }
    }
  }

  SEXP means = PROTECT(allocMatrix(REALSXP, p, np));
  SEXP drawn = PROTECT(allocMatrix(INTSXP, m, n_boots));
  double *out = (double *) R_alloc((size_t) p * np * n_boots + 1, sizeof(double));

  // min_cl: LEfSe's minimum per class in a subsample.
  ld_X       = X;
  ld_cls     = cls;
  ld_pa      = INTEGER(sexp_pa);
  ld_pb      = INTEGER(sexp_pb);
  ld_boots   = isNull(sexp_boots) ? NULL : INTEGER(sexp_boots);
  ld_n       = n;
  ld_p       = p;
  ld_k       = k;
  ld_m       = m;
  ld_np      = np;
  ld_min_cl  = asInteger(sexp_min_cl);
  ld_n_boots = n_boots;
  ld_seed    = (uint64_t) (uint32_t) asInteger(sexp_seed);
  ld_drawn   = INTEGER(drawn);
  ld_out     = out;

  // Work per bootstrap: about the fit's m p min(m, p).
  double c = (m < p) ? m : p;
  int nt = threads_for(n_thr, (double) n_boots * m * p * (c + 10) * 5, 1e6);
  if (nt > n_boots) nt = n_boots;
  if (nt < 1) nt = 1;

  ld_dbuf = (double *) R_alloc((size_t) nt * LD_DOUBLES, sizeof(double));
  ld_ibuf = (int *)    R_alloc((size_t) nt * LD_INTS,    sizeof(int));
  ld_acc  = (long double *) R_alloc((size_t) nt * k, sizeof(long double));

  // (run_parallel() stays serial under 100 tasks; threads_for() has already
  // judged the work, and each thread takes the bootstraps i = t, t + nt, ...)
  if (n_boots > 0 && p > 0) run_parallel(ld_worker, nt, n_boots < 100 ? 100 : n_boots);

  // Averages over the bootstraps, in their order.
  double *mo = REAL(means);
  for (int pr = 0; pr < np; pr++)
    for (int f = 0; f < p; f++) {
      double s = 0;
      for (int b = 0; b < n_boots; b++) s += out[((size_t) b * np + pr) * p + f];
      mo[f + (size_t) pr * p] = s / n_boots;
    }

  SEXP res = PROTECT(allocVector(VECSXP, 3));
  SEXP nms = PROTECT(allocVector(STRSXP, 3));
  SET_VECTOR_ELT(res, 0, means); SET_STRING_ELT(nms, 0, mkChar("means"));
  SET_VECTOR_ELT(res, 1, x);     SET_STRING_ELT(nms, 1, mkChar("x"));
  SET_VECTOR_ELT(res, 2, drawn); SET_STRING_ELT(nms, 2, mkChar("boots"));
  setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(5);
  return res;
}
