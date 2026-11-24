// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Closed forms of rarefaction, from each sample's frequency counts f_k (the
// number of features with exactly k observations):
//
//   * Expected richness of a subsample of d of the sample's N observations,
//     drawn without replacement (Hurlbert 1971), and its variance (Heck,
//     van Belle & Simberloff 1975).
//   * Expected sample coverage of that subsample (Chao & Jost 2012; Chao et
//     al. 2014, Table 1), and at d = N the reference sample's estimate
//     (Chao et al. 2014, Eq. 12 and Appendix G).
//
// Written from the papers. See plans/rarefaction_curves.md.
//
//   C_rarefaction_curve()   value, sd and coverage at given depths
//   C_coverage_depth()      the smallest depth reaching a coverage
//
// rc_freq(), rc_absence(), rc_coverage() and rc_coverage_full() are
// declared in ecodive.h, for hill_extrapolation() to share.

#include "ecodive.h"


// The frequency counts of one sample's `nnz` non-zero values `vals` (whole
// numbers), written to `k` (the distinct abundances, increasing) and `f`
// (how many features have each). `buf` is scratch room for `nnz` doubles.
// Returns the number of distinct abundances, and sets *N (the total) and
// *S (the features present).

static int cmp_dbl (const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

int rc_freq (const double *vals, int nnz, double *buf, double *k, double *f,
             double *N, double *S) {
  int m = 0;
  for (int i = 0; i < nnz; i++) if (vals[i] > 0) buf[m++] = vals[i];
  qsort(buf, (size_t) m, sizeof(double), cmp_dbl);

  int    K   = 0;
  double tot = 0;
  for (int i = 0; i < m; i++) {
    tot += buf[i];
    if (K > 0 && k[K - 1] == buf[i]) { f[K - 1] += 1; }
    else                             { k[K] = buf[i]; f[K] = 1; K++; }
  }
  *N = tot;
  *S = (double) m;
  return K;
}


// The probabilities that a feature with U[i] of a sample's N observations
// is absent from a subsample of d of them (0 <= d <= N), drawn without
// replacement: A[i] = C(N - U[i], d) / C(N, d), the product over j = 0 ..
// U[i] - 1 of (N - d - j) / (N - j), and 0 when U[i] > N - d. `U` holds `u`
// positive whole numbers, increasing. One running product serves them all,
// which avoids differences of lchoose() of large numbers; once it
// underflows to 0, the rest are 0.
void rc_absence (double N, double d, const double *U, int u, double *A) {
  double p = 1, k = 0; // p = A(k)
  int    i = 0;
  for (; i < u; i++) {
    double target = U[i];
    if (target > N - d) break;
    while (k < target && p > 0) { p *= (N - d - k) / (N - k); k += 1; }
    if (!(p > 0)) break;
    A[i] = p;
  }
  for (; i < u; i++) A[i] = 0;
}


// Expected coverage of a subsample of d < N observations (Chao & Jost 2012;
// Chao et al. 2014, Table 1): 1 - sum_k f_k (k / N) C(N - k, d) / C(N - 1,
// d), from the absence probabilities `a` of the K abundances `k`
// (rc_absence()), since C(N - k, d) / C(N - 1, d) = a_k N / (N - d).
double rc_coverage (double N, double d, int K, const double *k,
                    const double *f, const double *a) {
  double s = 0;
  for (int i = 0; i < K; i++) s += f[i] * k[i] * a[i];
  return 1 - s / (N - d);
}


// The reference sample's estimated coverage, at d = N (Chao et al. 2014,
// Eq. 12): 1 - (f1 / N) A, with A = (N - 1) f1 / ((N - 1) f1 + 2 f2), or
// with f2 = 0, (N - 1)(f1 - 1) / ((N - 1)(f1 - 1) + 2) (their Appendix G).
double rc_coverage_full (double N, int K, const double *k, const double *f) {
  double f1 = 0, f2 = 0;
  for (int i = 0; i < K && k[i] <= 2; i++) {
    if (k[i] == 1) f1 = f[i];
    else           f2 = f[i];
  }
  if (f1 == 0) return 1;
  double A = (f2 > 0)
    ? (N - 1) * f1 / ((N - 1) * f1 + 2 * f2)
    : (N - 1) * (f1 - 1) / ((N - 1) * (f1 - 1) + 2);
  return 1 - f1 / N * A;
}




//======================================================
// Shared setup: each sample's frequency counts.
//======================================================

static int     n_samples;
static int     n_threads_;
static int    *pos_vec;      // ecomatrix: sample s's values are
static double *val_vec;      //   val_vec[pos_vec[s] .. pos_vec[s + 1] - 1]
static double *frq_k;        // Sample s's abundances: frq_k[pos_vec[s] ..
static double *frq_f;        //   pos_vec[s] + frq_K[s] - 1], and counts.
static int    *frq_K;
static double *frq_N;
static double *frq_S;

// The sorted values go in frq_f's slots (rc_freq()'s `buf`), and are
// compressed in place: rc_freq() writes the j-th distinct value's count only
// after reading the j-th sorted value, and j never passes the position read.
static void *freq_worker (void *arg) {
  int i = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n;
  for (int s = i; s < n_samples; s += n) {
    int p = pos_vec[s], nnz = pos_vec[s + 1] - p;
    frq_K[s] = rc_freq(val_vec + p, nnz, frq_f + p, frq_k + p, frq_f + p, frq_N + s, frq_S + s);
  }
  return NULL;
}

static void freq_setup (SEXP sexp_counts, SEXP sexp_margin, int n_threads) {
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  n_samples  = em->n_samples;
  pos_vec    = em->pos_vec;
  val_vec    = em->val_vec;
  n_threads_ = n_threads;
  int nnz    = em->nnz;
  frq_k = (double *) safe_malloc((size_t) (nnz + 1) * sizeof(double));
  frq_f = (double *) safe_malloc((size_t) (nnz + 1) * sizeof(double));
  frq_K = (int *)    safe_malloc((size_t) (n_samples + 1) * sizeof(int));
  frq_N = (double *) safe_malloc((size_t) (n_samples + 1) * sizeof(double));
  frq_S = (double *) safe_malloc((size_t) (n_samples + 1) * sizeof(double));
  run_parallel(freq_worker, n_threads, n_samples);
}

static int max_K (void) {
  int m = 1;
  for (int s = 0; s < n_samples; s++) if (frq_K[s] > m) m = frq_K[s];
  return m;
}




//======================================================
// Rarefaction curves.
//======================================================

static int    *row_pos;   // Sample s's rows: row_pos[s] .. row_pos[s + 1] - 1
static int    *depth_vec; // Each row's depth (1 <= d <= N).
static double *out_value, *out_sd, *out_cov;
static int     buf_K;     // Per-thread scratch, for buf_K abundances.
static double *buf_A;


// Heck et al.'s variance of the richness of a depth-d subsample:
//
//   Var = sum_i f_i a_i (1 - a_i) + sum_i f_i (f_i - 1) D(k_i, k_i)
//       + 2 sum_{j < i} f_i f_j D(k_i, k_j),
//
// where a_i = A[i] is the absence probability of k_i (rc_absence()), and
// D(k, l) = a_{k+l} - a_k a_l (two features with k and l observations are
// both absent, less the product of their absences). Computed as a
// difference, D loses most of its digits when d is small next to N. So
// instead, with g_m = a_{k+m} / a_k (the absence of m observations among
// the N - k left) and b_m = a_m, D(k, l) = a_k (g_l - b_l), and
// delta_m = g_m - b_m follows
//
//   delta_{m+1} = delta_m s_m + g_m (r_m - s_m),
//   r_m = (N - d - k - m) / (N - k - m),  s_m = (N - d - m) / (N - m),
//   r_m - s_m = -d k / ((N - k - m) (N - m)),
//
// whose terms all have the same sign, so nothing cancels. For abundance
// k_i, delta runs up to the largest k_j needed (k_i itself, or k_{i-1}
// when f_i = 1), so a depth costs at most the sum of the distinct
// abundances (at most N) steps.
static double heck_variance (double N, double d, int K, const double *k,
                             const double *f, const double *A) {
  double v1 = 0, vp = 0;
  for (int i = 0; i < K; i++) v1 += f[i] * A[i] * (1 - A[i]);

  for (int i = 0; i < K; i++) {
    int jmax = (f[i] >= 2) ? i : i - 1;
    if (jmax < 0 || !(A[i] > 0)) continue;
    double ki = k[i], g = 1, delta = 0, m = 0, pairs = 0;
    for (int j = 0; j <= jmax; j++) {
      for (; m < k[j]; m += 1) {
        double t = (g > 0) ? -(d * ki) / ((N - ki - m) * (N - m)) : 0;
        delta    = delta * ((N - d - m) / (N - m)) + g * t;
        double r = N - d - ki - m;
        g        = (r > 0) ? g * (r / (N - ki - m)) : 0;
      }
      pairs += ((j < i) ? 2 * f[j] : (f[i] - 1)) * delta;
    }
    vp += f[i] * A[i] * pairs;
  }

  return v1 + vp;
}


// Rows are shared out one by one (row r to thread r % n), so that a few
// samples with many depths use every thread too. Each row depends only on
// its sample and depth, so results are the same for any number of threads.
static int *row_sam;      // Each row's sample.

static void *curve_worker (void *arg) {

  int     t = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n;
  double *A = buf_A + (size_t) t * buf_K;
  int     n_rows = row_pos[n_samples];

  for (int r = t; r < n_rows; r += n) {

    int     s = row_sam[r];
    int     K = frq_K[s];
    double *k = frq_k + pos_vec[s], *f = frq_f + pos_vec[s];
    double  N = frq_N[s], S = frq_S[s];
    double  d = (double) depth_vec[r];

    rc_absence(N, d, k, K, A);

    // Hurlbert's expected richness (exactly 1 at d = 1).
    double e = 0;
    for (int i = 0; i < K; i++) e += f[i] * A[i];
    out_value[r] = (d > 1) ? S - e : 1;

    // Heck et al.'s standard deviation: 0 at d = 1 (one observation is
    // always one feature) and at d = N, exactly.
    double var = (d > 1 && d < N) ? heck_variance(N, d, K, k, f, A) : 0;
    out_sd[r]  = (var > 0) ? sqrt(var) : 0;

    // Expected coverage; the whole sample's estimate at d = N.
    out_cov[r] = (d < N) ? rc_coverage(N, d, K, k, f, A) : rc_coverage_full(N, K, k, f);
  }

  return NULL;
}


// counts, margin: the table. row_pos: n + 1 offsets into depths (whole
// numbers from 1 to the sample's total), sample by sample. Returns
// list(value, sd, coverage), one value per depth.
SEXP C_rarefaction_curve (SEXP sexp_counts, SEXP sexp_margin,
                          SEXP sexp_row_pos, SEXP sexp_depths,
                          SEXP sexp_n_threads) {

  init_n_ptrs(20);
  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  freq_setup(sexp_counts, sexp_margin, n_threads);

  row_pos   = INTEGER(sexp_row_pos);
  depth_vec = INTEGER(sexp_depths);
  int n_rows = LENGTH(sexp_depths);

  if (LENGTH(sexp_row_pos) != n_samples + 1 || row_pos[n_samples] != n_rows) {
    free_all();                                 // # nocov
    error("Rows don't match the samples.");     // # nocov
  }
  row_sam = (int *) safe_malloc((size_t) (n_rows + 1) * sizeof(int));
  for (int s = 0; s < n_samples; s++)
    for (int r = row_pos[s]; r < row_pos[s + 1]; r++) {
      row_sam[r] = s;
      if (depth_vec[r] < 1 || (double) depth_vec[r] > frq_N[s]) {
        free_all();                             // # nocov
        error("A depth is outside 1 .. the sample's total."); // # nocov
      }
    }

  // Scratch: only the samples with rows matter.
  buf_K = 1;
  for (int s = 0; s < n_samples; s++)
    if (row_pos[s + 1] > row_pos[s] && frq_K[s] > buf_K) buf_K = frq_K[s];
  buf_A = (double *) safe_malloc((size_t) n_threads * buf_K * sizeof(double));

  SEXP sexp_value = PROTECT(allocVector(REALSXP, n_rows));
  SEXP sexp_sd    = PROTECT(allocVector(REALSXP, n_rows));
  SEXP sexp_cov   = PROTECT(allocVector(REALSXP, n_rows));
  out_value = REAL(sexp_value);
  out_sd    = REAL(sexp_sd);
  out_cov   = REAL(sexp_cov);

  run_parallel(curve_worker, n_threads, n_rows);

  SEXP sexp_result = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(sexp_result, 0, sexp_value);
  SET_VECTOR_ELT(sexp_result, 1, sexp_sd);
  SET_VECTOR_ELT(sexp_result, 2, sexp_cov);
  SEXP sexp_names = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(sexp_names, 0, mkChar("value"));
  SET_STRING_ELT(sexp_names, 1, mkChar("sd"));
  SET_STRING_ELT(sexp_names, 2, mkChar("coverage"));
  setAttrib(sexp_result, R_NamesSymbol, sexp_names);

  free_all();
  UNPROTECT(5);
  return sexp_result;
}




//======================================================
// The depth at which each sample reaches a coverage.
//======================================================

static double  target_cov;
static int    *out_depth;

// Expected coverage at depth d (1 <= d <= N) of a sample.
static double coverage_at (double N, double d, int K, const double *k,
                           const double *f, double *A) {
  if (d >= N) return rc_coverage_full(N, K, k, f);
  rc_absence(N, d, k, K, A);
  return rc_coverage(N, d, K, k, f, A);
}

static void *depth_worker (void *arg) {
  int     t = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n;
  double *A = buf_A + (size_t) t * buf_K;
  for (int s = t; s < n_samples; s += n) {
    int     K = frq_K[s];
    double *k = frq_k + pos_vec[s], *f = frq_f + pos_vec[s];
    double  N = frq_N[s];
    // Unreachable (or an empty sample): one more than the total.
    if (N < 1 || coverage_at(N, N, K, k, f, A) < target_cov) {
      out_depth[s] = (int) N + 1;
      continue;
    }
    // Expected coverage increases with depth: the smallest d with
    // coverage >= target, by bisection over 1 .. N.
    double lo = 1, hi = N;
    while (lo < hi) {
      double mid = floor((lo + hi) / 2);
      if (coverage_at(N, mid, K, k, f, A) >= target_cov) hi = mid;
      else                                               lo = mid + 1;
    }
    out_depth[s] = (int) lo;
  }
  return NULL;
}


// The smallest depth (1 .. N) whose expected coverage is at least
// `coverage`, for each sample; the sample's total plus one where even the
// whole sample's estimated coverage falls short.
SEXP C_coverage_depth (SEXP sexp_counts, SEXP sexp_margin,
                       SEXP sexp_coverage, SEXP sexp_n_threads) {

  init_n_ptrs(20);
  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;
  target_cov = asReal(sexp_coverage);

  freq_setup(sexp_counts, sexp_margin, n_threads);

  buf_K = max_K();
  buf_A = (double *) safe_malloc((size_t) n_threads * buf_K * sizeof(double));

  SEXP sexp_result = PROTECT(allocVector(INTSXP, n_samples));
  out_depth = INTEGER(sexp_result);

  run_parallel(depth_worker, n_threads, n_samples);

  free_all();
  UNPROTECT(1);
  return sexp_result;
}
