// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

// Helpers for beta_div_rarefied() and alpha_div_rarefied() (R/rarefied.r):
//
//   C_row_medians()   the median of each row of a matrix, as
//                     apply(M, 1, median) gives it.
//   C_medoid_accum()  adds, to each cell (i, j) of a times x times matrix,
//                     the squared differences between columns i and j of a
//                     block of rows: the medoid's squared Euclidean
//                     distances between iterations (q2-boots' _medoid()),
//                     accumulated block by block.
//
// Each row (medians) or each pair of columns (medoid) is computed on its own,
// in a fixed order, so results are the same for any number of threads, and
// the medoid's sums don't depend on how the rows are split into blocks.

#include "ecodive.h"


/*
 * Row medians
 */

static const double *rm_mtx;
static double       *rm_res;
static double       *rm_buf;
static int           rm_nrow;
static int           rm_ncol;


// Rearranges x[lo..hi] so that x[k] holds the value it would hold if
// x[lo..hi] were sorted, with no larger value before it and no smaller value
// after it (Hoare's selection, with a median-of-three pivot).
static void select_kth (double *x, int lo, int hi, int k) {

  while (hi > lo) {

    int    mid = lo + (hi - lo) / 2;
    double a = x[lo], b = x[mid], c = x[hi];
    double pivot = (a < b) ? ((b < c) ? b : ((a < c) ? c : a))
                           : ((a < c) ? a : ((b < c) ? c : b));

    int i = lo, j = hi;
    while (i <= j) {
      while (x[i] < pivot) i++;
      while (x[j] > pivot) j--;
      if (i <= j) {
        double t = x[i]; x[i] = x[j]; x[j] = t;
        i++; j--;
      }
    }

    if      (k <= j) hi = j;
    else if (k >= i) lo = i;
    else             return;
  }
}


static void *row_medians_worker (void *arg) {

  int thread_i  = ((worker_t *)arg)->i;
  int n_threads = ((worker_t *)arg)->n;

  int     n   = rm_ncol;
  double *buf = rm_buf + (size_t) thread_i * n;

  for (int r = thread_i; r < rm_nrow; r += n_threads) {

    // NA (or NaN) anywhere gives NA, as median() without na.rm does.
    int has_na = 0;
    for (int c = 0; c < n; c++) {
      double v = rm_mtx[r + (size_t) c * rm_nrow];
      if (ISNAN(v)) { has_na = 1; break; }
      buf[c] = v;
    }
    if (has_na) { rm_res[r] = NA_REAL; continue; }

    int half = (n - 1) / 2;
    select_kth(buf, 0, n - 1, half);
    double m = buf[half];

    if (n % 2 == 0) {
      // The next value up: the smallest of those after `half`.
      double m2 = buf[half + 1];
      for (int c = half + 2; c < n; c++)
        if (buf[c] < m2) m2 = buf[c];
      // The same as R's mean() of the two: halving a rounded sum is exact.
      m = (m + m2) / 2;
    }

    rm_res[r] = m;
  }

  return NULL;
}


SEXP C_row_medians (SEXP sexp_mtx, SEXP sexp_n_threads) {

  rm_nrow = nrows(sexp_mtx);
  rm_ncol = ncols(sexp_mtx);
  rm_mtx  = REAL(sexp_mtx);

  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  SEXP sexp_res = PROTECT(allocVector(REALSXP, rm_nrow));
  rm_res = REAL(sexp_res);

  if (rm_ncol == 0) {
    for (int r = 0; r < rm_nrow; r++) rm_res[r] = NA_REAL;
    UNPROTECT(1);
    return sexp_res;
  }

  init_n_ptrs(1);
  rm_buf = (double*) safe_malloc((size_t) n_threads * rm_ncol * sizeof(double));

  run_parallel(row_medians_worker, n_threads, rm_nrow);

  free_all();
  UNPROTECT(1);
  return sexp_res;
}



/*
 * Medoid: summed squared differences between iterations
 */

static const double *ma_mtx;
static double       *ma_acc;
static int           ma_nrow;
static int           ma_ncol;


static void *medoid_accum_worker (void *arg) {

  int thread_i  = ((worker_t *)arg)->i;
  int n_threads = ((worker_t *)arg)->n;

  int    n_pairs = ma_ncol * (ma_ncol - 1) / 2;
  int    i = 0, j = 1; // Column pair of index p (i < j), in order.
  int    p_next = 0;

  for (int p = thread_i; p < n_pairs; p += n_threads) {

    // Walk forward from the last pair to pair p.
    for (; p_next < p; p_next++) {
      if (++j == ma_ncol) { i++; j = i + 1; }
    }

    const double *x = ma_mtx + (size_t) i * ma_nrow;
    const double *y = ma_mtx + (size_t) j * ma_nrow;

    // Continue the running sum from earlier blocks, row by row, so that the
    // total is the same however the rows are split into blocks. Rows with
    // NA or NaN in either column are skipped.
    double acc = ma_acc[i + (size_t) j * ma_ncol];
    for (int r = 0; r < ma_nrow; r++) {
      double d = x[r] - y[r];
      if (!ISNAN(d)) acc += d * d;
    }

    ma_acc[i + (size_t) j * ma_ncol] = acc;
    ma_acc[j + (size_t) i * ma_ncol] = acc;
  }

  return NULL;
}


SEXP C_medoid_accum (SEXP sexp_mtx, SEXP sexp_acc, SEXP sexp_n_threads) {

  ma_nrow = nrows(sexp_mtx);
  ma_ncol = ncols(sexp_mtx);
  ma_mtx  = REAL(sexp_mtx);

  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  SEXP sexp_res = PROTECT(duplicate(sexp_acc));
  ma_acc = REAL(sexp_res);

  int n_pairs = ma_ncol * (ma_ncol - 1) / 2;
  run_parallel(medoid_accum_worker, n_threads, n_pairs);

  UNPROTECT(1);
  return sexp_res;
}
