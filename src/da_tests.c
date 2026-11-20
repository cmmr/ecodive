// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The per-feature statistics of da_wilcoxon() and da_ttest(), and the
// normalized table by feature that they and da_lm() (feature_lm.c) share.
//
// da_table_init() normalizes the table per sample (normalize.c, without a
// dense copy), then transposes the non-zero values to feature order: each
// task fills one feature's column, of every sample, from its stored values
// and each sample's value for a zero (0, the clr's or log's value of the
// pseudocount, or NaN for rclr, whose zeros are missing).
//
// The statistics are what base R's wilcox.test() and t.test() compute from
// them: Wilcoxon's rank sum W (average ranks for ties by exact equality),
// the ties' sum of t^3 - t, and each side's median; or each side's mean and
// variance, by R's own algorithms (mean() and var(), with long double
// sums). The p-values are computed in R, on the main thread. Each feature
// depends only on its own column, so results are identical for any number
// of threads.

#include "ecodive.h"

#define NORM_NONE 0
#define NORM_CLR  2
#define NORM_RCLR 5
#define NORM_LOG  7


//======================================================
// The normalized table by feature.
//======================================================
void da_table_init(
    SEXP sexp_counts, SEXP sexp_margin, int norm, double pc, int n_threads,
    da_table_t *t ) {

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  if (norm != NORM_NONE) normalize(em, norm, n_threads, pc);

  int n   = em->n_samples;
  int D   = em->n_otus;
  int nnz = em->nnz;

  int    *fpos = (int *)    R_alloc((size_t) D + 1,   sizeof(int));
  int    *fsam = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  double *fval = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  double *z0   = (double *) R_alloc((size_t) n + 1,   sizeof(double));

  // A stable counting sort by feature keeps each feature's
  // samples in order.
  memset(fpos, 0, ((size_t) D + 1) * sizeof(int));
  for (int e = 0; e < nnz; e++) fpos[em->otu_vec[e] + 1]++;
  for (int f = 0; f < D; f++)   fpos[f + 1] += fpos[f];

  int *fill = (int *) R_alloc((size_t) D + 1, sizeof(int));
  memcpy(fill, fpos, ((size_t) D + 1) * sizeof(int));
  for (int s = 0; s < n; s++)
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++) {
      int d = fill[em->otu_vec[e]]++;
      fsam[d] = s;
      fval[d] = em->val_vec[e];
    }

  for (int s = 0; s < n; s++) {
    if      (norm == NORM_RCLR) z0[s] = R_NaN;
    else if ((norm == NORM_CLR || norm == NORM_LOG) && em->clr_vec)
      z0[s] = em->clr_vec[s];
    else    z0[s] = 0;
  }

  free_all();

  t->n    = n;
  t->D    = D;
  t->fpos = fpos;
  t->fsam = fsam;
  t->fval = fval;
  t->z0   = z0;
}


// Feature f's value in every sample, in table order.
void da_table_column(const da_table_t *t, int f, double *y) {
  memcpy(y, t->z0, (size_t) t->n * sizeof(double));
  for (int e = t->fpos[f]; e < t->fpos[f + 1]; e++) y[t->fsam[e]] = t->fval[e];
}



//======================================================
// The tests' shared state.
//======================================================
typedef struct { double v; int side; } dt_pair_t;

static const da_table_t *dt_tab;
static int           dt_test;   // 1 Wilcoxon, 2 t-test
static const int    *dt_side;   // per sample: 1, 2, or 0 (left out)
static double       *dt_buf;    // per thread: n + 2 * m doubles
static dt_pair_t     *dt_pairs;  // per thread: m + 1 (value, side) pairs
static int           dt_m;      // samples in the contrast
static double       *dt_out;    // D x 7, column-major


static int cmp_pair (const void *a, const void *b) {
  double x = ((const dt_pair_t *) a)->v, y = ((const dt_pair_t *) b)->v;
  return (x > y) - (x < y);
}


// R's mean(): the long double mean, refined by the mean of the residuals.
static double r_mean (const double *x, int n) {
  long double s = 0;
  for (int i = 0; i < n; i++) s += x[i];
  s /= n;
  if (R_FINITE((double) s)) {
    long double t = 0;
    for (int i = 0; i < n; i++) t += (x[i] - s);
    s += t / n;
  }
  return (double) s;
}


// R's var(), with its own (long double) mean.
static double r_var (const double *x, int n) {
  if (n < 2) return R_NaN;
  long double sum = 0;
  for (int i = 0; i < n; i++) sum += x[i];
  long double tmp = sum / n;
  if (R_FINITE((double) tmp)) {
    sum = 0;
    for (int i = 0; i < n; i++) sum += (x[i] - tmp);
    tmp = tmp + sum / n;
  }
  double xm = (double) tmp;
  sum = 0;
  for (int i = 0; i < n; i++) sum += (x[i] - xm) * (x[i] - xm);
  return (double) (sum / (n - 1));
}


// R's median() of sorted values.
static double sorted_median (const double *x, int n) {
  if (n == 0) return R_NaN;
  int half = (n + 1) / 2;
  if (n % 2 == 1) return x[half - 1];
  double a = x[half - 1], b = x[half];
  long double s = ((long double) a + b) / 2;
  if (R_FINITE((double) s)) s += ((a - s) + (b - s)) / 2;
  return (double) s;
}


static void *dt_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  const da_table_t *tab = dt_tab;
  int n = tab->n, D = tab->D, m = dt_m;

  double    *y    = dt_buf + (size_t) thread_i * ((size_t) n + 3 * (size_t) m + 2);
  double    *a    = y + n;          // treatment side's values
  double    *b    = a + m + 1;      // reference side's values
  dt_pair_t *pr   = dt_pairs + (size_t) thread_i * ((size_t) m + 1);

  for (int f = thread_i; f < D; f += n_threads) {

    da_table_column(tab, f, y);

    int n1 = 0, n2 = 0;
    for (int s = 0; s < n; s++) {
      int sd = dt_side[s];
      if (!sd || ISNAN(y[s])) continue;
      if (sd == 1) a[n1++] = y[s]; else b[n2++] = y[s];
    }

    double *out = dt_out + f;
    out[0]             = n1;
    out[(size_t) D]    = n2;

    if (dt_test == 1) {

      // Wilcoxon: ranks of the pooled values, ties averaged.
      int k = 0;
      for (int i = 0; i < n1; i++) { pr[k].v = a[i]; pr[k++].side = 1; }
      for (int i = 0; i < n2; i++) { pr[k].v = b[i]; pr[k++].side = 2; }
      qsort(pr, (size_t) k, sizeof(dt_pair_t), cmp_pair);

      double r1 = 0, ties = 0;
      int i1 = 0, i2 = 0;
      for (int i = 0; i < k; ) {
        int j = i + 1;
        while (j < k && pr[j].v == pr[i].v) j++;
        double rank = (i + 1 + j) / 2.0;
        double tt   = j - i;
        if (tt > 1) ties += tt * tt * tt - tt;
        for (int e = i; e < j; e++) {
          if (pr[e].side == 1) { r1 += rank; a[i1++] = pr[e].v; }
          else                 {             b[i2++] = pr[e].v; }
        }
        i = j;
      }

      out[2 * (size_t) D] = r1 - (double) n1 * (n1 + 1) / 2;
      out[3 * (size_t) D] = ties;
      out[4 * (size_t) D] = sorted_median(a, n1);
      out[5 * (size_t) D] = sorted_median(b, n2);
      out[6 * (size_t) D] = R_NaN;
    }
    else {

      out[2 * (size_t) D] = n1 ? r_mean(a, n1) : R_NaN;
      out[3 * (size_t) D] = n2 ? r_mean(b, n2) : R_NaN;
      out[4 * (size_t) D] = r_var(a, n1);
      out[5 * (size_t) D] = r_var(b, n2);

      // The mean of both sides' values together.
      memcpy(a + n1, b, (size_t) n2 * sizeof(double));
      out[6 * (size_t) D] = (n1 + n2) ? r_mean(a, n1 + n2) : R_NaN;
    }
  }

  return NULL;
}



//======================================================
// R interface: da_wilcoxon() (test 1) and da_ttest()
// (test 2). Returns a D x 7 matrix: n1, n2 (the values
// on each side), then
//   Wilcoxon: W, sum(t^3 - t) over ties, median1,
//             median2, NaN;
//   t-test:   mean1, mean2, var1, var2, mean of both.
//======================================================
SEXP C_da_tests(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_pseudocount,
    SEXP sexp_side, SEXP sexp_test, SEXP sexp_n_threads ) {

  int    norm  = asInteger(sexp_norm);
  double pc    = asReal(sexp_pseudocount);
  int    n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, norm, pc, n_thr, &tab);
  int n = tab.n, D = tab.D;

  if (LENGTH(sexp_side) != n) error("`side` must have one value per sample.");  // # nocov

  int m = 0;
  for (int s = 0; s < n; s++) if (INTEGER(sexp_side)[s]) m++;

  // Threads only for enough work: about 20 ns per value.
  int nt = threads_for(n_thr, (double) D * (n + m) * 20, 1e6);

  SEXP result = PROTECT(allocMatrix(REALSXP, D, 7));

  dt_tab  = &tab;
  dt_test = asInteger(sexp_test);
  dt_side = INTEGER(sexp_side);
  dt_m    = m;
  dt_out  = REAL(result);
  dt_buf  = (double *) R_alloc((size_t) nt * ((size_t) n + 3 * (size_t) m + 2), sizeof(double));
  dt_pairs = (dt_pair_t *) R_alloc((size_t) nt * ((size_t) m + 1), sizeof(dt_pair_t));

  run_parallel(dt_worker, nt, D < 100 ? 100 : D);

  UNPROTECT(1);
  return result;
}
