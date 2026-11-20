// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// LinDA (Zhou et al. 2022): linear models of each feature's centred log2
// ratios, with the compositional bias of every coefficient estimated as
// the mode of that coefficient over the features.
//
// Written from the paper and from the documentation and outputs of
// MicrobiomeStat's linda() (GPL-3, never read): see plans/linda.md for how
// each step below was settled from outputs.
//
// The table is read by feature (da_table_init(), norm 'none'). Per feature:
//
//   1. Winsorization of high values at the type-7 quantile 1 - outlier_pct
//      of the feature's values: of its proportions y / N (N: the sample
//      totals) for counts, capped and turned back into counts as
//      round(q N); of the values as given for proportions.
//   2. Zeros: a pseudocount added to every value; imputation, a zero in
//      sample s becoming M_s / max M_k over the samples k where the feature
//      is zero (M: the winsorized totals); or, for proportions, half the
//      feature's smallest non-zero value.
//   3. log2, centred per sample over all features (the clr in log2).
//   4. The least squares fit to the design (feature_lm.c).
//
// The per-sample means of step 3 are sums over features: each block of
// LD_BLOCK features adds its own partial sums, and the blocks are added in
// order, so results are identical for any number of threads.

#include "ecodive.h"

#define LD_BLOCK 64

#define LD_ZERO_NONE    0
#define LD_ZERO_PSEUDO  1
#define LD_ZERO_IMPUTE  2
#define LD_ZERO_HALFMIN 3


static int cmp_double (const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}


// The type-7 quantile of x[0..n-1] (sorted into `buf`) at probability p,
// (1 - g) x_j + g x_{j+1} with j + g = 1 + (n - 1) p.
static double quantile7 (const double *x, int n, double p, double *buf) {
  memcpy(buf, x, (size_t) n * sizeof(double));
  qsort(buf, (size_t) n, sizeof(double), cmp_double);
  double idx = 1 + (n - 1) * p;
  double lo  = floor(idx), hi = ceil(idx);
  double qs  = buf[(int) lo - 1];
  double g   = idx - lo;
  if (idx > lo && buf[(int) hi - 1] != qs) qs = (1 - g) * qs + g * buf[(int) hi - 1];
  return qs;
}



//======================================================
// The shared state of both passes.
//======================================================
typedef struct {
  const da_table_t *tab;
  int     is_count;     // 1: counts (winsorize proportions); 0: as given.
  int     winsor;       // 0: no winsorization.
  double  prob;         // 1 - outlier_pct.
  const double *depth;  // n: the sample totals (counts).
  const double *depth_w;// n: the winsorized totals (imputation), or NULL.
  int     zero_mode;
  double  pc;
} ld_spec_t;


// Feature f's values after steps 1 and 2 (`y`, n), with `buf` (n) scratch.
// Returns, through *min_nz, the smallest non-zero winsorized value (Inf if
// none).
static void ld_column (const ld_spec_t *sp, int f, double *y, double *buf, double *min_nz) {

  const da_table_t *tab = sp->tab;
  int n = tab->n;

  da_table_column(tab, f, y);

  if (sp->winsor) {
    if (sp->is_count) {
      for (int s = 0; s < n; s++) buf[n + s] = y[s] / sp->depth[s];
      double q = quantile7(buf + n, n, sp->prob, buf);
      for (int s = 0; s < n; s++)
        if (buf[n + s] > q) y[s] = nearbyint(q * sp->depth[s]);
    }
    else {
      double q = quantile7(y, n, sp->prob, buf);
      for (int s = 0; s < n; s++) if (y[s] > q) y[s] = q;
    }
  }

  double mn = R_PosInf;
  for (int s = 0; s < n; s++) if (y[s] != 0 && y[s] < mn) mn = y[s];
  if (min_nz) *min_nz = mn;

  switch (sp->zero_mode) {

    case LD_ZERO_PSEUDO:
      for (int s = 0; s < n; s++) y[s] += sp->pc;
      break;

    case LD_ZERO_IMPUTE: {
      double mx = 0;
      for (int s = 0; s < n; s++) if (y[s] == 0 && sp->depth_w[s] > mx) mx = sp->depth_w[s];
      for (int s = 0; s < n; s++) if (y[s] == 0) y[s] = sp->depth_w[s] / mx;
      break;
    }

    case LD_ZERO_HALFMIN:
      for (int s = 0; s < n; s++) if (y[s] == 0) y[s] = mn / 2;
      break;
  }
}



//======================================================
// Pass 0: the winsorized totals, and each feature's
// smallest non-zero winsorized value.
//======================================================
static const ld_spec_t *ld_sp;
static double          *ld_buf;       // per thread: 4 n + 2 P doubles
static double          *ld_part;      // per block: n partial sums
static double          *ld_minnz;     // D
static int             *ld_nzero;     // D


static void *ld_prep_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  const ld_spec_t *sp = ld_sp;
  int n = sp->tab->n, D = sp->tab->D;
  int nb = (D + LD_BLOCK - 1) / LD_BLOCK;

  double *y   = ld_buf + (size_t) thread_i * 3 * n;
  double *buf = y + n;

  for (int b = thread_i; b < nb; b += n_threads) {
    double *part = ld_part + (size_t) b * n;
    for (int s = 0; s < n; s++) part[s] = 0;
    int f1 = (b + 1) * LD_BLOCK; if (f1 > D) f1 = D;
    for (int f = b * LD_BLOCK; f < f1; f++) {
      ld_column(sp, f, y, buf, ld_minnz + f);
      int nz = 0;
      for (int s = 0; s < n; s++) { part[s] += y[s]; nz += y[s] == 0; }
      ld_nzero[f] = nz;
    }
  }
  return NULL;
}


// Returns list(depth (the totals), depth_w (winsorized totals), min_nz (per
// feature, Inf if none), n_zero (per feature, after winsorization)). Totals of counts are exact (whole numbers).
SEXP C_linda_prep(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_is_count, SEXP sexp_outlier_pct,
    SEXP sexp_n_threads ) {

  int    is_count = asLogical(sexp_is_count);
  double pct      = asReal(sexp_outlier_pct);
  int    n_thr    = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, 0, 0, n_thr, &tab);
  int n = tab.n, D = tab.D;

  const char *nms[] = { "depth", "depth_w", "min_nz", "n_zero", "" };
  SEXP result  = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_depth = allocVector(REALSXP, n); SET_VECTOR_ELT(result, 0, s_depth);
  SEXP s_dw    = allocVector(REALSXP, n); SET_VECTOR_ELT(result, 1, s_dw);
  SEXP s_mn    = allocVector(REALSXP, D); SET_VECTOR_ELT(result, 2, s_mn);
  SEXP s_nz    = allocVector(INTSXP, D);  SET_VECTOR_ELT(result, 3, s_nz);

  double *depth = REAL(s_depth);
  for (int s = 0; s < n; s++) depth[s] = 0;
  for (int f = 0; f < D; f++)
    for (int e = tab.fpos[f]; e < tab.fpos[f + 1]; e++) depth[tab.fsam[e]] += tab.fval[e];

  ld_spec_t sp = {
    &tab, is_count, pct > 0, 1.0 - pct, depth, NULL, LD_ZERO_NONE, 0 };

  int nb = (D + LD_BLOCK - 1) / LD_BLOCK;
  int nt = threads_for(n_thr, (double) D * n * 8, 1e6);
  if (nt > nb) nt = nb;

  ld_sp    = &sp;
  ld_buf   = (double *) R_alloc((size_t) nt * 3 * n, sizeof(double));
  ld_part  = (double *) R_alloc((size_t) nb * n, sizeof(double));
  ld_minnz = REAL(s_mn);
  ld_nzero = INTEGER(s_nz);

  run_parallel(ld_prep_worker, nt, nb < 100 ? 100 : nb);

  double *dw = REAL(s_dw);
  for (int s = 0; s < n; s++) {
    long double t = 0;
    for (int b = 0; b < nb; b++) t += ld_part[(size_t) b * n + s];
    dw[s] = (double) t;
  }

  UNPROTECT(1);
  return result;
}



//======================================================
// Passes 1 and 2: the per-sample means of the log2
// values, then each feature's fit.
//======================================================
static int           ld_P;
static const flm_t  *ld_full;
static double       *ld_mem;          // per thread: flm workspace
static size_t        ld_mem_sz;
static const double *ld_center;       // n: per-sample means
static double       *ld_coef, *ld_se, *ld_rss, *ld_mean;
static const lmm_design_t *ld_lmm;    // mixed models: the design, or NULL
static double       *ld_lmem;         // per thread: lmm workspace
static size_t        ld_lmem_sz;
static double       *ld_df, *ld_theta, *ld_vcov;


static void *ld_sum_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  const ld_spec_t *sp = ld_sp;
  int n = sp->tab->n, D = sp->tab->D, P = ld_P;
  int nb = (D + LD_BLOCK - 1) / LD_BLOCK;

  double *y   = ld_buf + (size_t) thread_i * (4 * (size_t) n + 2 * (size_t) P);
  double *buf = y + n;

  for (int b = thread_i; b < nb; b += n_threads) {
    double *part = ld_part + (size_t) b * n;
    for (int s = 0; s < n; s++) part[s] = 0;
    int f1 = (b + 1) * LD_BLOCK; if (f1 > D) f1 = D;
    for (int f = b * LD_BLOCK; f < f1; f++) {
      ld_column(sp, f, y, buf, NULL);
      for (int s = 0; s < n; s++) part[s] += log2(y[s]);
    }
  }
  return NULL;
}


static void *ld_fit_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  const ld_spec_t *sp = ld_sp;
  int n = sp->tab->n, D = sp->tab->D, P = ld_P;

  double *y     = ld_buf + (size_t) thread_i * (4 * (size_t) n + 2 * (size_t) P);
  double *buf   = y + n;          // 2 n (the proportions, for winsorizing)
  double *qty   = buf + 2 * (size_t) n;
  double *coef  = qty + n;
  double *se    = coef + P;

  // The shared QR, with this thread's own scratch.
  flm_t q = *ld_full;
  q.work  = ld_mem + (size_t) thread_i * ld_mem_sz;

  for (int f = thread_i; f < D; f += n_threads) {

    ld_column(sp, f, y, buf, NULL);

    long double m = 0;
    for (int s = 0; s < n; s++) { y[s] = log2(y[s]) - ld_center[s]; m += y[s]; }
    ld_mean[f] = (double) (m / n);

    if (ld_lmm) {
      double *lm  = ld_lmem + (size_t) thread_i * ld_lmem_sz;
      double *vc  = lm + lmm_doubles(n, P, ld_lmm->J);
      double *df  = vc + (size_t) P * P;
      double sigma, dev;
      ld_theta[f] = lmm_fit(ld_lmm, y, -1, lm, coef, vc, df, &sigma, &dev);
      for (int j = 0; j < P; j++) {
        ld_coef[f + (size_t) D * j] = coef[j];
        ld_se[f + (size_t) D * j]   = sqrt(vc[j + (size_t) P * j]);
        ld_df[f + (size_t) D * j]   = df[j];
        for (int l = 0; l < P; l++)
          ld_vcov[f + (size_t) D * (j + (size_t) P * l)] = vc[j + (size_t) P * l];
      }
      ld_rss[f] = sigma;
      continue;
    }

    double rss;
    flm_fit(&q, y, qty, coef, &rss);
    flm_se(&q, y, qty, rss, 0, se, NULL);

    for (int j = 0; j < P; j++) {
      ld_coef[f + (size_t) D * j] = coef[j];
      ld_se[f + (size_t) D * j]   = se[j];
    }
    ld_rss[f] = rss;
  }

  return NULL;
}


// Returns list(coef, se (D x P), rss, mean (each feature's mean centred
// log2 value), rank (of the design)). With `groups` (1-based group codes,
// or NULL), each feature gets a random-intercept model (feature_lmm.c)
// instead, and the list adds df (D x P, Satterthwaite's), theta and vcov
// (D x P x P); rss is then the residual standard deviation.
SEXP C_linda_fit(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_is_count, SEXP sexp_outlier_pct,
    SEXP sexp_depth, SEXP sexp_depth_w, SEXP sexp_zero_mode, SEXP sexp_pc,
    SEXP sexp_X, SEXP sexp_tol, SEXP sexp_groups, SEXP sexp_n_threads ) {

  int    is_count = asLogical(sexp_is_count);
  double pct      = asReal(sexp_outlier_pct);
  int    zmode    = asInteger(sexp_zero_mode);
  double pc       = asReal(sexp_pc);
  double tol      = asReal(sexp_tol);
  int    n_thr    = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, 0, 0, n_thr, &tab);
  int n = tab.n, D = tab.D;

  SEXP dims = getAttrib(sexp_X, R_DimSymbol);
  if (INTEGER(dims)[0] != n) error("The design must have one row per sample.");  // # nocov
  int P = INTEGER(dims)[1];

  ld_spec_t sp = {
    &tab, is_count, pct > 0, 1.0 - pct, REAL(sexp_depth), REAL(sexp_depth_w), zmode, pc };

  // The design's QR, on this thread.
  size_t  msz   = flm_doubles(n, P, 0);
  double *fmem  = (double *) R_alloc(msz, sizeof(double));
  int    *fimem = (int *)    R_alloc((size_t) P + 1, sizeof(int));
  flm_t   full;
  flm_init(&full, fmem, fimem, n, P, 0);
  flm_decompose(&full, REAL(sexp_X), n, NULL, n, P, tol);

  int mixed = !isNull(sexp_groups);
  const char *nms[] = { "coef", "se", "rss", "mean", "rank", "df", "theta", "vcov", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_coef = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 0, s_coef);
  SEXP s_se   = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 1, s_se);
  SEXP s_rss  = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 2, s_rss);
  SEXP s_mean = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 3, s_mean);
  SET_VECTOR_ELT(result, 4, ScalarInteger(full.rank));

  lmm_design_t lmm;
  ld_lmm = NULL;
  if (mixed) {
    int *grp = (int *) R_alloc((size_t) n + 1, sizeof(int));
    int J = 0;
    for (int i = 0; i < n; i++) { grp[i] = INTEGER(sexp_groups)[i] - 1; if (grp[i] + 1 > J) J = grp[i] + 1; }
    lmm_design_init(&lmm, REAL(sexp_X), grp, n, P, J);
    ld_lmm = &lmm;
    SEXP s_df = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 5, s_df);
    SEXP s_th = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 6, s_th);
    SEXP s_vc = PROTECT(allocVector(REALSXP, (R_xlen_t) D * P * P));
    SEXP dim  = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dim)[0] = D; INTEGER(dim)[1] = P; INTEGER(dim)[2] = P;
    setAttrib(s_vc, R_DimSymbol, dim);
    SET_VECTOR_ELT(result, 7, s_vc);
    UNPROTECT(2);
    ld_df    = REAL(s_df);
    ld_theta = REAL(s_th);
    ld_vcov  = REAL(s_vc);
  }

  int nb = (D + LD_BLOCK - 1) / LD_BLOCK;
  int nt  = threads_for(n_thr, (double) D * n * (P + 8), 1e6);
  int nt1 = nt > nb ? nb : nt;
  // A mixed model's fit iterates: worth threads for smaller tables.
  if (mixed) nt = threads_for(n_thr, (double) D * (n * P + ld_lmm->J * P * P) * 60, 1e6);
  if (nt < nt1) nt = nt1;

  size_t per = 4 * (size_t) n + 2 * (size_t) P;
  ld_sp     = &sp;
  ld_P      = P;
  ld_full   = &full;
  ld_buf    = (double *) R_alloc((size_t) nt * per, sizeof(double));
  ld_part   = (double *) R_alloc((size_t) nb * n, sizeof(double));
  ld_mem_sz = (size_t) n + 3 * (size_t) P;
  ld_mem    = (double *) R_alloc((size_t) nt * ld_mem_sz, sizeof(double));

  // Pass 1: the per-sample means of the log2 values.
  run_parallel(ld_sum_worker, nt1, nb < 100 ? 100 : nb);

  double *center = (double *) R_alloc((size_t) n + 1, sizeof(double));
  for (int s = 0; s < n; s++) {
    long double t = 0;
    for (int b = 0; b < nb; b++) t += ld_part[(size_t) b * n + s];
    center[s] = (double) (t / D);
  }

  // Pass 2: the fits.
  ld_center = center;
  ld_coef   = REAL(s_coef);
  ld_se     = REAL(s_se);
  ld_rss    = REAL(s_rss);
  ld_mean   = REAL(s_mean);

  if (mixed) {
    ld_lmem_sz = lmm_doubles(n, P, ld_lmm->J) + (size_t) P * P + P;
    ld_lmem    = (double *) R_alloc((size_t) nt * ld_lmem_sz, sizeof(double));
  }

  run_parallel(ld_fit_worker, nt, D < 100 ? 100 : D);

  UNPROTECT(1);
  return result;
}



//======================================================
// The mode of each column of B (D x P) by the mean shift
// algorithm with a Gaussian kernel: from start[j], with
// bandwidth h[j], x <- sum w b / sum w, w = exp(-((x -
// b) / h)^2 / 2), until a step is smaller than
// tol |x| (then the point before it is the mode) or
// max_iters steps. NaN coefficients are left out.
// Returns list(mode, n_iters).
//======================================================
SEXP C_linda_mode(SEXP sexp_B, SEXP sexp_start, SEXP sexp_h, SEXP sexp_tol, SEXP sexp_max_iters) {

  int D = nrows(sexp_B), P = ncols(sexp_B);
  double tol   = asReal(sexp_tol);
  int    max_i = asInteger(sexp_max_iters);

  const char *nms[] = { "mode", "n_iters", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_mode = allocVector(REALSXP, P); SET_VECTOR_ELT(result, 0, s_mode);
  SEXP s_it   = allocVector(INTSXP, P);  SET_VECTOR_ELT(result, 1, s_it);

  for (int j = 0; j < P; j++) {
    const double *b = REAL(sexp_B) + (size_t) D * j;
    double x = REAL(sexp_start)[j], h = REAL(sexp_h)[j];
    int it = 0;
    if (ISNAN(x) || !(h > 0)) { REAL(s_mode)[j] = R_NaN; INTEGER(s_it)[j] = 0; continue; }
    for (it = 1; it <= max_i; it++) {
      long double sw = 0, swb = 0;
      for (int k = 0; k < D; k++) {
        if (ISNAN(b[k])) continue;
        double z = (x - b[k]) / h;
        double w = exp(-0.5 * z * z);
        sw  += w;
        double wb = w * b[k];
        swb += wb;
      }
      double xn = (double) (swb / sw);
      if (fabs(xn - x) < tol * fabs(x)) break;
      x = xn;
    }
    REAL(s_mode)[j]   = x;
    INTEGER(s_it)[j]  = it > max_i ? max_i : it;
  }

  UNPROTECT(1);
  return result;
}
