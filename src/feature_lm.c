// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Least squares fits of one design to many responses (features), for
// da_lm(), and planned for linda(), fastancom() and dirmult_lm(): a pivoted
// Householder QR of the design, once, and per feature the coefficients, the
// residual sum of squares, and classical, HC0 or HC3 standard errors.
// Features with missing values (NaN) get their own QR of the design's rows
// that they were observed in (a "masked" design), as lm() with na.omit.
//
// The QR follows the algorithm of LINPACK's dqrdc2 with R's limited column
// pivoting, as described in ?qr: a column whose norm, after the
// reflections so far, is below `tol` (1e-7, as lm()) times its original
// norm is moved to the end and left out of the fit (aliased: lm()'s NA
// coefficients). Written from that description; no code is copied.
//
// No BLAS: plain loops in a fixed order, so that worker threads can call
// every function here, and results are the same for any number of them.

#include "ecodive.h"


//======================================================
// Workspace: doubles (and p ints) for an n x p design.
//======================================================
size_t flm_doubles(int n, int p, int want_q1) {
  size_t np = (size_t) n * p;
  size_t sz = np + (size_t) p + (size_t) p * p + (size_t) n + 3 * (size_t) p;
  if (want_q1) sz += np + (size_t) n;
  return sz;
}


void flm_init(flm_t *q, double *mem, int *imem, int n, int p, int want_q1) {
  size_t np = (size_t) n * p;
  q->n     = n;
  q->p     = p;
  q->rank  = 0;
  q->a     = mem;           mem += np;
  q->qraux = mem;           mem += p;
  q->rinv  = mem;           mem += (size_t) p * p;
  q->work  = mem;           mem += (size_t) n + 3 * (size_t) p;
  q->q1    = NULL;
  q->h     = NULL;
  if (want_q1) {
    q->q1 = mem;            mem += np;
    q->h  = mem;
  }
  q->pivot = imem;
}



//======================================================
// The QR of rows `rows` (n of them; NULL: the first n)
// of the design X (ldx rows, p columns, column-major).
//======================================================
void flm_decompose(
    flm_t *q, const double *X, int ldx, const int *rows, int n, int p, double tol ) {

  q->n = n;
  q->p = p;
  double *a     = q->a;
  double *qraux = q->qraux;
  int    *piv   = q->pivot;
  double *ref   = q->work + n;       // p: each column's original norm (or 1)

  for (int j = 0; j < p; j++) {
    const double *x = X + (size_t) ldx * j;
    double       *c = a + (size_t) n * j;
    if (rows) for (int i = 0; i < n; i++) c[i] = x[rows[i]];
    else      memcpy(c, x, (size_t) n * sizeof(double));
    long double ss = 0;
    for (int i = 0; i < n; i++) ss += (long double) c[i] * c[i];
    ref[j]   = sqrtl(ss);
    if (ref[j] == 0) ref[j] = 1;
    piv[j]   = j;
    qraux[j] = 0;
  }

  int lup  = p;
  int lmax = n < p ? n : p;
  int l    = 0;

  for (; l < lmax && l < lup; l++) {

    // Limited pivoting: move near-collinear columns to the end.
    for (;;) {
      double *c = a + (size_t) n * l;
      long double ss = 0;
      for (int i = l; i < n; i++) ss += (long double) c[i] * c[i];
      if ((double) sqrtl(ss) >= ref[l] * tol) break;

      double *tmp = q->work;      // n scratch
      memcpy(tmp, c, (size_t) n * sizeof(double));
      int    pv = piv[l];
      double rf = ref[l];
      for (int j = l; j < p - 1; j++) {
        memcpy(a + (size_t) n * j, a + (size_t) n * (j + 1), (size_t) n * sizeof(double));
        piv[j] = piv[j + 1];
        ref[j] = ref[j + 1];
      }
      memcpy(a + (size_t) n * (p - 1), tmp, (size_t) n * sizeof(double));
      piv[p - 1] = pv;
      ref[p - 1] = rf;

      if (--lup <= l) break;
    }
    if (l >= lup) break;

    // No reflection for the last row.
    if (l == n - 1) { qraux[l] = 0; continue; }

    double *c = a + (size_t) n * l;
    long double ss = 0;
    for (int i = l; i < n; i++) ss += (long double) c[i] * c[i];
    double nrmxl = (double) sqrtl(ss);
    if (nrmxl == 0) { qraux[l] = 0; continue; }
    if (c[l] != 0) nrmxl = copysign(nrmxl, c[l]);

    for (int i = l; i < n; i++) c[i] /= nrmxl;
    c[l] += 1;

    for (int j = l + 1; j < p; j++) {
      double *cj = a + (size_t) n * j;
      long double d = 0;
      for (int i = l; i < n; i++) d += (long double) c[i] * cj[i];
      double t = (double) (-d / c[l]);
      for (int i = l; i < n; i++) cj[i] += t * c[i];
    }

    qraux[l] = c[l];
    c[l]     = -nrmxl;
  }

  int rank = lup < n ? lup : n;
  q->rank  = rank;


  //--------------------------------------------------------
  // R's inverse (upper triangular, rank x rank, in p x p).
  //--------------------------------------------------------
  double *ri = q->rinv;
  for (int j = 0; j < rank; j++) {
    for (int i = 0; i < p; i++) ri[i + (size_t) p * j] = 0;
    ri[j + (size_t) p * j] = 1 / a[j + (size_t) n * j];
    for (int i = j - 1; i >= 0; i--) {
      long double s = 0;
      for (int k = i + 1; k <= j; k++)
        s += (long double) a[i + (size_t) n * k] * ri[k + (size_t) p * j];
      ri[i + (size_t) p * j] = (double) (-s / a[i + (size_t) n * i]);
    }
  }


  //--------------------------------------------------------
  // Q's first rank columns and the leverages, for HC
  // standard errors.
  //--------------------------------------------------------
  if (q->q1) {
    for (int k = 0; k < rank; k++) {
      double *v = q->q1 + (size_t) n * k;
      for (int i = 0; i < n; i++) v[i] = 0;
      v[k] = 1;
      for (int r = rank - 1; r >= 0; r--) {
        if (qraux[r] == 0) continue;
        const double *u = a + (size_t) n * r;
        long double d = (long double) qraux[r] * v[r];
        for (int i = r + 1; i < n; i++) d += (long double) u[i] * v[i];
        double t = (double) (-d / qraux[r]);
        v[r] += t * qraux[r];
        for (int i = r + 1; i < n; i++) v[i] += t * u[i];
      }
    }
    for (int i = 0; i < n; i++) {
      long double s = 0;
      for (int k = 0; k < rank; k++) { double v = q->q1[i + (size_t) n * k]; s += (long double) v * v; }
      q->h[i] = (double) s;
    }
  }
}



//======================================================
// Q'y (n), the coefficients (p, in the design's column
// order, NaN where aliased), and the residual sum of
// squares, for a response y of the decomposed rows.
//======================================================
void flm_fit(const flm_t *q, const double *y, double *qty, double *coef, double *rss) {

  int n = q->n, p = q->p, rank = q->rank;
  const double *a = q->a;

  memcpy(qty, y, (size_t) n * sizeof(double));
  for (int l = 0; l < rank; l++) {
    if (q->qraux[l] == 0) continue;
    const double *u = a + (size_t) n * l;
    long double d = (long double) q->qraux[l] * qty[l];
    for (int i = l + 1; i < n; i++) d += (long double) u[i] * qty[i];
    double t = (double) (-d / q->qraux[l]);
    qty[l] += t * q->qraux[l];
    for (int i = l + 1; i < n; i++) qty[i] += t * u[i];
  }

  // Back substitution: R b = (Q'y)[1:rank].
  double *b = q->work + n;   // p scratch (after flm_decompose's use)
  for (int i = rank - 1; i >= 0; i--) {
    long double s = qty[i];
    for (int k = i + 1; k < rank; k++) s -= (long double) a[i + (size_t) n * k] * b[k];
    b[i] = (double) (s / a[i + (size_t) n * i]);
  }

  for (int j = 0; j < p; j++) coef[j] = R_NaN;
  for (int j = 0; j < rank; j++) coef[q->pivot[j]] = b[j];

  long double ss = 0;
  for (int i = rank; i < n; i++) ss += (long double) qty[i] * qty[i];
  *rss = (double) ss;
}



//======================================================
// The residuals (n): Q (0, (Q'y)[rank + 1:n]), from flm_fit()'s
// Q'y. Uses no scratch, so worker threads can share `q`.
//======================================================
void flm_resid(const flm_t *q, const double *qty, double *resid) {

  int n = q->n, rank = q->rank;

  for (int i = 0; i < n; i++) resid[i] = i < rank ? 0 : qty[i];
  for (int r = rank - 1; r >= 0; r--) {
    if (q->qraux[r] == 0) continue;
    const double *u = q->a + (size_t) n * r;
    long double d = (long double) q->qraux[r] * resid[r];
    for (int i = r + 1; i < n; i++) d += (long double) u[i] * resid[i];
    double t = (double) (-d / q->qraux[r]);
    resid[r] += t * q->qraux[r];
    for (int i = r + 1; i < n; i++) resid[i] += t * u[i];
  }
}



//======================================================
// Standard errors (p, NaN where aliased or without
// residual degrees of freedom): se_type 0 classical,
// 1 HC0, 2 HC3 (these need q1). `resid` (n) gets the
// residuals for HC0 / HC3.
//======================================================
void flm_se(
    const flm_t *q, const double *y, const double *qty, double rss,
    int se_type, double *se, double *resid ) {

  (void) y;
  int n = q->n, p = q->p, rank = q->rank;
  const double *ri = q->rinv;

  for (int j = 0; j < p; j++) se[j] = R_NaN;
  if (n - rank <= 0 && se_type == 0) return;

  if (se_type == 0) {
    double s2 = rss / (n - rank);
    for (int k = 0; k < rank; k++) {
      long double s = 0;
      for (int m = k; m < rank; m++) { double v = ri[k + (size_t) p * m]; s += (long double) v * v; }
      se[q->pivot[k]] = sqrt(s2 * (double) s);
    }
    return;
  }

  flm_resid(q, qty, resid);

  // Var(b_k) = sum_j A_kj^2 w_j, A = R^-1 Q1'.
  for (int k = 0; k < rank; k++) {
    long double v = 0;
    for (int j = 0; j < n; j++) {
      long double akj = 0;
      for (int m = k; m < rank; m++)
        akj += (long double) ri[k + (size_t) p * m] * q->q1[j + (size_t) n * m];
      double e = resid[j];
      double w = e * e;
      if (se_type == 2) { double d = 1 - q->h[j]; w /= d * d; }
      v += akj * akj * w;
    }
    se[q->pivot[k]] = sqrt((double) v);
  }
}



//======================================================
// da_lm(): every feature's fit.
//======================================================
static const da_table_t *lm_tab;
static const double     *lm_X;
static int               lm_P;
static double            lm_tol;
static int               lm_se_type;
static const flm_t      *lm_full;     // the complete design's QR
static double           *lm_mem;      // per thread: flm workspace
static int              *lm_imem;     // per thread: p + n ints
static size_t            lm_mem_sz;
static double           *lm_buf;      // per thread: 4 n + 2 p doubles
static double           *lm_coef, *lm_se, *lm_rss, *lm_dof, *lm_mean, *lm_nobs;


static void *lm_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  const da_table_t *tab = lm_tab;
  int n = tab->n, D = tab->D, P = lm_P;

  double *y     = lm_buf + (size_t) thread_i * (4 * (size_t) n + 2 * (size_t) P);
  double *yo    = y  + n;
  double *qty   = yo + n;
  double *resid = qty + n;
  double *coef  = resid + n;
  double *se    = coef + P;
  int    *rows  = lm_imem + (size_t) thread_i * ((size_t) P + n) + P;

  flm_t mq;
  flm_init(&mq, lm_mem + (size_t) thread_i * lm_mem_sz,
           lm_imem + (size_t) thread_i * ((size_t) P + n), n, P, lm_se_type > 0);

  for (int f = thread_i; f < D; f += n_threads) {

    da_table_column(tab, f, y);

    int n_obs = 0;
    for (int s = 0; s < n; s++) if (!ISNAN(y[s])) { rows[n_obs] = s; yo[n_obs++] = y[s]; }

    double rss = R_NaN, dof = R_NaN;
    const flm_t *q = NULL;

    if (n_obs == n) {
      q = lm_full;
    }
    else if (n_obs > 0) {
      flm_decompose(&mq, lm_X, n, rows, n_obs, P, lm_tol);
      q = &mq;
    }

    if (q) {
      // lm_full is shared: its scratch isn't used by flm_fit()
      // beyond `work`, so fit with a copy of the struct whose
      // scratch is this thread's.
      flm_t qq = *q;
      qq.work  = mq.work;
      flm_fit(&qq, yo, qty, coef, &rss);
      flm_se(&qq, yo, qty, rss, lm_se_type, se, resid);
      dof = n_obs - qq.rank;
    }
    else {
      for (int j = 0; j < P; j++) coef[j] = se[j] = R_NaN;
    }

    for (int j = 0; j < P; j++) {
      lm_coef[f + (size_t) D * j] = coef[j];
      lm_se[f + (size_t) D * j]   = se[j];
    }
    lm_rss[f]  = rss;
    lm_dof[f]  = dof;
    lm_nobs[f] = n_obs;

    long double s = 0;
    for (int i = 0; i < n_obs; i++) s += yo[i];
    if (n_obs) {
      s /= n_obs;
      if (R_FINITE((double) s)) {
        long double t = 0;
        for (int i = 0; i < n_obs; i++) t += (yo[i] - s);
        s += t / n_obs;
      }
    }
    lm_mean[f] = n_obs ? (double) s : R_NaN;
  }

  return NULL;
}



//======================================================
// R interface, for da_lm(). Returns list(coef, se (D x P),
// rss, dof (residual degrees of freedom), n_obs, mean
// (of the observed normalized values), rank (of the
// complete design)).
//======================================================
SEXP C_da_lm(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_pseudocount,
    SEXP sexp_X, SEXP sexp_tol, SEXP sexp_se_type, SEXP sexp_n_threads ) {

  int    norm    = asInteger(sexp_norm);
  double pc      = asReal(sexp_pseudocount);
  double tol     = asReal(sexp_tol);
  int    se_type = asInteger(sexp_se_type);
  int    n_thr   = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, norm, pc, n_thr, &tab);
  int n = tab.n, D = tab.D;

  SEXP dims = getAttrib(sexp_X, R_DimSymbol);
  if (INTEGER(dims)[0] != n) error("The design must have one row per sample.");  // # nocov
  int P = INTEGER(dims)[1];

  // The complete design's QR, on this thread.
  size_t  msz  = flm_doubles(n, P, se_type > 0);
  double *fmem = (double *) R_alloc(msz, sizeof(double));
  int    *fimem = (int *)   R_alloc((size_t) P + 1, sizeof(int));
  flm_t   full;
  flm_init(&full, fmem, fimem, n, P, se_type > 0);
  flm_decompose(&full, REAL(sexp_X), n, NULL, n, P, tol);

  int nt = threads_for(n_thr, (double) D * n * (P + 2) * 4, 1e6);

  const char *nms[] = { "coef", "se", "rss", "dof", "n_obs", "mean", "rank", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_coef = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 0, s_coef);
  SEXP s_se   = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 1, s_se);
  SEXP s_rss  = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 2, s_rss);
  SEXP s_dof  = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 3, s_dof);
  SEXP s_nobs = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 4, s_nobs);
  SEXP s_mean = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 5, s_mean);
  SET_VECTOR_ELT(result, 6, ScalarInteger(full.rank));

  lm_tab     = &tab;
  lm_X       = REAL(sexp_X);
  lm_P       = P;
  lm_tol     = tol;
  lm_se_type = se_type;
  lm_full    = &full;
  lm_mem_sz  = msz;
  lm_mem     = (double *) R_alloc((size_t) nt * msz, sizeof(double));
  lm_imem    = (int *)    R_alloc((size_t) nt * ((size_t) P + n), sizeof(int));
  lm_buf     = (double *) R_alloc((size_t) nt * (4 * (size_t) n + 2 * (size_t) P), sizeof(double));
  lm_coef    = REAL(s_coef);
  lm_se      = REAL(s_se);
  lm_rss     = REAL(s_rss);
  lm_dof     = REAL(s_dof);
  lm_nobs    = REAL(s_nobs);
  lm_mean    = REAL(s_mean);

  run_parallel(lm_worker, nt, D < 100 ? 100 : D);

  UNPROTECT(1);
  return result;
}



//======================================================
// R interface, for the tests of the kernel: the fit of
// each column of Y (n x m, NaN for missing) to X (n x p).
// Returns list(coef, se (m x p), rss, rank, resid (n x m,
// NaN where missing; HC only)).
//======================================================
SEXP C_feature_lm(SEXP sexp_X, SEXP sexp_Y, SEXP sexp_tol, SEXP sexp_se_type) {

  int n = nrows(sexp_X), p = ncols(sexp_X), m = ncols(sexp_Y);
  double tol     = asReal(sexp_tol);
  int    se_type = asInteger(sexp_se_type);

  const char *nms[] = { "coef", "se", "rss", "rank", "resid", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_coef  = allocMatrix(REALSXP, m, p); SET_VECTOR_ELT(result, 0, s_coef);
  SEXP s_se    = allocMatrix(REALSXP, m, p); SET_VECTOR_ELT(result, 1, s_se);
  SEXP s_rss   = allocVector(REALSXP, m);    SET_VECTOR_ELT(result, 2, s_rss);
  SEXP s_rank  = allocVector(INTSXP, m);     SET_VECTOR_ELT(result, 3, s_rank);
  SEXP s_resid = allocMatrix(REALSXP, n, m); SET_VECTOR_ELT(result, 4, s_resid);

  double *mem  = (double *) R_alloc(flm_doubles(n, p, se_type > 0), sizeof(double));
  int    *imem = (int *)    R_alloc((size_t) p + 1, sizeof(int));
  int    *rows = (int *)    R_alloc((size_t) n + 1, sizeof(int));
  double *yo   = (double *) R_alloc((size_t) n + 1, sizeof(double));
  double *qty  = (double *) R_alloc((size_t) n + 1, sizeof(double));
  double *res  = (double *) R_alloc((size_t) n + 1, sizeof(double));
  double *coef = (double *) R_alloc((size_t) p + 1, sizeof(double));
  double *se   = (double *) R_alloc((size_t) p + 1, sizeof(double));

  flm_t q;
  flm_init(&q, mem, imem, n, p, se_type > 0);

  for (int c = 0; c < m; c++) {
    const double *y = REAL(sexp_Y) + (size_t) n * c;
    int n_obs = 0;
    for (int i = 0; i < n; i++) if (!ISNAN(y[i])) { rows[n_obs] = i; yo[n_obs++] = y[i]; }

    flm_decompose(&q, REAL(sexp_X), n, rows, n_obs, p, tol);
    double rss;
    flm_fit(&q, yo, qty, coef, &rss);
    flm_se(&q, yo, qty, rss, se_type, se, res);

    for (int j = 0; j < p; j++) {
      REAL(s_coef)[c + (size_t) m * j] = coef[j];
      REAL(s_se)[c + (size_t) m * j]   = se[j];
    }
    REAL(s_rss)[c]     = rss;
    INTEGER(s_rank)[c] = q.rank;
    double *r = REAL(s_resid) + (size_t) n * c;
    for (int i = 0; i < n; i++) r[i] = R_NaN;
    if (se_type > 0) for (int i = 0; i < n_obs; i++) r[rows[i]] = res[i];
  }

  UNPROTECT(1);
  return result;
}
