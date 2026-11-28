// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// fastancom(): ANCOM from K marginal models (Zhou et al. 2022), written
// from the paper's abstract, fastANCOM's documentation and its outputs
// (fastANCOM is GPL-3; its code was not read). How each step was settled
// from outputs is recorded in plans/fastancom.md.
//
// With one design X for every feature, the least squares fit of a log
// ratio L_k - L_j is the difference of the features' own fits: its
// coefficients are b_k - b_j and its residuals e_k - e_j. So every one of
// ANCOM's K (K - 1) / 2 log-ratio tests follows from K fits:
//
//   F_kj = d' A d / (q s2),   d = Lambda' (b_k - b_j),
//                             s2 = |e_k - e_j|^2 / (n - rank),
//
// for a comparison Lambda (P x q columns of contrasts) with A = (Lambda'
// (X'X)^-1 Lambda)^-1; for q = 1, F = t^2. W_k counts the j with F_kj above
// the critical value (p < alpha). The effects compare each feature with the
// mean of a reference set R (the features with the smallest W):
// d = Lambda' (b_k - mean_R b), s2 = |e_k - mean_R e|^2 / (n - rank).
//
// Phases, each split into tasks that don't depend on the thread count:
//   1. Per feature (dealt in turn): log2(x + pc), its fit, residuals and
//      contrasts.
//   2. Per row k of the triangle of pairs (dealt in turn): the tests of k
//      against every j > k. Each thread counts into its own integers, which
//      are summed afterwards (exact in any order).
//   3. Per comparison, on this thread: the reference set (the nref smallest
//      W, ties in the features' order) and its mean residuals; then per
//      feature (dealt in turn): the effect, its standard error and statistic.

#include "ecodive.h"


//======================================================
// Shared state.
//======================================================
static const da_table_t *fa_tab;
static const flm_t      *fa_qr;      // The design's QR (shared, read-only).
static double           *fa_work;    // Per thread: n + 3 P doubles (flm scratch).
static double           *fa_buf;     // Per thread: 3 n + P doubles.
static int               fa_n, fa_D, fa_P, fa_icol, fa_m, fa_qsum;
static double            fa_pc, fa_df;

static const int        *fa_q;       // m: each comparison's columns.
static const int        *fa_qoff;    // m: its offset among the qsum.
static const double    **fa_L;       // m: P x q contrasts (column-major).
static const double    **fa_A;       // m: q x q inverse covariances.
static const double     *fa_crit;    // m: F critical values.

static double *fa_coef;    // D x P (column-major, D rows).
static double *fa_resid;   // n per feature (feature-major).
static double *fa_b;       // qsum per feature: each comparison's contrasts.
static double *fa_mean;    // D: each feature's mean log2 value.
static int    *fa_Wt;      // Per thread: D x m counts.

// Phase 3's comparison.
static int     fa_c;
static double *fa_ebar;    // n: the reference set's mean residuals.
static double *fa_bbar;    // q: its mean contrasts.
static double *fa_est, *fa_se, *fa_stat;   // D x m (column-major).



//======================================================
// d' A d / q, for d (q) and A (q x q, column-major).
//======================================================
static double fa_quad(const double *d, const double *A, int q) {
  if (q == 1) return d[0] * d[0] * A[0];
  double s = 0;
  for (int j = 0; j < q; j++) {
    double r = 0;
    for (int i = 0; i < q; i++) r += A[i + (size_t) q * j] * d[i];
    s += d[j] * r;
  }
  return s / q;
}



//======================================================
// Phase 1: every feature's fit.
//======================================================
static void *fa_fit_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  int n = fa_n, D = fa_D, P = fa_P;

  double *y    = fa_buf + (size_t) thread_i * (3 * (size_t) n + P);
  double *qty  = y + n;
  double *e    = qty + n;
  double *coef = e + n;

  flm_t qq = *fa_qr;
  qq.work  = fa_work + (size_t) thread_i * ((size_t) n + 3 * (size_t) P);

  for (int f = thread_i; f < D; f += n_threads) {

    da_table_column(fa_tab, f, y);
    for (int s = 0; s < n; s++) y[s] = log2(y[s] + fa_pc);

    int constant = 1;
    for (int s = 1; s < n && constant; s++) if (y[s] != y[0]) constant = 0;

    if (constant) {
      // Exactly: the intercept alone, with no residuals.
      for (int j = 0; j < P; j++) coef[j] = 0;
      coef[fa_icol] = y[0];
      for (int s = 0; s < n; s++) e[s] = 0;
    }
    else {
      double rss;
      flm_fit(&qq, y, qty, coef, &rss);
      flm_resid(&qq, qty, e);
    }

    for (int j = 0; j < P; j++) fa_coef[f + (size_t) D * j] = coef[j];
    memcpy(fa_resid + (size_t) n * f, e, (size_t) n * sizeof(double));

    for (int c = 0; c < fa_m; c++) {
      const double *L = fa_L[c];
      double       *b = fa_b + (size_t) fa_qsum * f + fa_qoff[c];
      for (int i = 0; i < fa_q[c]; i++) {
        double s = 0;
        for (int j = 0; j < P; j++) s += L[j + (size_t) P * i] * coef[j];
        b[i] = s;
      }
    }

    long double s = 0;
    for (int i = 0; i < n; i++) s += y[i];
    fa_mean[f] = (double) (s / n);
  }

  return NULL;
}



//======================================================
// Phase 2: the log-ratio tests, row k against j > k.
//======================================================
static void *fa_pair_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  int n = fa_n, D = fa_D, m = fa_m, qsum = fa_qsum;
  int    *W = fa_Wt + (size_t) thread_i * D * m;
  double *d = fa_buf + (size_t) thread_i * (3 * (size_t) n + fa_P);

  for (int k = thread_i; k < D - 1; k += n_threads) {

    const double *ek = fa_resid + (size_t) n * k;
    const double *bk = fa_b + (size_t) qsum * k;

    for (int j = k + 1; j < D; j++) {

      const double *ej = fa_resid + (size_t) n * j;
      const double *bj = fa_b + (size_t) qsum * j;

      double rss = 0;
      for (int s = 0; s < n; s++) { double r = ek[s] - ej[s]; rss += r * r; }
      double s2 = rss / fa_df;

      for (int c = 0; c < m; c++) {
        int q = fa_q[c], o = fa_qoff[c];
        for (int i = 0; i < q; i++) d[i] = bk[o + i] - bj[o + i];
        double num = fa_quad(d, fa_A[c], q);
        if (num == 0 && s2 == 0) continue;        // No test: NaN.
        if (num / s2 > fa_crit[c]) {               // Inf when s2 == 0.
          W[k + (size_t) D * c]++;
          W[j + (size_t) D * c]++;
        }
      }
    }
  }

  return NULL;
}



//======================================================
// Phase 3: each feature against the reference set.
//======================================================
static void *fa_effect_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  int n = fa_n, D = fa_D, c = fa_c, q = fa_q[c], o = fa_qoff[c];
  double *d = fa_buf + (size_t) thread_i * (3 * (size_t) n + fa_P);

  for (int k = thread_i; k < D; k += n_threads) {

    const double *ek = fa_resid + (size_t) n * k;
    const double *bk = fa_b + (size_t) fa_qsum * k + o;

    double rss = 0;
    for (int s = 0; s < n; s++) { double r = ek[s] - fa_ebar[s]; rss += r * r; }
    double s2 = rss / fa_df;
    for (int i = 0; i < q; i++) d[i] = bk[i] - fa_bbar[i];

    size_t at = k + (size_t) D * c;
    if (q == 1) {
      double se = sqrt(s2 / fa_A[c][0]);
      fa_est[at]  = d[0];
      fa_se[at]   = se;
      fa_stat[at] = d[0] / se;                     // NaN for 0 / 0.
    }
    else {
      double num = fa_quad(d, fa_A[c], q);
      fa_est[at]  = R_NaN;
      fa_se[at]   = R_NaN;
      fa_stat[at] = (num == 0 && s2 == 0) ? R_NaN : num / s2;
    }
  }

  return NULL;
}



//======================================================
// R interface. X (n x P, full rank) with its intercept in
// column icol (0-based); cmp a list of comparisons, each
// list(L (P x q contrasts), A (q x q), crit); df = n - P;
// nref the size of the reference sets.
//
// Returns list(coef (D x P), W (D x m integer), est, se,
// stat (D x m; est and se NaN for q > 1, stat t for q = 1
// and F for q > 1), ref (m integer vectors, 1-based, in
// the order chosen), mean (D)).
//======================================================
SEXP C_fastancom(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_pseudocount, SEXP sexp_X,
    SEXP sexp_icol, SEXP sexp_cmp, SEXP sexp_nref, SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, 0, 0, n_thr, &tab);
  int n = tab.n, D = tab.D;

  SEXP dims = getAttrib(sexp_X, R_DimSymbol);
  if (INTEGER(dims)[0] != n) error("The design must have one row per sample.");  // # nocov
  int P = INTEGER(dims)[1];

  int m = length(sexp_cmp);
  int    *q    = (int *)    R_alloc((size_t) m, sizeof(int));
  int    *qoff = (int *)    R_alloc((size_t) m, sizeof(int));
  const double **Ls = (const double **) R_alloc((size_t) m, sizeof(double *));
  const double **As = (const double **) R_alloc((size_t) m, sizeof(double *));
  double *crit = (double *) R_alloc((size_t) m, sizeof(double));
  int qsum = 0, qmax = 1;
  for (int c = 0; c < m; c++) {
    SEXP cm = VECTOR_ELT(sexp_cmp, c);
    q[c]    = ncols(VECTOR_ELT(cm, 0));
    qoff[c] = qsum;
    qsum   += q[c];
    if (q[c] > qmax) qmax = q[c];
    Ls[c]   = REAL(VECTOR_ELT(cm, 0));
    As[c]   = REAL(VECTOR_ELT(cm, 1));
    crit[c] = asReal(VECTOR_ELT(cm, 2));
  }

  // The design's QR, on this thread.
  size_t  msz   = flm_doubles(n, P, 0);
  double *fmem  = (double *) R_alloc(msz, sizeof(double));
  int    *fimem = (int *)    R_alloc((size_t) P + 1, sizeof(int));
  flm_t   full;
  flm_init(&full, fmem, fimem, n, P, 0);
  flm_decompose(&full, REAL(sexp_X), n, NULL, n, P, 1e-7);
  if (full.rank < P) error("The design is rank deficient.");  // # nocov: checked in R

  const char *nms[] = { "coef", "W", "est", "se", "stat", "ref", "mean", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_coef = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 0, s_coef);
  SEXP s_W    = allocMatrix(INTSXP,  D, m); SET_VECTOR_ELT(result, 1, s_W);
  SEXP s_est  = allocMatrix(REALSXP, D, m); SET_VECTOR_ELT(result, 2, s_est);
  SEXP s_se   = allocMatrix(REALSXP, D, m); SET_VECTOR_ELT(result, 3, s_se);
  SEXP s_stat = allocMatrix(REALSXP, D, m); SET_VECTOR_ELT(result, 4, s_stat);
  SEXP s_ref  = allocVector(VECSXP, m);     SET_VECTOR_ELT(result, 5, s_ref);
  SEXP s_mean = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 6, s_mean);

  // Each thread's buffer: 3 n + P doubles (phase 1: the feature, Q'y,
  // the residuals and the coefficients; phases 2 and 3: q <= P of them).
  int nt = threads_for(n_thr, (double) D * D * (n + qsum) / 2, 1e6);
  size_t bsz = 3 * (size_t) n + P;

  fa_tab   = &tab;
  fa_qr    = &full;
  fa_n     = n;
  fa_D     = D;
  fa_P     = P;
  fa_icol  = asInteger(sexp_icol);
  fa_m     = m;
  fa_qsum  = qsum;
  fa_pc    = asReal(sexp_pseudocount);
  fa_df    = (double) (n - P);
  fa_q     = q;
  fa_qoff  = qoff;
  fa_L     = Ls;
  fa_A     = As;
  fa_crit  = crit;
  fa_work  = (double *) R_alloc((size_t) nt * ((size_t) n + 3 * (size_t) P), sizeof(double));
  fa_buf   = (double *) R_alloc((size_t) nt * bsz, sizeof(double));
  fa_coef  = REAL(s_coef);
  fa_resid = (double *) R_alloc((size_t) n * D, sizeof(double));
  fa_b     = (double *) R_alloc((size_t) qsum * D, sizeof(double));
  fa_mean  = REAL(s_mean);
  fa_Wt    = (int *)    R_alloc((size_t) nt * D * m, sizeof(int));
  memset(fa_Wt, 0, (size_t) nt * D * m * sizeof(int));

  #define FA_PHASE(func, tasks) run_parallel(func, nt, tasks)

  FA_PHASE(fa_fit_worker, D < 100 ? 100 : D);
  if (user_interrupted()) { UNPROTECT(1); error("Interrupted."); }  // # nocov

  FA_PHASE(fa_pair_worker, D < 100 ? 100 : D);
  if (user_interrupted()) { UNPROTECT(1); error("Interrupted."); }  // # nocov

  int *W = INTEGER(s_W);
  memset(W, 0, (size_t) D * m * sizeof(int));
  for (int t = 0; t < nt; t++)
    for (size_t i = 0; i < (size_t) D * m; i++) W[i] += fa_Wt[(size_t) t * D * m + i];


  // Phase 3, per comparison.
  int nref = asInteger(sexp_nref);
  if (nref < 1) nref = 1;
  if (nref > D) nref = D;

  int    *cnt  = (int *)    R_alloc((size_t) D + 1, sizeof(int));
  int    *ord  = (int *)    R_alloc((size_t) D, sizeof(int));
  fa_ebar      = (double *) R_alloc((size_t) n, sizeof(double));
  fa_bbar      = (double *) R_alloc((size_t) qmax, sizeof(double));
  fa_est       = REAL(s_est);
  fa_se        = REAL(s_se);
  fa_stat      = REAL(s_stat);

  for (int c = 0; c < m; c++) {

    // A stable counting sort of the features by W.
    const int *Wc = W + (size_t) D * c;
    memset(cnt, 0, ((size_t) D + 1) * sizeof(int));
    for (int k = 0; k < D; k++) cnt[Wc[k] + 1]++;
    for (int w = 0; w < D; w++) cnt[w + 1] += cnt[w];
    for (int k = 0; k < D; k++) ord[cnt[Wc[k]]++] = k;

    SEXP s_r = allocVector(INTSXP, nref);
    SET_VECTOR_ELT(s_ref, c, s_r);
    for (int r = 0; r < nref; r++) INTEGER(s_r)[r] = ord[r] + 1;

    for (int s = 0; s < n; s++) {
      double sum = 0;
      for (int r = 0; r < nref; r++) sum += fa_resid[(size_t) n * ord[r] + s];
      fa_ebar[s] = sum / nref;
    }
    for (int i = 0; i < q[c]; i++) {
      double sum = 0;
      for (int r = 0; r < nref; r++) sum += fa_b[(size_t) qsum * ord[r] + qoff[c] + i];
      fa_bbar[i] = sum / nref;
    }

    fa_c = c;
    FA_PHASE(fa_effect_worker, D < 100 ? 100 : D);
  }

  #undef FA_PHASE

  UNPROTECT(1);
  return result;
}
