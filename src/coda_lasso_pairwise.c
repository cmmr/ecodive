// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// coda_lasso_pairwise(): the elastic net on all pairwise log ratios
// log(x_j / x_k), j < k (coda4microbiome's coda_glmnet() model), without
// forming the n x p(p-1)/2 matrix of log ratios. Column (j,k) is
// L_j - L_k (L = log x), so the gradient of pair (j,k) is G_j - G_k with
// G = L'Wr, and its variance under the weights is C_jj + C_kk - 2 C_jk
// with C the weighted covariance of L.
//
// The coordinate update (soft thresholding) and its convergence measure
// (the largest weighted squared change of a coefficient) are those of
// Friedman, Hastie & Tibshirani (2010); the screening uses the sequential
// strong rules (Tibshirani et al. 2012) with a KKT check over all pairs; the
// binomial family is fitted by iteratively reweighted least squares around
// it (Friedman et al. 2010). No code is ported. glmnet's conventions (the objective's scaling, the
// standardization of each column and, for gaussian y, of y, the lambda
// path and when it stops early) come from its documentation and Friedman
// et al. (2010), and were checked against its output.
//
// For samples i with weights w_i summing to 1 and pairs c = (j,k) with
// penalty weights a_c (sd of the column with standardize, else 1):
//
//   gaussian:  sum_i w_i r_i^2 / 2
//              + lambda sum_c [ alpha a_c |t_c| + (1 - alpha) a_c^2 t_c^2 / (2 s_y) ]
//   binomial:  -sum_i w_i (y_i eta_i - log(1 + exp(eta_i)))
//              + lambda sum_c [ alpha a_c |t_c| + (1 - alpha) a_c^2 t_c^2 / 2 ]
//
// with s_y the weighted sd of y (glmnet fits y / s_y with lambda / s_y).

#include "ecodive.h"
#include <float.h>  // DBL_EPSILON

#define PW_GAUSSIAN 0
#define PW_BINOMIAL 1
#define PW_PMIN 1e-5      // glmnet.control()'s pmin
#define PW_FDEV 1e-5      // glmnet.control()'s fdev
#define PW_DEVMAX 0.999   // glmnet.control()'s devmax
#define PW_MNLAM 5        // glmnet.control()'s mnlam


typedef struct {
  int     n, p, K, L, family, standardize, early_stop, maxit;
  double  alpha, tol;
  const double *x;          // n x p log abundances
  const double *y, *off;    // n
  const int    *foldid;     // n, 1..K
  const double *lambdas;    // L
  // Full-data fit:
  double *beta_path;        // p x L
  double *intercept_path;   // L
  double *dev_ratio;        // L
  // Pairs of the full-data fit, per lambda: grown by the full fit's thread.
  int    *pr_l, *pr_j, *pr_k; double *pr_t; size_t pr_n, pr_cap;
  // CV:
  double *cv_link;          // n x L
  int    *n_done;           // K + 1
  int    *n_iter_hit;       // K + 1: maxit reached
  int     next;
} pw_shared_t;

static pw_shared_t pw;


typedef struct {
  int j, k;
  double t;     // the coefficient
  double xv;    // sum_i w_i x_c^2 (weights of the current IRLS step)
  double a;     // penalty weight
} pw_pair_t;


typedef struct {
  int     f;
  double *w;        // observation weights (0 out of the training set)
  double *W;        // IRLS weights (gaussian: w)
  double *z;        // working response (without the offset)
  double *r;        // residuals z - b0 - sum_c t_c x_c
  double *eta;      // linear predictor (with the offset)
  double *mu;       // p: weighted means of L
  double *C;        // p x p weighted covariance (standardize only)
  double *G;        // p: L'Wr (centered columns)
  double *xc;       // n: one column
  double  b0, ys, nulldev;
  pw_pair_t *S; int nS, capS;   // the working set, sorted by (j, k)
  char   *ever;                 // p x p (upper): ever in the working set
  int     passes;
} pw_fit_t;


static inline double pw_a (const pw_fit_t *F, int j, int k) {
  if (!pw.standardize) return 1.0;
  int p = pw.p;
  double v = F->C[j + (size_t) j * p] + F->C[k + (size_t) k * p] - 2.0 * F->C[j + (size_t) k * p];
  return v > 0 ? sqrt(v) : 0.0;
}

static inline void pw_column (const pw_fit_t *F, int j, int k, double *xc) {
  int n = pw.n;
  const double *xj = pw.x + (size_t) j * n, *xk = pw.x + (size_t) k * n;
  double m = F->mu[j] - F->mu[k];
  for (int i = 0; i < n; i++) xc[i] = (xj[i] - xk[i]) - m;
}

static int pw_add (pw_fit_t *F, int j, int k) {
  // Insert (j,k) in order; 1 if new.
  int lo = 0, hi = F->nS;
  while (lo < hi) {
    int m = (lo + hi) >> 1;
    pw_pair_t *q = F->S + m;
    if (q->j < j || (q->j == j && q->k < k)) lo = m + 1; else hi = m;
  }
  if (lo < F->nS && F->S[lo].j == j && F->S[lo].k == k) return 0;
  if (F->nS == F->capS) {
    int cap = F->capS ? 2 * F->capS : 64;
    pw_pair_t *S = (pw_pair_t *) realloc(F->S, (size_t) cap * sizeof(pw_pair_t));
    if (!S) return -1;
    F->S = S; F->capS = cap;
  }
  memmove(F->S + lo + 1, F->S + lo, (size_t) (F->nS - lo) * sizeof(pw_pair_t));
  pw_pair_t *q = F->S + lo;
  q->j = j; q->k = k; q->t = 0.0; q->a = pw_a(F, j, k);
  pw_column(F, j, k, F->xc);
  double xv = 0.0;
  for (int i = 0; i < pw.n; i++) xv += F->W[i] * F->xc[i] * F->xc[i];
  q->xv = xv;
  F->nS++;
  F->ever[j + (size_t) k * pw.p] = 1;
  return 1;
}

// G = L_c' W r, with L centered by the observation-weighted means.
static void pw_gradient (pw_fit_t *F) {
  int n = pw.n, p = pw.p;
  for (int j = 0; j < p; j++) {
    const double *xj = pw.x + (size_t) j * n;
    double m = F->mu[j], g = 0.0;
    for (int i = 0; i < n; i++) g += F->W[i] * (xj[i] - m) * F->r[i];
    F->G[j] = g;
  }
}

// The working response, weights and residuals at the current fit.
static void pw_irls (pw_fit_t *F) {
  int n = pw.n;
  for (int i = 0; i < n; i++) {
    if (F->w[i] == 0) { F->W[i] = 0; F->z[i] = 0; continue; }
    double p = 1.0 / (1.0 + exp(-F->eta[i]));
    if (p < PW_PMIN) p = PW_PMIN; else if (p > 1 - PW_PMIN) p = 1 - PW_PMIN;
    double v = p * (1 - p);
    F->W[i] = F->w[i] * v;
    F->z[i] = (F->eta[i] - pw.off[i]) + (pw.y[i] - p) / v;
  }
  // r = z - (eta - off)
  for (int i = 0; i < n; i++) F->r[i] = F->z[i] - (F->eta[i] - pw.off[i]);
  for (int c = 0; c < F->nS; c++) {
    pw_pair_t *q = F->S + c;
    pw_column(F, q->j, q->k, F->xc);
    double xv = 0.0;
    for (int i = 0; i < n; i++) xv += F->W[i] * F->xc[i] * F->xc[i];
    q->xv = xv;
  }
}

static double pw_deviance (const pw_fit_t *F) {
  int n = pw.n;
  double d = 0.0;
  if (pw.family == PW_GAUSSIAN) {
    for (int i = 0; i < n; i++) {
      if (F->w[i] == 0) continue;
      double e = pw.y[i] - F->eta[i];
      d += F->w[i] * e * e;
    }
    return d;
  }
  for (int i = 0; i < n; i++) {
    if (F->w[i] == 0) continue;
    double e = F->eta[i];
    // -2 loglik = 2 (log(1 + exp(e)) - y e)
    double l1pe = e > 0 ? e + log1p(exp(-e)) : log1p(exp(e));
    d += 2.0 * F->w[i] * (l1pe - pw.y[i] * e);
  }
  return d;
}

// Coordinate descent on the working set until the largest weighted squared
// change falls below tol * nulldev; returns the passes made, or -1 at maxit.
static int pw_cd (pw_fit_t *F, double lambda) {
  int n = pw.n;
  double l1 = lambda * pw.alpha;
  double l2 = lambda * (1.0 - pw.alpha) / (pw.family == PW_GAUSSIAN ? F->ys : 1.0);
  double thr = pw.tol * F->nulldev;
  int binom = pw.family == PW_BINOMIAL;
  for (int pass = 0; ; pass++) {
    if (F->passes++ >= pw.maxit) return -1;
    double dmax = 0.0;
    if (binom) {
      // The intercept.
      double sw = 0.0, sr = 0.0;
      for (int i = 0; i < n; i++) { sw += F->W[i]; sr += F->W[i] * F->r[i]; }
      if (sw > 0) {
        double d = sr / sw;
        F->b0 += d;
        for (int i = 0; i < n; i++) F->r[i] -= d;
        if (sw * d * d > dmax) dmax = sw * d * d;
      }
    }
    for (int c = 0; c < F->nS; c++) {
      pw_pair_t *q = F->S + c;
      if (q->xv <= 0) continue;
      pw_column(F, q->j, q->k, F->xc);
      double g = 0.0;
      for (int i = 0; i < n; i++) g += F->W[i] * F->xc[i] * F->r[i];
      double u = g + q->xv * q->t;
      double thrs = l1 * q->a;
      double tn = fabs(u) <= thrs ? 0.0 : (u > 0 ? u - thrs : u + thrs) / (q->xv + l2 * q->a * q->a);
      double d = tn - q->t;
      if (d != 0.0) {
        q->t = tn;
        for (int i = 0; i < n; i++) F->r[i] -= d * F->xc[i];
        double ch = q->xv * d * d;
        if (ch > dmax) dmax = ch;
      }
    }
    if (dmax < thr) return pass + 1;
  }
}

static void pw_eta (pw_fit_t *F) {
  // eta = off + b0 + sum_c t_c x_c  (from the working set)
  int n = pw.n;
  for (int i = 0; i < n; i++) F->eta[i] = pw.off[i] + F->b0;
  for (int c = 0; c < F->nS; c++) {
    pw_pair_t *q = F->S + c;
    if (q->t == 0.0) continue;
    pw_column(F, q->j, q->k, F->xc);
    for (int i = 0; i < n; i++) F->eta[i] += q->t * F->xc[i];
  }
}

// Fits one lambda; returns 0, 1 at maxit, -1 on allocation failure.
static int pw_fit_lambda (pw_fit_t *F, double lambda, double lambda_prev) {
  int p = pw.p;
  int binom = pw.family == PW_BINOMIAL;

  // The strong set, from the gradient at the previous solution.
  if (binom) pw_irls(F);
  pw_gradient(F);
  double cut = pw.alpha * (2.0 * lambda - lambda_prev);
  for (int j = 0; j < p - 1; j++)
    for (int k = j + 1; k < p; k++) {
      if (F->ever[j + (size_t) k * p]) continue;
      double a = pw_a(F, j, k);
      if (fabs(F->G[j] - F->G[k]) >= cut * a && a > 0)
        if (pw_add(F, j, k) < 0) return -1;
    }

  for (;;) {
    if (binom) {
      // IRLS around coordinate descent.
      for (int it = 0; ; it++) {
        pw_irls(F);
        int np = pw_cd(F, lambda);
        if (np < 0) return 1;
        pw_eta(F);
        if (np == 1) break;
      }
      pw_irls(F);
    }
    else {
      if (pw_cd(F, lambda) < 0) return 1;
      pw_eta(F);
    }

    // KKT over the pairs outside the working set.
    pw_gradient(F);
    int added = 0;
    double l1 = lambda * pw.alpha;
    for (int j = 0; j < p - 1; j++)
      for (int k = j + 1; k < p; k++) {
        if (F->ever[j + (size_t) k * p]) continue;
        double a = pw_a(F, j, k);
        if (a > 0 && fabs(F->G[j] - F->G[k]) > l1 * a) {
          int s = pw_add(F, j, k);
          if (s < 0) return -1;
          added += s;
        }
      }
    if (!added) break;
  }
  return 0;
}


static int pw_path (int f) {
  int n = pw.n, p = pw.p, K = pw.K, L = pw.L;
  pw_fit_t F;
  memset(&F, 0, sizeof(F));
  F.f = f;
  double *mem = (double *) calloc((size_t) 7 * n + 2 * (size_t) p + 1, sizeof(double));
  F.ever = (char *) calloc((size_t) p * p, 1);
  if (pw.standardize) F.C = (double *) calloc((size_t) p * p, sizeof(double));
  if (!mem || !F.ever || (pw.standardize && !F.C)) { free(mem); free(F.ever); free(F.C); return -1; }
  F.w = mem; F.W = F.w + n; F.z = F.W + n; F.r = F.z + n; F.eta = F.r + n;
  F.xc = F.eta + n; F.mu = F.xc + 2 * (size_t) n; F.G = F.mu + p;

  int nt = 0;
  for (int i = 0; i < n; i++) if (f == K || pw.foldid[i] - 1 != f) nt++;
  for (int i = 0; i < n; i++) F.w[i] = (f == K || pw.foldid[i] - 1 != f) ? 1.0 / nt : 0.0;

  for (int j = 0; j < p; j++) {
    const double *xj = pw.x + (size_t) j * n;
    double m = 0.0;
    for (int i = 0; i < n; i++) m += F.w[i] * xj[i];
    F.mu[j] = m;
  }
  if (pw.standardize) {
    for (int j = 0; j < p; j++) {
      const double *xj = pw.x + (size_t) j * n;
      for (int k = j; k < p; k++) {
        const double *xk = pw.x + (size_t) k * n;
        double s = 0.0;
        for (int i = 0; i < n; i++) s += F.w[i] * (xj[i] - F.mu[j]) * (xk[i] - F.mu[k]);
        F.C[j + (size_t) k * p] = F.C[k + (size_t) j * p] = s;
      }
    }
  }

  // The null model (intercept with the offset).
  if (pw.family == PW_GAUSSIAN) {
    double ym = 0.0;
    for (int i = 0; i < n; i++) ym += F.w[i] * (pw.y[i] - pw.off[i]);
    double v = 0.0;
    for (int i = 0; i < n; i++) if (F.w[i] > 0) { double e = pw.y[i] - pw.off[i] - ym; v += F.w[i] * e * e; }
    F.b0 = ym; F.ys = sqrt(v); F.nulldev = v;
    for (int i = 0; i < n; i++) {
      F.W[i] = F.w[i];
      F.z[i] = pw.y[i] - pw.off[i];
      F.eta[i] = pw.off[i] + ym;
      F.r[i] = F.z[i] - ym;
    }
  }
  else {
    double b = 0.0, ybar = 0.0;
    for (int i = 0; i < n; i++) ybar += F.w[i] * pw.y[i];
    b = log(ybar / (1 - ybar));
    for (int it = 0; it < 100; it++) {
      double g = 0.0, h = 0.0;
      for (int i = 0; i < n; i++) {
        if (F.w[i] == 0) continue;
        double q = 1.0 / (1.0 + exp(-(pw.off[i] + b)));
        g += F.w[i] * (pw.y[i] - q); h += F.w[i] * q * (1 - q);
      }
      double d = g / h;
      b += d;
      if (fabs(d) < 1e-12) break;
    }
    F.b0 = b;
    for (int i = 0; i < n; i++) F.eta[i] = pw.off[i] + b;
    F.nulldev = pw_deviance(&F);
    pw_irls(&F);
  }
  F.ys = pw.family == PW_GAUSSIAN ? F.ys : 1.0;

  int done = L, hit = 0;
  double dr_prev = 0.0;
  double *beta = (double *) malloc((size_t) p * sizeof(double));
  if (!beta) { free(mem); free(F.ever); free(F.C); free(F.S); return -1; }
  for (int l = 0; l < L; l++) {
    double lam = pw.lambdas[l];
    double lprev = l == 0 ? lam : pw.lambdas[l - 1];
    int st = pw_fit_lambda(&F, lam, lprev);
    if (st < 0) { free(mem); free(F.ever); free(F.C); free(F.S); free(beta); return -1; }
    if (st == 1) { hit = 1; done = l; break; }

    if (f < K) {
      for (int i = 0; i < n; i++)
        if (pw.foldid[i] - 1 == f) pw.cv_link[i + (size_t) l * n] = F.eta[i];
    }
    else {
      for (int j = 0; j < p; j++) beta[j] = 0.0;
      double b0 = F.b0;
      for (int c = 0; c < F.nS; c++) {
        pw_pair_t *q = F.S + c;
        if (q->t == 0.0) continue;
        beta[q->j] += q->t; beta[q->k] -= q->t;
        b0 -= q->t * (F.mu[q->j] - F.mu[q->k]);
        if (pw.pr_n == pw.pr_cap) {
          size_t cap = pw.pr_cap ? 2 * pw.pr_cap : 1024;
          int *a1 = (int *) realloc(pw.pr_l, cap * sizeof(int));
          if (a1) pw.pr_l = a1;
          int *a2 = (int *) realloc(pw.pr_j, cap * sizeof(int));
          if (a2) pw.pr_j = a2;
          int *a3 = (int *) realloc(pw.pr_k, cap * sizeof(int));
          if (a3) pw.pr_k = a3;
          double *a4 = (double *) realloc(pw.pr_t, cap * sizeof(double));
          if (a4) pw.pr_t = a4;
          if (!a1 || !a2 || !a3 || !a4) { free(mem); free(F.ever); free(F.C); free(F.S); free(beta); return -1; }
          pw.pr_cap = cap;
        }
        pw.pr_l[pw.pr_n] = l; pw.pr_j[pw.pr_n] = q->j; pw.pr_k[pw.pr_n] = q->k;
        pw.pr_t[pw.pr_n] = q->t; pw.pr_n++;
      }
      memcpy(pw.beta_path + (size_t) l * p, beta, (size_t) p * sizeof(double));
      pw.intercept_path[l] = b0;
    }

    // glmnet's early stop of a default path.
    double dr = 1.0 - pw_deviance(&F) / F.nulldev;
    if (f == K) pw.dev_ratio[l] = dr;
    if (pw.early_stop && l + 1 >= PW_MNLAM) {
      int stop = dr > PW_DEVMAX;
      if (pw.family == PW_GAUSSIAN) stop = stop || (dr - dr_prev < PW_FDEV * dr);
      else                          stop = stop || (dr - dr_prev < PW_FDEV);
      if (stop) { done = l + 1; break; }
    }
    dr_prev = dr;
  }

  // Beyond a fold's stop, its last fit predicts (as glmnet's predict()).
  if (f < K && done > 0)
    for (int l = done; l < L; l++)
      for (int i = 0; i < n; i++)
        if (pw.foldid[i] - 1 == f)
          pw.cv_link[i + (size_t) l * n] = pw.cv_link[i + (size_t) (done - 1) * n];

  pw.n_done[f] = done;
  pw.n_iter_hit[f] = hit;
  free(mem); free(F.ever); free(F.C); free(F.S); free(beta);
  return 0;
}


static void *pw_worker (void *arg) {
  (void) arg;
  for (;;) {
    int f = __atomic_fetch_add(&pw.next, 1, __ATOMIC_RELAXED);
    if (f > pw.K) break;
    if (pw_path(f)) pw.n_done[f] = -1;
  }
  return NULL;
}


//======================================================
// R interface.
//
// Returns list(beta (p x L), intercept (L), dev_ratio (L),
// cv_link (n x L), n_done (K+1: the full fit last),
// maxit (K+1), pairs = data.frame-like list(l, j, k, t),
// 1-based).
//======================================================
SEXP C_coda_lasso_pairwise (
    SEXP sexp_x, SEXP sexp_y, SEXP sexp_offset, SEXP sexp_foldid, SEXP sexp_K,
    SEXP sexp_family, SEXP sexp_alpha, SEXP sexp_lambdas, SEXP sexp_standardize,
    SEXP sexp_tol, SEXP sexp_maxit, SEXP sexp_early_stop, SEXP sexp_n_threads ) {

  int n = nrows(sexp_x), p = ncols(sexp_x);
  int K = asInteger(sexp_K), L = length(sexp_lambdas);
  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  SEXP sexp_beta = PROTECT(allocMatrix(REALSXP, p, L));
  SEXP sexp_icpt = PROTECT(allocVector(REALSXP, L));
  SEXP sexp_dr   = PROTECT(allocVector(REALSXP, L));
  SEXP sexp_link = PROTECT(allocMatrix(REALSXP, n, L));
  SEXP sexp_done = PROTECT(allocVector(INTSXP, K + 1));
  SEXP sexp_hit  = PROTECT(allocVector(INTSXP, K + 1));
  memset(REAL(sexp_beta), 0, (size_t) p * L * sizeof(double));
  for (int l = 0; l < L; l++) { REAL(sexp_icpt)[l] = NA_REAL; REAL(sexp_dr)[l] = NA_REAL; }
  for (size_t i = 0; i < (size_t) n * L; i++) REAL(sexp_link)[i] = NA_REAL;

  memset(&pw, 0, sizeof(pw));
  pw.n = n; pw.p = p; pw.K = K; pw.L = L;
  pw.family = asInteger(sexp_family);
  pw.standardize = asInteger(sexp_standardize);
  pw.early_stop = asInteger(sexp_early_stop);
  pw.maxit = asInteger(sexp_maxit);
  pw.alpha = asReal(sexp_alpha);
  pw.tol = asReal(sexp_tol);
  pw.x = REAL(sexp_x); pw.y = REAL(sexp_y); pw.off = REAL(sexp_offset);
  pw.foldid = INTEGER(sexp_foldid); pw.lambdas = REAL(sexp_lambdas);
  pw.beta_path = REAL(sexp_beta); pw.intercept_path = REAL(sexp_icpt);
  pw.dev_ratio = REAL(sexp_dr); pw.cv_link = REAL(sexp_link);
  pw.n_done = INTEGER(sexp_done); pw.n_iter_hit = INTEGER(sexp_hit);

  int nt = n_threads < K + 1 ? n_threads : K + 1;
  run_parallel(pw_worker, nt, nt > 1 ? 100 : 1);

  int failed = 0;
  for (int f = 0; f <= K; f++) if (pw.n_done[f] < 0) failed = 1;
  if (failed) {
    free(pw.pr_l); free(pw.pr_j); free(pw.pr_k); free(pw.pr_t);
    error("Memory allocation failed.");  // # nocov
  }

  // Move the full fit's done count first.
  SEXP sexp_pairs = PROTECT(allocVector(VECSXP, 4));
  SEXP pl = PROTECT(allocVector(INTSXP, (R_xlen_t) pw.pr_n));
  SEXP pj = PROTECT(allocVector(INTSXP, (R_xlen_t) pw.pr_n));
  SEXP pk = PROTECT(allocVector(INTSXP, (R_xlen_t) pw.pr_n));
  SEXP pt = PROTECT(allocVector(REALSXP, (R_xlen_t) pw.pr_n));
  for (size_t i = 0; i < pw.pr_n; i++) {
    INTEGER(pl)[i] = pw.pr_l[i] + 1; INTEGER(pj)[i] = pw.pr_j[i] + 1;
    INTEGER(pk)[i] = pw.pr_k[i] + 1; REAL(pt)[i] = pw.pr_t[i];
  }
  free(pw.pr_l); free(pw.pr_j); free(pw.pr_k); free(pw.pr_t);
  SET_VECTOR_ELT(sexp_pairs, 0, pl); SET_VECTOR_ELT(sexp_pairs, 1, pj);
  SET_VECTOR_ELT(sexp_pairs, 2, pk); SET_VECTOR_ELT(sexp_pairs, 3, pt);
  SEXP pn = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(pn, 0, mkChar("l")); SET_STRING_ELT(pn, 1, mkChar("j"));
  SET_STRING_ELT(pn, 2, mkChar("k")); SET_STRING_ELT(pn, 3, mkChar("t"));
  setAttrib(sexp_pairs, R_NamesSymbol, pn);

  SEXP res = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(res, 0, sexp_beta);
  SET_VECTOR_ELT(res, 1, sexp_icpt);
  SET_VECTOR_ELT(res, 2, sexp_dr);
  SET_VECTOR_ELT(res, 3, sexp_link);
  SET_VECTOR_ELT(res, 4, sexp_done);
  SET_VECTOR_ELT(res, 5, sexp_hit);
  SET_VECTOR_ELT(res, 6, sexp_pairs);
  SEXP nms = PROTECT(allocVector(STRSXP, 7));
  const char *nm[] = { "beta", "intercept", "dev_ratio", "cv_link", "n_done", "maxit", "pairs" };
  for (int i = 0; i < 7; i++) SET_STRING_ELT(nms, i, mkChar(nm[i]));
  setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(14);
  return res;
}
