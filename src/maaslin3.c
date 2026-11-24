// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// maaslin3() is a C port of the fitting steps of MaAsLin 3 (maaslin3 1.4.0,
// Bioconductor: R/utility_scripts.R TSSnorm(), CLRnorm(), NONEnorm(),
// LOG(); R/fit.R augment_data(), check_for_zero_one_obs(),
// check_missing_first_factor_level(), fit.model() and
// run_median_comparison_general()), restructured for ecodive (see
// plans/maaslin3.md), under this license:
//
//   Copyright (c) 2025 Harvard School of Public Health
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
// The rest of this file is Copyright 2026 ecodive authors, under the MIT
// license above. maaslin3 fits one R model per feature (lm(), glm(),
// lmerTest::lmer(), lme4::glmer(); those packages are GPL and their code
// was not read); here:
//
//   * Abundance, fixed effects: least squares on the feature's non-zero
//     samples, with feature_lm.c's pivoted QR (lm()'s tolerance, 1e-7).
//   * Abundance, a random intercept: REML at the exact optimum with
//     Satterthwaite's degrees of freedom (feature_lmm.c), on the feature's
//     non-zero samples and the design columns of full rank among them.
//   * Prevalence, fixed effects: iteratively reweighted least squares for
//     a logistic model, as glm() documents it (convergence when |dev -
//     dev_old| / (|dev| + 0.1) < 1e-8, at most 25 iterations, starting
//     values (w y + 1/2) / (w + 1)), with maaslin3's augmentation: every
//     sample twice more, once present and once absent, each with weight a =
//     p / (2 n). The copies share the sample's covariates, so each
//     iteration's weighted least squares collapses them into one row per
//     sample (the same normal equations). Standard errors come from the
//     last iteration's weighted QR (dispersion 1), as summary.glm()'s.
//     The logit's inverse and derivative saturate beyond |eta| = 30, as R's
//     binomial() family's do (settled from their outputs).
//   * Prevalence, a random intercept: the Laplace approximation of lme4's
//     glmer() (nAGQ = 1), minimised exactly over (theta, beta) by Newton's
//     method with analytic gradients. With the augmentation's fractional
//     weights, glmer's criterion (found from its devFunOnly outputs) is
//     -2 log L of the original samples only (the copies' binomial
//     likelihood has round(a) = 0 trials), plus |u|^2 and log det(I +
//     theta^2 Z'WZ), with the conditional modes u and the working weights
//     W of the augmented data. Standard errors: the fixed-effects block of
//     2 Hess^-1 of that criterion in (theta, beta), as glmer's (also
//     settled from outputs).
//   * The median comparison's Monte Carlo: maaslin3 draws 10,000 sets of
//     coefficients from R's rnorm(); here each draw has its own PCG stream
//     (rng_dist.c), keyed by (seed, covariate, draw), so the results are
//     the same for any number of threads. Its covariances are summed in
//     blocks of draws, in a fixed order.

#include "ecodive.h"

#define M3_NORM_NONE    0
#define M3_NORM_PERCENT 1
#define M3_NORM_RCLR    2

// Feature statuses (all of the feature's coefficients).
#define M3_OK          0
#define M3_SAME        1   // all values the same
#define M3_BASELINE    2   // no non-zero sample at a factor's reference level
#define M3_ONE_GROUP   3   // the random intercept: fewer than two groups
#define M3_MANY_GROUPS 4   //   as many groups as samples, or more
#define M3_FAILED      5   // the fit failed
#define M3_NO_CONV     6   // the mixed logistic model did not converge
#define M3_ONE_LEVEL   7   // a factor with one level among the samples fitted


//======================================================
// The table, by feature, and each sample's divisor:
// its total (percent) or the geometric mean of its
// non-zero values (rclr).
//======================================================
static da_table_t m3_tab;
static double    *m3_div;
static int        m3_norm, m3_log;


static void m3_prepare (SEXP sexp_counts, SEXP sexp_margin, int norm, int n_threads) {

  da_table_init(sexp_counts, sexp_margin, 0, 0, n_threads, &m3_tab);
  int n = m3_tab.n, D = m3_tab.D;
  m3_div = (double *) R_alloc((size_t) n + 1, sizeof(double));

  if (norm == M3_NORM_PERCENT) {
    // TSSnorm(): x / max(eps, rowSums(x)), summed in feature order.
    long double *tot = (long double *) R_alloc((size_t) n + 1, sizeof(long double));
    for (int s = 0; s < n; s++) tot[s] = 0;
    for (int f = 0; f < D; f++)
      for (int e = m3_tab.fpos[f]; e < m3_tab.fpos[f + 1]; e++)
        tot[m3_tab.fsam[e]] += m3_tab.fval[e];
    for (int s = 0; s < n; s++) {
      double t = (double) tot[s];
      m3_div[s] = t > DBL_EPSILON ? t : DBL_EPSILON;
    }
  }

  else if (norm == M3_NORM_RCLR) {
    // CLRnorm(): x / exp(1)^mean(log(x[x > 0])), the mean as R's mean().
    long double *sum = (long double *) R_alloc((size_t) n + 1, sizeof(long double));
    int         *cnt = (int *)         R_alloc((size_t) n + 1, sizeof(int));
    for (int s = 0; s < n; s++) { sum[s] = 0; cnt[s] = 0; }
    for (int f = 0; f < D; f++)
      for (int e = m3_tab.fpos[f]; e < m3_tab.fpos[f + 1]; e++)
        if (m3_tab.fval[e] > 0) { sum[m3_tab.fsam[e]] += log(m3_tab.fval[e]); cnt[m3_tab.fsam[e]]++; }
    long double *m = (long double *) R_alloc((size_t) n + 1, sizeof(long double));
    for (int s = 0; s < n; s++) { m[s] = cnt[s] ? sum[s] / cnt[s] : 0; sum[s] = 0; }
    for (int f = 0; f < D; f++)
      for (int e = m3_tab.fpos[f]; e < m3_tab.fpos[f + 1]; e++)
        if (m3_tab.fval[e] > 0) sum[m3_tab.fsam[e]] += log(m3_tab.fval[e]) - m[m3_tab.fsam[e]];
    for (int s = 0; s < n; s++) {
      double gm = cnt[s] ? (double) (m[s] + sum[s] / cnt[s]) : R_NaN;
      m3_div[s] = pow(M_E, gm);
    }
  }

  else {
    for (int s = 0; s < n; s++) m3_div[s] = 1;
  }

  m3_norm = norm;
}


// Feature f in every sample: its normalized (and log2) value where it is
// present, NaN where not (abundance); or 1 / 0 (prevalence).
static void m3_column (int f, double *y, int prevalence) {
  int n = m3_tab.n;
  for (int s = 0; s < n; s++) y[s] = prevalence ? 0 : R_NaN;
  for (int e = m3_tab.fpos[f]; e < m3_tab.fpos[f + 1]; e++) {
    double v = m3_tab.fval[e];
    if (!(v > 0)) continue;
    int s = m3_tab.fsam[e];
    if (prevalence) { y[s] = 1; continue; }
    if      (m3_norm == M3_NORM_PERCENT) v = v / m3_div[s];
    else if (m3_norm == M3_NORM_RCLR)    v = log(v / m3_div[s]);
    if (m3_log) v = log2(v);
    y[s] = v;
  }
}



//======================================================
// The logistic model's pieces, as R's binomial()
// family's (its outputs): the inverse link and its
// derivative saturate beyond |eta| = 30.
//======================================================
static inline double lg_inv (double eta) {
  double e = eta < -30 ? DBL_EPSILON : (eta > 30 ? 1 / DBL_EPSILON : exp(eta));
  return e / (1 + e);
}

static inline double lg_mueta (double eta) {
  if (eta > 30 || eta < -30) return DBL_EPSILON;
  double e = exp(eta);
  return e / ((1 + e) * (1 + e));
}

// The binomial deviance residual, 2 w (y log(y / mu) + (1 - y) log((1 - y) / (1 - mu))).
static inline double bin_dev (double y, double mu, double w) {
  double r = 0;
  if (y != 0) r += y * log(y / mu);
  if (y != 1) r += (1 - y) * log((1 - y) / (1 - mu));
  return 2 * w * r;
}



//======================================================
// Shared state of the workers.
//======================================================
static const double *m3_X;          // n x P design (all samples)
static int           m3_P;
static const int    *m3_grp;        // n: random intercept groups (0-based), or NULL
static int           m3_J;
static const int    *m3_base;       // n x K: each sample's level of factor k (1: the reference)
static int           m3_K;
static double        m3_a;          // the augmentation's weight (0: none)
static int           m3_prev;       // 1: prevalence model

static double *m3_coef, *m3_se, *m3_df, *m3_xsd;   // D x P
static double *m3_nobs, *m3_mean, *m3_rsd, *m3_theta;
static int    *m3_status, *m3_warn;

static double *m3_mem;    static size_t m3_mem_sz;    // per thread doubles
static int    *m3_imem;   static size_t m3_imem_sz;   // per thread ints
static long double *m3_lmem;                          // per thread J long doubles


// The standard deviation of x[rows] (or x, rows NULL), as sd(); NaN when
// every value is 0 or 1 (maaslin3 then leaves its threshold unscaled).
// `copies` > 1: the values repeated that many times (the augmented data).
static double m3_xsd_of (const double *x, const int *rows, int n, int copies) {
  int binary = 1;
  long double s = 0;
  for (int i = 0; i < n; i++) {
    double v = rows ? x[rows[i]] : x[i];
    if (v != 0 && v != 1) binary = 0;
    s += v;
  }
  if (binary || n < 2) return R_NaN;
  long double m = s / n, ss = 0;
  for (int i = 0; i < n; i++) { double v = (rows ? x[rows[i]] : x[i]) - (double) m; ss += (long double) v * v; }
  return sqrt((double) (ss * copies / ((long double) n * copies - 1)));
}



//======================================================
// One feature's logistic regression (glm(), with
// maaslin3's augmentation when a > 0). Writes coef, se
// (P), returns the warnings: 1 not converged, 2 fitted
// probabilities numerically 0 or 1.
//======================================================
typedef struct {
  double *xw, *zw, *qty, *eta, *mu, *b;    // n*P, n, n, 3n, 3n, P
  flm_t   q;
} glm_work_t;


static int m3_glm (
    const double *X, int n, int P, const double *y, double a, glm_work_t *w,
    double *coef, double *se ) {

  int    nc = a > 0 ? 3 : 1;          // copies of each sample
  double wr[3] = { 1, a, a };
  double *eta = w->eta, *mu = w->mu;  // per copy (first iteration), else copy 0

  #define YR(r, s) ((r) == 0 ? y[s] : ((r) == 1 ? 1.0 : 0.0))

  // Starting values: mustart = (w y + 1/2) / (w + 1), eta = logit(mustart).
  long double dold = 0;
  for (int r = 0; r < nc; r++)
    for (int s = 0; s < n; s++) {
      double ms = (wr[r] * YR(r, s) + 0.5) / (wr[r] + 1);
      double e  = log(ms / (1 - ms));
      eta[s + (size_t) n * r] = e;
      mu[s + (size_t) n * r]  = lg_inv(e);
      dold += bin_dev(YR(r, s), mu[s + (size_t) n * r], wr[r]);
    }
  double devold = (double) dold;

  int conv = 0, shared = 0;
  for (int it = 1; it <= 25; it++) {

    for (int s = 0; s < n; s++) {
      double W = 0, WZ = 0;
      for (int r = 0; r < nc; r++) {
        size_t i  = s + (size_t) n * (shared ? 0 : r);
        double me = lg_mueta(eta[i]);
        double v  = mu[i] * (1 - mu[i]);
        double z  = eta[i] + (YR(r, s) - mu[i]) / me;
        double w2 = wr[r] * me * me / v;
        W  += w2;
        WZ += w2 * z;
      }
      double sw = sqrt(W);
      for (int j = 0; j < P; j++) w->xw[s + (size_t) n * j] = sw * X[s + (size_t) n * j];
      w->zw[s] = sw * (WZ / W);
    }

    double rss;
    flm_decompose(&w->q, w->xw, n, NULL, n, P, 1e-11);
    flm_fit(&w->q, w->zw, w->qty, w->b, &rss);

    long double dv = 0;
    for (int s = 0; s < n; s++) {
      double e = 0;
      for (int j = 0; j < P; j++) if (!ISNAN(w->b[j])) e += X[s + (size_t) n * j] * w->b[j];
      eta[s] = e;
      mu[s]  = lg_inv(e);
    }
    shared = 1;
    for (int r = 0; r < nc; r++)
      for (int s = 0; s < n; s++) dv += bin_dev(YR(r, s), mu[s], wr[r]);
    double dev = (double) dv;

    if (fabs(dev - devold) / (fabs(dev) + 0.1) < 1e-8) { conv = 1; break; }
    devold = dev;
  }
  #undef YR

  // Standard errors from the last weighted QR, dispersion 1.
  int rank = w->q.rank;
  for (int j = 0; j < P; j++) { coef[j] = w->b[j]; se[j] = R_NaN; }
  for (int k = 0; k < rank; k++) {
    long double ss = 0;
    for (int m = k; m < rank; m++) { double v = w->q.rinv[k + (size_t) P * m]; ss += (long double) v * v; }
    se[w->q.pivot[k]] = sqrt((double) ss);
  }

  int warn = conv ? 0 : 1;
  double eps10 = 10 * DBL_EPSILON;
  for (int s = 0; s < n; s++) if (mu[s] > 1 - eps10 || mu[s] < eps10) { warn |= 2; break; }
  return warn;
}



//======================================================
// One feature's logistic model with a random intercept
// (glmer(), Laplace). par = (theta, beta).
//======================================================
typedef struct {
  const double *X;     // n x P
  const double *y;     // n: 0 / 1
  const int    *gidx;  // samples ordered by group
  const int    *gpos;  // J + 1: group j is gidx[gpos[j] .. gpos[j + 1] - 1]
  int     n, P, J;
  double  a;           // augmentation weight
  double *u;           // J: conditional modes (warm starts)
  double *xb;          // n
  double *g2;          // P + 1 scratch
  double *ub;          // P scratch
} glmm_t;


// The criterion at par; its gradient in g (if not NULL).
static double glmm_eval (glmm_t *m, const double *par, double *g) {

  int n = m->n, P = m->P, J = m->J;
  double th = par[0];
  const double *beta = par + 1;
  double W = 1 + 2 * m->a;

  for (int s = 0; s < n; s++) {
    double e = 0;
    for (int k = 0; k < P; k++) e += m->X[s + (size_t) n * k] * beta[k];
    m->xb[s] = e;
  }
  if (g) for (int k = 0; k <= P; k++) g[k] = 0;

  long double D = 0;

  for (int j = 0; j < J; j++) {
    int lo = m->gpos[j], hi = m->gpos[j + 1];

    // The conditional mode: u = theta sum W (yt - mu), by Newton.
    double u = m->u[j];
    for (int it = 0; it < 100; it++) {
      double F = u, Fu = 1;
      for (int i = lo; i < hi; i++) {
        int s = m->gidx[i];
        double mu = lg_inv(m->xb[s] + th * u);
        F  -= th * ((m->y[s] + m->a) - W * mu);
        Fu += th * th * W * mu * (1 - mu);
      }
      double st = F / Fu;
      if (st > 5) st = 5; else if (st < -5) st = -5;
      u -= st;
      if (fabs(st) <= 1e-15 * (1 + fabs(u))) break;
    }
    m->u[j] = u;

    // The criterion and the pieces of its gradient.
    long double ll = 0, S = 0, Swr = 0;
    for (int k = 0; k < P; k++) m->ub[k] = 0;
    for (int i = lo; i < hi; i++) {
      int s = m->gidx[i];
      double mu = lg_inv(m->xb[s] + th * u);
      ll += m->y[s] ? log(mu) : log(1 - mu);
      double v = mu * (1 - mu);
      S   += W * v;
      Swr += (m->y[s] + m->a) - W * mu;
      if (g) for (int k = 0; k < P; k++) m->ub[k] += W * v * m->X[s + (size_t) n * k];
    }
    double Sd = (double) S;
    double L  = 1 + th * th * Sd;
    D += -2 * ll + u * u + log(L);

    if (g) {
      double Fu  = L;
      double c   = th * th / L;
      double u_t = ((double) Swr - th * u * Sd) / Fu;
      // u_beta = -theta ub / Fu
      for (int k = 0; k < P; k++) m->ub[k] *= -th / Fu;

      // sum_s r_s (x_s + theta u_beta) and sum_s r_s (u + theta u_theta),
      // with r_s = -2 (y - mu) + c W v (1 - 2 mu).
      long double gt = 0;
      for (int k = 0; k < P; k++) m->g2[k] = 0;
      double rsum = 0;
      for (int i = lo; i < hi; i++) {
        int s = m->gidx[i];
        double mu = lg_inv(m->xb[s] + th * u);
        double v  = mu * (1 - mu);
        double r  = -2 * (m->y[s] - mu) + c * W * v * (1 - 2 * mu);
        rsum += r;
        for (int k = 0; k < P; k++) m->g2[k] += r * m->X[s + (size_t) n * k];
      }
      for (int k = 0; k < P; k++)
        g[k + 1] += m->g2[k] + rsum * th * m->ub[k] + 2 * u * m->ub[k];
      gt = rsum * (u + th * u_t) + 2 * u * u_t + 2 * th * Sd / L;
      g[0] += (double) gt;
    }
  }

  return (double) D;
}


// Cholesky solve of the k x k system H x = b (H overwritten by its factor);
// 0 if H isn't positive definite.
static int m3_chol_solve (double *H, int k, double *b) {
  for (int j = 0; j < k; j++) {
    double d = H[j + (size_t) k * j];
    for (int m = 0; m < j; m++) d -= H[j + (size_t) k * m] * H[j + (size_t) k * m];
    if (!(d > 0)) return 0;
    d = sqrt(d);
    H[j + (size_t) k * j] = d;
    for (int i = j + 1; i < k; i++) {
      double v = H[i + (size_t) k * j];
      for (int m = 0; m < j; m++) v -= H[i + (size_t) k * m] * H[j + (size_t) k * m];
      H[i + (size_t) k * j] = v / d;
    }
  }
  if (b) {
    for (int i = 0; i < k; i++) {
      double v = b[i];
      for (int m = 0; m < i; m++) v -= H[i + (size_t) k * m] * b[m];
      b[i] = v / H[i + (size_t) k * i];
    }
    for (int i = k - 1; i >= 0; i--) {
      double v = b[i];
      for (int m = i + 1; m < k; m++) v -= H[m + (size_t) k * i] * b[m];
      b[i] = v / H[i + (size_t) k * i];
    }
  }
  return 1;
}


// The Hessian by central differences of the analytic gradient.
static void glmm_hess (glmm_t *m, const double *par, double *H, double *p2, double *gp, double *gm) {
  int k = m->P + 1;
  for (int c = 0; c < k; c++) {
    double h = 1e-5 * (1 + fabs(par[c]));
    memcpy(p2, par, (size_t) k * sizeof(double));
    p2[c] = par[c] + h; glmm_eval(m, p2, gp);
    p2[c] = par[c] - h; glmm_eval(m, p2, gm);
    for (int r = 0; r < k; r++) H[r + (size_t) k * c] = (gp[r] - gm[r]) / (2 * h);
  }
  for (int r = 0; r < k; r++)
    for (int c = 0; c < r; c++) {
      double v = 0.5 * (H[r + (size_t) k * c] + H[c + (size_t) k * r]);
      H[r + (size_t) k * c] = H[c + (size_t) k * r] = v;
    }
}


// Fits; writes coef, se (P) and *theta; returns a status.
static int m3_glmm (glmm_t *m, double *par, double *work, double *coef, double *se, double *theta) {

  int P = m->P, k = P + 1;
  double *g  = work;  work += k;
  double *H  = work;  work += (size_t) k * k;
  double *d  = work;  work += k;
  double *p2 = work;  work += k;
  double *gp = work;  work += k;
  double *gm = work;  work += k;
  double *pn = work;  work += k;

  for (int j = 0; j < m->J; j++) m->u[j] = 0;
  for (int c = 0; c < P; c++) { coef[c] = se[c] = R_NaN; }
  *theta = R_NaN;

  int conv = 0;
  double D = glmm_eval(m, par, g);
  if (!R_FINITE(D)) return M3_FAILED;

  for (int it = 0; it < 200; it++) {

    double gmax = 0;
    for (int c = 0; c < k; c++) if (fabs(g[c]) > gmax) gmax = fabs(g[c]);
    if (gmax < 1e-9) { conv = 1; break; }

    glmm_hess(m, par, H, p2, gp, gm);
    glmm_eval(m, par, g);   // restore the modes at par

    // Newton's direction, regularized until H is positive definite.
    double lam = 0, dmax = 0;
    for (int c = 0; c < k; c++) if (fabs(H[c + (size_t) k * c]) > dmax) dmax = fabs(H[c + (size_t) k * c]);
    for (int tries = 0; tries < 60; tries++) {
      double *A = work;
      memcpy(A, H, (size_t) k * k * sizeof(double));
      for (int c = 0; c < k; c++) { A[c + (size_t) k * c] += lam; d[c] = -g[c]; }
      if (m3_chol_solve(A, k, d)) break;
      lam = lam ? lam * 10 : 1e-8 * (dmax > 0 ? dmax : 1);
      if (tries == 59) return M3_NO_CONV;
    }

    // Backtracking line search.
    double slope = 0;
    for (int c = 0; c < k; c++) slope += g[c] * d[c];
    double t = 1, Dn = R_NaN;
    int ok = 0;
    for (int ls = 0; ls < 60; ls++) {
      for (int c = 0; c < k; c++) pn[c] = par[c] + t * d[c];
      pn[0] = fabs(pn[0]);
      Dn = glmm_eval(m, pn, NULL);
      if (R_FINITE(Dn) && Dn <= D + 1e-4 * t * slope) { ok = 1; break; }
      t /= 2;
    }
    double smax = 0;
    for (int c = 0; c < k; c++) { double v = fabs(t * d[c]) / (1 + fabs(par[c])); if (v > smax) smax = v; }
    if (!ok) {
      // No decrease along the direction: at the optimum to rounding.
      glmm_eval(m, par, g);
      conv = smax < 1e-6 || fabs(slope) < 1e-12 * (1 + fabs(D));
      break;
    }
    memcpy(par, pn, (size_t) k * sizeof(double));
    D = glmm_eval(m, par, g);
    if (smax < 1e-12) { conv = 1; break; }
  }

  for (int c = 1; c < k; c++) if (!(fabs(par[c]) < 1e3)) conv = 0;

  double *A = work;

  // A singular fit (theta < 1e-4, lme4's isSingular()): glmer computes no
  // Hessian, and its standard errors are those of the penalized least
  // squares at the conditional modes, the fixed effects' Schur complement
  // X'WX - sum_j theta^2 / (1 + theta^2 S_j) (X'W 1_j)(X'W 1_j)'.
  if (par[0] < 1e-4) {
    glmm_eval(m, par, g);
    int n = m->n, J = m->J;
    double th = par[0], W = 1 + 2 * m->a;
    for (int r = 0; r < P * P; r++) H[r] = 0;
    for (int j = 0; j < J; j++) {
      double S = 0;
      for (int c = 0; c < P; c++) p2[c] = 0;
      for (int i = m->gpos[j]; i < m->gpos[j + 1]; i++) {
        int s = m->gidx[i];
        double mu = lg_inv(m->xb[s] + th * m->u[j]);
        double w  = W * mu * (1 - mu);
        S += w;
        for (int c = 0; c < P; c++) {
          p2[c] += w * m->X[s + (size_t) n * c];
          for (int r = c; r < P; r++)
            H[r + (size_t) P * c] += w * m->X[s + (size_t) n * c] * m->X[s + (size_t) n * r];
        }
      }
      double f = th * th / (1 + th * th * S);
      for (int c = 0; c < P; c++) for (int r = c; r < P; r++) H[r + (size_t) P * c] -= f * p2[c] * p2[r];
    }
    for (int c = 0; c < P; c++) for (int r = 0; r < c; r++) H[r + (size_t) P * c] = H[c + (size_t) P * r];
    memcpy(A, H, (size_t) P * P * sizeof(double));
    if (m3_chol_solve(A, P, NULL)) {
      for (int c = 0; c < P; c++) {
        for (int r = 0; r < P; r++) d[r] = r == c ? 1 : 0;
        memcpy(A, H, (size_t) P * P * sizeof(double));
        m3_chol_solve(A, P, d);
        se[c] = sqrt(d[c]);
      }
    }
    for (int c = 1; c < k; c++) coef[c - 1] = par[c];
    *theta = par[0];
    return conv ? M3_OK : M3_NO_CONV;
  }

  // Otherwise the fixed-effects block of 2 H^-1.
  glmm_hess(m, par, H, p2, gp, gm);
  glmm_eval(m, par, g);
  memcpy(A, H, (size_t) k * k * sizeof(double));
  if (m3_chol_solve(A, k, NULL)) {
    for (int c = 1; c < k; c++) {
      for (int r = 0; r < k; r++) d[r] = r == c ? 1 : 0;
      memcpy(A, H, (size_t) k * k * sizeof(double));
      m3_chol_solve(A, k, d);
      se[c - 1] = sqrt(2 * d[c]);
    }
  }
  for (int c = 1; c < k; c++) coef[c - 1] = par[c];
  *theta = par[0];

  return conv ? M3_OK : M3_NO_CONV;
}



//======================================================
// Each feature's fit.
//======================================================
static void *m3_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  int n = m3_tab.n, D = m3_tab.D, P = m3_P, J = m3_J;
  double      *mem  = m3_mem  + (size_t) thread_i * m3_mem_sz;
  int         *imem = m3_imem + (size_t) thread_i * m3_imem_sz;
  long double *lmem = m3_lmem ? m3_lmem + (size_t) thread_i * ((size_t) J + 1) : NULL;

  double *y    = mem;  mem += n;
  double *yo   = mem;  mem += n;
  double *qty  = mem;  mem += n;
  double *coef = mem;  mem += P;
  double *se   = mem;  mem += P;
  double *df   = mem;  mem += P;
  double *vc   = mem;  mem += (size_t) P * P;
  double *Xr   = mem;  mem += (size_t) n * P;
  double *fmem = mem;  mem += flm_doubles(n, P, 0);
  int    *rows = imem; imem += n;
  int    *gr   = imem; imem += n;
  int    *gmap = imem; imem += J + 1;
  int    *keep = imem; imem += P;
  int    *piv  = imem; imem += P;

  flm_t q;
  flm_init(&q, fmem, piv, n, P, 0);

  // Prevalence, fixed effects.
  glm_work_t gw;
  if (m3_prev) {
    gw.xw  = mem; mem += (size_t) n * P;
    gw.zw  = mem; mem += n;
    gw.qty = mem; mem += n;
    gw.eta = mem; mem += 3 * (size_t) n;
    gw.mu  = mem; mem += 3 * (size_t) n;
    gw.b   = mem; mem += P;
    flm_init(&gw.q, mem, piv, n, P, 0);  mem += flm_doubles(n, P, 0);
  }

  // Prevalence, a random intercept.
  glmm_t gm;
  double *gpar = NULL, *gwork = NULL;
  int    *gidx = NULL, *gpos = NULL;
  if (m3_prev && m3_grp) {
    gidx = imem; imem += n;
    gpos = imem; imem += J + 1;
    for (int j = 0; j <= J; j++) gpos[j] = 0;
    for (int s = 0; s < n; s++) gpos[m3_grp[s] + 1]++;
    for (int j = 0; j < J; j++) gpos[j + 1] += gpos[j];
    for (int j = 0; j < J; j++) gmap[j] = gpos[j];
    for (int s = 0; s < n; s++) gidx[gmap[m3_grp[s]]++] = s;
    gm.X = m3_X; gm.gidx = gidx; gm.gpos = gpos; gm.n = n; gm.P = P; gm.J = J; gm.a = m3_a;
    gm.u  = mem; mem += J;
    gm.xb = mem; mem += n;
    gm.g2 = mem; mem += P + 1;
    gm.ub = mem; mem += P;
    gpar  = mem; mem += P + 1;
    gwork = mem; mem += 8 * ((size_t) P + 1) + 2 * ((size_t) P + 1) * (P + 1);
  }

  // Abundance, a random intercept.
  double *lxd = NULL, *lmm = NULL;
  if (!m3_prev && m3_grp) {
    lxd = mem; mem += lmm_design_doubles(n, P, J);
    lmm = mem; mem += lmm_doubles(n, P, J);
  }

  for (int f = thread_i; f < D; f += n_threads) {

    m3_column(f, y, m3_prev);

    int no = 0;
    for (int s = 0; s < n; s++) if (!ISNAN(y[s])) { rows[no] = s; yo[no++] = y[s]; }

    for (int j = 0; j < P; j++) {
      coef[j] = se[j] = df[j] = R_NaN;
      m3_xsd[f + (size_t) D * j] = m3_xsd_of(
        m3_X + (size_t) n * j, m3_prev ? NULL : rows, no, m3_prev && m3_a > 0 ? 3 : 1 );
    }
    double rsd = R_NaN, th = R_NaN;
    int status = M3_OK, warn = 0;

    // All values the same (none, one, or equal).
    int same = 1;
    for (int i = 1; i < no; i++) if (yo[i] != yo[0]) { same = 0; break; }
    if (same) status = M3_SAME;

    // A factor's reference level without a sample (maaslin3's check), or
    // a factor with a single level among the samples (lm()'s contrasts
    // fail).
    if (!status && !m3_prev)
      for (int k = 0; k < m3_K && !status; k++) {
        const int *lv = m3_base + (size_t) n * k;
        int any = 0, two = 0;
        for (int i = 0; i < no; i++) {
          if (lv[rows[i]] == 1) any = 1;
          if (lv[rows[i]] != lv[rows[0]]) two = 1;
        }
        if (!any) status = M3_BASELINE;
        else if (!two) status = M3_ONE_LEVEL;
      }

    if (!status && !m3_prev && !m3_grp) {
      double rss;
      flm_decompose(&q, m3_X, n, rows, no, P, 1e-7);
      flm_fit(&q, yo, qty, coef, &rss);
      flm_se(&q, yo, qty, rss, 0, se, NULL);
      double dof = no - q.rank;
      for (int j = 0; j < P; j++) df[j] = dof;
      if (dof > 0) rsd = sqrt(rss / dof);
    }

    else if (!status && !m3_prev) {
      // The groups among the feature's samples.
      for (int j = 0; j < J; j++) gmap[j] = -1;
      int Jo = 0;
      for (int i = 0; i < no; i++) {
        int g0 = m3_grp[rows[i]];
        if (gmap[g0] < 0) gmap[g0] = Jo++;
      }
      // Number the groups in their original order.
      for (int j = 0, c = 0; j < J; j++) if (gmap[j] >= 0) gmap[j] = c++;
      for (int i = 0; i < no; i++) gr[i] = gmap[m3_grp[rows[i]]];

      if (Jo < 2)        status = M3_ONE_GROUP;
      else if (Jo >= no) status = M3_MANY_GROUPS;
      else {
        // The design columns of full rank among these samples.
        flm_decompose(&q, m3_X, n, rows, no, P, 1e-7);
        int r = q.rank;
        for (int j = 0; j < P; j++) keep[j] = 0;
        for (int j = 0; j < r; j++) keep[q.pivot[j]] = 1;
        int c = 0;
        for (int j = 0; j < P; j++) if (keep[j]) {
          for (int i = 0; i < no; i++) Xr[i + (size_t) no * c] = m3_X[rows[i] + (size_t) n * j];
          c++;
        }
        if (no - r <= 0) status = M3_FAILED;
        else {
          lmm_design_t des;
          lmm_design_init_mem(&des, Xr, gr, no, r, Jo, lxd, lmem);
          double sigma, dev;
          double t = lmm_fit(&des, yo, -1, lmm, coef, vc, df, &sigma, &dev);
          if (ISNAN(t)) status = M3_FAILED;
          else {
            // Back to the design's columns.
            for (int j = P - 1, cc = r - 1; j >= 0; j--) {
              if (keep[j]) {
                double cf = coef[cc], v = vc[cc + (size_t) r * cc], dd = df[cc];
                coef[j] = cf; se[j] = sqrt(v); df[j] = dd; cc--;
              }
              else coef[j] = se[j] = df[j] = R_NaN;
            }
            rsd = sigma;
            th  = t;
          }
        }
      }
      if (status) for (int j = 0; j < P; j++) coef[j] = se[j] = df[j] = R_NaN;
    }

    else if (!status && !m3_grp) {
      warn = m3_glm(m3_X, n, P, y, m3_a, &gw, coef, se);
      for (int j = 0; j < P; j++) df[j] = R_PosInf;
    }

    else if (!status) {
      // Start at the fixed-effects fit, theta = 1.
      warn = m3_glm(m3_X, n, P, y, m3_a, &gw, coef, se);
      gpar[0] = 1;
      for (int j = 0; j < P; j++) gpar[j + 1] = ISNAN(coef[j]) || !R_FINITE(coef[j]) ? 0 : coef[j];
      warn = 0;
      gm.y = y;
      status = m3_glmm(&gm, gpar, gwork, coef, se, &th);
      if (status == M3_FAILED) for (int j = 0; j < P; j++) coef[j] = se[j] = R_NaN;
      for (int j = 0; j < P; j++) df[j] = R_PosInf;
    }

    for (int j = 0; j < P; j++) {
      m3_coef[f + (size_t) D * j] = coef[j];
      m3_se[f + (size_t) D * j]   = se[j];
      m3_df[f + (size_t) D * j]   = df[j];
    }
    m3_status[f] = status;
    m3_warn[f]   = warn;
    m3_nobs[f]   = no;
    m3_rsd[f]    = rsd;
    m3_theta[f]  = th;

    long double s = 0;
    for (int i = 0; i < no; i++) s += yo[i];
    m3_mean[f] = (no && !m3_prev) ? (double) (s / no) : (m3_prev ? (double) (s / n) : R_NaN);
  }

  return NULL;
}



//======================================================
// R interface. Fits every feature's abundance (model 0)
// or prevalence (model 1) model.
//
//   norm       0 none, 1 percent (TSS), 2 rclr (CLR).
//   log        1: log2 of the normalized values.
//   X          n x P design.
//   groups     1-based random intercept groups, or NULL.
//   base       n x K integer matrix: each sample's level
//              of factor k (1: the reference; abundance),
//              or NULL.
//   a          The augmentation's weight (0: none).
//
// Returns list(coef, se, df, xsd (D x P), status, warn,
// n_obs, mean, resid_sd, theta (D)). df: residual (lm),
// Satterthwaite (lmer) or Inf (logistic). xsd: the design
// column's standard deviation over the samples fitted,
// NaN where it is 0 / 1. mean: the mean log2 abundance
// where present, or the prevalence.
//======================================================
SEXP C_maaslin3_fit(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_log,
    SEXP sexp_X, SEXP sexp_model, SEXP sexp_groups, SEXP sexp_base,
    SEXP sexp_a, SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  m3_prepare(sexp_counts, sexp_margin, asInteger(sexp_norm), n_thr);
  m3_log  = asInteger(sexp_log);
  m3_prev = asInteger(sexp_model);
  m3_a    = asReal(sexp_a);

  int n = m3_tab.n, D = m3_tab.D;
  if (nrows(sexp_X) != n) error("The design must have one row per sample.");  // # nocov
  int P = ncols(sexp_X);
  m3_X = REAL(sexp_X);
  m3_P = P;

  int *grp = NULL, J = 0;
  if (!isNull(sexp_groups)) {
    grp = (int *) R_alloc((size_t) n + 1, sizeof(int));
    for (int s = 0; s < n; s++) { grp[s] = INTEGER(sexp_groups)[s] - 1; if (grp[s] + 1 > J) J = grp[s] + 1; }
  }
  m3_grp = grp;
  m3_J   = J;

  m3_base = isNull(sexp_base) ? NULL : INTEGER(sexp_base);
  m3_K    = isNull(sexp_base) ? 0 : ncols(sexp_base);

  int nt = threads_for(n_thr, (double) D * n * (P + 2) * (grp ? 40 : 4) * (m3_prev ? 10 : 1), 1e6);

  m3_mem_sz = 3 * (size_t) n + 3 * (size_t) P + (size_t) P * P + (size_t) n * P + flm_doubles(n, P, 0) + 16;
  if (m3_prev) m3_mem_sz += (size_t) n * P + 8 * (size_t) n + P + flm_doubles(n, P, 0);
  if (m3_prev && grp) m3_mem_sz += J + n + 3 * ((size_t) P + 1) + 8 * ((size_t) P + 1) + 2 * ((size_t) P + 1) * (P + 1);
  if (!m3_prev && grp) m3_mem_sz += lmm_design_doubles(n, P, J) + lmm_doubles(n, P, J);
  m3_imem_sz = 4 * (size_t) n + 2 * ((size_t) J + 1) + 2 * (size_t) P + 8;

  m3_mem  = (double *) R_alloc((size_t) nt * m3_mem_sz, sizeof(double));
  m3_imem = (int *)    R_alloc((size_t) nt * m3_imem_sz, sizeof(int));
  m3_lmem = grp ? (long double *) R_alloc((size_t) nt * ((size_t) J + 1), sizeof(long double)) : NULL;

  const char *nms[] = { "coef", "se", "df", "xsd", "status", "warn", "n_obs", "mean", "resid_sd", "theta", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_coef = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 0, s_coef);
  SEXP s_se   = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 1, s_se);
  SEXP s_df   = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 2, s_df);
  SEXP s_xsd  = allocMatrix(REALSXP, D, P); SET_VECTOR_ELT(result, 3, s_xsd);
  SEXP s_st   = allocVector(INTSXP, D);     SET_VECTOR_ELT(result, 4, s_st);
  SEXP s_wn   = allocVector(INTSXP, D);     SET_VECTOR_ELT(result, 5, s_wn);
  SEXP s_no   = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 6, s_no);
  SEXP s_mn   = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 7, s_mn);
  SEXP s_rsd  = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 8, s_rsd);
  SEXP s_th   = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 9, s_th);

  m3_coef = REAL(s_coef); m3_se = REAL(s_se); m3_df = REAL(s_df); m3_xsd = REAL(s_xsd);
  m3_status = INTEGER(s_st); m3_warn = INTEGER(s_wn);
  m3_nobs = REAL(s_no); m3_mean = REAL(s_mn); m3_rsd = REAL(s_rsd); m3_theta = REAL(s_th);

  run_parallel(m3_worker, nt, D < 100 ? 100 : D);

  UNPROTECT(1);
  return result;
}



//======================================================
// The median comparison's covariances.
//======================================================
#define M3_BLOCK 256

static const double *mc_coef, *mc_se, *mc_sims;
static const int    *mc_use;
static int           mc_nc, mc_nsims, mc_nblocks;
static uint64_t      mc_seed, mc_key;
static double        mc_center;
static double       *mc_med, *mc_bz, *mc_bzm, *mc_buf;


// The k-th smallest (0-based) of x[0..n-1], rearranging x.
static double m3_select (double *x, int n, int k) {
  int lo = 0, hi = n - 1;
  while (hi > lo) {
    int mid = lo + (hi - lo) / 2;
    // Median of three as the pivot.
    double a = x[lo], b = x[mid], c = x[hi];
    double pv = a < b ? (b < c ? b : (a < c ? c : a)) : (a < c ? a : (b < c ? c : b));
    int i = lo, j = hi;
    while (i <= j) {
      while (x[i] < pv) i++;
      while (x[j] > pv) j--;
      if (i <= j) { double t = x[i]; x[i] = x[j]; x[j] = t; i++; j--; }
    }
    if (k <= j) hi = j;
    else if (k >= i) lo = i;
    else return x[k];
  }
  return x[k];
}


static void *mc_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int nc = mc_nc;

  double *z   = mc_buf + (size_t) thread_i * 2 * ((size_t) nc + 1);
  double *buf = z + nc + 1;

  for (int b = thread_i; b < mc_nblocks; b += n_threads) {

    double *bz  = mc_bz  + (size_t) b * nc;
    double *bzm = mc_bzm + (size_t) b * nc;
    for (int j = 0; j < nc; j++) bz[j] = bzm[j] = 0;

    int s1 = (b + 1) * M3_BLOCK;
    if (s1 > mc_nsims) s1 = mc_nsims;

    for (int s = b * M3_BLOCK; s < s1; s++) {
      if (mc_sims) memcpy(z, mc_sims + (size_t) nc * s, (size_t) nc * sizeof(double));
      else {
        rng_t rng;
        rng_seed(&rng, mc_seed, (mc_key << 32) + (uint64_t) s);
        for (int j = 0; j < nc; j++) z[j] = rng_norm(&rng);
      }
      int nu = 0;
      for (int j = 0; j < nc; j++) if (mc_use[j]) buf[nu++] = mc_coef[j] + mc_se[j] * z[j];
      double med;
      if (nu == 0) med = R_NaN;
      else if (nu % 2) med = m3_select(buf, nu, nu / 2);
      else {
        double hi = m3_select(buf, nu, nu / 2);
        double lo = buf[0];
        for (int i = 1; i < nu / 2; i++) if (buf[i] > lo) lo = buf[i];
        med = (double) (((long double) lo + hi) / 2);
      }
      mc_med[s] = med;
      double dm = med - mc_center;
      for (int j = 0; j < nc; j++) { bz[j] += z[j]; bzm[j] += z[j] * dm; }
    }
  }
  return NULL;
}


//======================================================
// R interface: for features with coefficients `coef`
// and standard errors `se`, n_sims draws of them all
// (coef + se z, z standard normal: from PCG streams
// keyed by (seed, key, draw), or the columns of `sims`),
// the median of the `use` ones in each draw, and each
// feature's sample covariance with it. `center` (near
// the medians) is subtracted from them for accuracy.
// Returns list(cov (nc), medians (n_sims)).
//======================================================
SEXP C_maaslin3_median(
    SEXP sexp_coef, SEXP sexp_se, SEXP sexp_use, SEXP sexp_n_sims, SEXP sexp_seed,
    SEXP sexp_key, SEXP sexp_sims, SEXP sexp_center, SEXP sexp_n_threads ) {

  mc_nc     = length(sexp_coef);
  mc_nsims  = asInteger(sexp_n_sims);
  mc_coef   = REAL(sexp_coef);
  mc_se     = REAL(sexp_se);
  mc_use    = LOGICAL(sexp_use);
  mc_seed   = (uint64_t) (uint32_t) asInteger(sexp_seed);
  mc_key    = (uint64_t) asInteger(sexp_key);
  mc_sims   = isNull(sexp_sims) ? NULL : REAL(sexp_sims);
  mc_center = asReal(sexp_center);
  mc_nblocks = (mc_nsims + M3_BLOCK - 1) / M3_BLOCK;

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;
  int nt = threads_for(n_thr, (double) mc_nsims * mc_nc * 20, 1e6);
  if (nt > mc_nblocks) nt = mc_nblocks;
  if (nt < 1) nt = 1;

  const char *nms[] = { "cov", "medians", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_cov  = allocVector(REALSXP, mc_nc);    SET_VECTOR_ELT(result, 0, s_cov);
  SEXP s_med  = allocVector(REALSXP, mc_nsims); SET_VECTOR_ELT(result, 1, s_med);

  mc_med = REAL(s_med);
  mc_bz  = (double *) R_alloc((size_t) mc_nblocks * mc_nc + 1, sizeof(double));
  mc_bzm = (double *) R_alloc((size_t) mc_nblocks * mc_nc + 1, sizeof(double));
  mc_buf = (double *) R_alloc((size_t) nt * 2 * ((size_t) mc_nc + 1), sizeof(double));

  run_parallel(mc_worker, nt, mc_nblocks < 100 ? 100 : mc_nblocks);

  // Sums over the blocks, in order.
  long double sm = 0;
  for (int s = 0; s < mc_nsims; s++) sm += mc_med[s] - mc_center;
  double N = mc_nsims;
  for (int j = 0; j < mc_nc; j++) {
    long double sz = 0, szm = 0;
    for (int b = 0; b < mc_nblocks; b++) {
      sz  += mc_bz[(size_t) b * mc_nc + j];
      szm += mc_bzm[(size_t) b * mc_nc + j];
    }
    REAL(s_cov)[j] = mc_se[j] * (double) ((szm - sz * sm / N) / (N - 1));
  }

  UNPROTECT(1);
  return result;
}
