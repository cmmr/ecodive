// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// coda_lasso(): the elastic net under a zero-sum constraint, a C port of the
// coordinate descent of zeroSum 2.0.7 (src/zeroSumCoordinateDescent.cpp,
// zeroSumMoves.cpp, zeroSumApproximation.cpp, zeroSumCostFunction.cpp,
// zeroSumUtils.cpp and zeroSumCV.cpp), restricted to the gaussian and
// binomial families without fusion or rotated moves, and restructured for
// ecodive: each training set (each CV fold's, and the full data) runs its
// whole lambda path on its own thread, and pairs of zero coefficients are
// tested against the soft threshold before their full update (see
// zs_pair()). zeroSum is covered by this license:
//
//   MIT License
//
//   Copyright (c) 2019 Thorsten Rehberg
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
// The model, for samples i with weights w_i (1/n), features j and the
// zero-sum weights u_j (1 for taxa, 0 for covariates):
//
//   minimize   -loglik(b0, beta) + lambda * sum_j [ (1 - alpha)/2 (v_j beta_j)^2
//                                                   + alpha v_j |beta_j| ]
//   subject to sum_j u_j beta_j = 0,
//
// with -loglik = sum_i w_i (y_i - eta_i)^2 / 2 (gaussian) or
// -sum_i w_i (y_i eta_i - log(1 + exp(eta_i))) (binomial). As in zeroSum,
// the gaussian fit scales lambda by the (weighted) sd of y and v by its
// inverse, which leaves the lasso term alone and divides the ridge term by
// sd(y), as glmnet does.

#include "ecodive.h"
#include <float.h>  // DBL_EPSILON

#define ZS_GAUSSIAN 0
#define ZS_BINOMIAL 1
#define ZS_BETA_CHANGE 1e-14   // zeroSum's BETA_CHANGE_PRECISION


//======================================================
// Shared inputs (read only while the threads run).
//======================================================
typedef struct {
  int     n, P, K, L;          // samples, features, CV folds, lambdas
  int     family;
  int     zero_sum, polish, standardize, max_steps;
  double  alpha, precision;
  const double *x;             // n x P, centered by the full-data means
  const double *y;             // n
  const int    *foldid;        // n, 1..K
  const double *u, *v;         // P
  const double *lambdas;       // L
  double *beta_path;           // P x L (full-data fit)
  double *intercept_path;      // L     (full-data fit, centered x)
  double *fold_err;            // (K + 1) x L: -2 loglik, wCV weights
  double *cv_link;             // n x L: held-out linear predictors
  int    *n_done;              // K + 1: lambdas fitted before a NaN
  int     next;                // the next fit to claim
} zs_shared_t;

static zs_shared_t zs;


//======================================================
// One training set's state (zeroSum's per-fold slices).
//======================================================
typedef struct {
  int     f;                   // 0..K-1 a CV fold, K the full data
  double  lambda, ySD;
  double *w, *wOrg, *wCV;      // IRLS (or gaussian) weights; training; held-out
  double *yw;                  // the working response (gaussian: y)
  double *xb, *xbTmp;
  double *beta, *last_beta, *v;
  double  intercept, last_intercept;
  double  loglik, lasso, ridge, cost;
  int     approxFailed;
  int    *act; int n_act;      // the active set, sorted (zeroSum's std::set)
  double *G; uint64_t *Gstamp; uint64_t version;
  // Gaussian covariance updates (zs_pair_gram()): the weighted Gram matrix
  // of the features ever active, by slot.
  int     gram;                // in use
  int    *slot;                // P: each feature's slot, or -1
  int    *sfeat;               // capE: each slot's feature
  int     nE, capE;
  double *Q;                   // capE x capE: sum_i w_i x_ij x_ik
  double *m;                   // capE: sum_i w_i x_ij
  double *Gs;                  // capE: sum_i w_i x_ij r_i (active slots)
  double  R, sumW;             // sum_i w_i r_i; sum_i w_i
  int     gstale;
  int     lazy_xb;             // xb is left stale by the Gram moves (no covariates)
} zs_fit_t;


//======================================================
// The active set: a sorted array, iterated as std::set
// is (the next element above the current one, so that
// elements inserted during a loop are visited).
//======================================================
static int as_find (const zs_fit_t *F, int k) {
  int lo = 0, hi = F->n_act;
  while (lo < hi) { int m = (lo + hi) >> 1; if (F->act[m] < k) lo = m + 1; else hi = m; }
  return lo;
}
static int as_contains (const zs_fit_t *F, int k) {
  int i = as_find(F, k);
  return i < F->n_act && F->act[i] == k;
}
static int as_insert (zs_fit_t *F, int k) {
  int i = as_find(F, k);
  if (i < F->n_act && F->act[i] == k) return 0;
  memmove(F->act + i + 1, F->act + i, (size_t) (F->n_act - i) * sizeof(int));
  F->act[i] = k; F->n_act++;
  return 1;
}
// The first element above k (k = -1: the first), or -1.
static int as_next (const zs_fit_t *F, int k) {
  int lo = 0, hi = F->n_act;
  while (lo < hi) { int m = (lo + hi) >> 1; if (F->act[m] <= k) lo = m + 1; else hi = m; }
  return lo < F->n_act ? F->act[lo] : -1;
}
static void as_remove_zeros (zs_fit_t *F) {
  int m = 0;
  for (int i = 0; i < F->n_act; i++)
    if (F->beta[F->act[i]] != 0.0) F->act[m++] = F->act[i];
  F->n_act = m;
}


//======================================================
// Sums over the samples, in four interleaved partial
// sums (a fixed order, so results don't depend on the
// threads; zeroSum's AVX kernels also sum in lanes).
//======================================================

// sum_i w_i (a_i - b_i)^2
static double zs_wss (const double *restrict w, const double *restrict a,
                      const double *restrict b, int n) {
  double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    double d0 = a[i] - b[i], d1 = a[i+1] - b[i+1], d2 = a[i+2] - b[i+2], d3 = a[i+3] - b[i+3];
    s0 += w[i] * d0 * d0; s1 += w[i+1] * d1 * d1; s2 += w[i+2] * d2 * d2; s3 += w[i+3] * d3 * d3;
  }
  for (; i < n; i++) { double d = a[i] - b[i]; s0 += w[i] * d * d; }
  return (s0 + s1) + (s2 + s3);
}

// sum_i w_i (a_i - b_i - dk xk_i - ds xs_i)^2: the residual sum of squares
// after a pair move, without forming it.
static double zs_wss_move (const double *restrict w, const double *restrict a,
                           const double *restrict b, const double *restrict xk, double dk,
                           const double *restrict xs, double ds, int n) {
  double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    double d0 = a[i]   - ((b[i]   + dk * xk[i])   + ds * xs[i]);
    double d1 = a[i+1] - ((b[i+1] + dk * xk[i+1]) + ds * xs[i+1]);
    double d2 = a[i+2] - ((b[i+2] + dk * xk[i+2]) + ds * xs[i+2]);
    double d3 = a[i+3] - ((b[i+3] + dk * xk[i+3]) + ds * xs[i+3]);
    s0 += w[i] * d0 * d0; s1 += w[i+1] * d1 * d1; s2 += w[i+2] * d2 * d2; s3 += w[i+3] * d3 * d3;
  }
  for (; i < n; i++) { double d = a[i] - ((b[i] + dk * xk[i]) + ds * xs[i]); s0 += w[i] * d * d; }
  return (s0 + s1) + (s2 + s3);
}


//======================================================
// The likelihood and cost (zeroSumCostFunction.cpp).
//======================================================
static double zs_loglik (const zs_fit_t *F, const double *wf, const double *xb) {
  int n = zs.n;
  const double *y = zs.y;
  if (zs.family == ZS_GAUSSIAN) return -zs_wss(wf, y, xb, n) * 0.5;
  double t = 0.0;
  for (int i = 0; i < n; i++)
    t = fma(wf[i], fma(y[i], xb[i], -log1p(exp(xb[i]))), t);
  (void) F;
  return t;
}

static double zs_approx_loglik (const zs_fit_t *F, const double *xb) {
  return -0.5 * zs_wss(F->w, F->yw, xb, zs.n);
}

static void zs_cost (zs_fit_t *F) {
  int n = zs.n, P = zs.P;
  double alpha = zs.alpha;
  F->ridge = 0.0; F->lasso = 0.0;
  if (alpha != 1.0)
    for (int j = 0; j < P; j++) { double t = F->beta[j] * F->v[j]; F->ridge += pow(t, 2); }
  if (alpha != 0.0)
    for (int j = 0; j < P; j++) F->lasso += F->v[j] * fabs(F->beta[j]);

  double *xb = F->xb;
  for (int i = 0; i < n; i++) xb[i] = F->intercept;
  for (int j = 0; j < P; j++) {
    double b = F->beta[j];
    if (b == 0.0) continue;
    const double *xj = zs.x + (size_t) j * n;
    for (int i = 0; i < n; i++) xb[i] += b * xj[i];
  }
  F->version++;

  F->loglik = zs_loglik(F, F->wOrg, xb);
  F->cost   = -F->loglik + F->lambda * ((1.0 - alpha) * F->ridge / 2.0 + alpha * F->lasso);
}

static void zs_update_cost (zs_fit_t *F) {
  F->cost  += F->loglik;
  F->loglik = zs_approx_loglik(F, F->xb);   // only used with the approximation
  F->cost  -= F->loglik;
}


//======================================================
// The quadratic approximation of the binomial
// likelihood (zeroSumApproximation.cpp).
//======================================================
static void zs_refresh (zs_fit_t *F, int update_cost) {
  int n = zs.n;
  const double *restrict y = zs.y, *restrict wOrg = F->wOrg, *restrict xb = F->xb;
  double *restrict wl = F->w, *restrict yl = F->yw;
  int lowW = 0, nan = 0;
  for (int i = 0; i < n; i++) {
    // Samples outside the training set get weight 0 (their working
    // response is never used), without the exp().
    if (wOrg[i] == 0.0) { wl[i] = 0.0; yl[i] = xb[i]; lowW++; continue; }
    double e = exp(-fabs(xb[i]));
    double p = xb[i] >= 0.0 ? 1.0 / (e + 1.0) : e / (e + 1.0);
    double v = p * (1.0 - p);
    if (v < DBL_EPSILON * 100) { yl[i] = 0.0; wl[i] = 0.0; }
    else { yl[i] = xb[i] + (y[i] - p) / v; wl[i] = v * wOrg[i]; }
    if (wl[i] < 1e-6) lowW++;
    nan |= isnan(wl[i]);
  }
  F->approxFailed = nan || (double) lowW / (double) n > 0.7;
  F->version++;
  if (update_cost) zs_update_cost(F);
}

static inline int zs_approx (void) { return zs.family != ZS_GAUSSIAN; }


//======================================================
// Moves (zeroSumMoves.cpp).
//======================================================
static void zs_intercept_move (zs_fit_t *F) {
  if (F->approxFailed) return;
  int n = zs.n;
  const double *restrict yw = F->yw, *restrict w = F->w;
  double *restrict xb = F->xb;
  double b0 = 0, b1 = 0, b2 = 0, b3 = 0, w0 = 0, w1 = 0, w2 = 0, w3 = 0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    b0 += (yw[i]   - xb[i])   * w[i];   w0 += w[i];
    b1 += (yw[i+1] - xb[i+1]) * w[i+1]; w1 += w[i+1];
    b2 += (yw[i+2] - xb[i+2]) * w[i+2]; w2 += w[i+2];
    b3 += (yw[i+3] - xb[i+3]) * w[i+3]; w3 += w[i+3];
  }
  for (; i < n; i++) { b0 += (yw[i] - xb[i]) * w[i]; w0 += w[i]; }
  double bk = (b0 + b1) + (b2 + b3), sumW = (w0 + w1) + (w2 + w3);
  if (fabs(sumW) < 1000 * DBL_EPSILON) return;
  double old = F->intercept;
  F->intercept = bk / sumW + old;
  double diff = F->intercept - old;
  for (i = 0; i < n; i++) xb[i] += diff;
  F->version++;
}

// One coefficient outside the constraint (u_k = 0).
static int zs_cd_move (zs_fit_t *F, int k) {
  if (F->approxFailed) return 0;
  int n = zs.n;
  const double *xk = zs.x + (size_t) k * n;
  double *betak = F->beta + k;
  double ak = 0.0, bk = 0.0;
  for (int i = 0; i < n; i++) {
    double t = F->w[i] * xk[i];
    ak += t * xk[i];
    bk += t * (F->yw[i] - F->xb[i] + *betak * xk[i]);
  }
  double vk = F->v[k];
  ak += F->lambda * (1.0 - zs.alpha) * pow(vk, 2);
  double t  = F->lambda * zs.alpha * vk;
  double b1 = bk + t, b2 = bk - t;
  if (b1 < 0.0)      t = b1 / ak;
  else if (b2 > 0.0) t = b2 / ak;
  else               t = 0.0;
  double diff = t - *betak;
  *betak = t;
  if (fabs(diff) < ZS_BETA_CHANGE) return 0;
  for (int i = 0; i < n; i++) F->xb[i] += diff * xk[i];
  F->version++;
  if (zs_approx()) zs_refresh(F, 0);
  zs_intercept_move(F);
  return 1;
}

// G_j = sum_i w_i x_ij (yw_i - xb_i), computed when needed and kept until
// the residuals or weights change.
static double zs_grad (zs_fit_t *F, int j) {
  if (F->Gstamp[j] == F->version) return F->G[j];
  int n = zs.n;
  const double *xj = zs.x + (size_t) j * n;
  const double *restrict w = F->w, *restrict yw = F->yw, *restrict xb = F->xb;
  double g0 = 0, g1 = 0, g2 = 0, g3 = 0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    g0 += w[i]   * xj[i]   * (yw[i]   - xb[i]);
    g1 += w[i+1] * xj[i+1] * (yw[i+1] - xb[i+1]);
    g2 += w[i+2] * xj[i+2] * (yw[i+2] - xb[i+2]);
    g3 += w[i+3] * xj[i+3] * (yw[i+3] - xb[i+3]);
  }
  for (; i < n; i++) g0 += w[i] * xj[i] * (yw[i] - xb[i]);
  double g = (g0 + g1) + (g2 + g3);
  F->G[j] = g; F->Gstamp[j] = F->version;
  return g;
}

// A pair move: beta_k changes by d and beta_s by -(u_k/u_s) d, to the
// minimum of the 1-d problem (zeroSum tries the four sign patterns). When
// one of them is zero, whether the minimum is at d = 0 (no change; none of
// the four patterns is consistent then) is tested first from cached
// gradients, and the O(n) update is skipped if so.
static int zs_pair (zs_fit_t *F, int k, int s) {
  if (F->approxFailed) return 0;
  int n = zs.n;
  double *betak = F->beta + k, *betas = F->beta + s;
  const double *u = zs.u;
  double ukus = u[k] / u[s];
  double vk = F->v[k], vs = F->v[s];

  // The move leaves both unchanged when d = 0 minimizes the 1-d problem
  // f(d) (beta_k + d, beta_s - ukus d); when one of them is zero, that is a
  // subgradient condition on D = f'(0) of the smooth part, from the cached
  // gradients. (With both non-zero, d = 0 is optimal only by chance.)
  if (*betak == 0.0 || *betas == 0.0) {
    double D  = ukus * zs_grad(F, s) - zs_grad(F, k);
    double la = F->lambda * zs.alpha, lr = F->lambda * (1.0 - zs.alpha);
    if (*betak == 0.0 && *betas == 0.0) {
      if (!(fabs(D) > la * (vk + vs * ukus))) return 0;
    }
    else if (*betak == 0.0) {
      double c = D - lr * ukus * vs * vs * *betas - la * vs * ukus * (*betas > 0 ? 1.0 : -1.0);
      if (!(fabs(c) > la * vk)) return 0;
    }
    else {
      double c = D + lr * vk * vk * *betak + la * vk * (*betak > 0 ? 1.0 : -1.0);
      if (!(fabs(c) > la * vs * ukus)) return 0;
    }
  }

  const double *xk = zs.x + (size_t) k * n, *xs = zs.x + (size_t) s * n;
  double l1 = F->lambda * (1.0 - zs.alpha);
  double ak = l1 * (pow(vk, 2) + pow(vs, 2) * pow(ukus, 2));
  double bk = l1 * pow(vs, 2) * ukus * (*betas + ukus * *betak);
  {
    const double *restrict w = F->w, *restrict yw = F->yw, *restrict xb = F->xb;
    double bkk = *betak;
    double a0 = 0, a1 = 0, a2 = 0, a3 = 0, c0 = 0, c1 = 0, c2 = 0, c3 = 0;
    int i = 0;
    for (; i + 3 < n; i += 4) {
      double t0 = xs[i]   * ukus - xk[i],   t1 = xs[i+1] * ukus - xk[i+1];
      double t2 = xs[i+2] * ukus - xk[i+2], t3 = xs[i+3] * ukus - xk[i+3];
      a0 += w[i]   * t0 * t0; c0 += w[i]   * t0 * (yw[i]   - xb[i]   - bkk * t0);
      a1 += w[i+1] * t1 * t1; c1 += w[i+1] * t1 * (yw[i+1] - xb[i+1] - bkk * t1);
      a2 += w[i+2] * t2 * t2; c2 += w[i+2] * t2 * (yw[i+2] - xb[i+2] - bkk * t2);
      a3 += w[i+3] * t3 * t3; c3 += w[i+3] * t3 * (yw[i+3] - xb[i+3] - bkk * t3);
    }
    for (; i < n; i++) {
      double t = xs[i] * ukus - xk[i];
      a0 += w[i] * t * t; c0 += w[i] * t * (yw[i] - xb[i] - bkk * t);
    }
    ak += (a0 + a1) + (a2 + a3);
    bk -= (c0 + c1) + (c2 + c3);
  }

  double tmp3 = F->lambda * zs.alpha / ak;
  vs *= ukus;
  double case14 = tmp3 * (vk - vs), case23 = tmp3 * (vk + vs);
  bk /= ak;
  double bk1 = bk - case14, bk4 = bk + case14, bk2 = bk - case23, bk3 = bk + case23;
  tmp3 = u[s] * *betas + u[k] * *betak;
  double bs1 = tmp3 - u[k] * bk1, bs2 = tmp3 - u[k] * bk2;
  double bs3 = tmp3 - u[k] * bk3, bs4 = tmp3 - u[k] * bk4;

  double diffk = *betak, diffs = *betas;
  if      (bk1 > 0 && bs1 > 0) { *betak = bk1; *betas = bs1 / u[s]; }
  else if (bk2 > 0 && bs2 < 0) { *betak = bk2; *betas = bs2 / u[s]; }
  else if (bk3 < 0 && bs3 > 0) { *betak = bk3; *betas = bs3 / u[s]; }
  else if (bk4 < 0 && bs4 < 0) { *betak = bk4; *betas = bs4 / u[s]; }
  else return 0;

  diffk -= *betak;
  diffs -= *betas;
  for (int i = 0; i < n; i++) F->xb[i] -= xk[i] * diffk + xs[i] * diffs;
  F->version++;
  if (zs_approx()) zs_refresh(F, 0);
  zs_intercept_move(F);
  return 1;
}

// The polish's local-search move: beta_k += delta_k, beta_s compensates;
// kept if the cost falls.
static int zs_ls_move (zs_fit_t *F, int k, int s, double delta_k) {
  const double *u = zs.u;
  if (u[s] == 0) return 0;
  int n = zs.n;
  double betak = F->beta[k], betakNew = betak + delta_k;
  double betas = 0.0, betasNew = 0.0, tmp_lasso, tmp_ridge;
  double vk = F->v[k], vs = F->v[s];
  const double *xk = zs.x + (size_t) k * n, *xs = zs.x + (size_t) s * n;
  double delta_s = 0.0;

  if (zs.zero_sum && u[k] != 0) {
    delta_s = -delta_k * u[k] / u[s];
    betas = F->beta[s]; betasNew = betas + delta_s;
    tmp_lasso = F->lasso + vk * (fabs(betakNew) - fabs(betak)) +
                vs * (fabs(betasNew) - fabs(betas));
    tmp_ridge = F->ridge + vk * vk * (betakNew * betakNew - betak * betak) +
                vs * vs * (betasNew * betasNew - betas * betas);
  }
  else {
    xs = xk;    // delta_s = 0
    tmp_lasso = F->lasso + vk * (fabs(betakNew) - fabs(betak));
    tmp_ridge = F->ridge + vk * vk * (betakNew * betakNew - betak * betak);
  }

  // The moved fit's log-likelihood (the approximation for binomial).
  double tmp_ll = zs_approx()
    ? -0.5 * zs_wss_move(F->w, F->yw, F->xb, xk, delta_k, xs, delta_s, n)
    : -0.5 * zs_wss_move(F->wOrg, zs.y, F->xb, xk, delta_k, xs, delta_s, n);
  double tmp_cost = -tmp_ll + F->lambda * ((1.0 - zs.alpha) * tmp_ridge / 2.0 + zs.alpha * tmp_lasso);

  if (tmp_cost - F->cost < 0.0) {
    double *restrict xb = F->xb;
    for (int i = 0; i < n; i++) xb[i] = (xb[i] + delta_k * xk[i]) + delta_s * xs[i];
    F->ridge = tmp_ridge; F->lasso = tmp_lasso; F->loglik = tmp_ll; F->cost = tmp_cost;
    F->beta[k] = betakNew;
    if (zs.zero_sum && u[k] != 0) F->beta[s] = betasNew;
    F->version++;
    if (zs_approx()) zs_refresh(F, 1);
    return 1;
  }
  return 0;
}


//======================================================
// Gaussian covariance updates. With fixed weights, a
// pair's move needs only the Gram entries of the two
// features and their gradients G_j = sum_i w_i x_ij r_i;
// a successful move then updates the residuals once and
// the active features' gradients from the Gram matrix,
// instead of zeroSum's several passes over the samples.
// (The same moves in the same order; sums differ only in
// rounding.)
//======================================================
#define ZS_MAX_SLOTS 4096

// The slot of feature j, adding it (and its Gram row) if new; -1 when full.
static int zs_slot (zs_fit_t *F, int j) {
  if (F->slot[j] >= 0) return F->slot[j];
  int n = zs.n;
  if (F->nE == F->capE) {
    if (F->capE >= ZS_MAX_SLOTS) return -1;
    int cap = F->capE ? 2 * F->capE : 32;
    double *Q = (double *) calloc((size_t) cap * cap, sizeof(double));
    double *m = (double *) realloc(F->m, (size_t) cap * sizeof(double));
    if (m) F->m = m;
    double *G = (double *) realloc(F->Gs, (size_t) cap * sizeof(double));
    if (G) F->Gs = G;
    int *sf = (int *) realloc(F->sfeat, (size_t) cap * sizeof(int));
    if (sf) F->sfeat = sf;
    if (!Q || !m || !G || !sf) { free(Q); return -1; }
    for (int a = 0; a < F->nE; a++)
      memcpy(Q + (size_t) a * cap, F->Q + (size_t) a * F->capE, (size_t) F->nE * sizeof(double));
    free(F->Q);
    F->Q = Q; F->capE = cap;
  }
  int e = F->nE++;
  F->slot[j] = e;
  F->sfeat[e] = j;
  const double *restrict w = F->w, *restrict xj = zs.x + (size_t) j * n;
  double s0 = 0;
  for (int i = 0; i < n; i++) s0 += w[i] * xj[i];
  F->m[e] = s0;
  for (int a = 0; a <= e; a++) {
    const double *restrict xa = zs.x + (size_t) F->sfeat[a] * n;
    double q0 = 0, q1 = 0, q2 = 0, q3 = 0;
    int i = 0;
    for (; i + 3 < n; i += 4) {
      q0 += w[i] * xj[i] * xa[i];       q1 += w[i+1] * xj[i+1] * xa[i+1];
      q2 += w[i+2] * xj[i+2] * xa[i+2]; q3 += w[i+3] * xj[i+3] * xa[i+3];
    }
    for (; i < n; i++) q0 += w[i] * xj[i] * xa[i];
    double q = (q0 + q1) + (q2 + q3);
    F->Q[a + (size_t) e * F->capE] = F->Q[e + (size_t) a * F->capE] = q;
  }
  return e;
}

// The active features' gradients and the residuals' sum, from the samples.
// Returns 0 if a slot couldn't be made (then the direct moves are used).
static int zs_gram_prepare (zs_fit_t *F) {
  int n = zs.n;
  for (int a = 0; a < F->n_act; a++) if (zs_slot(F, F->act[a]) < 0) return 0;
  const double *restrict w = F->w, *restrict yw = F->yw, *restrict xb = F->xb;
  double r0 = 0, s0 = 0;
  for (int i = 0; i < n; i++) { r0 += w[i] * (yw[i] - xb[i]); s0 += w[i]; }
  F->R = r0; F->sumW = s0;
  for (int a = 0; a < F->n_act; a++) {
    int j = F->act[a];
    F->Gs[F->slot[j]] = zs_grad(F, j);
  }
  F->gstale = 0;
  return 1;
}

static int zs_pair_gram (zs_fit_t *F, int k, int s) {
  if (F->gstale && !zs_gram_prepare(F)) return zs_pair(F, k, s);   // # nocov
  int n = zs.n, cap = F->capE;
  int ek = F->slot[k], es = F->slot[s];
  double *betak = F->beta + k, *betas = F->beta + s;
  const double *u = zs.u;
  double ukus = u[k] / u[s];
  double vk = F->v[k], vs = F->v[s];
  double l1 = F->lambda * (1.0 - zs.alpha);
  double T2 = ukus * ukus * F->Q[es + (size_t) es * cap] - 2.0 * ukus * F->Q[es + (size_t) ek * cap]
              + F->Q[ek + (size_t) ek * cap];
  double ak = l1 * (pow(vk, 2) + pow(vs, 2) * pow(ukus, 2)) + T2;
  double bk = l1 * pow(vs, 2) * ukus * (*betas + ukus * *betak)
              - ((ukus * F->Gs[es] - F->Gs[ek]) - *betak * T2);

  double tmp3 = F->lambda * zs.alpha / ak;
  double vsu = vs * ukus;
  double case14 = tmp3 * (vk - vsu), case23 = tmp3 * (vk + vsu);
  bk /= ak;
  double bk1 = bk - case14, bk4 = bk + case14, bk2 = bk - case23, bk3 = bk + case23;
  tmp3 = u[s] * *betas + u[k] * *betak;
  double bs1 = tmp3 - u[k] * bk1, bs2 = tmp3 - u[k] * bk2;
  double bs3 = tmp3 - u[k] * bk3, bs4 = tmp3 - u[k] * bk4;

  double oldk = *betak, olds = *betas;
  if      (bk1 > 0 && bs1 > 0) { *betak = bk1; *betas = bs1 / u[s]; }
  else if (bk2 > 0 && bs2 < 0) { *betak = bk2; *betas = bs2 / u[s]; }
  else if (bk3 < 0 && bs3 > 0) { *betak = bk3; *betas = bs3 / u[s]; }
  else if (bk4 < 0 && bs4 < 0) { *betak = bk4; *betas = bs4 / u[s]; }
  else return 0;

  double dk = *betak - oldk, ds = *betas - olds;
  // The intercept's move (zeroSum's interceptMove()), folded in.
  double R  = F->R - (dk * F->m[ek] + ds * F->m[es]);
  double d0 = fabs(F->sumW) < 1000 * DBL_EPSILON ? 0.0 : R / F->sumW;
  F->intercept += d0;
  F->R = R - d0 * F->sumW;
  // xb is recomputed by zs_cost() after the sweep; only covariates' moves
  // (zs_cd_move()) read it within one.
  if (!F->lazy_xb) {
    const double *restrict xk = zs.x + (size_t) k * n, *restrict xs = zs.x + (size_t) s * n;
    double *restrict xb = F->xb;
    for (int i = 0; i < n; i++) xb[i] += (xk[i] * dk + xs[i] * ds) + d0;
  }
  for (int a = 0; a < F->n_act; a++) {
    int e = F->slot[F->act[a]];
    F->Gs[e] -= dk * F->Q[e + (size_t) ek * cap] + ds * F->Q[e + (size_t) es * cap] + d0 * F->m[e];
  }
  F->version++;
  return 1;
}


// The polish's move (zs_ls_move()), gaussian, from the Gram matrix: the
// residual sum of squares after moving beta_k by dk and beta_s by ds is
// RSS - 2 (dk G_k + ds G_s) + dk^2 Q_kk + 2 dk ds Q_ks + ds^2 Q_ss.
static int zs_ls_move_gram (zs_fit_t *F, int k, int s, double delta_k) {
  const double *u = zs.u;
  if (u[s] == 0) return 0;  // # nocov
  int cap = F->capE, ek = F->slot[k], es = F->slot[s];
  double betak = F->beta[k], betakNew = betak + delta_k;
  double delta_s = -delta_k * u[k] / u[s];
  double betas = F->beta[s], betasNew = betas + delta_s;
  double vk = F->v[k], vs = F->v[s];
  double tmp_lasso = F->lasso + vk * (fabs(betakNew) - fabs(betak)) +
                     vs * (fabs(betasNew) - fabs(betas));
  double tmp_ridge = F->ridge + vk * vk * (betakNew * betakNew - betak * betak) +
                     vs * vs * (betasNew * betasNew - betas * betas);
  double rss = -2.0 * F->loglik
    - 2.0 * (delta_k * F->Gs[ek] + delta_s * F->Gs[es])
    + delta_k * delta_k * F->Q[ek + (size_t) ek * cap]
    + 2.0 * delta_k * delta_s * F->Q[ek + (size_t) es * cap]
    + delta_s * delta_s * F->Q[es + (size_t) es * cap];
  double tmp_ll = -0.5 * rss;
  double tmp_cost = -tmp_ll + F->lambda * ((1.0 - zs.alpha) * tmp_ridge / 2.0 + zs.alpha * tmp_lasso);
  if (!(tmp_cost - F->cost < 0.0)) return 0;

  F->ridge = tmp_ridge; F->lasso = tmp_lasso; F->loglik = tmp_ll; F->cost = tmp_cost;
  F->beta[k] = betakNew; F->beta[s] = betasNew;
  F->R -= delta_k * F->m[ek] + delta_s * F->m[es];
  for (int a = 0; a < F->n_act; a++) {
    int e = F->slot[F->act[a]];
    F->Gs[e] -= delta_k * F->Q[e + (size_t) ek * cap] + delta_s * F->Q[e + (size_t) es * cap];
  }
  if (!F->lazy_xb) {
    int n = zs.n;
    const double *restrict xk = zs.x + (size_t) k * n, *restrict xs = zs.x + (size_t) s * n;
    double *restrict xb = F->xb;
    for (int i = 0; i < n; i++) xb[i] = (xb[i] + delta_k * xk[i]) + delta_s * xs[i];
  }
  F->version++;
  return 1;
}


//======================================================
// One lambda (zeroSum::doFitUsingCoordinateDescent).
//======================================================
static void zs_save (zs_fit_t *F) {
  memcpy(F->last_beta, F->beta, (size_t) zs.P * sizeof(double));
  F->last_intercept = F->intercept;
}
static void zs_restore (zs_fit_t *F) {
  memcpy(F->beta, F->last_beta, (size_t) zs.P * sizeof(double));
  F->intercept = F->last_intercept;
}

static void zs_fit_lambda (zs_fit_t *F) {
  int P = zs.P;
  const double *u = zs.u;
  zs_cost(F);
  F->n_act = 0;

  int gettingWorse = 0;
  for (int steps = 0; steps < zs.max_steps && !gettingWorse; steps++) {
    int asc = 0;
    double costBefore = F->cost;
    zs_save(F);

    F->approxFailed = 0;
    if (zs_approx()) zs_refresh(F, 0);
    zs_intercept_move(F);

    if (zs.zero_sum) {
      for (int s = 0; s < P; s++) {
        if (u[s] == 0) {
          if (zs_cd_move(F, s)) asc += as_insert(F, s);
          continue;
        }
        if (steps > 0) {
          if (as_contains(F, s)) continue;
          for (int k = as_next(F, -1); k >= 0; k = as_next(F, k)) {
            if (u[k] == 0) continue;
            if (zs_pair(F, s, k)) { asc += as_insert(F, k); asc += as_insert(F, s); }
          }
        }
        else {
          for (int k = s + 1; k < P; k++) {
            if (u[k] == 0) continue;
            if (zs_pair(F, s, k)) { asc += as_insert(F, k); asc += as_insert(F, s); }
          }
        }
      }
    }
    else {
      for (int s = 0; s < P; s++)
        if (zs_cd_move(F, s)) asc += as_insert(F, s);
    }

    if (zs.polish) {
      zs_cost(F);
      if (zs_approx()) zs_refresh(F, 1);
      int gram = F->gram && zs.zero_sum && zs_gram_prepare(F);
      for (int s = as_next(F, -1); s >= 0; s = as_next(F, s)) {
        if (fabs(F->beta[s]) < DBL_EPSILON * 100 || u[s] == 0) continue;
        for (int k = as_next(F, -1); k >= 0; k = as_next(F, k)) {
          if (s == k || fabs(F->beta[k]) < DBL_EPSILON * 100 || u[k] == 0) continue;
          if (gram) zs_ls_move_gram(F, k, s, -F->beta[k]);
          else      zs_ls_move(F, k, s, -F->beta[k]);
        }
      }
      as_remove_zeros(F);
    }

    zs_cost(F);
    double costAfter = F->cost;
    if (isnan(costAfter) || isinf(costAfter) || costBefore < costAfter) {
      zs_restore(F); zs_cost(F);
      break;
    }
    if (asc == 0 || F->n_act == 0) break;

    for (;;) {
      costBefore = F->cost;
      zs_save(F);
      int success = 0;
      F->approxFailed = 0;
      if (zs_approx()) zs_refresh(F, 0);
      zs_intercept_move(F);

      if (zs.zero_sum) {
        int gram = F->gram && zs_gram_prepare(F);
        for (int s = as_next(F, -1); s >= 0; s = as_next(F, s)) {
          if (u[s] == 0) {
            int c = zs_cd_move(F, s);
            success += c;
            if (c) F->gstale = 1;
            continue;
          }
          for (int k = as_next(F, -1); k >= 0; k = as_next(F, k)) {
            if (s == k || u[k] == 0) continue;
            success += gram ? zs_pair_gram(F, s, k) : zs_pair(F, s, k);
          }
        }
      }
      else {
        for (int s = as_next(F, -1); s >= 0; s = as_next(F, s))
          success += zs_cd_move(F, s);
      }

      zs_cost(F);
      costAfter = F->cost;
      if (isnan(costAfter) || costBefore < costAfter) {
        gettingWorse = 1;
        zs_restore(F); zs_cost(F);
        break;
      }
      if (success == 0 || fabs(costAfter - costBefore) < zs.precision) break;
    }
  }

  for (int j = 0; j < P; j++)
    if (fabs(F->beta[j]) < 100 * DBL_EPSILON) F->beta[j] = 0.0;
}


//======================================================
// One training set's whole path, with zeroSum's weights
// (adjustWeights, standardizeResponse, standardizeData).
//======================================================
static int zs_path (int f) {

  int n = zs.n, P = zs.P, K = zs.K, L = zs.L;
  zs_fit_t F;
  memset(&F, 0, sizeof(F));
  F.f = f;

  double *mem = (double *) calloc((size_t) 7 * n + (size_t) 4 * P, sizeof(double));
  int    *act = (int *) malloc(((size_t) P + 1) * sizeof(int));
  uint64_t *stamp = (uint64_t *) malloc(((size_t) P + 1) * sizeof(uint64_t));
  if (!mem || !act || !stamp) { free(mem); free(act); free(stamp); return 1; }
  F.w = mem; F.wOrg = F.w + n; F.wCV = F.wOrg + n; F.yw = F.wCV + n;
  F.xb = F.yw + n; F.xbTmp = F.xb + n;
  F.beta = F.xbTmp + 2 * (size_t) n; F.last_beta = F.beta + P; F.v = F.last_beta + P;
  F.G = F.v + P;
  F.act = act; F.Gstamp = stamp;
  for (int j = 0; j < P; j++) stamp[j] = 0;
  F.version = 1;
  F.gram = zs.family == ZS_GAUSSIAN && zs.zero_sum;
  F.slot = F.gram ? (int *) malloc((size_t) P * sizeof(int)) : NULL;
  if (F.gram && !F.slot) F.gram = 0;   // # nocov
  if (F.slot) for (int j = 0; j < P; j++) F.slot[j] = -1;
  F.lazy_xb = F.gram;
  for (int j = 0; j < P; j++) if (zs.u[j] == 0) F.lazy_xb = 0;

  // Weights: 1/n; a CV fold rescales its training samples to 1/(n - fold
  // size) and its held-out ones (wCV) to 1/(fold size).
  for (int i = 0; i < n; i++) F.w[i] = F.wOrg[i] = F.wCV[i] = 1.0 / n;
  if (f < K) {
    int fs = 0;
    for (int i = 0; i < n; i++) if (zs.foldid[i] - 1 == f) fs++;
    double sc1 = (double) n / (double) (n - fs), sc2 = (double) n / (double) fs;
    for (int i = 0; i < n; i++) {
      if (zs.foldid[i] - 1 == f) { F.w[i] = 0.0; F.wOrg[i] = 0.0; F.wCV[i] *= sc2; }
      else { F.w[i] *= sc1; F.wOrg[i] *= sc1; F.wCV[i] *= 0.0; }
    }
  }
  for (int i = 0; i < n; i++) F.yw[i] = zs.y[i];
  for (int j = 0; j < P; j++) F.v[j] = zs.v[j];

  F.ySD = 1.0;
  if (zs.family == ZS_GAUSSIAN) {
    double sw = 0.0;
    for (int i = 0; i < n; i++) sw += F.wOrg[i];
    sw = 1.0 / sw;
    double ym = 0.0;
    for (int i = 0; i < n; i++) ym += F.wOrg[i] * zs.y[i];
    ym *= sw;
    double sd = 0.0;
    for (int i = 0; i < n; i++) sd += F.wOrg[i] * pow(zs.y[i] - ym, 2);
    F.ySD = sqrt(sd);
    for (int j = 0; j < P; j++) F.v[j] /= F.ySD;
  }
  if (zs.standardize) {
    double sw = 0.0;
    for (int i = 0; i < n; i++) sw += F.wOrg[i];
    sw = 1.0 / sw;
    for (int j = 0; j < P; j++) {
      const double *xj = zs.x + (size_t) j * n;
      double m = 0.0;
      for (int i = 0; i < n; i++) m += xj[i] * F.wOrg[i];
      m *= sw;
      double t = 0.0;
      for (int i = 0; i < n; i++) t += F.wOrg[i] * pow(xj[i] - m, 2);
      F.v[j] *= sqrt(t * sw);
    }
  }

  int done = L;
  for (int l = 0; l < L; l++) {
    F.lambda = zs.lambdas[l] * F.ySD;
    zs_fit_lambda(&F);

    // The error with the held-out weights (zeroSum's cv_tmp).
    zs_cost(&F);
    double err = -2.0 * zs_loglik(&F, F.wCV, F.xb);
    zs.fold_err[f + (size_t) l * (K + 1)] = err;

    if (f < K) {
      for (int i = 0; i < n; i++)
        if (zs.foldid[i] - 1 == f) zs.cv_link[i + (size_t) l * n] = F.xb[i];
    }
    else {
      memcpy(zs.beta_path + (size_t) l * P, F.beta, (size_t) P * sizeof(double));
      zs.intercept_path[l] = F.intercept;
    }
    if (isnan(err) || isinf(err)) { done = l; break; }
  }
  zs.n_done[f] = done;

  free(mem); free(act); free(stamp);
  free(F.slot); free(F.sfeat); free(F.Q); free(F.m); free(F.Gs);
  return 0;
}


static void *zs_worker (void *arg) {
  (void) arg;
  for (;;) {
    int f = __atomic_fetch_add(&zs.next, 1, __ATOMIC_RELAXED);
    if (f > zs.K) break;
    if (zs_path(f)) zs.n_done[f] = -1;
  }
  return NULL;
}


//======================================================
// R interface.
//
// x        n x P (log abundances, then covariates)
// y        n (gaussian), or 0/1 (binomial)
// foldid   n, 1..K (ignored when K = 0)
// u, v     P: zero-sum weights and penalty factors
// lambdas  decreasing
//
// Returns list(beta (P x L), intercept (L), cv_error (L),
// cv_sd (L), train_error (L), fold_error ((K+1) x L),
// cv_link (n x L), n_fit), on the uncentered x.
//======================================================
SEXP C_coda_lasso (
    SEXP sexp_x, SEXP sexp_y, SEXP sexp_foldid, SEXP sexp_K, SEXP sexp_family,
    SEXP sexp_alpha, SEXP sexp_lambdas, SEXP sexp_u, SEXP sexp_v,
    SEXP sexp_standardize, SEXP sexp_precision, SEXP sexp_max_steps,
    SEXP sexp_polish, SEXP sexp_zero_sum, SEXP sexp_n_threads ) {

  int n = nrows(sexp_x), P = ncols(sexp_x);
  int K = asInteger(sexp_K), L = length(sexp_lambdas);
  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  // Center x by the full-data means (zeroSum's meanCentering).
  double *x = (double *) R_alloc((size_t) n * P, sizeof(double));
  double *mean = (double *) R_alloc((size_t) P, sizeof(double));
  const double *x0 = REAL(sexp_x);
  double w0 = 1.0 / n, sw0 = 0.0;
  for (int i = 0; i < n; i++) sw0 += w0;
  sw0 = 1.0 / sw0;
  for (int j = 0; j < P; j++) {
    const double *xj = x0 + (size_t) j * n;
    double m = 0.0;
    for (int i = 0; i < n; i++) m += xj[i] * w0;
    m *= sw0;
    mean[j] = m;
    for (int i = 0; i < n; i++) x[i + (size_t) j * n] = xj[i] - m;
  }

  SEXP sexp_beta  = PROTECT(allocMatrix(REALSXP, P, L));
  SEXP sexp_icpt  = PROTECT(allocVector(REALSXP, L));
  SEXP sexp_cve   = PROTECT(allocVector(REALSXP, L));
  SEXP sexp_cvsd  = PROTECT(allocVector(REALSXP, L));
  SEXP sexp_train = PROTECT(allocVector(REALSXP, L));
  SEXP sexp_ferr  = PROTECT(allocMatrix(REALSXP, K + 1, L));
  SEXP sexp_link  = PROTECT(allocMatrix(REALSXP, n, L));
  int *n_done = (int *) R_alloc((size_t) K + 1, sizeof(int));

  memset(REAL(sexp_beta), 0, (size_t) P * L * sizeof(double));
  for (size_t i = 0; i < (size_t) n * L; i++) REAL(sexp_link)[i] = NA_REAL;
  for (size_t i = 0; i < (size_t) (K + 1) * L; i++) REAL(sexp_ferr)[i] = NA_REAL;
  for (int l = 0; l < L; l++) REAL(sexp_icpt)[l] = NA_REAL;

  zs.n = n; zs.P = P; zs.K = K; zs.L = L;
  zs.family      = asInteger(sexp_family);
  zs.alpha       = asReal(sexp_alpha);
  zs.precision   = asReal(sexp_precision);
  zs.max_steps   = asInteger(sexp_max_steps);
  zs.polish      = asInteger(sexp_polish);
  zs.zero_sum    = asInteger(sexp_zero_sum);
  zs.standardize = asInteger(sexp_standardize);
  zs.x = x; zs.y = REAL(sexp_y); zs.foldid = INTEGER(sexp_foldid);
  zs.u = REAL(sexp_u); zs.v = REAL(sexp_v); zs.lambdas = REAL(sexp_lambdas);
  zs.beta_path = REAL(sexp_beta); zs.intercept_path = REAL(sexp_icpt);
  zs.fold_err = REAL(sexp_ferr); zs.cv_link = REAL(sexp_link);
  zs.n_done = n_done;
  zs.next = 0;

  int nt = n_threads < K + 1 ? n_threads : K + 1;
  run_parallel(zs_worker, nt, nt > 1 ? 100 : 1);

  int n_fit = L;
  for (int f = 0; f <= K; f++) {
    if (n_done[f] < 0) error("Memory allocation failed.");  // # nocov
    if (n_done[f] < n_fit) n_fit = n_done[f];
  }

  // CV error and its SD, as zeroSum (in fold order).
  double *ferr = REAL(sexp_ferr);
  for (int l = 0; l < L; l++) {
    const double *e = ferr + (size_t) l * (K + 1);
    REAL(sexp_train)[l] = e[K];
    if (K == 0) { REAL(sexp_cve)[l] = NA_REAL; REAL(sexp_cvsd)[l] = NA_REAL; continue; }
    double cv = 0.0;
    for (int f = 0; f < K; f++) cv += e[f];
    cv /= K;
    double sd = 0.0;
    for (int f = 0; f < K; f++) sd += e[f] * e[f];
    sd /= (double) K;
    sd -= cv * cv;
    sd /= (double) K - 1.0;
    sd = sqrt(sd);
    REAL(sexp_cve)[l] = cv; REAL(sexp_cvsd)[l] = sd;
  }

  // Back to the uncentered x: the intercept absorbs the means.
  for (int l = 0; l < L; l++) {
    const double *b = REAL(sexp_beta) + (size_t) l * P;
    double t = REAL(sexp_icpt)[l];
    for (int j = 0; j < P; j++) t -= b[j] * mean[j];
    REAL(sexp_icpt)[l] = t;
  }

  SEXP res = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(res, 0, sexp_beta);
  SET_VECTOR_ELT(res, 1, sexp_icpt);
  SET_VECTOR_ELT(res, 2, sexp_cve);
  SET_VECTOR_ELT(res, 3, sexp_cvsd);
  SET_VECTOR_ELT(res, 4, sexp_train);
  SET_VECTOR_ELT(res, 5, sexp_ferr);
  SET_VECTOR_ELT(res, 6, sexp_link);
  SET_VECTOR_ELT(res, 7, ScalarInteger(n_fit));
  SEXP nms = PROTECT(allocVector(STRSXP, 8));
  const char *nm[] = { "beta", "intercept", "cv_error", "cv_sd", "train_error",
                       "fold_error", "cv_link", "n_fit" };
  for (int i = 0; i < 8; i++) SET_STRING_ELT(nms, i, mkChar(nm[i]));
  setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(9);
  return res;
}
