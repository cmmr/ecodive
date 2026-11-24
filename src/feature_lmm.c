// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Linear mixed models with one random intercept, fitted by REML, for one
// design and many responses (features): for linda(), dirmult_lm() and
// maaslin3().
//
//   y = X b + Z u + e,  u ~ N(0, s^2 t^2 I_J),  e ~ N(0, s^2 I_n),
//
// with Z the indicators of J groups (sizes n_j) and t the ratio of the
// random intercept's standard deviation to the residual's (lme4's theta).
// Written from Bates, Maechler, Bolker & Walker (2015, J. Stat. Softw.
// 67(1), the profiled REML criterion) and Kuznetsova, Brockhoff &
// Christensen (2017, J. Stat. Softw. 82(13), Satterthwaite's degrees of
// freedom from the variance parameters (t, s)); no lme4 or lmerTest code
// was read.
//
// With H = I + t^2 Z Z', block j of H^-1 is I - (t^2 / (1 + n_j t^2)) 1 1',
// so with c_j = n_j / (1 + n_j t^2), group means xbar_j, ybar_j and
// within-group centred Xc, yc:
//
//   M(t)  = X' H^-1 X = Xc'Xc + sum_j c_j xbar_j xbar_j'
//   X'H^-1 y          = Xc'yc + sum_j c_j ybar_j xbar_j
//   r2(t) = (y - X b)' H^-1 (y - X b)
//         = |yc - Xc b|^2 + sum_j c_j (ybar_j - xbar_j' b)^2
//
// (no subtraction of nearly equal terms), b = M^-1 X'H^-1 y, and the
// profiled REML deviance is
//
//   d(t) = sum_j log(1 + n_j t^2) + log|M| + (n - p)(1 + log(2 pi r2 / (n - p))).
//
// d is even in t, so its derivative is t h(t) with h analytic: t = 0 is a
// minimum when h(0) >= 0 and nothing smaller lies below the first root.
// The optimum is found by Brent's minimisation on [0, t_hi], then made
// exact by bracketing the root of h. Then s^2 = r2 / (n - p), Var(b) =
// s^2 M^-1, and for each coefficient k, Satterthwaite's
//
//   df_k = 2 V_kk^2 / (g' A g),  g = dV_kk / d(t, s),  A = 2 Hess^-1,
//
// with Hess the Hessian of the (unprofiled) REML deviance in (t, s), all
// derivatives analytic.

#include "ecodive.h"


//======================================================
// The design, once for all responses.
//======================================================
void lmm_design_init(lmm_design_t *d, const double *X, const int *grp, int n, int p, int J) {
  double      *mem  = (double *) R_alloc(lmm_design_doubles(n, p, J), sizeof(double));
  long double *lmem = (long double *) R_alloc((size_t) J + 1, sizeof(long double));
  lmm_design_init_mem(d, X, grp, n, p, J, mem, lmem);
}


// The same, in the caller's memory (lmm_design_doubles(n, p, J) doubles and
// J long doubles), so that worker threads can set up a design per response
// (maaslin3(), whose responses are observed in different samples).
size_t lmm_design_doubles(int n, int p, int J) {
  return (size_t) J + (size_t) J * p + (size_t) n * p + (size_t) p * p + 4;
}


void lmm_design_init_mem(
    lmm_design_t *d, const double *X, const int *grp, int n, int p, int J,
    double *mem, long double *acc ) {

  d->n = n; d->p = p; d->J = J; d->grp = grp;

  double *nj   = mem;  mem += J;
  double *xbar = mem;  mem += (size_t) J * p;   // J x p
  double *xc   = mem;  mem += (size_t) n * p;   // n x p
  double *wxx  = mem;

  for (int j = 0; j < J; j++) nj[j] = 0;
  for (int i = 0; i < n; i++) nj[grp[i]]++;

  for (int k = 0; k < p; k++) {
    const double *x = X + (size_t) n * k;
    for (int j = 0; j < J; j++) xbar[j + (size_t) J * k] = 0;
    for (int j = 0; j < J; j++) acc[j] = 0;
    for (int i = 0; i < n; i++) acc[grp[i]] += x[i];
    for (int j = 0; j < J; j++) xbar[j + (size_t) J * k] = (double) (acc[j] / nj[j]);
    for (int i = 0; i < n; i++) xc[i + (size_t) n * k] = x[i] - xbar[grp[i] + (size_t) J * k];
  }

  for (int a = 0; a < p; a++)
    for (int b = 0; b <= a; b++) {
      long double s = 0;
      const double *u = xc + (size_t) n * a, *v = xc + (size_t) n * b;
      for (int i = 0; i < n; i++) s += (long double) u[i] * v[i];
      wxx[a + (size_t) p * b] = wxx[b + (size_t) p * a] = (double) s;
    }

  d->nj = nj; d->xbar = xbar; d->xc = xc; d->wxx = wxx;
}


size_t lmm_doubles(int n, int p, int J) {
  return (size_t) n + (size_t) J * 3 + (size_t) p * p * 8 + (size_t) p * 8 + 16;
}



//======================================================
// One response.
//======================================================
typedef struct {
  const lmm_design_t *d;
  double *yc;      // n: within-group centred y
  double *ybar;    // J
  double *wxy;     // p
  double *e;       // J: ybar_j - xbar_j' b
  double *M;       // p x p, then its Cholesky factor (lower)
  double *Mi;      // p x p inverse
  double *b;       // p
  double *t1, *t2, *t3, *t4;   // p x p scratch
  double *v;       // p scratch
  double  r2, logdetM;
  int     ok;
} lmm_work_t;


// Cholesky (lower, in place) of a p x p matrix; log determinant. 0 if not
// positive definite.
static int chol (double *A, int p, double *logdet) {
  long double ld = 0;
  for (int j = 0; j < p; j++) {
    long double s = A[j + (size_t) p * j];
    for (int k = 0; k < j; k++) s -= (long double) A[j + (size_t) p * k] * A[j + (size_t) p * k];
    if (!(s > 0)) return 0;
    double ljj = sqrt((double) s);
    A[j + (size_t) p * j] = ljj;
    ld += 2 * logl(ljj);
    for (int i = j + 1; i < p; i++) {
      long double t = A[i + (size_t) p * j];
      for (int k = 0; k < j; k++) t -= (long double) A[i + (size_t) p * k] * A[j + (size_t) p * k];
      A[i + (size_t) p * j] = (double) (t / ljj);
    }
  }
  *logdet = (double) ld;
  return 1;
}


// Solves L L' x = b in place.
static void chol_solve (const double *L, int p, double *x) {
  for (int i = 0; i < p; i++) {
    long double s = x[i];
    for (int k = 0; k < i; k++) s -= (long double) L[i + (size_t) p * k] * x[k];
    x[i] = (double) (s / L[i + (size_t) p * i]);
  }
  for (int i = p - 1; i >= 0; i--) {
    long double s = x[i];
    for (int k = i + 1; k < p; k++) s -= (long double) L[k + (size_t) p * i] * x[k];
    x[i] = (double) (s / L[i + (size_t) p * i]);
  }
}


// M, b, e, r2 and log|M| at t. Returns 0 if M isn't positive definite.
static int lmm_eval (lmm_work_t *w, double t) {

  const lmm_design_t *d = w->d;
  int n = d->n, p = d->p, J = d->J;
  double t2 = t * t;

  for (int a = 0; a < p; a++) {
    w->b[a] = w->wxy[a];
    for (int c = 0; c <= a; c++) w->M[a + (size_t) p * c] = d->wxx[a + (size_t) p * c];
  }
  for (int j = 0; j < J; j++) {
    double cj = d->nj[j] / (1 + d->nj[j] * t2);
    for (int a = 0; a < p; a++) {
      double xa = d->xbar[j + (size_t) J * a];
      w->b[a] += cj * w->ybar[j] * xa;
      for (int c = 0; c <= a; c++) w->M[a + (size_t) p * c] += cj * xa * d->xbar[j + (size_t) J * c];
    }
  }
  for (int a = 0; a < p; a++) for (int c = 0; c < a; c++) w->M[c + (size_t) p * a] = w->M[a + (size_t) p * c];

  if (!chol(w->M, p, &w->logdetM)) return (w->ok = 0);
  chol_solve(w->M, p, w->b);

  long double ss = 0;
  for (int i = 0; i < n; i++) {
    long double r = w->yc[i];
    for (int a = 0; a < p; a++) r -= (long double) d->xc[i + (size_t) n * a] * w->b[a];
    ss += r * r;
  }
  for (int j = 0; j < J; j++) {
    long double ej = w->ybar[j];
    for (int a = 0; a < p; a++) ej -= (long double) d->xbar[j + (size_t) J * a] * w->b[a];
    w->e[j] = (double) ej;
    ss += (d->nj[j] / (1 + d->nj[j] * t2)) * ej * ej;
  }
  w->r2 = (double) ss;
  return (w->ok = w->r2 > 0 && R_FINITE(w->r2));
}


static double lmm_dev (lmm_work_t *w, double t) {
  const lmm_design_t *d = w->d;
  if (!lmm_eval(w, t)) return R_PosInf;
  long double s = 0;
  for (int j = 0; j < d->J; j++) s += log1p(d->nj[j] * t * t);
  int np = d->n - d->p;
  return (double) (s + w->logdetM + np * (1 + log(2 * M_PI * w->r2 / np)));
}


// The inverse of M from its Cholesky factor, into Mi.
static void lmm_inverse (lmm_work_t *w) {
  int p = w->d->p;
  for (int c = 0; c < p; c++) {
    for (int a = 0; a < p; a++) w->v[a] = a == c;
    chol_solve(w->M, p, w->v);
    for (int a = 0; a < p; a++) w->Mi[a + (size_t) p * c] = w->v[a];
  }
}


// h(t) = d'(t) / t, at t (after lmm_eval(w, t)).
static double lmm_h (lmm_work_t *w, double t) {

  const lmm_design_t *d = w->d;
  int p = d->p, J = d->J, np = d->n - d->p;
  double t2 = t * t;

  lmm_inverse(w);

  long double s = 0, tr = 0, rr = 0;
  for (int j = 0; j < J; j++) {
    double nj = d->nj[j], aj = 1 + nj * t2;
    double vj = 2 * nj * nj / (aj * aj);           // -c_j' / t
    s += 2 * nj / aj;
    long double q = 0;                             // xbar_j' Mi xbar_j
    for (int a = 0; a < p; a++) {
      long double r = 0;
      for (int c = 0; c < p; c++) r += (long double) w->Mi[a + (size_t) p * c] * d->xbar[j + (size_t) J * c];
      q += r * d->xbar[j + (size_t) J * a];
    }
    tr -= vj * q;
    rr -= vj * (long double) w->e[j] * w->e[j];
  }
  return (double) (s + tr + np * rr / w->r2);
}


// Brent's minimisation of the deviance on [a, b] (Brent 1973, ch. 5:
// golden sections with parabolic steps).
static double lmm_brent (lmm_work_t *w, double a, double b, double tol) {
  const double cg = 0.3819660112501051;
  double x, v, wv, fx, fv, fw, e = 0, dd = 0;
  x = v = wv = a + cg * (b - a);
  fx = fv = fw = lmm_dev(w, x);
  for (int it = 0; it < 200; it++) {
    double xm = 0.5 * (a + b), tol1 = tol * fabs(x) + 1e-12, tol2 = 2 * tol1;
    if (fabs(x - xm) <= tol2 - 0.5 * (b - a)) break;
    int golden = 1;
    if (fabs(e) > tol1) {
      double r = (x - wv) * (fx - fv), q = (x - v) * (fx - fw), pp = (x - v) * q - (x - wv) * r;
      q = 2 * (q - r);
      if (q > 0) pp = -pp; else q = -q;
      double etemp = e;
      e = dd;
      if (!(fabs(pp) >= fabs(0.5 * q * etemp) || pp <= q * (a - x) || pp >= q * (b - x))) {
        dd = pp / q;
        double u = x + dd;
        if (u - a < tol2 || b - u < tol2) dd = xm >= x ? tol1 : -tol1;
        golden = 0;
      }
    }
    if (golden) { e = (x >= xm ? a : b) - x; dd = cg * e; }
    double u  = fabs(dd) >= tol1 ? x + dd : x + (dd > 0 ? tol1 : -tol1);
    double fu = lmm_dev(w, u);
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = wv; fv = fw; wv = x; fw = fx; x = u; fx = fu;
    }
    else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || wv == x) { v = wv; fv = fw; wv = u; fw = fu; }
      else if (fu <= fv || v == x || v == wv) { v = u; fv = fu; }
    }
  }
  return x;
}


// The REML estimate of t.
static double lmm_theta (lmm_work_t *w) {

  // An upper end at which the deviance increases.
  double hi = 2;
  for (int k = 0; k < 30; k++) {
    if (!lmm_eval(w, hi)) return R_NaN;
    if (lmm_h(w, hi) > 0) break;
    hi *= 4;
  }

  double t = lmm_brent(w, 0, hi, 1e-8);

  // Bracket the root of h at the minimum: h(lo) < 0 < h(up).
  if (!(t > 0) || !lmm_eval(w, t)) return 0;
  double ht = lmm_h(w, t);
  double lo, up, hlo, hup;
  if (ht == 0) return t;
  if (ht > 0) {
    up = t; hup = ht;
    int found = 0;
    double x = t;
    for (int k = 0; k < 60; k++) {
      x *= 0.5;
      if (!lmm_eval(w, x)) return 0;
      double hx = lmm_h(w, x);
      if (hx < 0) { lo = x; hlo = hx; found = 1; break; }
      up = x; hup = hx;
    }
    if (!found) return 0;   // increasing all the way down: the boundary
  }
  else {
    lo = t; hlo = ht;
    int found = 0;
    double x = t;
    for (int k = 0; k < 200; k++) {
      x = x * 1.5 + 1e-8;
      if (!lmm_eval(w, x)) return R_NaN;
      double hx = lmm_h(w, x);
      if (hx > 0) { up = x; hup = hx; found = 1; break; }
      lo = x; hlo = hx;
    }
    if (!found) return lo;
  }

  // Illinois method, then bisection to the last bit.
  int side = 0;
  for (int k = 0; k < 300; k++) {
    if (up - lo <= 4 * DBL_EPSILON * up) break;
    double m = (lo * hup - up * hlo) / (hup - hlo);
    if (!(m > lo && m < up) || k > 100) m = 0.5 * (lo + up);
    if (!lmm_eval(w, m)) return R_NaN;
    double hm = lmm_h(w, m);
    if (hm == 0) return m;
    if (hm < 0) { lo = m; hlo = hm; if (side == -1) hup /= 2; side = -1; }
    else        { up = m; hup = hm; if (side == +1) hlo /= 2; side = +1; }
  }
  return 0.5 * (lo + up);
}


// Fits response y (n), at the REML estimate of t, or at `theta` when it is
// >= 0. Writes coef (p), vcov (p x p), df (p), and returns t; *sigma gets
// the residual standard deviation and *dev the REML criterion (lme4's
// REMLcrit()). NaN where it fails.
double lmm_fit(
    const lmm_design_t *d, const double *y, double theta, double *mem,
    double *coef, double *vcov, double *df, double *sigma, double *dev ) {

  int n = d->n, p = d->p, J = d->J, np = n - p;
  lmm_work_t w;
  w.d    = d;
  w.yc   = mem;            mem += n;
  w.ybar = mem;            mem += J;
  w.e    = mem;            mem += J;
  w.wxy  = mem;            mem += p;
  w.b    = mem;            mem += p;
  w.v    = mem;            mem += p;
  w.M    = mem;            mem += (size_t) p * p;
  w.Mi   = mem;            mem += (size_t) p * p;
  w.t1   = mem;            mem += (size_t) p * p;
  w.t2   = mem;            mem += (size_t) p * p;
  w.t3   = mem;            mem += (size_t) p * p;
  w.t4   = mem;

  for (int k = 0; k < p; k++) { coef[k] = df[k] = R_NaN; for (int c = 0; c < p; c++) vcov[k + (size_t) p * c] = R_NaN; }
  *sigma = R_NaN;
  *dev   = R_NaN;

  // Sums by group, compensated (Kahan; e holds the compensations).
  for (int j = 0; j < J; j++) { w.ybar[j] = 0; w.e[j] = 0; }
  for (int i = 0; i < n; i++) {
    int g = d->grp[i];
    double yv = y[i] - w.e[g];
    double t  = w.ybar[g] + yv;
    w.e[g] = (t - w.ybar[g]) - yv;
    w.ybar[g] = t;
  }
  for (int j = 0; j < J; j++) w.ybar[j] /= d->nj[j];
  for (int i = 0; i < n; i++) w.yc[i] = y[i] - w.ybar[d->grp[i]];
  for (int a = 0; a < p; a++) {
    long double s = 0;
    for (int i = 0; i < n; i++) s += (long double) d->xc[i + (size_t) n * a] * w.yc[i];
    w.wxy[a] = (double) s;
  }

  double t = theta >= 0 ? theta : lmm_theta(&w);
  if (ISNAN(t)) return R_NaN;
  *dev = lmm_dev(&w, t);
  if (!w.ok) return R_NaN;
  lmm_inverse(&w);

  double s2 = w.r2 / np, s = sqrt(s2);
  *sigma = s;
  for (int k = 0; k < p; k++) {
    coef[k] = w.b[k];
    for (int c = 0; c < p; c++) vcov[k + (size_t) p * c] = s2 * w.Mi[k + (size_t) p * c];
  }

  //--------------------------------------------------------
  // Satterthwaite. M' and M'' (t1, t2), Mi M' Mi (t3).
  //--------------------------------------------------------
  double tt = t * t;
  long double dHH = 0, r2p = 0, r2pp = 0;
  for (int a = 0; a < p * p; a++) w.t1[a] = w.t2[a] = 0;
  for (int a = 0; a < p; a++) w.v[a] = 0;   // u = sum c_j' e_j xbar_j
  for (int j = 0; j < J; j++) {
    double nj = d->nj[j], aj = 1 + nj * tt;
    double c1 = -2 * t * nj * nj / (aj * aj);
    double c2 = 2 * nj * nj * (3 * nj * tt - 1) / (aj * aj * aj);
    dHH  += 2 * nj * (1 - nj * tt) / (aj * aj);
    r2p  += c1 * w.e[j] * w.e[j];
    r2pp += c2 * w.e[j] * w.e[j];
    for (int a = 0; a < p; a++) {
      double xa = d->xbar[j + (size_t) J * a];
      w.v[a] += c1 * w.e[j] * xa;
      for (int c = 0; c < p; c++) {
        double xx = xa * d->xbar[j + (size_t) J * c];
        w.t1[a + (size_t) p * c] += c1 * xx;
        w.t2[a + (size_t) p * c] += c2 * xx;
      }
    }
  }
  // u' Mi u
  long double uMu = 0;
  for (int a = 0; a < p; a++) for (int c = 0; c < p; c++)
    uMu += (long double) w.v[a] * w.Mi[a + (size_t) p * c] * w.v[c];
  r2pp -= 2 * uMu;

  // t4 = Mi M'; t3 = Mi M' Mi.
  for (int a = 0; a < p; a++) for (int c = 0; c < p; c++) {
    long double q = 0;
    for (int k = 0; k < p; k++) q += (long double) w.Mi[a + (size_t) p * k] * w.t1[k + (size_t) p * c];
    w.t4[a + (size_t) p * c] = (double) q;
  }
  for (int a = 0; a < p; a++) for (int c = 0; c < p; c++) {
    long double q = 0;
    for (int k = 0; k < p; k++) q += (long double) w.t4[a + (size_t) p * k] * w.Mi[k + (size_t) p * c];
    w.t3[a + (size_t) p * c] = (double) q;
  }
  long double trMM2 = 0, trMiM1sq = 0;
  for (int a = 0; a < p; a++) for (int k = 0; k < p; k++) {
    trMM2    += (long double) w.Mi[a + (size_t) p * k] * w.t2[k + (size_t) p * a];
    trMiM1sq += (long double) w.t4[a + (size_t) p * k] * w.t4[k + (size_t) p * a];
  }

  double Dtt = (double) (dHH + trMM2 - trMiM1sq + r2pp / s2);
  double Dts = (double) (-2 * r2p / (s2 * s));
  double Dss = -2.0 * np / s2 + 6 * w.r2 / (s2 * s2);
  double det = Dtt * Dss - Dts * Dts;
  if (!(det != 0) || !R_FINITE(det)) return t;
  // A = 2 Hess^-1
  double Att =  2 * Dss / det, Ats = -2 * Dts / det, Ass = 2 * Dtt / det;

  for (int k = 0; k < p; k++) {
    double V  = vcov[k + (size_t) p * k];
    double gt = -s2 * w.t3[k + (size_t) p * k];
    double gs = 2 * s * w.Mi[k + (size_t) p * k];
    double q  = gt * gt * Att + 2 * gt * gs * Ats + gs * gs * Ass;
    df[k] = q > 0 ? 2 * V * V / q : R_NaN;
  }

  return t;
}



//======================================================
// R interface, for the tests and linda()'s adaptive
// choice: each column of Y (n x m) fitted to X (n x p)
// with random intercepts by `groups` (1-based codes),
// at the REML estimates, or at `theta` (one per column)
// if not NULL. Returns list(coef, se, df (m x p),
// theta, sigma, vcov (m x p x p), deviance).
//======================================================
SEXP C_feature_lmm(SEXP sexp_X, SEXP sexp_Y, SEXP sexp_groups, SEXP sexp_theta) {

  int n = nrows(sexp_X), p = ncols(sexp_X), m = ncols(sexp_Y);
  int *grp = (int *) R_alloc((size_t) n, sizeof(int));
  int J = 0;
  for (int i = 0; i < n; i++) { grp[i] = INTEGER(sexp_groups)[i] - 1; if (grp[i] + 1 > J) J = grp[i] + 1; }

  lmm_design_t d;
  lmm_design_init(&d, REAL(sexp_X), grp, n, p, J);

  const char *nms[] = { "coef", "se", "df", "theta", "sigma", "vcov", "deviance", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_coef = allocMatrix(REALSXP, m, p); SET_VECTOR_ELT(result, 0, s_coef);
  SEXP s_se   = allocMatrix(REALSXP, m, p); SET_VECTOR_ELT(result, 1, s_se);
  SEXP s_df   = allocMatrix(REALSXP, m, p); SET_VECTOR_ELT(result, 2, s_df);
  SEXP s_th   = allocVector(REALSXP, m);    SET_VECTOR_ELT(result, 3, s_th);
  SEXP s_sg   = allocVector(REALSXP, m);    SET_VECTOR_ELT(result, 4, s_sg);
  SEXP s_vc   = PROTECT(allocVector(REALSXP, (R_xlen_t) m * p * p));
  SEXP dim    = PROTECT(allocVector(INTSXP, 3));
  INTEGER(dim)[0] = m; INTEGER(dim)[1] = p; INTEGER(dim)[2] = p;
  setAttrib(s_vc, R_DimSymbol, dim);
  SET_VECTOR_ELT(result, 5, s_vc);
  SEXP s_dev  = allocVector(REALSXP, m);    SET_VECTOR_ELT(result, 6, s_dev);

  double *mem  = (double *) R_alloc(lmm_doubles(n, p, J), sizeof(double));
  double *coef = (double *) R_alloc((size_t) p, sizeof(double));
  double *df   = (double *) R_alloc((size_t) p, sizeof(double));
  double *vc   = (double *) R_alloc((size_t) p * p, sizeof(double));

  for (int c = 0; c < m; c++) {
    double sigma, dev;
    double th = isNull(sexp_theta) ? -1 : REAL(sexp_theta)[c];
    double t  = lmm_fit(&d, REAL(sexp_Y) + (size_t) n * c, th, mem, coef, vc, df, &sigma, &dev);
    REAL(s_dev)[c] = dev;
    REAL(s_th)[c] = t;
    REAL(s_sg)[c] = sigma;
    for (int k = 0; k < p; k++) {
      REAL(s_coef)[c + (size_t) m * k] = coef[k];
      REAL(s_se)[c + (size_t) m * k]   = sqrt(vc[k + (size_t) p * k]);
      REAL(s_df)[c + (size_t) m * k]   = df[k];
      for (int l = 0; l < p; l++)
        REAL(s_vc)[c + (size_t) m * (k + (size_t) p * l)] = vc[k + (size_t) p * l];
    }
  }

  UNPROTECT(3);
  return result;
}
