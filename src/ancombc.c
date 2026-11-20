// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// C port of the core of ANCOM-BC and ANCOM-BC2 from scikit-bio 0.7.4
// (skbio/stats/composition/_ancombc.py: _ancombc_core(), _transform_data(),
// _estimate_params_dense(), _estimate_params_sparse(),
// _lstsq_sparse_batch(), _apply_pinv(), _solve_sparse_iter(),
// _r_fit_info(), _r_rebase_categorical(), _invert_singular(),
// _calc_grouped_var_cov(), _estimate_bias_em(), _init_bias_params(),
// _estimate_bias_var() and _sample_fractions()), which is itself derived
// from the R package ANCOMBC. Modified for ecodive: restructured as
// threaded C over blocks of features, with ecodive's own singular value
// decomposition (one-sided Jacobi) in place of LAPACK's, and with the
// design matrix, its pseudoinverse and the per-covariate statistics
// computed in R (R/ancombc.r). See plans/ancombc.md for the differences
// from scikit-bio and ANCOMBC. scikit-bio's license:
//
//   Copyright (c) 2013--, scikit-bio development team.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   * Redistributions of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//
//   * Neither the names scikit-bio, skbio, or biocore nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
//   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
//   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// scikit-bio's implementation "is based on an analysis of the source code
// of the R package ANCOMBC", whose ANCOM-BC and ANCOM-BC2 code it follows,
// and ANCOMBC's source (version 2.8.1 and 2.14.0) was read to check this
// port. ANCOMBC's notice:
//
//   ANCOMBC: Copyright (c) Huang Lin.
//   Licensed under the Artistic License 2.0. The full text of the license
//   is in inst/COPYRIGHTS, and at
//   https://opensource.org/license/artistic-2-0.
//
//   This file is a modified version of that code (Artistic License 2.0,
//   section 4): a C port, restructured for ecodive, of scikit-bio 0.7.4's
//   Python port of it, distributed as part of the package ecodive, which
//   is not ANCOMBC and doesn't prevent ANCOMBC from being installed.
//
// The bounded one-dimensional Nelder-Mead minimiser, nm_minimize(), follows
// _minimize_neldermead() from SciPy 1.18.1 (scipy/optimize/_optimize.py),
// which scikit-bio calls with method="Nelder-Mead", bounds [(0, None)] and
// xatol = fatol = 1e-8. SciPy's license:
//
//   Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//
//   3. Neither the name of the copyright holder nor the names of its
//      contributors may be used to endorse or promote products derived
//      from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
// Outline (n samples, D features, P design columns; X the n x P design):
//
//   1. L = log(x + pseudocount), an n x D table, column per feature. With
//      ANCOM-BC2 and no pseudocount, zeros are missing (NaN), and each
//      feature is centred on the mean of its observed logs.
//   2. The initial fit (biased): per feature, beta_f = X+ L_f and the
//      residuals; each sample's mean residual over the features (theta) is
//      the starting estimate of its sampling fraction; the coefficients'
//      variances (and covariances) are the sandwich sums of the squared,
//      theta-centred residuals with H = (X+)'. With missing values, each
//      feature has its own masked design (rows of missing values zeroed)
//      and pseudoinverse, beta and theta are refined alternately (up to 20
//      rounds, to 1e-2), R's whole-fit failures are found (a factor with
//      fewer than two observed levels), and a factor whose reference level
//      is unobserved is re-expressed against its first observed level.
//   3. The bias EM, per covariate: a three-component Gaussian mixture of the
//      coefficients over the features, whose null component's mean is the
//      bias (delta). Covariates are independent tasks.
//   4. ANCOM-BC2 only: each sample's sampling fraction from the initial fit
//      and the biases, subtracted from L, and a final, unbiased fit.
//
// Threads: features are split into blocks of AB_BLK; a block's per-sample
// sums (for theta and the sampling fractions) are kept separately and added
// in block order, and every per-feature quantity depends on that feature
// alone, so results are identical for any number of threads.

#include "ecodive.h"
#include <float.h> // DBL_EPSILON, DBL_MAX

#define AB_BLK      64      // Features per block.
#define AB_MAX_COND 1e4     // _lstsq_sparse_batch(max_cond).
#define AB_SP_TOL   1e-2    // _estimate_params_sparse(tol).
#define AB_SP_ITERS 20      // _estimate_params_sparse(max_iter).
#define NM_XATOL    1e-8    // OPTARG
#define NM_FATOL    1e-8
#define NM_MAXFUN   200     // SciPy's default for one variable: 200 N.


//======================================================
// Shared state for the workers.
//======================================================
typedef struct {
  int           n, D, P, nb;
  const double *X;        // n x P design.
  const double *H;        // n x P: (X+)', the global influence matrix.
  const double *Xinv;     // P x n: X+, the global pseudoinverse.
  const double *gsum;     // P: (X'X)+ 1, row sums of the inverse Gram.
  double       *data;     // n x D logs; NaN = missing.
  int           sparse;   // Any missing values.

  // Sparse: each feature's masked-design SVD.
  double       *V;        // D blocks of P x P: right singular vectors (columns).
  double       *Sinv;     // D x P: inverted singular values, 0 when dropped.
  int          *rank;     // D
  double      **pinv;     // D: P x n full pseudoinverse of an ill-conditioned
                          //    masked design, else NULL.
  const int    *valid;    // D (1/0), or NULL for all valid.
  const int    *nmiss;    // D: missing values per feature.

  // The current fit.
  double       *beta;     // D blocks of P (feature f's at beta + f * P).
  double       *beta_new;
  const double *theta;    // n, or NULL for 0.
  double       *part;     // nb x n: per-block sums over features.
  double       *part2;    // nb: per-block sums.
  double       *scratch;  // Per thread: AB_SCR doubles.
  int           scr_len;
  int           failed;   // Set by a worker whose malloc() failed.

  // Variances and covariances.
  double       *var;      // D x P
  double       *cov;      // D x P x P, or NULL.
  double       *rss;      // D: the residual sums of squares (NaN if failed).

  // Sampling fractions.
  const double *delta;    // P: the EM biases.
} ab_state_t;

static ab_state_t ab;



//======================================================
// One-sided (Hestenes) Jacobi SVD of the n x P column-
// major matrix A (n >= P), in place: on return A's
// columns are U S (unnormalised left singular vectors),
// V (P x P) holds the right singular vectors, and S the
// singular values, sorted in decreasing order with A and
// V's columns.
//======================================================
static void jacobi_svd (double *A, int n, int P, double *V, double *S) {

  for (int i = 0; i < P * P; i++) V[i] = 0;
  for (int i = 0; i < P; i++)     V[i + i * P] = 1;

  for (int sweep = 0; sweep < 100; sweep++) {
    int rotated = 0;

    for (int i = 0; i < P - 1; i++) {
      double *ai = A + (size_t) i * n;
      for (int j = i + 1; j < P; j++) {
        double *aj = A + (size_t) j * n;

        double alpha = 0, beta = 0, gamma = 0;
        for (int s = 0; s < n; s++) {
          alpha += ai[s] * ai[s];
          beta  += aj[s] * aj[s];
          gamma += ai[s] * aj[s];
        }
        if (gamma == 0 || fabs(gamma) <= DBL_EPSILON * sqrt(alpha * beta))
          continue;

        double zeta = (beta - alpha) / (2 * gamma);
        double t    = (zeta >= 0 ? 1.0 : -1.0) / (fabs(zeta) + sqrt(1 + zeta * zeta));
        double c    = 1 / sqrt(1 + t * t);
        double sn   = c * t;

        for (int s = 0; s < n; s++) {
          double x = ai[s], y = aj[s];
          ai[s] = c * x - sn * y;
          aj[s] = sn * x + c * y;
        }
        double *vi = V + (size_t) i * P, *vj = V + (size_t) j * P;
        for (int p = 0; p < P; p++) {
          double x = vi[p], y = vj[p];
          vi[p] = c * x - sn * y;
          vj[p] = sn * x + c * y;
        }
        rotated = 1;
      }
    }
    if (!rotated) break;
  }

  for (int j = 0; j < P; j++) {
    double ss = 0, *aj = A + (size_t) j * n;
    for (int s = 0; s < n; s++) ss += aj[s] * aj[s];
    S[j] = sqrt(ss);
  }

  // Sort decreasing (selection sort; P is small).
  for (int i = 0; i < P - 1; i++) {
    int m = i;
    for (int j = i + 1; j < P; j++) if (S[j] > S[m]) m = j;
    if (m == i) continue;
    double t = S[i]; S[i] = S[m]; S[m] = t;
    double *ai = A + (size_t) i * n, *am = A + (size_t) m * n;
    for (int s = 0; s < n; s++) { t = ai[s]; ai[s] = am[s]; am[s] = t; }
    double *vi = V + (size_t) i * P, *vm = V + (size_t) m * P;
    for (int p = 0; p < P; p++) { t = vi[p]; vi[p] = vm[p]; vm[p] = t; }
  }
}



//======================================================
// numpy.quantile(x, q) (method 'linear') of m values,
// sorted in place.
//======================================================
static int cmp_dbl (const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

static double np_quantile_sorted (const double *x, int m, double q) {
  double vi = (m - 1) * q;
  double lo = floor(vi);
  int    i0, i1;
  if (vi >= m - 1) { i0 = i1 = m - 1; lo = -1; }
  else if (vi < 0) { i0 = i1 = 0; }
  else             { i0 = (int) lo; i1 = i0 + 1; }
  double t = vi - lo;
  double a = x[i0], b = x[i1], d = b - a;
  return (t >= 0.5) ? b - d * (1 - t) : a + d * t;
}



//======================================================
// SciPy's Nelder-Mead for one variable bounded below by
// 0, minimising
//   f(x) = 0.5 * sum_i resp_i (sq_i / (var_i + x) + log(var_i + x)).
//======================================================
typedef struct {
  const double *var, *sq, *resp;
  int m;
} nm_obj_t;

static double nm_f (const nm_obj_t *o, double x) {
  double s = 0;
  for (int i = 0; i < o->m; i++) {
    double v = o->var[i] + x;
    s += o->resp[i] * (o->sq[i] / v + log(v));
  }
  return 0.5 * s;
}

static inline double nm_clip (double x) { return (x < 0) ? 0 : x; }

// numpy.argsort of two values (NaN last; ties keep their order).
static inline void nm_sort (double *sim, double *fsim) {
  if (fsim[1] < fsim[0] || (ISNAN(fsim[0]) && !ISNAN(fsim[1]))) {
    double t = sim[0];  sim[0]  = sim[1];  sim[1]  = t;
    t = fsim[0];        fsim[0] = fsim[1]; fsim[1] = t;
  }
}

static double nm_minimize (const nm_obj_t *o, double x0) {

  double sim[2], fsim[2];
  int    fcalls = 0, iterations;

  x0     = nm_clip(x0);
  sim[0] = x0;
  sim[1] = (x0 != 0) ? (1 + 0.05) * x0 : 0.00025;
  sim[1] = nm_clip(sim[1]);

  fsim[0] = nm_f(o, sim[0]); fcalls++;
  fsim[1] = nm_f(o, sim[1]); fcalls++;
  nm_sort(sim, fsim);

  #define NM_EVAL(out, x) { if (fcalls >= NM_MAXFUN) goto step_done; fcalls++; out = nm_f(o, x); }

  iterations = 1;
  while (fcalls < NM_MAXFUN && iterations < NM_MAXFUN) {

    if (fabs(sim[1] - sim[0]) <= NM_XATOL && fabs(fsim[0] - fsim[1]) <= NM_FATOL)
      break;

    {
      double xbar = sim[0];               // (sum of all but the worst) / N
      double xr   = nm_clip(2 * xbar - 1 * sim[1]);
      double fxr, fxe, fxc, fxcc;
      int    doshrink = 0;
      NM_EVAL(fxr, xr);

      if (fxr < fsim[0]) {
        double xe = nm_clip(3 * xbar - 2 * sim[1]);
        NM_EVAL(fxe, xe);
        if (fxe < fxr) { sim[1] = xe; fsim[1] = fxe; }
        else           { sim[1] = xr; fsim[1] = fxr; }
      }
      else if (fxr < fsim[0]) {           // fsim[-2] is fsim[0] for N = 1.
        sim[1] = xr; fsim[1] = fxr;
      }
      else {
        if (fxr < fsim[1]) {
          double xc = nm_clip(1.5 * xbar - 0.5 * sim[1]);
          NM_EVAL(fxc, xc);
          if (fxc <= fxr) { sim[1] = xc; fsim[1] = fxc; }
          else            doshrink = 1;
        }
        else {
          double xcc = nm_clip(0.5 * xbar + 0.5 * sim[1]);
          NM_EVAL(fxcc, xcc);
          if (fxcc < fsim[1]) { sim[1] = xcc; fsim[1] = fxcc; }
          else                doshrink = 1;
        }
        if (doshrink) {
          sim[1] = nm_clip(sim[0] + 0.5 * (sim[1] - sim[0]));
          NM_EVAL(fsim[1], sim[1]);
        }
      }
      iterations++;
    }

    step_done:
    nm_sort(sim, fsim);
  }
  #undef NM_EVAL

  return sim[0];
}



//======================================================
// The bias EM for one covariate's m coefficients (beta)
// and variances (var). w: 14 m doubles of scratch.
// out: delta_em, delta_wls, var_delta; returns the
// number of EM iterations run.
//======================================================
static int bias_em (
    const double *beta, const double *var, int m, double tol, int max_iter,
    double *w, double *out ) {

  double *nu_inv = w;           // 3 x m
  double *sd_inv = w + 3 * m;   // 3 x m
  double *ratios = w + 6 * m;   // 3 x m
  double *resp   = w + 9 * m;   // 3 x m
  double *srt    = w + 12 * m;  // m: sorted copy
  double *tmp    = w + 13 * m;  // m

  for (int i = 0; i < m; i++) srt[i] = beta[i];
  qsort(srt, (size_t) m, sizeof(double), cmp_dbl);


  //--------------------------------------------------------
  // _init_bias_params()
  //--------------------------------------------------------
  double e0 = np_quantile_sorted(srt, m, 0.125);
  double e1 = np_quantile_sorted(srt, m, 0.25);
  double e2 = np_quantile_sorted(srt, m, 0.75);
  double e3 = np_quantile_sorted(srt, m, 0.875);

  double delta, l1, l2, kappa1, kappa2, s;
  int    c;

  s = 0; c = 0;
  for (int i = 0; i < m; i++) if (beta[i] >= e1 && beta[i] <= e2) { s += beta[i]; c++; }
  if (c) delta = s / c;
  else {
    s = 0; for (int i = 0; i < m; i++) s += beta[i];
    delta = s / m;
  }

  s = 0; c = 0;
  for (int i = 0; i < m; i++) if (beta[i] < e0) { s += beta[i]; c++; }
  if (c) {
    l1 = s / c;
    if (c > 1) {
      double v = 0;
      for (int i = 0; i < m; i++) if (beta[i] < e0) v += (beta[i] - l1) * (beta[i] - l1);
      v /= (c - 1);
      kappa1 = (v == 0) ? 1.0 : v;
    }
    else kappa1 = 1.0;
  }
  else { l1 = srt[0]; kappa1 = 1.0; }

  s = 0; c = 0;
  for (int i = 0; i < m; i++) if (beta[i] > e3) { s += beta[i]; c++; }
  if (c) {
    l2 = s / c;
    if (c > 1) {
      double v = 0;
      for (int i = 0; i < m; i++) if (beta[i] > e3) v += (beta[i] - l2) * (beta[i] - l2);
      v /= (c - 1);
      kappa2 = (v == 0) ? 1.0 : v;
    }
    else kappa2 = 1.0;
  }
  else { l2 = srt[m - 1]; kappa2 = 1.0; }

  double params[8]  = { 0.75, 0.125, 0.125, delta, l1, l2, kappa1, kappa2 };
  double updated[8];


  //--------------------------------------------------------
  // _estimate_bias_em()
  //--------------------------------------------------------
  for (int i = 0; i < m; i++) {
    nu_inv[i] = 1.0 / var[i];
    sd_inv[i] = sqrt(nu_inv[i]);
    ratios[i] = beta[i] / var[i];
  }

  double loss = R_PosInf;
  int    it   = 0;

  while (loss > tol && it < max_iter) {

    for (int r = 1; r <= 2; r++) {
      double kap = params[5 + r], l = params[3 + r];
      double *nu = nu_inv + r * m, *sd = sd_inv + r * m, *ra = ratios + r * m;
      for (int i = 0; i < m; i++) {
        nu[i] = 1.0 / (var[i] + kap);
        sd[i] = sqrt(nu[i]);
        ra[i] = (beta[i] - l) * nu[i];
      }
    }

    // E-step
    double means[3];
    delta    = params[3];
    means[0] = delta;
    means[1] = params[4] + delta;
    means[2] = params[5] + delta;

    for (int r = 0; r < 3; r++) {
      double *re = resp + r * m, *nu = nu_inv + r * m, *sd = sd_inv + r * m;
      for (int i = 0; i < m; i++) {
        double d = beta[i] - means[r];
        d  = d * d;
        d *= nu[i];
        d *= -0.5;
        re[i] = exp(d) * sd[i] * params[r];
      }
    }
    for (int i = 0; i < m; i++) {
      double tot = resp[i] + resp[m + i] + resp[2 * m + i];
      resp[i] /= tot; resp[m + i] /= tot; resp[2 * m + i] /= tot;
    }

    // M-step
    for (int r = 0; r < 3; r++) {
      double t = 0, *re = resp + r * m;
      for (int i = 0; i < m; i++) t += re[i];
      updated[r] = t / m;
    }

    double num = 0, den = 0;
    for (int i = 0; i < 3 * m; i++) { num += resp[i] * ratios[i]; den += resp[i] * nu_inv[i]; }
    updated[3] = num / den;

    double ln[2];
    for (int r = 1; r <= 2; r++) {
      double dn = 0, nm = 0, *re = resp + r * m, *nu = nu_inv + r * m;
      for (int i = 0; i < m; i++) {
        double x = re[i] * nu[i];
        dn += x;
        nm += x * (beta[i] - delta);
      }
      ln[r - 1] = nm / dn;
    }
    updated[4] = ISNAN(ln[0]) ? ln[0] : (ln[0] < 0 ? ln[0] : 0);
    updated[5] = ISNAN(ln[1]) ? ln[1] : (ln[1] > 0 ? ln[1] : 0);

    // kappa1, kappa2: squared residuals about the iteration's means.
    for (int r = 1; r <= 2; r++) {
      for (int i = 0; i < m; i++) {
        double d = beta[i] - means[r];
        tmp[i] = d * d;
      }
      nm_obj_t o = { .var = var, .sq = tmp, .resp = resp + r * m, .m = m };
      updated[5 + r] = nm_minimize(&o, params[5 + r]);
    }

    double ss = 0;
    for (int k = 0; k < 8; k++) ss += (updated[k] - params[k]) * (updated[k] - params[k]);
    loss = sqrt(ss);

    for (int k = 0; k < 8; k++) params[k] = updated[k];
    it++;
  }


  //--------------------------------------------------------
  // _estimate_bias_var()
  //--------------------------------------------------------
  double q1 = np_quantile_sorted(srt, m, params[1]);
  double q2 = np_quantile_sorted(srt, m, 1.0 - params[2]);

  double wden = 0, wnum = 0;
  for (int i = 0; i < m; i++) {
    double nu = var[i], b = beta[i];
    if (beta[i] <  q1) { nu += params[6]; b -= params[4]; }
    if (beta[i] >= q2) { nu += params[7]; b -= params[5]; }
    double ni = 1.0 / nu;
    wden += ni;
    wnum += ni * b;
  }

  double dinv = 1.0 / wden;
  out[0] = params[3];
  out[1] = wnum * dinv;
  out[2] = ISNAN(dinv) ? 0 : (dinv == R_PosInf ? DBL_MAX : (dinv == R_NegInf ? -DBL_MAX : dinv));

  return it;
}



//======================================================
// Applies feature f's masked-design pseudoinverse to
// resp (n; 0 where missing): out = X_f+ resp.
// w: 2 P doubles.
//======================================================
static void apply_pinv (int f, const double *resp, double *out, double *w) {

  int n = ab.n, P = ab.P;

  if (ab.pinv[f]) {
    const double *pi = ab.pinv[f];
    for (int p = 0; p < P; p++) out[p] = 0;
    for (int s = 0; s < n; s++) {
      double r = resp[s];
      for (int p = 0; p < P; p++) out[p] += pi[p + (size_t) s * P] * r;
    }
    return;
  }

  double *rhs = w, *tmp = w + P;
  const double *V  = ab.V + (size_t) f * P * P;
  const double *Si = ab.Sinv + (size_t) f * P;

  for (int p = 0; p < P; p++) {
    const double *x = ab.X + (size_t) p * n;
    double t = 0;
    for (int s = 0; s < n; s++) t += x[s] * resp[s];
    rhs[p] = t;
  }
  for (int k = 0; k < P; k++) {
    const double *v = V + (size_t) k * P;
    double t = 0;
    for (int p = 0; p < P; p++) t += v[p] * rhs[p];
    t *= Si[k];
    t *= Si[k];
    tmp[k] = t;
  }
  for (int p = 0; p < P; p++) {
    double t = 0;
    for (int k = 0; k < P; k++) t += V[p + (size_t) k * P] * tmp[k];
    out[p] = t;
  }
}



//======================================================
// Worker: each feature's masked-design SVD (sparse).
// Scratch: n P + P P + P doubles, plus n P for U.
//======================================================
static void *svd_worker (void *arg) {

  worker_t *wk = (worker_t *) arg;
  int n = ab.n, P = ab.P;
  double *A  = ab.scratch + (size_t) wk->i * ab.scr_len;
  double *S  = A + (size_t) n * P;

  for (int b = wk->i; b < ab.nb; b += wk->n) {
    int f0 = b * AB_BLK, f1 = f0 + AB_BLK < ab.D ? f0 + AB_BLK : ab.D;

    for (int f = f0; f < f1; f++) {
      const double *y = ab.data + (size_t) f * n;
      double *V = ab.V + (size_t) f * P * P;

      for (int p = 0; p < P; p++)
        for (int s = 0; s < n; s++)
          A[s + (size_t) p * n] = ISNAN(y[s]) ? 0 : ab.X[s + (size_t) p * n];

      jacobi_svd(A, n, P, V, S);

      // _invert_singular()
      double *Si = ab.Sinv + (size_t) f * P;
      double  cut = 1e-15 * S[0];
      int     rk = 0;
      double  min_ret = R_PosInf;
      for (int k = 0; k < P; k++) {
        if (S[k] > cut) { Si[k] = 1.0 / S[k]; rk++; if (S[k] < min_ret) min_ret = S[k]; }
        else            Si[k] = 0;
      }
      ab.rank[f] = rk;

      double cond = (R_FINITE(min_ret) && min_ret > 0) ? S[0] / min_ret : 1;
      ab.pinv[f]  = NULL;

      if (cond > AB_MAX_COND) {
        double *pi = (double *) malloc((size_t) P * n * sizeof(double));
        if (!pi) { ab.failed = 1; continue; }
        // pinv = V diag(Sinv) U', with U = A / S.
        for (int s = 0; s < n; s++)
          for (int p = 0; p < P; p++) {
            double t = 0;
            for (int k = 0; k < P; k++)
              if (Si[k] != 0) t += V[p + (size_t) k * P] * Si[k] * (A[s + (size_t) k * n] / S[k]);
            pi[p + (size_t) s * P] = t;
          }
        ab.pinv[f] = pi;
      }
    }
  }

  return NULL;
}



//======================================================
// Worker: coefficients and per-block residual sums.
//
// Dense:  beta_new_f = X+ y_f; part[b, s] = sum of the
//         block's residuals at s (theta == NULL here).
// Sparse: beta_new_f = X_f+ (y_f - theta), 0 where
//         missing; part[b, s] = sum over the block's valid
//         features observed at s of y - X beta_new;
//         part2[b] = sum over its (valid) features of
//         |beta_new - beta|^2.
//======================================================
static void *fit_worker (void *arg) {

  worker_t *wk = (worker_t *) arg;
  int n = ab.n, P = ab.P;
  double *resp = ab.scratch + (size_t) wk->i * ab.scr_len;
  double *w    = resp + n;

  for (int b = wk->i; b < ab.nb; b += wk->n) {
    int f0 = b * AB_BLK, f1 = f0 + AB_BLK < ab.D ? f0 + AB_BLK : ab.D;
    double *part = ab.part + (size_t) b * n;
    double  dsum = 0;
    for (int s = 0; s < n; s++) part[s] = 0;

    for (int f = f0; f < f1; f++) {
      const double *y  = ab.data + (size_t) f * n;
      double       *bn = ab.beta_new + (size_t) f * P;
      int           ok = ab.valid ? ab.valid[f] : 1;

      if (!ab.sparse) {
        for (int p = 0; p < P; p++) {
          double t = 0;
          for (int s = 0; s < n; s++) t += ab.Xinv[p + (size_t) s * P] * y[s];
          bn[p] = t;
        }
      }
      else {
        for (int s = 0; s < n; s++)
          resp[s] = ISNAN(y[s]) ? 0 : (ab.theta ? y[s] - ab.theta[s] : y[s]);
        apply_pinv(f, resp, bn, w);
      }

      if (!ok) continue;

      for (int s = 0; s < n; s++) {
        if (ISNAN(y[s])) continue;
        double fit = 0;
        for (int p = 0; p < P; p++) fit += ab.X[s + (size_t) p * n] * bn[p];
        part[s] += y[s] - fit;
      }

      if (ab.sparse) {
        const double *bo = ab.beta + (size_t) f * P;
        double t = 0;
        for (int p = 0; p < P; p++) {
          double d = bn[p] - bo[p];
          d = d * d;
          if (ab.valid || !ISNAN(d)) t += d;
        }
        dsum += t;
      }
    }
    ab.part2[b] = dsum;
  }

  return NULL;
}



//======================================================
// Worker: variances and covariances of the coefficients
// from the squared residuals r2 = (y - X beta - theta)^2
// (0 where missing): var_fp = sum_s r2_s H_sp^2, cov_fpq
// = sum_s r2_s H_sp H_sq, plus R's corrections for
// missing values and failed fits (sparse); and the
// residual sum of squares, sum_s r2_s.
//======================================================
static void *var_worker (void *arg) {

  worker_t *wk = (worker_t *) arg;
  int n = ab.n, P = ab.P, D = ab.D;
  double *r2  = ab.scratch + (size_t) wk->i * ab.scr_len;
  double *acc = r2 + n;            // P x P

  for (int b = wk->i; b < ab.nb; b += wk->n) {
    int f0 = b * AB_BLK, f1 = f0 + AB_BLK < D ? f0 + AB_BLK : D;

    for (int f = f0; f < f1; f++) {
      const double *y  = ab.data + (size_t) f * n;
      const double *bt = ab.beta + (size_t) f * P;

      for (int s = 0; s < n; s++) {
        if (ISNAN(y[s])) { r2[s] = 0; continue; }
        double fit = 0;
        for (int p = 0; p < P; p++) fit += ab.X[s + (size_t) p * n] * bt[p];
        double r = y[s] - fit;
        if (ab.theta) r -= ab.theta[s];
        r2[s] = r * r;
      }

      double rss = 0;
      for (int s = 0; s < n; s++) rss += r2[s];
      ab.rss[f] = (ab.valid && !ab.valid[f]) ? R_NaN : rss;

      for (int i = 0; i < P * P; i++) acc[i] = 0;
      for (int s = 0; s < n; s++) {
        double rs = r2[s];
        if (rs == 0) continue;
        for (int p = 0; p < P; p++) {
          double hp = ab.H[s + (size_t) p * n] * rs;
          for (int q = p; q < P; q++) acc[p + q * P] += hp * ab.H[s + (size_t) q * n];
        }
      }

      if (ab.sparse) {
        int failed = ab.valid && !ab.valid[f];
        double c = failed ? 0.1 * n : 0.1 * ab.nmiss[f];
        if (failed) for (int i = 0; i < P * P; i++) acc[i] = 0;
        if (c != 0)
          for (int p = 0; p < P; p++)
            for (int q = p; q < P; q++)
              acc[p + q * P] += c * (ab.gsum[p] * ab.gsum[q]);
      }

      for (int p = 0; p < P; p++) ab.var[f + (size_t) p * D] = acc[p + p * P];

      if (ab.cov)
        for (int p = 0; p < P; p++)
          for (int q = 0; q < P; q++) {
            double v = (p <= q) ? acc[p + q * P] : acc[q + p * P];
            ab.cov[f + (size_t) p * D + (size_t) q * D * P] = v;
          }
    }
  }

  return NULL;
}



//======================================================
// Worker: per-block sums for the sampling fractions
// (sparse): part[b, s] = sum over the block's features
// observed at s of y - X beta + (X delta) valid, and
// part[nb + b, s] the count.
//======================================================
static void *sf_worker (void *arg) {

  worker_t *wk = (worker_t *) arg;
  int n = ab.n, P = ab.P;
  double *xd = ab.scratch + (size_t) wk->i * ab.scr_len;

  for (int s = 0; s < n; s++) {
    double t = 0;
    for (int p = 0; p < P; p++) t += ab.X[s + (size_t) p * n] * ab.delta[p];
    xd[s] = t;
  }

  for (int b = wk->i; b < ab.nb; b += wk->n) {
    int f0 = b * AB_BLK, f1 = f0 + AB_BLK < ab.D ? f0 + AB_BLK : ab.D;
    double *sum = ab.part + (size_t) b * n;
    double *cnt = ab.part + (size_t) (ab.nb + b) * n;
    for (int s = 0; s < n; s++) { sum[s] = 0; cnt[s] = 0; }

    for (int f = f0; f < f1; f++) {
      const double *y  = ab.data + (size_t) f * n;
      const double *bt = ab.beta + (size_t) f * P;
      int ok = ab.valid ? ab.valid[f] : 1;

      for (int s = 0; s < n; s++) {
        if (ISNAN(y[s])) continue;
        double fit = 0;
        if (ok) for (int p = 0; p < P; p++) fit += ab.X[s + (size_t) p * n] * bt[p];
        double v = fit;
        if (ok) v -= xd[s];
        v -= y[s];
        sum[s] += -v;
        cnt[s] += 1;
      }
    }
  }

  return NULL;
}



//======================================================
// Worker: the bias EM, one covariate per task.
//======================================================
static const double *em_beta_all;  // D blocks of P (ab.beta)
static const double *em_var_all;   // D x P
static double       *em_out;       // P x 3
static int          *em_iters;     // P
static int          *em_m;         // P: features used
static int           em_max_iter;
static double        em_tol;

static void *em_worker (void *arg) {

  worker_t *wk = (worker_t *) arg;
  int D = ab.D, P = ab.P;

  double *buf = (double *) malloc((size_t) 16 * (D > 0 ? D : 1) * sizeof(double));
  if (!buf) { ab.failed = 1; return NULL; }
  double *b = buf, *v = buf + D, *w = buf + 2 * (size_t) D;

  for (int i = wk->i; i < P; i += wk->n) {
    int m = 0;
    for (int f = 0; f < D; f++) {
      double vv = em_var_all[f + (size_t) i * D];
      if (ab.valid && !ab.valid[f]) continue;
      if (!R_FINITE(vv)) continue;
      b[m] = em_beta_all[i + (size_t) f * P];
      v[m] = vv;
      m++;
    }
    em_m[i] = m;
    if (m == 0) { em_iters[i] = 0; continue; }
    em_iters[i] = bias_em(b, v, m, em_tol, em_max_iter, w, em_out + 3 * i);
  }

  free(buf);
  return NULL;
}



//======================================================
// Helpers run on the main thread.
//======================================================

// One fit pass: beta_new from (data - theta), and the
// per-sample residual sums over blocks, added in block
// order into sums (n).
static void fit_pass (int nt, double *sums) {
  run_parallel(fit_worker, nt, ab.D);
  for (int s = 0; s < ab.n; s++) sums[s] = 0;
  for (int b = 0; b < ab.nb; b++) {
    const double *part = ab.part + (size_t) b * ab.n;
    for (int s = 0; s < ab.n; s++) sums[s] += part[s];
  }
}

static void swap_beta (void) {
  double *t = ab.beta; ab.beta = ab.beta_new; ab.beta_new = t;
}



//======================================================
// R interface.
//
//   counts, margin  The table.
//   pseudocount     Added to every value; 0 = none (zeros are
//                   then missing, which R allows only for
//                   ANCOM-BC2).
//   v2              ANCOM-BC2: centre the logs, correct for the
//                   sampling fractions, refit.
//   X               n x P design (with the intercept, if any).
//   Xinv            P x n: its pseudoinverse.
//   gsum            P: (X'X)+ 1.
//   factors         A list, one element per factor with a main-
//                   effect term: integer c(standard, rebase, L,
//                   level_col[1..L], codes[1..n]), with standard
//                   = 1 for treatment coding (level 0 the
//                   reference), rebase = 1 when it may be re-
//                   expressed against a feature-local reference,
//                   level_col the 0-based design column of each
//                   level (-1 for the reference), codes each
//                   sample's 0-based level.
//   intercept       0-based column of the intercept, or -1.
//   max_iters, tol  The bias EM.
//   n_threads
//
// Returns a list:
//   beta       D x P coefficients before bias correction (initial fit
//              for ANCOM-BC; final fit for ANCOM-BC2).
//   var        D x P variances (unadjusted).
//   vcov       D x P x P covariances.
//   valid      D logical, or NULL when every fit succeeded.
//   rank       D integer masked-design ranks, or NULL (no missing values).
//   n_missing  D integer, or NULL.
//   bias       P x 3: delta_em, delta_wls, var_delta.
//   n_iters    P integer: EM iterations.
//   n_used     P integer: features in each covariate's EM.
//   theta      n sampling fractions (ANCOM-BC2), else NULL.
//   rss        D residual sums of squares of the fit whose variances
//              these are (NaN for failed fits).
//   mean       D mean clr values (rclr with missing values).
//======================================================
SEXP C_ancombc(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_pc, SEXP sexp_v2,
    SEXP sexp_X, SEXP sexp_Xinv, SEXP sexp_gsum, SEXP sexp_factors,
    SEXP sexp_intercept, SEXP sexp_max_iters, SEXP sexp_tol,
    SEXP sexp_n_threads ) {

  double pc        = asReal(sexp_pc);
  int    v2        = asLogical(sexp_v2);
  int    intercept = asInteger(sexp_intercept);
  int    n_thr     = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  init_n_ptrs(40);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n = em->n_samples, D = em->n_otus, P = ncols(sexp_X);
  int nb = (D + AB_BLK - 1) / AB_BLK;

  memset(&ab, 0, sizeof(ab));
  ab.n = n; ab.D = D; ab.P = P; ab.nb = nb;
  ab.X    = REAL(sexp_X);
  ab.Xinv = REAL(sexp_Xinv);
  ab.gsum = REAL(sexp_gsum);

  double *H = (double *) R_alloc((size_t) n * P, sizeof(double));
  for (int s = 0; s < n; s++)
    for (int p = 0; p < P; p++)
      H[s + (size_t) p * n] = ab.Xinv[p + (size_t) s * P];
  ab.H = H;


  //--------------------------------------------------------
  // 1. The logs (_transform_data()).
  //--------------------------------------------------------
  double *data = (double *) safe_malloc((size_t) n * D * sizeof(double));
  ab.data = data;

  int miss_mode = (pc == 0);   // R refuses zeros with pc 0 for ANCOM-BC.
  double lz = miss_mode ? R_NaN : log(pc);
  for (size_t i = 0; i < (size_t) n * D; i++) data[i] = lz;
  int any_missing = 0;
  for (int s = 0; s < n; s++)
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++) {
      double v = em->val_vec[e];
      data[s + (size_t) em->otu_vec[e] * n] = (miss_mode && v == 0) ? R_NaN : log(v + pc);
    }

  int *nmiss = (int *) R_alloc((size_t) D, sizeof(int));
  for (int f = 0; f < D; f++) {
    const double *y = data + (size_t) f * n;
    int c = 0;
    for (int s = 0; s < n; s++) if (ISNAN(y[s])) c++;
    nmiss[f] = c;
    if (c) any_missing = 1;
  }
  ab.nmiss  = nmiss;
  ab.sparse = any_missing;

  if (any_missing && !v2) {
    free_all();
    errorcall(R_NilValue, "`counts` has zeros, and `pseudocount = 0` leaves their logs undefined.");
  }
  if (any_missing)
    for (int s = 0; s < n; s++) {
      int obs = 0;
      for (int f = 0; f < D && !obs; f++) if (!ISNAN(data[s + (size_t) f * n])) obs = 1;
      if (!obs) {
        free_all();
        errorcall(R_NilValue,
          "A sample has only zeros: with zeros as missing values (`pseudocount = 0`), "
          "every sample needs a non-zero value.");
      }
    }

  // Each feature's mean clr (rclr with missing values): the mean over its
  // observed samples of its log minus the sample's mean log over the
  // features observed in it.
  SEXP sexp_mclr = PROTECT(allocVector(REALSXP, D));
  {
    double *smean = (double *) R_alloc((size_t) n + 1, sizeof(double));
    for (int s = 0; s < n; s++) {
      double t = 0; int c = 0;
      for (int f = 0; f < D; f++) {
        double v = data[s + (size_t) f * n];
        if (!ISNAN(v)) { t += v; c++; }
      }
      smean[s] = c ? t / c : R_NaN;
    }
    for (int f = 0; f < D; f++) {
      const double *y = data + (size_t) f * n;
      double t = 0; int c = 0;
      for (int s = 0; s < n; s++) if (!ISNAN(y[s])) { t += y[s] - smean[s]; c++; }
      REAL(sexp_mclr)[f] = c ? t / c : R_NaN;
    }
  }

  if (v2)
    for (int f = 0; f < D; f++) {
      double *y = data + (size_t) f * n, sum = 0;
      int obs = n - nmiss[f];
      for (int s = 0; s < n; s++) if (!ISNAN(y[s])) sum += y[s];
      double mean = sum / (obs > 0 ? obs : 1);
      for (int s = 0; s < n; s++) y[s] -= mean;
    }


  //--------------------------------------------------------
  // Buffers.
  //--------------------------------------------------------
  int nt = threads_for(n_thr, (double) n * D * P * (any_missing ? 4 : 1), 2e5);
  if (nt > nb) nt = nb > 0 ? nb : 1;

  int scr = n + 2 * P;
  if (scr < n + P * P) scr = n + P * P;
  if (scr < n * P + P) scr = n * P + P;
  ab.scr_len = scr;
  ab.scratch = (double *) safe_malloc((size_t) nt * scr * sizeof(double));

  ab.beta     = (double *) safe_malloc((size_t) D * P * sizeof(double) + 1);
  ab.beta_new = (double *) safe_malloc((size_t) D * P * sizeof(double) + 1);
  ab.part     = (double *) safe_malloc((size_t) 2 * nb * n * sizeof(double) + 1);
  ab.part2    = (double *) safe_malloc((size_t) nb * sizeof(double) + 1);
  double *sums  = (double *) R_alloc((size_t) n + 1, sizeof(double));
  double *theta = (double *) R_alloc((size_t) n + 1, sizeof(double));


  //--------------------------------------------------------
  // Sparse: masked-design SVDs, R's failed fits and
  // feature-local references (_r_fit_info()).
  //--------------------------------------------------------
  int  *valid    = NULL;
  int   n_fac    = length(sexp_factors);
  int **rb_col   = NULL;      // per factor: D local reference columns (-1 none)
  unsigned char **present = NULL;  // per factor: D x L

  if (any_missing) {

    ab.V    = (double *)  safe_malloc((size_t) D * P * P * sizeof(double) + 1);
    ab.Sinv = (double *)  safe_malloc((size_t) D * P * sizeof(double) + 1);
    ab.rank = (int *)     R_alloc((size_t) D, sizeof(int));
    ab.pinv = (double **) safe_malloc((size_t) D * sizeof(double *) + 1);
    for (int f = 0; f < D; f++) ab.pinv[f] = NULL;

    run_parallel(svd_worker, nt, D);

    valid = (int *) R_alloc((size_t) D, sizeof(int));
    for (int f = 0; f < D; f++) valid[f] = nmiss[f] < n;

    rb_col  = (int **)           R_alloc((size_t) n_fac + 1, sizeof(int *));
    present = (unsigned char **) R_alloc((size_t) n_fac + 1, sizeof(unsigned char *));

    for (int k = 0; k < n_fac; k++) {
      const int *fi       = INTEGER(VECTOR_ELT(sexp_factors, k));
      int        standard = fi[0], rebase = fi[1], L = fi[2];
      const int *lcol     = fi + 3;
      const int *codes    = fi + 3 + L;

      unsigned char *pr = (unsigned char *) R_alloc((size_t) D * L + 1, 1);
      memset(pr, 0, (size_t) D * L);
      for (int f = 0; f < D; f++) {
        const double *y = data + (size_t) f * n;
        for (int s = 0; s < n; s++) if (!ISNAN(y[s])) pr[f * (size_t) L + codes[s]] = 1;
      }
      present[k] = pr;

      int *rc = (int *) R_alloc((size_t) D, sizeof(int));
      for (int f = 0; f < D; f++) {
        int np = 0;
        for (int l = 0; l < L; l++) np += pr[f * (size_t) L + l];
        if (np < 2) valid[f] = 0;
        rc[f] = -1;
        if (standard && rebase && intercept >= 0 && !pr[f * (size_t) L + 0])
          for (int l = 1; l < L; l++)
            if (pr[f * (size_t) L + l]) { rc[f] = lcol[l]; break; }
      }
      rb_col[k] = rc;
    }

    int all_valid = 1;
    for (int f = 0; f < D; f++) if (!valid[f]) { all_valid = 0; break; }
    ab.valid = all_valid ? NULL : valid;
  }

  if (ab.failed) { free_all(); errorcall(R_NilValue, "Out of memory."); }


  //--------------------------------------------------------
  // A fit (initial: biased = 1; ANCOM-BC2's final: 0).
  //--------------------------------------------------------
  #define REBASE() if (any_missing) {                                          \
    for (int k = 0; k < n_fac; k++) {                                          \
      const int *fi = INTEGER(VECTOR_ELT(sexp_factors, k));                    \
      int L = fi[2]; const int *lcol = fi + 3;                                 \
      for (int f = 0; f < D; f++) {                                            \
        int c = rb_col[k][f];                                                  \
        if (c < 0 || (ab.valid && !ab.valid[f])) continue;                     \
        double *bt = ab.beta + (size_t) f * P, shift = bt[c];                  \
        bt[intercept] += shift;                                                \
        for (int l = 1; l < L; l++) {                                          \
          if (present[k][f * (size_t) L + l]) bt[lcol[l]] -= shift;            \
          else                                bt[lcol[l]]  = 0;                \
        }                                                                      \
        bt[c] = 0;                                                             \
      }                                                                        \
    }                                                                          \
  }

  #define FIT(biased) {                                                        \
    ab.theta = NULL;                                                           \
    fit_pass(nt, sums);                                                        \
    swap_beta();                                                               \
    if (!(biased)) {                                                           \
      for (int s = 0; s < n; s++) theta[s] = 0;                                \
    }                                                                          \
    else if (!any_missing) {                                                   \
      for (int s = 0; s < n; s++) theta[s] = sums[s] / D;                      \
    }                                                                          \
    else {                                                                     \
      int *nobs = (int *) R_alloc((size_t) n, sizeof(int));                   \
      for (int s = 0; s < n; s++) {                                            \
        int c = 0;                                                             \
        for (int f = 0; f < D; f++)                                            \
          if ((!ab.valid || ab.valid[f]) && !ISNAN(data[s + (size_t) f * n])) c++; \
        nobs[s] = c;                                                           \
        theta[s] = 0;                                                          \
      }                                                                        \
      double eps = 100.0;                                                      \
      int    it  = 0;                                                          \
      while (eps > AB_SP_TOL && it < AB_SP_ITERS) {                            \
        ab.theta = theta;                                                      \
        fit_pass(nt, sums);                                                    \
        double bd = 0, td = 0;                                                 \
        for (int b = 0; b < nb; b++) bd += ab.part2[b];                        \
        for (int s = 0; s < n; s++) {                                          \
          double tn = (nobs[s] == 0) ? R_NaN : sums[s] / nobs[s];              \
          double d  = tn - theta[s];                                           \
          d = d * d;                                                           \
          if (!ISNAN(d)) td += d;                                              \
          theta[s] = tn;                                                       \
        }                                                                      \
        eps = sqrt(bd + td);                                                   \
        swap_beta();                                                           \
        it++;                                                                  \
      }                                                                        \
    }                                                                          \
    REBASE();                                                                  \
    ab.theta = theta;                                                          \
    run_parallel(var_worker, nt, D);                                           \
  }


  SEXP sexp_var  = PROTECT(allocMatrix(REALSXP, D, P));
  SEXP sexp_cov  = PROTECT(alloc3DArray(REALSXP, D, P, P));
  SEXP sexp_rss  = PROTECT(allocVector(REALSXP, D));
  ab.var = REAL(sexp_var);
  ab.rss = REAL(sexp_rss);
  for (size_t i = 0; i < (size_t) D * P * P; i++) REAL(sexp_cov)[i] = R_NaN;


  //--------------------------------------------------------
  // 2. The initial fit. ANCOM-BC keeps its covariances.
  //--------------------------------------------------------
  ab.cov = v2 ? NULL : REAL(sexp_cov);
  FIT(1);


  //--------------------------------------------------------
  // 3. The bias EM, one covariate per task.
  //--------------------------------------------------------
  SEXP sexp_bias  = PROTECT(allocMatrix(REALSXP, P, 3));
  SEXP sexp_iters = PROTECT(allocVector(INTSXP, P));
  SEXP sexp_used  = PROTECT(allocVector(INTSXP, P));
  double *bias3 = (double *) R_alloc((size_t) 3 * P, sizeof(double));
  for (int i = 0; i < 3 * P; i++) bias3[i] = R_NaN;

  em_beta_all = ab.beta;
  em_var_all  = ab.var;
  em_out      = bias3;
  em_iters    = INTEGER(sexp_iters);
  em_m        = INTEGER(sexp_used);
  em_max_iter = asInteger(sexp_max_iters);
  em_tol      = asReal(sexp_tol);

  int nte = threads_for(n_thr < P ? n_thr : P, (double) D * P * 1e3, 1e6);
  run_parallel(em_worker, nte, nte > 1 ? 100 : 1);
  if (ab.failed) { free_all(); errorcall(R_NilValue, "Out of memory."); }

  for (int i = 0; i < P; i++)
    for (int j = 0; j < 3; j++)
      REAL(sexp_bias)[i + (size_t) j * P] = bias3[3 * i + j];


  //--------------------------------------------------------
  // 4. ANCOM-BC2: sampling fractions, then the final fit.
  //--------------------------------------------------------
  SEXP sexp_theta = R_NilValue;
  int  em_ok = 1;
  for (int i = 0; i < P; i++) if (em_m[i] == 0) em_ok = 0;

  if (v2 && em_ok) {

    double *delta = (double *) R_alloc((size_t) P, sizeof(double));
    for (int i = 0; i < P; i++) delta[i] = bias3[3 * i];
    ab.delta = delta;

    sexp_theta = PROTECT(allocVector(REALSXP, n));
    double *th = REAL(sexp_theta);

    if (!any_missing) {
      // mean(data) - X mean(beta) + X delta
      double *mb = (double *) R_alloc((size_t) P, sizeof(double));
      for (int p = 0; p < P; p++) {
        double t = 0;
        for (int f = 0; f < D; f++) t += ab.beta[p + (size_t) f * P];
        mb[p] = t / D;
      }
      for (int s = 0; s < n; s++) {
        double t = 0;
        for (int f = 0; f < D; f++) t += data[s + (size_t) f * n];
        double a = t / D, xm = 0, xd = 0;
        for (int p = 0; p < P; p++) {
          xm += ab.X[s + (size_t) p * n] * mb[p];
          xd += ab.X[s + (size_t) p * n] * delta[p];
        }
        th[s] = a - xm + xd;
      }
    }
    else {
      run_parallel(sf_worker, nt, D);
      for (int s = 0; s < n; s++) {
        double t = 0, c = 0;
        for (int b = 0; b < nb; b++) {
          t += ab.part[s + (size_t) b * n];
          c += ab.part[s + (size_t) (nb + b) * n];
        }
        th[s] = (c > 0) ? t / c : R_NaN;
      }
    }

    for (int f = 0; f < D; f++) {
      double *y = data + (size_t) f * n;
      for (int s = 0; s < n; s++) y[s] -= th[s];
    }

    ab.cov = REAL(sexp_cov);
    FIT(0);
  }
  else {
    sexp_theta = PROTECT(R_NilValue);
  }

  #undef FIT
  #undef REBASE


  //--------------------------------------------------------
  // Results.
  //--------------------------------------------------------
  SEXP sexp_beta = PROTECT(allocMatrix(REALSXP, D, P));
  for (int f = 0; f < D; f++)
    for (int p = 0; p < P; p++)
      REAL(sexp_beta)[f + (size_t) p * D] = ab.beta[p + (size_t) f * P];

  SEXP sexp_valid = R_NilValue, sexp_rank = R_NilValue, sexp_nmiss = R_NilValue;
  if (any_missing) {
    if (ab.valid) {
      sexp_valid = allocVector(LGLSXP, D);
      for (int f = 0; f < D; f++) LOGICAL(sexp_valid)[f] = valid[f];
    }
    PROTECT(sexp_valid);
    sexp_rank = PROTECT(allocVector(INTSXP, D));
    for (int f = 0; f < D; f++) INTEGER(sexp_rank)[f] = ab.rank[f];
    sexp_nmiss = PROTECT(allocVector(INTSXP, D));
    for (int f = 0; f < D; f++) INTEGER(sexp_nmiss)[f] = nmiss[f];
    for (int f = 0; f < D; f++) if (ab.pinv[f]) free(ab.pinv[f]);
  }
  else {
    PROTECT(sexp_valid); PROTECT(sexp_rank); PROTECT(sexp_nmiss);
  }

  const char *nms[] = {
    "beta", "var", "vcov", "valid", "rank", "n_missing", "bias", "n_iters",
    "n_used", "theta", "rss", "mean" };
  SEXP vals[] = {
    sexp_beta, sexp_var, sexp_cov, sexp_valid, sexp_rank, sexp_nmiss,
    sexp_bias, sexp_iters, sexp_used, sexp_theta, sexp_rss, sexp_mclr };

  SEXP result = PROTECT(allocVector(VECSXP, 12));
  SEXP names  = PROTECT(allocVector(STRSXP, 12));
  for (int i = 0; i < 12; i++) {
    SET_VECTOR_ELT(result, i, vals[i]);
    SET_STRING_ELT(names, i, mkChar(nms[i]));
  }
  setAttrib(result, R_NamesSymbol, names);

  free_all();
  UNPROTECT(14);
  return result;
}



//======================================================================
// Post-hoc tests (ancombc_global(), ancombc_dunnett(), ancombc_trend())
// and structural_zeros(): scikit-bio's _global_stats(), _dunn_global(),
// _trend_test(), _constrain_est(), _constrain_est_identity(),
// _safe_inverse_spd(), _l_infty() and struc_zero(). The p-values,
// adjustments and mdFDR rules are in R/ancombc_posthoc.r.
//
// Each feature is one task, with its own PCG stream (seed, feature), so
// results are identical for any number of threads, and a feature's
// bootstrap doesn't depend on the other features.
//======================================================================

// Moore-Penrose inverse of a k x k matrix A (numpy.linalg.pinv: SVD, cutoff
// 1e-15 times the largest singular value). w: 2 k k + k doubles.
static void small_pinv (const double *A, int k, double *out, double *w) {
  double *U = w, *V = w + k * k, *S = w + 2 * k * k;
  for (int i = 0; i < k * k; i++) U[i] = A[i];
  jacobi_svd(U, k, k, V, S);
  double cut = 1e-15 * S[0];
  for (int p = 0; p < k; p++)
    for (int q = 0; q < k; q++) {
      double t = 0;
      for (int j = 0; j < k; j++)
        if (S[j] > cut) t += V[p + j * k] * (U[q + j * k] / S[j]) / S[j];
      out[p + q * k] = t;
    }
}


// Cholesky factorisation of the symmetric positive definite k x k matrix A
// (lower triangle into L). Returns 0 if A isn't positive definite.
static int cholesky (const double *A, int k, double *L) {
  for (int i = 0; i < k * k; i++) L[i] = 0;
  for (int j = 0; j < k; j++) {
    double s = A[j + j * k];
    for (int m = 0; m < j; m++) s -= L[j + m * k] * L[j + m * k];
    if (!(s > 0)) return 0;
    double d = sqrt(s);
    L[j + j * k] = d;
    for (int i = j + 1; i < k; i++) {
      double t = A[i + j * k];
      for (int m = 0; m < j; m++) t -= L[i + m * k] * L[j + m * k];
      L[i + j * k] = t / d;
    }
  }
  return 1;
}

// Inverse from a Cholesky factor. w: k doubles.
static void chol_inverse (const double *L, int k, double *out, double *w) {
  for (int c = 0; c < k; c++) {
    // Solve L y = e_c, then L' x = y.
    for (int i = 0; i < k; i++) {
      double t = (i == c) ? 1 : 0;
      for (int m = 0; m < i; m++) t -= L[i + m * k] * w[m];
      w[i] = t / L[i + i * k];
    }
    for (int i = k - 1; i >= 0; i--) {
      double t = w[i];
      for (int m = i + 1; m < k; m++) t -= L[m + i * k] * out[m + c * k];
      out[i + c * k] = t / L[i + i * k];
    }
  }
}

// Gauss-Jordan inverse with partial pivoting; returns 0 if singular.
// w: k k doubles.
static int gj_inverse (const double *A, int k, double *out, double *w) {
  for (int i = 0; i < k * k; i++) { w[i] = A[i]; out[i] = 0; }
  for (int i = 0; i < k; i++) out[i + i * k] = 1;
  for (int c = 0; c < k; c++) {
    int piv = c;
    for (int r = c + 1; r < k; r++) if (fabs(w[r + c * k]) > fabs(w[piv + c * k])) piv = r;
    if (w[piv + c * k] == 0) return 0;
    if (piv != c)
      for (int j = 0; j < k; j++) {
        double t = w[c + j * k]; w[c + j * k] = w[piv + j * k]; w[piv + j * k] = t;
        t = out[c + j * k]; out[c + j * k] = out[piv + j * k]; out[piv + j * k] = t;
      }
    double d = w[c + c * k];
    for (int j = 0; j < k; j++) { w[c + j * k] /= d; out[c + j * k] /= d; }
    for (int r = 0; r < k; r++) {
      if (r == c) continue;
      double f = w[r + c * k];
      if (f == 0) continue;
      for (int j = 0; j < k; j++) { w[r + j * k] -= f * w[c + j * k]; out[r + j * k] -= f * out[c + j * k]; }
    }
  }
  return 1;
}


//======================================================
// A small dense convex QP with homogeneous inequality
// constraints, by a primal active-set method (Nocedal &
// Wright 2006, Algorithm 16.3), starting from x = 0,
// which is always feasible:
//
//   minimise 0.5 x' G x + c' x   subject to  C x >= 0
//
// G (k x k, positive definite) is given with its inverse
// Gi; C is m x k (column-major). The minimiser is unique,
// so it is the solution SLSQP (scikit-bio) and
// quadprog::solve.QP (ANCOMBC) approximate. Returns 1 on
// success (x set), 0 if it didn't converge.
// w: qp_wlen(m, k) doubles, iw: 2 m ints.
//======================================================
static int qp_solve (
    const double *G, const double *Gi, const double *c, const double *C,
    int m, int k, double *x, double *w, int *iw ) {

  double *g   = w;            // k
  double *p   = g + k;        // k
  double *gig = p + k;        // k: Gi g
  double *t   = gig + k;      // k
  double *lam = t + k;        // m
  double *r   = lam + m;      // m
  double *S   = r + m;        // m x m
  double *Si  = S + m * m;    // m x m
  double *sw  = Si + m * m;   // 2 m m + m: small_pinv() scratch
  int    *act = iw;           // working set, nW entries
  int    *inW = iw + m;       // m flags
  int     nW  = 0;

  // Scale for the tolerances: the unconstrained solution's size.
  double scale = 0;
  for (int i = 0; i < k; i++) {
    double s = 0;
    for (int j = 0; j < k; j++) s += Gi[i + j * k] * c[j];
    scale += s * s;
  }
  scale = sqrt(scale) + 1e-300;

  for (int i = 0; i < k; i++) x[i] = 0;
  for (int i = 0; i < m; i++) inW[i] = 0;

  for (int iter = 0; iter < 50 + 10 * (m + k); iter++) {

    for (int i = 0; i < k; i++) {
      double s = c[i];
      for (int j = 0; j < k; j++) s += G[i + j * k] * x[j];
      g[i] = s;
    }
    for (int i = 0; i < k; i++) {
      double s = 0;
      for (int j = 0; j < k; j++) s += Gi[i + j * k] * g[j];
      gig[i] = s;
    }

    // The equality-constrained step: p = Gi (A' lam - g), with
    // (A Gi A') lam = A Gi g, for the working set's rows A.
    if (nW == 0) {
      for (int i = 0; i < k; i++) p[i] = -gig[i];
    }
    else {
      for (int a = 0; a < nW; a++) {
        const int ra = act[a];
        double s = 0;
        for (int j = 0; j < k; j++) s += C[ra + j * m] * gig[j];
        r[a] = s;
        for (int b = 0; b < nW; b++) {
          const int rb = act[b];
          double v = 0;
          for (int i = 0; i < k; i++) {
            double gi_cb = 0;
            for (int j = 0; j < k; j++) gi_cb += Gi[i + j * k] * C[rb + j * m];
            v += C[ra + i * m] * gi_cb;
          }
          S[a + b * nW] = v;
        }
      }
      small_pinv(S, nW, Si, sw);
      for (int a = 0; a < nW; a++) {
        double s = 0;
        for (int b = 0; b < nW; b++) s += Si[a + b * nW] * r[b];
        lam[a] = s;
      }
      for (int i = 0; i < k; i++) {
        double s = -g[i];
        for (int a = 0; a < nW; a++) s += C[act[a] + i * m] * lam[a];
        t[i] = s;
      }
      for (int i = 0; i < k; i++) {
        double s = 0;
        for (int j = 0; j < k; j++) s += Gi[i + j * k] * t[j];
        p[i] = s;
      }
    }

    double pn = 0, xn = 0;
    for (int i = 0; i < k; i++) { pn += p[i] * p[i]; xn += x[i] * x[i]; }
    pn = sqrt(pn);

    if (pn <= 1e-12 * (scale + sqrt(xn))) {
      if (nW == 0) return 1;
      int    worst = -1;
      double lmin  = -1e-12 * (1 + fabs(lam[0]));
      for (int a = 0; a < nW; a++) if (lam[a] < lmin) { lmin = lam[a]; worst = a; }
      if (worst < 0) return 1;
      inW[act[worst]] = 0;
      for (int a = worst; a < nW - 1; a++) act[a] = act[a + 1];
      nW--;
      continue;
    }

    double alpha = 1;
    int    block = -1;
    for (int i = 0; i < m; i++) {
      if (inW[i]) continue;
      double ap = 0, ax = 0;
      for (int j = 0; j < k; j++) { ap += C[i + j * m] * p[j]; ax += C[i + j * m] * x[j]; }
      if (ap < -1e-14 * pn) {
        double st = -ax / ap;
        if (st < 0) st = 0;
        if (st < alpha) { alpha = st; block = i; }
      }
    }
    for (int i = 0; i < k; i++) x[i] += alpha * p[i];
    if (block >= 0) { inW[block] = 1; act[nW++] = block; }
  }

  return 0;
}

static int qp_wlen (int m, int k) { return 4 * k + 3 * m + 4 * m * m + 4; }



//======================================================
// R interface: the global test's quadratic forms,
// W_f = b_f' pinv(V_f) b_f, over the k grouping
// coefficients (beta: D x k; vcov: D x k x k). NaN where
// `valid` (logical D, or NULL) is FALSE.
//======================================================
static const double *ph_beta, *ph_vcov, *ph_var, *ph_W, *ph_dof, *ph_draws;
static const int    *ph_valid;
static int           ph_D, ph_k, ph_B, ph_npat;
static uint64_t      ph_seed;
static double       *ph_out, *ph_out2;
static int          *ph_iout;
static const double **ph_pmat;   // each pattern's m_i x k matrix
static const int    *ph_pm;      // its m_i
static const int    *ph_nodes;
static double       *ph_scr;
static int           ph_scr_len;

static void *global_worker (void *arg) {
  worker_t *wk = (worker_t *) arg;
  int D = ph_D, k = ph_k;
  double *V  = ph_scr + (size_t) wk->i * ph_scr_len;
  double *Vi = V + k * k, *w = Vi + k * k;
  for (int f = wk->i; f < D; f += wk->n) {
    if (ph_valid && !ph_valid[f]) { ph_out[f] = R_NaN; continue; }
    for (int p = 0; p < k; p++)
      for (int q = 0; q < k; q++)
        V[p + q * k] = ph_vcov[f + (size_t) p * D + (size_t) q * D * k];
    small_pinv(V, k, Vi, w);
    double s = 0;
    for (int p = 0; p < k; p++) {
      double t = 0;
      for (int q = 0; q < k; q++) t += Vi[p + q * k] * ph_beta[f + (size_t) q * D];
      s += ph_beta[f + (size_t) p * D] * t;
    }
    ph_out[f] = s;
  }
  return NULL;
}

SEXP C_ancombc_global (SEXP sexp_beta, SEXP sexp_vcov, SEXP sexp_valid, SEXP sexp_n_threads) {
  ph_D = nrows(sexp_beta); ph_k = ncols(sexp_beta);
  ph_beta = REAL(sexp_beta); ph_vcov = REAL(sexp_vcov);
  ph_valid = isNull(sexp_valid) ? NULL : LOGICAL(sexp_valid);
  SEXP res = PROTECT(allocVector(REALSXP, ph_D));
  ph_out = REAL(res);
  int nt = threads_for(asInteger(sexp_n_threads), (double) ph_D * ph_k * ph_k * ph_k * 20, 1e5);
  ph_scr_len = 4 * ph_k * ph_k + ph_k;
  ph_scr = (double *) R_alloc((size_t) nt * ph_scr_len, sizeof(double));
  run_parallel(global_worker, nt, ph_D);
  UNPROTECT(1);
  return res;
}



//======================================================
// R interface: Dunnett's global screen (_dunn_global()),
// the number of bootstraps whose max |W_null| over the k
// comparisons exceeds the observed max |W|, per feature.
// W_null is t(dof_f) (dof NaN: 999), or N(0, 1) when dof
// is NULL; draws (k x B values per feature, as a D x k x B
// array) replace the sampler when given (tests).
//======================================================
static void *dunnett_worker (void *arg) {
  worker_t *wk = (worker_t *) arg;
  int D = ph_D, k = ph_k, B = ph_B;
  for (int f = wk->i; f < D; f += wk->n) {
    double obs = R_NegInf;
    int    nan = 0;
    for (int j = 0; j < k; j++) {
      double a = fabs(ph_W[f + (size_t) j * D]);
      if (ISNAN(a)) nan = 1; else if (a > obs) obs = a;
    }
    rng_t rng;
    rng_seed(&rng, ph_seed, (uint64_t) f);
    double nu = ph_dof ? (ISNAN(ph_dof[f]) ? 999.0 : ph_dof[f]) : 0;
    int exceed = 0;
    for (int b = 0; b < B; b++) {
      double mx = R_NegInf;
      for (int j = 0; j < k; j++) {
        double z;
        if (ph_draws) z = ph_draws[f + (size_t) j * D + (size_t) b * D * k];
        else if (ph_dof) z = rng_norm(&rng) / sqrt(2 * rng_gamma(&rng, nu / 2) / nu);
        else z = rng_norm(&rng);
        z = fabs(z);
        if (ISNAN(z) || z > mx) mx = z;
      }
      // numpy: max propagates NaN, and a NaN comparison is FALSE.
      if (!nan && !ISNAN(mx) && mx > obs) exceed++;
    }
    ph_iout[f] = exceed;
  }
  return NULL;
}

SEXP C_ancombc_dunnett (
    SEXP sexp_W, SEXP sexp_dof, SEXP sexp_B, SEXP sexp_seed, SEXP sexp_draws,
    SEXP sexp_n_threads ) {
  ph_D = nrows(sexp_W); ph_k = ncols(sexp_W);
  ph_W = REAL(sexp_W);
  ph_dof = isNull(sexp_dof) ? NULL : REAL(sexp_dof);
  ph_B = asInteger(sexp_B);
  ph_seed = (uint64_t) (uint32_t) asInteger(sexp_seed);
  ph_draws = isNull(sexp_draws) ? NULL : REAL(sexp_draws);
  SEXP res = PROTECT(allocVector(INTSXP, ph_D));
  ph_iout = INTEGER(res);
  int nt = threads_for(asInteger(sexp_n_threads), (double) ph_D * ph_k * ph_B * 50, 1e5);
  run_parallel(dunnett_worker, nt, ph_D);
  UNPROTECT(1);
  return res;
}



//======================================================
// R interface: the trend test (_trend_test()).
//
//   beta, var  D x k; vcov D x k x k (the grouping block).
//   valid      logical D, or NULL.
//   patterns   list of m_i x k constraint matrices (C x >= 0).
//   nodes      0-based node of each pattern.
//   draws      NULL, or a D x k x B array of N(0, 1) null draws.
//
// Returns list(W (D), exceed (D, integer), beta (D x k: the
// constrained estimates of the winning pattern), pattern (D,
// 1-based)). For each valid feature: P = V^-1 (Cholesky, or
// (V + 1e-8 I)^-1), G = (2P + 2P')/2 + 1e-10 I, and for each
// pattern x = argmin 0.5 x'Gx - (2Pb)'x s.t. C x >= 0 (zeros
// if the solver fails); l = max(|x_node|, |x_node - x_k|);
// W = max l (first on ties). Each bootstrap draws z ~ N(0, I_k),
// projects it onto each pattern's cone (identity metric),
// scales it by sqrt(var) (var NaN: 1), takes the same l, and
// counts max l > W.
//======================================================
static void *trend_worker (void *arg) {

  worker_t *wk = (worker_t *) arg;
  int D = ph_D, k = ph_k, B = ph_B, np = ph_npat;
  int mmax = 0;
  for (int t = 0; t < np; t++) if (ph_pm[t] > mmax) mmax = ph_pm[t];

  double *w   = ph_scr + (size_t) wk->i * ph_scr_len;
  double *V   = w;             // k k
  double *P   = V + k * k;     // k k
  double *G   = P + k * k;     // k k
  double *Gi  = G + k * k;     // k k
  double *L   = Gi + k * k;    // k k
  double *I   = L + k * k;     // k k (identity)
  double *c   = I + k * k;     // k
  double *x   = c + k;         // k
  double *xb  = x + k;         // k: best pattern's x
  double *z   = xb + k;        // k
  double *sd  = z + k;         // k
  double *qw  = sd + k;        // qp work
  int    *iw  = (int *) (qw + qp_wlen(mmax, k) + 4 * k * k + k);

  for (int i = 0; i < k * k; i++) I[i] = 0;
  for (int i = 0; i < k; i++) I[i + i * k] = 1;

  for (int f = wk->i; f < D; f += wk->n) {

    if (ph_valid && !ph_valid[f]) {
      ph_out[f] = R_NaN; ph_iout[f] = 0; ph_iout[D + f] = NA_INTEGER;
      for (int j = 0; j < k; j++) ph_out2[f + (size_t) j * D] = R_NaN;
      continue;
    }

    // _safe_inverse_spd()
    for (int p = 0; p < k; p++)
      for (int q = 0; q < k; q++) {
        double a = ph_vcov[f + (size_t) p * D + (size_t) q * D * k];
        double b = ph_vcov[f + (size_t) q * D + (size_t) p * D * k];
        V[p + q * k] = (a + b) / 2;
      }
    if (cholesky(V, k, L)) chol_inverse(L, k, P, qw);
    else {
      for (int i = 0; i < k; i++) V[i + i * k] += 1e-8;
      if (!gj_inverse(V, k, P, qw)) for (int i = 0; i < k * k; i++) P[i] = R_NaN;
    }

    // Dmat, dvec
    for (int p = 0; p < k; p++)
      for (int q = 0; q < k; q++)
        G[p + q * k] = (2 * P[p + q * k] + 2 * P[q + p * k]) / 2 + (p == q ? 1e-10 : 0);
    for (int p = 0; p < k; p++) {
      double s = 0;
      for (int q = 0; q < k; q++) s += 2 * P[p + q * k] * ph_beta[f + (size_t) q * D];
      c[p] = -s;
    }
    int gok = cholesky(G, k, L);
    if (gok) chol_inverse(L, k, Gi, qw);

    // numpy.argmax(): the first largest (a NaN wins, as the first NaN).
    double best = R_NegInf;
    int    bpat = 0;
    for (int t = 0; t < np; t++) {
      if (!gok || !qp_solve(G, Gi, c, ph_pmat[t], ph_pm[t], k, x, qw, iw))
        for (int j = 0; j < k; j++) x[j] = 0;
      int nd = ph_nodes[t];
      double l1 = fabs(x[nd]), l2 = fabs(x[nd] - x[k - 1]);
      double l  = (l2 > l1) ? l2 : l1;     // Python's max(l1, l2)
      if (t == 0 || (!ISNAN(best) && (l > best || ISNAN(l)))) {
        best = l; bpat = t;
        for (int j = 0; j < k; j++) xb[j] = x[j];
      }
    }
    ph_out[f] = best;
    ph_iout[D + f] = bpat + 1;
    for (int j = 0; j < k; j++) ph_out2[f + (size_t) j * D] = xb[j];

    // Bootstrap.
    for (int j = 0; j < k; j++) {
      double v = ph_var[f + (size_t) j * D];
      if (ISNAN(v)) v = 1;
      sd[j] = sqrt(v > 0 ? v : 0);
    }
    rng_t rng;
    rng_seed(&rng, ph_seed, (uint64_t) f);
    int exceed = 0;
    for (int b = 0; b < B; b++) {
      for (int j = 0; j < k; j++)
        z[j] = ph_draws ? ph_draws[f + (size_t) j * D + (size_t) b * D * k] : rng_norm(&rng);
      for (int j = 0; j < k; j++) c[j] = -z[j];
      double lmax = R_NegInf;
      for (int t = 0; t < np; t++) {
        if (!qp_solve(I, I, c, ph_pmat[t], ph_pm[t], k, x, qw, iw))
          for (int j = 0; j < k; j++) x[j] = 0;
        int nd = ph_nodes[t];
        double a = x[nd] * sd[nd], e = x[k - 1] * sd[k - 1];
        double l1 = fabs(a), l2 = fabs(a - e);
        double l  = (ISNAN(l1) || ISNAN(l2)) ? R_NaN : (l1 >= l2 ? l1 : l2);
        if (ISNAN(l) || l > lmax) lmax = l;
      }
      if (lmax > best) exceed++;            // FALSE for NaN
    }
    ph_iout[f] = exceed;
  }

  return NULL;
}

SEXP C_ancombc_trend (
    SEXP sexp_beta, SEXP sexp_var, SEXP sexp_vcov, SEXP sexp_valid,
    SEXP sexp_patterns, SEXP sexp_nodes, SEXP sexp_B, SEXP sexp_seed,
    SEXP sexp_draws, SEXP sexp_n_threads ) {

  ph_D = nrows(sexp_beta); ph_k = ncols(sexp_beta);
  ph_beta = REAL(sexp_beta); ph_var = REAL(sexp_var); ph_vcov = REAL(sexp_vcov);
  ph_valid = isNull(sexp_valid) ? NULL : LOGICAL(sexp_valid);
  ph_npat = length(sexp_patterns);
  ph_pmat = (const double **) R_alloc((size_t) ph_npat + 1, sizeof(double *));
  int *pm = (int *) R_alloc((size_t) ph_npat + 1, sizeof(int));
  for (int t = 0; t < ph_npat; t++) {
    ph_pmat[t] = REAL(VECTOR_ELT(sexp_patterns, t));
    pm[t]      = nrows(VECTOR_ELT(sexp_patterns, t));
  }
  ph_pm = pm;
  ph_nodes = INTEGER(sexp_nodes);
  ph_B = asInteger(sexp_B);
  ph_seed = (uint64_t) (uint32_t) asInteger(sexp_seed);
  ph_draws = isNull(sexp_draws) ? NULL : REAL(sexp_draws);

  int mmax = 0;
  for (int t = 0; t < ph_npat; t++) if (ph_pm[t] > mmax) mmax = ph_pm[t];

  SEXP sexp_W    = PROTECT(allocVector(REALSXP, ph_D));
  SEXP sexp_ex   = PROTECT(allocVector(INTSXP, 2 * (R_xlen_t) ph_D));
  SEXP sexp_beta2 = PROTECT(allocMatrix(REALSXP, ph_D, ph_k));
  ph_out = REAL(sexp_W); ph_iout = INTEGER(sexp_ex); ph_out2 = REAL(sexp_beta2);

  int k = ph_k;
  int nt = threads_for(asInteger(sexp_n_threads),
                       (double) ph_D * (ph_B + 1) * ph_npat * (k + mmax) * k * k * 4, 1e5);
  ph_scr_len = 6 * k * k + 5 * k + qp_wlen(mmax, k) + 4 * k * k + k + 2 * mmax + 8;
  ph_scr = (double *) R_alloc((size_t) nt * ph_scr_len, sizeof(double));
  run_parallel(trend_worker, nt, ph_D);

  SEXP res = PROTECT(allocVector(VECSXP, 4));
  SEXP nms = PROTECT(allocVector(STRSXP, 4));
  SET_VECTOR_ELT(res, 0, sexp_W);   SET_STRING_ELT(nms, 0, mkChar("W"));
  SET_VECTOR_ELT(res, 1, sexp_ex);  SET_STRING_ELT(nms, 1, mkChar("counts"));
  SET_VECTOR_ELT(res, 2, sexp_beta2); SET_STRING_ELT(nms, 2, mkChar("beta"));
  setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(5);
  return res;
}



//======================================================
// R interface: structural_zeros()'s counts: the number of
// non-zero values of each feature in each group (D x n_groups,
// integer); groups are 0-based codes per sample.
//======================================================
SEXP C_ancombc_nonzero (SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_groups, SEXP sexp_n_groups) {
  init_n_ptrs(20);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int D = em->n_otus, G = asInteger(sexp_n_groups);
  const int *g = INTEGER(sexp_groups);
  SEXP res = PROTECT(allocMatrix(INTSXP, D, G));
  int *out = INTEGER(res);
  for (size_t i = 0; i < (size_t) D * G; i++) out[i] = 0;
  for (int s = 0; s < em->n_samples; s++)
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++)
      if (em->val_vec[e] != 0) out[em->otu_vec[e] + (size_t) g[s] * D]++;
  free_all();
  UNPROTECT(1);
  return res;
}
