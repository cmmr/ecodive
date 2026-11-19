// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The discriminant analysis in this file (the within-group scaling, the
// SVDs of the sphered residuals and of the weighted group means, the rank
// tolerances, the weighted centre and the scores) is ported from
// scikit-learn 1.9.1: `LinearDiscriminantAnalysis._solve_svd()`, `fit()`
// and `transform()` in sklearn/discriminant_analysis.py. That file is
// covered by this license:
//
//
//   BSD 3-Clause License
//
//   Copyright (c) 2007-2026 The scikit-learn developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


/*
 * Linear discriminant analysis (lda()): Fisher's discriminant axes for
 * known groups, as scikit-learn's LinearDiscriminantAnalysis(solver =
 * "svd") computes them, on a dense n x p table X (samples x features) with
 * groups 0..g-1 and prior probabilities pi_k:
 *
 *   1. The group means m_k, and the within-group residuals
 *      Xc = X - m_grp (n x p).
 *   2. Each feature's within-group standard deviation (divisor n), s_j, or
 *      1 if it is 0. This only conditions the next step.
 *   3. The SVD of Z = Xc S^-1 / sqrt(n - g) = U D V'. With rank(Z) = p,
 *      W = S^-1 V D^-1 sphers the pooled within-group covariance:
 *      W' (Xc'Xc / (n - g)) W = I.
 *   4. The weighted centre xbar = sum_k pi_k m_k, and the g x p matrix
 *      B = diag(sqrt(n pi_k / (g - 1))) (M - 1 xbar') W of the sphered
 *      group means, and its SVD, B = U_b D_b V_b'.
 *   5. The discriminant coefficients ("scalings") W V_b, keeping the axes
 *      whose singular value is above tol times the first, and the scores
 *      (X - 1 xbar') W V_b.
 *
 * D_b^2 is the ratio of the between-group to the within-group mean square
 * on each axis.
 *
 * Step 3 divides by n - g, as MASS::lda() does, where scikit-learn divides
 * by n: so the scores have pooled within-group variance 1 (scikit-learn's
 * are larger by sqrt(n / (n - g))), and D_b^2 is the canonical F statistic
 * (scikit-learn's are larger by n / (n - g)). The proportions D_b^2 /
 * sum(D_b^2) are the same.
 *
 * The SVDs are computed as a Householder QR of Z, Z = Q R (parallel over
 * the columns), then one-sided (Hestenes) Jacobi rotations of the p x p R,
 * which give D and V to full relative accuracy; and one-sided Jacobi
 * rotations of B's p columns, which give D_b and V_b.
 *
 * lda_loo() classifies each sample by the fit to the other n - 1 samples
 * (as MASS::lda(CV = TRUE), with the priors held fixed), for the first
 * m = 1, ..., p columns of X at once, for capda(). See there.
 */

#include "ecodive.h"
#include <float.h> // DBL_EPSILON

#define LDA_TOL        1e-4  // scikit-learn's (and MASS's) default tol.
#define LDA_MAX_SWEEPS 100   // Jacobi sweeps; 10 or so are usual.



//======================================================
// Householder QR of the n x p A (column-major, n >= p),
// in place: R (upper triangle) is left in R (p x p).
// The reflection of each column onto later columns is
// done in parallel, each column on its own.
//======================================================

static double       *qr_A;
static const double *qr_v;
static int           qr_n, qr_p, qr_a;

static void *qr_worker(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = qr_n, a = qr_a;
  const double *v = qr_v;          // length n - a, v'v = 2.

  for (int j = a + 1 + thread_i; j < qr_p; j += n_threads) {
    double *col = qr_A + (size_t) j * n + a;
    double  s   = 0;
    for (int i = 0; i < n - a; i++) s += v[i] * col[i];
    for (int i = 0; i < n - a; i++) col[i] -= s * v[i];
  }
  return NULL;
}

static void lda_qr(double *A, int n, int p, double *R, double *v, int n_threads) {

  qr_A = A; qr_n = n; qr_p = p; qr_v = v;
  memset(R, 0, (size_t) p * p * sizeof(double));

  for (int a = 0; a < p; a++) {

    double *col = A + (size_t) a * n + a;
    int     len = n - a;

    double norm = 0;
    for (int i = 0; i < len; i++) norm += col[i] * col[i];
    norm = sqrt(norm);

    double r = (col[0] > 0) ? -norm : norm;  // R[a, a]
    double vn = 0;
    if (norm > 0) {
      v[0] = col[0] - r;
      for (int i = 1; i < len; i++) v[i] = col[i];
      for (int i = 0; i < len; i++) vn += v[i] * v[i];
    }

    if (vn > 0) {
      double f = sqrt(2 / vn);
      for (int i = 0; i < len; i++) v[i] *= f;
      qr_a = a;
      run_parallel(qr_worker, threads_for(n_threads, (double) len * (p - a - 1), 1e5), p - a - 1);
    }
    else {
      r = 0;
    }

    R[a + (size_t) a * p] = r;
    for (int j = a + 1; j < p; j++) R[a + (size_t) j * p] = A[a + (size_t) j * n];
  }
}



//======================================================
// One-sided Jacobi SVD of the m x c A (column-major), in
// place: on return A V = U D, with V (c x c) orthogonal,
// the singular values s (decreasing) the norms of A's
// columns, and A's and V's columns in that order. `tmp`
// holds m * c + c doubles. Returns 0 if converged.
//======================================================

static int lda_jacobi_svd(double *A, int m, int c, double *V, double *s, double *tmp) {

  memset(V, 0, (size_t) c * c * sizeof(double));
  for (int j = 0; j < c; j++) V[j + (size_t) j * c] = 1;

  // Columns at rounding level (A's rank is often below c) are left alone.
  double tiny = 0;
  for (size_t e = 0; e < (size_t) m * c; e++) tiny += A[e] * A[e];
  tiny *= DBL_EPSILON * DBL_EPSILON;

  int converged = 0;
  for (int sweep = 0; sweep < LDA_MAX_SWEEPS && !converged; sweep++) {
    converged = 1;
    for (int i = 0; i < c - 1; i++) {
      double *ai = A + (size_t) i * m;
      double *vi = V + (size_t) i * c;
      for (int j = i + 1; j < c; j++) {
        double *aj = A + (size_t) j * m;
        double *vj = V + (size_t) j * c;
        double alpha = 0, beta = 0, gamma = 0;
        for (int r = 0; r < m; r++) {
          alpha += ai[r] * ai[r];
          beta  += aj[r] * aj[r];
          gamma += ai[r] * aj[r];
        }
        if (alpha <= tiny || beta <= tiny) continue;
        if (fabs(gamma) <= DBL_EPSILON * sqrt(alpha * beta)) continue;
        converged = 0;
        double zeta = (beta - alpha) / (2 * gamma);
        double t    = ((zeta >= 0) ? 1 : -1) / (fabs(zeta) + sqrt(1 + zeta * zeta));
        double cs   = 1 / sqrt(1 + t * t);
        double sn   = cs * t;
        for (int r = 0; r < m; r++) {
          double x = ai[r], y = aj[r];
          ai[r] = cs * x - sn * y;
          aj[r] = sn * x + cs * y;
        }
        for (int r = 0; r < c; r++) {
          double x = vi[r], y = vj[r];
          vi[r] = cs * x - sn * y;
          vj[r] = sn * x + cs * y;
        }
      }
    }
  }

  // Singular values, and a stable sort into decreasing order
  // (selection: c is small).
  for (int j = 0; j < c; j++) s[j] = sqrt(lz_dot(A + (size_t) j * m, A + (size_t) j * m, m));

  double *ta = tmp, *tv = tmp + m;
  for (int a = 0; a < c; a++) {
    int best = a;
    for (int j = a + 1; j < c; j++) if (s[j] > s[best]) best = j;
    if (best == a) continue;
    // Rotate a..best right by one, so that ties keep their order.
    double sb = s[best];
    memcpy(ta, A + (size_t) best * m, (size_t) m * sizeof(double));
    memcpy(tv, V + (size_t) best * c, (size_t) c * sizeof(double));
    for (int j = best; j > a; j--) {
      s[j] = s[j - 1];
      memcpy(A + (size_t) j * m, A + (size_t) (j - 1) * m, (size_t) m * sizeof(double));
      memcpy(V + (size_t) j * c, V + (size_t) (j - 1) * c, (size_t) c * sizeof(double));
    }
    s[a] = sb;
    memcpy(A + (size_t) a * m, ta, (size_t) m * sizeof(double));
    memcpy(V + (size_t) a * c, tv, (size_t) c * sizeof(double));
  }

  return !converged;
}



//======================================================
// Group means (g x p) and counts, and the within-group
// residuals Xc = X - M[grp, ] (n x p), each scaled by
// 1 / sd (divisor n; 1 for a constant column), with the
// sds in `sd` (p). Sums run over samples in order.
//======================================================

static void lda_center(
    const double *X, int n, int p, const int *grp, int g,
    double *M, int *cnt, double *Xc, double *sd ) {

  memset(cnt, 0, (size_t) g * sizeof(int));
  for (int i = 0; i < n; i++) cnt[grp[i]]++;

  memset(M, 0, (size_t) g * p * sizeof(double));
  for (int j = 0; j < p; j++) {
    const double *x  = X  + (size_t) j * n;
    double       *mj = M  + (size_t) j * g;
    double       *xc = Xc + (size_t) j * n;
    for (int i = 0; i < n; i++) mj[grp[i]] += x[i];
    for (int k = 0; k < g; k++) mj[k] /= cnt[k];
    double ss = 0;
    for (int i = 0; i < n; i++) {
      xc[i] = x[i] - mj[grp[i]];
      ss   += xc[i] * xc[i];
    }
    double s = sqrt(ss / n);
    if (s == 0) s = 1;
    sd[j] = s;
    for (int i = 0; i < n; i++) xc[i] /= s;
  }
}



//======================================================
// The fit (steps 1 to 5 above). X is n x p, grp in
// 0..g-1 (every group non-empty), n > g >= 2, p >= 1.
// Allocates with R_alloc. Returns the rank of Z: if it
// is below p, the within-group covariance is singular
// and nothing else is set but `w_sv`.
//======================================================

typedef struct {
  double *means;    // g x p group means.
  double *xbar;     // p: the weighted centre.
  double *scaling;  // p x r coefficients.
  double *sv;       // r between-group singular values D_b.
  double  sv_sum2;  // sum of all D_b^2.
  int     r;        // Discriminant axes kept.
  double  w_min_sv; // The smallest singular value of Z.
} lda_fit_t;

static int lda_fit(
    const double *X, int n, int p, const int *grp, int g,
    const double *prior, int n_threads, lda_fit_t *fit ) {

  double *M   = (double *) R_alloc((size_t) g * p, sizeof(double));
  int    *cnt = (int *)    R_alloc(g, sizeof(int));
  double *Xc  = (double *) R_alloc((size_t) n * p, sizeof(double));
  double *sd  = (double *) R_alloc(p, sizeof(double));

  lda_center(X, n, p, grp, g, M, cnt, Xc, sd);
  fit->means = M;

  // Z = Xc / sqrt(n - g): its SVD from QR, then Jacobi on R.
  double f = 1 / sqrt((double) (n - g));
  for (size_t e = 0; e < (size_t) n * p; e++) Xc[e] *= f;

  double *R   = (double *) R_alloc((size_t) p * p, sizeof(double));
  double *V   = (double *) R_alloc((size_t) p * p, sizeof(double));
  double *D   = (double *) R_alloc(p, sizeof(double));
  int     mx  = (n > g) ? n : g;
  double *tmp = (double *) R_alloc((size_t) mx + p + 1, sizeof(double));

  lda_qr(Xc, n, p, R, tmp, n_threads);
  if (lda_jacobi_svd(R, p, p, V, D, tmp))
    warning("The SVD in lda() did not converge."); // # nocov

  int rank = 0;
  for (int j = 0; j < p; j++) if (D[j] > LDA_TOL) rank++;
  fit->w_min_sv = D[p - 1];
  if (rank < p) return rank;

  // W = S^-1 V D^-1 (p x p).
  double *W = V;
  for (int a = 0; a < p; a++)
    for (int j = 0; j < p; j++)
      W[j + (size_t) a * p] /= sd[j] * D[a];

  // xbar = prior' M.
  double *xbar = (double *) R_alloc(p, sizeof(double));
  for (int j = 0; j < p; j++) {
    double s = 0;
    for (int k = 0; k < g; k++) s += prior[k] * M[k + (size_t) j * g];
    xbar[j] = s;
  }
  fit->xbar = xbar;

  // B = diag(sqrt(n prior / (g - 1))) (M - xbar) W  (g x p).
  double *B  = (double *) R_alloc((size_t) g * p, sizeof(double));
  double *Mc = (double *) R_alloc((size_t) g * p, sizeof(double));
  for (int j = 0; j < p; j++)
    for (int k = 0; k < g; k++)
      Mc[k + (size_t) j * g] = sqrt(n * prior[k] / (g - 1)) * (M[k + (size_t) j * g] - xbar[j]);
  for (int a = 0; a < p; a++)
    for (int k = 0; k < g; k++) {
      double s = 0;
      for (int j = 0; j < p; j++) s += Mc[k + (size_t) j * g] * W[j + (size_t) a * p];
      B[k + (size_t) a * g] = s;
    }

  double *Vb  = (double *) R_alloc((size_t) p * p, sizeof(double));
  double *Db  = (double *) R_alloc(p, sizeof(double));
  if (lda_jacobi_svd(B, g, p, Vb, Db, tmp))
    warning("The SVD in lda() did not converge."); // # nocov

  // At most min(g, p) singular values are not 0 (g - 1, as the rows of
  // B are dependent, but rounding leaves the g-th nonzero).
  int nb = (g < p) ? g : p;
  double ss = 0;
  for (int a = 0; a < nb; a++) ss += Db[a] * Db[a];
  fit->sv_sum2 = ss;

  int r = 0;
  // (And none if the group means are the same, up to rounding.)
  for (int a = 0; a < nb; a++) if (Db[a] > LDA_TOL * Db[0] && Db[a] > 1e-8) r++;
  fit->r  = r;
  fit->sv = Db;

  // Scaling = W Vb[, 1:r].
  double *sc = (double *) R_alloc((size_t) p * (r > 0 ? r : 1), sizeof(double));
  for (int a = 0; a < r; a++)
    for (int j = 0; j < p; j++) {
      double s = 0;
      for (int l = 0; l < p; l++) s += W[j + (size_t) l * p] * Vb[l + (size_t) a * p];
      sc[j + (size_t) a * p] = s;
    }
  fit->scaling = sc;

  return rank;
}



//======================================================
// Scores (X - xbar) scaling, n x r, in parallel over
// samples.
//======================================================

static const double *sc_X, *sc_xbar, *sc_scaling;
static double       *sc_out;
static int           sc_n, sc_p, sc_r;

#define LDA_ROWS 64

static void *score_worker(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = sc_n, p = sc_p, r = sc_r;
  int n_tasks = (n + LDA_ROWS - 1) / LDA_ROWS;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int i1 = (task + 1) * LDA_ROWS;
    if (i1 > n) i1 = n;
    for (int i = task * LDA_ROWS; i < i1; i++)
      for (int a = 0; a < r; a++) {
        const double *w = sc_scaling + (size_t) a * p;
        double s = 0;
        for (int j = 0; j < p; j++) s += (sc_X[i + (size_t) j * n] - sc_xbar[j]) * w[j];
        sc_out[i + (size_t) a * n] = s;
      }
  }
  return NULL;
}



//======================================================
// .Call entry: lda() on the dense n x p `x`, with groups
// `grp` (1..g) and `prior` (g). Returns list(scores =
// n x r, scaling = p x r, means = g x p, sv = r between
// singular values, sv_sum2, rank = rank of Z, min_sv).
// With rank < p, only `rank` and `min_sv` are set.
//======================================================

SEXP C_lda(SEXP sexp_x, SEXP sexp_grp, SEXP sexp_g, SEXP sexp_prior, SEXP sexp_n_threads) {

  int n = nrows(sexp_x), p = ncols(sexp_x);
  int g = asInteger(sexp_g);
  int n_threads = asInteger(sexp_n_threads);

  if (!isReal(sexp_x) || !isInteger(sexp_grp) || LENGTH(sexp_grp) != n ||
      !isReal(sexp_prior) || LENGTH(sexp_prior) != g || g < 2 || n <= g ||
      p < 1 || n_threads < 1)
    error("Invalid arguments to C_lda()."); // # nocov

  const double *X = REAL(sexp_x);
  int *grp = (int *) R_alloc(n, sizeof(int));
  for (int i = 0; i < n; i++) {
    grp[i] = INTEGER(sexp_grp)[i] - 1;
    if (grp[i] < 0 || grp[i] >= g) error("Invalid arguments to C_lda()."); // # nocov
  }

  lda_fit_t fit;
  int rank = lda_fit(X, n, p, grp, g, REAL(sexp_prior), n_threads, &fit);

  const char *names[] = { "scores", "scaling", "means", "sv", "sv_sum2", "rank", "min_sv", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 5, ScalarInteger(rank));
  SET_VECTOR_ELT(res, 6, ScalarReal(fit.w_min_sv));

  if (rank == p) {
    int r = fit.r;
    SEXP scores  = allocMatrix(REALSXP, n, r); SET_VECTOR_ELT(res, 0, scores);
    SEXP scaling = allocMatrix(REALSXP, p, r); SET_VECTOR_ELT(res, 1, scaling);
    SEXP means   = allocMatrix(REALSXP, g, p); SET_VECTOR_ELT(res, 2, means);
    SEXP sv      = allocVector(REALSXP, r);    SET_VECTOR_ELT(res, 3, sv);
    SET_VECTOR_ELT(res, 4, ScalarReal(fit.sv_sum2));

    if (r > 0) memcpy(REAL(scaling), fit.scaling, (size_t) p * r * sizeof(double));
    memcpy(REAL(means), fit.means, (size_t) g * p * sizeof(double));
    if (r > 0) memcpy(REAL(sv), fit.sv, (size_t) r * sizeof(double));

    sc_X = X; sc_xbar = fit.xbar; sc_scaling = fit.scaling; sc_out = REAL(scores);
    sc_n = n; sc_p = p; sc_r = r;
    if (r > 0)
      run_parallel(score_worker, threads_for(n_threads, (double) n * p * r, 1e5),
                   n);
  }

  UNPROTECT(1);
  return res;
}



/*
 * Leave-one-out classification, for every number of leading columns
 * m = 1, ..., p at once (capda() chooses m by it).
 *
 * The fit without sample i classifies x_i to the group k that maximizes
 * log(pi_k) - d_k / 2, with d_k the Mahalanobis distance from x_i to the
 * group's mean without x_i, under the pooled within-group covariance
 * without x_i (divisor n - 1 - g). That is the full-rank form of the
 * classifier that lda()'s axes give; MASS's predict() and CV = TRUE use
 * it too (and the same priors in every fit).
 *
 * With E = Xc'Xc = R'R (R from the QR of the scaled residuals, as in the
 * fit), sphere every point by z = R^-T x. For sample i in group h, with
 * n_h > 1 and c = n_h / (n_h - 1), and u = z_i - mu_h:
 *
 *   - the other samples' scatter is E_-i = I - c u u' (sphered), with
 *     inverse I + c u u' / (1 - c u'u);
 *   - group h's mean moves to mu_h - u / (n_h - 1), so z_i - it = c u;
 *   - so, with e_k = z_i - mu_k for k != h, and df = n - 1 - g,
 *       d_k = df (|e_k|^2 + c (u'e_k)^2 / (1 - c u'u)),
 *       d_h = df c^2 u'u / (1 - c u'u).
 *
 * R^-T is lower triangular, so z's first m coordinates come from x's first
 * m columns alone, and the sums |e_k|^2, u'e_k and u'u for m columns are
 * running sums over coordinates: every m costs O(g) more per sample.
 *
 * A sample alone in its group is NA, as in MASS (leaving it out leaves its
 * group without a mean). A fit is singular (NA) from
 * the first m where R[m, m] is below LDA_TOL times sqrt(n - g) (the
 * column is collinear with the earlier ones within groups), or where
 * 1 - c u'u <= sqrt(eps) (x_i alone spans a direction), or df < 1.
 */

static const double *lo_X, *lo_R, *lo_sd, *lo_mu, *lo_lp;
static const int    *lo_grp, *lo_cnt;
static int           lo_n, lo_p, lo_g, lo_m_ok, lo_want_post;
static int          *lo_cls;
static double       *lo_post, *lo_work;

static void *loo_worker(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = lo_n, p = lo_p, g = lo_g;

  // Per thread: z (p), e2 (g), ue (g).
  double *z  = lo_work + (size_t) thread_i * (p + 2 * g);
  double *e2 = z + p;
  double *ue = e2 + g;

  for (int i = thread_i; i < n; i += n_threads) {

    int    h  = lo_grp[i];
    int    nh = lo_cnt[h];
    double c  = (nh > 1) ? (double) nh / (nh - 1) : 0;
    int    df = n - 1 - g;

    // z = R^-T (x_i / sd), by forward substitution.
    for (int a = 0; a < p; a++) {
      double s = lo_X[i + (size_t) a * n] / lo_sd[a];
      for (int b = 0; b < a; b++) s -= lo_R[b + (size_t) a * p] * z[b];
      z[a] = s / lo_R[a + (size_t) a * p];
    }

    for (int k = 0; k < g; k++) e2[k] = ue[k] = 0;
    double uu   = 0;
    int    dead = (df < 1 || nh == 1);

    for (int a = 0; a < p; a++) {

      int *cls = lo_cls + i + (size_t) a * n;
      if (dead || a >= lo_m_ok) { *cls = NA_INTEGER; dead = 1; continue; }

      double ua = z[a] - lo_mu[h + (size_t) a * g];
      uu += ua * ua;
      for (int k = 0; k < g; k++) {
        double ek = z[a] - lo_mu[k + (size_t) a * g];
        e2[k] += ek * ek;
        ue[k] += ua * ek;
      }

      double den = 1 - c * uu;
      if (den <= 1.4901161193847656e-08) { *cls = NA_INTEGER; dead = 1; continue; }

      // Scores, and the best (the first, on ties).
      int    best = -1;
      double top  = -INFINITY;
      for (int k = 0; k < g; k++) {
        double d;
        if (k == h) {
          d = c * c * uu / den;
        }
        else {
          d = e2[k] + c * ue[k] * ue[k] / den;
        }
        double sk = lo_lp[k] - df * d / 2;
        if (sk > top) { top = sk; best = k; }
      }
      *cls = best + 1;

      if (lo_want_post && a == p - 1) {
        double tot = 0;
        double *post = lo_post + i;
        for (int k = 0; k < g; k++) {
          double d  = (k == h) ? c * c * uu / den : e2[k] + c * ue[k] * ue[k] / den;
          double pk = exp(lo_lp[k] - df * d / 2 - top);
          post[(size_t) k * n] = pk;
          tot += pk;
        }
        for (int k = 0; k < g; k++) post[(size_t) k * n] /= tot;
      }
    }

    if (lo_want_post && (dead || p > lo_m_ok))
      for (int k = 0; k < g; k++) lo_post[i + (size_t) k * n] = NA_REAL;
  }
  return NULL;
}



//======================================================
// .Call entry: leave-one-out classes for the dense n x p
// `x`, groups `grp` (1..g) and `prior` (g). Returns
// list(class = n x p integer matrix, column m for the
// first m columns of x, NA where that fit is singular;
// posterior = n x g for all p columns).
//======================================================

SEXP C_lda_loo(SEXP sexp_x, SEXP sexp_grp, SEXP sexp_g, SEXP sexp_prior, SEXP sexp_n_threads) {

  int n = nrows(sexp_x), p = ncols(sexp_x);
  int g = asInteger(sexp_g);
  int n_threads = asInteger(sexp_n_threads);

  if (!isReal(sexp_x) || !isInteger(sexp_grp) || LENGTH(sexp_grp) != n ||
      !isReal(sexp_prior) || LENGTH(sexp_prior) != g || g < 2 || n <= g ||
      p < 1 || n_threads < 1)
    error("Invalid arguments to C_lda_loo()."); // # nocov

  const double *X = REAL(sexp_x);
  int *grp = (int *) R_alloc(n, sizeof(int));
  for (int i = 0; i < n; i++) {
    grp[i] = INTEGER(sexp_grp)[i] - 1;
    if (grp[i] < 0 || grp[i] >= g) error("Invalid arguments to C_lda_loo()."); // # nocov
  }

  double *M   = (double *) R_alloc((size_t) g * p, sizeof(double));
  int    *cnt = (int *)    R_alloc(g, sizeof(int));
  double *Xc  = (double *) R_alloc((size_t) n * p, sizeof(double));
  double *sd  = (double *) R_alloc(p, sizeof(double));
  double *R   = (double *) R_alloc((size_t) p * p, sizeof(double));
  double *v   = (double *) R_alloc(n, sizeof(double));

  lda_center(X, n, p, grp, g, M, cnt, Xc, sd);
  lda_qr(Xc, n, p, R, v, n_threads);

  // The first column at which the fit is singular.
  int m_ok = p;
  double lim = LDA_TOL * sqrt((double) (n - g));
  for (int a = 0; a < p; a++)
    if (!(fabs(R[a + (size_t) a * p]) > lim)) { m_ok = a; break; }

  // Sphered group means mu = R^-T (M / sd)' (g x p), up to m_ok.
  double *mu = (double *) R_alloc((size_t) g * p, sizeof(double));
  for (int k = 0; k < g; k++)
    for (int a = 0; a < m_ok; a++) {
      double s = M[k + (size_t) a * g] / sd[a];
      for (int b = 0; b < a; b++) s -= R[b + (size_t) a * p] * mu[k + (size_t) b * g];
      mu[k + (size_t) a * g] = s / R[a + (size_t) a * p];
    }

  double *lp = (double *) R_alloc(g, sizeof(double));
  for (int k = 0; k < g; k++) lp[k] = log(REAL(sexp_prior)[k]);

  const char *names[] = { "class", "posterior", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  SEXP cls  = allocMatrix(INTSXP,  n, p); SET_VECTOR_ELT(res, 0, cls);
  SEXP post = allocMatrix(REALSXP, n, g); SET_VECTOR_ELT(res, 1, post);

  int nt = threads_for(n_threads, (double) n * p * (p + g), 1e5);
  if (nt > n) nt = n;
  lo_work = (double *) R_alloc((size_t) nt * (p + 2 * g), sizeof(double));
  lo_X = X; lo_R = R; lo_sd = sd; lo_mu = mu; lo_lp = lp;
  lo_grp = grp; lo_cnt = cnt; lo_n = n; lo_p = p; lo_g = g; lo_m_ok = m_ok;
  lo_cls = INTEGER(cls); lo_post = REAL(post); lo_want_post = 1;

  // The forward substitution would divide by a zero R[a, a] past m_ok.
  if (m_ok < p) for (int a = m_ok; a < p; a++) if (R[a + (size_t) a * p] == 0) R[a + (size_t) a * p] = 1;

  run_parallel(loo_worker, nt, n);

  UNPROTECT(1);
  return res;
}
