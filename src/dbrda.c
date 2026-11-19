// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Distance-based redundancy analysis (db-RDA), as in
 * vegan::dbrda() (McArdle & Anderson 2001).
 *
 * With A = D^2 (elementwise) and J = I - 11'/n, Gower's
 * centred matrix is
 *
 *     G = -1/2 * J * A * J
 *
 * The design's columns (centred) are orthonormalized by
 * Gram-Schmidt into Q = [Qz, Qx]: Qz spans the conditioning
 * variables, and Qx the constraints after removing Qz. Then
 *
 *   total inertia        trace(G)
 *   conditional inertia  trace(Qz' G Qz)
 *   constrained axes     eigenpairs (lambda, y) of Qx' G Qx,
 *                        with u = Qx y
 *   unconstrained axes   eigenpairs of (I-H) G (I-H), with
 *                        H = 11'/n + Q Q'
 *
 * vegan builds G and H G H in full (n x n) and decomposes
 * them. Here the constrained part only needs G * Q (a few
 * products with G, computed straight from the packed lower
 * triangle of D) and a p x p eigenproblem. Unconstrained axes
 * are only computed when there are fewer positive
 * constrained axes than requested, by block Lanczos with
 * thick restarts on (I-H) G (I-H), keeping every basis
 * vector orthogonal to Q and the ones vector.
 *
 * Site scores for constrained axes are vegan's "wa" scores,
 * (I - Qz Qz') G u / lambda; for unconstrained axes they are
 * the unit eigenvectors (vegan's CA$u). The R side applies
 * vegan's scaling constant.
 *
 * The Lanczos iteration and the small dense helpers are
 * shared with pcoa() and pca() (lanczos.c), and the products
 * with G with pcoa() (gower.c).
 *
 * Determinism: results are bit-identical for any number of
 * threads. G * X is computed (gower.c) in fixed blocks of
 * rows, each summed by one thread in a fixed order. All else
 * runs on the main thread. Lanczos starting vectors come from
 * a fixed PCG stream.
 *
 * Sign convention: each axis is flipped so that its score
 * with the largest absolute value is positive (the first such
 * sample, if tied).
 *
 *
 * capscale() is the same analysis of G's positive part, G+:
 * the principal coordinates with positive eigenvalues, as
 * vegan::capscale() analyzes them. When G has no eigenvalue
 * below zero (beyond vegan's tolerance; Euclidean
 * distances), G+ = G and nothing changes. A Lanczos run on
 * -G, as products straight from D, settles that: it stops as
 * soon as a Ritz value of -G proves a negative eigenvalue
 * (a Ritz value never exceeds the largest eigenvalue), and
 * otherwise runs until -G's largest eigenvalue has
 * converged.
 *
 * Otherwise G+ needs all of G's eigenvalues: G (scaled) is
 * formed, reduced to tridiagonal form T by Householder
 * reflections, and T's eigenvalues are found by the QL
 * algorithm. Eigenvectors are only computed for one side:
 * for the positive eigenvalues, so that G+ = Y Y', or for
 * the others, so that G+ = G - Y S Y' (S their signs),
 * whichever has fewer. They come from inverse iteration on
 * T, then the reflections. That needs about 4n^3/3 + 2n^2 r
 * operations, for r eigenvectors, against several times as
 * many for all of them.
 *
 * rda() (C_rda(), at the end) is the same analysis of
 * B = A A' for the centred (scaled) count table A, whose
 * products pca.c computes from the table's non-zeros: the
 * constrained and unconstrained solvers below never see
 * where B comes from (bop_t). Its species scores are A' u
 * for each axis's unit vector u.
 */

#include "ecodive.h"

#include <float.h>  // DBL_EPSILON
#include <stdlib.h> // qsort

// A design column is aliased (redundant) when Gram-Schmidt
// leaves less than this fraction of its centred norm: the
// same rule and tolerance as R's qr(), which vegan uses.
#define DBRDA_QR_TOL 1e-7

// capscale()'s probe for a negative eigenvalue: Lanczos
// with a large basis, and Rayleigh-Ritz every few steps
// (cheaper per step than pcoa()'s 200 and 4, for about as
// many steps).
#define CAP_BASIS    120
#define CAP_RR_EVERY 6

// Inverse iteration: at most this many solves for an
// eigenvector, and this many more once it has converged.
#define CAP_INV_ITERS 10
#define CAP_INV_EXTRA 2

// Eigenvalues of T closer than this (times T's norm) are a
// cluster, whose eigenvectors are kept orthogonal to each
// other (LAPACK's dstein uses the same).
#define CAP_ORTOL 1e-3

// Columns of a dense product per block.
#define CAP_BLOCK 16

// Eigenvectors per back-transformation task.
#define CAP_VEC_TASK 8


//------------------------------------------------------
// Small dense helpers, shared with lanczos.c. Vectors
// have length n; matrices are column-major.
//------------------------------------------------------
#define dot       lz_dot
#define center    lz_center
#define sym_eigen lz_sym_eigen
#define combine   lz_combine

#define project_out lz_project_out



//------------------------------------------------------
// Flip x (and y, if given) so that x's element with the
// largest absolute value is positive.
//------------------------------------------------------
static void fix_sign(double *x, int n, double *y, int ny) {
  int imax = lz_sign_index(x, n);
  if (x[imax] < 0) {
    for (int i = 0; i < n;  i++) x[i] = -x[i];
    for (int i = 0; i < ny; i++) y[i] = -y[i];
  }
}



//======================================================
// Products with B, the matrix analyzed: G itself, from
// D (gower.c), or capscale()'s G+ (see pos_part()):
//
//   r == 0           B = G
//   Y, rest == 0     B = Y Y'
//   Y, rest == 1     B = G - Y S Y', S = diag(sign)
//   g == NULL        B = A A', rda()'s centred table
//
// Each product is centred and has F projected out, as
// gower_mult() does for G (which reads F from g).
//======================================================
typedef struct {
  int           n;
  gower_t      *g;     // G from D, or NULL for rda()'s
  lanczos_op_t  tmult; // ... B = A A' from a table: tmult
  void         *tctx;  // ... returns B V, centred.
  const double *Y;     // n x r, or NULL.
  const double *sign;  // r signs (rest == 1).
  int           r;
  int           rest;
  double       *t;     // r x CAP_BLOCK scratch.
  double       *w;     // n x CAP_BLOCK scratch.
  const double *F;     // n x f orthonormal columns to
  int           f;     // project out, or NULL / 0.
} bop_t;


// Point both B and G's products at F.
static void bop_set_basis(bop_t *b, const double *F, int f) {
  b->F = F;
  b->f = f;
  if (b->g) {
    b->g->F = F;
    b->g->f = f;
  }
}


static void bop_mult(const double *V, double *W, int m, void *ctx) {

  bop_t *b = (bop_t *) ctx;
  int    n = b->n, r = b->r;

  if (b->g == NULL) {
    b->tmult(V, W, m, b->tctx);
    for (int v = 0; v < m; v++) project_out(b->F, b->f, W + (size_t) v * n, n);
    return;
  }

  if (b->Y == NULL || b->rest) gower_mult(V, W, m, b->g);
  if (b->Y == NULL) return;

  for (int v0 = 0; v0 < m; v0 += CAP_BLOCK) {

    int mc = (m - v0 < CAP_BLOCK) ? m - v0 : CAP_BLOCK;

    // t = Y' V, then w = Y t (times the signs).
    for (int v = 0; v < mc; v++) {
      const double *x = V + (size_t) (v0 + v) * n;
      double       *t = b->t + (size_t) v * r;
      for (int c = 0; c < r; c++) {
        t[c] = dot(b->Y + (size_t) c * n, x, n);
        if (b->rest) t[c] *= b->sign[c];
      }
    }

    for (int v = 0; v < mc; v++) {
      double *w = (b->rest) ? b->w : W + (size_t) (v0 + v) * n;
      combine(b->Y, n, r, b->t + (size_t) v * r, 1, w);
      center(w, n);
      project_out(b->F, b->f, w, n);
      if (b->rest) {
        double *out = W + (size_t) (v0 + v) * n;
        for (int i = 0; i < n; i++) out[i] -= w[i];
      }
    }
  }
}



//======================================================
// capscale(): is there an eigenvalue of G below -limit?
// The largest eigenvalue of -G, by Lanczos on products
// straight from D. A Ritz value above `limit` proves one,
// and stops the run; otherwise it runs until -G's
// largest eigenvalue has converged, which can take many
// steps when G's smallest eigenvalues are close together
// (full-rank Euclidean distances). After CAP_PROBE_MIN
// steps, or n / CAP_PROBE_DIV if more, by when the
// eigen-decomposition would cost about as much, it gives
// up and answers yes: pos_part() is right either way.
//======================================================
#define CAP_PROBE_MIN 30
#define CAP_PROBE_DIV 10

static void neg_mult(const double *V, double *W, int m, void *ctx) {
  gower_mult(V, W, m, ctx);
  size_t len = (size_t) m * ((gower_t *) ctx)->n;
  for (size_t i = 0; i < len; i++) W[i] = -W[i];
}

static int has_negative(gower_t *g, int n, double trace, double limit, int *steps) {

  double  val = 0, res = 0;
  double *vec = (double *) R_alloc(n, sizeof(double));
  int     max = (n / CAP_PROBE_DIV > CAP_PROBE_MIN) ? n / CAP_PROBE_DIV : CAP_PROBE_MIN;

  lanczos_t lz = {
    .n = n, .k = 1, .dims = n - 1, .deflate = 1, .scale = trace,
    .op = neg_mult, .ctx = g, .vecs = vec, .vals = &val, .resid = &res,
    .max_basis = CAP_BASIS, .rr_every = CAP_RR_EVERY,
    .early_stop = 1, .stop_above = limit, .stop_tol = INFINITY,
    .stream = 61u, .max_steps = max };

  lanczos(&lz);

  *steps += lz.steps;

  return val > limit || !lz.converged;
}



//======================================================
// Dense symmetric eigenvalues, for capscale().
//======================================================

// Householder reduction of the symmetric n x n matrix in
// A's lower triangle (column-major) to tridiagonal T =
// H' A H, H = H_0 H_1 ... H_{n-3}, H_k = I - tau_k v v'.
// T's diagonal goes to d, its sub-diagonal to e (n - 1),
// and each v (v[0] = 1, at row k + 1) to column k of A,
// below the diagonal.
//
// Step k needs p = B v for the trailing block B, then
// updates B -= v w' + w v'. Both are passes over B, so
// the update is put off: step k + 1 applies it to each
// column of its (smaller) trailing block just before
// reading that column for its own product, so that each
// step makes one pass. Only the column that step k + 1
// takes its reflection from is updated at once.
//
// A pass is split into tasks of TD_BLOCK columns, run in
// parallel. A column's own sum goes to sdot; its
// contributions to rows below go to its task's buffer,
// and the buffers are added up in task order. The tasks
// don't depend on the number of threads, so neither do
// the results.
#define TD_BLOCK 32

static double       *td_B;      // the trailing block, leading dimension td_ld
static int           td_ld, td_m;
static const double *td_pv;     // pending update: B -= pv pw' + pw pv',
static const double *td_pw;     // ... or NULL
static const double *td_v;      // product vector, or NULL
static double       *td_part;   // tasks x m
static double       *td_sdot;   // m

static void *td_pass(void *arg) {

  int share = ((worker_t *) arg)->i, shares = ((worker_t *) arg)->n;
  int m = td_m, tasks = (m + TD_BLOCK - 1) / TD_BLOCK;

  const double *pv = td_pv, *pw = td_pw, *v = td_v;

  for (int t = share; t < tasks; t += shares) {

    int     c0   = t * TD_BLOCK;
    int     c1   = (c0 + TD_BLOCK < m) ? c0 + TD_BLOCK : m;
    double *part = td_part + (size_t) t * m;

    if (v) for (int i = c0; i < m; i++) part[i] = 0;

    for (int j = c0; j < c1; j++) {

      double *col = td_B + (size_t) j * td_ld;

      if (pv && !v) {
        double pvj = pv[j], pwj = pw[j];
        for (int i = j; i < m; i++) col[i] -= pv[i] * pwj + pw[i] * pvj;
        continue;
      }

      if (pv) col[j] -= pv[j] * pw[j] + pw[j] * pv[j];

      double vj = v[j], s0 = col[j] * vj, s1 = 0;
      int    i  = j + 1;

      if (pv) {
        double pvj = pv[j], pwj = pw[j];
        for (; i + 1 < m; i += 2) {
          double a0 = col[i]     - (pv[i]     * pwj + pw[i]     * pvj);
          double a1 = col[i + 1] - (pv[i + 1] * pwj + pw[i + 1] * pvj);
          col[i] = a0; col[i + 1] = a1;
          part[i] += a0 * vj; part[i + 1] += a1 * vj;
          s0 += a0 * v[i];    s1 += a1 * v[i + 1];
        }
        if (i < m) {
          double a0 = col[i] - (pv[i] * pwj + pw[i] * pvj);
          col[i] = a0;
          part[i] += a0 * vj;
          s0 += a0 * v[i];
        }
      }
      else {
        for (; i + 1 < m; i += 2) {
          double a0 = col[i], a1 = col[i + 1];
          part[i] += a0 * vj; part[i + 1] += a1 * vj;
          s0 += a0 * v[i];    s1 += a1 * v[i + 1];
        }
        if (i < m) {
          part[i] += col[i] * vj;
          s0 += col[i] * v[i];
        }
      }

      td_sdot[j] = s0 + s1;
    }
  }

  return NULL;
}

// One pass over the m x m block B (leading dimension ld):
// the pending update (pv, pw), if any, then p = B v, if v.
static void td_run(
    double *B, int ld, int m, const double *pv, const double *pw,
    const double *v, double *p, int n_threads ) {

  td_B = B; td_ld = ld; td_m = m; td_pv = pv; td_pw = pw; td_v = v;

  int tasks = (m + TD_BLOCK - 1) / TD_BLOCK;
  int nt    = threads_for(n_threads, (double) m * m, 2.5e5);
  run_parallel(td_pass, (nt < tasks) ? nt : tasks, m);

  if (!v) return;
  for (int i = 0; i < m; i++) {
    double s = td_sdot[i];
    for (int t = 0; t < tasks && t * TD_BLOCK < i; t++)
      s += td_part[(size_t) t * m + i];
    p[i] = s;
  }
}

static void tridiagonalize(
    double *A, int n, double *d, double *e, double *tau, int n_threads ) {

  int tasks = (n + TD_BLOCK - 1) / TD_BLOCK;
  td_part   = (double *) R_alloc((size_t) tasks * n + 1, sizeof(double));
  td_sdot   = (double *) R_alloc(n, sizeof(double));

  double *p  = (double *) R_alloc(n, sizeof(double));
  double *w0 = (double *) R_alloc(n, sizeof(double));
  double *w1 = (double *) R_alloc(n, sizeof(double));

  const double *pv = NULL, *pw = NULL; // pending update of the trailing block

  for (int k = 0; k < n - 2; k++) {

    if (k % 64 == 0 && user_interrupted()) error("Interrupted."); // # nocov

    int     m = n - k - 1;
    double *v = A + (size_t) k * n + k + 1;
    double *B = A + (size_t) (k + 1) * n + k + 1; // trailing m x m block
    double *w = (k % 2) ? w1 : w0;

    // Column k is up to date: its reflection.
    d[k] = A[(size_t) k * n + k];

    double alpha = v[0], sig = 0;
    for (int i = 1; i < m; i++) sig += v[i] * v[i];

    tau[k] = 0;
    e[k]   = alpha;
    if (sig > 0) {
      double beta = -copysign(sqrt(alpha * alpha + sig), alpha);
      double sc   = 1 / (alpha - beta);
      for (int i = 1; i < m; i++) v[i] *= sc;
      e[k]   = beta;
      tau[k] = (beta - alpha) / beta;
    }
    v[0] = 1;

    // The pending update, and p = B v.
    td_run(B, n, m, pv, pw, (tau[k] != 0) ? v : NULL, p, n_threads);
    pv = pw = NULL;
    if (tau[k] == 0) continue;

    // w = tau p - (tau^2 / 2) (p'v) v.
    double h = 0;
    for (int i = 0; i < m; i++) { w[i] = tau[k] * p[i]; h += w[i] * v[i]; }
    h *= -0.5 * tau[k];
    for (int i = 0; i < m; i++) w[i] += h * v[i];

    // B's first column now; the rest at the next step.
    for (int i = 0; i < m; i++) B[i] -= v[i] * w[0] + w[i] * v[0];
    pv = v + 1;
    pw = w + 1;
  }

  if (pv) td_run(A + (size_t) (n - 1) * n + n - 1, n, 1, pv, pw, NULL, NULL, 1);

  if (n >= 2) {
    d[n - 2] = A[(size_t) (n - 2) * n + n - 2];
    e[n - 2] = A[(size_t) (n - 2) * n + n - 1];
  }
  d[n - 1] = A[(size_t) (n - 1) * n + n - 1];
}

// Eigenvalues of the symmetric tridiagonal matrix with
// diagonal d and sub-diagonal e (n - 1; both overwritten),
// by the QL algorithm with implicit shifts, as in the
// public-domain JAMA routine tql2() (lz_sym_eigen()),
// without the eigenvectors. Returns 1 if an eigenvalue
// failed to converge.
static int tri_eigenvalues(double *d, double *e, int n) {

  e[n - 1] = 0;

  double f = 0, tst1 = 0;
  int    failed = 0;

  for (int l = 0; l < n; l++) {

    double t = fabs(d[l]) + fabs(e[l]);
    if (t > tst1) tst1 = t;
    int m = l;
    while (m < n - 1 && fabs(e[m]) > DBL_EPSILON * tst1) m++;

    if (m > l) {
      int done = 0;
      for (int iter = 0; iter < 100; iter++) {

        double g = d[l];
        double p = (d[l + 1] - g) / (2 * e[l]);
        double r = hypot(p, 1);
        if (p < 0) r = -r;
        d[l]       = e[l] / (p + r);
        d[l + 1]   = e[l] * (p + r);
        double dl1 = d[l + 1];
        double h   = g - d[l];
        for (int i = l + 2; i < n; i++) d[i] -= h;
        f += h;

        p = d[m];
        double c = 1, c2 = 1, c3 = 1, sn = 0, s2 = 0;
        double el1 = e[l + 1];
        for (int i = m - 1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = sn;
          g  = c * e[i];
          h  = c * p;
          r  = hypot(p, e[i]);
          e[i + 1] = sn * r;
          sn = e[i] / r;
          c  = p / r;
          p  = c * d[i] - sn * g;
          d[i + 1] = h + sn * (c * g + sn * d[i]);
        }
        p    = -sn * s2 * c3 * el1 * e[l] / dl1;
        e[l] = sn * p;
        d[l] = c * p;

        if (fabs(e[l]) <= DBL_EPSILON * tst1) { done = 1; break; }
      }
      if (!done) failed = 1; // # nocov
    }
    d[l] = d[l] + f;
    e[l] = 0;
  }

  return failed;
}


static int cmp_double(const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}


// Eigenvectors of the symmetric tridiagonal T (diagonal d,
// sub-diagonal e) for its eigenvalues w[0..r-1], ascending,
// by inverse iteration, to the columns of Z (n x r).
// Eigenvalues within CAP_ORTOL * tnorm of the one before
// are a cluster: each eigenvector in it is kept orthogonal
// to the cluster's earlier ones, and its shift at least
// 10 eps * tnorm above the earlier one's, as in LAPACK's
// dstein. Starting vectors come from a fixed PCG stream.
static void tri_vectors(
    const double *d, const double *e, int n, const double *w, int r,
    double *Z, double tnorm ) {

  double *u1  = (double *) R_alloc(n, sizeof(double)); // U's diagonal,
  double *u2  = (double *) R_alloc(n, sizeof(double)); // ... first and
  double *u3  = (double *) R_alloc(n, sizeof(double)); // ... second super-diagonal.
  double *l   = (double *) R_alloc(n, sizeof(double)); // multipliers
  int    *piv = (int *)    R_alloc(n, sizeof(int));    // rows swapped?

  double ortol = CAP_ORTOL * tnorm;
  double pert  = 10 * DBL_EPSILON * tnorm;
  double tiny  = DBL_EPSILON * tnorm;
  double grow  = sqrt(0.1 / n) / (n * DBL_EPSILON * tnorm);
  if (!(tnorm > 0)) { pert = tiny = DBL_MIN; grow = 1 / DBL_EPSILON; } // # nocov

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 20261117u, 62u);

  int    c0   = 0;  // first eigenvector of the current cluster
  double prev = 0;  // the previous shift

  for (int j = 0; j < r; j++) {

    if (j % 64 == 0 && user_interrupted()) error("Interrupted."); // # nocov

    double lam = w[j];
    if (j == 0 || lam - w[j - 1] > ortol) c0 = j;
    else if (lam - prev < pert) lam = prev + pert;
    prev = lam;


    //------------------------------------------------
    // T - lam I = P L U, with partial pivoting. U has
    // two super-diagonals; tiny pivots are replaced by
    // +/- tiny, so that every solve is defined.
    //------------------------------------------------
    double a = d[0] - lam, b = (n > 1) ? e[0] : 0;
    for (int i = 0; i < n - 1; i++) {
      double c  = e[i];                               // row i + 1, column i
      double dn = d[i + 1] - lam;                     // ... column i + 1
      double en = (i + 1 < n - 1) ? e[i + 1] : 0;     // ... column i + 2
      if (fabs(c) > fabs(a)) {
        piv[i] = 1;
        l[i]   = a / c;
        u1[i]  = c;
        u2[i]  = dn;
        u3[i]  = en;
        a      = b - l[i] * dn;
        b      = -l[i] * en;
      }
      else {
        piv[i] = 0;
        if (fabs(a) < tiny) a = copysign(tiny, a);
        l[i]   = c / a;
        u1[i]  = a;
        u2[i]  = b;
        u3[i]  = 0;
        a      = dn - l[i] * b;
        b      = en;
      }
    }
    if (fabs(a) < tiny) a = copysign(tiny, a);
    u1[n - 1] = a;
    u2[n - 1] = u3[n - 1] = 0;
    if (n >= 2) u3[n - 2] = 0;


    //------------------------------------------------
    // Iterate from a random start.
    //------------------------------------------------
    double *x = Z + (size_t) j * n;
    lz_fill_random(x, n, &rng);
    double nrm = sqrt(dot(x, x, n));
    for (int i = 0; i < n; i++) x[i] /= nrm;

    int conv = 0;
    for (int it = 0; it < CAP_INV_ITERS; it++) {

      // Solve (T - lam I) y = x in place.
      for (int i = 0; i < n - 1; i++) {
        if (piv[i]) { double t = x[i]; x[i] = x[i + 1]; x[i + 1] = t; }
        x[i + 1] -= l[i] * x[i];
      }
      x[n - 1] /= u1[n - 1];
      if (n >= 2) x[n - 2] = (x[n - 2] - u2[n - 2] * x[n - 1]) / u1[n - 2];
      for (int i = n - 3; i >= 0; i--)
        x[i] = (x[i] - u2[i] * x[i + 1] - u3[i] * x[i + 2]) / u1[i];

      // Orthogonal to the cluster's earlier eigenvectors.
      double y_nrm = sqrt(dot(x, x, n));
      for (int c = c0; c < j; c++) {
        const double *z = Z + (size_t) c * n;
        double s = dot(z, x, n);
        for (int i = 0; i < n; i++) x[i] -= s * z[i];
      }

      nrm = sqrt(dot(x, x, n));
      if (!(nrm > 0)) { // # nocov start: a fresh start
        lz_fill_random(x, n, &rng);
        nrm = sqrt(dot(x, x, n));
      } // # nocov end
      for (int i = 0; i < n; i++) x[i] /= nrm;

      if (y_nrm >= grow && ++conv > CAP_INV_EXTRA) break;
    }
  }
}


// Back-transformation: Z <- H Z for the reflections that
// tridiagonalize() left in A and tau. Each task takes
// CAP_VEC_TASK columns of Z, applying the reflections
// from last to first.
static const double *bt_A, *bt_tau;
static double       *bt_Z;
static int           bt_n, bt_r;

static void *back_transform(void *arg) {

  int share = ((worker_t *) arg)->i, shares = ((worker_t *) arg)->n;
  int n = bt_n, tasks = (bt_r + CAP_VEC_TASK - 1) / CAP_VEC_TASK;

  for (int task = share; task < tasks; task += shares) {
    int j0 = task * CAP_VEC_TASK;
    int j1 = (j0 + CAP_VEC_TASK < bt_r) ? j0 + CAP_VEC_TASK : bt_r;
    for (int k = n - 3; k >= 0; k--) {
      double t = bt_tau[k];
      if (t == 0) continue;
      int           m = n - k - 1;
      const double *v = bt_A + (size_t) k * n + k + 1;
      int           j = j0;

      // Four columns at a time, each summed as on its own.
      for (; j + 3 < j1; j += 4) {
        double *z0 = bt_Z + (size_t) j * n + k + 1;
        double *z1 = z0 + n, *z2 = z1 + n, *z3 = z2 + n;
        double  s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (int i = 0; i < m; i++) {
          s0 += v[i] * z0[i]; s1 += v[i] * z1[i];
          s2 += v[i] * z2[i]; s3 += v[i] * z3[i];
        }
        s0 *= t; s1 *= t; s2 *= t; s3 *= t;
        for (int i = 0; i < m; i++) {
          z0[i] -= s0 * v[i]; z1[i] -= s1 * v[i];
          z2[i] -= s2 * v[i]; z3[i] -= s3 * v[i];
        }
      }
      for (; j < j1; j++) {
        double *z = bt_Z + (size_t) j * n + k + 1;
        double  s = 0;
        for (int i = 0; i < m; i++) s += v[i] * z[i];
        s *= t;
        for (int i = 0; i < m; i++) z[i] -= s * v[i];
      }
    }
  }

  return NULL;
}


//------------------------------------------------------
// capscale(): sets up `b` for products with G+, the
// positive part of G (with distances `dist`), and returns
// trace(G+) and the sum of G's negative eigenvalues in
// *real and *imag, on the scale of trace_g = trace(G).
// Eigenvalues within vegan's tolerance of 0,
// max(sqrt(eps), sqrt(eps) * lambda_1) on the scale that
// `adj` gives, count as 0.
//------------------------------------------------------
static void pos_part(
    const double *dist, int n, double trace_g, double adj, int n_threads,
    bop_t *b, double *real, double *imag ) {

  //--------------------------------------------------
  // A = G / trace(G), from distances scaled by their
  // largest value (so that no square overflows).
  //--------------------------------------------------
  size_t nd  = (size_t) n * (n - 1) / 2;
  double top = 0;
  for (size_t i = 0; i < nd; i++) if (dist[i] > top) top = dist[i];

  double *A  = (double *) R_alloc((size_t) n * n, sizeof(double));
  double *rs = (double *) R_alloc(n, sizeof(double));
  for (int i = 0; i < n; i++) rs[i] = 0;

  double tot = 0;
  for (int j = 0; j < n; j++) {
    const double *dj = dist + dist_col_start(j, n);
    for (int i = j + 1; i < n; i++) {
      double a = dj[i - j - 1] / top;
      a *= a;
      rs[i] += a;
      rs[j] += a;
    }
  }
  for (int i = 0; i < n; i++) { tot += rs[i]; rs[i] /= n; }
  double mean = tot / ((double) n * n);
  double tr   = tot / (2.0 * n); // trace of the scaled G

  for (int j = 0; j < n; j++) {
    const double *dj = dist + dist_col_start(j, n);
    double       *aj = A + (size_t) j * n;
    aj[j] = (rs[j] - 0.5 * mean) / tr;
    for (int i = j + 1; i < n; i++) {
      double a = dj[i - j - 1] / top;
      aj[i] = -0.5 * (a * a - rs[i] - rs[j] + mean) / tr;
    }
  }


  //--------------------------------------------------
  // T, and its eigenvalues (ascending).
  //--------------------------------------------------
  double *d   = (double *) R_alloc(n, sizeof(double));
  double *e   = (double *) R_alloc(n, sizeof(double));
  double *tau = (double *) R_alloc(n, sizeof(double));
  double *w   = (double *) R_alloc(n, sizeof(double));
  double *we  = (double *) R_alloc(n, sizeof(double));

  for (int i = 0; i < n; i++) tau[i] = e[i] = 0;
  tridiagonalize(A, n, d, e, tau, n_threads);

  double tnorm = 0;
  for (int i = 0; i < n; i++) {
    double t = fabs(d[i]) + ((i > 0) ? fabs(e[i - 1]) : 0) + ((i < n - 1) ? fabs(e[i]) : 0);
    if (t > tnorm) tnorm = t;
  }

  memcpy(w,  d, n * sizeof(double));
  memcpy(we, e, n * sizeof(double));
  if (tri_eigenvalues(w, we, n))
    error("The eigenvalues of the distances' Gower matrix did not converge."); // # nocov
  qsort(w, n, sizeof(double), cmp_double);


  //--------------------------------------------------
  // vegan's tolerance, on A's scale: G = trace_g * A,
  // and vegan's eigenvalues are G's times adj.
  //--------------------------------------------------
  double zero = sqrt(DBL_EPSILON);
  double tol  = zero * w[n - 1];
  if (tol < zero / (trace_g * adj)) tol = zero / (trace_g * adj);

  int    n_pos = 0;
  double pos = 0, neg = 0;
  for (int i = 0; i < n; i++) {
    if (w[i] >  tol) { pos += w[i]; n_pos++; }
    if (w[i] < -tol) { neg += w[i]; }
  }

  *real = pos * trace_g;
  *imag = neg * trace_g;


  //--------------------------------------------------
  // Eigenvectors for the side with fewer: the top n_pos
  // eigenvalues, or the other n - n_pos.
  //--------------------------------------------------
  int           rest = (n - n_pos < n_pos);
  int           r    = rest ? n - n_pos : n_pos;
  const double *wr   = rest ? w : w + (n - n_pos);

  double *Z    = (double *) R_alloc((size_t) n * r + 1, sizeof(double));
  double *sign = (double *) R_alloc(r + 1, sizeof(double));

  if (r > 0) {

    tri_vectors(d, e, n, wr, r, Z, tnorm);

    bt_A = A; bt_tau = tau; bt_Z = Z; bt_n = n; bt_r = r;
    // (run_parallel() threads only for 100 or more tasks; the
    // shares split the column groups among themselves.)
    int tasks = (r + CAP_VEC_TASK - 1) / CAP_VEC_TASK;
    int nt    = threads_for(n_threads, (double) n * n * r, 1e6);
    run_parallel(back_transform, (nt < tasks) ? nt : tasks, (r < 100) ? r : 100);
  }

  // Y = Z |w|^(1/2), on G's scale.
  for (int c = 0; c < r; c++) {
    double  s = sqrt(fabs(wr[c]) * trace_g);
    double *z = Z + (size_t) c * n;
    for (int i = 0; i < n; i++) z[i] *= s;
    sign[c] = (wr[c] < 0) ? -1 : 1;
  }

  b->Y    = Z;
  b->sign = sign;
  b->r    = r;
  b->rest = rest;
  b->t    = (double *) R_alloc((size_t) r * CAP_BLOCK + 1, sizeof(double));
  b->w    = (double *) R_alloc((size_t) n * CAP_BLOCK, sizeof(double));
}


//======================================================
// The design: its columns centred (Xc), and orthonormalized
// in order by Gram-Schmidt (two passes) into Q = [Qz, Qx],
// dropping the aliased ones. The first nz of the p columns
// are conditions, the rest constraints.
//======================================================
typedef struct {
  int     n, p, nz;
  int     rz, rx, rq;  // Kept condition / constraint columns.
  double *Xc, *Q;      // n x p, n x rq.
  int    *kept;        // p: 1 if not aliased.
  const double *w;     // cca(): sqrt(r), the square roots of the
                       // row weights (summing to 1), or NULL.
} cord_design_t;

// With w (cca()), columns are centred with the weights w^2
// and multiplied by w, so they are orthogonal to w, which
// takes the place of the ones vector.
static void cord_design(
    SEXP sexp_x, int n, int nz, const double *w, cord_design_t *d ) {

  int p = (n > 0) ? (int) (XLENGTH(sexp_x) / n) : 0;

  if (nz < 0 || nz > p || !isReal(sexp_x) || XLENGTH(sexp_x) != (R_xlen_t) n * p)
    error("Invalid design."); // # nocov

  double *Xc   = (double *) R_alloc((size_t) n * (p + 1), sizeof(double));
  double *Q    = (double *) R_alloc((size_t) n * (p + 1), sizeof(double));
  int    *kept = (int *)    R_alloc(p + 1, sizeof(int));
  int     rz = 0, rx = 0;

  // No columns (e.g. only Condition() terms, or ~ 1): an empty vector's data
  // pointer may be a sentinel, which even a 0-byte memcpy() may not read.
  if (p > 0) memcpy(Xc, REAL(sexp_x), (size_t) n * p * sizeof(double));

  for (int c = 0; c < p; c++) {

    double *xc = Xc + (size_t) c * n;
    double *q  = Q  + (size_t) (rz + rx) * n;

    for (int i = 0; i < n; i++)
      if (!R_FINITE(xc[i]))
        error("The design must not contain NA, NaN, or infinite values.");

    if (w) {
      double mu = 0;
      for (int i = 0; i < n; i++) mu += w[i] * w[i] * xc[i];
      for (int i = 0; i < n; i++) xc[i] = w[i] * (xc[i] - mu);
    }
    else {
      center(xc, n);
    }
    memcpy(q, xc, (size_t) n * sizeof(double));

    double norm0 = sqrt(dot(q, q, n));
    project_out(Q, rz + rx, q, n);
    project_out(Q, rz + rx, q, n);
    double norm = sqrt(dot(q, q, n));

    kept[c] = norm0 > 0 && norm >= DBRDA_QR_TOL * norm0;
    if (!kept[c]) continue;

    for (int i = 0; i < n; i++) q[i] /= norm;
    if (c < nz) rz++; else rx++;
  }

  d->n = n; d->p = p; d->nz = nz;
  d->rz = rz; d->rx = rx; d->rq = rz + rx;
  d->Xc = Xc; d->Q = Q; d->kept = kept; d->w = w;
}



//======================================================
// The constrained ordination of B (see bop_t), on the
// design d. trace_g is trace(B), real_g trace(B+) (the
// same, but for capscale()), imag_g the sum of B's
// negative eigenvalues (capscale()), all unscaled; adj
// scales every eigenvalue and inertia. `steps` and
// `dense` come from capscale()'s probe. Returns the list
// that C_dbrda() documents.
//======================================================
static SEXP cord_solve(
    bop_t *B, const cord_design_t *d, int k, double adj, double trace_g,
    double real_g, double imag_g, int steps, int dense ) {

  const double ZERO = sqrt(DBL_EPSILON); // vegan's threshold

  int     n  = d->n, p = d->p, nz = d->nz;
  int     rz = d->rz, rx = d->rx, rq = d->rq;
  int    *kept = d->kept;
  double *Xc   = d->Xc;
  const double *Q  = d->Q;
  const double *Qz = Q;
  const double *Qx = Q + (size_t) n * rz;
  double  total = trace_g * adj;
  int     converged = 1;


  //--------------------------------------------------
  // GQ = B * Q, then the conditional inertia and the
  // constrained matrix C = Qx' B Qx.
  //--------------------------------------------------
  double *GQ = (double *) R_alloc((size_t) n * (rq + 1), sizeof(double));

  bop_mult(Q, GQ, rq, B);
  const double *GQx = GQ + (size_t) n * rz;

  double conditional = 0;
  for (int j = 0; j < rz; j++)
    conditional += dot(Qz + (size_t) j * n, GQ + (size_t) j * n, n);
  conditional *= adj;
  if (conditional < ZERO) conditional = 0;

  double *C    = (double *) R_alloc((size_t) rx * rx + 1, sizeof(double));
  double *Y    = (double *) R_alloc((size_t) rx * rx + 1, sizeof(double));
  double *lam  = (double *) R_alloc(rx + 1, sizeof(double));
  double *tmp  = (double *) R_alloc(2 * (size_t) rx + 2, sizeof(double));

  for (int b = 0; b < rx; b++)
    for (int a = b; a < rx; a++) {
      double h = (a == b)
        ? dot(Qx + (size_t) n * a, GQx + (size_t) n * b, n)
        : (dot(Qx + (size_t) n * a, GQx + (size_t) n * b, n) +
           dot(Qx + (size_t) n * b, GQx + (size_t) n * a, n)) / 2;
      C[(size_t) b * rx + a] = h;
      C[(size_t) a * rx + b] = h;
    }

  if (rx > 0 && sym_eigen(C, rx, lam, Y, tmp))
    error("The constrained eigenproblem did not converge."); // # nocov


  //--------------------------------------------------
  // vegan's rule for zero eigenvalues: |lambda| below
  // max(ZERO, ZERO * largest). The rest, positive and
  // negative, make up the constrained inertia.
  //--------------------------------------------------
  int    n_nonzero = 0, n_pos = 0;
  double constrained = 0;
  double con_zero    = (rx > 0 && ZERO * lam[0] * adj > ZERO) ? ZERO * lam[0] * adj : ZERO;

  int *con_idx = (int *) R_alloc(rx + 1, sizeof(int)); // nonzero eigenvalues' indices
  for (int j = 0; j < rx; j++) {
    double l = lam[j] * adj;
    if (fabs(l) < con_zero) continue;
    con_idx[n_nonzero++] = j;
    constrained += l;
    if (l > 0) n_pos++;
  }

  int kc = (k < n_pos) ? k : n_pos;


  //--------------------------------------------------
  // Unconstrained axes, if more are wanted: the top
  // eigenpairs of (I-H) G (I-H), in the space orthogonal
  // to Q and the ones vector.
  //--------------------------------------------------
  int     kr = 0;
  int     dims_r = n - 1 - rq;
  double *r_vecs = NULL, *r_vals = NULL;

  if (kc < k && dims_r > 0 && real_g > 0) {

    kr     = (k - kc < dims_r) ? k - kc : dims_r;
    r_vecs = (double *) R_alloc((size_t) n * kr, sizeof(double));
    r_vals = (double *) R_alloc(kr, sizeof(double));

    bop_set_basis(B, Q, rq);

    // cca(): deflate w, as the first basis column, in
    // place of the ones vector.
    const double *basis = Q;
    if (d->w) {
      double *wq = (double *) R_alloc((size_t) n * (rq + 1), sizeof(double));
      memcpy(wq, d->w, (size_t) n * sizeof(double));
      if (rq > 0) memcpy(wq + n, Q, (size_t) n * rq * sizeof(double));
      basis = wq;
    }

    lanczos_t lz = {
      .n = n, .k = kr, .dims = dims_r, .deflate = (d->w == NULL),
      .basis = basis, .n_basis = rq + (d->w != NULL), .scale = real_g,
      .op = bop_mult, .ctx = B,
      .vecs = r_vecs, .vals = r_vals,
      .resid = (double *) R_alloc(kr, sizeof(double)) };

    lanczos(&lz);
    steps    += lz.steps;
    converged = converged && lz.converged;

    // Keep the positive, non-zero ones (vegan's ordResid).
    double res_zero = (ZERO * r_vals[0] * adj > ZERO) ? ZERO * r_vals[0] * adj : ZERO;
    int    keep     = 0;
    while (keep < kr && r_vals[keep] * adj >= res_zero) keep++;
    kr = keep;
  }


  //--------------------------------------------------
  // Output.
  //--------------------------------------------------
  int m = kc + kr;

  int n_bp = 0;
  for (int c = nz; c < p; c++) n_bp += kept[c];

  SEXP sexp_scores  = PROTECT(allocMatrix(REALSXP, n, m));
  SEXP sexp_eig     = PROTECT(allocVector(REALSXP, m));
  SEXP sexp_con_eig = PROTECT(allocVector(REALSXP, n_nonzero));
  SEXP sexp_kept    = PROTECT(allocVector(LGLSXP,  p));
  SEXP sexp_biplot  = PROTECT(allocMatrix(REALSXP, n_bp, kc));
  SEXP sexp_u       = PROTECT(allocMatrix(REALSXP, n, n_pos));

  double *scores = REAL(sexp_scores);
  double *eig    = REAL(sexp_eig);
  double *biplot = REAL(sexp_biplot);
  double *u      = REAL(sexp_u);

  for (int j = 0; j < n_nonzero; j++)
    REAL(sexp_con_eig)[j] = lam[con_idx[j]] * adj;

  for (int c = 0; c < p; c++)
    LOGICAL(sexp_kept)[c] = kept[c];

  // Constrained axes: u = Qx y; wa = (I - Qz Qz') G u / lambda.
  // u is returned for every positive axis (vegan's centroids
  // use them all), wa only for the first kc.
  for (int j = 0; j < n_pos; j++) {

    int           e  = con_idx[j];
    const double *y  = Y + (size_t) e * rx;
    double       *uj = u + (size_t) j * n;

    combine(Qx, n, rx, y, 1, uj);
    if (j >= kc) continue;

    double *wa = scores + (size_t) j * n;
    combine(GQx, n, rx, y, 1, wa);
    project_out(Qz, rz, wa, n);
    for (int i = 0; i < n; i++) wa[i] /= lam[e];
    if (d->w) for (int i = 0; i < n; i++) wa[i] /= d->w[i];

    fix_sign(wa, n, uj, n);
    eig[j] = lam[e] * adj;

    // Biplot scores: correlations of the centred constraint
    // columns with u.
    int row = 0;
    for (int c = nz; c < p; c++) {
      if (!kept[c]) continue;
      const double *xc = Xc + (size_t) c * n;
      biplot[(size_t) j * n_bp + row++] = dot(xc, uj, n) / sqrt(dot(xc, xc, n));
    }
  }

  // Unconstrained axes.
  for (int j = 0; j < kr; j++) {
    double *x = scores + (size_t) (kc + j) * n;
    memcpy(x, r_vecs + (size_t) j * n, (size_t) n * sizeof(double));
    if (d->w) for (int i = 0; i < n; i++) x[i] /= d->w[i];
    fix_sign(x, n, NULL, 0);
    eig[kc + j] = r_vals[j] * adj;
  }

  const char *names[] = {
    "scores", "eig", "n_con", "con_eig", "total", "conditional",
    "constrained", "kept", "biplot", "u", "steps", "converged",
    "real", "imaginary", "dense", "" };

  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result,  0, sexp_scores);
  SET_VECTOR_ELT(sexp_result,  1, sexp_eig);
  SET_VECTOR_ELT(sexp_result,  2, ScalarInteger(kc));
  SET_VECTOR_ELT(sexp_result,  3, sexp_con_eig);
  SET_VECTOR_ELT(sexp_result,  4, ScalarReal(total));
  SET_VECTOR_ELT(sexp_result,  5, ScalarReal(conditional));
  SET_VECTOR_ELT(sexp_result,  6, ScalarReal(constrained));
  SET_VECTOR_ELT(sexp_result,  7, sexp_kept);
  SET_VECTOR_ELT(sexp_result,  8, sexp_biplot);
  SET_VECTOR_ELT(sexp_result,  9, sexp_u);
  SET_VECTOR_ELT(sexp_result, 10, ScalarInteger(steps));
  SET_VECTOR_ELT(sexp_result, 11, ScalarLogical(converged));
  SET_VECTOR_ELT(sexp_result, 12, ScalarReal(real_g * adj));
  SET_VECTOR_ELT(sexp_result, 13, ScalarReal(imag_g * adj));
  SET_VECTOR_ELT(sexp_result, 14, ScalarLogical(dense));

  UNPROTECT(7);
  return sexp_result;
}



//======================================================
// R interface.
//
//   sexp_dist   Packed lower triangle of D (double),
//               for n >= 2 samples.
//   sexp_x      n x p design matrix (double): the first
//               `nz` columns are conditioning variables,
//               the rest constraints. Centred here.
//   sexp_k      Axes wanted (>= 1).
//   sexp_adj    Multiplier for G (vegan's 1/(n-1) for
//               "mean squared distances", else 1).
//   sexp_pos    1 for capscale(): analyze G+, not G.
//
// Returns list(
//   scores        n x m: wa scores for constrained axes,
//                 unit eigenvectors for unconstrained ones,
//   eig           m eigenvalues, times adj,
//   n_con         how many of the m axes are constrained,
//   con_eig       all non-zero constrained eigenvalues,
//   total, conditional, constrained   inertia,
//   kept          per design column, 1 if not aliased,
//   biplot        (kept constraint columns) x n_con,
//   u             n x (positive constrained axes) unit
//                 eigenvectors (lc scores, unscaled),
//   steps, converged,
//   real          trace(G+) times adj (= total for G),
//   imaginary     the sum of G's negative eigenvalues
//                 times adj, for capscale() (else 0),
//   dense         TRUE if G+ was found by a full eigen-
//                 decomposition (G had a negative one) )
//======================================================
SEXP C_dbrda(
    SEXP sexp_dist, SEXP sexp_n, SEXP sexp_x, SEXP sexp_nz,
    SEXP sexp_k, SEXP sexp_adj, SEXP sexp_pos, SEXP sexp_n_threads ) {

  int    n         = asInteger(sexp_n);
  int    nz        = asInteger(sexp_nz);
  int    k         = asInteger(sexp_k);
  double adj       = asReal(sexp_adj);
  int    positive  = asLogical(sexp_pos) == TRUE;
  int    n_threads = asInteger(sexp_n_threads);
  int    p         = (n > 0) ? (int) (XLENGTH(sexp_x) / n) : 0;

  if (n < 2 || k < 1 || n_threads < 1 || nz < 0 || nz > p ||
      !isReal(sexp_x) || XLENGTH(sexp_x) != (R_xlen_t) n * p ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2 || !(adj > 0))
    error("Invalid arguments to C_dbrda()."); // # nocov

  const double ZERO = sqrt(DBL_EPSILON); // vegan's threshold

  //--------------------------------------------------
  // trace(G) = (1/n) * sum_{i<j} d_ij^2, and products
  // with G (gower.c, shared with pcoa()).
  //--------------------------------------------------
  gower_t g;
  // Up to max(16, k, p) vectors per pass over D (k axes; p model columns),
  // as pcoa() does, rather than several passes for more than 16.
  int    max_m   = (k > p) ? k : p;
  if (max_m < 16) max_m = 16;
  if (max_m > n)  max_m = n;
  double trace_g = gower_init(&g, REAL(sexp_dist), n, max_m, n_threads); // unscaled


  cord_design_t des;
  cord_design(sexp_x, n, nz, NULL, &des);


  //--------------------------------------------------
  // B: G, or for capscale() G+. G+ = G unless G has an
  // eigenvalue below vegan's tolerance, -max(ZERO, ZERO
  // * lambda_1) after adj. lambda_1 is at least
  // trace(G) / (n - 1), so a probe with that bound in
  // its place misses none (and may find ones that
  // pos_part() then counts as 0).
  //--------------------------------------------------
  int     steps = 0, dense = 0;
  double  real_g = trace_g, imag_g = 0;
  bop_t   B = { .n = n, .g = &g };

  if (positive && trace_g > 0) {

    double limit = ZERO * trace_g * adj / (n - 1);
    if (limit < ZERO) limit = ZERO;

    if (has_negative(&g, n, trace_g, limit / adj, &steps)) {
      dense = 1;
      pos_part(REAL(sexp_dist), n, trace_g, adj, n_threads, &B, &real_g, &imag_g);
    }
  }


  return cord_solve(&B, &des, k, adj, trace_g, real_g, imag_g, steps, dense);
}




//======================================================
// rda(): the same analysis of B = A A', where A = J Z D^-1
// is the centred (and, with `scale`, standardized) count
// table after normalization (pca.c), with adj = 1/(n - 1)
// so that eigenvalues are variances, as in vegan::rda().
// Products with B take two passes over the table's
// non-zeros, never forming A or B.
//
//   sexp_counts, sexp_margin, sexp_norm, sexp_pseudocount,
//   sexp_scale   as for C_pca();
//   sexp_x, sexp_nz, sexp_k   as for C_dbrda().
//
// Returns C_dbrda()'s list, plus
//   loadings   p x m: A' u for each axis's unit
//              eigenvector u (lc scores, unscaled), in
//              the axes' (fixed) signs.
//======================================================
SEXP C_rda(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm,
    SEXP sexp_pseudocount, SEXP sexp_scale, SEXP sexp_x, SEXP sexp_nz,
    SEXP sexp_k, SEXP sexp_n_threads ) {

  int    norm        = asInteger(sexp_norm);
  double pseudocount = asReal(sexp_pseudocount);
  int    scale       = asLogical(sexp_scale) == TRUE;
  int    nz          = asInteger(sexp_nz);
  int    k           = asInteger(sexp_k);
  int    n_threads   = asInteger(sexp_n_threads);

  if (k < 1 || n_threads < 1)
    error("Invalid arguments to C_rda()."); // # nocov

  init_n_ptrs(10);

  pca_table_t tbl;
  double trace_g = pca_table_init(
    sexp_counts, sexp_margin, norm, pseudocount, scale, n_threads, &tbl );

  int n = tbl.n, p = tbl.p;

  cord_design_t des;
  cord_design(sexp_x, n, nz, NULL, &des);

  bop_t B = { .n = n, .g = NULL, .tmult = pca_gram_mult, .tctx = tbl.ctx };

  SEXP sexp_res = PROTECT(cord_solve(
    &B, &des, k, 1.0 / (n - 1), trace_g, trace_g, 0, 0, 0 ));


  //--------------------------------------------------
  // Loadings: A' u, for the constrained axes' u (the
  // first n_con columns of `u`) and the unconstrained
  // axes' eigenvectors (columns n_con.. of `scores`).
  //--------------------------------------------------
  int     kc     = asInteger(VECTOR_ELT(sexp_res, 2));
  SEXP    sexp_s = VECTOR_ELT(sexp_res, 0);
  int     m      = ncols(sexp_s);
  double *U      = (double *) R_alloc((size_t) n * m + 1, sizeof(double));

  if (kc > 0)
    memcpy(U, REAL(VECTOR_ELT(sexp_res, 9)), (size_t) n * kc * sizeof(double));
  if (m > kc)
    memcpy(U + (size_t) n * kc, REAL(sexp_s) + (size_t) n * kc, (size_t) n * (m - kc) * sizeof(double));

  SEXP sexp_load = PROTECT(allocMatrix(REALSXP, p, m));
  if (m > 0) pca_cross(U, m, REAL(sexp_load), n_threads);

  free_all();

  // The result, with `loadings` appended.
  int  len   = length(sexp_res);
  SEXP out   = PROTECT(allocVector(VECSXP, len + 1));
  SEXP names = PROTECT(allocVector(STRSXP, len + 1));
  SEXP old   = getAttrib(sexp_res, R_NamesSymbol);
  for (int i = 0; i < len; i++) {
    SET_VECTOR_ELT(out, i, VECTOR_ELT(sexp_res, i));
    SET_STRING_ELT(names, i, STRING_ELT(old, i));
  }
  SET_VECTOR_ELT(out, len, sexp_load);
  SET_STRING_ELT(names, len, mkChar("loadings"));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(4);
  return out;
}




//======================================================
// cca(): the same analysis of B = Qbar Qbar', with Qbar
// the chi-square residuals of the table (cca.c), the
// design centred and weighted by the row weights r, and
// sqrt(r) in place of the ones vector, as in ter Braak
// (1986) and vegan::cca(). Scores are divided by sqrt(r)
// (vegan's scaling 2, without its rda() constant).
//
//   sexp_counts, sexp_margin   the table, as for C_pca();
//   sexp_x, sexp_nz, sexp_k    as for C_dbrda();
//   sexp_perms   n x n_perm 0-based permutations, or NULL;
//   sexp_n_perm  permutations for the adjusted R-squared
//                (0 for none);
//   sexp_seed    their seed, when sexp_perms is NULL;
//   sexp_fixed   TRUE to permute Qbar's rows against the
//                fixed design even without conditions
//                (varpart(chisquare = TRUE), as vegan's).
//
// Returns C_dbrda()'s list, plus
//   loadings      p x m: vegan's species scores, scaling 2;
//   w             sqrt(r), n;
//   perm_inertia  the constrained inertia for each
//                 permutation (n_perm): with conditions
//                 (or sexp_fixed), that of conditions and constraints
//                 together, with Qbar's rows permuted
//                 against the fixed design;
//   perm_cond     ... and of the conditions alone (only
//                 with conditions; else 0s).
// or, if some samples have no observations, list(empty =
// their 1-based indices).
//======================================================
SEXP C_cca(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_x, SEXP sexp_nz,
    SEXP sexp_k, SEXP sexp_perms, SEXP sexp_n_perm, SEXP sexp_seed,
    SEXP sexp_fixed, SEXP sexp_n_threads ) {

  int fixed     = asLogical(sexp_fixed) == TRUE;
  int nz        = asInteger(sexp_nz);
  int k         = asInteger(sexp_k);
  int n_perm    = asInteger(sexp_n_perm);
  int n_threads = asInteger(sexp_n_threads);

  if (k < 1 || n_threads < 1 || n_perm < 0)
    error("Invalid arguments to C_cca()."); // # nocov

  init_n_ptrs(10);

  int    n, p, *empty;
  double trace_g = cca_table_init(sexp_counts, sexp_margin, n_threads, &n, &p, &empty);

  if (trace_g < 0) {
    int ne = 0;
    for (int i = 0; i < n; i++) ne += empty[i];
    SEXP idx = PROTECT(allocVector(INTSXP, ne));
    for (int i = 0, a = 0; i < n; i++) if (empty[i]) INTEGER(idx)[a++] = i + 1;
    const char *names[] = { "empty", "" };
    SEXP out = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(out, 0, idx);
    UNPROTECT(2);
    return out;
  }

  if (n < 2) error("Invalid arguments to C_cca()."); // # nocov

  const double *w = cca_sqrt_r();

  cord_design_t des;
  cord_design(sexp_x, n, nz, w, &des);

  bop_t B = { .n = n, .g = NULL, .tmult = cca_gram_mult, .tctx = NULL };

  SEXP sexp_res = PROTECT(cord_solve(&B, &des, k, 1.0, trace_g, trace_g, 0, 0, 0));


  //--------------------------------------------------
  // Loadings: Dc^-1/2 Qbar' u, for the constrained axes'
  // u and the unconstrained axes' unit eigenvectors
  // (their scores times w).
  //--------------------------------------------------
  int     kc     = asInteger(VECTOR_ELT(sexp_res, 2));
  SEXP    sexp_s = VECTOR_ELT(sexp_res, 0);
  int     m      = ncols(sexp_s);
  double *U      = (double *) R_alloc((size_t) n * m + 1, sizeof(double));

  if (kc > 0)
    memcpy(U, REAL(VECTOR_ELT(sexp_res, 9)), (size_t) n * kc * sizeof(double));
  for (int j = kc; j < m; j++)
    for (int i = 0; i < n; i++)
      U[(size_t) j * n + i] = REAL(sexp_s)[(size_t) j * n + i] * w[i];

  SEXP sexp_load = PROTECT(allocMatrix(REALSXP, p, m));
  if (m > 0) cca_cross(U, m, REAL(sexp_load));

  SEXP sexp_w = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(sexp_w), w, (size_t) n * sizeof(double));


  //--------------------------------------------------
  // Permuted constrained inertia, for the adjusted
  // R-squared.
  //--------------------------------------------------
  SEXP sexp_pi = PROTECT(allocVector(REALSXP, n_perm));
  SEXP sexp_pc = PROTECT(allocVector(REALSXP, n_perm));
  int  px      = des.p - nz;
  for (int i = 0; i < n_perm; i++) REAL(sexp_pi)[i] = REAL(sexp_pc)[i] = 0;
  if (n_perm > 0 && des.rx > 0) {
    perm_source_t ps;
    perm_source_init(&ps, n, sexp_perms, R_NilValue, asInteger(sexp_seed));
    if (des.rz > 0 || fixed)
      cca_perm_inertia_fixed(des.Q, des.rq, des.rz, &ps, n_perm, REAL(sexp_pi), REAL(sexp_pc));
    else
      cca_perm_inertia(REAL(sexp_x) + (size_t) n * nz, px, &ps, n_perm, REAL(sexp_pi));
  }


  // The result, with loadings, w and perm_inertia appended.
  int  len   = length(sexp_res);
  SEXP out   = PROTECT(allocVector(VECSXP, len + 4));
  SEXP names = PROTECT(allocVector(STRSXP, len + 4));
  SEXP old   = getAttrib(sexp_res, R_NamesSymbol);
  for (int i = 0; i < len; i++) {
    SET_VECTOR_ELT(out, i, VECTOR_ELT(sexp_res, i));
    SET_STRING_ELT(names, i, STRING_ELT(old, i));
  }
  SET_VECTOR_ELT(out, len,     sexp_load);
  SET_VECTOR_ELT(out, len + 1, sexp_w);
  SET_VECTOR_ELT(out, len + 2, sexp_pi);
  SET_STRING_ELT(names, len,     mkChar("loadings"));
  SET_STRING_ELT(names, len + 1, mkChar("w"));
  SET_STRING_ELT(names, len + 2, mkChar("perm_inertia"));
  SET_VECTOR_ELT(out, len + 3, sexp_pc);
  SET_STRING_ELT(names, len + 3, mkChar("perm_cond"));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(7);
  return out;
}




//======================================================
// varpart() (and forward_select()): the inertia that any
// set of design columns explains, from one pass over the
// data. With Xc the design's columns centred, returns
//
//   total   trace(B), unscaled (vegan's SS.Y), and
//   cross   Xc' B Xc (m x m),
//
// for B = G (Gower's matrix, from distances) or B = A A'
// (the centred, normalized count table, as in rda()). The
// inertia explained by any subset S of the columns is then
// trace(R^-T Xc_S' B Xc_S R^-1), for the QR factor R of
// Xc_S, without another pass over the data.
//
//   sexp_mode        0: sexp_y is a packed lower triangle
//                    of distances, for sexp_n samples;
//                    1: sexp_y is a count table, read with
//                    sexp_margin, sexp_norm, sexp_pseudocount
//                    as for C_pca() (no scaling).
//   sexp_x           n x m design (double), centred here.
//
// Deterministic for any number of threads: products with
// B are (gower.c, pca.c), and the cross products are
// summed here on one thread in a fixed order.
//======================================================
SEXP C_cord_cross(
    SEXP sexp_mode, SEXP sexp_y, SEXP sexp_n, SEXP sexp_margin,
    SEXP sexp_norm, SEXP sexp_pseudocount, SEXP sexp_x,
    SEXP sexp_n_threads ) {

  int mode      = asInteger(sexp_mode);
  int n_threads = asInteger(sexp_n_threads);

  if ((mode != 0 && mode != 1) || n_threads < 1 || !isReal(sexp_x) || !isMatrix(sexp_x))
    error("Invalid arguments to C_cord_cross()."); // # nocov

  int n = nrows(sexp_x), m = ncols(sexp_x);

  if (n < 2) error("Invalid arguments to C_cord_cross()."); // # nocov

  double *Xc    = (double *) R_alloc((size_t) n * m + 1, sizeof(double));
  SEXP    sexp_cross = PROTECT(allocMatrix(REALSXP, m, m));
  double *cross = REAL(sexp_cross);
  double  total;

  if (m > 0) memcpy(Xc, REAL(sexp_x), (size_t) n * m * sizeof(double));
  for (size_t i = 0; i < (size_t) n * m; i++)
    if (!R_FINITE(Xc[i]))
      error("The design must not contain NA, NaN, or infinite values.");
  for (int c = 0; c < m; c++) center(Xc + (size_t) c * n, n);


  if (mode == 0) {

    if (asInteger(sexp_n) != n || XLENGTH(sexp_y) != (R_xlen_t) n * (n - 1) / 2 || !isReal(sexp_y))
      error("Invalid arguments to C_cord_cross()."); // # nocov

    gower_t g;
    int     max_m = (m < 16) ? 16 : m;
    if (max_m > n) max_m = n;
    total = gower_init(&g, REAL(sexp_y), n, max_m, n_threads);

    double *GX = (double *) R_alloc((size_t) n * m + 1, sizeof(double));
    if (m > 0) gower_mult(Xc, GX, m, &g);

    for (int b = 0; b < m; b++)
      for (int a = b; a < m; a++) {
        double h = (a == b)
          ? dot(Xc + (size_t) n * a, GX + (size_t) n * b, n)
          : (dot(Xc + (size_t) n * a, GX + (size_t) n * b, n) +
             dot(Xc + (size_t) n * b, GX + (size_t) n * a, n)) / 2;
        cross[(size_t) b * m + a] = h;
        cross[(size_t) a * m + b] = h;
      }
  }

  else {

    init_n_ptrs(10);

    pca_table_t tbl;
    total = pca_table_init(
      sexp_y, sexp_margin, asInteger(sexp_norm), asReal(sexp_pseudocount),
      0, n_threads, &tbl );

    if (tbl.n != n) {
      free_all();
      error("Invalid arguments to C_cord_cross()."); // # nocov
    }

    // L = A' Xc (p x m), then cross = L' L.
    int     p = tbl.p;
    double *L = (double *) R_alloc((size_t) p * m + 1, sizeof(double));
    if (m > 0) pca_cross(Xc, m, L, n_threads);

    free_all();

    for (int b = 0; b < m; b++)
      for (int a = b; a < m; a++) {
        double h = dot(L + (size_t) p * a, L + (size_t) p * b, p);
        cross[(size_t) b * m + a] = h;
        cross[(size_t) a * m + b] = h;
      }
  }

  const char *names[] = { "total", "cross", "" };
  SEXP out = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, ScalarReal(total));
  SET_VECTOR_ELT(out, 1, sexp_cross);

  UNPROTECT(2);
  return out;
}
