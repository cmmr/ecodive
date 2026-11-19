// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Principal Coordinates Analysis (PCoA), a.k.a. classical
 * multidimensional scaling.
 *
 * With A = D^2 (elementwise) and J = I - 11'/n, the
 * coordinates come from the top eigenpairs of
 *
 *     B = -1/2 * J * A * J
 *
 * Only the k largest (algebraic) eigenpairs are wanted, with
 * k small. So instead of forming B (n x n) and decomposing it
 * in full like cmdscale(), they are found by block Lanczos
 * with thick restarts (lanczos.c), which only needs products
 * B * X. Those are computed straight from the packed lower
 * triangle of D that R's `dist` objects hold (gower.c, shared
 * with dbrda()).
 *
 * Every basis vector is kept orthogonal to the ones vector, so
 * J * X = X, and B * X = -1/2 * J * (A * X): one product with
 * A, then subtracting each column's mean.
 *
 * Determinism: results are bit-identical for any number of
 * threads (see gower.c). The starting vectors come from fixed
 * PCG streams.
 *
 * Sign convention: each axis is flipped so that its coordinate
 * with the largest absolute value is positive (the first such
 * sample, allowing for rounding; see lz_sign_index()).
 *
 *
 * Negative eigenvalue corrections (C_pcoa_correct(); Gower &
 * Legendre 1986; Legendre & Legendre 2012, sec. 9.3):
 *
 *   Lingoes (1971): d'^2 = d^2 + 2 c1 (i != j), so A' = A +
 *     2 c1 (11' - I) and B' = B + c1 J. Every eigenvector of B
 *     orthogonal to 1 keeps its eigenvalue plus c1; the ones
 *     vector keeps 0. With c1 = -lambda_min(B), every
 *     eigenvalue becomes >= 0, and the eigenvectors are B's.
 *     So it takes B's top k eigenpairs, as without a
 *     correction, plus B's smallest eigenvalue: the largest
 *     of -B, from the same solver.
 *
 *   Cailliez (1983): d' = d + c2 (i != j), so B' = B + 2 c2
 *     B2 + c2^2 / 2 J, with B2 = -1/2 J D J. c2 is the largest
 *     real eigenvalue of the 2n x 2n matrix [[0, 2B], [-I,
 *     -4 B2]]. For y orthogonal to 1, that matrix has
 *     eigenvalue c exactly when Q(c) y = 0, with
 *
 *         Q(c) = B + 2 c B2 + c^2 / 2 I,
 *
 *     i.e. when B' (for that c) is singular. Q(c) is positive
 *     definite for every c above c2, so c2 is the largest
 *     root of lambda_min(Q(c)). It is found without forming
 *     either matrix, by Rayleigh-Ritz on a growing subspace V
 *     (cailliez_constant()): the projected problem's largest
 *     real root, from the eigenvalues of its own small 2s x 2s
 *     companion matrix (lz_gen_eigenvalues()), is a lower
 *     bound on c2 that rises to it as V grows, and V grows by
 *     the residuals of Q(c)'s Ritz pairs. Products with B and
 *     B2 come from one pass over D (gower_mult2()). Then B''s
 *     top k eigenpairs come from Lanczos, reading d + c2 for
 *     each d.
 *
 *     The iteration can settle on a root that isn't the
 *     largest (with exactly symmetric designs, where V can be
 *     invariant under B and B2), so c2 is checked: B''s
 *     smallest eigenvalue, from a fresh random start, must
 *     not be below -sqrt(eps) * trace(B'). If it is, its
 *     eigenvector restarts the iteration, up to PC_CHECKS
 *     rounds.
 *
 *   Negative eigenvalues above -sqrt(eps) * trace(B) are
 *   taken as rounding error: then no constant is added, and
 *   the results are those of C_pcoa().
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON

// As LZ_TOL, LZ_TOL_STALL and LZ_STALL in lanczos.c, and its
// limit on products with blocks of vectors.
#define PC_TOL       1e-14
#define PC_TOL_STALL 1e-10
#define PC_STALL     10
#define PC_MAX_STEPS 10000

// The smallest eigenvalue of B is usually at the end of a
// dense cluster of small negative ones, where Lanczos
// converges far more slowly than for the largest. Its runs
// keep a basis of PC_BASIS_MIN vectors, and do the
// Rayleigh-Ritz step (O(basis^3)) only every PC_RR_EVERY
// steps and before each restart.
#define PC_BASIS_MIN 200
#define PC_RR_EVERY  4

// Ritz vectors from that run that start the Cailliez
// iteration, and that iteration's basis size and block
// size. Each of its steps solves a projected problem of
// the basis's size several times, so its basis is smaller.
#define PC_WARM           30
#define PC_BASIS_CAILLIEZ 120
#define PC_BLOCK          4

// The Cailliez run for B's smallest eigenvalue stops once
// that is certainly below -sqrt(eps) * trace and its
// residual is at most PC_TOL_WARM times the trace, so that
// its Ritz vectors make a good start. The check of B' stops
// as soon as a Ritz value shows an eigenvalue below the
// threshold, and otherwise when its residual is at most
// PC_TOL_CHECK times the trace. At most PC_CHECKS rounds.
#define PC_TOL_WARM  1e-6
#define PC_TOL_CHECK 1e-9
#define PC_CHECKS    10

#define CORRECT_LINGOES  1
#define CORRECT_CAILLIEZ 2


//------------------------------------------------------
// The top k eigenpairs of the Gower matrix that `g`
// reads, with unit eigenvectors to `vecs` (n x k).
//------------------------------------------------------
static void top_axes(
    gower_t *g, int n, int k, double trace,
    double *vecs, double *eig, double *res, int *steps, int *converged ) {

  lanczos_t lz = {
    .n = n, .k = k, .dims = n - 1, .deflate = 1, .scale = trace,
    .op = gower_mult, .ctx = g, .n_threads = g->n_req,
    .vecs = vecs, .vals = eig, .resid = res };

  lanczos(&lz);

  *steps     += lz.steps;
  *converged &= lz.converged;
}


//------------------------------------------------------
// Unit eigenvectors to coordinates: each column scaled by
// sqrt(eig), or 0 where eig <= 0, with its sign fixed.
//------------------------------------------------------
static void scale_axes(double *points, const double *eig, int n, int k) {

  for (int j = 0; j < k; j++) {

    double *x = points + (size_t) n * j;

    int imax = lz_sign_index(x, n);

    double mult = (eig[j] > 0) ? sqrt(eig[j]) : 0;
    if (x[imax] < 0) mult = -mult;

    for (int i = 0; i < n; i++) x[i] *= mult;
  }
}


//======================================================
// R interface. `sexp_dist` is a double vector holding a
// dist object's packed lower triangle for `sexp_n`
// samples (n >= 2). Returns list(points, eig, trace,
// steps, converged, resid), with points an n x k matrix of
// eigenvectors scaled by sqrt(eig), or 0 where eig <= 0.
//======================================================
SEXP C_pcoa(SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_n_threads) {

  int n         = asInteger(sexp_n);
  int k         = asInteger(sexp_k);
  int n_threads = asInteger(sexp_n_threads);
  int dims      = n - 1; // B's rank is at most n - 1.

  if (n < 2 || k < 1 || k > dims || n_threads < 1 ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_pcoa()."); // # nocov


  //--------------------------------------------------
  // trace(B) = (1/n) * sum_{i<j} d_ij^2, and products
  // with B for up to Lanczos's block size of vectors at
  // a time.
  //--------------------------------------------------
  gower_t g;
  double trace = gower_init(&g, REAL(sexp_dist), n, (k > 16) ? k : 16, n_threads);


  //--------------------------------------------------
  // The top k eigenpairs of B, as coordinates.
  //--------------------------------------------------
  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_res    = PROTECT(allocVector(REALSXP, k));
  double *points   = REAL(sexp_points);
  double *eig      = REAL(sexp_eig);

  int steps = 0, converged = 1;
  top_axes(&g, n, k, trace, points, eig, REAL(sexp_res), &steps, &converged);
  scale_axes(points, eig, n, k);

  const char *names[] = { "points", "eig", "trace", "steps", "converged", "resid", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 2, ScalarReal(trace));
  SET_VECTOR_ELT(sexp_result, 3, ScalarInteger(steps));
  SET_VECTOR_ELT(sexp_result, 4, ScalarLogical(converged));
  SET_VECTOR_ELT(sexp_result, 5, sexp_res);

  UNPROTECT(4);
  return sexp_result;
}



//------------------------------------------------------
// W = -B * V: a lanczos_op_t for B's smallest eigenvalues,
// as -B's largest. `ctx` is a gower_t.
//------------------------------------------------------
static void neg_gower_mult(const double *V, double *W, int m, void *ctx) {
  gower_mult(V, W, m, ctx);
  size_t len = (size_t) m * ((gower_t *) ctx)->n;
  for (size_t i = 0; i < len; i++) W[i] = -W[i];
}


//------------------------------------------------------
// The largest eigenvalue of -B (B as `g` reads it), from
// a run with a large basis (see PC_BASIS_MIN). Optionally
// stops early (see lanczos_t's early_stop), and writes up
// to n_ritz Ritz vectors to `ritz`. The eigenvector goes
// to `vec`.
//------------------------------------------------------
static double neg_largest(
    gower_t *g, int n, double trace, uint64_t stream, double tol,
    int early_stop, double stop_above, double stop_tol,
    double *vec, double *ritz, int n_ritz, int *n_ritz_out,
    int *steps, int *converged ) {

  double val = 0, res = 0;

  lanczos_t lz = {
    .n = n, .k = 1, .dims = n - 1, .deflate = 1, .scale = trace,
    .op = neg_gower_mult, .ctx = g, .n_threads = g->n_req,
    .vecs = vec, .vals = &val, .resid = &res,
    .max_basis = PC_BASIS_MIN, .rr_every = PC_RR_EVERY, .tol = tol,
    .early_stop = early_stop, .stop_above = stop_above, .stop_tol = stop_tol,
    .stream = stream, .ritz = ritz, .n_ritz = n_ritz };

  lanczos(&lz);

  if (n_ritz_out) *n_ritz_out = lz.n_ritz_out;
  *steps     += lz.steps;
  *converged &= lz.converged;

  return val;
}



//------------------------------------------------------
// Small helpers for the projected quadratic eigenvalue
// problem (H1 + 2 c H2 + c^2 / 2 I) z = 0, with H1 and H2
// s x s symmetric (leading dimension ld).
//------------------------------------------------------

// G = H1 + 2 c H2 + c^2 / 2 I (s x s, compact), then its
// eigen-decomposition (eigenvalues decreasing).
static void qep_matrix(
    const double *H1, const double *H2, int s, int ld, double c,
    double *G, double *gv, double *Y, double *tmp ) {

  for (int j = 0; j < s; j++)
    for (int i = 0; i < s; i++)
      G[(size_t) j * s + i] = H1[(size_t) j * ld + i] + 2 * c * H2[(size_t) j * ld + i] +
        ((i == j) ? c * c / 2 : 0);

  lz_sym_eigen(G, s, gv, Y, tmp);
}

// z' H z for a unit vector z.
static double quad_form(const double *Hm, int s, int ld, const double *z) {
  double total = 0;
  for (int j = 0; j < s; j++) {
    double hz = 0;
    for (int i = 0; i < s; i++) hz += Hm[(size_t) j * ld + i] * z[i];
    total += z[j] * hz;
  }
  return total;
}

// The larger root of c^2 / 2 + 2 b2 c + b1 = 0, i.e. of
// c^2 + 4 b2 c + 2 b1 = 0, without cancellation, or NaN if
// the roots aren't real.
static double larger_root(double b1, double b2) {
  double disc = 4 * b2 * b2 - 2 * b1;
  if (disc < 0) return NAN; // # nocov
  double sq = sqrt(disc);
  return (b2 > 0) ? -2 * b1 / (2 * b2 + sq) : -2 * b2 + sq;
}


//------------------------------------------------------
// The largest real root c of the projected problem: the
// largest real eigenvalue of [[0, 2 H1], [-I, -4 H2]].
// `tau`, a power of two near the distances' scale, keeps
// that matrix's blocks of similar size: its eigenvalues
// are computed for c / tau. A few Newton steps on the
// smallest eigenvalue of G(c) = H1 + 2 c H2 + c^2 / 2 I,
// which is 0 at the root and has derivative 2 z' H2 z + c,
// polish it. Returns 0 if there is no real root.
//------------------------------------------------------
static int qep_largest(
    const double *H1, const double *H2, int s, int ld, double tau,
    double *C, double *wr, double *wi, double *ort,
    double *G, double *gv, double *Y, double *tmp, double *c_out ) {

  int    s2   = 2 * s;
  double tau2 = tau * tau;

  memset(C, 0, (size_t) s2 * s2 * sizeof(double));
  for (int j = 0; j < s; j++) {
    for (int i = 0; i < s; i++) {
      C[(size_t) (s + j) * s2 + i]       =  2 * H1[(size_t) j * ld + i] / tau2;
      C[(size_t) (s + j) * s2 + (s + i)] = -4 * H2[(size_t) j * ld + i] / tau;
    }
    C[(size_t) j * s2 + (s + j)] = -1;
  }

  if (lz_gen_eigenvalues(C, s2, wr, wi, ort)) return 0; // # nocov

  // A repeated real eigenvalue can come back as a complex
  // pair with a tiny imaginary part: count those as real.
  double big = 1;
  for (int i = 0; i < s2; i++) {
    if (fabs(wr[i]) > big) big = fabs(wr[i]);
    if (fabs(wi[i]) > big) big = fabs(wi[i]);
  }

  int    found = 0;
  double best  = 0;
  for (int i = 0; i < s2; i++)
    if (fabs(wi[i]) <= sqrt(DBL_EPSILON) * big && (!found || wr[i] > best)) {
      best  = wr[i];
      found = 1;
    }
  if (!found) return 0; // # nocov

  double c = best * tau;

  for (int iter = 0; iter < 4; iter++) {

    qep_matrix(H1, H2, s, ld, c, G, gv, Y, tmp);

    double slope = 2 * quad_form(H2, s, ld, Y + (size_t) (s - 1) * s) + c;
    if (!(slope > 0)) break; // # nocov

    // Only small corrections: a large one would mean the
    // root is ill-conditioned, so keep the QR value.
    double dc = -gv[s - 1] / slope;
    if (!(fabs(dc) <= 1e-6 * tau)) break; // # nocov

    c += dc;
    if (fabs(dc) <= 4 * DBL_EPSILON * fabs(c)) break;
  }

  *c_out = c;
  return 1;
}


//------------------------------------------------------
// Follow the largest real root from c0, a value with G(c0)
// not positive definite (so that the root is at or above
// c0), by c <- the larger root of z' G(c) z = 0, with z
// G(c)'s eigenvector for its smallest eigenvalue. That is
// a lower bound on the root (each z gives one, and the
// root is their maximum), at least c, and it converges
// quadratically. Leaves G's decomposition at the final c
// in gv and Y.
//------------------------------------------------------
static double qep_follow(
    const double *H1, const double *H2, int s, int ld, double c0,
    double *G, double *gv, double *Y, double *tmp ) {

  double c = c0;

  for (int iter = 0; iter < 50; iter++) {

    qep_matrix(H1, H2, s, ld, c, G, gv, Y, tmp);

    const double *z  = Y + (size_t) (s - 1) * s;
    double        cn = larger_root(quad_form(H1, s, ld, z), quad_form(H2, s, ld, z));

    if (!(cn > c)) return c;
    double dc = cn - c;
    c = cn;
    if (dc <= 4 * DBL_EPSILON * fabs(c)) return c;
  }

  qep_matrix(H1, H2, s, ld, c, G, gv, Y, tmp); // # nocov
  return c;                                    // # nocov
}


//------------------------------------------------------
// Convergence bookkeeping, as in lanczos(): true once
// `worst` (a relative residual norm) is at most PC_TOL, or
// at most PC_TOL_STALL without halving for PC_STALL steps.
//------------------------------------------------------
typedef struct {
  double best;
  int    stalled;
} stall_t;

static int check_converged(stall_t *st, double worst) {
  if (worst < st->best / 2) { st->best = worst; st->stalled = 0; }
  else                      { st->stalled++;                    }
  return worst <= PC_TOL || (worst <= PC_TOL_STALL && st->stalled >= PC_STALL);
}


//------------------------------------------------------
// The Cailliez constant: the largest c with Q(c) = B + 2 c
// B2 + c^2 / 2 I singular on the space orthogonal to 1 (B
// and B2 as `g` reads distances, with no shift), by
// Rayleigh-Ritz on a subspace V that grows by the
// residuals Q(c) y - theta y of the Ritz pairs with the
// smallest theta.
//
// `start` holds n_start vectors orthogonal to 1, the first
// with Q(c) not positive definite on it for some c > 0
// (e.g. a negative Rayleigh quotient of B), which start V.
// The projected problem's root is followed from step to
// step by qep_follow(), and found afresh by qep_largest()
// at the first step, at each restart, and before accepting
// convergence, so that it is the largest root within V.
// Convergence shows that c is a root, but not that it is
// the largest: C_pcoa_correct() checks that afterwards.
// Sets *conv to whether it converged.
//------------------------------------------------------
static double cailliez_constant(
    gower_t *g, int n, double trace, const double *start, int n_start,
    int *steps_out, int *conv_out ) {

  int dims = n - 1;
  int p    = (PC_BLOCK < dims) ? PC_BLOCK : dims;
  int M    = (PC_BASIS_CAILLIEZ < dims) ? PC_BASIS_CAILLIEZ : dims;
  int q    = (M - p + 1) / 2;
  if (n_start > M - p) n_start = M - p;
  if (n_start < 1)     n_start = 1;

  int e;
  frexp(g->max_d, &e);
  double tau = ldexp(1.0, e);

  double *V   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *W1  = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *W2  = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *T   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *R   = (double *) R_alloc((size_t) n * p, sizeof(double));
  double *u   = (double *) R_alloc(n, sizeof(double));
  double *bu  = (double *) R_alloc(n, sizeof(double));
  double *u2  = (double *) R_alloc(n, sizeof(double));
  double *H1  = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *H2  = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *G   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Y   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Yk  = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *C   = (double *) R_alloc((size_t) 4 * M * M, sizeof(double));
  double *wr  = (double *) R_alloc((size_t) 2 * M, sizeof(double));
  double *wi  = (double *) R_alloc((size_t) 2 * M, sizeof(double));
  double *ort = (double *) R_alloc((size_t) 2 * M, sizeof(double));
  double *gv  = (double *) R_alloc(M, sizeof(double));
  double *tmp = (double *) R_alloc((size_t) 2 * M, sizeof(double));

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 20260116u, 57u);


  //--------------------------------------------------
  // Starting basis.
  //--------------------------------------------------
  memcpy(T, start, (size_t) n * n_start * sizeof(double));
  int s_v = lz_orthonormalize(V, 0, T, n_start, n, 1, NULL, 0, &rng);
  int s_w = 0;
  memcpy(V, T, (size_t) n * s_v * sizeof(double));

  int     steps = 0, converged = 0, fresh = 1, jumps = 0;
  double  c     = 0;
  stall_t st    = { INFINITY, 0 };

  while (1) {

    //------------------------------------------------
    // Products of the newest vectors with B and B2.
    //------------------------------------------------
    int m = s_v - s_w;
    gower_mult2(V + (size_t) n * s_w, W1 + (size_t) n * s_w, W2 + (size_t) n * s_w, m, g);
    steps++;

    lz_extend_h(V, W1, H1, M, s_w, s_v, n);
    lz_extend_h(V, W2, H2, M, s_w, s_v, n);
    s_w = s_v;


    //------------------------------------------------
    // The projected problem's root c, and Q_V(c)'s
    // eigenvectors; the last (eigenvalue 0) gives the
    // Ritz vector y = V z.
    //------------------------------------------------
    if (fresh) {
      if (!qep_largest(H1, H2, s_w, M, tau, C, wr, wi, ort, G, gv, Y, tmp, &c))
        break; // # nocov
      fresh = 0;
    }
    c = qep_follow(H1, H2, s_w, M, c, G, gv, Y, tmp);

    const double *z = Y + (size_t) (s_w - 1) * s_w;
    lz_combine(V,  n, s_w, z, 1, u);
    lz_combine(W1, n, s_w, z, 1, bu);
    lz_combine(W2, n, s_w, z, 1, u2);
    for (int i = 0; i < n; i++) bu[i] += 2 * c * u2[i] + c * c / 2 * u[i];

    double scale = trace + c * c * (n - 1) / 2;
    for (int j = 0; j < s_w; j++)
      if (fabs(gv[j]) > scale) scale = fabs(gv[j]);

    converged = check_converged(&st, sqrt(lz_dot(bu, bu, n)) / scale);

    // Before accepting, make sure c is the largest root.
    if (converged && jumps < 10) {
      double ch;
      if (qep_largest(H1, H2, s_w, M, tau, C, wr, wi, ort, G, gv, Y, tmp, &ch) &&
          ch > c + 1e-10 * tau) { // # nocov start
        c          = ch;
        converged  = 0;
        st.best    = INFINITY;
        st.stalled = 0;
        jumps++;
        continue; // Recompute the Ritz pair from the new c.
      } // # nocov end
    }

    if (converged || s_w == dims || steps >= PC_MAX_STEPS) break;

    if (user_interrupted()) error("Interrupted."); // # nocov


    //------------------------------------------------
    // Next block: the residuals Q(c) y - theta y of the
    // p Ritz pairs with the smallest theta (the first is
    // bu, above), orthonormalized against the basis. As c
    // changes between steps, the basis isn't a Krylov
    // space of any one Q(c), so these, rather than
    // products with the newest block, are what it lacks.
    //------------------------------------------------
    int nr = (p < s_w) ? p : s_w;
    memcpy(R, bu, (size_t) n * sizeof(double));
    for (int b = 1; b < nr; b++) {
      const double *zb = Y + (size_t) (s_w - 1 - b) * s_w;
      double       *rb = R + (size_t) n * b;
      lz_combine(V,  n, s_w, zb, 1, u);
      lz_combine(W1, n, s_w, zb, 1, rb);
      lz_combine(W2, n, s_w, zb, 1, u2);
      double th = gv[s_w - 1 - b];
      for (int i = 0; i < n; i++)
        rb[i] += 2 * c * u2[i] + (c * c / 2 - th) * u[i];
    }
    int r = lz_orthonormalize(V, s_v, R, nr, n, 1, NULL, 0, &rng);
    if (r == 0) break; // # nocov (only if the space is exhausted)


    //------------------------------------------------
    // Thick restart: keep the q Ritz vectors with the
    // smallest eigenvalues of Q_V(c), smallest first.
    //------------------------------------------------
    if (s_v + r > M) {

      int keep = (q + r <= M) ? q : M - r;

      for (int j = 0; j < keep; j++)
        memcpy(Yk + (size_t) j * s_w, Y + (size_t) (s_w - 1 - j) * s_w,
               (size_t) s_w * sizeof(double));

      lz_combine(V, n, s_w, Yk, keep, T);
      memcpy(V, T, (size_t) n * keep * sizeof(double));
      lz_combine(W1, n, s_w, Yk, keep, T);
      memcpy(W1, T, (size_t) n * keep * sizeof(double));
      lz_combine(W2, n, s_w, Yk, keep, T);
      memcpy(W2, T, (size_t) n * keep * sizeof(double));

      lz_restart_h(V, W1, H1, M, keep, n);
      lz_restart_h(V, W2, H2, M, keep, n);

      s_v   = s_w = keep;
      fresh = 1;
    }

    memcpy(V + (size_t) n * s_v, R, (size_t) n * r * sizeof(double));
    s_v += r;
  }

  *steps_out += steps;
  *conv_out   = converged;
  return c;
}


//======================================================
// R interface. As C_pcoa(), with a negative eigenvalue
// correction: `sexp_correction` 1 (Lingoes) or 2
// (Cailliez). Returns list(points, eig, trace, steps,
// converged, resid, constant, checked): eig and trace are
// those of the corrected matrix; constant is c1 (Lingoes)
// or c2 (Cailliez), or 0 if B had no eigenvalue below
// -sqrt(eps) * trace(B), in which case the rest is exactly
// as from C_pcoa(). `checked` is FALSE if c2 failed its
// check after PC_CHECKS rounds.
//======================================================
SEXP C_pcoa_correct(
    SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_correction,
    SEXP sexp_n_threads ) {

  int n          = asInteger(sexp_n);
  int k          = asInteger(sexp_k);
  int correction = asInteger(sexp_correction);
  int n_threads  = asInteger(sexp_n_threads);
  int dims       = n - 1;

  if (n < 2 || k < 1 || k > dims || n_threads < 1 ||
      (correction != CORRECT_LINGOES && correction != CORRECT_CAILLIEZ) ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_pcoa_correct()."); // # nocov

  gower_t g;
  double trace = gower_init(&g, REAL(sexp_dist), n, (k > 16) ? k : 16, n_threads);

  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_res    = PROTECT(allocVector(REALSXP, k));
  double *points   = REAL(sexp_points);
  double *eig      = REAL(sexp_eig);
  double *res      = REAL(sexp_res);

  int    steps = 0, converged = 1, checked = 1, fix = 0, n_warm = 0;
  double constant = 0;
  double *warm    = NULL;


  //--------------------------------------------------
  // B's smallest eigenvalue, the largest of -B.
  // Eigenvalues above -sqrt(eps) * trace count as zero,
  // and need no correction. Cailliez only needs to know
  // whether it is below that, which is certain as soon
  // as a Ritz value of -B is above sqrt(eps) * trace
  // (each is at most -B's largest eigenvalue); Lingoes
  // needs its value.
  //--------------------------------------------------
  double lmin = 0;

  if (trace > 0) {

    double *vmin = (double *) R_alloc(n, sizeof(double));
    int     cai  = (correction == CORRECT_CAILLIEZ);
    if (cai) warm = (double *) R_alloc((size_t) n * (PC_WARM + 1), sizeof(double));

    lmin = -neg_largest(
      &g, n, trace, 55u, 0, cai, sqrt(DBL_EPSILON) * trace, PC_TOL_WARM,
      vmin, cai ? warm + n : NULL, cai ? PC_WARM : 0, &n_warm, &steps, &converged );

    fix = lmin < -sqrt(DBL_EPSILON) * trace;
  }


  if (fix && correction == CORRECT_CAILLIEZ) {

    //------------------------------------------------
    // c2, checked: B' (distances d + c2) must have no
    // eigenvalue below -sqrt(eps) * trace(B'). If it has,
    // c2 was a root but not the largest, and B''s
    // eigenvector w for that eigenvalue has w' Q(c2) w < 0,
    // so a subspace that holds w has a larger root. Each
    // round starts afresh from w and the warm vectors; the
    // first from the warm vectors alone.
    //------------------------------------------------
    double *w      = warm;
    double  tracec = trace;

    for (int round = 0; ; round++) {

      int conv = 0;
      gower_set_shift(&g, 0);
      constant = cailliez_constant(
        &g, n, trace, (round == 0) ? warm + n : warm,
        (round == 0) ? n_warm : n_warm + 1, &steps, &conv );

      tracec = gower_set_shift(&g, constant);

      double limit = sqrt(DBL_EPSILON) * tracec;
      int    conv_check = 1;
      double lc = neg_largest(
        &g, n, tracec, 58u + (uint64_t) round, PC_TOL_CHECK, 1, limit, INFINITY,
        w, NULL, 0, NULL, &steps, &conv_check );

      if (lc <= limit) { converged &= conv && conv_check; break; }

      if (round + 1 >= PC_CHECKS) { converged = 0; checked = 0; break; } // # nocov
    }

    trace = tracec;
    top_axes(&g, n, k, trace, points, eig, res, &steps, &converged);
  }

  else {

    //------------------------------------------------
    // As C_pcoa(); for Lingoes, every eigenvalue gains
    // c1, and the trace (n - 1) c1.
    //------------------------------------------------
    top_axes(&g, n, k, trace, points, eig, res, &steps, &converged);

    if (fix) {
      constant = -lmin;
      for (int j = 0; j < k; j++) eig[j] += constant;
      trace += dims * constant;
    }
  }

  scale_axes(points, eig, n, k);

  const char *names[] = {
    "points", "eig", "trace", "steps", "converged", "resid", "constant", "checked", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 2, ScalarReal(trace));
  SET_VECTOR_ELT(sexp_result, 3, ScalarInteger(steps));
  SET_VECTOR_ELT(sexp_result, 4, ScalarLogical(converged));
  SET_VECTOR_ELT(sexp_result, 5, sexp_res);
  SET_VECTOR_ELT(sexp_result, 6, ScalarReal(constant));
  SET_VECTOR_ELT(sexp_result, 7, ScalarLogical(checked));

  UNPROTECT(4);
  return sexp_result;
}
