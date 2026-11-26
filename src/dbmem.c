// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * dbmem(): distance-based Moran's eigenvector maps (Borcard & Legendre
 * 2002; Dray, Legendre & Peres-Neto 2006). Written from those papers;
 * the results were checked against adespatial's outputs only.
 *
 * 1. The threshold t defaults to the longest edge of the minimum
 *    spanning tree of the distances (Prim's algorithm on the packed
 *    triangle, O(n^2)).
 *
 * 2. The spatial weights are w_ij = 1 - (d_ij / 4t)^2 for d_ij <= t
 *    (i != j), and 0 otherwise. Up to a constant, these are the
 *    truncated distances of PCNM: with d* = d (d <= t) or 4t (d > t),
 *    -1/2 d*^2 + 8 t^2 = 8 t^2 w, so PCoA of d* has the same
 *    eigenvectors as the doubly centred weights A = C W C (C = I - 11'/n),
 *    with eigenvalues 8 t^2 (lambda + 1) for A's lambda.
 *
 * 3. The MEMs are A's eigenvectors, other than the ones vector (A 1 = 0).
 *    The ones vector is removed first: H, the Householder reflection that
 *    takes 1 / sqrt(n) to -e1, gives H A H a zero first row and column,
 *    and the eigenvectors of its trailing (n - 1) x (n - 1) block B,
 *    padded with a leading 0 and reflected back by H, are A's
 *    eigenvectors orthogonal to 1 (all n - 1 of them, however many of
 *    A's eigenvalues are 0).
 *
 * 4. B's eigenpairs come from a dense solver: Householder reduction to
 *    tridiagonal form, accumulation of the reflections, then the QL
 *    algorithm with implicit shifts (as the public-domain JAMA / EISPACK
 *    tred2() and tql2(), which lanczos.c's lz_sym_eigen() follows). The
 *    O(n^3) steps run in parallel: the matrix-vector products and rank-2
 *    updates of the reduction and the accumulation by column, and the
 *    plane rotations of each QL sweep (recorded, then applied) by row of
 *    the eigenvector matrix.
 *
 * 5. Moran's I of each MEM v, with the weights W as given, is
 *    n / S0 * v'Wv / v'v (S0 = sum of W), from the list of neighbour
 *    pairs. It equals n / S0 * lambda.
 *
 * Determinism: every element is computed by one thread, in the same
 * order as on one thread, so results are bit-identical for any number
 * of threads.
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON

// Neighbours are pairs within t * (1 + DBM_EPS), so that rounding in
// computed distances doesn't drop pairs at the threshold.
#define DBM_EPS 1e-6

// Work (in flops) that makes a thread worth waking.
#define DBM_MIN_WORK 5e5

#define AUTOCOR_POSITIVE 0
#define AUTOCOR_NEGATIVE 1
#define AUTOCOR_NONNULL  2
#define AUTOCOR_ALL      3


// Index of pair (i, j), i > j, in R's packed lower triangle.
static R_INLINE size_t pk(int n, int i, int j) {
  return (size_t) j * (2 * (size_t) n - j - 1) / 2 + (size_t) (i - j - 1);
}


//------------------------------------------------------
// The longest edge of the minimum spanning tree (Prim).
//------------------------------------------------------
static double mst_longest(const double *d, int n, double *key, int *done) {

  for (int i = 0; i < n; i++) { key[i] = R_PosInf; done[i] = 0; }
  key[0] = 0;
  double longest = 0;

  for (int it = 0; it < n; it++) {

    int u = -1;
    for (int i = 0; i < n; i++)
      if (!done[i] && (u < 0 || key[i] < key[u])) u = i;

    done[u] = 1;
    if (key[u] > longest) longest = key[u];

    for (int v = 0; v < n; v++) {
      if (done[v]) continue;
      double dv = (v > u) ? d[pk(n, v, u)] : d[pk(n, u, v)];
      if (dv < key[v]) key[v] = dv;
    }
  }

  return longest;
}



//======================================================
// Shared state of the parallel steps. Each step is one
// of the jobs below; workers take every n-th column (or
// row, or MEM), starting at their share's index.
//======================================================
enum {
  JOB_BUILD,    // B from the distances.
  JOB_MATVEC,   // p = beta * B22 v.
  JOB_RANK2,    // B22 -= v w' + w v'.
  JOB_ACCUM,    // Q22 = H_k Q22.
  JOB_ROTATE,   // A QL sweep's rotations on blocks of rows of V.
  JOB_MEMS      // Eigenvectors to MEMs, and Moran's I.
};

static struct {
  int job;

  // The input.
  const double *d;
  int    n;
  double tt, w4;     // Neighbour cutoff, and 1 / (4t)^2.
  const double *rm;  // Row means of W.
  double g;          // Grand mean of W.
  const double *wv;  // Deflation: w for B = A - u w' - w u'.
  double ui;         // u_i = 1 / sqrt(n), i >= 1.

  // The eigenproblem, m = n - 1.
  int    m;
  double *B, *Q, *V;
  int    k;          // Step: the column of the reflection.
  const double *v;   // Its vector (length m - k - 1).
  double beta;
  double *p;         // Matvec result / w (length m - k - 1).
  const double *rc, *rs;
  int    l, mm;      // Rotations mm - 1 down to l.

  // The MEMs.
  const int    *sel; // Columns of V, one per MEM.
  int    n_mem;
  double *out;       // n x n_mem.
  double *moran;
  double h, u0, s0;
  const int    *pi, *pj;
  const double *pw;
  size_t n_pairs;
} dbm;


static R_INLINE double weight(const double *d, int n, int i, int j, double tt, double w4) {
  if (i == j) return 0;
  double x = (i > j) ? d[pk(n, i, j)] : d[pk(n, j, i)];
  return (x <= tt) ? 1 - x * x * w4 : 0;
}


static void build_col(int j) {
  // B column j (original index J = j + 1), rows i (I = i + 1).
  int    n  = dbm.n, m = dbm.m, J = j + 1;
  double *b = dbm.B + (size_t) j * m;
  for (int i = 0; i < m; i++) {
    int I = i + 1;
    double a = weight(dbm.d, n, I, J, dbm.tt, dbm.w4) - dbm.rm[I] - dbm.rm[J] + dbm.g;
    b[i] = a - dbm.ui * dbm.wv[J] - dbm.wv[I] * dbm.ui;
  }
}


// A dot product in four interleaved partial sums (a fixed
// order, whatever the threads), which pipelines better
// than one running sum.
static R_INLINE double dot4(const double *x, const double *y, int L) {
  double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  int i = 0;
  for (; i + 3 < L; i += 4) {
    s0 += x[i]     * y[i];
    s1 += x[i + 1] * y[i + 1];
    s2 += x[i + 2] * y[i + 2];
    s3 += x[i + 3] * y[i + 3];
  }
  for (; i < L; i++) s0 += x[i] * y[i];
  return (s0 + s1) + (s2 + s3);
}


static void matvec_col(int j) {
  int L = dbm.m - dbm.k - 1, off = dbm.k + 1;
  const double *col = dbm.B + (size_t) (off + j) * dbm.m + off;
  dbm.p[j] = dbm.beta * dot4(col, dbm.v, L);
}


static void rank2_col(int j) {
  int L = dbm.m - dbm.k - 1, off = dbm.k + 1;
  double *col = dbm.B + (size_t) (off + j) * dbm.m + off;
  double vj = dbm.v[j], wj = dbm.p[j];
  for (int i = 0; i < L; i++) col[i] -= dbm.v[i] * wj + dbm.p[i] * vj;
}


static void accum_col(int j) {
  int L = dbm.m - dbm.k - 1, off = dbm.k + 1;
  double *col = dbm.Q + (size_t) (off + j) * dbm.m + off;
  double s = dbm.beta * dot4(dbm.v, col, L);
  if (s != 0) for (int i = 0; i < L; i++) col[i] -= s * dbm.v[i];
}


// A QL sweep's rotations mm - 1 down to l, on the rows of
// V in block b: rows DBM_ROWS * b onwards. Rotation i
// takes (V[r, i], V[r, i + 1]) to (c V[r, i] - s V[r, i +
// 1], s V[r, i] + c V[r, i + 1]); the new V[r, i] is the
// next rotation's V[r, i + 1], so it is carried. Each row
// is a chain of dependent steps, so a block's rows are
// interleaved. Blocks are fixed, whatever the threads.
#define DBM_ROWS 4

static void rotate_block(int b) {

  int m = dbm.m, l = dbm.l, mm = dbm.mm;
  int r0 = b * DBM_ROWS;
  const double *rc = dbm.rc, *rs = dbm.rs;

  if (r0 + DBM_ROWS <= m) {
    double *x0 = dbm.V + (size_t) r0 * m;
    double *x1 = x0 + m, *x2 = x1 + m, *x3 = x2 + m;
    double h0 = x0[mm], h1 = x1[mm], h2 = x2[mm], h3 = x3[mm];
    for (int i = mm - 1; i >= l; i--) {
      double c = rc[i], s = rs[i];
      double a0 = x0[i], a1 = x1[i], a2 = x2[i], a3 = x3[i];
      x0[i + 1] = s * a0 + c * h0;  h0 = c * a0 - s * h0;
      x1[i + 1] = s * a1 + c * h1;  h1 = c * a1 - s * h1;
      x2[i + 1] = s * a2 + c * h2;  h2 = c * a2 - s * h2;
      x3[i + 1] = s * a3 + c * h3;  h3 = c * a3 - s * h3;
    }
    x0[l] = h0; x1[l] = h1; x2[l] = h2; x3[l] = h3;
    return;
  }

  for (int r = r0; r < m; r++) {
    double *x = dbm.V + (size_t) r * m;
    double  h = x[mm];
    for (int i = mm - 1; i >= l; i--) {
      double c = rc[i], s = rs[i];
      double a = x[i];
      x[i + 1] = s * a + c * h;  h = c * a - s * h;
    }
    x[l] = h;
  }
}


static void mem_col(int c) {

  int n = dbm.n, m = dbm.m, col = dbm.sel[c];
  double *x = dbm.out + (size_t) c * n;

  // z = [0; y], then x = H z = z - h u (u'z).
  double s = 0;
  x[0] = 0;
  for (int r = 0; r < m; r++) {
    double y = dbm.V[(size_t) r * m + col];
    x[r + 1] = y;
    s += y;
  }
  s *= dbm.ui * dbm.h;
  x[0] -= dbm.u0 * s;
  for (int i = 1; i < n; i++) x[i] -= dbm.ui * s;

  // Unit length, then sqrt(n), with the sign fixed.
  double ss = 0;
  for (int i = 0; i < n; i++) ss += x[i] * x[i];
  int    imax = lz_sign_index(x, n);
  double mult = sqrt((double) n / ss);
  if (x[imax] < 0) mult = -mult;
  for (int i = 0; i < n; i++) x[i] *= mult;

  // Moran's I: n / S0 * x'Wx / x'x, with x'x = n and
  // x'Wx = 2 q over the pairs.
  double q = 0;
  for (size_t e = 0; e < dbm.n_pairs; e++)
    q += dbm.pw[e] * x[dbm.pi[e]] * x[dbm.pj[e]];
  dbm.moran[c] = 2 * q / dbm.s0;
}


static void *dbm_worker(void *arg) {

  int t  = ((worker_t *) arg)->i;
  int nt = ((worker_t *) arg)->n;
  int L  = dbm.m - dbm.k - 1;

  switch (dbm.job) {
    case JOB_BUILD:  for (int j = t; j < dbm.m;     j += nt) build_col(j);  break;
    case JOB_MATVEC: for (int j = t; j < L;         j += nt) matvec_col(j); break;
    case JOB_RANK2:  for (int j = t; j < L;         j += nt) rank2_col(j);  break;
    case JOB_ACCUM:  for (int j = t; j < L;         j += nt) accum_col(j);  break;
    case JOB_ROTATE:
      for (int b = t; b * DBM_ROWS < dbm.m; b += nt) rotate_block(b);
      break;
    case JOB_MEMS:   for (int c = t; c < dbm.n_mem; c += nt) mem_col(c);    break;
  }

  return NULL;
}


static void run_job(int job, int n_threads, double work) {
  dbm.job = job;
  int nt  = threads_for(n_threads, work, DBM_MIN_WORK);
  run_parallel(dbm_worker, nt, nt > 1 ? 100 : 1);
}



//------------------------------------------------------
// Eigenpairs of the symmetric m x m matrix dbm.B
// (column-major, full): eigenvalues to `vals` (in the
// QL algorithm's order, not sorted), and eigenvectors to
// dbm.V, row-major: V[r * m + c] is element r of
// eigenvector c. dbm.B is overwritten (V takes its
// storage); dbm.Q is workspace. Returns 1 if the QL
// iteration failed to converge.
//------------------------------------------------------
static int dbm_eigen(double *vals, double *e, double *betas, double *rc, double *rs, int n_threads) {

  int     m = dbm.m;
  double *B = dbm.B;
  double *d = vals;

  #define BB(i, j) B[(size_t) (j) * m + (i)]


  //--------------------------------------------------
  // Householder reduction to tridiagonal form. Step k
  // reflects rows / columns k + 1 .. m - 1, and keeps
  // its vector in column k below the subdiagonal.
  //--------------------------------------------------
  for (int k = 0; k < m - 2; k++) {

    int     L = m - k - 1;
    double *x = B + (size_t) k * m + k + 1;

    d[k] = BB(k, k);

    double tail = 0;
    for (int i = 1; i < L; i++) tail += x[i] * x[i];

    if (tail == 0) {
      betas[k] = 0;
      e[k]     = x[0];
      continue;
    }

    double norm  = sqrt(x[0] * x[0] + tail);
    double alpha = (x[0] > 0) ? -norm : norm;
    e[k]     = alpha;
    betas[k] = 1 / (norm * (norm + fabs(x[0])));
    x[0]    -= alpha;

    dbm.k = k; dbm.v = x; dbm.beta = betas[k]; dbm.p = rc; // rc as scratch (length m)
    run_job(JOB_MATVEC, n_threads, 2.0 * L * L);

    double K = 0;
    for (int i = 0; i < L; i++) K += rc[i] * x[i];
    K *= betas[k] / 2;
    for (int i = 0; i < L; i++) rc[i] -= K * x[i];

    run_job(JOB_RANK2, n_threads, 4.0 * L * L);

    if (k % 64 == 63) R_CheckUserInterrupt();
  }

  if (m >= 2) {
    d[m - 2] = BB(m - 2, m - 2);
    e[m - 2] = BB(m - 1, m - 2);
  }
  d[m - 1] = BB(m - 1, m - 1);
  e[m - 1] = 0;


  //--------------------------------------------------
  // Q = H_0 H_1 ... H_{m-3}, accumulated backwards.
  //--------------------------------------------------
  double *Q = dbm.Q;
  memset(Q, 0, sizeof(double) * (size_t) m * m);
  for (int i = 0; i < m; i++) Q[(size_t) i * m + i] = 1;

  for (int k = m - 3; k >= 0; k--) {
    if (betas[k] == 0) continue;
    int L = m - k - 1;
    dbm.k = k; dbm.v = B + (size_t) k * m + k + 1; dbm.beta = betas[k];
    run_job(JOB_ACCUM, n_threads, 4.0 * L * L);
    if (k % 64 == 0) R_CheckUserInterrupt();
  }

  // V (row-major) = Q, in B's storage.
  double *V = B;
  for (int c = 0; c < m; c++)
    for (int r = 0; r < m; r++)
      V[(size_t) r * m + c] = Q[(size_t) c * m + r];
  dbm.V = V;

  #undef BB


  //--------------------------------------------------
  // QL with implicit shifts (tql2), its rotations
  // applied to the rows of V once per sweep.
  //--------------------------------------------------
  double f = 0, tst1 = 0;
  int    failed = 0, sweeps = 0;

  for (int l = 0; l < m; l++) {

    double t = fabs(d[l]) + fabs(e[l]);
    if (t > tst1) tst1 = t;
    int mm = l;
    while (mm < m - 1 && fabs(e[mm]) > DBL_EPSILON * tst1) mm++;

    if (mm > l) {
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
        for (int i = l + 2; i < m; i++) d[i] -= h;
        f += h;

        p = d[mm];
        double c = 1, c2 = 1, c3 = 1, sn = 0, s2 = 0;
        double el1 = e[l + 1];
        for (int i = mm - 1; i >= l; i--) {
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
          rc[i] = c;
          rs[i] = sn;
        }

        dbm.rc = rc; dbm.rs = rs; dbm.l = l; dbm.mm = mm;
        run_job(JOB_ROTATE, n_threads, 6.0 * m * (mm - l));
        if (++sweeps % 64 == 0) R_CheckUserInterrupt();

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


// Decreasing eigenvalue, then increasing index.
static const double *sort_vals;
static int cmp_desc(const void *a, const void *b) {
  int i = *(const int *) a, j = *(const int *) b;
  if (sort_vals[i] > sort_vals[j]) return -1;
  if (sort_vals[i] < sort_vals[j]) return  1;
  return (i > j) - (i < j);
}



//======================================================
// R interface.
//
// d:         packed lower triangle of the distances.
// n:         samples (n >= 3).
// threshold: the truncation threshold, or NA for the
//            minimum spanning tree's longest edge.
// autocor:   0 positive, 1 negative, 2 non-null, 3 all.
// n_threads: threads.
//
// Returns list(vectors, values, moran, threshold, mst,
// n_pairs, islands, converged). With no pair of
// neighbours, only threshold, mst and n_pairs (0).
//======================================================
SEXP C_dbmem(SEXP sexp_d, SEXP sexp_n, SEXP sexp_threshold, SEXP sexp_autocor, SEXP sexp_n_threads) {

  int    n         = asInteger(sexp_n);
  double threshold = asReal(sexp_threshold);
  int    autocor   = asInteger(sexp_autocor);
  int    n_threads = asInteger(sexp_n_threads);

  if (!isReal(sexp_d) || n < 3 || n_threads < 1 || autocor < 0 || autocor > 3 ||
      XLENGTH(sexp_d) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_dbmem()."); // # nocov

  const double *d = REAL(sexp_d);
  int m = n - 1;


  //--------------------------------------------------
  // The threshold, and the pairs of neighbours.
  //--------------------------------------------------
  double *key  = (double *) R_alloc(n, sizeof(double));
  int    *flag = (int *)    R_alloc(n, sizeof(int));
  double  mst  = mst_longest(d, n, key, flag);
  if (ISNAN(threshold)) threshold = mst;

  double tt = threshold * (1 + DBM_EPS);
  double w4 = 1 / (16 * threshold * threshold);

  size_t n_dist = (size_t) n * (n - 1) / 2, n_pairs = 0;
  for (size_t e = 0; e < n_dist; e++) if (d[e] <= tt) n_pairs++;

  // All samples at one place (threshold 0) gives no
  // weights either.
  if (!(threshold > 0)) n_pairs = 0;

  if (n_pairs == 0) {
    const char *names[] = { "threshold", "mst", "n_pairs", "" };
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, ScalarReal(threshold));
    SET_VECTOR_ELT(res, 1, ScalarReal(mst));
    SET_VECTOR_ELT(res, 2, ScalarReal(0));
    UNPROTECT(1);
    return res;
  }

  int    *pi = (int *)    R_alloc(n_pairs, sizeof(int));
  int    *pj = (int *)    R_alloc(n_pairs, sizeof(int));
  double *pw = (double *) R_alloc(n_pairs, sizeof(double));
  double *rm = (double *) R_alloc(n, sizeof(double));
  for (int i = 0; i < n; i++) rm[i] = 0;

  size_t e = 0;
  for (int j = 0; j < n - 1; j++) {
    const double *col = d + pk(n, j + 1, j);
    for (int i = j + 1; i < n; i++) {
      double x = col[i - j - 1];
      if (x <= tt) {
        double w = 1 - x * x * w4;
        pi[e] = i; pj[e] = j; pw[e] = w; e++;
        rm[i] += w; rm[j] += w;
      }
    }
  }

  double s0 = 0;
  int islands = 0;
  for (int i = 0; i < n; i++) {
    s0 += rm[i];
    if (rm[i] == 0) islands++;
  }
  for (int i = 0; i < n; i++) rm[i] /= n;
  double g = s0 / ((double) n * n);


  //--------------------------------------------------
  // B: the trailing block of H A H. With u = 1/sqrt(n)
  // + e1, h = 2 / u'u and a = A e1 (A 1 = 0, so u'a =
  // a_0), H A H = A - u w' - w u' for w = h a - h^2 a_0
  // / 2 u.
  //--------------------------------------------------
  double ui = 1 / sqrt((double) n);
  double u0 = ui + 1;
  double h  = 2 / (u0 * u0 + (n - 1) * ui * ui);

  double *a  = (double *) R_alloc(n, sizeof(double));
  double *wv = (double *) R_alloc(n, sizeof(double));
  for (int i = 0; i < n; i++)
    a[i] = weight(d, n, i, 0, tt, w4) - rm[i] - rm[0] + g;
  for (int i = 0; i < n; i++)
    wv[i] = h * a[i] - h * h * a[0] / 2 * (i == 0 ? u0 : ui);

  dbm.d  = d;  dbm.n  = n;  dbm.m  = m;
  dbm.tt = tt; dbm.w4 = w4; dbm.rm = rm; dbm.g = g;
  dbm.wv = wv; dbm.ui = ui; dbm.k  = 0;

  dbm.B = (double *) R_alloc((size_t) m * m, sizeof(double));
  dbm.Q = (double *) R_alloc((size_t) m * m, sizeof(double));
  run_job(JOB_BUILD, n_threads, 6.0 * m * m);


  //--------------------------------------------------
  // Eigenpairs, sorted by decreasing eigenvalue.
  //--------------------------------------------------
  double *vals  = (double *) R_alloc(m, sizeof(double));
  double *ework = (double *) R_alloc(m, sizeof(double));
  double *betas = (double *) R_alloc(m, sizeof(double));
  double *rc    = (double *) R_alloc(m, sizeof(double));
  double *rs    = (double *) R_alloc(m, sizeof(double));

  int failed = dbm_eigen(vals, ework, betas, rc, rs, n_threads);

  int *ord = (int *) R_alloc(m, sizeof(int));
  for (int i = 0; i < m; i++) ord[i] = i;
  sort_vals = vals;
  qsort(ord, m, sizeof(int), cmp_desc);


  //--------------------------------------------------
  // Which to keep. Eigenvalues within sqrt(eps) of the
  // largest in absolute value are 0.
  //--------------------------------------------------
  double big = 0;
  for (int i = 0; i < m; i++) if (fabs(vals[i]) > big) big = fabs(vals[i]);
  double tol = sqrt(DBL_EPSILON) * big;

  int *sel = (int *) R_alloc(m, sizeof(int));
  int  k   = 0;
  for (int r = 0; r < m; r++) {
    double lam  = vals[ord[r]];
    int    keep = 0;
    switch (autocor) {
      case AUTOCOR_POSITIVE: keep = lam >  tol;       break;
      case AUTOCOR_NEGATIVE: keep = lam < -tol;       break;
      case AUTOCOR_NONNULL:  keep = fabs(lam) > tol;  break;
      case AUTOCOR_ALL:      keep = 1;                break;
    }
    if (keep) sel[k++] = ord[r];
  }

  SEXP sexp_vec    = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_values = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_moran  = PROTECT(allocVector(REALSXP, k));

  for (int c = 0; c < k; c++) {
    double lam = vals[sel[c]];
    REAL(sexp_values)[c] = (fabs(lam) > tol) ? lam / n : 0;
  }

  dbm.sel = sel;  dbm.n_mem = k;
  dbm.out = REAL(sexp_vec); dbm.moran = REAL(sexp_moran);
  dbm.h   = h;    dbm.u0 = u0;
  dbm.s0  = s0;
  dbm.pi  = pi;   dbm.pj = pj; dbm.pw = pw; dbm.n_pairs = n_pairs;
  if (k > 0) run_job(JOB_MEMS, n_threads < k ? n_threads : k, (double) k * (2.0 * n_pairs + 6.0 * n));

  const char *names[] = {
    "vectors", "values", "moran", "threshold", "mst", "n_pairs", "islands", "converged", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 0, sexp_vec);
  SET_VECTOR_ELT(res, 1, sexp_values);
  SET_VECTOR_ELT(res, 2, sexp_moran);
  SET_VECTOR_ELT(res, 3, ScalarReal(threshold));
  SET_VECTOR_ELT(res, 4, ScalarReal(mst));
  SET_VECTOR_ELT(res, 5, ScalarReal((double) n_pairs));
  SET_VECTOR_ELT(res, 6, ScalarInteger(islands));
  SET_VECTOR_ELT(res, 7, ScalarLogical(!failed));

  UNPROTECT(4);
  return res;
}
