// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Uniform Manifold Approximation and Projection (UMAP), from
 * a distance matrix, following the uwot R package's defaults.
 *
 * Stages, each a separate entry point so that R can check it:
 *
 *   C_umap_graph()     Exact k nearest neighbours (including
 *                      self) from the packed dist triangle,
 *                      smooth-kNN calibration (rho, sigma),
 *                      membership strengths, and their fuzzy
 *                      union, as a symmetric sparse matrix.
 *   C_umap_spectral()  Eigenvectors of the graph's normalized
 *                      Laplacian, by thick-restart Lanczos.
 *   C_umap_sgd()       The layout optimization: uwot's batch
 *                      mode (gradients from each epoch's start,
 *                      then Adam), with the same arithmetic
 *                      (single precision coordinates) and the
 *                      same random number scheme, so that given
 *                      the same seeds the result is the same as
 *                      uwot's batch = TRUE.
 *   C_umap_rand()      Uniform or normal deviates from PCG,
 *                      for random starts and jitter.
 *
 * Determinism: the kNN, smooth-kNN, graph and Laplacian
 * products are split into fixed tasks (per point or block of
 * points), each written by one thread with no shared
 * accumulators; reductions run on the main thread in a fixed
 * order. The SGD sums each vertex's gradient within the task
 * that owns the vertex. Results are thus identical for any
 * number of threads.
 */

#include "ecodive.h"

#include <float.h>  // FLT_EPSILON, DBL_MAX
#include <limits.h> // INT_MAX


// Points per kNN task.
#define UMAP_KNN_ROWS 64

// Rows per task in Laplacian products.
#define UMAP_LZ_ROWS 256

// PCG streams, so that different uses of `seed` are
// independent.
#define UMAP_STREAM_UNIF   1
#define UMAP_STREAM_NORM   2
#define UMAP_STREAM_EPOCH  3
#define UMAP_STREAM_LZ     4

// pcg-cpp's default stream, which uwot's pcg32 uses:
// increment 1442695040888963407 = (this << 1) | 1.
#define UWOT_PCG_STREAM 721347520444481703ULL

// smooth-kNN settings, as in uwot::smooth_knn().
#define SKNN_ITER   64
#define SKNN_TOL    1e-5
#define SKNN_MINK   1e-3


//------------------------------------------------------
// Offset of element (j + 1, j) in a dist vector. Column
// j's elements, rows j + 1 to n - 1, follow it in order.
//------------------------------------------------------
static inline size_t col_start(int j, int n) {
  return (size_t) j * (size_t) n - (size_t) j * (size_t) (j + 1) / 2;
}


//======================================================
// Stage 1: k nearest neighbours.
//======================================================

static int           kn_n;     // Points.
static int           kn_k;     // Neighbours, including self.
static const double *kn_d;     // Packed lower triangle.
static int          *kn_idx;   // n x k, row-major, 0-based.
static double       *kn_dist;  // n x k, row-major.


// Is (d1, j1) closer than (d2, j2)? Ties go to the lower
// index, as R's order() does in uwot's dist_nn().
static inline int nbr_less(double d1, int j1, double d2, int j2) {
  return d1 < d2 || (d1 == d2 && j1 < j2);
}


// Offer (d, j) to a max-heap of the k closest so far.
static inline void heap_offer(double *hd, int *hj, int *size, int k, double d, int j) {

  int i;

  if (*size < k) {
    i = (*size)++;
    while (i > 0) {
      int p = (i - 1) / 2;
      if (!nbr_less(hd[p], hj[p], d, j)) break;
      hd[i] = hd[p]; hj[i] = hj[p];
      i = p;
    }
    hd[i] = d; hj[i] = j;
    return;
  }

  if (!nbr_less(d, j, hd[0], hj[0])) return;

  // Replace the root and sift down.
  i = 0;
  for (;;) {
    int c = 2 * i + 1;
    if (c >= k) break;
    if (c + 1 < k && nbr_less(hd[c], hj[c], hd[c + 1], hj[c + 1])) c++;
    if (!nbr_less(d, j, hd[c], hj[c])) break;
    hd[i] = hd[c]; hj[i] = hj[c];
    i = c;
  }
  hd[i] = d; hj[i] = j;
}


// Sort a full max-heap into ascending order, in place.
static void heap_sort(double *hd, int *hj, int k) {

  for (int end = k - 1; end > 0; end--) {

    double d = hd[end]; int j = hj[end];
    hd[end] = hd[0];    hj[end] = hj[0];

    int i = 0;
    for (;;) {
      int c = 2 * i + 1;
      if (c >= end) break;
      if (c + 1 < end && nbr_less(hd[c], hj[c], hd[c + 1], hj[c + 1])) c++;
      if (!nbr_less(d, j, hd[c], hj[c])) break;
      hd[i] = hd[c]; hj[i] = hj[c];
      i = c;
    }
    hd[i] = d; hj[i] = j;
  }
}


//------------------------------------------------------
// Each task finds the neighbours of UMAP_KNN_ROWS points.
// Distances to earlier points come from other columns of
// the packed triangle; they are read as contiguous runs,
// one per column, for all of the task's points at once.
//------------------------------------------------------
static void *knn_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = kn_n;
  int k         = kn_k;
  int n_tasks   = (n + UMAP_KNN_ROWS - 1) / UMAP_KNN_ROWS;

  int sizes[UMAP_KNN_ROWS];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int i0 = task * UMAP_KNN_ROWS;
    int i1 = (i0 + UMAP_KNN_ROWS < n) ? i0 + UMAP_KNN_ROWS : n;
    int nr = i1 - i0;

    // Heaps live in the output rows.
    for (int r = 0; r < nr; r++) sizes[r] = 0;

    #define HEAP(i) kn_dist + (size_t) (i) * k, kn_idx + (size_t) (i) * k, sizes + ((i) - i0), k

    // Columns left of the block.
    for (int j = 0; j < i0; j++) {
      const double *d = kn_d + col_start(j, n) + (i0 - j - 1);
      for (int r = 0; r < nr; r++)
        heap_offer(HEAP(i0 + r), d[r], j);
    }

    for (int i = i0; i < i1; i++) {

      // Earlier points within the block.
      for (int j = i0; j < i; j++)
        heap_offer(HEAP(i), kn_d[col_start(j, n) + (i - j - 1)], j);

      // Self.
      heap_offer(HEAP(i), 0.0, i);

      // Later points: column i, contiguous.
      const double *d = kn_d + col_start(i, n) - (i + 1);
      for (int j = i + 1; j < n; j++)
        heap_offer(HEAP(i), d[j], j);

      heap_sort(kn_dist + (size_t) i * k, kn_idx + (size_t) i * k, k);
    }

    #undef HEAP
  }

  return NULL;
}


//======================================================
// Stage 2: smooth-kNN calibration and membership
// strengths, as in uwot's smooth_knn.h, with uwot's
// defaults: skip_first = TRUE, n_iter = 64, tol = 1e-5,
// min_k_dist_scale = 1e-3.
//======================================================

static double  sk_target;   // log2(k) * bandwidth
static double  sk_local;    // local_connectivity
static double  sk_mean;     // Mean of all neighbour distances.
static double *sk_sigma;    // n
static double *sk_rho;      // n
static double *sk_w;        // n x k, row-major.
static int    *sk_fails;    // n: 1 where the search failed.


// Welford mean in long double, as uwot's mean_average().
static double mean_average(const double *v, size_t len) {
  long double mean = 0.0;
  for (size_t i = 0; i < len; i++)
    mean += (v[i] - mean) / (long double) (i + 1);
  return (double) mean;
}


static void *sknn_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = kn_n;
  int k         = kn_k;

  for (int i = thread_i; i < n; i += n_threads) {

    const double *nd = kn_dist + (size_t) i * k;
    double       *w  = sk_w    + (size_t) i * k;

    // First non-zero distance.
    int nz = k;
    for (int j = 0; j < k; j++)
      if (nd[j] > 0.0) { nz = j; break; }

    // rho
    double rho    = 0.0;
    int    nnzero = k - nz;
    if ((double) nnzero >= sk_local) {
      int    index  = (int) floor(sk_local);
      double interp = sk_local - index;
      if (index > 0) {
        rho = nd[nz + index - 1];
        if (interp >= SKNN_TOL)
          rho += interp * (nd[nz + index] - rho);
      } else if (nnzero > 0) {
        rho = interp * nd[nz];
      }
    } else if (nnzero > 0) {
      rho = nd[k - 1];
    }

    // sigma, by binary search, skipping the first (self).
    double sigma = 1.0, sigma_best = 1.0, adiff_min = DBL_MAX;
    double lo = 0.0, hi = DBL_MAX;
    int    converged = 0;

    for (int iter = 0; iter < SKNN_ITER; iter++) {

      double val = 0.0;
      for (int j = 1; j < k; j++) {
        double rj = nd[j] - rho;
        val += (rj <= 0.0) ? 1.0 : exp(-rj / sigma);
      }

      double adiff = fabs(val - sk_target);
      if (adiff < SKNN_TOL) { converged = 1; break; }

      if (adiff < adiff_min) { adiff_min = adiff; sigma_best = sigma; }

      if (val > sk_target) {
        hi    = sigma;
        sigma = 0.5 * (lo + hi);
      } else {
        lo = sigma;
        if (hi == DBL_MAX) { sigma *= 2;                }
        else               { sigma = 0.5 * (lo + hi);   }
      }
    }

    if (!converged) sigma = sigma_best;
    sk_fails[i] = !converged;

    // Safeguard sigma.
    double floor_sigma = (rho > 0.0)
      ? SKNN_MINK * mean_average(nd, (size_t) k)
      : SKNN_MINK * sk_mean;
    if (floor_sigma > sigma) sigma = floor_sigma;

    for (int j = 0; j < k; j++) {
      double rj = nd[j] - rho;
      w[j] = (rj <= 0.0) ? 1.0 : exp(-rj / sigma);
    }

    sk_sigma[i] = sigma;
    sk_rho[i]   = rho;
  }

  return NULL;
}


//======================================================
// Stage 3: fuzzy union of the membership strengths,
//   mix * (P + P' - P o P') + (1 - mix) * (P o P'),
// as uwot's fuzzy_set_union(), dropping the diagonal and
// zeros. The result is symmetric, stored by rows with
// ascending column indices (the same as uwot's
// column-compressed matrix, read by columns).
//======================================================

static double  un_mix;
static int    *un_rp;     // Reverse lists: n + 1 offsets.
static int    *un_rj;     // Sources, ascending within a list.
static double *un_rw;     // Their weights.
static int    *un_sorted; // n x k: each row's neighbour slots,
                          // ordered by neighbour index.
static int    *un_cnt;    // n: entries in each output row.
static int    *un_p;      // n + 1 output offsets (pass 2).
static int    *un_j;      // Output column indices.
static double *un_x;      // Output values.


// Union value from p_ij = a and p_ji = b, with Matrix's
// order of operations.
static inline double union_val(double a, double b) {
  double xx = a * b;
  if (un_mix == 1) return (a + b) - xx;
  if (un_mix == 0) return xx;
  return un_mix * ((a + b) - xx) + (1 - un_mix) * xx;
}


// Pass 1 (un_p == NULL) counts each row's entries; pass 2
// writes them.
static void *union_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = kn_n;
  int k         = kn_k;

  for (int i = thread_i; i < n; i += n_threads) {

    const int    *idx = kn_idx    + (size_t) i * k;
    const double *w   = sk_w      + (size_t) i * k;
    const int    *ord = un_sorted + (size_t) i * k;

    int a = 0, b = un_rp[i], b_end = un_rp[i + 1];
    int cnt = 0, out = un_p ? un_p[i] : 0;

    while (a < k || b < b_end) {

      int    ja = INT_MAX, jb = INT_MAX;
      double wa = 0,       wb = 0;

      // Skip self and zero weights in the forward list.
      while (a < k && (idx[ord[a]] == i || w[ord[a]] == 0)) a++;
      if (a < k)     { ja = idx[ord[a]]; wa = w[ord[a]]; }
      if (b < b_end) { jb = un_rj[b];    wb = un_rw[b];  }
      if (ja == INT_MAX && jb == INT_MAX) break;

      int    j;
      double pa = 0, pb = 0;
      if (ja <= jb) { j = ja; pa = wa; a++; }
      else          { j = jb;             }
      if (jb == j)  { pb = wb; b++;       }

      double v = union_val(pa, pb);
      if (v != 0) {
        if (un_p) { un_j[out] = j; un_x[out] = v; out++; }
        cnt++;
      }
    }

    un_cnt[i] = cnt;
  }

  return NULL;
}


//======================================================
// R interface. Returns list(idx, dist, sigma, rho,
// n_failures, p, j, x, n_components): the n x k
// neighbour indices (1-based) and distances, each point's
// sigma and rho, and the symmetric fuzzy graph by rows
// (0-based p and j). `sexp_dbl` is c(local_connectivity,
// bandwidth, set_op_mix_ratio).
//======================================================
SEXP C_umap_graph(SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_dbl, SEXP sexp_n_threads) {

  int n         = asInteger(sexp_n);
  int k         = asInteger(sexp_k);
  int n_threads = asInteger(sexp_n_threads);

  if (n < 2 || k < 2 || k > n || n_threads < 1 ||
      (size_t) XLENGTH(sexp_dist) != (size_t) n * (n - 1) / 2 ||
      (size_t) n * k > INT_MAX)
    error("Invalid arguments to C_umap_graph()."); // # nocov

  size_t nk = (size_t) n * k;

  SEXP sexp_idx   = PROTECT(allocMatrix(INTSXP,  n, k));
  SEXP sexp_ndist = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_sigma = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_rho   = PROTECT(allocVector(REALSXP, n));

  kn_n    = n;
  kn_k    = k;
  kn_d    = REAL(sexp_dist);
  kn_idx  = (int *)    R_alloc(nk, sizeof(int));
  kn_dist = (double *) R_alloc(nk, sizeof(double));

  run_parallel(knn_worker, n_threads, n);


  //--------------------------------------------------
  // Smooth kNN.
  //--------------------------------------------------
  sk_local  = REAL(sexp_dbl)[0];
  sk_target = log2((double) k) * REAL(sexp_dbl)[1];
  sk_mean   = mean_average(kn_dist, nk);
  sk_sigma  = REAL(sexp_sigma);
  sk_rho    = REAL(sexp_rho);
  sk_w      = (double *) R_alloc(nk, sizeof(double));
  sk_fails  = (int *)    R_alloc(n,  sizeof(int));

  run_parallel(sknn_worker, n_threads, n);

  int n_fails = 0;
  for (int i = 0; i < n; i++) n_fails += sk_fails[i];


  //--------------------------------------------------
  // Fuzzy union.
  //--------------------------------------------------
  un_mix = REAL(sexp_dbl)[2];

  // Reverse lists: (i -> j, w) filed under j, by i.
  un_rp = (int *) R_alloc(n + 1, sizeof(int));
  memset(un_rp, 0, (n + 1) * sizeof(int));
  for (size_t e = 0; e < nk; e++) {
    int i = (int) (e / k), j = kn_idx[e];
    if (j != i && sk_w[e] != 0) un_rp[j + 1]++;
  }
  for (int i = 0; i < n; i++) un_rp[i + 1] += un_rp[i];

  int nrev = un_rp[n];
  int *fill = (int *) R_alloc(n, sizeof(int));
  memcpy(fill, un_rp, n * sizeof(int));
  un_rj = (int *)    R_alloc(nrev > 0 ? nrev : 1, sizeof(int));
  un_rw = (double *) R_alloc(nrev > 0 ? nrev : 1, sizeof(double));
  for (size_t e = 0; e < nk; e++) {
    int i = (int) (e / k), j = kn_idx[e];
    if (j != i && sk_w[e] != 0) {
      un_rj[fill[j]] = i;
      un_rw[fill[j]] = sk_w[e];
      fill[j]++;
    }
  }

  // Each row's neighbour slots, by neighbour index
  // (insertion sort; k is small).
  un_sorted = (int *) R_alloc(nk, sizeof(int));
  for (int i = 0; i < n; i++) {
    const int *idx = kn_idx    + (size_t) i * k;
    int       *ord = un_sorted + (size_t) i * k;
    for (int a = 0; a < k; a++) {
      int b = a;
      while (b > 0 && idx[ord[b - 1]] > idx[a]) { ord[b] = ord[b - 1]; b--; }
      ord[b] = a;
    }
  }

  un_cnt = (int *) R_alloc(n, sizeof(int));
  un_p   = NULL;
  run_parallel(union_worker, n_threads, n);

  SEXP sexp_p = PROTECT(allocVector(INTSXP, n + 1));
  int *p = INTEGER(sexp_p);
  p[0] = 0;
  for (int i = 0; i < n; i++) p[i + 1] = p[i] + un_cnt[i];

  SEXP sexp_j = PROTECT(allocVector(INTSXP,  p[n]));
  SEXP sexp_x = PROTECT(allocVector(REALSXP, p[n]));
  un_p = p;
  un_j = INTEGER(sexp_j);
  un_x = REAL(sexp_x);
  run_parallel(union_worker, n_threads, n);
  un_p = NULL;


  //--------------------------------------------------
  // Connected components, by union-find.
  //--------------------------------------------------
  int *parent = fill; // reuse
  for (int i = 0; i < n; i++) parent[i] = i;
  for (int i = 0; i < n; i++) {
    for (int e = p[i]; e < p[i + 1]; e++) {
      int a = i, b = un_j[e];
      while (parent[a] != a) a = parent[a] = parent[parent[a]];
      while (parent[b] != b) b = parent[b] = parent[parent[b]];
      if (a != b) { if (a < b) parent[b] = a; else parent[a] = b; }
    }
  }
  int n_comp = 0;
  for (int i = 0; i < n; i++) if (parent[i] == i) n_comp++;


  //--------------------------------------------------
  // Neighbours to R's column-major layout, 1-based.
  //--------------------------------------------------
  int    *ri = INTEGER(sexp_idx);
  double *rd = REAL(sexp_ndist);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < k; j++) {
      ri[i + (size_t) j * n] = kn_idx [(size_t) i * k + j] + 1;
      rd[i + (size_t) j * n] = kn_dist[(size_t) i * k + j];
    }

  const char *names[] = {
    "idx", "dist", "sigma", "rho", "n_failures", "p", "j", "x", "n_components", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_idx);
  SET_VECTOR_ELT(sexp_result, 1, sexp_ndist);
  SET_VECTOR_ELT(sexp_result, 2, sexp_sigma);
  SET_VECTOR_ELT(sexp_result, 3, sexp_rho);
  SET_VECTOR_ELT(sexp_result, 4, ScalarInteger(n_fails));
  SET_VECTOR_ELT(sexp_result, 5, sexp_p);
  SET_VECTOR_ELT(sexp_result, 6, sexp_j);
  SET_VECTOR_ELT(sexp_result, 7, sexp_x);
  SET_VECTOR_ELT(sexp_result, 8, ScalarInteger(n_comp));

  UNPROTECT(8);
  return sexp_result;
}



//======================================================
// Stage 4: spectral initialization.
//
// uwot takes the eigenvectors of the k + 1 smallest
// eigenvalues of L = I - D^-1/2 A D^-1/2 (D = A's column
// sums), and drops the first, constant one (its
// eigenvector is D^1/2 1). Those are the eigenvectors of
// the k largest eigenvalues of M = D^-1/2 A D^-1/2 in the
// complement of D^1/2 1, which are found here by Lanczos
// iteration with thick restarts (Wu & Simon 2000) and full
// reorthogonalization.
//======================================================

#define LZ_TOL       1e-10
#define LZ_MAX_ITER  5000   // Restarts.

static int           lz_n;
static const int    *lz_p;
static const int    *lz_j;
static const double *lz_x;
static const double *lz_s;    // D^-1/2
static const double *lz_in;
static double       *lz_out;


static void *lz_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = lz_n;
  int n_tasks   = (n + UMAP_LZ_ROWS - 1) / UMAP_LZ_ROWS;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int i1 = (task + 1) * UMAP_LZ_ROWS;
    if (i1 > n) i1 = n;
    for (int i = task * UMAP_LZ_ROWS; i < i1; i++) {
      double s = 0;
      for (int e = lz_p[i]; e < lz_p[i + 1]; e++)
        s += lz_x[e] * lz_s[lz_j[e]] * lz_in[lz_j[e]];
      lz_out[i] = lz_s[i] * s;
    }
  }

  return NULL;
}


static void lz_mult(const double *in, double *out, int n_threads) {
  lz_in  = in;
  lz_out = out;
  int n_tasks = (lz_n + UMAP_LZ_ROWS - 1) / UMAP_LZ_ROWS;
  run_parallel(lz_worker, n_tasks >= 4 ? n_threads : 1, 100);
}


static double dot(const double *a, const double *b, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) s += a[i] * b[i];
  return s;
}


// w -= (v . w) v for unit vectors v[0..nv-1] (twice, for
// numerical orthogonality). Adds the coefficients to h, if
// given.
static void orthogonalize(double *w, double **v, int nv, double *h, int n) {
  for (int pass = 0; pass < 2; pass++)
    for (int i = 0; i < nv; i++) {
      double c = dot(v[i], w, n);
      for (int r = 0; r < n; r++) w[r] -= c * v[i][r];
      if (h) h[i] += c;
    }
}


// Eigen-decomposition of a symmetric m x m matrix `a`
// (destroyed) by cyclic Jacobi rotations. Eigenvalues to
// `val`, unit eigenvectors to the columns of `vec`, both
// sorted by decreasing eigenvalue.
static void jacobi_eigen(double *a, int m, double *val, double *vec) {

  for (int i = 0; i < m * m; i++) vec[i] = 0;
  for (int i = 0; i < m; i++) vec[i * m + i] = 1;

  for (int sweep = 0; sweep < 100; sweep++) {

    double off = 0, tot = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        tot += a[i * m + j] * a[i * m + j];
        if (i != j) off += a[i * m + j] * a[i * m + j];
      }
    if (off <= 1e-30 * tot || off == 0) break;

    for (int p = 0; p < m - 1; p++)
      for (int q = p + 1; q < m; q++) {

        double apq = a[p * m + q];
        if (apq == 0) continue;

        double app   = a[p * m + p], aqq = a[q * m + q];
        double theta = (aqq - app) / (2 * apq);
        double t     = (theta >= 0 ? 1 : -1) / (fabs(theta) + sqrt(theta * theta + 1));
        double c     = 1 / sqrt(t * t + 1), s = t * c;

        for (int r = 0; r < m; r++) {
          double arp = a[r * m + p], arq = a[r * m + q];
          a[r * m + p] = c * arp - s * arq;
          a[r * m + q] = s * arp + c * arq;
        }
        for (int r = 0; r < m; r++) {
          double apr = a[p * m + r], aqr = a[q * m + r];
          a[p * m + r] = c * apr - s * aqr;
          a[q * m + r] = s * apr + c * aqr;
        }
        for (int r = 0; r < m; r++) {
          double vrp = vec[r * m + p], vrq = vec[r * m + q];
          vec[r * m + p] = c * vrp - s * vrq;
          vec[r * m + q] = s * vrp + c * vrq;
        }
      }
  }

  for (int i = 0; i < m; i++) val[i] = a[i * m + i];

  // Selection sort, decreasing, with stable ties.
  for (int i = 0; i < m; i++) {
    int best = i;
    for (int j = i + 1; j < m; j++) if (val[j] > val[best]) best = j;
    if (best != i) {
      double t = val[i]; val[i] = val[best]; val[best] = t;
      for (int r = 0; r < m; r++) {
        t = vec[r * m + i]; vec[r * m + i] = vec[r * m + best]; vec[r * m + best] = t;
      }
    }
  }
}


// A random unit vector orthogonal to v[0..nv-1].
static void random_start(double *w, double **v, int nv, int n, pcg32_random_t *rng) {
  for (int tries = 0; tries < 10; tries++) {
    for (int r = 0; r < n; r++) w[r] = pcg32_random_r(rng) / 4294967296.0 - 0.5;
    orthogonalize(w, v, nv, NULL, n);
    double nrm = sqrt(dot(w, w, n));
    if (nrm > 1e-8) {
      for (int r = 0; r < n; r++) w[r] /= nrm;
      return;
    }
  }
  memset(w, 0, n * sizeof(double)); // # nocov
}


//======================================================
// R interface. `sexp_graph` is C_umap_graph()'s result.
// Returns list(vectors, values, converged): the n x k
// unit eigenvectors of L for its 2nd to (k + 1)th
// smallest eigenvalues, in that order, and those
// eigenvalues. Each vector's largest coordinate in
// absolute value is positive.
//======================================================
SEXP C_umap_spectral(SEXP sexp_graph, SEXP sexp_k, SEXP sexp_n_threads) {

  SEXP sexp_p = VECTOR_ELT(sexp_graph, 5);
  int  n         = (int) XLENGTH(sexp_p) - 1;
  int  nev       = asInteger(sexp_k);
  int  n_threads = asInteger(sexp_n_threads);

  if (n < 3 || nev < 1 || nev > n - 2 || n_threads < 1)
    error("Invalid arguments to C_umap_spectral()."); // # nocov

  lz_n = n;
  lz_p = INTEGER(sexp_p);
  lz_j = INTEGER(VECTOR_ELT(sexp_graph, 6));
  lz_x = REAL(VECTOR_ELT(sexp_graph, 7));

  // Degrees (column sums = row sums, as A is symmetric).
  double *s  = (double *) R_alloc(n, sizeof(double));
  double *v0 = (double *) R_alloc(n, sizeof(double));
  for (int i = 0; i < n; i++) {
    double d = 0;
    for (int e = lz_p[i]; e < lz_p[i + 1]; e++) d += lz_x[e];
    s[i]  = (d > 0) ? 1 / sqrt(d) : 0;
    v0[i] = sqrt(d);
  }
  lz_s = s;
  double nrm0 = sqrt(dot(v0, v0, n));
  for (int i = 0; i < n; i++) v0[i] /= nrm0;

  // Basis size.
  int m = 2 * nev + 20;
  if (m > n - 1) m = n - 1;

  // basis[0] = v0 (deflated), basis[1..m+1] = Krylov
  // vectors V[0..m].
  double **basis = (double **) R_alloc(m + 2, sizeof(double *));
  basis[0] = v0;
  for (int i = 1; i < m + 2; i++)
    basis[i] = (double *) R_alloc(n, sizeof(double));
  double **V = basis + 1;

  double *w    = (double *) R_alloc(n, sizeof(double));
  double *T    = (double *) R_alloc((size_t) m * m, sizeof(double));
  double *Tc   = (double *) R_alloc((size_t) m * m, sizeof(double));
  double *S    = (double *) R_alloc((size_t) m * m, sizeof(double));
  double *th   = (double *) R_alloc(m, sizeof(double));
  double *h    = (double *) R_alloc(m + 1, sizeof(double));
  double *tmp  = (double *) R_alloc((size_t) n * m, sizeof(double));

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 0, UMAP_STREAM_LZ);

  memset(T, 0, (size_t) m * m * sizeof(double));
  random_start(V[0], basis, 1, n, &rng);

  int    p         = 0; // Kept Ritz vectors.
  int    converged = 0;
  double beta      = 0;

  for (int iter = 0; iter < LZ_MAX_ITER; iter++) {

    for (int j = p; j < m; j++) {

      lz_mult(V[j], w, n_threads);

      // Coefficients against v0 are dropped; those against
      // V[0..j] are column j of T = V' M V.
      memset(h, 0, (j + 1) * sizeof(double));
      orthogonalize(w, basis, 1, NULL, n);
      orthogonalize(w, V, j + 1, h, n);
      orthogonalize(w, basis, 1, NULL, n);

      for (int i = 0; i <= j; i++)
        T[i * m + j] = T[j * m + i] = h[i];

      beta = sqrt(dot(w, w, n));

      if (beta > 1e-12) {
        for (int r = 0; r < n; r++) V[j + 1][r] = w[r] / beta;
      } else {
        // Invariant subspace: continue from a fresh vector.
        beta = 0;
        random_start(V[j + 1], basis, j + 2, n, &rng);
      }
    }

    memcpy(Tc, T, (size_t) m * m * sizeof(double));
    jacobi_eigen(Tc, m, th, S);

    // Converged when each wanted residual, |beta * S[m-1, i]|,
    // is small.
    converged = 1;
    for (int i = 0; i < nev; i++)
      if (fabs(beta * S[(m - 1) * m + i]) > LZ_TOL) converged = 0;

    if (converged || m >= n - 1 || iter == LZ_MAX_ITER - 1) break;

    if (iter % 10 == 9) R_CheckUserInterrupt();

    // Thick restart: keep the best p Ritz vectors.
    p = nev + (m - nev) / 2;
    if (p > m - 1) p = m - 1;

    memset(tmp, 0, (size_t) n * p * sizeof(double));
    for (int i = 0; i < p; i++)
      for (int j = 0; j < m; j++) {
        double c = S[j * m + i];
        double *t = tmp + (size_t) i * n;
        for (int r = 0; r < n; r++) t[r] += c * V[j][r];
      }
    for (int i = 0; i < p; i++) memcpy(V[i], tmp + (size_t) i * n, n * sizeof(double));
    memcpy(V[p], V[m], n * sizeof(double));

    memset(T, 0, (size_t) m * m * sizeof(double));
    for (int i = 0; i < p; i++) T[i * m + i] = th[i];
  }

  if (m >= n - 1) converged = 1;


  //--------------------------------------------------
  // Ritz vectors, with the sign convention.
  //--------------------------------------------------
  SEXP sexp_vec = PROTECT(allocMatrix(REALSXP, n, nev));
  SEXP sexp_val = PROTECT(allocVector(REALSXP, nev));
  double *vec = REAL(sexp_vec);

  for (int i = 0; i < nev; i++) {

    double *out = vec + (size_t) i * n;
    memset(out, 0, n * sizeof(double));
    for (int j = 0; j < m; j++) {
      double c = S[j * m + i];
      for (int r = 0; r < n; r++) out[r] += c * V[j][r];
    }

    double nrm = sqrt(dot(out, out, n));
    int    big = 0;
    for (int r = 0; r < n; r++) {
      out[r] /= nrm;
      if (fabs(out[r]) > fabs(out[big])) big = r;
    }
    if (out[big] < 0)
      for (int r = 0; r < n; r++) out[r] = -out[r];

    REAL(sexp_val)[i] = 1 - th[i];
  }

  const char *names[] = { "vectors", "values", "converged", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_vec);
  SET_VECTOR_ELT(sexp_result, 1, sexp_val);
  SET_VECTOR_ELT(sexp_result, 2, ScalarLogical(converged));

  UNPROTECT(3);
  return sexp_result;
}



//======================================================
// Stage 5: layout optimization.
//
// uwot's optimize_layout() in batch mode (uwot::umap(...,
// batch = TRUE)), with its defaults. Each epoch, every
// vertex sums the gradient of its own edges, computed from
// the coordinates as they were at the start of the epoch:
// attraction along each sampled edge (twice its gradient,
// as the graph is symmetric and each edge is listed under
// both of its ends) and repulsion from negative samples.
// Gradients are clipped to [-4, 4]. At the end of the
// epoch, Adam moves every coordinate, with a learning rate
// that decays linearly. The arithmetic follows uwot's:
// coordinates, gradients and sampling schedules in single
// precision, and the same mixed float/double expressions.
//
// Random numbers, as uwot's batch mode: each epoch draws
// two 32-bit seeds per vertex, and each vertex draws its
// negative samples from its own pcg32 stream, seeded with
// ((seed[2i] << 32) | seed[2i + 1]) on pcg-cpp's default
// stream. uwot draws the seeds from R's RNG; here they come
// from `seed` via PCG, or are given.
//
// Determinism: a vertex's gradient depends only on the
// start-of-epoch coordinates and its own seeds, and is
// written only by the task that owns the vertex, so the
// result is the same for any number of threads. Tasks are
// idempotent (they read the sampling schedule from one
// buffer and write it to another), as run_parallel() may
// rerun the work on one thread.
//======================================================

// Vertices per SGD task.
#define UMAP_SGD_ROWS 64

static inline uint32_t pcg_bounded(pcg32_random_t *rng, uint32_t bound) {
  uint32_t threshold = (uint32_t) (-bound) % bound;
  for (;;) {
    uint32_t r = pcg32_random_r(rng);
    if (r >= threshold) return r % bound;
  }
}

static inline float clamp4(float g) {
  return (g < -4.0f) ? -4.0f : (4.0f < g) ? 4.0f : g;
}


static int             sg_n;
static int             sg_ndim;
static const int      *sg_ptr;    // n + 1 edge offsets, by head.
static const int      *sg_tail;   // Each edge's tail.
static const float    *sg_eps;    // Epochs per sample.
static const float    *sg_epns;   // Epochs per negative sample.
static const float    *sg_eons0;  // Epoch of next sample (in).
static const float    *sg_eonn0;  // Epoch of next negative sample (in).
static float          *sg_eons1;  // (out)
static float          *sg_eonn1;  // (out)
static const uint32_t *sg_seeds;  // 2 per vertex.
static const float    *sg_x;      // n x ndim, point-major.
static float          *sg_grad;   // n x ndim, point-major.
static float          *sg_disp;   // ndim per thread.
static float           sg_epoch;
static float           sg_a, sg_b, sg_a_b_m2, sg_gamma_b_2;


static void *sgd_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = sg_n;
  int ndim      = sg_ndim;
  int n_tasks   = (n + UMAP_SGD_ROWS - 1) / UMAP_SGD_ROWS;

  const float *x    = sg_x;
  float       *disp = sg_disp + (size_t) thread_i * ndim;
  float        a    = sg_a, b = sg_b, fepoch = sg_epoch;

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int p1 = (task + 1) * UMAP_SGD_ROWS;
    if (p1 > n) p1 = n;

    for (int p = task * UMAP_SGD_ROWS; p < p1; p++) {

      size_t dj   = (size_t) ndim * p;
      float *grad = sg_grad + dj;
      for (int d = 0; d < ndim; d++) grad[d] = 0.0f;

      pcg32_random_t rng;
      pcg32_srandom_r(
        &rng, ((uint64_t) sg_seeds[2 * (size_t) p] << 32) | sg_seeds[2 * (size_t) p + 1],
        UWOT_PCG_STREAM );

      for (int e = sg_ptr[p]; e < sg_ptr[p + 1]; e++) {

        if (!(sg_eons0[e] <= fepoch)) {
          sg_eons1[e] = sg_eons0[e];
          sg_eonn1[e] = sg_eonn0[e];
          continue;
        }

        size_t dk = (size_t) ndim * sg_tail[e];

        // Attraction.
        float d2 = 0.0f;
        for (int d = 0; d < ndim; d++) {
          float diff = x[dj + d] - x[dk + d];
          disp[d] = diff;
          d2 += diff * diff;
        }
        if (!(FLT_EPSILON < d2)) d2 = FLT_EPSILON;

        float pd2b  = powf(d2, b);
        float coeff = (float) ((sg_a_b_m2 * pd2b) / (d2 * (a * pd2b + 1.0)));

        for (int d = 0; d < ndim; d++)
          grad[d] = (float) (grad[d] + 2.0 * clamp4(coeff * disp[d]));

        // Repulsion from negative samples.
        size_t n_neg = (size_t) ((fepoch - sg_eonn0[e]) / sg_epns[e]);

        for (size_t q = 0; q < n_neg; q++) {

          size_t dkn = (size_t) pcg_bounded(&rng, (uint32_t) n) * ndim;
          if (dj == dkn) continue;

          d2 = 0.0f;
          for (int d = 0; d < ndim; d++) {
            float diff = x[dj + d] - x[dkn + d];
            disp[d] = diff;
            d2 += diff * diff;
          }
          if (!(FLT_EPSILON < d2)) d2 = FLT_EPSILON;

          coeff = (float) (sg_gamma_b_2 / ((0.001 + d2) * (a * powf(d2, b) + 1.0)));

          for (int d = 0; d < ndim; d++)
            grad[d] += clamp4(coeff * disp[d]);
        }

        sg_eons1[e] = sg_eons0[e] + sg_eps[e];
        sg_eonn1[e] = sg_eonn0[e] + (float) n_neg * sg_epns[e];
      }
    }
  }

  return NULL;
}


//======================================================
// R interface. `sexp_graph` is C_umap_graph()'s result,
// `sexp_init` the n x k starting coordinates, `sexp_dbl`
// c(a, b, gamma, learning_rate, negative_sample_rate,
// beta1, beta2, eps) (the last three for Adam),
// `sexp_seeds` either an integer seed or the
// 2 * n * n_epochs per-vertex seeds (doubles in
// [0, 2^32)), epoch by epoch. Returns the n x k
// coordinates (not centered).
//======================================================
SEXP C_umap_sgd(
    SEXP sexp_graph, SEXP sexp_init, SEXP sexp_n_epochs, SEXP sexp_dbl,
    SEXP sexp_seeds, SEXP sexp_n_threads ) {

  SEXP sexp_p = VECTOR_ELT(sexp_graph, 5);
  const int    *gp = INTEGER(sexp_p);
  const int    *gj = INTEGER(VECTOR_ELT(sexp_graph, 6));
  const double *gx = REAL(VECTOR_ELT(sexp_graph, 7));

  int n         = (int) XLENGTH(sexp_p) - 1;
  int n_epochs  = asInteger(sexp_n_epochs);
  int n_threads = asInteger(sexp_n_threads);
  int ndim      = (n > 0) ? (int) (XLENGTH(sexp_init) / n) : 0;

  if (n < 1 || ndim < 1 || n_epochs < 0 || n_threads < 1 ||
      (size_t) XLENGTH(sexp_init) != (size_t) n * ndim ||
      (isReal(sexp_seeds) && (size_t) XLENGTH(sexp_seeds) != 2 * (size_t) n * n_epochs))
    error("Invalid arguments to C_umap_sgd()."); // # nocov

  float a      = (float) REAL(sexp_dbl)[0];
  float b      = (float) REAL(sexp_dbl)[1];
  float gamma  = (float) REAL(sexp_dbl)[2];
  float alpha0 = (float) REAL(sexp_dbl)[3];
  float nsr_in = (float) REAL(sexp_dbl)[4];
  float beta1  = (float) REAL(sexp_dbl)[5];
  float beta2  = (float) REAL(sexp_dbl)[6];
  float eps    = (float) REAL(sexp_dbl)[7];


  //--------------------------------------------------
  // Edges, by head, as uwot: drop weights below
  // max / n_epochs; epochs_per_sample =
  // n_epochs / (n_epochs * w / max).
  //--------------------------------------------------
  int nnz = gp[n];
  double wmax = 0;
  for (int e = 0; e < nnz; e++) if (gx[e] > wmax) wmax = gx[e];

  size_t n_edges = 0;
  for (int e = 0; e < nnz; e++)
    if (!(gx[e] < wmax / n_epochs)) n_edges++;

  size_t ea = n_edges > 0 ? n_edges : 1;
  int   *ptr  = (int *)   R_alloc(n + 1, sizeof(int));
  int   *tail = (int *)   R_alloc(ea, sizeof(int));
  float *eps_ = (float *) R_alloc(ea, sizeof(float));
  float *epns = (float *) R_alloc(ea, sizeof(float));
  float *eons[2], *eonn[2];
  for (int t = 0; t < 2; t++) {
    eons[t] = (float *) R_alloc(ea, sizeof(float));
    eonn[t] = (float *) R_alloc(ea, sizeof(float));
  }

  float nsr = (float) (1.0 / nsr_in);

  int ne = 0;
  ptr[0] = 0;
  for (int r = 0; r < n; r++) {
    for (int e = gp[r]; e < gp[r + 1]; e++) {
      if (gx[e] < wmax / n_epochs) continue;
      double n_samples = n_epochs * (gx[e] / wmax);
      tail[ne]    = gj[e];
      eps_[ne]    = (float) (n_epochs / n_samples);
      epns[ne]    = eps_[ne] * nsr;
      eons[0][ne] = eps_[ne];
      eonn[0][ne] = epns[ne];
      ne++;
    }
    ptr[r + 1] = ne;
  }


  //--------------------------------------------------
  // Coordinates and Adam's state: point-major, single
  // precision.
  //--------------------------------------------------
  size_t nk = (size_t) n * ndim;

  const double *init = REAL(sexp_init);
  float *x    = (float *) R_alloc(nk, sizeof(float));
  float *grad = (float *) R_alloc(nk, sizeof(float));
  float *mt   = (float *) R_alloc(nk, sizeof(float));
  float *vt   = (float *) R_alloc(nk, sizeof(float));
  for (int i = 0; i < n; i++)
    for (int d = 0; d < ndim; d++)
      x[(size_t) i * ndim + d] = (float) init[i + (size_t) d * n];
  memset(mt, 0, nk * sizeof(float));
  memset(vt, 0, nk * sizeof(float));

  uint32_t *seeds = (uint32_t *) R_alloc(2 * (size_t) n, sizeof(uint32_t));
  const double *given = isReal(sexp_seeds) ? REAL(sexp_seeds) : NULL;

  pcg32_random_t master;
  if (!given)
    pcg32_srandom_r(&master, (uint64_t) (uint32_t) asInteger(sexp_seeds), UMAP_STREAM_EPOCH);

  // Adam, as uwot's: debiasing folded into ad_scale and epsc.
  float beta11   = (float) (1.0 - beta1);
  float beta21   = (float) (1.0 - beta2);
  float beta1t   = beta1;
  float beta2t   = beta2;
  float epsc     = eps * sqrtf(beta21);
  float ad_scale = alpha0 * sqrtf(beta21) / beta11;

  sg_n         = n;
  sg_ndim      = ndim;
  sg_ptr       = ptr;
  sg_tail      = tail;
  sg_eps       = eps_;
  sg_epns      = epns;
  sg_seeds     = seeds;
  sg_x         = x;
  sg_grad      = grad;
  sg_disp      = (float *) R_alloc((size_t) n_threads * ndim, sizeof(float));
  sg_a         = a;
  sg_b         = b;
  sg_a_b_m2    = (float) (-2.0 * a * b);
  sg_gamma_b_2 = (float) (2.0 * gamma * b);

  for (int epoch = 0; epoch < n_epochs; epoch++) {

    for (size_t i = 0; i < 2 * (size_t) n; i++)
      seeds[i] = given
        ? (uint32_t) given[(size_t) epoch * 2 * n + i]
        : pcg32_random_r(&master);

    sg_epoch = (float) epoch;
    sg_eons0 = eons[epoch & 1]; sg_eons1 = eons[(epoch + 1) & 1];
    sg_eonn0 = eonn[epoch & 1]; sg_eonn1 = eonn[(epoch + 1) & 1];

    run_parallel(sgd_worker, n_threads, n);

    for (size_t i = 0; i < nk; i++) {
      vt[i] += beta21 * (grad[i] * grad[i] - vt[i]);
      mt[i] += beta11 * (grad[i] - mt[i]);
      x[i]  += ad_scale * mt[i] / (sqrtf(vt[i]) + epsc);
    }

    float alpha = (float) (alpha0 * (1.0 - ((float) epoch / (float) n_epochs)));
    beta1t *= beta1;
    beta2t *= beta2;
    float sqrt_b2t1 = (float) sqrt(1.0 - beta2t);
    ad_scale = (float) (alpha * sqrt_b2t1 / (1.0 - beta1t));
    epsc     = sqrt_b2t1 * eps;

    if (epoch % 10 == 9) R_CheckUserInterrupt();
  }

  SEXP sexp_result = PROTECT(allocMatrix(REALSXP, n, ndim));
  double *res = REAL(sexp_result);
  for (int i = 0; i < n; i++)
    for (int d = 0; d < ndim; d++)
      res[i + (size_t) d * n] = x[(size_t) i * ndim + d];

  UNPROTECT(1);
  return sexp_result;
}



// A 53-bit uniform deviate on [0, 1).
static inline double unif53(pcg32_random_t *rng) {
  uint32_t hi = pcg32_random_r(rng) >> 5;
  uint32_t lo = pcg32_random_r(rng) >> 6;
  return (hi * 67108864.0 + lo) / 9007199254740992.0;
}


//======================================================
// R interface. `len` deviates from PCG stream `seed`:
// uniform on [-10, 10) for kind 0, standard normal (by
// Box-Muller) for kind 1.
//======================================================
SEXP C_umap_rand(SEXP sexp_len, SEXP sexp_kind, SEXP sexp_seed) {

  R_xlen_t len  = (R_xlen_t) asReal(sexp_len);
  int      kind = asInteger(sexp_kind);

  SEXP sexp_result = PROTECT(allocVector(REALSXP, len));
  double *res = REAL(sexp_result);

  pcg32_random_t rng;
  pcg32_srandom_r(
    &rng, (uint64_t) (uint32_t) asInteger(sexp_seed),
    kind == 0 ? UMAP_STREAM_UNIF : UMAP_STREAM_NORM );

  for (R_xlen_t i = 0; i < len; i++) {
    if (kind == 0) {
      res[i] = 20 * unif53(&rng) - 10;
    } else {
      double u1 = 1 - unif53(&rng); // (0, 1]
      double u2 = unif53(&rng);
      double r  = sqrt(-2 * log(u1));
      res[i] = r * cos(2 * M_PI * u2);
      if (++i < len) res[i] = r * sin(2 * M_PI * u2);
    }
  }

  UNPROTECT(1);
  return sexp_result;
}
