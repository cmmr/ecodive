// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Parts of this file are a port of scikit-learn 1.9.1's t-SNE
// (sklearn/manifold/_t_sne.py, sklearn/manifold/_utils.pyx,
// sklearn/manifold/_barnes_hut_tsne.pyx and
// sklearn/neighbors/_quad_tree.pyx), which is distributed under the
// following license:
//
//   BSD 3-Clause License
//
//   Copyright (c) 2007-2026 The scikit-learn developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in
//     the documentation and/or other materials provided with the
//     distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//   POSSIBILITY OF SUCH DAMAGE.

/*
 * t-distributed Stochastic Neighbor Embedding (t-SNE).
 *
 * Input affinities (C_tsne_p). Each sample's conditional
 * probabilities p_{j|i} are Gaussian in the squared distance,
 * with a precision found by the same binary search on the
 * entropy as scikit-learn's _binary_search_perplexity(). Then
 * P = (P_cond + P_cond') / sum.
 *   - Dense ("exact"): over all other samples. Stored as a
 *     packed lower triangle, like a dist object. Rows of the
 *     distance matrix are read from the dist object's packed
 *     lower triangle a block at a time, as contiguous runs.
 *   - Sparse: over each sample's K nearest neighbors only,
 *     found exactly from the distances by knn_from_dist()
 *     (knn.c; ties go to the lower index). Stored as CSR.
 * Distances are multiplied by a power of 2 (from R) so the
 * largest is about 1, and a row whose search doesn't converge
 * in scikit-learn's 100 steps is searched again (search_row()).
 *
 * Gradient (gradient()). With w_ij = 1 / (1 + |y_i - y_j|^2),
 * Z = sum_{i != j} w_ij and q_ij = w_ij / Z,
 *
 *   dC/dy_i = 4 sum_j (a p_ij - q_ij) w_ij (y_i - y_j)
 *
 * where a is the early exaggeration factor (1 after the early
 * phase).
 *   - Exact: every pair, from the dense P.
 *   - Barnes-Hut: the attractive (p) part over P's nonzeros; the
 *     repulsive (q) part from a quadtree (k = 2), octree (k = 3)
 *     or binary tree (k = 1), summarizing a cell by its center
 *     of mass when (cell width)^2 / (distance)^2 < theta^2.
 * Each gradient (with its KL divergence terms) is one parallel
 * pass over the samples.
 *
 * Optimization (C_tsne_run). scikit-learn's _gradient_descent()
 * schedule: gradient descent with momentum and per-parameter
 * gains; momentum 0.5 with exaggerated P for the first
 * iterations, then 0.8; progress checked every 50 iterations.
 *
 * Threads come from the persistent pool (run_parallel() in
 * parallel.c), so are cheap to wake: a step uses them when it
 * has about TS_MIN_WORK ns of work per thread.
 *
 * Determinism. Results are bit-identical for any number of
 * threads. Every per-sample quantity is computed by one thread,
 * in a fixed order, into its own output slot; sums over samples
 * (Z, the KL divergence, the sum of P, the gradient norm) are
 * taken afterwards on the main thread, in sample order. The
 * tree is built on the main thread. Ties between equal
 * distances are broken by sample index.
 */

#include "ecodive.h"

#include <float.h>  // DBL_EPSILON, DBL_MAX
#include <limits.h> // INT_MAX

#if defined(__GNUC__) || defined(__clang__)
#  define TS_INLINE static inline __attribute__((always_inline))
#else
#  define TS_INLINE static inline
#endif


// Rows of the distance matrix gathered at once by a thread.
#define TS_ROWS 32

// Samples per parallel task in the gradient.
#define TS_BLOCK 64

// Most halvings of a tree cell. Points still together after
// this many are kept in one leaf and handled pairwise.
#define BH_MAX_DEPTH 200

// scikit-learn's constants. The perplexity search's are
// single-precision floats there.
#define TS_PERP_TOL   ((double) 1e-5f)
#define TS_PERP_EPS   ((double) 1e-8f)
#define TS_PERP_STEPS 100
#define TS_MACHINE_EPSILON DBL_EPSILON


// Least work per thread (roughly in nanoseconds on one
// thread) for a step to use threads. Pool threads are cheap
// to wake: with this gate, threads are never slower than
// one, even for 60 samples.
#define TS_MIN_WORK 1e4

// Threads to use for `work`, at most one per task.
static int use_threads(int n_threads, int n_tasks, double work) {
  int t = threads_for(n_threads, work, TS_MIN_WORK);
  return (t < n_tasks) ? t : (n_tasks > 1 ? n_tasks : 1);
}

// run_parallel() starts threads only when told of at least
// 100 tasks; the workers here split their own tasks.
static void run_tasks(pthread_func_t func, int n_threads) {
  run_parallel(func, n_threads, (n_threads > 1) ? 100 : 1);
}



//======================================================
// Input affinities
//======================================================

static int           pp_n;        // Samples.
static int           pp_K;        // Neighbors, or 0 for all (dense).
static int           pp_f32;      // Round squared distances to float.
static double        pp_scale;    // Multiplies each distance (a power of 2).
static const double *pp_dist;     // Packed lower triangle of D.
static double        pp_target;   // Desired entropy: log(perplexity).
static double       *pp_rowbuf;   // Per thread: TS_ROWS x n (dense).
static double       *pp_sqbuf;    // Per thread: n (dense) or K.
static double       *pp_pbuf;     // Per thread: n (dense).
static double       *pp_beta;     // n: precision used for each row.
static double       *pp_sum;      // n: normalizer used for each row.
static double       *pp_shift;    // n: shift used for each row.
static double       *pp_rowsum;   // n: sum of each row's p_{j|i}.
static int          *pp_nbr;      // n x K: neighbors, ascending.
static double       *pp_cond;     // n x K: their p_{j|i}.
static int           pp_threads;
static int           pp_task0;    // Tasks (blocks of rows or columns)
static int           pp_task1;    // ... in this batch.


//------------------------------------------------------
// buf[r * n + j] = d(i0 + r, j) for r < nr, and every j.
//------------------------------------------------------
static void gather_rows(int i0, int nr, double *buf) {

  int n  = pp_n;
  int i1 = i0 + nr;

  // Columns left of the block: contiguous runs.
  for (int j = 0; j < i0; j++) {
    const double *d = pp_dist + dist_col_start(j, n) + (i0 - j - 1);
    for (int r = 0; r < nr; r++)
      buf[(size_t) r * n + j] = d[r];
  }

  // The block's own columns.
  for (int j = i0; j < i1; j++) {
    const double *d  = pp_dist + dist_col_start(j, n);  // d[t]: row j + 1 + t
    double       *bj = buf + (size_t) (j - i0) * n;
    bj[j] = 0;
    for (int i = j + 1; i < n; i++) {
      bj[i] = d[i - j - 1];
      if (i < i1) buf[(size_t) (i - i0) * n + j] = bj[i];
    }
  }
}


//------------------------------------------------------
// Binary search for the precision beta of row i's
// Gaussian, as in scikit-learn's
// _binary_search_perplexity(). `sq` holds m squared
// distances; `self` is the entry to leave out (or -1).
// Leaves the normalized probabilities in p, and the
// beta, normalizer and shift they were computed with in
// *beta_out, *sum_out and *shift_out:
//
//   p[j] = exp(-(sq[j] - shift) * beta) / sum
//
// Returns the sum of p.
//
// scikit-learn's search starts at beta = 1 and has 100
// steps, so it can't reach a precision far from 1, e.g.
// for rows whose nearest neighbors are many orders of
// magnitude closer than the farthest sample (distances are
// scaled so that the largest is about 1). When it doesn't
// converge, the search starts again from a better guess,
// with shifted squared distances (so the largest term is
// 1, and the sum can't underflow), bisecting beta
// geometrically for as many steps as the range of doubles
// needs.
//------------------------------------------------------
static double entropy_at(
    const double *sq, int m, int self, double shift, double beta,
    double *p, double *sum_out ) {

  double sum_p = 0;
  for (int j = 0; j < m; j++) {
    if (j == self) continue;
    p[j]   = exp(-(sq[j] - shift) * beta);
    sum_p += p[j];
  }
  if (sum_p == 0.0) sum_p = TS_PERP_EPS;

  double sum_dp = 0;
  for (int j = 0; j < m; j++) {
    p[j]   /= sum_p;
    sum_dp += (sq[j] - shift) * p[j];
  }

  *sum_out = sum_p;
  return log(sum_p) + beta * sum_dp;
}


static double search_row(
    const double *sq, int m, int self, double target,
    double *p, double *beta_out, double *sum_out, double *shift_out ) {

  double beta_min = -INFINITY, beta_max = INFINITY, beta = 1.0;
  double sum_p = 0, used_beta = beta;
  int    done = 0;

  if (self >= 0) p[self] = 0;

  for (int step = 0; step < TS_PERP_STEPS && !done; step++) {

    double diff = entropy_at(sq, m, self, 0, beta, p, &sum_p) - target;
    used_beta   = beta;

    if (fabs(diff) <= TS_PERP_TOL) { done = 1; break; }

    if (diff > 0.0) {
      beta_min = beta;
      beta     = (beta_max == INFINITY) ? beta * 2.0 : (beta + beta_max) / 2.0;
    }
    else {
      beta_max = beta;
      beta     = (beta_min == -INFINITY) ? beta / 2.0 : (beta + beta_min) / 2.0;
    }
  }

  double shift = 0;

  if (!done) {

    // Start from 1 / (mean squared distance beyond the
    // nearest), with the nearest's squared distance as
    // the shift.
    shift = INFINITY;
    for (int j = 0; j < m; j++)
      if (j != self && sq[j] < shift) shift = sq[j];

    double mean = 0;
    int    cnt  = 0;
    for (int j = 0; j < m; j++)
      if (j != self && sq[j] > shift) { mean += sq[j] - shift; cnt++; }

    beta     = cnt ? cnt / mean : 1.0;
    beta_min = 0;
    beta_max = INFINITY;

    // Doubling or halving covers the doubles' range in
    // about 2100 steps; then each bisection halves
    // log(beta_max / beta_min).
    for (int step = 0; step < 2400; step++) {

      double diff = entropy_at(sq, m, self, shift, beta, p, &sum_p) - target;
      used_beta   = beta;

      if (fabs(diff) <= TS_PERP_TOL) break;

      double next;
      if (diff > 0.0) {
        beta_min = beta;
        next = (beta_max == INFINITY) ? beta * 2.0 : sqrt(beta) * sqrt(beta_max);
      }
      else {
        beta_max = beta;
        next = (beta_min == 0) ? beta / 2.0 : sqrt(beta) * sqrt(beta_min);
      }

      // The target can't be reached (e.g. it is below the
      // entropy of several equally near neighbors).
      if (!(next > 0 && next <= DBL_MAX)) break;
      beta = next;
    }
  }

  *beta_out  = used_beta;
  *sum_out   = sum_p;
  *shift_out = shift;

  double total = 0;
  for (int j = 0; j < m; j++) total += p[j];
  return total;
}


// A squared distance, after scaling the distance by
// pp_scale (a power of 2, so exactly), and optionally
// rounding to single precision as scikit-learn does.
static inline double sqd(double d) {
  double s = d * pp_scale;
  s = s * s;
  return pp_f32 ? (double) (float) s : s;
}


//------------------------------------------------------
// Dense: rows of D a block at a time, each row searched
// over all other samples.
//------------------------------------------------------
static void *p_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pp_n;

  double *buf = pp_rowbuf + (size_t) thread_i * TS_ROWS * n;
  double *sq  = pp_sqbuf  + (size_t) thread_i * n;
  double *pr  = pp_pbuf   + (size_t) thread_i * n;

  for (int task = pp_task0 + thread_i; task < pp_task1; task += n_threads) {

    int i0 = task * TS_ROWS;
    int nr = (i0 + TS_ROWS < n) ? TS_ROWS : n - i0;

    gather_rows(i0, nr, buf);

    for (int r = 0; r < nr; r++) {
      int           i   = i0 + r;
      const double *row = buf + (size_t) r * n;
      for (int j = 0; j < n; j++) sq[j] = sqd(row[j]);
      pp_rowsum[i] = search_row(
        sq, n, i, pp_target, pr, pp_beta + i, pp_sum + i, pp_shift + i );
    }
  }

  return NULL;
}


//------------------------------------------------------
// Sparse: each row's K nearest neighbors (from
// knn_from_dist(), by distance and then index) put in
// ascending index order, as scikit-learn's sorted CSR
// graph has them, then searched. On entry, pp_nbr and
// pp_cond hold the neighbors and their distances; on exit,
// the neighbors (sorted) and their p_{j|i}.
//------------------------------------------------------
typedef struct { int j; double d; } ts_nbr_t;

static int nbr_by_index(const void *a, const void *b) {
  int ja = ((const ts_nbr_t *) a)->j, jb = ((const ts_nbr_t *) b)->j;
  return (ja > jb) - (ja < jb);
}

static ts_nbr_t *pp_pairs;  // Per thread: K.

static void *nn_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pp_n;
  int K         = pp_K;

  ts_nbr_t *pr = pp_pairs + (size_t) thread_i * K;
  double   *sq = pp_sqbuf + (size_t) thread_i * K;

  for (int task = pp_task0 + thread_i; task < pp_task1; task += n_threads) {

    int i1 = (task + 1) * TS_ROWS;
    if (i1 > n) i1 = n;

    for (int i = task * TS_ROWS; i < i1; i++) {

      int    *nb = pp_nbr  + (size_t) i * K;
      double *pc = pp_cond + (size_t) i * K;

      for (int t = 0; t < K; t++) { pr[t].j = nb[t]; pr[t].d = pc[t]; }
      qsort(pr, (size_t) K, sizeof(ts_nbr_t), nbr_by_index);
      for (int t = 0; t < K; t++) { nb[t] = pr[t].j; sq[t] = sqd(pr[t].d); }

      pp_rowsum[i] = search_row(
        sq, K, -1, pp_target, pc, pp_beta + i, pp_sum + i, pp_shift + i );
    }
  }

  return NULL;
}


//------------------------------------------------------
// Dense P: packed lower triangle, each column filled by
// one thread.
//------------------------------------------------------
static double *pp_out;
static double  pp_total;

static void *dense_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pp_n;

  for (int j = pp_task0 + thread_i; j < pp_task1; j += n_threads) {
    size_t        off = dist_col_start(j, n);  // Element (j + 1, j).
    const double *d   = pp_dist + off;
    double       *out = pp_out  + off;
    double        bj  = pp_beta[j], sj = pp_sum[j], hj = pp_shift[j];
    for (int t = 0; t < n - j - 1; t++) {
      int    i = j + 1 + t;
      double s = sqd(d[t]);
      double v = exp(-(s - pp_shift[i]) * pp_beta[i]) / pp_sum[i] + exp(-(s - hj) * bj) / sj;
      v /= pp_total;
      out[t] = (v > TS_MACHINE_EPSILON) ? v : TS_MACHINE_EPSILON;
    }
  }

  return NULL;
}


//------------------------------------------------------
// Sparse P: row i of P_cond + P_cond', merging row i of
// P_cond with row i of its transpose (both ascending).
//------------------------------------------------------
static const int    *sp_tptr;  // n + 1
static const int    *sp_tidx;  // n * K
static const double *sp_tval;  // n * K
static int          *sp_ptr;   // n + 1
static int          *sp_col;
static double       *sp_val;
static int          *sp_cnt;   // n
static int           sp_fill;  // 0 = count, 1 = fill.

static void *sparse_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pp_n;
  int K         = pp_K;

  for (int i = thread_i; i < n; i += n_threads) {

    const int    *a  = pp_nbr  + (size_t) i * K;
    const double *av = pp_cond + (size_t) i * K;
    const int    *b  = sp_tidx + sp_tptr[i];
    const double *bv = sp_tval + sp_tptr[i];
    int na = K, nb = sp_tptr[i + 1] - sp_tptr[i];
    int x = 0, y = 0, c = 0;

    int    *col = sp_fill ? sp_col + sp_ptr[i] : NULL;
    double *val = sp_fill ? sp_val + sp_ptr[i] : NULL;

    while (x < na || y < nb) {
      int    j;
      double v;
      if (y >= nb || (x < na && a[x] < b[y])) { j = a[x]; v = av[x]; x++; }
      else if (x >= na || b[y] < a[x])        { j = b[y]; v = bv[y]; y++; }
      else { j = a[x]; v = av[x] + bv[y]; x++; y++; }
      if (sp_fill) { col[c] = j; val[c] = v; }
      c++;
    }

    if (!sp_fill) sp_cnt[i] = c;
  }

  return NULL;
}


static void *scale_worker(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  for (int i = thread_i; i < pp_n; i += n_threads)
    for (int e = sp_ptr[i]; e < sp_ptr[i + 1]; e++)
      sp_val[e] /= pp_total;
  return NULL;
}


//======================================================
// R interface. Joint probabilities P from a dist vector
// for n >= 2 samples. `sexp_K` is the number of nearest
// neighbors (1 to n - 1), or 0 for all samples (dense).
// `sexp_f32` rounds squared distances to single precision
// as scikit-learn does; tsne() passes FALSE, and only the
// tests, which compare P with scikit-learn's, pass TRUE.
// Each distance is first
// multiplied by `sexp_scale`, a power of 2, so that the
// largest is about 1. Returns the packed lower triangle of
// P (dense), or list(ptr, col, val), a CSR matrix with
// 0-based indices (sparse).
//======================================================

// Runs `func` over tasks 0 to n_tasks - 1 (blocks of rows
// or columns) in batches of about 50 ms, checking for a
// user interrupt between batches.
static void run_batches(pthread_func_t func, int n_threads, int n_tasks, double task_work) {
  int batch = (int) (5e7 / (task_work > 1 ? task_work : 1)) * n_threads;
  if (batch < n_threads) batch = n_threads;
  for (int t0 = 0; t0 < n_tasks; t0 += batch) {
    pp_task0 = t0;
    pp_task1 = (n_tasks - t0 > batch) ? t0 + batch : n_tasks;
    run_tasks(func, n_threads);
    R_CheckUserInterrupt();
  }
}

SEXP C_tsne_p(SEXP sexp_dist, SEXP sexp_n, SEXP sexp_perplexity,
              SEXP sexp_K, SEXP sexp_f32, SEXP sexp_scale, SEXP sexp_n_threads) {

  int    n          = asInteger(sexp_n);
  int    K          = asInteger(sexp_K);
  double perplexity = asReal(sexp_perplexity);
  int    n_threads  = asInteger(sexp_n_threads);

  if (n < 2 || K < 0 || K > n - 1 || !(perplexity > 0) || n_threads < 1 ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_tsne_p()."); // # nocov

  if (K > 0 && (double) n * K * 2 > INT_MAX)
    error("Too many nearest neighbors: reduce `perplexity`."); // # nocov

  pp_n      = n;
  pp_K      = K;
  pp_f32    = asLogical(sexp_f32) == 1;
  pp_scale  = asReal(sexp_scale);
  pp_dist   = REAL(sexp_dist);
  pp_target = log((double) (float) perplexity);

  // Work per row: reading n distances, then the binary
  // search over K (or n) of them.
  int    n_tasks  = (n + TS_ROWS - 1) / TS_ROWS;
  double row_work = (K ? 10.0 * K : (double) n) + 150.0 * (K ? K : n);
  int    all_threads = n_threads;
  n_threads  = use_threads(n_threads, n_tasks, row_work * n);
  pp_threads = n_threads;

  pp_beta   = (double *) R_alloc(n, sizeof(double));
  pp_sum    = (double *) R_alloc(n, sizeof(double));
  pp_rowsum = (double *) R_alloc(n, sizeof(double));
  pp_shift  = (double *) R_alloc(n, sizeof(double));

  if (K == 0) {
    pp_rowbuf = (double *) R_alloc((size_t) n_threads * TS_ROWS * n, sizeof(double));
    pp_sqbuf  = (double *) R_alloc((size_t) n_threads * n, sizeof(double));
    pp_pbuf   = (double *) R_alloc((size_t) n_threads * n, sizeof(double));
    run_batches(p_worker, n_threads, n_tasks, row_work * TS_ROWS);
  }
  else {
    pp_nbr   = (int *)      R_alloc((size_t) n * K, sizeof(int));
    pp_cond  = (double *)   R_alloc((size_t) n * K, sizeof(double));
    pp_sqbuf = (double *)   R_alloc((size_t) n_threads * K, sizeof(double));
    pp_pairs = (ts_nbr_t *) R_alloc((size_t) n_threads * K, sizeof(ts_nbr_t));

    // Each point's K nearest: about one comparison per
    // distance, read twice.
    knn_from_dist(
      pp_dist, n, K, 0, pp_nbr, pp_cond,
      use_threads(all_threads, n, 2.0 * n * n) );
    R_CheckUserInterrupt();

    run_batches(nn_worker, n_threads, n_tasks, row_work * TS_ROWS);
  }


  // Sum of P_cond + P_cond', in row order.
  double total = 0;
  for (int i = 0; i < n; i++) total += pp_rowsum[i];
  total *= 2;
  if (total < TS_MACHINE_EPSILON) total = TS_MACHINE_EPSILON; // # nocov
  pp_total = total;

  if (K == 0) {
    int col_threads = use_threads(n_threads, n - 1, (double) n * n * 20);
    SEXP sexp_p = PROTECT(allocVector(REALSXP, XLENGTH(sexp_dist)));
    pp_out = REAL(sexp_p);
    run_batches(dense_worker, col_threads, n - 1, (double) n * 20);
    UNPROTECT(1);
    return sexp_p;
  }


  //--------------------------------------------------
  // Transpose of P_cond: counting sort by column. Rows
  // are visited in order, so each transposed row is
  // ascending.
  //--------------------------------------------------
  size_t nnz  = (size_t) n * K;
  int   *tptr = (int *)    R_alloc((size_t) n + 1, sizeof(int));
  int   *tidx = (int *)    R_alloc(nnz, sizeof(int));
  double *tval = (double *) R_alloc(nnz, sizeof(double));
  int   *fill = (int *)    R_alloc(n, sizeof(int));

  memset(tptr, 0, ((size_t) n + 1) * sizeof(int));
  for (size_t e = 0; e < nnz; e++) tptr[pp_nbr[e] + 1]++;
  for (int i = 0; i < n; i++) tptr[i + 1] += tptr[i];
  memcpy(fill, tptr, (size_t) n * sizeof(int));
  for (int i = 0; i < n; i++)
    for (int t = 0; t < K; t++) {
      size_t e = (size_t) i * K + t;
      int    p = fill[pp_nbr[e]]++;
      tidx[p]  = i;
      tval[p]  = pp_cond[e];
    }

  sp_tptr = tptr;
  sp_tidx = tidx;
  sp_tval = tval;
  sp_cnt  = (int *) R_alloc(n, sizeof(int));

  int sym_threads = use_threads(n_threads, n, (double) nnz * 10);

  sp_fill = 0;
  run_tasks(sparse_worker, sym_threads);

  SEXP sexp_ptr = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
  sp_ptr    = INTEGER(sexp_ptr);
  sp_ptr[0] = 0;
  for (int i = 0; i < n; i++) sp_ptr[i + 1] = sp_ptr[i] + sp_cnt[i];

  SEXP sexp_col = PROTECT(allocVector(INTSXP,  sp_ptr[n]));
  SEXP sexp_val = PROTECT(allocVector(REALSXP, sp_ptr[n]));
  sp_col = INTEGER(sexp_col);
  sp_val = REAL(sexp_val);

  sp_fill = 1;
  run_tasks(sparse_worker, sym_threads);

  // Normalize by the sum of all entries, in order.
  total = 0;
  for (int e = 0; e < sp_ptr[n]; e++) total += sp_val[e];
  if (total < TS_MACHINE_EPSILON) total = TS_MACHINE_EPSILON; // # nocov
  pp_total = total;
  run_tasks(scale_worker, sym_threads);

  const char *names[] = { "ptr", "col", "val", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_ptr);
  SET_VECTOR_ELT(sexp_result, 1, sexp_col);
  SET_VECTOR_ELT(sexp_result, 2, sexp_val);

  UNPROTECT(4);
  return sexp_result;
}



//======================================================
// Gradient
//======================================================

static int           gr_n, gr_k;
static int           gr_exact;    // 1 = exact, 0 = Barnes-Hut.
static const double *gr_y;        // n x k, row-major.
static const double *gr_pv;       // Exact: packed lower triangle of P.
static const int    *gr_ptr;      // Barnes-Hut: CSR P.
static const int    *gr_col;
static const double *gr_val;
static double        gr_theta2;
static double        gr_exag;
static int           gr_kl;       // Compute the KL divergence?
static double        gr_zsum;     // Z
static double       *gr_attr;     // n x k: sum_j p_ij w_ij (y_i - y_j)
static double       *gr_rep;      // n x k: sum_j w_ij^2 (y_i - y_j)
static double       *gr_z;        // n: sum_j w_ij
static double       *gr_klv;      // n: per-sample KL terms (see kl_term())
static double       *gr_ps;       // n: per-sample sums of a * p_ij
static int           gr_threads;


//------------------------------------------------------
// KL divergence terms. With pe = a * p_ij (a the
// exaggeration) and q_ij = w_ij / Z,
//
//   KL = sum pe log(pe / q) = sum pe log(pe / w) + log(Z) sum pe
//
// so each sample's terms can be summed before Z is known.
//------------------------------------------------------
static inline double kl_term(double pe, double w) {
  return (pe > 0) ? pe * log(pe / w) : 0;
}


//------------------------------------------------------
// Exact: rows i0 to i1 - 1 of attr, rep and z (and, if
// `kl`, the KL terms and sum of pe), over j in ascending
// order. Pairs are read from P's packed lower triangle as
// contiguous runs; pairs within the block are computed
// once, for both rows.
//------------------------------------------------------
TS_INLINE void exact_rows(const int k, const int kl, int i0, int i1) {

  int           n  = gr_n;
  const double *Y  = gr_y;
  const double *P  = gr_pv;
  double       *A  = gr_attr;
  double       *B  = gr_rep;
  double       *Zv = gr_z;
  double       *Kv = gr_klv;
  double       *Sv = gr_ps;
  double        ex = gr_exag;

  memset(A + (size_t) i0 * k, 0, (size_t) (i1 - i0) * k * sizeof(double));
  memset(B + (size_t) i0 * k, 0, (size_t) (i1 - i0) * k * sizeof(double));
  memset(Zv + i0, 0, (size_t) (i1 - i0) * sizeof(double));
  memset(Kv + i0, 0, (size_t) (i1 - i0) * sizeof(double));
  memset(Sv + i0, 0, (size_t) (i1 - i0) * sizeof(double));

  for (int j = 0; j < i0; j++) {
    const double *p  = P + dist_col_start(j, n) + (i0 - j - 1);
    const double *yj = Y + (size_t) j * k;
    for (int i = i0; i < i1; i++) {
      const double *yi = Y + (size_t) i * k;
      double d2 = 0;
      for (int d = 0; d < k; d++) { double t = yi[d] - yj[d]; d2 += t * t; }
      double w = 1 / (1 + d2), a = p[i - i0] * w, b = w * w;
      double *ai = A + (size_t) i * k, *bi = B + (size_t) i * k;
      for (int d = 0; d < k; d++) {
        double t = yi[d] - yj[d];
        ai[d] += a * t;
        bi[d] += b * t;
      }
      Zv[i] += w;
      if (kl) {
        double pe = ex * p[i - i0];
        Kv[i] += kl_term(pe, w);
        Sv[i] += pe;
      }
    }
  }

  for (int j = i0; j < i1; j++) {
    const double *p  = P + dist_col_start(j, n);  // p[t]: row j + 1 + t
    const double *yj = Y + (size_t) j * k;
    double *aj = A + (size_t) j * k, *bj = B + (size_t) j * k;
    for (int i = j + 1; i < n; i++) {
      const double *yi = Y + (size_t) i * k;
      double d2 = 0;
      for (int d = 0; d < k; d++) { double t = yi[d] - yj[d]; d2 += t * t; }
      double pij = p[i - j - 1];
      double w = 1 / (1 + d2), a = pij * w, b = w * w;
      for (int d = 0; d < k; d++) {
        double t = yj[d] - yi[d];
        aj[d] += a * t;
        bj[d] += b * t;
      }
      Zv[j] += w;
      double term = 0, pe = 0;
      if (kl) {
        pe   = ex * pij;
        term = kl_term(pe, w);
        Kv[j] += term;
        Sv[j] += pe;
      }
      if (i < i1) {
        double *ai = A + (size_t) i * k, *bi = B + (size_t) i * k;
        for (int d = 0; d < k; d++) {
          double t = yi[d] - yj[d];
          ai[d] += a * t;
          bi[d] += b * t;
        }
        Zv[i] += w;
        if (kl) { Kv[i] += term; Sv[i] += pe; }
      }
    }
  }
}


static void *exact_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = gr_n;
  int n_tasks   = (n + TS_BLOCK - 1) / TS_BLOCK;
  int kl        = gr_kl;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int i0 = task * TS_BLOCK;
    int i1 = (i0 + TS_BLOCK < n) ? i0 + TS_BLOCK : n;
    switch (gr_k * 2 + kl) {
      case 2:  exact_rows(1, 0, i0, i1); break;
      case 3:  exact_rows(1, 1, i0, i1); break;
      case 4:  exact_rows(2, 0, i0, i1); break;
      case 5:  exact_rows(2, 1, i0, i1); break;
      case 6:  exact_rows(3, 0, i0, i1); break;
      case 7:  exact_rows(3, 1, i0, i1); break;
      default: exact_rows(gr_k, kl, i0, i1);
    }
  }

  return NULL;
}



//------------------------------------------------------
// Barnes-Hut tree: a quadtree (k = 2), octree (k = 3) or
// binary tree (k = 1), as scikit-learn's _QuadTree: the
// root box is the points' bounding box, slightly enlarged;
// cells split at their centers; and a child's index has
// one bit per axis (axis 0 most significant), set when the
// point is >= the center. Unlike _QuadTree, a cell whose
// points all fall in one child is not kept (it summarizes
// whenever that child would, so this changes nothing);
// identical points share a leaf, which acts on other
// points as one point of their combined weight; and points
// still together after BH_MAX_DEPTH halvings share a leaf
// and interact pairwise.
//------------------------------------------------------
typedef struct {
  double com[3];  // Center of mass.
  double w2;      // Squared width of the cell (largest axis).
  int    count;   // Points in the cell.
  int    begin;   // Their first position in bh_perm.
  int    child;   // First child node, or -1 for a leaf.
  int    nchild;  // Number of (nonempty) children.
  int    same;    // Leaf whose points are all identical?
} bh_node_t;

static bh_node_t     *bh_nodes;
static int            bh_count;    // Nodes in use.
static int           *bh_perm;     // Points, grouped by cell.
static int           *bh_inv;      // Position of each point in bh_perm.
static int           *bh_tmp;      // n
static unsigned char *bh_code;     // n
static int            bh_depth;    // Deepest node.
static int           *bh_stack;    // Per thread: BH_STACK.

#define BH_STACK (8 * (BH_MAX_DEPTH + 2))


static void bh_build(const double *Y, int D, int node, double *lo, double *hi,
                     int begin, int count, int depth, int level, double *sum) {

  bh_node_t *nd = bh_nodes + node;
  int *perm = bh_perm + begin;

  if (level > bh_depth) bh_depth = level;

  nd->count = count;
  nd->begin = begin;
  nd->child = -1;
  nd->nchild = 0;
  nd->same   = 0;

  while (1) {

    // Leaf: all points identical (or just one), or too deep.
    const double *y0 = Y + (size_t) perm[0] * D;
    int same = 1;
    for (int p = 1; p < count && same; p++) {
      const double *yp = Y + (size_t) perm[p] * D;
      for (int d = 0; d < D; d++)
        if (yp[d] != y0[d]) { same = 0; break; }
    }
    int leaf = same || depth >= BH_MAX_DEPTH;

    double w2 = 0;
    for (int d = 0; d < D; d++) {
      double w = hi[d] - lo[d];
      if (w * w > w2) w2 = w * w;
    }
    nd->w2 = w2;

    if (same) {
      nd->same = 1;
      for (int d = 0; d < D; d++) { nd->com[d] = y0[d]; sum[d] = y0[d] * count; }
      return;
    }

    if (leaf) {
      for (int d = 0; d < D; d++) sum[d] = 0;
      for (int p = 0; p < count; p++)
        for (int d = 0; d < D; d++)
          sum[d] += Y[(size_t) perm[p] * D + d];
      break;
    }

    // Children of each point.
    double ctr[3];
    for (int d = 0; d < D; d++) ctr[d] = (lo[d] + hi[d]) / 2;

    int counts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int p = 0; p < count; p++) {
      const double *yp = Y + (size_t) perm[p] * D;
      int c = 0;
      for (int d = 0; d < D; d++) c = 2 * c + (yp[d] >= ctr[d]);
      bh_code[begin + p] = (unsigned char) c;
      counts[c]++;
    }

    int nc = 1 << D, nonempty = 0, last = 0;
    for (int c = 0; c < nc; c++)
      if (counts[c]) { nonempty++; last = c; }

    // All in one child: that child is this cell.
    if (nonempty == 1) {
      for (int d = 0; d < D; d++) {
        if ((last >> (D - 1 - d)) & 1) lo[d] = ctr[d];
        else                           hi[d] = ctr[d];
      }
      depth++;
      continue;
    }

    // Stable counting sort of the points by child.
    int start[8], pos[8];
    start[0] = 0;
    for (int c = 1; c < nc; c++) start[c] = start[c - 1] + counts[c - 1];
    memcpy(pos, start, sizeof(pos));
    for (int p = 0; p < count; p++)
      bh_tmp[pos[bh_code[begin + p]]++] = perm[p];
    memcpy(perm, bh_tmp, (size_t) count * sizeof(int));

    int first = bh_count;
    bh_count += nonempty;
    nd->child  = first;
    nd->nchild = nonempty;

    for (int d = 0; d < D; d++) sum[d] = 0;

    int ci = first;
    for (int c = 0; c < nc; c++) {
      if (!counts[c]) continue;
      double clo[3], chi[3], csum[3];
      for (int d = 0; d < D; d++) {
        if ((c >> (D - 1 - d)) & 1) { clo[d] = ctr[d]; chi[d] = hi[d];  }
        else                        { clo[d] = lo[d];  chi[d] = ctr[d]; }
      }
      bh_build(Y, D, ci, clo, chi, begin + start[c], counts[c],
               depth + 1, level + 1, csum);
      for (int d = 0; d < D; d++) sum[d] += csum[d];
      ci++;
    }

    nd = bh_nodes + node; // (Unchanged: nodes are preallocated.)
    break;
  }

  for (int d = 0; d < D; d++) nd->com[d] = sum[d] / count;
}


static void bh_make_tree(const double *Y, int n, int D) {

  double lo[3], hi[3], sum[3];

  for (int d = 0; d < D; d++) {
    lo[d] = hi[d] = Y[d];
    for (int i = 1; i < n; i++) {
      double v = Y[(size_t) i * D + d];
      if (v < lo[d]) lo[d] = v;
      if (v > hi[d]) hi[d] = v;
    }
    // As _QuadTree: every point strictly inside the box.
    double m1 = hi[d] * (1 + 1e-3 * ((hi[d] > 0) - (hi[d] < 0)));
    double m2 = hi[d] + 1e-3;
    hi[d] = (m1 > m2) ? m1 : m2;
  }

  for (int i = 0; i < n; i++) bh_perm[i] = i;

  bh_count = 1;
  bh_depth = 0;
  bh_build(Y, D, 0, lo, hi, 0, n, 0, 0, sum);

  for (int p = 0; p < n; p++) bh_inv[bh_perm[p]] = p;
}


//------------------------------------------------------
// Repulsive forces on point i: rep_i and z_i, walking the
// tree depth-first in child order.
//------------------------------------------------------
TS_INLINE void bh_point(const int D, int i, int *stack) {

  const double *yi = gr_y + (size_t) i * D;
  int    pi = bh_inv[i];
  double rep[3] = {0, 0, 0}, z = 0;
  double theta2 = gr_theta2;
  int    top = 0;

  stack[top++] = 0;

  while (top > 0) {

    const bh_node_t *nd = bh_nodes + stack[--top];

    double diff[3], d2 = 0;
    for (int d = 0; d < D; d++) {
      diff[d] = yi[d] - nd->com[d];
      d2     += diff[d] * diff[d];
    }

    int holds = pi >= nd->begin && pi < nd->begin + nd->count;

    if (nd->child < 0) {

      if (nd->same) {
        // Identical points: exact, as one weighted point.
        int m = nd->count - holds;
        double w = 1 / (1 + d2), b = m * w * w;
        z += m * w;
        for (int d = 0; d < D; d++) rep[d] += b * diff[d];
      }
      else {
        // Points too close to separate, pairwise.
        for (int p = nd->begin; p < nd->begin + nd->count; p++) {
          int j = bh_perm[p];
          if (j == i) continue;
          const double *yj = gr_y + (size_t) j * D;
          double t[3], e2 = 0;
          for (int d = 0; d < D; d++) { t[d] = yi[d] - yj[d]; e2 += t[d] * t[d]; }
          double w = 1 / (1 + e2), b = w * w;
          z += w;
          for (int d = 0; d < D; d++) rep[d] += b * t[d];
        }
      }
      continue;
    }

    // Summarize a cell that doesn't hold point i and looks
    // small enough from it.
    if (!holds && nd->w2 / d2 < theta2) {
      double w = 1 / (1 + d2), b = nd->count * w * w;
      z += nd->count * w;
      for (int d = 0; d < D; d++) rep[d] += b * diff[d];
      continue;
    }

    for (int c = nd->nchild - 1; c >= 0; c--)
      stack[top++] = nd->child + c;
  }

  for (int d = 0; d < D; d++) gr_rep[(size_t) i * D + d] = rep[d];
  gr_z[i] = z;
}


//------------------------------------------------------
// Attractive forces on point i from the sparse P, and its
// KL divergence terms (see kl_term()).
//------------------------------------------------------
TS_INLINE void bh_attract(const int k, int i) {

  const double *yi = gr_y + (size_t) i * k;
  double       *ai = gr_attr + (size_t) i * k;
  double        kl = 0, ps = 0;

  for (int d = 0; d < k; d++) ai[d] = 0;

  for (int e = gr_ptr[i]; e < gr_ptr[i + 1]; e++) {
    const double *yj = gr_y + (size_t) gr_col[e] * k;
    double d2 = 0;
    for (int d = 0; d < k; d++) { double t = yi[d] - yj[d]; d2 += t * t; }
    double w = 1 / (1 + d2), a = gr_val[e] * w;
    for (int d = 0; d < k; d++) ai[d] += a * (yi[d] - yj[d]);
    if (gr_kl) {
      double pe = gr_exag * gr_val[e];
      kl += kl_term(pe, w);
      ps += pe;
    }
  }

  gr_klv[i] = kl;
  gr_ps[i]  = ps;
}


//------------------------------------------------------
// Both forces on each point, in one parallel pass (the
// attractive force doesn't need Z). Points are taken in
// tree order, for locality.
//------------------------------------------------------
static void *bh_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = gr_n;
  int n_tasks   = (n + TS_BLOCK - 1) / TS_BLOCK;
  int *stack    = bh_stack + (size_t) thread_i * BH_STACK;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int p1 = (task + 1) * TS_BLOCK;
    if (p1 > n) p1 = n;
    for (int p = task * TS_BLOCK; p < p1; p++) {
      int i = bh_perm[p];
      switch (gr_k) {
        case 1:  bh_point(1, i, stack); bh_attract(1, i); break;
        case 2:  bh_point(2, i, stack); bh_attract(2, i); break;
        default: bh_point(3, i, stack); bh_attract(3, i);
      }
    }
  }

  return NULL;
}


//------------------------------------------------------
// The gradient at gr_y, into `grad` (n x k, row-major).
// Returns the KL divergence if gr_kl, else NaN. One
// parallel pass either way.
//------------------------------------------------------
static double gradient(double *grad) {

  int n = gr_n, k = gr_k;

  if (gr_exact) {
    run_tasks(exact_worker, gr_threads);
  }
  else {
    bh_make_tree(gr_y, n, k);
    run_tasks(bh_worker, gr_threads);
  }

  double Z = 0;
  for (int i = 0; i < n; i++) Z += gr_z[i];
  if (!gr_exact && Z < DBL_EPSILON) Z = DBL_EPSILON; // # nocov
  gr_zsum = Z;

  for (size_t e = 0; e < (size_t) n * k; e++)
    grad[e] = 4 * (gr_exag * gr_attr[e] - gr_rep[e] / Z);

  if (!gr_kl) return NAN;

  double kl = 0, ps = 0;
  for (int i = 0; i < n; i++) { kl += gr_klv[i]; ps += gr_ps[i]; }
  return kl + ps * log(Z);
}


//------------------------------------------------------
// Set up the gradient's inputs and work space. `sexp_p`
// is a dense P (numeric vector) for the exact gradient,
// or a CSR list for Barnes-Hut.
//------------------------------------------------------
static void gradient_setup(SEXP sexp_p, int n, int k, double theta, int n_threads) {

  gr_n     = n;
  gr_k     = k;
  gr_exact = TYPEOF(sexp_p) == REALSXP;

  if (gr_exact) {
    if (XLENGTH(sexp_p) != (R_xlen_t) n * (n - 1) / 2)
      error("Invalid P for C_tsne_*()."); // # nocov
    gr_pv = REAL(sexp_p);
  }
  else {
    if (TYPEOF(sexp_p) != VECSXP || XLENGTH(sexp_p) != 3 || k > 3 ||
        XLENGTH(VECTOR_ELT(sexp_p, 0)) != (R_xlen_t) n + 1)
      error("Invalid P for C_tsne_*()."); // # nocov
    gr_ptr = INTEGER(VECTOR_ELT(sexp_p, 0));
    gr_col = INTEGER(VECTOR_ELT(sexp_p, 1));
    gr_val = REAL(VECTOR_ELT(sexp_p, 2));
  }

  gr_theta2 = theta * theta;

  double work = gr_exact ? (double) n * n * k : (double) n * 1500 * k;
  gr_threads  = use_threads(n_threads, (n + TS_BLOCK - 1) / TS_BLOCK, work);

  gr_attr = (double *) R_alloc((size_t) n * k, sizeof(double));
  gr_rep  = (double *) R_alloc((size_t) n * k, sizeof(double));
  gr_z    = (double *) R_alloc(n, sizeof(double));
  gr_klv  = (double *) R_alloc(n, sizeof(double));
  gr_ps   = (double *) R_alloc(n, sizeof(double));

  if (!gr_exact) {
    bh_nodes = (bh_node_t *)     R_alloc((size_t) 2 * n, sizeof(bh_node_t));
    bh_perm  = (int *)           R_alloc(n, sizeof(int));
    bh_inv   = (int *)           R_alloc(n, sizeof(int));
    bh_tmp   = (int *)           R_alloc(n, sizeof(int));
    bh_code  = (unsigned char *) R_alloc(n, sizeof(unsigned char));
    bh_stack = (int *)           R_alloc((size_t) gr_threads * BH_STACK, sizeof(int));
  }
}


// R's column-major n x k matrix to row-major, and back.
static void to_rows(const double *x, int n, int k, double *y) {
  for (int i = 0; i < n; i++)
    for (int d = 0; d < k; d++)
      y[(size_t) i * k + d] = x[(size_t) d * n + i];
}

static void to_cols(const double *y, int n, int k, double *x) {
  for (int i = 0; i < n; i++)
    for (int d = 0; d < k; d++)
      x[(size_t) d * n + i] = y[(size_t) i * k + d];
}





//------------------------------------------------------
// Starting layouts (row-major y, n x k).
//------------------------------------------------------

// PCG stream for random starts.
#define TS_STREAM 20261117u

// Normal deviates with standard deviation 1e-4 from the
// PCG stream for `seed` (Box-Muller), filled in R's
// column-major order.
static void random_start(double *y, int n, int k, int seed) {

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, (uint64_t) seed /* sign-extended, as elsewhere */, TS_STREAM);

  size_t nk = (size_t) n * k;
  for (size_t e = 0; e < nk; e += 2) {
    double u[2];
    for (int t = 0; t < 2; t++) {
      // 53 random bits, as a double in (0, 1).
      uint64_t a = pcg32_random_r(&rng) >> 5;
      uint64_t b = pcg32_random_r(&rng) >> 6;
      u[t] = ((double) a * 67108864.0 + (double) b + 0.5) / 9007199254740992.0;
    }
    double r = sqrt(-2 * log(u[0]));
    for (int t = 0; t < 2 && e + t < nk; t++) {
      size_t c = e + t, i = c % (size_t) n, d = c / (size_t) n;
      y[i * k + d] = r * (t ? sin(2 * M_PI * u[1]) : cos(2 * M_PI * u[1])) * 1e-4;
    }
  }
}

// Scales y so that its first column's population standard
// deviation is 1e-4, as scikit-learn scales its PCA start,
// in double precision and a fixed order (R's mean() and
// sum() accumulate in long double, whose precision varies
// between platforms).
static void scale_start(double *y, int n, int k) {

  double mean = 0, ss = 0;
  for (int i = 0; i < n; i++) mean += y[(size_t) i * k];
  mean /= n;
  for (int i = 0; i < n; i++) {
    double t = y[(size_t) i * k] - mean;
    ss += t * t;
  }

  double sc = 1e-4 / sqrt(ss / n);
  if (!R_FINITE(sc)) error("Invalid start for C_tsne_run()."); // # nocov
  for (size_t e = 0; e < (size_t) n * k; e++) y[e] *= sc;
}



//======================================================
// R interface. Runs the optimization, with P from
// C_tsne_p() (dense: exact gradient; sparse: Barnes-Hut).
//
//   sexp_init: the n x k starting layout, or NULL for a
//              random one.
//   sexp_iopt: n, k, max_iters, exaggeration iterations,
//              threads, iterations between checks,
//              iterations without progress, whether to
//              center the result, how to start, and the
//              seed. Starts:
//                0 = sexp_init as given;
//                1 = sexp_init scaled so that its first
//                    column has standard deviation 1e-4
//                    (the PCoA start);
//                2 = random, from the seed (see
//                    random_start()).
//   sexp_dopt: theta, early exaggeration, learning rate,
//              min gradient norm, momentum (early phase),
//              momentum (late phase).
//
// Returns list(points, kl, n_iter, diverged, grad): the
// final embedding, its KL divergence, the number of
// iterations run, whether the layout diverged (a
// coordinate became NaN or beyond +-1e100; the run stops
// then), and the gradient at the final embedding (without
// exaggeration). With max_iters = 0, they are the start's
// KL divergence and gradient, which the tests compare with
// scikit-learn's.
//======================================================

// Real layouts are nowhere near this; beyond it, squared
// distances lose all meaning, then overflow.
#define TS_MAX_COORD 1e100

static int layout_ok(const double *y, size_t nk) {
  for (size_t e = 0; e < nk; e++)
    if (!(fabs(y[e]) <= TS_MAX_COORD)) return 0;
  return 1;
}

SEXP C_tsne_run(SEXP sexp_p, SEXP sexp_init, SEXP sexp_iopt, SEXP sexp_dopt) {

  if (TYPEOF(sexp_iopt) != INTSXP || XLENGTH(sexp_iopt) != 10 ||
      TYPEOF(sexp_dopt) != REALSXP || XLENGTH(sexp_dopt) != 6)
    error("Invalid arguments to C_tsne_run()."); // # nocov

  int    *iopt         = INTEGER(sexp_iopt);
  double *dopt         = REAL(sexp_dopt);
  int     n            = iopt[0];
  int     k            = iopt[1];
  int     max_iters     = iopt[2];
  int     exag_iter    = iopt[3];
  int     n_iter_check = iopt[5];
  int     no_progress  = iopt[6];
  int     center       = iopt[7];
  int     start        = iopt[8];
  double  exaggeration = dopt[1];
  double  lr           = dopt[2];
  double  min_grad     = dopt[3];

  if (n < 2 || k < 1 || max_iters < 0 || n_iter_check < 1 || start < 0 || start > 2 ||
      (double) n * k > INT_MAX)
    error("Invalid arguments to C_tsne_run()."); // # nocov

  if (start < 2 && (TYPEOF(sexp_init) != REALSXP || !isMatrix(sexp_init) ||
                    nrows(sexp_init) != n || ncols(sexp_init) != k))
    error("Invalid start for C_tsne_run()."); // # nocov

  gradient_setup(sexp_p, n, k, dopt[0], iopt[4]);

  size_t  nk     = (size_t) n * k;
  double *y      = (double *) R_alloc(nk, sizeof(double));
  double *grad   = (double *) R_alloc(nk, sizeof(double));
  double *update = (double *) R_alloc(nk, sizeof(double));
  double *gains  = (double *) R_alloc(nk, sizeof(double));

  if (start == 2) {
    random_start(y, n, k, iopt[9]);
  }
  else {
    to_rows(REAL(sexp_init), n, k, y);
    if (start == 1) scale_start(y, n, k);
  }
  gr_y = y;

  int it = 0, iters = 0, diverged = 0;

  for (int phase = 0; phase < 2 && !diverged; phase++) {

    int    stop     = phase ? max_iters : ((exag_iter < max_iters) ? exag_iter : max_iters);
    double momentum = dopt[4 + phase];
    double best     = DBL_MAX;
    int    best_it  = it;

    gr_exag = phase ? 1 : exaggeration;

    for (size_t e = 0; e < nk; e++) { update[e] = 0; gains[e] = 1; }

    int i;
    for (i = it; i < stop; i++) {

      int check = (i + 1) % n_iter_check == 0;
      gr_kl     = check || i == stop - 1;

      double error = gradient(grad);

      for (size_t e = 0; e < nk; e++) {
        if (update[e] * grad[e] < 0.0) gains[e] += 0.2;
        else                           gains[e] *= 0.8;
        if (gains[e] < 0.01) gains[e] = 0.01;
        grad[e]  *= gains[e];
        update[e] = momentum * update[e] - lr * grad[e];
        y[e]     += update[e];
      }
      iters++;

      if (check) {
        double norm = 0;
        for (size_t e = 0; e < nk; e++) norm += grad[e] * grad[e];
        norm = sqrt(norm);

        // A runaway layout won't recover.
        if (!R_FINITE(norm) || !layout_ok(y, nk)) { diverged = 1; i++; break; }

        if (error < best)                    { best = error; best_it = i; }
        else if (i - best_it > no_progress)  { i++; break; }
        if (norm <= min_grad)                { i++; break; }
      }

      R_CheckUserInterrupt();
    }

    it = i;
  }

  if (!diverged && !layout_ok(y, nk)) diverged = 1;

  // KL divergence and gradient of the final embedding.
  double kl = NAN;
  if (!diverged) {
    gr_exag = 1;
    gr_kl   = 1;
    kl = gradient(grad);
  }

  // Center (which leaves the KL divergence and gradient as
  // they are), in plain double precision and a fixed
  // order. Signs are fixed in R, by ord_flip().
  if (center) {
    for (int d = 0; d < k; d++) {
      double mean = 0;
      for (int i = 0; i < n; i++) mean += y[(size_t) i * k + d];
      mean /= n;
      for (int i = 0; i < n; i++) y[(size_t) i * k + d] -= mean;
    }
  }

  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_grad   = PROTECT(allocMatrix(REALSXP, n, k));
  to_cols(y, n, k, REAL(sexp_points));
  to_cols(grad, n, k, REAL(sexp_grad));

  const char *names[] = { "points", "kl", "n_iter", "diverged", "grad", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, ScalarReal(kl));
  SET_VECTOR_ELT(sexp_result, 2, ScalarInteger(iters));
  SET_VECTOR_ELT(sexp_result, 3, ScalarLogical(diverged));
  SET_VECTOR_ELT(sexp_result, 4, sexp_grad);

  UNPROTECT(3);
  return sexp_result;
}
