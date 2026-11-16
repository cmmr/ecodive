// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

//======================================================
// t-distributed stochastic neighbour embedding (t-SNE)
// from a packed distance triangle.
//
// Follows van der Maaten's bhtsne as wrapped by the
// Rtsne package (its `is_distance = TRUE` path), with
// the same arithmetic in the same order wherever that
// costs neither speed nor determinism. From the same
// starting coordinates, Barnes-Hut results are
// identical to Rtsne's. In exact mode, the sum that
// normalizes Q runs in parallel, so it differs from
// Rtsne's serial one in rounding, which t-SNE then
// amplifies over the iterations.
//
//   theta == 0: exact t-SNE. P is dense, kept in the
//               same packed layout as the distances.
//   theta  > 0: Barnes-Hut t-SNE. P is sparse, from
//               each point's 3 * perplexity nearest
//               neighbours; repulsion is approximated
//               with a space-partitioning tree.
//
// Points are processed in parallel; each point's work
// (calibration, gradient, cost) is done by one thread,
// which writes only that point's outputs. Sums across
// points run afterwards on the main thread, in point
// order. Results are therefore identical for any
// number of threads.
//
// Coordinates are kept row-major (n x k), as in bhtsne.
//======================================================

#include "ecodive.h"
#include <float.h> // DBL_MIN, DBL_MAX, FLT_MIN

#define TS_CHUNK   16 // Points per work unit.
#define TS_MAXTREE 3  // Highest k for Barnes-Hut.


//------------------------------------------------------
// Shared state for the worker threads.
//------------------------------------------------------
static int           ts_n;      // Points.
static int           ts_k;      // Output dimensions.
static int           ts_K;      // Neighbours per point (Barnes-Hut).
static const double *ts_dist;   // Packed lower triangle, as R's dist.
static double        ts_perp;   // Perplexity.
static double        ts_theta;  // Barnes-Hut accuracy.
static double       *ts_buf;    // n doubles of work space per thread.
static int          *ts_ibuf;   // K ints of work space per thread.

static double *ts_beta;  // n: each point's Gaussian precision ...
static double *ts_zsum;  // n: ... and its kernel's row sum.

// Exact P: packed like ts_dist.
static double *ts_P;

// Sparse P (Barnes-Hut): ts_K columns per point before
// symmetrizing, then rows of varying length.
static size_t *ts_row;   // n + 1
static int    *ts_col;
static double *ts_val;
static int    *ts_mate;  // Position of the reverse edge, or -1.

static const double *ts_Y;    // n x k, row-major.
static double       *ts_sq;   // n: squared norms of Y's rows.
static double       *ts_Q;    // Exact: q_ij, packed like ts_dist.
static double       *ts_out;  // n: per-point sums.
static double       *ts_dC;   // n x k: gradient.
static double       *ts_pos;  // n x k: attractive forces.
static double       *ts_neg;  // n x k: repulsive forces.
static double        ts_sumQ;

static int ts_next; // Next chunk of points to process.

#ifdef HAVE_PTHREAD
  static pthread_mutex_t ts_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif


//------------------------------------------------------
// Claim the next chunk of points: [*lo, *hi). Returns 0
// when none are left. Which thread takes which chunk
// varies, but never what is computed for each point.
//------------------------------------------------------
static int next_chunk(int *lo, int *hi) {

  int c;

  #ifdef HAVE_PTHREAD
    pthread_mutex_lock(&ts_mutex);
  #endif
  c = ts_next++;
  #ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&ts_mutex);
  #endif

  if ((size_t) c * TS_CHUNK >= (size_t) ts_n) return 0;

  *lo = c * TS_CHUNK;
  *hi = (*lo + TS_CHUNK < ts_n) ? *lo + TS_CHUNK : ts_n;
  return 1;
}

static void run_points(pthread_func_t func, int n_threads) {
  ts_next = 0;
  run_parallel(func, n_threads, ts_n);
}


//------------------------------------------------------
// Distances from point i to every point, into out[0..n-1].
//------------------------------------------------------
static inline size_t col_start(int i) {
  return dist_col_start(i, ts_n);
}

static void dist_row(int i, double *out) {

  int n = ts_n;

  for (int m = 0; m < i; m++)
    out[m] = ts_dist[col_start(m) + (size_t) (i - m - 1)];

  out[i] = 0;

  const double *d = ts_dist + col_start(i);
  for (int m = i + 1; m < n; m++)
    out[m] = d[m - i - 1];
}

// The packed index of pair (i, j), i < j.
static inline size_t pair_idx(int i, int j) {
  return col_start(i) + (size_t) (j - i - 1);
}


//------------------------------------------------------
// Binary search for the Gaussian precision (beta) that
// gives a point's conditional distribution the target
// perplexity, as bhtsne does it: start at 1, double or
// halve until bracketed, then bisect, until the entropy
// is within 1e-5 of log(perplexity) or 200 steps.
//
// Exact mode: DD holds squared distances to all n
// points, with the point itself at `self`; its kernel
// value is DBL_MIN, as in bhtsne.
//
// Barnes-Hut: d holds (unsquared) distances to the K
// nearest neighbours; self = -1.
//
// Stores the beta that was last evaluated and its row
// sum; the kernel values are recomputed from them.
//------------------------------------------------------
static void calibrate(
    const double *x, int len, int self, int squared,
    double *beta_out, double *sum_out ) {

  int    found    = 0;
  double beta     = 1.0;
  double min_beta = -DBL_MAX;
  double max_beta =  DBL_MAX;
  double tol      = 1e-5;
  double sum_P    = DBL_MIN;
  double used     = beta;
  int    iter     = 0;

  while (!found && iter < 200) {

    double H = 0.0;
    sum_P    = DBL_MIN;

    if (squared) {
      for (int m = 0; m < len; m++) {
        double p = (m == self) ? DBL_MIN : exp(-beta * x[m]);
        sum_P += p;
        H     += beta * (x[m] * p);
      }
    }
    else {
      for (int m = 0; m < len; m++) {
        double p = exp(-beta * x[m] * x[m]);
        sum_P += p;
        H     += beta * (x[m] * x[m] * p);
      }
    }

    H    = (H / sum_P) + log(sum_P);
    used = beta;

    double Hdiff = H - log(ts_perp);

    if (Hdiff < tol && -Hdiff < tol) {
      found = 1;
    }
    else if (Hdiff > 0) {
      min_beta = beta;
      if (max_beta == DBL_MAX || max_beta == -DBL_MAX) beta *= 2.0;
      else                                             beta  = (beta + max_beta) / 2.0;
    }
    else {
      max_beta = beta;
      if (min_beta == -DBL_MAX || min_beta == DBL_MAX) beta /= 2.0;
      else                                             beta  = (beta + min_beta) / 2.0;
    }

    iter++;
  }

  *beta_out = used;
  *sum_out  = sum_P;
}


//======================================================
// Input similarities (P), exact mode.
//======================================================

// Calibrate each point against all others.
static void *exact_calib_worker(void *arg) {

  double *DD = ts_buf + (size_t) ((worker_t *) arg)->i * ts_n;
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {
      dist_row(i, DD);
      for (int m = 0; m < ts_n; m++) DD[m] = DD[m] * DD[m];
      calibrate(DD, ts_n, i, 1, ts_beta + i, ts_zsum + i);
    }
  }

  return NULL;
}

// P[i,j] = p(j|i) + p(i|j), for i < j, by column of the
// packed triangle.
static void *exact_pairs_worker(void *arg) {

  (void) arg;
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {
      size_t        s  = col_start(i);
      const double *d  = ts_dist + s;
      double       *P  = ts_P + s;
      double        bi = ts_beta[i], zi = ts_zsum[i];
      for (int j = i + 1; j < ts_n; j++) {
        double DD = d[j - i - 1] * d[j - i - 1];
        P[j - i - 1] = exp(-bi * DD) / zi + exp(-ts_beta[j] * DD) / ts_zsum[j];
      }
    }
  }

  return NULL;
}

static void exact_affinities(int n_threads) {

  int    n  = ts_n;
  size_t nd = (size_t) n * (n - 1) / 2;

  run_points(exact_calib_worker, n_threads);
  run_points(exact_pairs_worker, n_threads);

  // Normalize to sum 1. The total runs over the full
  // n x n matrix, row by row, as bhtsne sums it.
  double total = 0;
  for (int i = 0; i < n; i++) {
    for (int m = 0; m < n; m++) {
      if      (m < i) total += ts_P[pair_idx(m, i)];
      else if (m > i) total += ts_P[pair_idx(i, m)];
      else            total += DBL_MIN / ts_zsum[i];
    }
  }

  for (size_t x = 0; x < nd; x++) ts_P[x] /= total;
}


//======================================================
// Input similarities (P), Barnes-Hut mode.
//======================================================

// Max-heap on (distance, index), for the K nearest.
static inline int nn_less(double da, int ia, double db, int ib) {
  return da < db || (da == db && ia < ib);
}

static void heap_down(double *d, int *idx, int len, int at) {
  while (1) {
    int big = at, l = 2 * at + 1, r = l + 1;
    if (l < len && nn_less(d[idx[big]], idx[big], d[idx[l]], idx[l])) big = l;
    if (r < len && nn_less(d[idx[big]], idx[big], d[idx[r]], idx[r])) big = r;
    if (big == at) return;
    int t = idx[at]; idx[at] = idx[big]; idx[big] = t;
    at = big;
  }
}

// Each point's K nearest neighbours, nearest first, and
// their conditional probabilities p(j|i).
static void *sparse_calib_worker(void *arg) {

  int     t   = ((worker_t *) arg)->i;
  double *row = ts_buf  + (size_t) t * ts_n;
  int    *idx = ts_ibuf + (size_t) t * ts_K;
  int     K   = ts_K;
  int     lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {

      dist_row(i, row);

      // Select the K nearest, excluding i itself.
      int len = 0;
      for (int m = 0; m < ts_n; m++) {
        if (m == i) continue;
        if (len < K) {
          idx[len++] = m;
          if (len == K)
            for (int a = K / 2 - 1; a >= 0; a--) heap_down(row, idx, K, a);
        }
        else if (nn_less(row[m], m, row[idx[0]], idx[0])) {
          idx[0] = m;
          heap_down(row, idx, K, 0);
        }
      }

      // Heap sort: nearest first.
      for (int end = K - 1; end > 0; end--) {
        int tmp = idx[0]; idx[0] = idx[end]; idx[end] = tmp;
        heap_down(row, idx, end, 0);
      }

      size_t  off = (size_t) i * K;
      double *val = ts_val + off;
      int    *col = ts_col + off;

      for (int m = 0; m < K; m++) {
        col[m] = idx[m];
        val[m] = row[idx[m]];
      }

      double beta, sum_P;
      calibrate(val, K, -1, 0, &beta, &sum_P);
      ts_beta[i] = beta;
      ts_zsum[i] = sum_P;

      for (int m = 0; m < K; m++)
        val[m] = exp(-beta * val[m] * val[m]) / sum_P;
    }
  }

  return NULL;
}

// For each edge (i, col), where i appears in col's row.
static void *sparse_mate_worker(void *arg) {

  (void) arg;
  int K = ts_K;
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {
      for (size_t e = (size_t) i * K; e < (size_t) (i + 1) * K; e++) {
        int c = ts_col[e], at = -1;
        for (size_t f = (size_t) c * K; f < (size_t) (c + 1) * K; f++)
          if (ts_col[f] == i) at = (int) (f - (size_t) c * K);
        ts_mate[e] = at;
      }
    }
  }

  return NULL;
}

// Symmetrize: P = (P + t(P)) / 2, in bhtsne's order of
// entries within each row, then normalize to sum 1.
// Returns the number of entries.
static size_t sparse_affinities(int n_threads) {

  int    n   = ts_n;
  int    K   = ts_K;
  size_t nnz = (size_t) n * K;

  run_points(sparse_calib_worker, n_threads);
  run_points(sparse_mate_worker,  n_threads);

  const int    *col  = ts_col;
  const double *val  = ts_val;
  const int    *mate = ts_mate;

  // Row lengths.
  size_t *cnt = (size_t *) R_alloc(n, sizeof(size_t));
  memset(cnt, 0, n * sizeof(size_t));
  for (int i = 0; i < n; i++) {
    for (size_t e = (size_t) i * K; e < (size_t) (i + 1) * K; e++) {
      cnt[i]++;
      if (mate[e] < 0) cnt[col[e]]++;
    }
  }

  size_t *row = (size_t *) R_alloc((size_t) n + 1, sizeof(size_t));
  row[0] = 0;
  for (int i = 0; i < n; i++) row[i + 1] = row[i] + cnt[i];
  size_t total = row[n];

  int    *scol = (int *)    R_alloc(total, sizeof(int));
  double *sval = (double *) R_alloc(total, sizeof(double));

  memset(cnt, 0, n * sizeof(size_t)); // Now the fill offsets.
  for (int i = 0; i < n; i++) {
    for (size_t e = (size_t) i * K; e < (size_t) (i + 1) * K; e++) {
      int c = col[e];
      if (mate[e] >= 0 && i > c) continue; // Added from c's row.
      double v = (mate[e] >= 0) ? val[e] + val[(size_t) c * K + mate[e]] : val[e];
      scol[row[i] + cnt[i]] = c;
      scol[row[c] + cnt[c]] = i;
      sval[row[i] + cnt[i]] = v;
      sval[row[c] + cnt[c]] = v;
      cnt[i]++;
      cnt[c]++;
    }
  }

  for (size_t x = 0; x < total; x++) sval[x] /= 2.0;

  double sum_P = 0;
  for (size_t x = 0; x < total; x++) sum_P += sval[x];
  for (size_t x = 0; x < total; x++) sval[x] /= sum_P;

  (void) nnz;
  ts_row = row;
  ts_col = scol;
  ts_val = sval;

  return total;
}


//======================================================
// Space-partitioning tree (quadtree for k = 2, octree
// for k = 3, binary for k = 1), as bhtsne's SPTree,
// with node capacity 1. Built serially on the main
// thread; read concurrently by the workers.
//======================================================

typedef struct {
  double       corner[TS_MAXTREE];
  double       width[TS_MAXTREE];
  double       com[TS_MAXTREE];   // Center of mass.
  unsigned int cum;               // Points within.
  int          size;              // Points stored here: 0 or 1.
  int          index;             // The stored point.
  int          leaf;
  int          child;             // First of 2^k children.
} tnode_t;

static tnode_t *tr_pool = NULL;
static size_t   tr_len  = 0;
static size_t   tr_cap  = 0;
static int      tr_nch;           // Children per node: 2^k.

static void tree_free(void) {
  free(tr_pool);
  tr_pool = NULL;
  tr_len  = 0;
  tr_cap  = 0;
}

// Reserve m new nodes; returns the first one's index.
static int tree_new(int m) {

  if (tr_len + m > tr_cap) {
    size_t   cap = (tr_cap < 1024) ? 1024 : tr_cap * 2;
    while (cap < tr_len + m) cap *= 2;
    tnode_t *p = (tnode_t *) realloc(tr_pool, cap * sizeof(tnode_t));
    if (!p) { tree_free(); error("Insufficient memory."); } // # nocov
    tr_pool = p;
    tr_cap  = cap;
  }

  int first = (int) tr_len;
  tr_len += m;
  return first;
}

static void tree_init(int ni, const double *corner, const double *width) {
  tnode_t *nd = tr_pool + ni;
  for (int d = 0; d < ts_k; d++) {
    nd->corner[d] = corner[d];
    nd->width[d]  = width[d];
    nd->com[d]    = 0.0;
  }
  nd->cum   = 0;
  nd->size  = 0;
  nd->index = -1;
  nd->leaf  = 1;
  nd->child = -1;
}

static int tree_insert(int ni, int pi);

static void tree_subdivide(int ni) {

  int k  = ts_k;
  int ch = tree_new(tr_nch); // May move the pool.

  tnode_t *nd = tr_pool + ni;
  double corner[TS_MAXTREE], width[TS_MAXTREE];

  for (int i = 0; i < tr_nch; i++) {
    unsigned int div = 1;
    for (int d = 0; d < k; d++) {
      width[d] = .5 * nd->width[d];
      if ((i / div) % 2 == 1) corner[d] = nd->corner[d] - .5 * nd->width[d];
      else                    corner[d] = nd->corner[d] + .5 * nd->width[d];
      div *= 2;
    }
    tree_init(ch + i, corner, width);
  }

  nd->child = ch;

  // Move the stored point down. Children never need
  // subdividing for it, so the pool doesn't move.
  if (nd->size == 1) {
    int pi = nd->index, ok = 0;
    for (int j = 0; j < tr_nch; j++)
      if (!ok) ok = tree_insert(ch + j, pi);
  }

  nd         = tr_pool + ni;
  nd->index  = -1;
  nd->size   = 0;
  nd->leaf   = 0;
}

static int tree_insert(int ni, int pi) {

  int           k  = ts_k;
  tnode_t      *nd = tr_pool + ni;
  const double *pt = ts_Y + (size_t) pi * k;

  for (int d = 0; d < k; d++) {
    if (nd->corner[d] - nd->width[d] > pt[d]) return 0;
    if (nd->corner[d] + nd->width[d] < pt[d]) return 0;
  }

  // Running update of the size and center of mass.
  nd->cum++;
  double mult1 = (double) (nd->cum - 1) / (double) nd->cum;
  double mult2 = 1.0 / (double) nd->cum;
  for (int d = 0; d < k; d++)
    nd->com[d] = nd->com[d] * mult1 + mult2 * pt[d];

  if (nd->leaf && nd->size < 1) {
    nd->index = pi;
    nd->size  = 1;
    return 1;
  }

  // A duplicate of the stored point stays here, counted
  // only in `cum`.
  if (nd->size == 1) {
    const double *q = ts_Y + (size_t) nd->index * k;
    int dup = 1;
    for (int d = 0; d < k; d++)
      if (pt[d] != q[d]) { dup = 0; break; }
    if (dup) return 1;
  }

  if (nd->leaf) tree_subdivide(ni);

  int ch = tr_pool[ni].child;
  for (int i = 0; i < tr_nch; i++)
    if (tree_insert(ch + i, pi)) return 1;

  return 0; // # nocov (every point fits in some child)
}

// Build the tree on ts_Y.
static void tree_build(void) {

  int n = ts_n, k = ts_k;
  // Zero-filled so gcc can see tree_init() never reads unset entries.
  double mean[TS_MAXTREE] = {0}, mn[TS_MAXTREE], mx[TS_MAXTREE], width[TS_MAXTREE] = {0};

  for (int d = 0; d < k; d++) {
    mean[d] = 0.0;
    mn[d]   =  DBL_MAX;
    mx[d]   = -DBL_MAX;
  }

  for (int i = 0; i < n; i++) {
    const double *y = ts_Y + (size_t) i * k;
    for (int d = 0; d < k; d++) {
      mean[d] += y[d];
      if (y[d] < mn[d]) mn[d] = y[d];
      if (y[d] > mx[d]) mx[d] = y[d];
    }
  }

  for (int d = 0; d < k; d++) {
    mean[d] /= (double) n;
    double a = mx[d] - mean[d], b = mean[d] - mn[d];
    width[d] = ((a <= b) ? b : a) + 1e-5;
  }

  tr_nch = 1 << k;
  tr_len = 0;
  int root = tree_new(1);
  tree_init(root, mean, width);

  for (int i = 0; i < n; i++) tree_insert(root, i);
}

// Repulsive force on point pi, into neg[0..k-1]; returns
// its share of the normalization sum Z.
static double tree_repulse(int ni, int pi, double *neg) {

  const tnode_t *nd = tr_pool + ni;
  int    k  = ts_k;
  double rs = 0;

  if (nd->cum == 0 || (nd->leaf && nd->size == 1 && nd->index == pi))
    return rs;

  const double *pt = ts_Y + (size_t) pi * k;
  double buff[TS_MAXTREE], sq = .0;
  for (int d = 0; d < k; d++) {
    buff[d] = pt[d] - nd->com[d];
    sq     += buff[d] * buff[d];
  }

  double max_width = 0.0;
  for (int d = 0; d < k; d++)
    max_width = (max_width > nd->width[d]) ? max_width : nd->width[d];

  if (nd->leaf || max_width / sqrt(sq) < ts_theta) {
    sq = 1.0 / (1.0 + sq);
    double mult = nd->cum * sq;
    rs   += mult;
    mult *= sq;
    for (int d = 0; d < k; d++) neg[d] += mult * buff[d];
  }
  else {
    for (int i = 0; i < tr_nch; i++)
      rs += tree_repulse(nd->child + i, pi, neg);
  }

  return rs;
}


//======================================================
// Gradients. Like bhtsne, without the constant factor 4
// of the true gradient (it is absorbed by the learning
// rate).
//======================================================

// Exact mode. Squared distances between rows of Y are
// computed as bhtsne computes them (with BLAS's dgemm):
// |y_i|^2 + |y_j|^2 - 2 y_i . y_j. The unnormalized
// q_ij = 1 / (1 + that) are kept in ts_Q, packed like
// `dm`.
static void *exact_q_worker(void *arg) {

  (void) arg;
  int k = ts_k;
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {
      const double *a = ts_Y + (size_t) i * k;
      double       *Q = ts_Q + col_start(i);
      for (int j = i + 1; j < ts_n; j++) {
        const double *b = ts_Y + (size_t) j * k;
        double dot = 0;
        for (int d = 0; d < k; d++) dot += b[d] * a[d];
        Q[j - i - 1] = 1 / (1 + (-2.0 * dot + (ts_sq[i] + ts_sq[j])));
      }
    }
  }

  return NULL;
}

// Sum of a packed symmetric array over the full n x n
// matrix without its diagonal, starting from `s`. Each
// row is added up by one thread, in the order of its
// columns, reading the packed triangle in contiguous
// runs (as exact_grad_worker does); the row sums are
// then added in row order. The result is the same for
// any number of threads. (bhtsne adds up the whole
// matrix in one serial pass; this differs from it only
// in rounding.)
static const double *ps_A;

static void *packed_rowsum_worker(void *arg) {

  (void) arg;
  int    n = ts_n;
  int    lo, hi;
  size_t xs[TS_CHUNK];
  double rs[TS_CHUNK];

  while (next_chunk(&lo, &hi)) {

    for (int j = lo; j < hi; j++) {
      xs[j - lo] = col_start(j);
      rs[j - lo] = 0;
    }

    for (int m = 0; m < n; m++) {
      size_t base = (m < lo) ? col_start(m) + (size_t) (lo - m - 1) : 0;
      for (int j = lo; j < hi; j++) {
        size_t x;
        if      (m <  lo) x = base + (size_t) (j - lo);
        else if (m == j)  continue;
        else if (m <  j)  x = pair_idx(m, j);
        else              x = xs[j - lo]++;
        rs[j - lo] += ps_A[x];
      }
    }

    for (int j = lo; j < hi; j++) ts_out[j] = rs[j - lo];
  }

  return NULL;
}

static double packed_sum(const double *A, double s, int n_threads) {
  ps_A = A;
  run_points(packed_rowsum_worker, n_threads);
  for (int i = 0; i < ts_n; i++) s += ts_out[i];
  return s;
}

static void *exact_grad_worker(void *arg) {

  (void) arg;
  int    k    = ts_k;
  int    n    = ts_n;
  double sumQ = ts_sumQ;
  int    lo, hi;

  // Rows lo..hi-1 together, each adding up its terms in the
  // order of m, as bhtsne does. Going through m in the outer
  // loop reads P and Q in contiguous runs: for m < lo, from
  // column m of the packed triangle; for m > j, from column
  // j (one stream per row).
  size_t xs[TS_CHUNK]; // Row j's next position in column j.

  while (next_chunk(&lo, &hi)) {

    for (int j = lo; j < hi; j++) {
      for (int d = 0; d < k; d++) ts_dC[(size_t) j * k + d] = 0.0;
      xs[j - lo] = col_start(j);
    }

    for (int m = 0; m < n; m++) {

      const double *ym   = ts_Y + (size_t) m * k;
      size_t        base = (m < lo) ? col_start(m) + (size_t) (lo - m - 1) : 0;

      for (int j = lo; j < hi; j++) {

        size_t x;
        if      (m <  lo) x = base + (size_t) (j - lo);
        else if (m == j)  continue;
        else if (m <  j)  x = pair_idx(m, j);
        else              x = xs[j - lo]++;

        double       *dC   = ts_dC + (size_t) j * k;
        const double *yj   = ts_Y  + (size_t) j * k;
        double        q    = ts_Q[x];
        double        mult = (ts_P[x] - (q / sumQ)) * q;
        for (int d = 0; d < k; d++) dC[d] += (yj[d] - ym[d]) * mult;
      }
    }
  }

  return NULL;
}

static void exact_q(int n_threads) {

  for (int i = 0; i < ts_n; i++) {
    const double *y = ts_Y + (size_t) i * ts_k;
    double s = 0;
    for (int d = 0; d < ts_k; d++) s += y[d] * y[d];
    ts_sq[i] = s;
  }

  run_points(exact_q_worker, n_threads);
}

static void exact_gradient(int n_threads) {
  exact_q(n_threads);
  ts_sumQ = packed_sum(ts_Q, .0, n_threads);
  run_points(exact_grad_worker, n_threads);
}


// Barnes-Hut: each point's attractive (edge) and
// repulsive (tree) forces.
static void *bh_force_worker(void *arg) {

  (void) arg;
  int k = ts_k;
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {

      double       *pos = ts_pos + (size_t) i * k;
      double       *neg = ts_neg + (size_t) i * k;
      const double *yi  = ts_Y   + (size_t) i * k;

      for (int d = 0; d < k; d++) { pos[d] = 0.0; neg[d] = 0.0; }

      for (size_t e = ts_row[i]; e < ts_row[i + 1]; e++) {
        const double *ym = ts_Y + (size_t) ts_col[e] * k;
        double buff[TS_MAXTREE], sq = 1.0;
        for (int d = 0; d < k; d++) {
          buff[d] = yi[d] - ym[d];
          sq     += buff[d] * buff[d];
        }
        sq = ts_val[e] / sq;
        for (int d = 0; d < k; d++) pos[d] += sq * buff[d];
      }

      ts_out[i] = tree_repulse(0, i, neg);
    }
  }

  return NULL;
}

static void bh_gradient(int n_threads) {

  tree_build();
  run_points(bh_force_worker, n_threads);

  double sumQ = .0;
  for (int i = 0; i < ts_n; i++) sumQ += ts_out[i];

  size_t nk = (size_t) ts_n * ts_k;
  for (size_t x = 0; x < nk; x++) ts_dC[x] = ts_pos[x] - (ts_neg[x] / sumQ);
}


//======================================================
// Final cost: the KL divergence of Q from P, computed as
// bhtsne computes it. Each term is computed in parallel.
// Barnes-Hut terms are added up serially, in bhtsne's
// order; exact ones by packed_sum().
//======================================================

// Exact: ts_Q holds the unnormalized q_ij; replaced by
// the terms.
static void *exact_cost_worker(void *arg) {

  (void) arg;
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {
      size_t  s = col_start(i);
      double *Q = ts_Q + s;
      double *P = ts_P + s;
      for (int j = 0; j < ts_n - i - 1; j++) {
        double q = Q[j] / ts_sumQ;
        Q[j] = P[j] * log((P[j] + 1e-9) / (q + 1e-9));
      }
    }
  }

  return NULL;
}

// Barnes-Hut: each edge's term, into ts_neg's place
// (ts_term).
static double *ts_term;

static void *bh_cost_worker(void *arg) {

  (void) arg;
  int k = ts_k;
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {
      const double *yi = ts_Y + (size_t) i * k;
      for (size_t e = ts_row[i]; e < ts_row[i + 1]; e++) {
        const double *ym = ts_Y + (size_t) ts_col[e] * k;
        double q = .0;
        for (int d = 0; d < k; d++) q += (yi[d] - ym[d]) * (yi[d] - ym[d]);
        q = (1.0 / (1.0 + q)) / ts_sumQ;
        ts_term[e] = ts_val[e] * log((ts_val[e] + FLT_MIN) / (q + FLT_MIN));
      }
    }
  }

  return NULL;
}

static void *bh_zsum_worker(void *arg) {

  (void) arg;
  double neg[TS_MAXTREE];
  int lo, hi;

  while (next_chunk(&lo, &hi)) {
    for (int i = lo; i < hi; i++) {
      for (int d = 0; d < ts_k; d++) neg[d] = 0.0;
      ts_out[i] = tree_repulse(0, i, neg);
    }
  }

  return NULL;
}

static double kl_cost(int exact, int n_threads) {

  int n = ts_n;

  if (exact) {
    exact_q(n_threads);
    ts_sumQ = packed_sum(ts_Q, DBL_MIN, n_threads);
    run_points(exact_cost_worker, n_threads);
    return packed_sum(ts_Q, .0, n_threads);
  }

  tree_build();
  run_points(bh_zsum_worker, n_threads);
  double sumQ = .0;
  for (int i = 0; i < n; i++) sumQ += ts_out[i];
  ts_sumQ = sumQ;

  ts_term = (double *) R_alloc(ts_row[n], sizeof(double));
  run_points(bh_cost_worker, n_threads);

  double C = .0;
  for (size_t e = 0; e < ts_row[n]; e++) C += ts_term[e];
  return C;
}


//------------------------------------------------------
// Set up the shared state and compute P.
//------------------------------------------------------
static void setup(SEXP sexp_dist, int n, double perp, int exact, int n_threads) {

  ts_n    = n;
  ts_dist = REAL(sexp_dist);
  ts_perp = perp;
  ts_K    = (int) (3 * perp);

  // Barnes-Hut needs at least one neighbor per point.
  if (!exact && ts_K < 1)
    error("Barnes-Hut t-SNE needs a perplexity of at least 1/3.");

  ts_beta = (double *) R_alloc(n, sizeof(double));
  ts_zsum = (double *) R_alloc(n, sizeof(double));
  ts_buf  = (double *) R_alloc((size_t) n_threads * n, sizeof(double));

  if (exact) {
    ts_P = (double *) R_alloc((size_t) n * (n - 1) / 2, sizeof(double));
    exact_affinities(n_threads);
  }
  else {
    size_t nnz = (size_t) n * ts_K;
    ts_ibuf = (int *)    R_alloc((size_t) n_threads * ts_K, sizeof(int));
    ts_col  = (int *)    R_alloc(nnz, sizeof(int));
    ts_val  = (double *) R_alloc(nnz, sizeof(double));
    ts_mate = (int *)    R_alloc(nnz, sizeof(int));
    sparse_affinities(n_threads);
  }
}


// Gaussian deviate, by the polar method as in bhtsne.
static double randn(pcg32_random_t *rng) {
  double x, y, r;
  do {
    x = 2 * (pcg32_random_r(rng) / 4294967296.0) - 1;
    y = 2 * (pcg32_random_r(rng) / 4294967296.0) - 1;
    r = (x * x) + (y * y);
  } while ((r >= 1.0) || (r == 0.0));
  r = sqrt(-2 * log(r) / r);
  return x * r;
}

static inline double sign_of(double x) {
  return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0));
}


//======================================================
// R interface: run t-SNE.
//
// sexp_dist:  packed lower triangle of distances.
// sexp_init:  n x k starting coordinates, or NULL for
//             random ones from `seed`.
// sexp_iopts: c(maxit, stop_lying_iter, mom_switch_iter,
//               n_threads, seed)
// sexp_dopts: c(perplexity, theta, learning_rate,
//               momentum, final_momentum,
//               exaggeration_factor)
//
// Returns list(points = n x k matrix, kl_divergence).
//======================================================
SEXP C_tsne(
    SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_init,
    SEXP sexp_iopts, SEXP sexp_dopts ) {

  int n = asInteger(sexp_n);
  int k = asInteger(sexp_k);

  int maxit           = INTEGER(sexp_iopts)[0];
  int stop_lying_iter = INTEGER(sexp_iopts)[1];
  int mom_switch_iter = INTEGER(sexp_iopts)[2];
  int n_threads       = INTEGER(sexp_iopts)[3];
  int seed            = INTEGER(sexp_iopts)[4];

  double perp           = REAL(sexp_dopts)[0];
  double theta          = REAL(sexp_dopts)[1];
  double eta            = REAL(sexp_dopts)[2];
  double momentum       = REAL(sexp_dopts)[3];
  double final_momentum = REAL(sexp_dopts)[4];
  double exaggeration   = REAL(sexp_dopts)[5];

  int    exact = (theta == 0);
  size_t nk    = (size_t) n * k;

  if (n < 2 || k < 1 || n_threads < 1 || maxit < 0 || !(perp > 0) ||
      n - 1 < 3 * perp || (!exact && k > TS_MAXTREE) ||
      (size_t) XLENGTH(sexp_dist) != (size_t) n * (n - 1) / 2 ||
      (!isNull(sexp_init) && (size_t) XLENGTH(sexp_init) != nk))
    error("Invalid arguments to C_tsne()."); // # nocov

  tree_free(); // In case an interrupt left one behind.

  ts_k     = k;
  ts_theta = theta;

  setup(sexp_dist, n, perp, exact, n_threads);

  // Each iteration runs short parallel phases, which start their
  // threads anew. Use no more threads than their work keeps busy
  // for about a millisecond each: 4e5 pairs, or (Barnes-Hut) 1500
  // points, per thread.
  double units = exact ? 0.5 * n * (double) n / 4e5 : n / 1500.0;
  if (units + 1 < n_threads) n_threads = (int) units + 1;

  double *Y     = (double *) R_alloc(nk, sizeof(double));
  double *uY    = (double *) R_alloc(nk, sizeof(double));
  double *gains = (double *) R_alloc(nk, sizeof(double));
  ts_dC  = (double *) R_alloc(nk, sizeof(double));
  ts_out = (double *) R_alloc(n,  sizeof(double));
  ts_Y   = Y;

  if (exact) {
    ts_sq = (double *) R_alloc(n, sizeof(double));
    ts_Q  = (double *) R_alloc((size_t) n * (n - 1) / 2, sizeof(double));
  }
  else {
    ts_pos = (double *) R_alloc(nk, sizeof(double));
    ts_neg = (double *) R_alloc(nk, sizeof(double));
  }

  for (size_t x = 0; x < nk; x++) { uY[x] = .0; gains[x] = 1.0; }


  // Exaggerate P.
  size_t np = exact ? (size_t) n * (n - 1) / 2 : ts_row[n];
  double *P = exact ? ts_P : ts_val;
  for (size_t x = 0; x < np; x++) P[x] *= exaggeration;


  // Starting coordinates.
  if (isNull(sexp_init)) {
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (uint64_t) (int64_t) seed, 0x74736E65u); // "tsne"
    for (size_t x = 0; x < nk; x++) Y[x] = randn(&rng) * .0001;
  }
  else {
    const double *init = REAL(sexp_init); // n x k, column-major
    for (int i = 0; i < n; i++)
      for (int d = 0; d < k; d++)
        Y[(size_t) i * k + d] = init[(size_t) d * n + i];
  }


  //--------------------------------------------------
  // Gradient descent with momentum and gains.
  //--------------------------------------------------
  double *mean = (double *) R_alloc(k, sizeof(double));

  for (int iter = 0; iter < maxit; iter++) {

    if (iter == stop_lying_iter)
      for (size_t x = 0; x < np; x++) P[x] /= exaggeration;

    if (iter == mom_switch_iter) momentum = final_momentum;

    if (exact) exact_gradient(n_threads);
    else       bh_gradient(n_threads);

    double *dY = ts_dC;

    for (size_t x = 0; x < nk; x++)
      gains[x] = (sign_of(dY[x]) != sign_of(uY[x])) ? (gains[x] + .2) : (gains[x] * .8);
    for (size_t x = 0; x < nk; x++)
      if (gains[x] < .01) gains[x] = .01;

    for (size_t x = 0; x < nk; x++) uY[x] = momentum * uY[x] - eta * gains[x] * dY[x];
    for (size_t x = 0; x < nk; x++) Y[x]  = Y[x] + uY[x];

    // Center.
    for (int d = 0; d < k; d++) mean[d] = 0.0;
    for (int i = 0; i < n; i++)
      for (int d = 0; d < k; d++) mean[d] += Y[(size_t) i * k + d];
    for (int d = 0; d < k; d++) mean[d] /= (double) n;

    // If a sum is not finite, some coordinate is not: the optimization
    // diverged. (NaN coordinates would also make the Barnes-Hut tree
    // subdivide without end.)
    for (int d = 0; d < k; d++)
      if (!R_FINITE(mean[d])) {
        tree_free();
        error(
          "t-SNE diverged: coordinates became infinite or NaN. Try a "
          "smaller `learning_rate`, `momentum`, or `exaggeration_factor`." );
      }

    for (int i = 0; i < n; i++)
      for (int d = 0; d < k; d++) Y[(size_t) i * k + d] -= mean[d];

    if (iter % 10 == 9) R_CheckUserInterrupt();
  }

  double cost = kl_cost(exact, n_threads);

  tree_free();


  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  double *out = REAL(sexp_points);
  for (int i = 0; i < n; i++)
    for (int d = 0; d < k; d++)
      out[(size_t) d * n + i] = Y[(size_t) i * k + d];

  const char *names[] = { "points", "kl_divergence", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, ScalarReal(cost));

  UNPROTECT(2);
  return sexp_result;
}


//======================================================
// R interface, for tests: the input similarities P and
// each point's calibrated beta.
//
// Exact (theta == 0): list(beta, P), P packed as `dm`.
// Barnes-Hut:         list(beta, row, col, val), P in
//                     compressed sparse rows (0-based).
//======================================================
SEXP C_tsne_affinity(
    SEXP sexp_dist, SEXP sexp_n, SEXP sexp_perp, SEXP sexp_exact,
    SEXP sexp_n_threads ) {

  int    n         = asInteger(sexp_n);
  double perp      = asReal(sexp_perp);
  int    exact     = asLogical(sexp_exact);
  int    n_threads = asInteger(sexp_n_threads);

  if (n < 2 || !(perp > 0) || n - 1 < 3 * perp || n_threads < 1 ||
      (size_t) XLENGTH(sexp_dist) != (size_t) n * (n - 1) / 2)
    error("Invalid arguments to C_tsne_affinity()."); // # nocov

  setup(sexp_dist, n, perp, exact, n_threads);

  SEXP sexp_beta = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(sexp_beta), ts_beta, n * sizeof(double));

  SEXP sexp_result;

  if (exact) {
    size_t nd = (size_t) n * (n - 1) / 2;
    SEXP sexp_P = PROTECT(allocVector(REALSXP, (R_xlen_t) nd));
    memcpy(REAL(sexp_P), ts_P, nd * sizeof(double));
    const char *names[] = { "beta", "P", "" };
    sexp_result = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(sexp_result, 0, sexp_beta);
    SET_VECTOR_ELT(sexp_result, 1, sexp_P);
    UNPROTECT(3);
  }
  else {
    size_t nnz = ts_row[n];
    SEXP sexp_row = PROTECT(allocVector(REALSXP, (R_xlen_t) n + 1));
    SEXP sexp_col = PROTECT(allocVector(INTSXP,  (R_xlen_t) nnz));
    SEXP sexp_val = PROTECT(allocVector(REALSXP, (R_xlen_t) nnz));
    for (int i = 0; i <= n; i++) REAL(sexp_row)[i] = (double) ts_row[i];
    memcpy(INTEGER(sexp_col), ts_col, nnz * sizeof(int));
    memcpy(REAL(sexp_val),    ts_val, nnz * sizeof(double));
    const char *names[] = { "beta", "row", "col", "val", "" };
    sexp_result = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(sexp_result, 0, sexp_beta);
    SET_VECTOR_ELT(sexp_result, 1, sexp_row);
    SET_VECTOR_ELT(sexp_result, 2, sexp_col);
    SET_VECTOR_ELT(sexp_result, 3, sexp_val);
    UNPROTECT(5);
  }

  return sexp_result;
}
