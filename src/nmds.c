// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The SMACOF iteration in this file (metric and non-metric) is ported from
// scikit-learn 1.9.1: `_smacof_single()` in sklearn/manifold/_mds.py, and
// the isotonic regression it calls (sklearn/isotonic.py and
// sklearn/_isotonic.pyx). Those files are covered by this license:
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
 * Non-metric (nmds()) and metric (mds()) multidimensional scaling by
 * SMACOF.
 *
 * Each start iterates, from a starting configuration X (n x k):
 *
 *   1. Disparities. The first iteration uses the dissimilarities
 *      themselves; later ones, the isotonic (monotone) regression of
 *      the current configuration distances on the dissimilarities.
 *      Tied dissimilarities get one disparity, the mean of their
 *      distances (Kruskal's secondary approach, scikit-learn's); or,
 *      with Kruskal's primary approach, each tie group is sorted by
 *      distance and its pairs are fitted one by one. Pairs with a
 *      dissimilarity of exactly 0 have a disparity of 0. Disparities
 *      are scaled so their sum of squares is the number of pairs.
 *
 *   2. Guttman transform: X <- (1/n) B(X) X, with
 *      B(X) X row i = sum_j r_ij (x_i - x_j), r_ij = dhat_ij / d_ij
 *      (0 where d_ij = 0).
 *
 *   3. Raw stress, sum over pairs of (d_ij - dhat_ij)^2, with the new
 *      distances and the disparities of step 1. Stop when
 *      (previous stress - stress) / sum(d_ij^2) <= eps, or after
 *      maxit iterations.
 *
 * With the secondary approach and plain steps, this is scikit-learn's
 * non-metric `_smacof_single()`, step for step. Metric MDS
 * (scikit-learn's metric `_smacof_single()`) differs only in step 1:
 * the disparities are the dissimilarities themselves, unscaled. The
 * disparity stage, fit_disparities(), is skipped, and the pairs are
 * visited in dm's own order, unsorted. As the distances that give an
 * iteration's stress are the ones the next Guttman transform needs,
 * one pass over the pairs does both.
 *
 * Relaxed steps. The Guttman transform is scale-free (Gamma(cX) =
 * Gamma(X)), so the relaxed update X <- 2 Gamma(X) - X, which about
 * doubles the progress per iteration in the slowly converging
 * directions, does not converge in the scale direction (or, for k = 1,
 * in any direction: there, Gamma only depends on the points' order),
 * where it flips back and forth. So, optionally, and only for k >= 2,
 * a start takes relaxed steps only until one makes no more than `eps`
 * progress, and then plain steps until a plain step does. The last
 * iteration allowed by maxit is always plain. The configuration that
 * is measured and returned therefore always comes from a plain step.
 *
 * After the last iteration, each start's stress is measured: for
 * non-metric MDS, Kruskal's stress-1 with disparities fitted to the
 * final configuration; for metric MDS, sqrt(raw stress / sum of
 * squared dissimilarities). Starts are ranked by it.
 *
 * Shortcuts, as k is small and n can be large:
 *
 *   - Pairs are sorted by dissimilarity once, in prep_sort(). The
 *     sorted order, a pair of 16-bit sample indices per pair, and a
 *     bit per pair marking where tie groups start, are shared by all
 *     starts. Nothing n x n is ever formed.
 *
 *   - Each non-metric start keeps one double per pair, `y`:
 *     configuration distances, then disparities, in sorted order.
 *     Distances are recomputed from the coordinates when needed,
 *     rather than stored, and B(X) is never formed: the Guttman
 *     transform accumulates r_ij (x_i - x_j) pair by pair.
 *
 *   - The isotonic regression is pool-adjacent-violators, run in
 *     place on `y`, with the pooled blocks encoded in `y` itself, so
 *     it needs no extra memory.
 *
 *   - dm is read as dm * dscale, where dscale is a power of 2 that
 *     makes its largest value about 1, so that squares neither
 *     overflow nor underflow. Multiplying by a power of 2 is exact.
 *
 * Threads and determinism. The pairs are cut into a fixed number of
 * slices, which depends only on the number of pairs. Each slice's
 * sums (and its own copy of the Guttman accumulators) are computed in
 * order by one thread, and then added up slice by slice in order. So a
 * start gives bit-identical results however many threads run it, and
 * whether its slices run in parallel (one start at a time, "mode B")
 * or one after another on one thread (several starts at once, "mode
 * A"). The sort is a total order (ties broken by position in `dm`) and
 * runs in fixed chunks too, and so does the re-sorting of tie groups
 * for the primary approach (whose result is unique). Random starts
 * come from PCG streams keyed by `seed` and the start's number.
 */

#include "ecodive.h"

// Largest number of samples for non-metric MDS: indices must fit
// in 16 bits.
#define NMDS_MAX_N 65536

// Slices of the pairs: at most NMDS_SLICES, of at least
// NMDS_SLICE_MIN pairs each.
#define NMDS_SLICES    32
#define NMDS_SLICE_MIN 4096

// Mode A threads advance their starts at most this many iterations
// between checks for user interrupts.
#define NMDS_ROUND 10

// Sorting: chunks of `dm` for counting and scattering, and the
// target size of the buckets that are sorted individually.
#define NMDS_CHUNKS      64
#define NMDS_BUCKET      4096
#define NMDS_MAX_BUCKETS 16384
#define NMDS_OVERSAMPLE  32

// Parts of the sorted pairs, each a run of whole tie groups, in
// which the primary approach re-sorts the groups.
#define NMDS_PARTS 64

// Insertion sort a tie group while it moves fewer than this many
// elements per element; then merge sort it, in runs of NMDS_RUN.
#define NMDS_BUDGET 8
#define NMDS_RUN    32

// A dissimilarity starts a new tie group when it is at least
// NMDS_TIE_TOL times the largest dissimilarity above the group's
// first one, as scikit-learn compares with np.finfo(float64)
// .resolution (1e-15, absolute).
#define NMDS_TIE_TOL 1e-15

// Start states.
#define ST_INIT   0
#define ST_ITER   1
#define ST_FINISH 2
#define ST_DONE   3

// Slice operations. The last three are for metric MDS.
#define OP_GATHER  0
#define OP_GUTTMAN 1
#define OP_REDUCE  2
#define OP_DIST    3
#define OP_EVAL    4
#define OP_RESID   5
#define OP_MACC    6
#define OP_MFUSED  7
#define OP_MEVAL   8

#define BIT(bits, p) (((bits)[(p) >> 6] >> ((p) & 63)) & 1u)


//------------------------------------------------------
// Shared by all starts. Read-only once prepared.
//------------------------------------------------------
typedef struct {
  int             n;         // Samples.
  int             k;         // Dimensions.
  size_t          n_pairs;   // n * (n - 1) / 2.
  size_t          n_zero;    // Pairs with dissimilarity 0 (sorted first).
  const double   *dist;      // dm's packed lower triangle.
  double          dscale;    // dm is read as dm * dscale (a power of 2).
  uint32_t       *order;     // Sorted pairs, (j << 16) | i, i > j.
  uint64_t       *ties;      // Bit p: sorted pair p starts a tie group.
  int             n_slices;
  size_t          slice[NMDS_SLICES + 1];
  size_t          part[NMDS_PARTS + 1];   // Primary approach: parts.
  size_t          gmax;      // Largest tie group.
  double          sumsq;     // Sum of squared (scaled) dissimilarities.
  double          norm0;     // sqrt(n_pairs / sumsq).
  int             maxit;
  double          eps;
  int             primary;   // Primary approach to ties?
  int             relax;     // Relaxed steps (when k >= 2)?
  int             metric;    // Metric MDS: the disparities are dm itself.
  double          check;     // Pairs of work between interrupt checks.
} prep_t;


//------------------------------------------------------
// One start's working memory and state.
//------------------------------------------------------
typedef struct {
  double   *x;      // n x k, row-major: the configuration.
  double   *y;      // n_pairs: distances / disparities, sorted order.
  uint32_t *order;  // The pairs, in the order of y: P->order, or, for
                    // the primary approach, this start's own copy,
                    // re-sorted within tie groups by distance.
  double   *sy;     // Primary approach: merge sort scratch, gmax...
  uint32_t *su;     // ... elements each.
  double   *acc;    // n_slices x n x k: Guttman accumulators.
  double   *part;   // 3 x n_slices: per-slice partial sums.
  double   *trace;  // maxit raw stresses, or NULL.
  double    scale;  // Disparities are scale * y.
  double    stress; // Raw stress of the last iteration.
  double    old;    // ... and of the one before.
  double    ss;     // Sum of squared distances.
  double    ss_nz;  // ... over pairs with nonzero dissimilarity.
  double    stress1;// The stress used to rank starts, when done.
  int       iter;
  int       state;
  int       converged;
  int       degenerate;
  int       relaxing;   // Still taking relaxed steps?
  int       relax_step; // Is the current step relaxed?
  int       main;       // Running on R's thread (may check interrupts)?
} run_t;


//------------------------------------------------------
// Offset of element (j + 1, j) in a dist vector.
//------------------------------------------------------
static inline size_t col_start(int j, int n) {
  return (size_t) j * (size_t) n - (size_t) j * (size_t) (j + 1) / 2;
}

static inline size_t tri_index(uint32_t u, int n) {
  int i = (int) (u & 0xFFFFu), j = (int) (u >> 16);
  return col_start(j, n) + (size_t) (i - j - 1);
}

// The column j holding element t of a dist vector.
static int tri_col(size_t t, int n) {
  int lo = 0, hi = n - 2;
  while (lo < hi) {
    int mid = lo + (hi - lo + 1) / 2;
    if (col_start(mid, n) <= t) lo = mid; else hi = mid - 1;
  }
  return lo;
}

// Every `every` units of work on R's thread, check for interrupts.
static inline void maybe_check(size_t *count, size_t add, double every, int main) {
  *count += add;
  if (main && *count >= every) {
    *count = 0;
    R_CheckUserInterrupt();
  }
}



/*======================================================
 * Sorting the pairs by dissimilarity: a sample sort.
 *
 * Splitters are drawn from a fixed PCG sample of `dm`. Each
 * value goes to a bucket: strictly between two splitters
 * ("L" buckets, sorted later) or equal to one ("E" buckets,
 * already in order). Buckets are counted per fixed chunk of
 * `dm`, then filled in `dm`'s order, so equal values stay
 * in `dm`'s order. Each L bucket is then sorted by
 * (value, position in dm), a total order.
 *====================================================*/

typedef struct {
  double   key;
  uint32_t u;
} keyed_t;

static const double *s_dist;
static int           s_n;
static size_t        s_n_pairs;
static double        s_dscale;
static const double *s_split;    // n_split splitters, sorted.
static int           s_n_split;
static int           s_n_buckets;
static size_t       *s_counts;   // NMDS_CHUNKS x n_buckets.
static size_t       *s_pos;      // ... a working copy, for sort_scatter().
static size_t       *s_bucket;   // n_buckets + 1 bucket starts.
static size_t        s_max_l;    // Largest L bucket.
static keyed_t      *s_tmp;      // Per-thread scratch, s_max_l each.
static uint32_t     *s_order;
static uint64_t     *s_ties;
static size_t        s_n_zero;
static double        s_tol;
static double        s_chunk_ss[NMDS_CHUNKS];
static double        s_chunk_max[NMDS_CHUNKS];
static size_t        s_chunk_zero[NMDS_CHUNKS];

// Tie groups: each task's first position after a gap of at least
// s_tol from the previous value, and the first value of its last
// group as the task found it.
static size_t       *s_gap;
static double       *s_last;

static inline size_t chunk_lo(int c) {
  return (size_t) ((double) s_n_pairs * c / NMDS_CHUNKS);
}

static inline int bucket_of(double x) {
  int lo = 0, hi = s_n_split; // first splitter >= x
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (s_split[mid] < x) lo = mid + 1; else hi = mid;
  }
  return (lo < s_n_split && s_split[lo] == x) ? 2 * lo + 1 : 2 * lo;
}

static void *sort_count(void *arg) {
  int ti = ((worker_t *) arg)->i, tn = ((worker_t *) arg)->n;
  for (int c = ti; c < NMDS_CHUNKS; c += tn) {
    size_t *cnt = s_counts + (size_t) c * s_n_buckets;
    size_t  hi  = chunk_lo(c + 1), nz = 0;
    double  ss  = 0, mx = 0;
    memset(cnt, 0, (size_t) s_n_buckets * sizeof(size_t));
    for (size_t t = chunk_lo(c); t < hi; t++) {
      double x = s_dist[t], xs = x * s_dscale;
      cnt[bucket_of(x)]++;
      ss += xs * xs;
      if (x > mx)  mx = x;
      if (x == 0) nz++;
    }
    s_chunk_ss[c]   = ss;
    s_chunk_max[c]  = mx;
    s_chunk_zero[c] = nz;
  }
  return NULL;
}

static void *sort_scatter(void *arg) {
  int ti = ((worker_t *) arg)->i, tn = ((worker_t *) arg)->n;
  int n  = s_n;
  for (int c = ti; c < NMDS_CHUNKS; c += tn) {
    // The chunk's write offsets (a working copy of its counts).
    size_t *pos = s_pos + (size_t) c * s_n_buckets;
    memcpy(pos, s_counts + (size_t) c * s_n_buckets, (size_t) s_n_buckets * sizeof(size_t));
    size_t  t   = chunk_lo(c), hi = chunk_lo(c + 1);
    if (t >= hi) continue;
    int j = tri_col(t, n);
    int i = (int) (t - col_start(j, n)) + j + 1;
    for (; t < hi; t++) {
      s_order[pos[bucket_of(s_dist[t])]++] = ((uint32_t) j << 16) | (uint32_t) i;
      if (++i == n) { j++; i = j + 1; }
    }
  }
  return NULL;
}

static int cmp_keyed(const void *a, const void *b) {
  const keyed_t *x = (const keyed_t *) a, *y = (const keyed_t *) b;
  if (x->key < y->key) return -1;
  if (x->key > y->key) return  1;
  return (x->u > y->u) - (x->u < y->u);
}

static void *sort_buckets(void *arg) {
  int ti = ((worker_t *) arg)->i, tn = ((worker_t *) arg)->n;
  keyed_t *tmp = s_tmp + (size_t) ti * s_max_l;
  for (int b = 2 * ti; b < s_n_buckets; b += 2 * tn) { // L buckets only
    size_t lo = s_bucket[b], m = s_bucket[b + 1] - lo;
    if (m < 2) continue;
    uint32_t *ord = s_order + lo;
    for (size_t t = 0; t < m; t++) {
      tmp[t].u   = ord[t];
      tmp[t].key = s_dist[tri_index(ord[t], s_n)];
    }
    qsort(tmp, m, sizeof(keyed_t), cmp_keyed);
    for (size_t t = 0; t < m; t++) ord[t] = tmp[t].u;
  }
  return NULL;
}

// Tie group starts, 64 bits per word, 1024 words per task. Each task
// starts a group at its first position; that is corrected later, up
// to the task's first gap of at least s_tol from the previous value,
// where a group starts whatever the group's first value was.
static inline double sorted_value(size_t p) {
  return s_dist[tri_index(s_order[p], s_n)];
}

static void *sort_ties(void *arg) {
  int    ti = ((worker_t *) arg)->i, tn = ((worker_t *) arg)->n;
  size_t n_words = (s_n_pairs + 63) / 64;
  size_t n_tasks = (n_words + 1023) / 1024;
  for (size_t task = (size_t) ti; task < n_tasks; task += (size_t) tn) {
    size_t w0 = task * 1024, w1 = w0 + 1024;
    if (w1 > n_words) w1 = n_words;
    size_t p     = w0 * 64, hi = w1 * 64;
    double prev  = (p > 0) ? sorted_value(p - 1) : 0;
    double first = 0;
    size_t gap   = (hi < s_n_pairs) ? hi : s_n_pairs;
    for (size_t w = w0; w < w1; w++) {
      uint64_t word = 0;
      for (int b = 0; b < 64 && p < s_n_pairs; b++, p++) {
        double x = sorted_value(p);
        if (p >= s_n_zero) {
          int is_gap = (p == s_n_zero) || (x > prev && x - prev >= s_tol);
          if (is_gap && p < gap) gap = p;
          if (p == w0 * 64 || is_gap || (x > first && x - first >= s_tol)) {
            word |= (uint64_t) 1 << b;
            first = x;
          }
        }
        prev = x;
      }
      s_ties[w] = word;
    }
    s_gap[task]  = gap;
    s_last[task] = first;
  }
  return NULL;
}

// Serially, redo each task's positions before its first gap, with
// the group first value carried over from the task before.
static void fix_ties(void) {
  size_t n_words = (s_n_pairs + 63) / 64;
  size_t n_tasks = (n_words + 1023) / 1024;
  double first   = 0;
  for (size_t task = 0; task < n_tasks; task++) {
    size_t lo = task * 1024 * 64;
    for (size_t p = (lo > s_n_zero) ? lo : s_n_zero; p < s_gap[task]; p++) {
      double   x   = sorted_value(p);
      uint64_t bit = (uint64_t) 1 << (p & 63);
      if (x > first && x - first >= s_tol) { s_ties[p >> 6] |=  bit; first = x; }
      else                    { s_ties[p >> 6] &= ~bit;            }
    }
    if (s_gap[task] < lo + 1024 * 64 && s_gap[task] < s_n_pairs) first = s_last[task];
  }
}

static void prep_sort(prep_t *P, int n_threads) {

  size_t N = P->n_pairs;

  // No pass has more than NMDS_CHUNKS parts worth a thread.
  if (n_threads > NMDS_CHUNKS) n_threads = NMDS_CHUNKS;

  s_dist    = P->dist;
  s_n       = P->n;
  s_n_pairs = N;
  s_dscale  = P->dscale;

  //--------------------------------------------------
  // Splitters from a fixed random sample. Metric MDS
  // only needs the sums: one bucket, and no sort.
  //--------------------------------------------------
  int B = P->metric ? 1 : (int) (N / NMDS_BUCKET);
  if (B < 1) B = 1;
  if (B > NMDS_MAX_BUCKETS) B = NMDS_MAX_BUCKETS;

  int     n_split = B - 1;
  double *split   = (double *) R_alloc(n_split > 0 ? n_split : 1, sizeof(double));

  if (n_split > 0) {
    int     n_sample = B * NMDS_OVERSAMPLE;
    double *sample   = (double *) R_alloc(n_sample, sizeof(double));
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, 20261117u, 7u);
    for (int s = 0; s < n_sample; s++)
      sample[s] = P->dist[((uint64_t) pcg32_random_r(&rng) * N) >> 32];
    R_qsort(sample, 1, n_sample);
    for (int b = 1; b < B; b++)
      split[b - 1] = sample[(size_t) b * n_sample / B];
  }

  s_split     = split;
  s_n_split   = n_split;
  s_n_buckets = 2 * B - 1;


  //--------------------------------------------------
  // Count each chunk's buckets; then turn the counts
  // into write offsets, by bucket, then chunk.
  //--------------------------------------------------
  s_counts = (size_t *) R_alloc((size_t) NMDS_CHUNKS * s_n_buckets, sizeof(size_t));
  s_bucket = (size_t *) R_alloc((size_t) s_n_buckets + 1, sizeof(size_t));
  run_parallel(sort_count, n_threads, (int) (N > 1000000 ? 1000000 : N));

  P->sumsq = 0;
  double mx = 0;
  size_t nz = 0;
  for (int c = 0; c < NMDS_CHUNKS; c++) {
    P->sumsq += s_chunk_ss[c];
    if (s_chunk_max[c] > mx) mx = s_chunk_max[c];
    nz += s_chunk_zero[c];
  }

  if (P->metric) return;
  R_CheckUserInterrupt();

  size_t pos = 0;
  s_max_l = 0;
  for (int b = 0; b < s_n_buckets; b++) {
    s_bucket[b] = pos;
    for (int c = 0; c < NMDS_CHUNKS; c++) {
      size_t *cnt = s_counts + (size_t) c * s_n_buckets + b;
      size_t  m   = *cnt;
      *cnt = pos;
      pos += m;
    }
    if (b % 2 == 0 && pos - s_bucket[b] > s_max_l) s_max_l = pos - s_bucket[b];
  }
  s_bucket[s_n_buckets] = pos;


  //--------------------------------------------------
  // Scatter, then sort the L buckets.
  //--------------------------------------------------
  P->order = (uint32_t *) R_alloc(N, sizeof(uint32_t));
  s_order  = P->order;
  s_pos    = (size_t *) R_alloc((size_t) NMDS_CHUNKS * s_n_buckets, sizeof(size_t));
  run_parallel(sort_scatter, n_threads, (int) (N > 1000000 ? 1000000 : N));
  R_CheckUserInterrupt();

  int nt = (n_threads < B) ? n_threads : B; // one L bucket per thread at most
  s_tmp  = (keyed_t *) R_alloc((size_t) nt * (s_max_l > 0 ? s_max_l : 1), sizeof(keyed_t));
  run_parallel(sort_buckets, nt, (int) (N > 1000000 ? 1000000 : N));
  R_CheckUserInterrupt();


  //--------------------------------------------------
  // Tie groups.
  //--------------------------------------------------
  size_t n_tasks = ((N + 63) / 64 + 1023) / 1024;
  s_n_zero = nz;
  s_tol    = NMDS_TIE_TOL * mx;
  s_gap    = (size_t *) R_alloc(n_tasks, sizeof(size_t));
  s_last   = (double *) R_alloc(n_tasks, sizeof(double));
  P->ties  = (uint64_t *) R_alloc((N + 63) / 64, sizeof(uint64_t));
  s_ties   = P->ties;
  run_parallel(sort_ties, n_threads, (int) (N > 1000000 ? 1000000 : N));
  fix_ties();

  P->n_zero = nz;


  //--------------------------------------------------
  // Parts of whole tie groups, of about equal size, and
  // the largest group.
  //--------------------------------------------------
  size_t gstart = nz;
  int    pt     = 0;
  P->part[0] = nz;
  P->gmax    = 0;
  for (size_t p = nz + 1; p <= N; p++) {
    if (p < N && !BIT(P->ties, p)) continue;
    size_t g = p - gstart; // the group [gstart, p)
    if (g > P->gmax) P->gmax = g;
    gstart = p;
    while (pt < NMDS_PARTS - 1 && p - nz >= (N - nz) * (double) (pt + 1) / NMDS_PARTS)
      P->part[++pt] = p;
  }
  P->part[NMDS_PARTS] = N;
}



/*======================================================
 * Per-slice operations on one start.
 *====================================================*/

// Copy each pair's (scaled) dissimilarity into y, in sorted order.
static void op_gather(const prep_t *P, run_t *R, int s) {
  size_t lo = P->slice[s], hi = P->slice[s + 1];
  if (R->order != P->order)
    memcpy(R->order + lo, P->order + lo, (hi - lo) * sizeof(uint32_t));
  if (lo < P->n_zero) lo = P->n_zero;
  for (size_t p = lo; p < hi; p++)
    R->y[p] = P->dist[tri_index(R->order[p], P->n)] * P->dscale;
}

// This slice's share of B(X) X.
static inline void guttman_k(const prep_t *P, run_t *R, int s, const int k) {

  int     n  = P->n;
  double *a  = R->acc + (size_t) s * n * k;
  size_t  lo = P->slice[s], hi = P->slice[s + 1];
  if (lo < P->n_zero) lo = P->n_zero;

  const double   *x   = R->x;
  const double   *y   = R->y;
  const uint32_t *ord = R->order;
  double          sc  = R->scale;

  memset(a, 0, (size_t) n * k * sizeof(double));

  for (size_t p = lo; p < hi; p++) {
    uint32_t      u  = ord[p];
    int           i  = (int) (u & 0xFFFFu), j = (int) (u >> 16);
    const double *xi = x + (size_t) i * k;
    const double *xj = x + (size_t) j * k;
    double        d2 = 0;
    for (int c = 0; c < k; c++) { double t = xi[c] - xj[c]; d2 += t * t; }
    if (d2 > 0) {
      double  r  = sc * y[p] / sqrt(d2);
      double *ai = a + (size_t) i * k;
      double *aj = a + (size_t) j * k;
      for (int c = 0; c < k; c++) {
        double t = r * (xi[c] - xj[c]);
        ai[c] += t;
        aj[c] -= t;
      }
    }
  }
}

static void op_guttman(const prep_t *P, run_t *R, int s) {
  switch (P->k) {
    case 1:  guttman_k(P, R, s, 1);    break;
    case 2:  guttman_k(P, R, s, 2);    break;
    case 3:  guttman_k(P, R, s, 3);    break;
    default: guttman_k(P, R, s, P->k); break;
  }
}

// x = (1/n) * (sum of the slices' accumulators), for this
// slice's share of the rows: the Guttman transform. A relaxed
// step moves x twice as far: x = 2 * (Guttman transform) - x.
static void op_reduce(const prep_t *P, run_t *R, int s) {
  int    n  = P->n, k = P->k, S = P->n_slices;
  size_t nk = (size_t) n * k;
  size_t lo = nk * s / S, hi = nk * (s + 1) / S;
  double inv_n = 1.0 / n;
  for (size_t e = lo; e < hi; e++) {
    double sum = 0;
    for (int t = 0; t < S; t++) sum += R->acc[(size_t) t * nk + e];
    double g = inv_n * sum;
    R->x[e] = R->relax_step ? 2.0 * g - R->x[e] : g;
  }
}

// Distances of the new configuration: raw stress against the
// current disparities (unless `eval`), sums of squares, and y = d.
static inline void dist_k(const prep_t *P, run_t *R, int s, const int k, int eval) {

  size_t lo = P->slice[s], hi = P->slice[s + 1], nz = P->n_zero;
  const double   *x   = R->x;
  double         *y   = R->y;
  const uint32_t *ord = R->order;
  double          sc  = R->scale;
  double          st = 0, ss = 0, ss_nz = 0;

  for (size_t p = lo; p < hi; p++) {
    uint32_t      u  = ord[p];
    const double *xi = x + (size_t) (u & 0xFFFFu) * k;
    const double *xj = x + (size_t) (u >> 16) * k;
    double        d2 = 0;
    for (int c = 0; c < k; c++) { double t = xi[c] - xj[c]; d2 += t * t; }
    double d = sqrt(d2);
    if (p < nz) {
      st += d2;
    } else {
      if (!eval) { double e = d - sc * y[p]; st += e * e; }
      y[p]   = d;
      ss_nz += d2;
    }
    ss += d2;
  }

  R->part[s]                   = st;
  R->part[P->n_slices + s]     = ss;
  R->part[2 * P->n_slices + s] = ss_nz;
}

static void op_dist(const prep_t *P, run_t *R, int s, int eval) {
  switch (P->k) {
    case 1:  dist_k(P, R, s, 1,    eval); break;
    case 2:  dist_k(P, R, s, 2,    eval); break;
    case 3:  dist_k(P, R, s, 3,    eval); break;
    default: dist_k(P, R, s, P->k, eval); break;
  }
}

// Residual sum of squares against fitted (unscaled) disparities.
static void op_resid(const prep_t *P, run_t *R, int s) {

  int    k  = P->k;
  size_t lo = P->slice[s], hi = P->slice[s + 1], nz = P->n_zero;
  double rss = 0;

  for (size_t p = lo; p < hi; p++) {
    uint32_t      u  = R->order[p];
    const double *xi = R->x + (size_t) (u & 0xFFFFu) * k;
    const double *xj = R->x + (size_t) (u >> 16) * k;
    double        d2 = 0;
    for (int c = 0; c < k; c++) { double t = xi[c] - xj[c]; d2 += t * t; }
    if (p < nz) { rss += d2; }
    else        { double e = sqrt(d2) - R->y[p]; rss += e * e; }
  }

  R->part[s] = rss;
}

// Metric MDS, over the slice's pairs in dm's order: the
// configuration's distances and their raw stress against dm
// (if `stress`), and this slice's share of B(X) X (if `gut`).
// The arithmetic is that of op_dist() and op_guttman().
static inline void metric_k(const prep_t *P, run_t *R, int s, const int k, int gut, int stress) {

  int     n  = P->n, S = P->n_slices;
  double *a  = R->acc + (size_t) s * n * k;
  size_t  lo = P->slice[s], hi = P->slice[s + 1];

  const double *x   = R->x;
  const double *dm  = P->dist;
  double        sc  = R->scale;
  double        st  = 0, ss = 0;

  if (gut) memset(a, 0, (size_t) n * k * sizeof(double));

  int j = tri_col(lo, n);
  int i = (int) (lo - col_start(j, n)) + j + 1;

  for (size_t t = lo; t < hi; t++) {
    const double *xi = x + (size_t) i * k;
    const double *xj = x + (size_t) j * k;
    double        d2 = 0, dh = sc * dm[t];
    for (int c = 0; c < k; c++) { double v = xi[c] - xj[c]; d2 += v * v; }
    double d = sqrt(d2);
    if (stress) {
      double e = d - dh;
      st += e * e;
      ss += d2;
    }
    if (gut && d2 > 0) {
      double  r  = dh / d;
      double *ai = a + (size_t) i * k;
      double *aj = a + (size_t) j * k;
      for (int c = 0; c < k; c++) {
        double v = r * (xi[c] - xj[c]);
        ai[c] += v;
        aj[c] -= v;
      }
    }
    if (++i == n) { j++; i = j + 1; }
  }

  if (stress) {
    R->part[s]         = st;
    R->part[S + s]     = ss;
    R->part[2 * S + s] = ss;
  }
}

static void op_metric(const prep_t *P, run_t *R, int s, int gut, int stress) {
  switch (P->k) {
    case 1:  metric_k(P, R, s, 1,    gut, stress); break;
    case 2:  metric_k(P, R, s, 2,    gut, stress); break;
    case 3:  metric_k(P, R, s, 3,    gut, stress); break;
    default: metric_k(P, R, s, P->k, gut, stress); break;
  }
}

static void do_op(int op, const prep_t *P, run_t *R, int s) {
  switch (op) {
    case OP_GATHER:  op_gather(P, R, s);        break;
    case OP_GUTTMAN: op_guttman(P, R, s);       break;
    case OP_REDUCE:  op_reduce(P, R, s);        break;
    case OP_DIST:    op_dist(P, R, s, 0);       break;
    case OP_EVAL:    op_dist(P, R, s, 1);       break;
    case OP_RESID:   op_resid(P, R, s);         break;
    case OP_MACC:    op_metric(P, R, s, 1, 0);  break;
    case OP_MFUSED:  op_metric(P, R, s, 1, 1);  break;
    default:         op_metric(P, R, s, 0, 1);  break; // OP_MEVAL
  }
}


// Mode B: one start's slices, spread over threads, in waves of
// slices between interrupt checks when slices are large. run_parallel()
// runs each share exactly once, so each slice is done once.
static const prep_t *g_P;
static run_t        *g_R;
static int           g_op;
static int           g_lo, g_hi;

static void *slice_worker(void *arg) {
  int ti = ((worker_t *) arg)->i, tn = ((worker_t *) arg)->n;
  for (int s = g_lo + ti; s < g_hi; s += tn) {
    do_op(g_op, g_P, g_R, s);
  }
  return NULL;
}

// n_threads: 0 on a mode A worker thread; 1 or more on R's thread.
static void each_slice(int op, const prep_t *P, run_t *R, int n_threads) {

  int    S    = P->n_slices;
  double each = (double) P->n_pairs / S;

  if (n_threads > 1) {
    // Enough slices per wave to keep every thread busy, and
    // about P->check pairs of work.
    int wave = (int) (P->check / each);
    if (wave < n_threads) wave = n_threads;
    g_op = op; g_P = P; g_R = R;
    for (g_lo = 0; g_lo < S; g_lo += wave) {
      g_hi = (g_lo + wave < S) ? g_lo + wave : S;
      run_parallel(slice_worker, n_threads, 100);
      if (g_hi < S) R_CheckUserInterrupt();
    }
  }
  else {
    size_t count = 0;
    for (int s = 0; s < S; s++) {
      do_op(op, P, R, s);
      maybe_check(&count, (size_t) each, P->check, n_threads);
    }
  }
}

static double sum_part(const prep_t *P, const run_t *R, int which) {
  double sum = 0;
  for (int s = 0; s < P->n_slices; s++) sum += R->part[which * P->n_slices + s];
  return sum;
}



/*======================================================
 * Primary approach to ties: sort each tie group's
 * distances (with their pairs) into increasing order,
 * ties broken by pair, so that the isotonic regression
 * may fit them with different disparities. The result is
 * the one order sorted by (distance, pair), however it is
 * reached: insertion sort while the group is nearly in
 * order, as it is once iterations settle, else a merge
 * sort. Parts of whole groups are sorted in parallel.
 *====================================================*/
#define LESS(ya, ua, yb, ub) ((ya) < (yb) || ((ya) == (yb) && (ua) < (ub)))

static void insertion_sort(double *y, uint32_t *u, size_t g) {
  for (size_t a = 1; a < g; a++) {
    double   v = y[a];
    uint32_t w = u[a];
    size_t   b = a;
    for (; b > 0 && LESS(v, w, y[b - 1], u[b - 1]); b--) {
      y[b] = y[b - 1];
      u[b] = u[b - 1];
    }
    y[b] = v;
    u[b] = w;
  }
}

// Natural merge sort, using sy and su (g elements each): blocks of
// NMDS_RUN are insertion sorted, and then each pass merges pairs of
// adjacent ascending runs, so runs already in order cost little.
static void merge_sort(double *y, uint32_t *u, size_t g, double *sy, uint32_t *su) {

  for (size_t a = 0; a < g; a += NMDS_RUN)
    insertion_sort(y + a, u + a, (g - a < NMDS_RUN) ? g - a : NMDS_RUN);

  double   *fy = y, *ty = sy;
  uint32_t *fu = u, *tu = su;

  for (;;) {
    size_t a = 0, runs = 0;
    while (a < g) {
      size_t m = a + 1;
      while (m < g && !LESS(fy[m], fu[m], fy[m - 1], fu[m - 1])) m++;
      size_t e = (m < g) ? m + 1 : g;
      while (e < g && !LESS(fy[e], fu[e], fy[e - 1], fu[e - 1])) e++;
      size_t i = a, j = m, o = a;
      while (i < m && j < e) {
        if (LESS(fy[j], fu[j], fy[i], fu[i])) { ty[o] = fy[j]; tu[o++] = fu[j++]; }
        else                                  { ty[o] = fy[i]; tu[o++] = fu[i++]; }
      }
      while (i < m) { ty[o] = fy[i]; tu[o++] = fu[i++]; }
      while (j < e) { ty[o] = fy[j]; tu[o++] = fu[j++]; }
      runs++;
      a = e;
    }
    double *t1 = fy; fy = ty; ty = t1;
    uint32_t *t2 = fu; fu = tu; tu = t2;
    if (runs == 1) break;
  }

  if (fy != y) {
    memcpy(y, fy, g * sizeof(double));
    memcpy(u, fu, g * sizeof(uint32_t));
  }
}

static void sort_group(double *y, uint32_t *u, size_t g, double *sy, uint32_t *su) {

  size_t moves = 0, budget = NMDS_BUDGET * g;

  for (size_t a = 1; a < g; a++) {
    double   v = y[a];
    uint32_t w = u[a];
    size_t   b = a;
    for (; b > 0 && LESS(v, w, y[b - 1], u[b - 1]); b--) {
      y[b] = y[b - 1];
      u[b] = u[b - 1];
    }
    y[b]   = v;
    u[b]   = w;
    moves += a - b;
    if (moves > budget) {
      merge_sort(y, u, g, sy, su);
      return;
    }
  }
}

#undef LESS

// Sort the groups of part `pt`; on R's thread if `main`.
static void sort_part(const prep_t *P, run_t *R, int pt, double *sy, uint32_t *su, int main) {
  size_t N = P->n_pairs, count = 0;
  for (size_t p = P->part[pt]; p < P->part[pt + 1]; ) {
    size_t q = p + 1;
    while (q < N && !BIT(P->ties, q)) q++;
    if (q - p > 1) sort_group(R->y + p, R->order + p, q - p, sy, su);
    maybe_check(&count, q - p, P->check, main);
    p = q;
  }
}

static double   *g_sy;
static uint32_t *g_su;

static void *part_worker(void *arg) {
  int ti = ((worker_t *) arg)->i, tn = ((worker_t *) arg)->n;
  for (int pt = ti; pt < NMDS_PARTS; pt += tn)
    sort_part(g_P, g_R, pt, g_sy + (size_t) ti * g_P->gmax, g_su + (size_t) ti * g_P->gmax, 0);
  return NULL;
}

// n_threads > 1: parallel, with scratch for that many threads
// in g_sy and g_su.
static void sort_groups(const prep_t *P, run_t *R, int n_threads) {
  if (n_threads > 1) {
    g_P = P; g_R = R;
    run_parallel(part_worker, n_threads, 100);
  }
  else {
    for (int pt = 0; pt < NMDS_PARTS; pt++) sort_part(P, R, pt, R->sy, R->su, R->main);
  }
}



/*======================================================
 * Isotonic regression of y[n_zero..n_pairs) on the
 * sorted dissimilarities, in place, by pool-adjacent-
 * violators. Each tie group (or, for the primary
 * approach, each pair) enters as one block, holding
 * its sum and count. While being built, a block [a, b]
 * is stored in y itself: y[a] holds its sum and, if
 * b > a, y[b] holds minus its length. (Distances are
 * never negative, so a negative y[b] marks a length.)
 * Returns the sum of squared fitted values.
 *====================================================*/
static double pava(const prep_t *P, run_t *R, int sort_threads) {

  size_t nz = P->n_zero, N = P->n_pairs, count = 0;
  double *y = R->y;
  const uint64_t *ties = P->ties;
  int primary = P->primary;

  if (primary) sort_groups(P, R, sort_threads);

  for (size_t p = nz; p < N; ) {

    // The next tie group, [p, q).
    double s = y[p];
    size_t q = p + 1;
    if (!primary) while (q < N && !BIT(ties, q)) s += y[q++];
    size_t len   = q - p;
    size_t start = p;

    // Pool with the blocks before it while they are not lower.
    while (start > nz) {
      size_t e  = start - 1;
      double v  = y[e];
      size_t lt = (v < 0) ? (size_t) (-v) : 1;
      double st = (v < 0) ? y[e - lt + 1] : v;
      if (st * (double) len < s * (double) lt) break; // st / lt < s / len
      s     += st;
      len   += lt;
      start  = e - lt + 1;
    }

    y[start] = s;
    if (len > 1) y[q - 1] = -(double) len;
    maybe_check(&count, q - p, P->check, R->main);
    p = q;
  }

  // Expand each block's mean, from the last block back.
  double sumsq = 0;
  for (size_t e = N; e > nz; ) {
    double v    = y[e - 1];
    size_t len  = (v < 0) ? (size_t) (-v) : 1;
    double s    = (v < 0) ? y[e - len] : v;
    double mean = s / (double) len;
    for (size_t t = e - len; t < e; t++) y[t] = mean;
    sumsq += (double) len * mean * mean;
    e     -= len;
  }

  return sumsq;
}



/*======================================================
 * The disparity stage. On entry, y holds the current
 * configuration distances in sorted order; on return, it
 * holds the disparities fitted to them (isotonic
 * regression), and `scale` makes their sum of squares the
 * number of pairs. Everything else in an iteration
 * (distances, stress, the Guttman transform, the starts,
 * and threading) only reads disparities as scale * y.
 *====================================================*/
static void fit_disparities(const prep_t *P, run_t *R, int sort_threads) {
  double sumsq = pava(P, R, sort_threads);
  R->scale = sqrt((double) P->n_pairs / sumsq);
}



/*======================================================
 * One step of a start: initialization, one SMACOF
 * iteration, or the final measurement. n_threads is 0 on
 * a mode A worker thread, else the threads for passes
 * over the pairs (mode B); sort_threads, for re-sorting
 * tie groups.
 *====================================================*/
static void step(const prep_t *P, run_t *R, int n_threads, int sort_threads) {

  R->main = (n_threads > 0);

  if (R->state == ST_INIT) {

    if (P->metric) {
      R->scale = P->dscale;
    }
    else {
      each_slice(OP_GATHER, P, R, n_threads);
      R->scale = P->norm0;
    }

    if (P->maxit == 0) {
      each_slice(P->metric ? OP_MEVAL : OP_EVAL, P, R, n_threads);
      R->stress = P->metric ? sum_part(P, R, 0) : NAN;
      R->ss     = sum_part(P, R, 1);
      R->ss_nz  = sum_part(P, R, 2);
      R->old    = NAN;
      R->degenerate = (R->ss_nz == 0);
      R->state  = ST_FINISH;
    }
    else {
      R->relaxing = P->relax && P->k > 1;
      R->state    = ST_ITER;
    }
  }

  else if (R->state == ST_ITER) {

    // Relaxed, unless it is the first or last iteration.
    R->relax_step = R->relaxing && R->iter > 0 && R->iter < P->maxit - 1;

    // For metric MDS after the first iteration, the previous
    // OP_MFUSED pass has already accumulated B(X) X.
    if (!P->metric)
      each_slice(OP_GUTTMAN, P, R, n_threads);
    else if (R->iter == 0)
      each_slice(OP_MACC, P, R, n_threads);
    each_slice(OP_REDUCE, P, R, n_threads);
    each_slice(P->metric ? OP_MFUSED : OP_DIST, P, R, n_threads);

    R->stress = sum_part(P, R, 0);
    R->ss     = sum_part(P, R, 1);
    R->ss_nz  = sum_part(P, R, 2);
    if (R->trace) R->trace[R->iter] = R->stress;
    R->iter++;

    int stop = 0;
    if (R->ss_nz == 0) {
      // Every pair with a nonzero dissimilarity has collapsed
      // to distance 0: no disparities can be scaled.
      R->degenerate = 1;
      stop = 1;
    }
    else if (R->iter > 1 && (R->old - R->stress) / R->ss <= P->eps) {
      // A relaxed step that makes no more progress hands over to
      // plain steps; only a plain step can stop the start.
      if (R->relax_step) { R->relaxing = 0; }
      else               { R->converged = 1; stop = 1; }
    }
    else if (R->iter >= P->maxit) {
      stop = 1;
    }
    R->old = R->stress;

    if (stop) {
      R->state = ST_FINISH;
    }
    else if (!P->metric) {
      fit_disparities(P, R, sort_threads);
    }
  }

  else { // ST_FINISH

    if (P->metric) {
      // The raw stress of the final configuration, over the sum
      // of squared (scaled) dissimilarities.
      R->stress1 = sqrt(R->stress / P->sumsq);
    }
    else if (R->degenerate) {
      R->stress1 = 1;
    }
    else {
      pava(P, R, sort_threads); // unscaled: stress-1 is scale-free
      each_slice(OP_RESID, P, R, n_threads);
      R->stress1 = sqrt(sum_part(P, R, 0) / R->ss);
    }
    R->state = ST_DONE;
  }
}


// Mode A: several starts at once, one share of run_parallel() each.
// A share may run on any of the pool's threads or on R's own, so
// step() is told n_threads = 0 and sort_threads = 1: it then neither
// calls R (no interrupt checks) nor run_parallel(), which must not be
// nested.
static run_t *g_batch;
static int    g_n_batch;
static int    g_round;

static void *start_worker(void *arg) {
  int ti = ((worker_t *) arg)->i, tn = ((worker_t *) arg)->n;
  for (int r = ti; r < g_n_batch; r += tn)
    for (int it = 0; it < g_round && g_batch[r].state != ST_DONE; it++)
      step(g_P, g_batch + r, 0, 1);
  return NULL;
}



/*======================================================
 * Starting configurations.
 *====================================================*/

// Uniform on [0, 1), with 53 random bits.
static double unif(pcg32_random_t *rng) {
  uint32_t a = pcg32_random_r(rng) >> 5, b = pcg32_random_r(rng) >> 6;
  return (a * 67108864.0 + b) / 9007199254740992.0;
}

// Start number `id`: a given configuration (column-major n x k),
// or uniform random coordinates if `init` is NULL or not finite. A
// given configuration's constant columns (e.g. PCoA axes without a
// positive eigenvalue) would stay constant under the Guttman
// transform, so they get small random values instead.
static void start_config(
    double *x, const double *init, int n, int k, int seed, int id ) {

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, (uint64_t) (uint32_t) seed, (uint64_t) id);

  for (size_t e = 0; init && e < (size_t) n * k; e++)
    if (!R_FINITE(init[e])) init = NULL;

  if (!init) {
    for (int i = 0; i < n; i++)
      for (int c = 0; c < k; c++)
        x[(size_t) i * k + c] = unif(&rng);
    return;
  }

  double ref = 0;
  for (int c = 0; c < k; c++) {
    const double *v = init + (size_t) c * n;
    double mean = 0, ss = 0;
    for (int i = 0; i < n; i++) mean += v[i];
    mean /= n;
    for (int i = 0; i < n; i++) ss += (v[i] - mean) * (v[i] - mean);
    double sd = sqrt(ss / n);
    if (sd > ref) ref = sd;
  }
  double fill = (ref > 0) ? 0.01 * ref : 1;

  for (int c = 0; c < k; c++) {
    const double *v = init + (size_t) c * n;
    int constant = 1;
    for (int i = 1; i < n; i++) if (v[i] != v[0]) { constant = 0; break; }
    for (int i = 0; i < n; i++)
      x[(size_t) i * k + c] = constant ? fill * (unif(&rng) - 0.5) : v[i];
  }
}



/*======================================================
 * Final configuration: centered, rotated to principal
 * axes (Jacobi eigendecomposition of X'X), multiplied by
 * `mult`, or, if rms > 0, scaled so that the root mean
 * square distance between points is rms; and each axis
 * flipped so that its largest absolute coordinate is
 * positive. All sums run in a fixed order.
 *====================================================*/
static void finalize(const double *x, int n, int k, double rms, double mult, double *out) {

  double *z = (double *) R_alloc((size_t) n * k, sizeof(double));
  double *A = (double *) R_alloc((size_t) k * k, sizeof(double));
  double *V = (double *) R_alloc((size_t) k * k, sizeof(double));

  // Center (z is column-major).
  for (int c = 0; c < k; c++) {
    double mean = 0;
    for (int i = 0; i < n; i++) mean += x[(size_t) i * k + c];
    mean /= n;
    for (int i = 0; i < n; i++) z[(size_t) c * n + i] = x[(size_t) i * k + c] - mean;
  }

  // A = Z'Z; V = I.
  double total = 0;
  for (int a = 0; a < k; a++)
    for (int b = 0; b < k; b++) {
      double s = 0;
      for (int i = 0; i < n; i++) s += z[(size_t) a * n + i] * z[(size_t) b * n + i];
      A[a * k + b] = s;
      V[a * k + b] = (a == b);
      if (a == b) total += s;
    }

  // Cyclic Jacobi rotations.
  for (int sweep = 0; sweep < 100; sweep++) {
    double off = 0;
    for (int a = 0; a < k; a++)
      for (int b = a + 1; b < k; b++) off += A[a * k + b] * A[a * k + b];
    if (off <= 1e-30 * total * total) break;
    for (int p = 0; p < k; p++)
      for (int q = p + 1; q < k; q++) {
        double apq = A[p * k + q];
        if (apq == 0) continue;
        double theta = (A[q * k + q] - A[p * k + p]) / (2 * apq);
        double t     = ((theta >= 0) ? 1 : -1) / (fabs(theta) + sqrt(theta * theta + 1));
        double cs    = 1 / sqrt(t * t + 1), sn = t * cs;
        for (int r = 0; r < k; r++) { // columns p, q
          double ap = A[r * k + p], aq = A[r * k + q];
          A[r * k + p] = cs * ap - sn * aq;
          A[r * k + q] = sn * ap + cs * aq;
        }
        for (int r = 0; r < k; r++) { // rows p, q
          double ap = A[p * k + r], aq = A[q * k + r];
          A[p * k + r] = cs * ap - sn * aq;
          A[q * k + r] = sn * ap + cs * aq;
        }
        for (int r = 0; r < k; r++) {
          double vp = V[r * k + p], vq = V[r * k + q];
          V[r * k + p] = cs * vp - sn * vq;
          V[r * k + q] = sn * vp + cs * vq;
        }
      }
  }

  // Order the axes by decreasing variance.
  int *ax = (int *) R_alloc(k, sizeof(int));
  for (int c = 0; c < k; c++) ax[c] = c;
  for (int a = 0; a < k; a++) {
    int best = a;
    for (int b = a + 1; b < k; b++)
      if (A[ax[b] * k + ax[b]] > A[ax[best] * k + ax[best]]) best = b;
    int t = ax[a]; ax[a] = ax[best]; ax[best] = t;
  }

  // Scale: the sum over pairs of squared distances is n * total,
  // and their mean is n * total / (n (n - 1) / 2).
  if (rms > 0) mult *= (total > 0) ? rms / sqrt(2 * total / (n - 1)) : 0;

  for (int c = 0; c < k; c++) {
    double *o = out + (size_t) c * n;
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int a = 0; a < k; a++) s += z[(size_t) a * n + i] * V[a * k + ax[c]];
      o[i] = s * mult;
    }
    int imax = 0;
    for (int i = 1; i < n; i++) if (fabs(o[i]) > fabs(o[imax])) imax = i;
    if (o[imax] < 0) for (int i = 0; i < n; i++) o[i] = -o[i];
  }
}



/*======================================================
 * R interface.
 *
 * sexp_dist:  dm's packed lower triangle (double), with no
 *             NA, negative, or infinite values, and at least
 *             one value > 0.
 * sexp_init:  NULL, or n * k * m0 doubles: m0 given starting
 *             configurations, column-major, in the units of
 *             dm * dscale.
 * sexp_iopts: c(tries, seed, maxit, cpus, detail, primary,
 *             MiB for starts run at once, relax, metric,
 *             pairs of work between interrupt checks (and the
 *             most for running several starts at once),
 *             fewest pairs for threads within a start).
 *             `tries` random starts follow the m0 given ones.
 * sexp_dopts: c(eps, dscale).
 *
 * Returns list(points, stress, best, iters, converged,
 * stress_all, raw_all, iters_all, converged_all) and, if
 * `detail`, also x_all (n x k x m, unrotated, in the units
 * of dm * dscale), trace (maxit x m, NA after a start
 * stops), n_zero, n_groups.
 *====================================================*/
SEXP C_nmds(SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_init, SEXP sexp_iopts, SEXP sexp_dopts) {

  int  n         = asInteger(sexp_n);
  int  k         = asInteger(sexp_k);
  int *iopts     = INTEGER(sexp_iopts);
  int  tries     = iopts[0];
  int  seed      = iopts[1];
  int  maxit     = iopts[2];
  int  n_threads = iopts[3];
  int  detail    = iopts[4];
  int  primary   = iopts[5];
  double par_mem = iopts[6] * 1048576.0;
  int  relax     = iopts[7];
  int  metric    = iopts[8];
  double check   = iopts[9];
  double par_min = iopts[10];
  int  m0        = isNull(sexp_init) ? 0 : (int) (XLENGTH(sexp_init) / ((R_xlen_t) n * k));
  int  m         = m0 + tries;

  if (n < 2 || (!metric && n > NMDS_MAX_N) || k < 1 || k > n - 1 || tries < 0 ||
      maxit < 0 || n_threads < 1 || m < 1 || XLENGTH(sexp_iopts) != 11 ||
      XLENGTH(sexp_dopts) != 2 || check < 1 ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_nmds()."); // # nocov


  //--------------------------------------------------
  // Shared preparation.
  //--------------------------------------------------
  prep_t P;
  memset(&P, 0, sizeof(prep_t));
  P.n       = n;
  P.k       = k;
  P.n_pairs = (size_t) n * (n - 1) / 2;
  P.dist    = REAL(sexp_dist);
  P.maxit   = maxit;
  P.eps     = REAL(sexp_dopts)[0];
  P.dscale  = REAL(sexp_dopts)[1];
  P.primary = primary;
  P.relax   = relax;
  P.metric  = metric;
  P.check   = check;

  size_t N = P.n_pairs;
  int    S = (int) ((N + NMDS_SLICE_MIN - 1) / NMDS_SLICE_MIN);
  if (S > NMDS_SLICES) S = NMDS_SLICES;
  P.n_slices = S;
  for (int s = 0; s <= S; s++) P.slice[s] = (size_t) ((double) N * s / S);
  P.slice[S] = N;

  prep_sort(&P, n_threads);
  P.norm0 = sqrt((double) N / P.sumsq);

  if (!(P.sumsq > 0) || !R_FINITE(P.sumsq))
    error("Invalid distances passed to C_nmds()."); // # nocov


  //--------------------------------------------------
  // How to run the starts: several at once (mode A),
  // while their memory allows, and while one iteration
  // of one is short enough to check for interrupts in
  // time; else one at a time (mode B), with threads for
  // its passes over the pairs when they are long enough.
  //--------------------------------------------------
  double gbytes    = primary ? 12.0 * P.gmax : 0;
  double run_bytes = (metric ? 0 : primary ? 12.0 : 8.0) * N + gbytes +
                     8.0 * S * n * k + 8.0 * n * k + 24.0 * S;
  int    n_batch   = (n_threads < m) ? n_threads : m;
  if (n_batch * run_bytes > par_mem) n_batch = (int) (par_mem / run_bytes);
  if ((double) N > check) n_batch = 1;
  if (n_batch < 2) n_batch = 1;

  int b_threads = ((double) N < par_min) ? 1 : (n_threads < S) ? n_threads : S;

  // Threads to re-sort tie groups (mode B), each with scratch for
  // the largest group, in all at most max(N / 4, gmax) elements.
  int s_threads = (b_threads < NMDS_PARTS) ? b_threads : NMDS_PARTS;
  if (primary && n_batch == 1 && s_threads > 1) {
    size_t cap = (P.gmax > 0) ? (N / 4) / P.gmax : 1;
    if ((size_t) s_threads > cap) s_threads = (cap > 1) ? (int) cap : 1;
    if (s_threads > 1) {
      g_sy = (double *)   R_alloc((size_t) s_threads * P.gmax, sizeof(double));
      g_su = (uint32_t *) R_alloc((size_t) s_threads * P.gmax, sizeof(uint32_t));
    }
  }
  else {
    s_threads = 1;
  }

  int round = (int) (check / (double) N);
  g_round   = (round < 1) ? 1 : (round > NMDS_ROUND) ? NMDS_ROUND : round;

  run_t *runs = (run_t *) R_alloc(n_batch, sizeof(run_t));
  for (int r = 0; r < n_batch; r++) {
    memset(runs + r, 0, sizeof(run_t));
    runs[r].x     = (double *) R_alloc((size_t) n * k, sizeof(double));
    runs[r].y     = metric  ? NULL : (double *) R_alloc(N, sizeof(double));
    runs[r].order = primary ? (uint32_t *) R_alloc(N, sizeof(uint32_t)) : P.order;
    if (primary) {
      size_t g = (P.gmax > 0) ? P.gmax : 1;
      runs[r].sy = (double *)   R_alloc(g, sizeof(double));
      runs[r].su = (uint32_t *) R_alloc(g, sizeof(uint32_t));
    }
    runs[r].acc   = (double *) R_alloc((size_t) S * n * k, sizeof(double));
    runs[r].part  = (double *) R_alloc(3 * (size_t) S, sizeof(double));
    runs[r].trace = detail ? (double *) R_alloc(maxit > 0 ? maxit : 1, sizeof(double)) : NULL;
  }


  //--------------------------------------------------
  // Outputs.
  //--------------------------------------------------
  SEXP sexp_stress_all = PROTECT(allocVector(REALSXP, m));
  SEXP sexp_raw_all    = PROTECT(allocVector(REALSXP, m));
  SEXP sexp_iters_all  = PROTECT(allocVector(INTSXP,  m));
  SEXP sexp_conv_all   = PROTECT(allocVector(LGLSXP,  m));
  SEXP sexp_x_all      = PROTECT(detail ? alloc3DArray(REALSXP, n, k, m) : R_NilValue);
  SEXP sexp_trace      = PROTECT(detail ? allocMatrix(REALSXP, maxit, m) : R_NilValue);

  double *best_x = (double *) R_alloc((size_t) n * k, sizeof(double));
  int     best   = -1;

  g_P = &P;


  //--------------------------------------------------
  // Run the starts, in batches; evaluate them in order.
  //--------------------------------------------------
  for (int first = 0; first < m; ) {

    int nb = m - first;
    if (nb > n_batch) nb = n_batch;

    for (int r = 0; r < nb; r++) {
      int id = first + r;
      start_config(
        runs[r].x, (id < m0) ? REAL(sexp_init) + (size_t) id * n * k : NULL,
        n, k, seed, id );
      runs[r].state      = ST_INIT;
      runs[r].iter       = 0;
      runs[r].converged  = 0;
      runs[r].degenerate = 0;
      runs[r].stress     = NAN;
      runs[r].old        = NAN;
    }

    if (nb > 1) {
      g_batch   = runs;
      g_n_batch = nb;
      for (;;) {
        int done = 1;
        for (int r = 0; r < nb; r++) if (runs[r].state != ST_DONE) done = 0;
        if (done) break;
        run_parallel(start_worker, nb, 100);
        R_CheckUserInterrupt();
      }
    }
    else {
      while (runs[0].state != ST_DONE) {
        step(&P, runs, b_threads, s_threads);
        R_CheckUserInterrupt();
      }
    }

    for (int r = 0; r < nb; r++) {

      int    id = first + r;
      run_t *R  = runs + r;

      REAL(sexp_stress_all)[id]   = R->stress1;
      REAL(sexp_raw_all)[id]      = R->stress;
      INTEGER(sexp_iters_all)[id] = R->iter;
      LOGICAL(sexp_conv_all)[id]  = R->converged;

      if (detail) {
        double *xa = REAL(sexp_x_all) + (size_t) id * n * k;
        for (int i = 0; i < n; i++)
          for (int c = 0; c < k; c++)
            xa[(size_t) c * n + i] = R->x[(size_t) i * k + c];
        double *tr = REAL(sexp_trace) + (size_t) id * maxit;
        for (int it = 0; it < maxit; it++)
          tr[it] = (it < R->iter) ? R->trace[it] : NA_REAL;
      }

      if (best < 0 || R->stress1 < REAL(sexp_stress_all)[best]) {
        best = id;
        memcpy(best_x, R->x, (size_t) n * k * sizeof(double));
      }
    }

    first += nb;
  }


  //--------------------------------------------------
  // The best start, finalized, in dm's units: metric MDS
  // keeps its fitted scale; non-metric MDS gets the root
  // mean square of dm.
  //--------------------------------------------------
  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  finalize(
    best_x, n, k, metric ? 0 : sqrt(P.sumsq / (double) N),
    1 / P.dscale, REAL(sexp_points) );

  size_t n_groups = 0;
  for (size_t w = 0; !metric && w < (N + 63) / 64; w++) {
    uint64_t v = P.ties[w];
    while (v) { v &= v - 1; n_groups++; }
  }

  const char *names[] = {
    "points", "stress", "best", "iters", "converged",
    "stress_all", "raw_all", "iters_all", "converged_all",
    "x_all", "trace", "n_zero", "n_groups", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result,  0, sexp_points);
  SET_VECTOR_ELT(sexp_result,  1, ScalarReal(REAL(sexp_stress_all)[best]));
  SET_VECTOR_ELT(sexp_result,  2, ScalarInteger(best));
  SET_VECTOR_ELT(sexp_result,  3, ScalarInteger(INTEGER(sexp_iters_all)[best]));
  SET_VECTOR_ELT(sexp_result,  4, ScalarLogical(LOGICAL(sexp_conv_all)[best]));
  SET_VECTOR_ELT(sexp_result,  5, sexp_stress_all);
  SET_VECTOR_ELT(sexp_result,  6, sexp_raw_all);
  SET_VECTOR_ELT(sexp_result,  7, sexp_iters_all);
  SET_VECTOR_ELT(sexp_result,  8, sexp_conv_all);
  SET_VECTOR_ELT(sexp_result,  9, sexp_x_all);
  SET_VECTOR_ELT(sexp_result, 10, sexp_trace);
  SET_VECTOR_ELT(sexp_result, 11, ScalarReal((double) P.n_zero));
  SET_VECTOR_ELT(sexp_result, 12, ScalarReal((double) n_groups));

  UNPROTECT(8);
  return sexp_result;
}
