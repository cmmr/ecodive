// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Permutation tests on a distance matrix: the sums of squares
 * behind permanova(), for the observed design and for each
 * permutation of the samples.
 *
 * With A = D^2 (elementwise) and J = I - 11'/n, Gower's
 * centred matrix is G = -1/2 * J * A * J. A model whose
 * centred design spans the orthonormal columns q_1 .. q_m
 * explains
 *
 *     SS = trace(Q' G Q) = sum_c q_c' G q_c
 *
 * of the total sum of squares, trace(G) = sum_{i<j} d_ij^2 / n
 * (McArdle & Anderson 2001). Each q_c sums to zero, so J q_c
 * = q_c, and
 *
 *     q' G q = -1/2 * q' A q = -sum_{i<j} q_i q_j d_ij^2,
 *
 * which is read straight from the packed lower triangle of D,
 * with no n x n matrix.
 *
 * Permutations: vegan permutes G's rows and columns together,
 * G_p = G[p, p]. Then q' G_p q = w' G w with w[p[a]] = q[a],
 * so each permutation scatters the design's rows instead.
 *
 * Partial models (conditions Z, with orthonormal basis Qz):
 * as vegan does, the permuted design columns are made
 * orthogonal to Qz and re-orthonormalized (modified Gram-
 * Schmidt, twice), within each set of columns. A column whose
 * norm falls below PERM_TOL times its original norm is aliased
 * and contributes 0, as in R's qr().
 *
 * One factor, with no conditions: the model's sum of squares
 * is trace(G) - SS_W, with SS_W = sum_g (1/n_g) times the sum
 * of d_ij^2 over pairs within group g (Anderson 2001). That
 * costs one comparison per pair and permutation, instead of
 * g - 1 multiply-adds.
 *
 * Determinism: permutation k (k >= 1) comes from its own PCG
 * stream, seeded with (seed, k), or from a supplied matrix.
 * Permutation 0 is the identity: the observed statistics. The
 * permutations are split into fixed blocks of PERM_BLOCK
 * vectors; each block is computed by one thread, and every
 * sum for one column runs in a fixed order. So results are
 * bit-identical for any number of threads.
 */

#include "ecodive.h"

// Aim for this many vectors (permutations x columns) per
// pass over D.
#define PERM_VECS 32

// Columns of a partial model with less than this fraction
// of their original norm left are aliased.
#define PERM_TOL 1e-7


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the parts of the output and workspace that each task or
// thread owns.
//------------------------------------------------------
static int           pm_n;       // Samples.
static const double *pm_d;       // Packed lower triangle of D.
static int           pm_m;       // Design columns (0 for one factor).
static const double *pm_x;       // n x m, column-major.
static const int    *pm_set;     // Each column's set, for Gram-Schmidt.
static const double *pm_qz;      // n x pz, column-major, or NULL.
static int           pm_pz;
static const int    *pm_grp;     // One factor: each sample's group, or NULL.
static const double *pm_inv_ng;  // ... 1 / group size.
static double        pm_sst;     // trace(G).
static const int    *pm_perms;   // n x n_perm supplied permutations, or NULL.
static const int    *pm_order;   // Samples, ordered by stratum.
static const int    *pm_sbeg;    // Stratum s is pm_order[pm_sbeg[s]] ..
static int           pm_ns;      // ... pm_order[pm_sbeg[s + 1] - 1].
static uint64_t      pm_seed;
static int           pm_bsize;   // Permutations per block.
static int           pm_blk0;    // Blocks for this round:
static int           pm_blk1;    // ... pm_blk0 to pm_blk1 - 1.
static int           pm_total;   // Permutations, including the identity.
static int           pm_mout;    // Output values per permutation.
static double       *pm_out;     // pm_mout x pm_total.
static double       *pm_work;    // Per-thread workspace ...
static size_t        pm_wlen;    // ... of this many doubles.
static int          *pm_iwork;   // Per-thread integer workspace ...
static size_t        pm_iwlen;   // ... of this many ints.


//------------------------------------------------------
// Offset of element (j + 1, j) in a dist vector. Column
// j's elements, rows j + 1 to n - 1, follow it in order.
//------------------------------------------------------
static inline size_t col_start(int j, int n) {
  return (size_t) j * (size_t) n - (size_t) j * (size_t) (j + 1) / 2;
}


//------------------------------------------------------
// Uniform integer in 0 .. bound - 1, without modulo bias.
//------------------------------------------------------
static inline uint32_t bounded_rand(pcg32_random_t *rng, uint32_t bound) {
  uint32_t threshold = (0u - bound) % bound;
  while (1) {
    uint32_t r = pcg32_random_r(rng);
    if (r >= threshold) return r % bound;
  }
}


//------------------------------------------------------
// Permutation k into p: permuted row a is original row
// p[a]. k = 0 is the identity.
//------------------------------------------------------
static void get_perm(int k, int *p) {

  int n = pm_n;

  if (k == 0) {
    for (int a = 0; a < n; a++) p[a] = a;
    return;
  }

  if (pm_perms) {
    memcpy(p, pm_perms + (size_t) (k - 1) * n, (size_t) n * sizeof(int));
    return;
  }

  // Fisher-Yates shuffle within each stratum.
  pcg32_random_t rng;
  pcg32_srandom_r(&rng, pm_seed, (uint64_t) k);

  for (int a = 0; a < n; a++) p[a] = a;

  for (int s = 0; s < pm_ns; s++) {
    const int *idx = pm_order + pm_sbeg[s];
    int        len = pm_sbeg[s + 1] - pm_sbeg[s];
    for (int t = len - 1; t > 0; t--) {
      int u   = (int) bounded_rand(&rng, (uint32_t) t + 1);
      int tmp = p[idx[t]];
      p[idx[t]] = p[idx[u]];
      p[idx[u]] = tmp;
    }
  }
}


//------------------------------------------------------
// Make column v (length n) orthogonal to the pz columns
// of qz and to the nb columns of basis, then normalize
// it. Two passes of modified Gram-Schmidt. Returns 0 (and
// zeroes v) if less than PERM_TOL of its norm is left.
//------------------------------------------------------
static int orthonormalize(double *v, const double *basis, int nb, int n) {

  double norm0 = 0;
  for (int i = 0; i < n; i++) norm0 += v[i] * v[i];
  norm0 = sqrt(norm0);

  for (int pass = 0; pass < 2; pass++) {
    for (int c = 0; c < pm_pz + nb; c++) {
      const double *q = (c < pm_pz)
        ? pm_qz + (size_t) c * n
        : basis + (size_t) (c - pm_pz) * n;
      double dot = 0;
      for (int i = 0; i < n; i++) dot += q[i] * v[i];
      for (int i = 0; i < n; i++) v[i] -= dot * q[i];
    }
  }

  double norm = 0;
  for (int i = 0; i < n; i++) norm += v[i] * v[i];
  norm = sqrt(norm);

  if (!(norm > PERM_TOL * norm0)) {
    memset(v, 0, (size_t) n * sizeof(double));
    return 0;
  }

  for (int i = 0; i < n; i++) v[i] /= norm;
  return 1;
}


//------------------------------------------------------
// General designs. Each task is one block of permutations.
// Their vectors (permutations x columns) go in chunks of
// PERM_W: chunk h is V + h * n * PERM_W, row-major
// (n x PERM_W), zero-padded. One pass over D's triangle
// computes the quadratic forms of every vector in the
// block, each summed in a fixed order.
//------------------------------------------------------

// Repeat F(k) for k = 0 to PERM_W - 1.
#define PERM_W 8
#define PERM_REP(F) F(0) F(1) F(2) F(3) F(4) F(5) F(6) F(7)


//------------------------------------------------------
// The inner loops, for one chunk of PERM_W vectors and one
// column j of D's triangle (nd = n - j - 1 rows below it).
// GCC and clang get explicit two-lane vectors: compilers
// don't vectorize these loops at R's default -O2. Each lane
// does the same operations in the same order as the scalar
// fallback, so results don't depend on which is used.
//
// quad_chunk:  s[k] += vj[k] * sum_t dsq[t] * vi[t, k]
// group_chunk: s[k] += inv_ng[lj[k]] *
//                      sum_t dsq[t] * (li[t, k] == lj[k])
//------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)

typedef double    v2d __attribute__((vector_size(16)));
typedef long long v2l __attribute__((vector_size(16)));

static inline v2d load2(const double *p) { v2d x; memcpy(&x, p, sizeof(x)); return x; }

#define V2_REP(F) F(0) F(1) F(2) F(3)

static inline void quad_chunk(
    const double *dsq, const double *vi, const double *vj, int nd, double *s) {

  #define Q_INIT(h) v2d acc##h = { 0, 0 };
  #define Q_STEP(h) acc##h += a * load2(vi + 2 * h);
  #define Q_DONE(h) s[2 * h] += vj[2 * h] * acc##h[0]; s[2 * h + 1] += vj[2 * h + 1] * acc##h[1];

  V2_REP(Q_INIT)
  for (int t = 0; t < nd; t++, vi += PERM_W) {
    v2d a = { dsq[t], dsq[t] };
    V2_REP(Q_STEP)
  }
  V2_REP(Q_DONE)
}

static inline void group_chunk(
    const double *dsq, const double *li, const double *lj, int nd, double *s,
    const double *inv_ng) {

  #define G_INIT(h) v2d acc##h = { 0, 0 }, g##h = load2(lj + 2 * h);
  #define G_STEP(h) acc##h += (v2d) ((v2l) a & (v2l) (load2(li + 2 * h) == g##h));
  #define G_DONE(h) s[2 * h]     += acc##h[0] * inv_ng[(int) lj[2 * h]]; \
                    s[2 * h + 1] += acc##h[1] * inv_ng[(int) lj[2 * h + 1]];

  V2_REP(G_INIT)
  for (int t = 0; t < nd; t++, li += PERM_W) {
    v2d a = { dsq[t], dsq[t] };
    V2_REP(G_STEP)
  }
  V2_REP(G_DONE)
}

// Two columns, j and j + 1, at once: each row's vectors are
// loaded once for both, and the eight sums are independent,
// so they don't wait on each other. Row j + 1 is column j's
// alone. Every sum runs in the same order as with one column
// at a time.
static inline void quad_chunk2(
    const double *dsq0, const double *dsq1, const double *vi,
    const double *vj0, const double *vj1, int nd, double *s) {

  #define Q2_INIT(h) v2d acc##h = { 0, 0 }, bcc##h = { 0, 0 };
  #define Q2_FIRST(h) acc##h += a * load2(vi + 2 * h);
  #define Q2_STEP(h) { v2d x = load2(vi + 2 * h); acc##h += a * x; bcc##h += b * x; }
  #define Q2_DONE(h) \
    s[2 * h]     += vj0[2 * h]     * acc##h[0]; s[2 * h + 1] += vj0[2 * h + 1] * acc##h[1]; \
    s[2 * h]     += vj1[2 * h]     * bcc##h[0]; s[2 * h + 1] += vj1[2 * h + 1] * bcc##h[1];

  V2_REP(Q2_INIT)
  { v2d a = { dsq0[0], dsq0[0] }; V2_REP(Q2_FIRST) }
  vi += PERM_W;
  for (int t = 1; t < nd; t++, vi += PERM_W) {
    v2d a = { dsq0[t],     dsq0[t]     };
    v2d b = { dsq1[t - 1], dsq1[t - 1] };
    V2_REP(Q2_STEP)
  }
  V2_REP(Q2_DONE)
}

static inline void group_chunk2(
    const double *dsq0, const double *dsq1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s,
    const double *inv_ng) {

  #define G2_INIT(h) v2d acc##h = { 0, 0 }, bcc##h = { 0, 0 }, \
                     g##h = load2(lj0 + 2 * h), k##h = load2(lj1 + 2 * h);
  #define G2_FIRST(h) acc##h += (v2d) ((v2l) a & (v2l) (load2(li + 2 * h) == g##h));
  #define G2_STEP(h) { v2d x = load2(li + 2 * h); \
                       acc##h += (v2d) ((v2l) a & (v2l) (x == g##h)); \
                       bcc##h += (v2d) ((v2l) b & (v2l) (x == k##h)); }
  #define G2_DONE(h) \
    s[2 * h]     += acc##h[0] * inv_ng[(int) lj0[2 * h]];     \
    s[2 * h + 1] += acc##h[1] * inv_ng[(int) lj0[2 * h + 1]]; \
    s[2 * h]     += bcc##h[0] * inv_ng[(int) lj1[2 * h]];     \
    s[2 * h + 1] += bcc##h[1] * inv_ng[(int) lj1[2 * h + 1]];

  V2_REP(G2_INIT)
  { v2d a = { dsq0[0], dsq0[0] }; V2_REP(G2_FIRST) }
  li += PERM_W;
  for (int t = 1; t < nd; t++, li += PERM_W) {
    v2d a = { dsq0[t],     dsq0[t]     };
    v2d b = { dsq1[t - 1], dsq1[t - 1] };
    V2_REP(G2_STEP)
  }
  V2_REP(G2_DONE)
}

#else

static inline void quad_chunk(
    const double *dsq, const double *vi, const double *vj, int nd, double *s) {

  #define Q_INIT(k) double acc##k = 0;
  #define Q_STEP(k) acc##k += a * vi[k];
  #define Q_DONE(k) s[k] += vj[k] * acc##k;

  PERM_REP(Q_INIT)
  for (int t = 0; t < nd; t++, vi += PERM_W) {
    double a = dsq[t];
    PERM_REP(Q_STEP)
  }
  PERM_REP(Q_DONE)
}

static inline void group_chunk(
    const double *dsq, const double *li, const double *lj, int nd, double *s,
    const double *inv_ng) {

  #define G_INIT(k) double acc##k = 0, g##k = lj[k];
  #define G_STEP(k) acc##k += (li[k] == g##k) ? a : 0.0;
  #define G_DONE(k) s[k] += acc##k * inv_ng[(int) g##k];

  PERM_REP(G_INIT)
  for (int t = 0; t < nd; t++, li += PERM_W) {
    double a = dsq[t];
    PERM_REP(G_STEP)
  }
  PERM_REP(G_DONE)
}

static inline void quad_chunk2(
    const double *dsq0, const double *dsq1, const double *vi,
    const double *vj0, const double *vj1, int nd, double *s) {
  quad_chunk(dsq0, vi, vj0, nd, s);
  quad_chunk(dsq1, vi + PERM_W, vj1, nd - 1, s);
}

static inline void group_chunk2(
    const double *dsq0, const double *dsq1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s,
    const double *inv_ng) {
  group_chunk(dsq0, li, lj0, nd, s, inv_ng);
  group_chunk(dsq1, li + PERM_W, lj1, nd - 1, s, inv_ng);
}

#endif

static void *quad_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pm_n;
  int m         = pm_m;
  int max_ch    = (pm_bsize * m + PERM_W - 1) / PERM_W;

  double *V    = pm_work + (size_t) thread_i * pm_wlen;   // max_ch x n x W
  double *s    = V + (size_t) max_ch * n * PERM_W;        // max_ch x W
  double *dsq  = s + (size_t) max_ch * PERM_W;            // 2 x n
  double *cols = dsq + 2 * (size_t) n;                                 // n x m (partial only)
  int    *p    = pm_iwork + (size_t) thread_i * pm_iwlen; // n

  for (int blk = pm_blk0 + thread_i; blk < pm_blk1; blk += n_threads) {

    int k0 = blk * pm_bsize;
    int nb = pm_total - k0;
    if (nb > pm_bsize) nb = pm_bsize;
    int M  = nb * m;
    int nc = (M + PERM_W - 1) / PERM_W;

    if (M % PERM_W)
      memset(V + (size_t) (nc - 1) * n * PERM_W, 0, (size_t) n * PERM_W * sizeof(double));

    //--------------------------------------------------
    // Scatter each permutation's design rows into V.
    // Vector v = b * m + c is lane v % W of chunk v / W.
    //--------------------------------------------------
    for (int b = 0; b < nb; b++) {

      get_perm(k0 + b, p);

      if (pm_qz) {

        // Partial model: permute, then orthonormalize
        // against Qz and the set's earlier columns.
        int set_start = 0;
        for (int c = 0; c < m; c++) {
          if (c > 0 && pm_set[c] != pm_set[c - 1]) set_start = c;
          double       *v = cols + (size_t) c * n;
          const double *x = pm_x + (size_t) c * n;
          for (int a = 0; a < n; a++) v[p[a]] = x[a];
          // Aliased columns are zeroed and still serve
          // (harmlessly) as basis vectors.
          orthonormalize(v, cols + (size_t) set_start * n, c - set_start, n);
        }
        for (int c = 0; c < m; c++) {
          int           v  = b * m + c;
          double       *vc = V + (size_t) (v / PERM_W) * n * PERM_W + v % PERM_W;
          const double *y  = cols + (size_t) c * n;
          for (int i = 0; i < n; i++) vc[(size_t) i * PERM_W] = y[i];
        }
      }
      else {
        for (int c = 0; c < m; c++) {
          int           v  = b * m + c;
          double       *vc = V + (size_t) (v / PERM_W) * n * PERM_W + v % PERM_W;
          const double *x  = pm_x + (size_t) c * n;
          for (int a = 0; a < n; a++) vc[(size_t) p[a] * PERM_W] = x[a];
        }
      }
    }

    //--------------------------------------------------
    // s[v] = sum_j V[j, v] * sum_{i > j} d_ij^2 V[i, v]
    //--------------------------------------------------
    for (int v = 0; v < nc * PERM_W; v++) s[v] = 0;

    for (int j = 0; j < n - 1; ) {

      int nd = n - j - 1;
      const double *d0 = pm_d + col_start(j, n);
      for (int t = 0; t < nd; t++) dsq[t] = d0[t] * d0[t];

      if (nd >= 2) {

        // Columns j and j + 1.
        const double *d1 = pm_d + col_start(j + 1, n);
        for (int t = 0; t < nd - 1; t++) dsq[n + t] = d1[t] * d1[t];

        for (int h = 0; h < nc; h++) {
          const double *vc = V + (size_t) h * n * PERM_W;
          quad_chunk2(
            dsq, dsq + n, vc + (size_t) (j + 1) * PERM_W,
            vc + (size_t) j * PERM_W, vc + (size_t) (j + 1) * PERM_W,
            nd, s + h * PERM_W );
        }
        j += 2;
      }
      else {
        for (int h = 0; h < nc; h++) {
          const double *vc = V + (size_t) h * n * PERM_W;
          quad_chunk(
            dsq, vc + (size_t) (j + 1) * PERM_W, vc + (size_t) j * PERM_W,
            nd, s + h * PERM_W );
        }
        j += 1;
      }
    }

    for (int b = 0; b < nb; b++)
      for (int c = 0; c < m; c++)
        pm_out[(size_t) (k0 + b) * pm_mout + c] = -s[b * m + c];
  }

  return NULL;
}


//------------------------------------------------------
// One factor. Each task is one block of PERM_FW
// permutations. L holds each sample's group (as a double,
// so comparisons vectorize alongside the sums) under each
// permutation: chunk h is L + h * n * PERM_W, row-major
// (n x PERM_W). Unused lanes are group 0, and ignored.
//------------------------------------------------------
#define PERM_FW (4 * PERM_W)

static void *factor_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pm_n;
  int nc        = PERM_FW / PERM_W;

  double *L   = pm_work + (size_t) thread_i * pm_wlen;   // nc x n x W
  double *s   = L + (size_t) n * PERM_FW;                // PERM_FW
  double *dsq = s + PERM_FW;                             // 2 x n
  int    *p   = pm_iwork + (size_t) thread_i * pm_iwlen; // n

  for (int blk = pm_blk0 + thread_i; blk < pm_blk1; blk += n_threads) {

    int k0 = blk * PERM_FW;
    int nb = pm_total - k0;
    if (nb > PERM_FW) nb = PERM_FW;

    memset(L, 0, (size_t) n * PERM_FW * sizeof(double));

    for (int b = 0; b < nb; b++) {
      double *lc = L + (size_t) (b / PERM_W) * n * PERM_W + b % PERM_W;
      get_perm(k0 + b, p);
      for (int a = 0; a < n; a++)
        lc[(size_t) p[a] * PERM_W] = pm_grp[a];
    }

    for (int b = 0; b < PERM_FW; b++) s[b] = 0;

    for (int j = 0; j < n - 1; ) {

      int nd = n - j - 1;
      const double *d0 = pm_d + col_start(j, n);
      for (int t = 0; t < nd; t++) dsq[t] = d0[t] * d0[t];

      if (nd >= 2) {

        // Columns j and j + 1.
        const double *d1 = pm_d + col_start(j + 1, n);
        for (int t = 0; t < nd - 1; t++) dsq[n + t] = d1[t] * d1[t];

        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * PERM_W;
          group_chunk2(
            dsq, dsq + n, vc + (size_t) (j + 1) * PERM_W,
            vc + (size_t) j * PERM_W, vc + (size_t) (j + 1) * PERM_W,
            nd, s + h * PERM_W, pm_inv_ng );
        }
        j += 2;
      }
      else {
        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * PERM_W;
          group_chunk(
            dsq, vc + (size_t) (j + 1) * PERM_W, vc + (size_t) j * PERM_W,
            nd, s + h * PERM_W, pm_inv_ng );
        }
        j += 1;
      }
    }

    for (int b = 0; b < nb; b++)
      pm_out[k0 + b] = pm_sst - s[b];
  }

  return NULL;
}


//------------------------------------------------------
// Sum of squares of each column of the lower triangle,
// and the smallest distance, per column.
//------------------------------------------------------
static double *ss_col, *min_col;

static void *sum_squares(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pm_n;

  for (int j = thread_i; j < n - 1; j += n_threads) {
    const double *d  = pm_d + col_start(j, n);
    int           nd = n - j - 1;
    double        ss = 0, mn = INFINITY;
    for (int i = 0; i < nd; i++) {
      ss += d[i] * d[i];
      if (d[i] < mn) mn = d[i];
    }
    ss_col[j]  = ss;
    min_col[j] = mn;
  }

  return NULL;
}


//======================================================
// R interface.
//
// dm:     packed lower triangle of D (double).
// n:      number of samples.
// x:      n x m design columns (double), or NULL. Without
//         qz, each set of columns must be orthonormal and
//         centred; with qz, centred.
// set:    integer set id of each column (non-decreasing).
// qz:     n x pz orthonormal, centred basis of the
//         conditions, or NULL.
// grp:    one factor instead of x: each sample's group,
//         0 .. g - 1, or NULL.
// perms:  n x n_perm integer matrix of 0-based permutations,
//         or NULL for random ones.
// strata: each sample's stratum, 0 .. n_strata - 1, or NULL.
// opts:   c(n_perm, seed, n_threads).
//
// Returns list(ss, total): ss is m x (n_perm + 1) (1 row for
// one factor), column 1 for the identity; total = trace(G).
//======================================================
SEXP C_permanova(
    SEXP sexp_dm, SEXP sexp_n, SEXP sexp_x, SEXP sexp_set, SEXP sexp_qz,
    SEXP sexp_grp, SEXP sexp_perms, SEXP sexp_strata, SEXP sexp_opts ) {

  int n         = asInteger(sexp_n);
  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];

  if (n < 2 || n_perm < 0 || n_threads < 1 ||
      XLENGTH(sexp_dm) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_permanova()."); // # nocov

  pm_n     = n;
  pm_d     = REAL(sexp_dm);
  pm_seed  = (uint64_t) seed;
  pm_total = n_perm + 1;


  //--------------------------------------------------
  // trace(G), and a check for bad distances.
  //--------------------------------------------------
  ss_col  = (double *) R_alloc(n, sizeof(double));
  min_col = (double *) R_alloc(n, sizeof(double));
  run_parallel(sum_squares, n_threads, n);

  double total = 0, dmin = INFINITY;
  for (int j = 0; j < n - 1; j++) {
    total += ss_col[j];
    if (min_col[j] < dmin) dmin = min_col[j];
  }
  if (!R_FINITE(total))
    error("`dm` must not contain NA, NaN, or infinite values.");
  if (dmin < -1e-7)
    error("`dm` must not contain negative distances.");
  pm_sst = total / n;


  //--------------------------------------------------
  // The design: one factor, or columns.
  //--------------------------------------------------
  pm_grp = NULL; pm_inv_ng = NULL;
  pm_x   = NULL; pm_set    = NULL; pm_m = 0;
  pm_qz  = NULL; pm_pz     = 0;

  if (!isNull(sexp_grp)) {

    pm_grp = INTEGER(sexp_grp);
    int g  = 0;
    for (int a = 0; a < n; a++)
      if (pm_grp[a] + 1 > g) g = pm_grp[a] + 1;

    double *inv = (double *) R_alloc(g, sizeof(double));
    for (int c = 0; c < g; c++) inv[c] = 0;
    for (int a = 0; a < n; a++) inv[pm_grp[a]] += 1;
    for (int c = 0; c < g; c++) inv[c] = (inv[c] > 0) ? 1 / inv[c] : 0;
    pm_inv_ng = inv;
    pm_mout   = 1;
  }
  else {
    pm_x    = REAL(sexp_x);
    pm_m    = ncols(sexp_x);
    pm_set  = INTEGER(sexp_set);
    pm_mout = pm_m;
    if (!isNull(sexp_qz)) {
      pm_qz = REAL(sexp_qz);
      pm_pz = ncols(sexp_qz);
    }
  }


  //--------------------------------------------------
  // Permutations: supplied, or shuffles within strata.
  //--------------------------------------------------
  pm_perms = isNull(sexp_perms) ? NULL : INTEGER(sexp_perms);

  int *order = (int *) R_alloc(n, sizeof(int));
  int  ns    = 1;
  int *sbeg;

  if (isNull(sexp_strata)) {
    for (int a = 0; a < n; a++) order[a] = a;
    sbeg    = (int *) R_alloc(2, sizeof(int));
    sbeg[0] = 0;
    sbeg[1] = n;
  }
  else {
    const int *strata = INTEGER(sexp_strata);
    ns = 0;
    for (int a = 0; a < n; a++)
      if (strata[a] + 1 > ns) ns = strata[a] + 1;
    sbeg = (int *) R_alloc(ns + 1, sizeof(int));
    for (int s = 0; s <= ns; s++) sbeg[s] = 0;
    for (int a = 0; a < n; a++) sbeg[strata[a] + 1]++;
    for (int s = 0; s < ns; s++) sbeg[s + 1] += sbeg[s];
    int *fill = (int *) R_alloc(ns, sizeof(int));
    memcpy(fill, sbeg, (size_t) ns * sizeof(int));
    for (int a = 0; a < n; a++) order[fill[strata[a]]++] = a;
  }
  pm_order = order;
  pm_sbeg  = sbeg;
  pm_ns    = ns;


  //--------------------------------------------------
  // Blocks, and workspace for each thread.
  //--------------------------------------------------
  int m    = (pm_m > 0) ? pm_m : 1;
  pm_bsize = PERM_VECS / m;
  if (pm_bsize < 1)  pm_bsize = 1;
  if (pm_grp)        pm_bsize = PERM_FW;

  int n_blocks = (pm_total + pm_bsize - 1) / pm_bsize;
  if (n_threads > n_blocks) n_threads = n_blocks;

  pm_iwlen = (size_t) n;
  if (pm_grp) {
    pm_wlen = (size_t) n * PERM_FW + PERM_FW + 2 * (size_t) n;
  }
  else {
    size_t W = (size_t) ((pm_bsize * pm_m + PERM_W - 1) / PERM_W) * PERM_W;
    pm_wlen  = (size_t) n * W + W + 2 * (size_t) n + (pm_qz ? (size_t) n * pm_m : 0);
  }
  pm_work  = (double *) R_alloc(pm_wlen  * n_threads, sizeof(double));
  pm_iwork = (int *)    R_alloc(pm_iwlen * n_threads, sizeof(int));

  SEXP sexp_ss = PROTECT(allocMatrix(REALSXP, pm_mout, pm_total));
  pm_out = REAL(sexp_ss);


  //--------------------------------------------------
  // Rounds of blocks, checking for interrupts between
  // them. Each round has at least ~1e8 multiply-adds
  // per thread, or every remaining block.
  //--------------------------------------------------
  double pairs     = (double) n * (n - 1) / 2;
  double blk_work  = pairs * pm_bsize * (pm_grp ? 1 : pm_m) + 1;
  int    per_round = (int) (1e8 * n_threads / blk_work);
  if (per_round < n_threads) per_round = n_threads;

  // run_parallel() starts threads only when told of at
  // least 100 tasks.
  int tasks_hint = (blk_work > 1e5) ? 100 : 1;

  for (pm_blk0 = 0; pm_blk0 < n_blocks; pm_blk0 = pm_blk1) {

    pm_blk1 = pm_blk0 + per_round;
    if (pm_blk1 > n_blocks) pm_blk1 = n_blocks;

    int nt = (pm_blk1 - pm_blk0 < n_threads) ? pm_blk1 - pm_blk0 : n_threads;
    run_parallel(pm_grp ? factor_worker : quad_worker, nt,
                 (pm_blk1 - pm_blk0 > 1) ? tasks_hint : 1);

    if (pm_blk1 < n_blocks) R_CheckUserInterrupt();
  }

  const char *names[] = { "ss", "total", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_ss);
  SET_VECTOR_ELT(sexp_result, 1, ScalarReal(pm_sst));

  UNPROTECT(2);
  return sexp_result;
}
