// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Multi-response permutation procedure (MRPP): the weighted
 * mean within-group distance, delta, for the observed groups
 * and for each permutation of them.
 *
 * Group g, with n_g samples, has P_g = n_g (n_g - 1) / 2
 * pairs, and mean within-group distance
 *
 *     delta_g = S_g / P_g,   S_g = sum of d_ij over its pairs.
 *
 * With weights w_g (summing to W over the groups with at
 * least two samples),
 *
 *     delta = sum_g (w_g / W) delta_g = sum_g c_g S_g,
 *
 * where c_g = w_g / (W P_g) is computed in R. So each
 * permutation costs one comparison per pair: add d_ij to
 * pair (i, j)'s sum when i and j are in the same group, and
 * scale each column's sum by its group's c_g. That is the
 * kernel of permanova()'s one-factor path, with d_ij in place
 * of d_ij^2.
 *
 * Permutations: as vegan::mrpp(), which takes grouping[p],
 * permuted sample a is in group g[p[a]].
 *
 * Determinism: permutation k (k >= 1) comes from its own PCG
 * stream, seeded with (seed, k), or from a supplied matrix.
 * Permutation 0 is the identity: the observed delta. The
 * permutations are split into fixed blocks of MRPP_FW; each
 * block is computed by one thread, and every sum runs in a
 * fixed order. So results are bit-identical for any number
 * of threads.
 *
 * TODO: use shared helpers after Batch C. The permutation
 * generator, bounded_rand(), col_start(), and the group
 * kernels are copies of permanova.c's; share them (and the
 * thread pool / parallel_for) once Batch C is merged.
 */

#include "ecodive.h"

#define MRPP_W  8              // Lanes per chunk.
#define MRPP_FW (4 * MRPP_W)   // Permutations per block.


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the parts of the output that each block owns.
//------------------------------------------------------
static int           mr_n;       // Samples.
static const double *mr_d;       // Packed lower triangle of D.
static const int    *mr_grp;     // Each sample's group, 0 .. g - 1.
static const double *mr_coef;    // c_g for each group.
static const int    *mr_perms;   // n x n_perm supplied permutations, or NULL.
static const int    *mr_order;   // Samples, ordered by stratum.
static const int    *mr_sbeg;    // Stratum s is mr_order[mr_sbeg[s]] ..
static int           mr_ns;      // ... mr_order[mr_sbeg[s + 1] - 1].
static uint64_t      mr_seed;
static int           mr_blk0;    // Blocks for this round:
static int           mr_blk1;    // ... mr_blk0 to mr_blk1 - 1.
static int           mr_total;   // Permutations, including the identity.
static double       *mr_out;     // mr_total deltas.
static double       *mr_work;    // Per-thread workspace ...
static size_t        mr_wlen;    // ... of this many doubles.
static int          *mr_iwork;   // Per-thread integer workspace (n each).


//------------------------------------------------------
// Offset of element (j + 1, j) in a dist vector. Column
// j's elements, rows j + 1 to n - 1, follow it in order.
// TODO: use shared helper after Batch C.
//------------------------------------------------------
static inline size_t mr_col_start(int j, int n) {
  return (size_t) j * (size_t) n - (size_t) j * (size_t) (j + 1) / 2;
}


//------------------------------------------------------
// Uniform integer in 0 .. bound - 1, without modulo bias.
// TODO: use shared helper after Batch C.
//------------------------------------------------------
static inline uint32_t mr_bounded_rand(pcg32_random_t *rng, uint32_t bound) {
  uint32_t threshold = (0u - bound) % bound;
  while (1) {
    uint32_t r = pcg32_random_r(rng);
    if (r >= threshold) return r % bound;
  }
}


//------------------------------------------------------
// Permutation k into p: permuted sample a takes original
// sample p[a]'s group. k = 0 is the identity.
// TODO: use shared helper after Batch C.
//------------------------------------------------------
static void mr_get_perm(int k, int *p) {

  int n = mr_n;

  if (k == 0) {
    for (int a = 0; a < n; a++) p[a] = a;
    return;
  }

  if (mr_perms) {
    memcpy(p, mr_perms + (size_t) (k - 1) * n, (size_t) n * sizeof(int));
    return;
  }

  // Fisher-Yates shuffle within each stratum.
  pcg32_random_t rng;
  pcg32_srandom_r(&rng, mr_seed, (uint64_t) k);

  for (int a = 0; a < n; a++) p[a] = a;

  for (int s = 0; s < mr_ns; s++) {
    const int *idx = mr_order + mr_sbeg[s];
    int        len = mr_sbeg[s + 1] - mr_sbeg[s];
    for (int t = len - 1; t > 0; t--) {
      int u   = (int) mr_bounded_rand(&rng, (uint32_t) t + 1);
      int tmp = p[idx[t]];
      p[idx[t]] = p[idx[u]];
      p[idx[u]] = tmp;
    }
  }
}


//------------------------------------------------------
// The inner loops, for one chunk of MRPP_W permutations
// and one column j of D's triangle (nd = n - j - 1 rows
// below it):
//
//   s[k] += coef[lj[k]] * sum_t d[t] * (li[t, k] == lj[k])
//
// GCC and clang get explicit two-lane vectors: compilers
// don't vectorize these loops at R's default -O2. Each lane
// does the same operations in the same order as the scalar
// fallback, so results don't depend on which is used.
// TODO: share with permanova.c after Batch C.
//------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)

typedef double    mr_v2d __attribute__((vector_size(16)));
typedef long long mr_v2l __attribute__((vector_size(16)));

static inline mr_v2d mr_load2(const double *p) { mr_v2d x; memcpy(&x, p, sizeof(x)); return x; }

#define MR_REP(F) F(0) F(1) F(2) F(3)

static inline void mr_chunk(
    const double *d, const double *li, const double *lj, int nd, double *s,
    const double *coef) {

  #define G_INIT(h) mr_v2d acc##h = { 0, 0 }, g##h = mr_load2(lj + 2 * h);
  #define G_STEP(h) acc##h += (mr_v2d) ((mr_v2l) a & (mr_v2l) (mr_load2(li + 2 * h) == g##h));
  #define G_DONE(h) s[2 * h]     += acc##h[0] * coef[(int) lj[2 * h]]; \
                    s[2 * h + 1] += acc##h[1] * coef[(int) lj[2 * h + 1]];

  MR_REP(G_INIT)
  for (int t = 0; t < nd; t++, li += MRPP_W) {
    mr_v2d a = { d[t], d[t] };
    MR_REP(G_STEP)
  }
  MR_REP(G_DONE)
}

// Two columns, j and j + 1, at once: each row's labels are
// loaded once for both. Row j + 1 is column j's alone.
// Every sum runs in the same order as with one column at a
// time.
static inline void mr_chunk2(
    const double *d0, const double *d1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s,
    const double *coef) {

  #define G2_INIT(h) mr_v2d acc##h = { 0, 0 }, bcc##h = { 0, 0 }, \
                     g##h = mr_load2(lj0 + 2 * h), k##h = mr_load2(lj1 + 2 * h);
  #define G2_FIRST(h) acc##h += (mr_v2d) ((mr_v2l) a & (mr_v2l) (mr_load2(li + 2 * h) == g##h));
  #define G2_STEP(h) { mr_v2d x = mr_load2(li + 2 * h); \
                       acc##h += (mr_v2d) ((mr_v2l) a & (mr_v2l) (x == g##h)); \
                       bcc##h += (mr_v2d) ((mr_v2l) b & (mr_v2l) (x == k##h)); }
  #define G2_DONE(h) \
    s[2 * h]     += acc##h[0] * coef[(int) lj0[2 * h]];     \
    s[2 * h + 1] += acc##h[1] * coef[(int) lj0[2 * h + 1]]; \
    s[2 * h]     += bcc##h[0] * coef[(int) lj1[2 * h]];     \
    s[2 * h + 1] += bcc##h[1] * coef[(int) lj1[2 * h + 1]];

  MR_REP(G2_INIT)
  { mr_v2d a = { d0[0], d0[0] }; MR_REP(G2_FIRST) }
  li += MRPP_W;
  for (int t = 1; t < nd; t++, li += MRPP_W) {
    mr_v2d a = { d0[t],     d0[t]     };
    mr_v2d b = { d1[t - 1], d1[t - 1] };
    MR_REP(G2_STEP)
  }
  MR_REP(G2_DONE)
}

#else

#define MR_REP8(F) F(0) F(1) F(2) F(3) F(4) F(5) F(6) F(7)

static inline void mr_chunk(
    const double *d, const double *li, const double *lj, int nd, double *s,
    const double *coef) {

  #define G_INIT(k) double acc##k = 0, g##k = lj[k];
  #define G_STEP(k) acc##k += (li[k] == g##k) ? a : 0.0;
  #define G_DONE(k) s[k] += acc##k * coef[(int) g##k];

  MR_REP8(G_INIT)
  for (int t = 0; t < nd; t++, li += MRPP_W) {
    double a = d[t];
    MR_REP8(G_STEP)
  }
  MR_REP8(G_DONE)
}

static inline void mr_chunk2(
    const double *d0, const double *d1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s,
    const double *coef) {
  mr_chunk(d0, li, lj0, nd, s, coef);
  mr_chunk(d1, li + MRPP_W, lj1, nd - 1, s, coef);
}

#endif


//------------------------------------------------------
// Each task is one block of MRPP_FW permutations. L holds
// each sample's group (as a double, so comparisons
// vectorize alongside the sums) under each permutation:
// chunk h is L + h * n * MRPP_W, row-major (n x MRPP_W).
// Unused lanes are group 0, and ignored.
//------------------------------------------------------
static void *mrpp_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = mr_n;
  int nc        = MRPP_FW / MRPP_W;

  double *L = mr_work + (size_t) thread_i * mr_wlen;     // nc x n x W
  double *s = L + (size_t) n * MRPP_FW;                  // MRPP_FW
  int    *p = mr_iwork + (size_t) thread_i * n;          // n

  for (int blk = mr_blk0 + thread_i; blk < mr_blk1; blk += n_threads) {

    int k0 = blk * MRPP_FW;
    int nb = mr_total - k0;
    if (nb > MRPP_FW) nb = MRPP_FW;

    memset(L, 0, (size_t) n * MRPP_FW * sizeof(double));

    for (int b = 0; b < nb; b++) {
      double *lc = L + (size_t) (b / MRPP_W) * n * MRPP_W + b % MRPP_W;
      mr_get_perm(k0 + b, p);
      for (int a = 0; a < n; a++)
        lc[(size_t) a * MRPP_W] = mr_grp[p[a]];
    }

    for (int b = 0; b < MRPP_FW; b++) s[b] = 0;

    for (int j = 0; j < n - 1; ) {

      int nd = n - j - 1;
      const double *d0 = mr_d + mr_col_start(j, n);

      if (nd >= 2) {
        const double *d1 = mr_d + mr_col_start(j + 1, n);
        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * MRPP_W;
          mr_chunk2(
            d0, d1, vc + (size_t) (j + 1) * MRPP_W,
            vc + (size_t) j * MRPP_W, vc + (size_t) (j + 1) * MRPP_W,
            nd, s + h * MRPP_W, mr_coef );
        }
        j += 2;
      }
      else {
        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * MRPP_W;
          mr_chunk(
            d0, vc + (size_t) (j + 1) * MRPP_W, vc + (size_t) j * MRPP_W,
            nd, s + h * MRPP_W, mr_coef );
        }
        j += 1;
      }
    }

    for (int b = 0; b < nb; b++)
      mr_out[k0 + b] = s[b];
  }

  return NULL;
}


//======================================================
// R interface.
//
// dm:     packed lower triangle of D (double), checked.
// n:      number of samples.
// grp:    each sample's group, 0 .. g - 1 (integer).
// coef:   c_g = w_g / (W P_g) for each group (double); 0
//         for groups left out of delta.
// perms:  n x n_perm integer matrix of 0-based permutations,
//         or NULL for random ones.
// strata: each sample's stratum, 0 .. n_strata - 1, or NULL.
// opts:   c(n_perm, seed, n_threads).
//
// Returns list(delta, class_sum): delta has n_perm + 1
// values, the observed one first; class_sum is S_g for the
// observed groups.
//======================================================
SEXP C_mrpp(
    SEXP sexp_dm, SEXP sexp_n, SEXP sexp_grp, SEXP sexp_coef,
    SEXP sexp_perms, SEXP sexp_strata, SEXP sexp_opts ) {

  int n         = asInteger(sexp_n);
  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];
  int g         = LENGTH(sexp_coef);

  if (n < 2 || n_perm < 0 || n_threads < 1 || g < 1 ||
      XLENGTH(sexp_dm) != (R_xlen_t) n * (n - 1) / 2 ||
      LENGTH(sexp_grp) != n)
    error("Invalid arguments to C_mrpp()."); // # nocov

  const int *grp = INTEGER(sexp_grp);
  for (int a = 0; a < n; a++)
    if (grp[a] < 0 || grp[a] >= g)
      error("Invalid arguments to C_mrpp()."); // # nocov

  mr_n     = n;
  mr_d     = REAL(sexp_dm);
  mr_grp   = grp;
  mr_coef  = REAL(sexp_coef);
  mr_seed  = (uint64_t) seed;
  mr_total = n_perm + 1;
  mr_perms = isNull(sexp_perms) ? NULL : INTEGER(sexp_perms);


  //--------------------------------------------------
  // Strata: samples ordered by stratum.
  //--------------------------------------------------
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
  mr_order = order;
  mr_sbeg  = sbeg;
  mr_ns    = ns;


  //--------------------------------------------------
  // Observed within-group sums, S_g.
  //--------------------------------------------------
  SEXP sexp_cs = PROTECT(allocVector(REALSXP, g));
  double *cs = REAL(sexp_cs);
  for (int c = 0; c < g; c++) cs[c] = 0;
  for (int j = 0; j < n - 1; j++) {
    const double *d  = mr_d + mr_col_start(j, n);
    int           gj = grp[j];
    for (int i = j + 1; i < n; i++)
      if (grp[i] == gj) cs[gj] += d[i - j - 1];
  }


  //--------------------------------------------------
  // Blocks, and workspace for each thread.
  //--------------------------------------------------
  int n_blocks = (mr_total + MRPP_FW - 1) / MRPP_FW;
  if (n_threads > n_blocks) n_threads = n_blocks;

  mr_wlen  = (size_t) n * MRPP_FW + MRPP_FW;
  mr_work  = (double *) R_alloc(mr_wlen * n_threads, sizeof(double));
  mr_iwork = (int *)    R_alloc((size_t) n * n_threads, sizeof(int));

  SEXP sexp_delta = PROTECT(allocVector(REALSXP, mr_total));
  mr_out = REAL(sexp_delta);


  //--------------------------------------------------
  // Rounds of blocks, checking for interrupts between
  // them. Each round has at least ~1e8 comparisons per
  // thread, or every remaining block.
  // TODO: use the thread pool / parallel_for after
  // Batch C.
  //--------------------------------------------------
  double pairs     = (double) n * (n - 1) / 2;
  double blk_work  = pairs * MRPP_FW + 1;
  int    per_round = (int) (1e8 * n_threads / blk_work);
  if (per_round < n_threads) per_round = n_threads;

  // run_parallel() starts threads only when told of at
  // least 100 tasks.
  int tasks_hint = (blk_work > 1e5) ? 100 : 1;

  for (mr_blk0 = 0; mr_blk0 < n_blocks; mr_blk0 = mr_blk1) {

    mr_blk1 = mr_blk0 + per_round;
    if (mr_blk1 > n_blocks) mr_blk1 = n_blocks;

    int nt = (mr_blk1 - mr_blk0 < n_threads) ? mr_blk1 - mr_blk0 : n_threads;
    run_parallel(mrpp_worker, nt, (mr_blk1 - mr_blk0 > 1) ? tasks_hint : 1);

    if (mr_blk1 < n_blocks) R_CheckUserInterrupt();
  }

  const char *names[] = { "delta", "class_sum", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_delta);
  SET_VECTOR_ELT(sexp_result, 1, sexp_cs);

  UNPROTECT(3);
  return sexp_result;
}
