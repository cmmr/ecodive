// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Analysis of similarities (ANOSIM; Clarke 1993): the rank
 * sums behind anosim(), for the observed grouping and for
 * each permutation of the samples.
 *
 * The M = n(n-1)/2 distances are ranked once, with tied
 * distances given their average rank, as R's rank() does.
 * Ranks are stored doubled, so that they and every sum of
 * them are whole numbers, exact in a double up to 2^53
 * (any n below about 13,000).
 *
 * The number of pairs within groups, N_w, is the same for
 * every permutation of the group labels, and the ranks sum
 * to M(M + 1)/2. So each permutation needs only the sum of
 * the ranks of the pairs within groups, S_w; the R side
 * turns it into
 *
 *     R = (r_B - r_W) / (M / 2),
 *     r_W = S_w / N_w,
 *     r_B = (M(M + 1)/2 - S_w) / (M - N_w).
 *
 * Each permutation reads the packed triangle of ranks once,
 * shared by a block of PERM_FW permutations, with no n x n
 * matrix.
 *
 * Permutations come from permute.c, shared with ecodive's
 * other permutation tests: permutation k (k >= 1) comes from
 * its own PCG stream, seeded with (seed, k), or from a
 * supplied matrix; permutation 0 is the identity, the
 * observed statistic. As in vegan::anosim(), sample a gets
 * the group of sample p[a] (vegan's grouping[take]).
 *
 * Determinism: each block of permutations is computed by one
 * thread, and every sum runs in a fixed order, so results
 * are bit-identical for any number of threads.
 */

#include "ecodive.h"


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the parts of the output and workspace that each task or
// thread owns.
//------------------------------------------------------
static int           an_n;       // Samples.
static const double *an_rk;      // Packed lower triangle of doubled ranks.
static const int    *an_grp;     // Each sample's group, 0 .. g - 1.
static perm_source_t an_ps;      // The permutations.
static int           an_blk0;    // Blocks for this round:
static int           an_blk1;    // ... an_blk0 to an_blk1 - 1.
static int           an_total;   // Permutations, including the identity.
static double       *an_out;     // Doubled S_w, one per permutation.
static double       *an_work;    // Per-thread workspace ...
static size_t        an_wlen;    // ... of this many doubles.
static int          *an_iwork;   // Per-thread integer workspace ...
static size_t        an_iwlen;   // ... of this many ints.


//------------------------------------------------------
// Doubled average ranks of the M values in d, into rk:
// a run of equal values in sorted positions a .. b
// (1-based) all get a + b.
//
// LSD radix sort of the values' bit patterns, 11 bits
// per pass, skipping digits that every value shares.
// Equal values have equal keys, so the order within ties
// doesn't matter. O(M) time; 24 bytes per value of
// scratch space, released before returning.
//------------------------------------------------------
#define RDX_BITS   11
#define RDX_SIZE   (1 << RDX_BITS)
#define RDX_PASSES 6   // ceil(64 / RDX_BITS)

static void rank_values(const double *d, size_t M, double *rk) {

  if (M == 0) return;

  const void *vmax = vmaxget();

  uint64_t *k0 = (uint64_t *) R_alloc(M, sizeof(uint64_t));
  uint64_t *k1 = (uint64_t *) R_alloc(M, sizeof(uint64_t));
  uint32_t *i0 = (uint32_t *) R_alloc(M, sizeof(uint32_t));
  uint32_t *i1 = (uint32_t *) R_alloc(M, sizeof(uint32_t));
  size_t   *h  = (size_t *)   R_alloc(RDX_PASSES * RDX_SIZE, sizeof(size_t));

  memset(h, 0, RDX_PASSES * RDX_SIZE * sizeof(size_t));

  // Keys that sort as the values do: flip every bit of
  // negatives, and only the sign bit of the rest. Adding
  // 0 turns -0 into +0.
  for (size_t t = 0; t < M; t++) {
    double   x = d[t] + 0.0;
    uint64_t u;
    memcpy(&u, &x, sizeof(u));
    u ^= (u >> 63) ? ~(uint64_t) 0 : ((uint64_t) 1 << 63);
    k0[t] = u;
    i0[t] = (uint32_t) t;
    for (int p = 0; p < RDX_PASSES; p++)
      h[p * RDX_SIZE + ((u >> (p * RDX_BITS)) & (RDX_SIZE - 1))]++;
  }

  for (int p = 0; p < RDX_PASSES; p++) {

    size_t *hp    = h + p * RDX_SIZE;
    int     shift = p * RDX_BITS;

    // Every key has the same digit: nothing to do.
    if (hp[(k0[0] >> shift) & (RDX_SIZE - 1)] == M) continue;

    size_t sum = 0;
    for (int b = 0; b < RDX_SIZE; b++) {
      size_t c = hp[b];
      hp[b]    = sum;
      sum     += c;
    }

    for (size_t t = 0; t < M; t++) {
      size_t dst = hp[(k0[t] >> shift) & (RDX_SIZE - 1)]++;
      k1[dst] = k0[t];
      i1[dst] = i0[t];
    }

    uint64_t *tk = k0; k0 = k1; k1 = tk;
    uint32_t *ti = i0; i0 = i1; i1 = ti;
  }

  for (size_t a = 0; a < M; ) {
    size_t b = a;
    while (b + 1 < M && k0[b + 1] == k0[a]) b++;
    double r2 = (double) (a + 1) + (double) (b + 1);
    for (size_t t = a; t <= b; t++) rk[i0[t]] = r2;
    a = b + 1;
  }

  vmaxset(vmax);
}


//------------------------------------------------------
// Each task is one block of PERM_FW permutations. L holds
// each sample's group (as a double, so comparisons
// vectorize alongside the sums) under each permutation:
// chunk h is L + h * n * PERM_W, row-major (n x PERM_W).
// Unused lanes are group 0, and ignored.
//
// For one chunk, column j of the triangle, and nd of the
// rows below it:
//
//     s[k] += sum_t rk[t] * (li[t, k] == lj[k])
//
// The rows go in tiles, so that their labels (TILE_ROWS x
// PERM_FW doubles: 128 KB) stay in cache while each column
// passes over them, instead of streaming all n rows' labels
// for every column.
//
// GCC and clang get explicit two-lane vectors: compilers
// don't vectorize these loops at R's default -O2. Each lane
// does the same operations in the same order as the scalar
// fallback, so results don't depend on which is used (and,
// the sums being of whole numbers, are exact).
//------------------------------------------------------
#define PERM_W  8
#define PERM_FW (4 * PERM_W)
#define TILE_ROWS 512
#define PERM_REP(F) F(0) F(1) F(2) F(3) F(4) F(5) F(6) F(7)

#if defined(__GNUC__) || defined(__clang__)

typedef double    v2d __attribute__((vector_size(16)));
typedef long long v2l __attribute__((vector_size(16)));

static inline v2d load2(const double *p) { v2d x; memcpy(&x, p, sizeof(x)); return x; }

#define V2_REP(F) F(0) F(1) F(2) F(3)

static inline void within_chunk(
    const double *rk, const double *li, const double *lj, int nd, double *s) {

  #define W_INIT(h) v2d acc##h = { 0, 0 }, g##h = load2(lj + 2 * h);
  #define W_STEP(h) acc##h += (v2d) ((v2l) a & (v2l) (load2(li + 2 * h) == g##h));
  #define W_DONE(h) s[2 * h] += acc##h[0]; s[2 * h + 1] += acc##h[1];

  V2_REP(W_INIT)
  for (int t = 0; t < nd; t++, li += PERM_W) {
    v2d a = { rk[t], rk[t] };
    V2_REP(W_STEP)
  }
  V2_REP(W_DONE)
}

// Two columns at once, over the same nd rows: each row's
// labels are loaded once for both, and the eight sums are
// independent. Every sum runs in the same order as with
// one column at a time.
static inline void within_chunk2(
    const double *rk0, const double *rk1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s) {

  #define W2_INIT(h) v2d acc##h = { 0, 0 }, bcc##h = { 0, 0 },                      g##h = load2(lj0 + 2 * h), k##h = load2(lj1 + 2 * h);
  #define W2_STEP(h) { v2d x = load2(li + 2 * h);                        acc##h += (v2d) ((v2l) a & (v2l) (x == g##h));                        bcc##h += (v2d) ((v2l) b & (v2l) (x == k##h)); }
  #define W2_DONE(h)     s[2 * h] += acc##h[0]; s[2 * h + 1] += acc##h[1];     s[2 * h] += bcc##h[0]; s[2 * h + 1] += bcc##h[1];

  V2_REP(W2_INIT)
  for (int t = 0; t < nd; t++, li += PERM_W) {
    v2d a = { rk0[t], rk0[t] };
    v2d b = { rk1[t], rk1[t] };
    V2_REP(W2_STEP)
  }
  V2_REP(W2_DONE)
}

#else

static inline void within_chunk(
    const double *rk, const double *li, const double *lj, int nd, double *s) {

  #define W_INIT(k) double acc##k = 0, g##k = lj[k];
  #define W_STEP(k) acc##k += (li[k] == g##k) ? a : 0.0;
  #define W_DONE(k) s[k] += acc##k;

  PERM_REP(W_INIT)
  for (int t = 0; t < nd; t++, li += PERM_W) {
    double a = rk[t];
    PERM_REP(W_STEP)
  }
  PERM_REP(W_DONE)
}

static inline void within_chunk2(
    const double *rk0, const double *rk1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s) {
  within_chunk(rk0, li, lj0, nd, s);
  within_chunk(rk1, li, lj1, nd, s);
}

#endif


static void *within_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = an_n;
  int nc        = PERM_FW / PERM_W;

  double *L = an_work + (size_t) thread_i * an_wlen;   // nc x n x W
  double *s = L + (size_t) n * PERM_FW;                // PERM_FW
  int    *p = an_iwork + (size_t) thread_i * an_iwlen; // n

  for (int blk = an_blk0 + thread_i; blk < an_blk1; blk += n_threads) {

    int k0 = blk * PERM_FW;
    int nb = an_total - k0;
    if (nb > PERM_FW) nb = PERM_FW;

    memset(L, 0, (size_t) n * PERM_FW * sizeof(double));

    // Sample a gets the group of sample p[a].
    for (int b = 0; b < nb; b++) {
      double *lc = L + (size_t) (b / PERM_W) * n * PERM_W + b % PERM_W;
      perm_source_get(&an_ps, k0 + b, p);
      for (int a = 0; a < n; a++)
        lc[(size_t) a * PERM_W] = an_grp[p[a]];
    }

    for (int b = 0; b < PERM_FW; b++) s[b] = 0;

    // Rows in tiles of TILE_ROWS, whose labels stay in
    // cache while every column above them goes by.
    for (int i0 = 1; i0 < n; i0 += TILE_ROWS) {

      int i1 = (n - i0 > TILE_ROWS) ? i0 + TILE_ROWS : n;
      int j  = 0;

      // Pairs of columns that span the whole tile.
      for (; j + 1 < i0; j += 2) {
        const double *r0 = an_rk + dist_col_start(j,     n) + (i0 - j - 1);
        const double *r1 = an_rk + dist_col_start(j + 1, n) + (i0 - j - 2);
        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * PERM_W;
          within_chunk2(
            r0, r1, vc + (size_t) i0 * PERM_W,
            vc + (size_t) j * PERM_W, vc + (size_t) (j + 1) * PERM_W,
            i1 - i0, s + h * PERM_W );
        }
      }

      // The rest: rows max(i0, j + 1) to i1 - 1 of column j.
      for (; j < i1 - 1; j++) {
        int start = (j + 1 > i0) ? j + 1 : i0;
        const double *r0 = an_rk + dist_col_start(j, n) + (start - j - 1);
        for (int h = 0; h < nc; h++) {
          const double *vc = L + (size_t) h * n * PERM_W;
          within_chunk(
            r0, vc + (size_t) start * PERM_W, vc + (size_t) j * PERM_W,
            i1 - start, s + h * PERM_W );
        }
      }
    }

    for (int b = 0; b < nb; b++)
      an_out[k0 + b] = s[b];
  }

  return NULL;
}


//======================================================
// R interface.
//
// dm:     packed lower triangle of D (double), with no NA.
// n:      number of samples.
// grp:    each sample's group, 0 .. g - 1 (integer).
// perms:  n x n_perm integer matrix of 0-based permutations,
//         or NULL for random ones.
// strata: each sample's stratum, 0 .. n_strata - 1, or NULL.
// opts:   c(n_perm, seed, n_threads).
//
// Returns list(sw, group_sw, group_pairs): sw is the doubled
// sum of within-group ranks for the identity, then each
// permutation; group_sw and group_pairs are the observed
// doubled rank sum and number of pairs within each group.
//======================================================
SEXP C_anosim(
    SEXP sexp_dm, SEXP sexp_n, SEXP sexp_grp, SEXP sexp_perms,
    SEXP sexp_strata, SEXP sexp_opts ) {

  int n         = asInteger(sexp_n);
  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];
  size_t M      = (size_t) n * (size_t) (n - 1) / 2;

  if (n < 2 || n_perm < 0 || n_threads < 1 ||
      XLENGTH(sexp_dm) != (R_xlen_t) M || XLENGTH(sexp_grp) != n)
    error("Invalid arguments to C_anosim()."); // # nocov

  if (M > (size_t) UINT32_MAX)
    error("`dm` has too many samples for anosim()."); // # nocov

  an_n     = n;
  an_grp   = INTEGER(sexp_grp);
  an_total = n_perm + 1;


  //--------------------------------------------------
  // Doubled ranks of the distances.
  //--------------------------------------------------
  double *rk = (double *) R_alloc(M, sizeof(double));
  rank_values(REAL(sexp_dm), M, rk);
  an_rk = rk;

  R_CheckUserInterrupt();


  //--------------------------------------------------
  // Observed rank sums and pairs within each group.
  //--------------------------------------------------
  int g = 0;
  for (int a = 0; a < n; a++)
    if (an_grp[a] < 0 || an_grp[a] >= n) error("Invalid arguments to C_anosim()."); // # nocov
    else if (an_grp[a] + 1 > g) g = an_grp[a] + 1;

  SEXP sexp_gsw = PROTECT(allocVector(REALSXP, g));
  SEXP sexp_gnp = PROTECT(allocVector(REALSXP, g));
  double *gsw = REAL(sexp_gsw), *gnp = REAL(sexp_gnp);
  for (int c = 0; c < g; c++) gsw[c] = gnp[c] = 0;

  for (int j = 0; j < n - 1; j++) {
    const double *r  = rk + dist_col_start(j, n);
    int           gj = an_grp[j];
    for (int i = j + 1; i < n; i++) {
      if (an_grp[i] == gj) {
        gsw[gj] += r[i - j - 1];
        gnp[gj] += 1;
      }
    }
  }


  //--------------------------------------------------
  // Permutations: supplied, or shuffles within strata.
  //--------------------------------------------------
  perm_source_init(&an_ps, n, sexp_perms, sexp_strata, seed);


  //--------------------------------------------------
  // Blocks, and workspace for each thread.
  //--------------------------------------------------
  int n_blocks = (int) (((size_t) an_total + PERM_FW - 1) / PERM_FW);
  if (n_threads > n_blocks) n_threads = n_blocks;

  an_wlen  = (size_t) n * PERM_FW + PERM_FW;
  an_iwlen = (size_t) n;
  an_work  = (double *) R_alloc(an_wlen  * n_threads, sizeof(double));
  an_iwork = (int *)    R_alloc(an_iwlen * n_threads, sizeof(int));

  SEXP sexp_sw = PROTECT(allocVector(REALSXP, an_total));
  an_out = REAL(sexp_sw);


  //--------------------------------------------------
  // Rounds of blocks, checking for interrupts between
  // them. Each round has at least ~1e8 pair visits per
  // thread, or every remaining block. Each share of a
  // round (one per thread) takes every nt-th block, and
  // has its own workspace.
  //--------------------------------------------------
  double blk_work  = (double) M * PERM_FW + 1;
  double rounds    = 1e8 * n_threads / blk_work;
  int    per_round = (rounds < n_blocks) ? (int) rounds : n_blocks; // no overflow
  if (per_round < n_threads) per_round = n_threads;

  for (an_blk0 = 0; an_blk0 < n_blocks; an_blk0 = an_blk1) {

    an_blk1 = (n_blocks - an_blk0 > per_round) ? an_blk0 + per_round : n_blocks;

    int nb = an_blk1 - an_blk0;
    int nt = threads_for(nb < n_threads ? nb : n_threads, nb * blk_work, 1e6);

    // run_parallel() starts threads only when told of at
    // least 100 tasks.
    run_parallel(within_worker, nt, nt > 1 ? 100 : 1);

    if (an_blk1 < n_blocks) R_CheckUserInterrupt();
  }

  //--------------------------------------------------
  // R = (r_B - r_W) / (M / 2). The mean ranks are
  // divided in long double, then rounded, as R's mean()
  // does (vegan's statistic): bit-identical with it where
  // the rounding agrees, within an ulp elsewhere.
  //--------------------------------------------------
  double n_w = 0;
  for (int c = 0; c < g; c++) n_w += gnp[c];
  double n_b = (double) M - n_w;
  double tot = (double) M * ((double) M + 1);  // Doubled.

  SEXP sexp_r = PROTECT(allocVector(REALSXP, an_total));
  double *r = REAL(sexp_r);
  for (int k = 0; k < an_total; k++) {
    double sw  = an_out[k];
    double r_w = (double) ((long double) sw / 2 / (long double) n_w);
    double r_b = (double) ((long double) (tot - sw) / 2 / (long double) n_b);
    r[k] = (r_b - r_w) / ((double) M / 2);
  }

  const char *names[] = { "r", "sw", "group_sw", "group_pairs", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_r);
  SET_VECTOR_ELT(sexp_result, 1, sexp_sw);
  SET_VECTOR_ELT(sexp_result, 2, sexp_gsw);
  SET_VECTOR_ELT(sexp_result, 3, sexp_gnp);

  UNPROTECT(5);
  return sexp_result;
}
