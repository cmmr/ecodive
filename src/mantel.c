// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"

#include <limits.h> // INT_MAX


//======================================================
// Mantel tests: the correlation of one distance matrix
// (x) with one or two others (targets), for the observed
// samples and for each permutation of x's samples.
//
// Permutation p relabels x: target pair (i, j) meets
// x[p[i], p[j]], as with vegan's xmat[take, take]. The
// targets stay fixed. Equivalently, x pair (a, b) meets
// target pair (q[a], q[b]), for q the inverse of p.
//
// Pearson (and Spearman, given ranks): x and the targets
// are centred and scaled to unit length once, so each
// correlation is a dot product over the pairs, with x
// read through the permutation. O(N) per permutation and
// target, for N = n(n - 1) / 2 pairs.
//
// Kendall: tau-b, as stats::cor(), after Knight (1966).
// The pairs are visited in x's order, once and for all;
// each permutation reads the target's integer codes in
// that order, through q. The discordant pairs are then
// the target's strict inversions between x's tie groups:
// counted with a Fenwick tree over the codes when there
// are few distinct target values, or else by sorting
// within x's tie groups and then counting inversions
// with a merge sort. O(N log N) per permutation and
// target.
//
// Each permutation is computed by one thread, start to
// finish, so results are identical for any thread count.
//======================================================


#define MT_PEARSON 0
#define MT_KENDALL 1

// Fenwick tree for targets with at most this many codes.
#define MT_SMALL_K 4096

// Merge sort: insertion-sorted runs of this length.
#define MT_RUN 32

static int           mt_n;       // Samples.
static size_t        mt_np;      // Pairs, n(n - 1) / 2.
static int           mt_m;       // Targets, 1 or 2.
static int           mt_method;
static int           mt_total;   // Permutations + 1.
static ptrdiff_t    *mt_cs1;     // Offset of pair (a, b), a > b: mt_cs1[b] + a.
static double       *mt_out;     // m x mt_total correlations.

// Pearson.
static const double *mt_x;       // Standardized x.
static const double *mt_t[2];    // Standardized targets.

// Kendall.
static int          *mt_xa;      // x's pairs, in x's order: samples a > b.
static int          *mt_xb;
static const int    *mt_gbeg;    // Start of each x tie group in that order.
static int           mt_ng;      // x tie groups.
static const int    *mt_tc[2];   // Targets' codes, 0 .. K - 1, in pair order.
static int           mt_k[2];    // Targets' distinct values, K.
static double        mt_n0;      // Pairs of pairs, N(N - 1) / 2.
static double        mt_n1;      // ... tied in x.
static double        mt_n2[2];   // ... tied in each target.

// Permutations: permuted sample a is original sample p[a],
// as vegan's xmat[take, take].
static perm_source_t mt_ps;

// Work.
static int          *mt_iwork;   // Per-thread integer workspace ...
static size_t        mt_iwlen;   // ... of this many ints.
static int           mt_k0;      // This round's permutations:
static int           mt_k1;      //   k0 .. k1 - 1.


//------------------------------------------------------
// Pearson: dot products of permuted x with the targets.
// p: the permutation; q: its inverse.
//
// Target pair (i, j) meets x[p[i], p[j]]. The columns j
// are visited in order of b = p[j], so that the x[a, b]
// with a < b, one in each column a of the packed
// triangle, are next to the previous column's.
//------------------------------------------------------
static void pearson_perm(const int *p, const int *q, double *out) {

  int              n   = mt_n;
  const double    *x   = mt_x;
  const double    *t0  = mt_t[0];
  const double    *t1  = mt_t[1];
  const ptrdiff_t *cs1 = mt_cs1;

  double s0 = 0, s1 = 0;

  for (int b = 0; b < n; b++) {

    int       j  = q[b];
    ptrdiff_t xb = cs1[b];
    ptrdiff_t tb = cs1[j];  // tb + i: target pair (i, j)

    if (mt_m == 1) {
      for (int i = j + 1; i < n; i++) {
        int       a  = p[i];
        ptrdiff_t lo = cs1[a] + b, hi = xb + a;
        s0 += x[(a > b) ? hi : lo] * t0[tb + i];
      }
    }
    else {
      for (int i = j + 1; i < n; i++) {
        int       a  = p[i];
        ptrdiff_t lo = cs1[a] + b, hi = xb + a;
        double    v  = x[(a > b) ? hi : lo];
        s0 += v * t0[tb + i];
        s1 += v * t1[tb + i];
      }
    }
  }

  out[0] = s0;
  if (mt_m == 2) out[1] = s1;
}


//------------------------------------------------------
// Sort a[0 .. len) with buf (len ints) as scratch.
// Returns the inversions: pairs with the earlier value
// strictly greater than the later. Also bioenv()'s.
//------------------------------------------------------
int64_t kendall_sort_count(int *a, int *buf, size_t len) {

  int64_t inv = 0;

  // Insertion sort runs of MT_RUN; each shift is an
  // inversion.
  for (size_t lo = 0; lo < len; lo += MT_RUN) {
    size_t hi = (lo + MT_RUN < len) ? lo + MT_RUN : len;
    for (size_t i = lo + 1; i < hi; i++) {
      int    v = a[i];
      size_t k = i;
      while (k > lo && a[k - 1] > v) { a[k] = a[k - 1]; k--; }
      inv += (int64_t) (i - k);
      a[k] = v;
    }
  }

  // Merge runs, bottom up. Taking from the right while
  // left values remain makes that many inversions.
  int *src = a, *dst = buf;
  for (size_t w = MT_RUN; w < len; w *= 2) {
    for (size_t lo = 0; lo < len; lo += 2 * w) {
      size_t mid = (lo + w < len)     ? lo + w     : len;
      size_t hi  = (lo + 2 * w < len) ? lo + 2 * w : len;
      size_t l = lo, r = mid, o = lo;
      while (l < mid && r < hi) {
        int al = src[l], ar = src[r];
        int right = ar < al;
        dst[o++] = right ? ar : al;
        inv += right ? (int64_t) (mid - l) : 0;
        r += (size_t) right;
        l += (size_t) !right;
      }
      while (l < mid) dst[o++] = src[l++];
      while (r < hi)  dst[o++] = src[r++];
    }
    int *tmp = src; src = dst; dst = tmp;
  }
  if (src != a) memcpy(a, src, len * sizeof(int));

  return inv;
}


//------------------------------------------------------
// Pairs with equal values in sorted a[0 .. len).
//------------------------------------------------------
double kendall_tied_sorted(const int *a, size_t len) {
  double ties = 0;
  size_t run  = 1;
  for (size_t t = 1; t < len; t++) {
    if (a[t] == a[t - 1]) { run++; continue; }
    ties += (double) run * (double) (run - 1) / 2;
    run   = 1;
  }
  return ties + (double) run * (double) (run - 1) / 2;
}


//------------------------------------------------------
// Kendall's tau-b of permuted x with target c.
// q: the inverse permutation. work: 2N + 2 * MT_SMALL_K
// + 2 ints.
//------------------------------------------------------
static double kendall_perm(const int *q, int c, int *work) {

  size_t           N    = mt_np;
  int              K    = mt_k[c];
  int              ng   = mt_ng;
  const int       *gbeg = mt_gbeg;
  const int       *xa   = mt_xa;
  const int       *xb   = mt_xb;
  const int       *tc   = mt_tc[c];
  const ptrdiff_t *cs1  = mt_cs1;

  int *seq = work;
  int *buf = work + N;

  // The target's codes, in x's order: x pair (a, b) meets
  // target pair (q[a], q[b]).
  for (size_t s = 0; s < N; s++) {
    int       i  = q[xa[s]], j = q[xb[s]];
    ptrdiff_t lo = cs1[i] + j, hi = cs1[j] + i;
    seq[s] = tc[(i > j) ? hi : lo];
  }

  double n3 = 0, disc = 0;

  if (K <= MT_SMALL_K) {

    // Fenwick tree of the codes seen in earlier x tie
    // groups; cnt counts codes within a group.
    int *tree = buf;
    int *cnt  = buf + K + 1;
    memset(tree, 0, (size_t) (K + 1) * sizeof(int));
    memset(cnt,  0, (size_t) K * sizeof(int));

    int64_t seen = 0, less_eq = 0, d = 0;

    for (int g = 0; g < ng; g++) {

      int lo = gbeg[g], hi = gbeg[g + 1];

      // Discordant: earlier codes strictly greater.
      for (int s = lo; s < hi; s++) {
        less_eq = 0;
        for (int v = seq[s] + 1; v > 0; v -= v & -v) less_eq += tree[v];
        d += seen - less_eq;
      }

      if (hi - lo > 1) {
        int64_t t3 = 0;
        for (int s = lo; s < hi; s++) t3 += cnt[seq[s]]++;
        for (int s = lo; s < hi; s++) cnt[seq[s]] = 0;
        n3 += (double) t3;
      }

      for (int s = lo; s < hi; s++)
        for (int v = seq[s] + 1; v <= K; v += v & -v) tree[v]++;
      seen += hi - lo;
    }

    disc = (double) d;
  }
  else {

    // Sort within x's tie groups, so that only pairs
    // between groups can be inversions.
    for (int g = 0; g < ng; g++) {
      size_t lo = (size_t) gbeg[g], len = (size_t) gbeg[g + 1] - lo;
      if (len < 2) continue;
      kendall_sort_count(seq + lo, buf, len);
      n3 += kendall_tied_sorted(seq + lo, len);
    }

    disc = (double) kendall_sort_count(seq, buf, N);
  }

  // As stats::cor(): the sum over ordered pairs of pairs
  // of sign products, over the square roots of the untied
  // ordered pairs in each.
  double n0 = mt_n0, n1 = mt_n1, n2 = mt_n2[c];
  double s  = 2 * (n0 - n1 - n2 + n3 - 2 * disc);
  double r  = s / (sqrt(2 * (n0 - n1)) * sqrt(2 * (n0 - n2)));

  return (r > 1) ? 1 : r;
}


//------------------------------------------------------
// Worker: permutations k0 + i, k0 + i + nt, ...
//------------------------------------------------------
static void *mantel_worker(void *arg) {

  worker_t *w    = (worker_t *) arg;
  int      *p    = mt_iwork + (size_t) w->i * mt_iwlen;
  int      *q    = p + mt_n;
  int      *work = q + mt_n;

  for (int k = mt_k0 + w->i; k < mt_k1; k += w->n) {

    double *out = mt_out + (size_t) k * mt_m;
    perm_source_get(&mt_ps, k, p);
    for (int a = 0; a < mt_n; a++) q[p[a]] = a;

    if (mt_method == MT_PEARSON) {
      pearson_perm(p, q, out);
    }
    else {
      for (int c = 0; c < mt_m; c++)
        out[c] = kendall_perm(q, c, work);
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// x:       Pearson: standardized packed x (double).
//          Kendall: x's pairs in x's order (0-based pair
//          indices, sorted by value; integer).
// targets: List of 1 or 2 targets, in pair order.
//          Pearson: standardized (double). Kendall: codes
//          0 .. K - 1 of the distinct values (integer).
// groups:  Kendall: x's tie group starts in x's order
//          (integer, groups + 1, ending in N); else NULL.
// ties:    Kendall: c(n1, n2[, n2']) tied pairs in x and
//          the targets; else NULL.
// perms:   n x n_perm integer matrix of 0-based
//          permutations, or NULL for random ones.
// strata:  Each sample's stratum, 0 .. n_strata - 1, or NULL.
// opts:    c(n, method, n_perm, seed, n_threads).
//
// Returns an m x (n_perm + 1) matrix of correlations,
// column 1 for the observed samples.
//======================================================
SEXP C_mantel(
    SEXP sexp_x, SEXP sexp_targets, SEXP sexp_groups, SEXP sexp_ties,
    SEXP sexp_perms, SEXP sexp_strata, SEXP sexp_opts ) {

  int n         = INTEGER(sexp_opts)[0];
  int method    = INTEGER(sexp_opts)[1];
  int n_perm    = INTEGER(sexp_opts)[2];
  int seed      = INTEGER(sexp_opts)[3];
  int n_threads = INTEGER(sexp_opts)[4];
  int m         = LENGTH(sexp_targets);

  size_t N = (size_t) n * (size_t) (n - 1) / 2;

  if (n < 2 || n_perm < 0 || n_threads < 1 || m < 1 || m > 2 ||
      (size_t) XLENGTH(sexp_x) != N || N > INT_MAX)
    error("Invalid arguments to C_mantel()."); // # nocov

  mt_n      = n;
  mt_np     = N;
  mt_m      = m;
  mt_method = method;
  mt_total  = n_perm + 1;

  mt_cs1 = (ptrdiff_t *) R_alloc(n, sizeof(ptrdiff_t));
  for (int b = 0; b < n; b++)
    mt_cs1[b] = (ptrdiff_t) dist_col_start(b, n) - b - 1;


  //--------------------------------------------------
  // The statistic's inputs.
  //--------------------------------------------------
  size_t wlen = 0;

  if (method == MT_PEARSON) {
    mt_x = REAL(sexp_x);
    for (int c = 0; c < m; c++) mt_t[c] = REAL(VECTOR_ELT(sexp_targets, c));
    if (m == 1) mt_t[1] = NULL;
  }
  else {

    // Each pair's samples, in x's order.
    const int *xord = INTEGER(sexp_x);
    int       *pos  = (int *) R_alloc(N, sizeof(int));
    for (size_t s = 0; s < N; s++) pos[xord[s]] = (int) s;

    mt_xa = (int *) R_alloc(N, sizeof(int));
    mt_xb = (int *) R_alloc(N, sizeof(int));
    size_t t = 0;
    for (int b = 0; b < n - 1; b++)
      for (int a = b + 1; a < n; a++, t++) {
        mt_xa[pos[t]] = a;
        mt_xb[pos[t]] = b;
      }

    mt_gbeg = INTEGER(sexp_groups);
    mt_ng   = LENGTH(sexp_groups) - 1;
    mt_n0   = (double) N * (double) (N - 1) / 2;
    mt_n1   = REAL(sexp_ties)[0];

    for (int c = 0; c < m; c++) {
      const int *tc = INTEGER(VECTOR_ELT(sexp_targets, c));
      int K = 0;
      for (size_t u = 0; u < N; u++) if (tc[u] + 1 > K) K = tc[u] + 1;
      mt_tc[c] = tc;
      mt_k[c]  = K;
      mt_n2[c] = REAL(sexp_ties)[1 + c];
    }

    wlen = 2 * N + 2 * (size_t) MT_SMALL_K + 2; // seq, buf (or tree + cnt)
  }


  //--------------------------------------------------
  // Permutations: supplied, or shuffles within strata.
  //--------------------------------------------------
  perm_source_init(&mt_ps, n, sexp_perms, sexp_strata, seed);


  //--------------------------------------------------
  // Workspace and output.
  //--------------------------------------------------
  if (n_threads > mt_total) n_threads = mt_total;

  mt_iwlen = 2 * (size_t) n + wlen;
  mt_iwork = (int *) R_alloc(mt_iwlen * n_threads, sizeof(int));

  SEXP sexp_r = PROTECT(allocMatrix(REALSXP, m, mt_total));
  mt_out = REAL(sexp_r);


  //--------------------------------------------------
  // Rounds of permutations, checking for interrupts
  // between them. Each round has at least ~1e8
  // operations per thread, or every remaining one.
  //--------------------------------------------------
  double per_perm = (double) N * m + 1;
  if (method == MT_KENDALL) per_perm *= 2 * log2((double) N + 2);

  int per_round = (int) fmin(1e8 * n_threads / per_perm, (double) mt_total);
  if (per_round < n_threads) per_round = n_threads;

  for (mt_k0 = 0; mt_k0 < mt_total; mt_k0 = mt_k1) {

    mt_k1 = (mt_total - mt_k0 > per_round) ? mt_k0 + per_round : mt_total;

    // Threads only for enough work: at least ~1e6
    // operations each. run_parallel() starts them only
    // when told of at least 100 tasks.
    int nt = threads_for(n_threads, per_perm * (mt_k1 - mt_k0), 1e6);
    if (nt > mt_k1 - mt_k0) nt = mt_k1 - mt_k0;
    run_parallel(mantel_worker, nt, (nt > 1) ? 100 : 1);

    if (mt_k1 < mt_total) R_CheckUserInterrupt();
  }

  UNPROTECT(1);
  return sexp_r;
}
