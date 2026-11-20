// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// SparCC (Friedman & Alm 2012), written from the paper alone (CC-BY). See
// plans/sparcc.md: no other implementation's source was read.
//
// For each Monte Carlo draw d = 1 .. n_draws:
//
//   1. Draw each sample's log fractions from the Dirichlet posterior
//      Dir(N_s + 1), as the logs of Gamma(N_sf + 1) draws: only log ratios
//      are used, so their total doesn't matter. Sample s of draw d of
//      dataset b (0 = the observed table; 1 .. B the null datasets of the
//      permutation test) has its own PCG stream, (b * n_draws + d) * (n + 1)
//      + s, and draws its features in table order, so results don't depend
//      on the number of threads, and sparse and dense tables agree.
//   2. The variation matrix t_ij = var_s(g_si - g_sj), divisor n - 1, from
//      the differences of the centred logs (as vlr() in proportionality.c;
//      never c_ii + c_jj - 2 c_ij).
//   3. The basis variances w_i, from t_i = sum_{j ~ i} t_ij = deg_i w_i +
//      sum_{j ~ i} w_j over the pairs not excluded: M w = t with M = B'B,
//      B the incidence matrix of the pairs that are kept. With no pair
//      excluded, M = (m - 2) I + 11', whose inverse is closed form; each
//      excluded pair (i, j) subtracts u u', u = e_i + e_j, so M = A - U U'
//      is solved with Woodbury's identity, a K x K Cholesky factorisation
//      (K <= max_exclusions), never an m x m one.
//   4. rho_ij = (w_i + w_j - t_ij) / (2 sqrt(w_i) sqrt(w_j)); NaN when either
//      variance isn't positive.
//   5. The exclusion loop, at most max_exclusions times: the pair with the
//      largest |rho| not yet excluded (the first in packed order on ties)
//      is excluded if |rho| > exclusion_threshold, else the loop stops. A
//      component left with no pair is dropped from the draw (NA). The loop
//      stops when fewer than 4 components would remain (keeping the last
//      estimates), or when the reduced system is singular (a bipartite set
//      of kept pairs: that exclusion is undone). Then step 3 again, and
//      every rho from the new variances.
//
//      Step 3 is the least-squares fit of t_ij ~ w_i + w_j over the kept
//      pairs (M = B'B, t_i = B't), and rho_ij is proportional to that
//      fit's residual. A pair whose exclusion would drop a component or
//      make M singular is one that no other pair's equation spans
//      (leverage 1), so its residual, and rho, is exactly 0: neither case
//      can happen for a threshold above rounding error. They are handled
//      for exclusion_threshold = 0, where rounding decides.
//
// The result for each pair is the median over draws of its non-missing
// values; each feature's variance likewise.
//
// Threads: one task per draw. Each draw's values are written to their own
// slots of a pair-major store, then each pair's median is taken alone, so
// results are identical for any number of threads. When the store of
// n_draws x pairs values would exceed SP_STORE_MAX bytes, pairs are stored
// in blocks, and the draws are made again for each later block (each
// draw's variances are kept from the first pass): the streams make the
// same numbers every time.
//
// Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from
// genomic survey data. PLoS Computational Biology, 8, e1002687.

#include "ecodive.h"

#define SP_PERM_STREAM ((uint64_t) 1 << 62)  // + (b - 1) * D + f: null datasets
#define SP_STORE_MAX   ((double) 256 * 1024 * 1024)
#define SP_PIVOT_MIN   1e-10

// Why a draw's exclusion loop stopped.
#define SP_STOP_MAX       1  // max_exclusions made
#define SP_STOP_THRESHOLD 2  // no pair left above the threshold
#define SP_STOP_THREE     3  // fewer than 4 components left
#define SP_STOP_SINGULAR  4  // the next exclusion made M singular (undone)

static int            sp_n, sp_D, sp_n_draws, sp_max_excl;
static double         sp_thresh;
static uint64_t       sp_seed, sp_dataset;
static size_t         sp_pairs;
static const double  *sp_cnt;     // D x n counts of the dataset, feature-major
static const double  *sp_given;   // n x D x n_draws log abundances, or NULL
static double        *sp_w2;      // n_draws x D: each draw's final variances
static int           *sp_n_excl;  // n_draws: exclusions in each draw
static int           *sp_stop;    // n_draws: why its loop stopped (SP_STOP_*), or NULL
static int           *sp_excl;    // n_draws x max_excl: their packed indices
static int            sp_solve_pass;  // 1: the whole draw; 0: t and rho only
static size_t         sp_p0, sp_p1;   // the block of pairs stored
static double        *sp_store;       // (p1 - p0) x n_draws, pair by pair
static double        *sp_med;         // pairs: medians (this block's part)
static int            sp_n_thr;

// Scratch, one set per thread (index i of run_parallel()'s shares).
static double        *sp_g, *sp_t, *sp_rho, *sp_w, *sp_wnew, *sp_ti, *sp_S, *sp_z;
static unsigned char *sp_ex, *sp_act;
static int           *sp_deg, *sp_er, *sp_ec;


// Pair (r, c), r < c, of D features: its position in the packed triangle.
static inline size_t pair_index (int r, int c, int D) {
  return (size_t) r * D - (size_t) r * (r + 1) / 2 + (size_t) (c - r - 1);
}

// The row r of the packed triangle that holds pair index p.
static int pair_row (size_t p, int D) {
  int r = 0;
  while (r < D - 1 && pair_index(r, D - 1, D) < p) r++;
  return r;
}

static inline double sp_rho_of (double wi, double wj, double t) {
  if (!(wi > 0 && wj > 0)) return R_NaN;
  return (wi + wj - t) / (2 * sqrt(wi) * sqrt(wj));
}


// Draw d's centred log abundances, feature by feature, into g (D x n).
static void sp_fill (int d, double *g) {

  int n = sp_n, D = sp_D;

  if (sp_given) {
    const double *x = sp_given + (size_t) d * n * D;
    for (size_t k = 0; k < (size_t) n * D; k++) g[k] = x[k];
  }
  else {
    rng_t rng;
    for (int s = 0; s < n; s++) {
      uint64_t stream = (sp_dataset * (uint64_t) sp_n_draws + (uint64_t) d) *
                        ((uint64_t) n + 1) + (uint64_t) s;
      rng_seed(&rng, sp_seed, stream);
      for (int f = 0; f < D; f++)
        g[(size_t) f * n + s] = rng_log_gamma(&rng, sp_cnt[(size_t) f * n + s] + 1);
    }
  }

  for (int f = 0; f < D; f++) {
    double *x = g + (size_t) f * n, sum = 0;
    for (int s = 0; s < n; s++) sum += x[s];
    double mean = sum / n;
    for (int s = 0; s < n; s++) x[s] -= mean;
  }
}


// t for the pairs p0 .. p1 - 1, into t[0 .. p1 - p0 - 1].
static void sp_variation (const double *g, double *t, size_t p0, size_t p1) {

  int    n  = sp_n, D = sp_D;
  double n1 = n - 1;
  size_t p  = p0;
  int    r  = pair_row(p0, D);
  int    c  = r + 1 + (int) (p0 - pair_index(r, r + 1, D));

  for (; p < p1; p++) {
    const double *a = g + (size_t) r * n, *b = g + (size_t) c * n;
    double ss = 0;
    for (int s = 0; s < n; s++) { double x = a[s] - b[s]; ss += x * x; }
    t[p - p0] = ss / n1;
    if (++c == D) { r++; c = r + 1; }
  }
}


// The basis variances of the m active components (act), given the pairs
// excluded among them (er[k], ec[k], k < K): solves M w = t_i. Returns 0
// (w unchanged) when M is singular.
static int sp_solve (const double *t, const unsigned char *act, int m,
                     const unsigned char *ex, const int *er, const int *ec,
                     int K, double *ti, double *S, double *z, double *w) {

  int D = sp_D;

  for (int i = 0; i < D; i++) ti[i] = 0;
  size_t p = 0;
  for (int r = 0; r < D - 1; r++)
    for (int c = r + 1; c < D; c++, p++)
      if (act[r] && act[c] && !ex[p]) { ti[r] += t[p]; ti[c] += t[p]; }

  double a = m - 2, h = 1.0 / (2.0 * m - 2), sum = 0;
  for (int i = 0; i < D; i++) if (act[i]) sum += ti[i];

  // y = A^-1 t, kept in ti.
  for (int i = 0; i < D; i++) if (act[i]) ti[i] = (ti[i] - h * sum) / a;

  if (K > 0) {

    // S = I - U' A^-1 U, with A^-1 u = (u - 2h 1) / a.
    for (int k = 0; k < K; k++)
      for (int l = 0; l <= k; l++) {
        int share = (er[k] == er[l]) + (er[k] == ec[l]) + (ec[k] == er[l]) + (ec[k] == ec[l]);
        S[k * K + l] = (k == l) - (share - 4 * h) / a;
      }

    // Cholesky, lower triangle.
    for (int j = 0; j < K; j++) {
      double v = S[j * K + j];
      for (int k = 0; k < j; k++) v -= S[j * K + k] * S[j * K + k];
      if (!(v > SP_PIVOT_MIN)) return 0;
      v = sqrt(v);
      S[j * K + j] = v;
      for (int i = j + 1; i < K; i++) {
        double x = S[i * K + j];
        for (int k = 0; k < j; k++) x -= S[i * K + k] * S[j * K + k];
        S[i * K + j] = x / v;
      }
    }

    // z = S^-1 U' y.
    for (int k = 0; k < K; k++) z[k] = ti[er[k]] + ti[ec[k]];
    for (int i = 0; i < K; i++) {
      double x = z[i];
      for (int k = 0; k < i; k++) x -= S[i * K + k] * z[k];
      z[i] = x / S[i * K + i];
    }
    for (int i = K - 1; i >= 0; i--) {
      double x = z[i];
      for (int k = i + 1; k < K; k++) x -= S[k * K + i] * z[k];
      z[i] = x / S[i * K + i];
    }

    // w = y + A^-1 U z = y + (U z - 2h (sum z) 1) / a.
    double zs = 0;
    for (int k = 0; k < K; k++) zs += z[k];
    for (int i = 0; i < D; i++) if (act[i]) w[i] = ti[i] - 2 * h * zs / a;
    for (int k = 0; k < K; k++) { w[er[k]] += z[k] / a; w[ec[k]] += z[k] / a; }
  }
  else {
    for (int i = 0; i < D; i++) if (act[i]) w[i] = ti[i];
  }

  return 1;
}


static void sp_all_rho (const double *t, const unsigned char *act,
                        const double *w, double *rho) {
  int D = sp_D;
  size_t p = 0;
  for (int r = 0; r < D - 1; r++)
    for (int c = r + 1; c < D; c++, p++)
      rho[p] = (act[r] && act[c]) ? sp_rho_of(w[r], w[c], t[p]) : NA_REAL;
}


// Draw d from scratch: the variances, the exclusions, and the stored rho.
static void sp_draw_full (int d, int ti_) {

  int     D = sp_D, Kmax = sp_max_excl;
  size_t  P = sp_pairs;
  double *g    = sp_g    + (size_t) ti_ * D * sp_n;
  double *t    = sp_t    + (size_t) ti_ * P;
  double *rho  = sp_rho  + (size_t) ti_ * P;
  double *w    = sp_w    + (size_t) ti_ * D;
  double *wn   = sp_wnew + (size_t) ti_ * D;
  double *tiv  = sp_ti   + (size_t) ti_ * D;
  double *S    = sp_S    + (size_t) ti_ * (Kmax * Kmax + 1);
  double *z    = sp_z    + (size_t) ti_ * (Kmax + 1);
  unsigned char *ex  = sp_ex  + (size_t) ti_ * P;
  unsigned char *act = sp_act + (size_t) ti_ * D;
  int    *deg  = sp_deg  + (size_t) ti_ * D;
  int    *er   = sp_er   + (size_t) ti_ * (Kmax + 1);
  int    *ec   = sp_ec   + (size_t) ti_ * (Kmax + 1);
  int    *excl = sp_excl ? sp_excl + (size_t) d * Kmax : NULL;

  sp_fill(d, g);
  sp_variation(g, t, 0, P);

  for (size_t p = 0; p < P; p++) ex[p] = 0;
  for (int i = 0; i < D; i++) { act[i] = 1; deg[i] = D - 1; w[i] = NA_REAL; }
  int m = D, K = 0, n_ex = 0, stop = SP_STOP_MAX;

  sp_solve(t, act, m, ex, er, ec, 0, tiv, S, z, w);   // never singular
  sp_all_rho(t, act, w, rho);

  for (int it = 0; it < Kmax; it++) {

    // The strongest pair not yet excluded.
    size_t best = P;
    double top  = sp_thresh;
    for (size_t p = 0; p < P; p++) {
      if (ex[p]) continue;
      double v = fabs(rho[p]);   // NA / NaN never compare greater
      if (v > top) { top = v; best = p; }
    }
    if (best == P) { stop = SP_STOP_THRESHOLD; break; }

    int r = pair_row(best, D);
    int c = r + 1 + (int) (best - pair_index(r, r + 1, D));

    ex[best] = 1;
    deg[r]--; deg[c]--;
    int drop_r = deg[r] == 0, drop_c = deg[c] == 0;
    if (drop_r) { act[r] = 0; m--; }
    if (drop_c) { act[c] = 0; m--; }
    if (excl) excl[n_ex] = (int) best;
    n_ex++;

    if (m < 4) { stop = SP_STOP_THREE; break; }

    // The excluded pairs among the components still active, in packed order.
    K = 0;
    {
      size_t p = 0;
      for (int a = 0; a < D - 1; a++)
        for (int b = a + 1; b < D; b++, p++)
          if (ex[p] && act[a] && act[b]) { er[K] = a; ec[K] = b; K++; }
    }

    for (int i = 0; i < D; i++) wn[i] = w[i];
    if (!sp_solve(t, act, m, ex, er, ec, K, tiv, S, z, wn)) {
      // Singular: undo this exclusion and stop.
      ex[best] = 0;
      deg[r]++; deg[c]++;
      if (drop_r) { act[r] = 1; m++; }
      if (drop_c) { act[c] = 1; m++; }
      n_ex--;
      if (excl) excl[n_ex] = NA_INTEGER;
      stop = SP_STOP_SINGULAR;
      break;
    }
    for (int i = 0; i < D; i++) w[i] = act[i] ? wn[i] : NA_REAL;
    sp_all_rho(t, act, w, rho);
  }

  sp_n_excl[d] = n_ex;
  if (sp_stop) sp_stop[d] = stop;
  double *w2 = sp_w2 + (size_t) d * D;
  for (int i = 0; i < D; i++) w2[i] = act[i] ? w[i] : NA_REAL;

  // Store this draw's rho for the pairs of the block (here, the first).
  int nd = sp_n_draws;
  // (NA for the pairs of a component dropped by the last exclusion.)
  int r = pair_row(sp_p0, D);
  int c = r + 1 + (int) (sp_p0 - pair_index(r, r + 1, D));
  for (size_t p = sp_p0; p < sp_p1; p++) {
    sp_store[(p - sp_p0) * nd + d] = (act[r] && act[c]) ? rho[p] : NA_REAL;
    if (++c == D) { r++; c = r + 1; }
  }
}


// Draw d made again, for a later block of pairs: t for those pairs only,
// and rho from the variances kept from the first pass.
static void sp_draw_block (int d, int ti_) {

  int     D  = sp_D, n = sp_n, nd = sp_n_draws;
  double *g  = sp_g + (size_t) ti_ * D * n;
  double *t  = sp_t + (size_t) ti_ * sp_pairs;
  const double *w2 = sp_w2 + (size_t) d * D;

  sp_fill(d, g);
  sp_variation(g, t, sp_p0, sp_p1);

  int r = pair_row(sp_p0, D);
  int c = r + 1 + (int) (sp_p0 - pair_index(r, r + 1, D));
  for (size_t p = sp_p0; p < sp_p1; p++) {
    double wi = w2[r], wj = w2[c];
    sp_store[(p - sp_p0) * nd + d] =
      (ISNA(wi) || ISNA(wj)) ? NA_REAL : sp_rho_of(wi, wj, t[p - sp_p0]);
    if (++c == D) { r++; c = r + 1; }
  }
}


static void *sp_draw_worker (void *arg) {
  int i = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n;
  for (int d = i; d < sp_n_draws; d += n) {
    if (sp_solve_pass) sp_draw_full(d, i);
    else               sp_draw_block(d, i);
  }
  return NULL;
}


// Quickselect: x[k] in place, x[0 .. k-1] <= x[k] <= x[k+1 ..].
static double select_k (double *x, int len, int k) {
  int lo = 0, hi = len - 1;
  while (lo < hi) {
    double pivot = x[(lo + hi) / 2];
    int i = lo, j = hi;
    while (i <= j) {
      while (x[i] < pivot) i++;
      while (x[j] > pivot) j--;
      if (i <= j) { double tmp = x[i]; x[i] = x[j]; x[j] = tmp; i++; j--; }
    }
    if (k <= j)      hi = j;
    else if (k >= i) lo = i;
    else             break;
  }
  return x[k];
}

// The median of the non-missing values of x[0 .. len-1] (R's median():
// the mean of the middle two of an even count), using buf; NA for none.
static double median_of (const double *x, int len, double *buf) {
  int k = 0;
  for (int i = 0; i < len; i++) if (!ISNAN(x[i])) buf[k++] = x[i];
  if (k == 0) return NA_REAL;
  int    h  = k / 2;
  double hi = select_k(buf, k, h);
  if (k % 2) return hi;
  double lo = buf[0];
  for (int i = 1; i < h; i++) if (buf[i] > lo) lo = buf[i];
  return (lo + hi) / 2;
}

static double *sp_buf;   // n_threads x n_draws

static void *sp_median_worker (void *arg) {
  int i = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n, nd = sp_n_draws;
  double *buf = sp_buf + (size_t) i * nd;
  size_t len = sp_p1 - sp_p0;
  for (size_t p = (size_t) i; p < len; p += (size_t) n)
    sp_med[sp_p0 + p] = median_of(sp_store + p * nd, nd, buf);
  return NULL;
}


// One whole SparCC run on sp_cnt (dataset sp_dataset): medians into med,
// each draw's variances into sp_w2. store_all: keep every draw's rho
// (pairs x n_draws, for the details) in `all` instead of blocks.
static void sp_run (double *med, double *store, size_t block, double *all) {

  size_t P  = sp_pairs;
  int    nd = sp_n_draws;
  double work = (double) nd * P * sp_n;
  int    thr  = threads_for(sp_n_thr, work, 1e6);
  if (thr > nd) thr = nd;
  int    thr_med = threads_for(sp_n_thr, (double) P * nd * 4, 1e6);

  sp_med = med;
  for (size_t p0 = 0; p0 < P; p0 += block) {
    sp_p0 = p0;
    sp_p1 = p0 + block < P ? p0 + block : P;
    sp_store = all ? all + p0 * nd : store;
    sp_solve_pass = p0 == 0;
    run_parallel(sp_draw_worker, thr, thr > 1 ? 100 : 1);
    run_parallel(sp_median_worker, thr_med, thr_med > 1 ? 100 : 1);
    if (p0 + block < P && user_interrupted()) { free_all(); error("Interrupted."); } // # nocov
  }
}


//======================================================
// R interface for sparcc().
//
// sexp_given:   n x D x n_draws log abundances, or NULL.
// sexp_block:   pairs per stored block (0: from SP_STORE_MAX).
// sexp_details: also return every draw's rho, exclusions
//               and variances, and the nulls' medians.
//======================================================
SEXP C_sparcc(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_n_draws, SEXP sexp_max_excl,
    SEXP sexp_thresh, SEXP sexp_perms, SEXP sexp_seed, SEXP sexp_given,
    SEXP sexp_block, SEXP sexp_details, SEXP sexp_n_threads ) {

  int n_draws = asInteger(sexp_n_draws);
  int Kmax    = asInteger(sexp_max_excl);
  int B       = asInteger(sexp_perms);
  int details = asLogical(sexp_details);
  int n_thr   = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  init_n_ptrs(30);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n = em->n_samples, D = em->n_otus;
  size_t P = (size_t) D * (D - 1) / 2;
  if ((size_t) Kmax > P) Kmax = (int) P;   // no more pairs to exclude

  sp_n = n; sp_D = D; sp_n_draws = n_draws; sp_max_excl = Kmax;
  sp_thresh  = asReal(sexp_thresh);
  sp_seed    = (uint64_t) (uint32_t) asInteger(sexp_seed);
  sp_pairs   = P;
  sp_given   = isNull(sexp_given) ? NULL : REAL(sexp_given);
  sp_n_thr   = n_thr;
  sp_dataset = 0;

  // The counts, densely, feature by feature.
  double *cnt = (double *) safe_malloc((size_t) D * n * sizeof(double));
  for (size_t k = 0; k < (size_t) D * n; k++) cnt[k] = 0;
  for (int s = 0; s < n; s++)
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++)
      cnt[(size_t) em->otu_vec[e] * n + s] = em->val_vec[e];
  sp_cnt = cnt;

  // Scratch for each thread.
  int T = n_thr < n_draws ? n_thr : n_draws;
  sp_g    = (double *) safe_malloc((size_t) T * D * n * sizeof(double));
  sp_t    = (double *) safe_malloc((size_t) T * P * sizeof(double));
  sp_rho  = (double *) safe_malloc((size_t) T * P * sizeof(double));
  sp_ex   = (unsigned char *) safe_malloc((size_t) T * P + 1);
  sp_w    = (double *) safe_malloc((size_t) T * D * sizeof(double));
  sp_wnew = (double *) safe_malloc((size_t) T * D * sizeof(double));
  sp_ti   = (double *) safe_malloc((size_t) T * D * sizeof(double));
  sp_S    = (double *) safe_malloc((size_t) T * (Kmax * Kmax + 1) * sizeof(double));
  sp_z    = (double *) safe_malloc((size_t) T * (Kmax + 1) * sizeof(double));
  sp_act  = (unsigned char *) safe_malloc((size_t) T * D);
  sp_deg  = (int *) safe_malloc((size_t) T * D * sizeof(int));
  sp_er   = (int *) safe_malloc((size_t) T * (Kmax + 1) * sizeof(int));
  sp_ec   = (int *) safe_malloc((size_t) T * (Kmax + 1) * sizeof(int));
  sp_buf  = (double *) safe_malloc((size_t) (n_thr > T ? n_thr : T) * n_draws * sizeof(double));

  // Blocks of pairs to store.
  size_t block = (size_t) asReal(sexp_block);
  if (block == 0) {
    double fit = SP_STORE_MAX / (8.0 * n_draws);
    block = fit < 1 ? 1 : (size_t) fit;
  }
  if (details || block > P) block = P;

  SEXP sexp_rho  = PROTECT(allocVector(REALSXP, (R_xlen_t) P));
  SEXP sexp_w2   = PROTECT(allocMatrix(REALSXP, D, n_draws));
  SEXP sexp_nex  = PROTECT(allocVector(INTSXP, n_draws));
  SEXP sexp_stop = PROTECT(allocVector(INTSXP, n_draws));
  SEXP sexp_all  = PROTECT(details ? allocMatrix(REALSXP, n_draws, (int) P) : allocVector(REALSXP, 0));
  SEXP sexp_exc  = PROTECT(details ? allocMatrix(INTSXP, Kmax > 0 ? Kmax : 1, n_draws) : allocVector(INTSXP, 0));
  SEXP sexp_pval = PROTECT(allocVector(REALSXP, B > 0 ? (R_xlen_t) P : 0));
  SEXP sexp_null = PROTECT(details && B > 0 ? allocMatrix(REALSXP, (int) P, B) : allocVector(REALSXP, 0));

  double *store = details ? NULL : (double *) safe_malloc((size_t) block * n_draws * sizeof(double) + 1);
  sp_w2     = REAL(sexp_w2);
  sp_n_excl = INTEGER(sexp_nex);
  sp_stop   = INTEGER(sexp_stop);
  sp_excl   = details ? INTEGER(sexp_exc) : NULL;
  if (details) for (R_xlen_t k = 0; k < XLENGTH(sexp_exc); k++) INTEGER(sexp_exc)[k] = NA_INTEGER;

  double *rho = REAL(sexp_rho);
  sp_run(rho, store, block, details ? REAL(sexp_all) : NULL);

  // Variances: each feature's median over the draws.
  SEXP sexp_var = PROTECT(allocVector(REALSXP, D));
  int  n_nonpos = 0, n_draws_nonpos = 0;
  {
    double *col = (double *) safe_malloc((size_t) n_draws * sizeof(double) + 1);
    for (int d = 0; d < n_draws; d++) {
      int any = 0;
      for (int f = 0; f < D; f++) {
        double v = sp_w2[(size_t) d * D + f];
        if (!ISNAN(v) && !(v > 0)) { n_nonpos++; any = 1; }
      }
      n_draws_nonpos += any;
    }
    for (int f = 0; f < D; f++) {
      for (int d = 0; d < n_draws; d++) col[d] = sp_w2[(size_t) d * D + f];
      REAL(sexp_var)[f] = median_of(col, n_draws, sp_buf);
    }
  }
  if (details) {
    // Exclusions as 1-based packed indices.
    int *x = INTEGER(sexp_exc);
    for (R_xlen_t k = 0; k < XLENGTH(sexp_exc); k++) if (x[k] != NA_INTEGER) x[k]++;
  }

  // Permutation test, as FastSpar's: each null dataset resamples every
  // feature's counts across the samples with replacement, independently of
  // the other features (stream SP_PERM_STREAM + (b - 1) D + f, one uniform
  // per sample), then runs SparCC with the same settings and its own draws.
  // (FastSpar's fastspar_bootstrap does this, judged from its outputs:
  // plans/sparcc.md.)
  if (B > 0) {

    double *pv   = REAL(sexp_pval);
    double *hits = (double *) safe_malloc(P * sizeof(double));
    double *medb = details ? NULL : (double *) safe_malloc(P * sizeof(double));
    double *cntb = (double *) safe_malloc((size_t) D * n * sizeof(double));
    double *w2b  = (double *) safe_malloc((size_t) D * n_draws * sizeof(double));
    int    *ex_b = (int *) safe_malloc((size_t) n_draws * sizeof(int));
    double  tol  = sqrt(DBL_EPSILON);
    for (size_t p = 0; p < P; p++) hits[p] = 0;

    // Null runs keep their variances apart from the observed ones.
    sp_w2     = w2b;
    sp_n_excl = ex_b;
    sp_stop   = NULL;
    sp_excl   = NULL;
    sp_given  = NULL;

    for (int b = 1; b <= B; b++) {

      rng_t rng;
      for (int f = 0; f < D; f++) {
        double *dst = cntb + (size_t) f * n;
        const double *src = cnt + (size_t) f * n;
        rng_seed(&rng, sp_seed, SP_PERM_STREAM + (uint64_t) (b - 1) * D + (uint64_t) f);
        for (int s = 0; s < n; s++) {
          int j = (int) (rng_unif(&rng) * n);
          if (j >= n) j = n - 1;
          dst[s] = src[j];
        }
      }

      sp_cnt     = cntb;
      sp_dataset = (uint64_t) b;
      double *mb = details ? REAL(sexp_null) + (size_t) (b - 1) * P : medb;
      // Details keep only the nulls' medians, stored in one block.
      if (!store) store = (double *) safe_malloc(P * n_draws * sizeof(double) + 1);
      sp_run(mb, store, details ? P : block, NULL);

      for (size_t p = 0; p < P; p++)
        if (!ISNAN(rho[p]) && !ISNAN(mb[p]) && fabs(mb[p]) >= fabs(rho[p]) - tol)
          hits[p] += 1;

      if (b < B && user_interrupted()) { free_all(); error("Interrupted."); } // # nocov
    }

    for (size_t p = 0; p < P; p++)
      pv[p] = ISNAN(rho[p]) ? NA_REAL : (hits[p] + 1) / (B + 1);
  }

  const char *nm[] = { "rho", "variances", "p_value", "n_nonpos", "n_draws_nonpos",
                       "draw_variances", "n_exclusions", "draw_rho", "exclusions", "null_rho",
                       "stop" };
  SEXP out   = PROTECT(allocVector(VECSXP, 11));
  SEXP names = PROTECT(allocVector(STRSXP, 11));
  SET_VECTOR_ELT(out, 0, sexp_rho);
  SET_VECTOR_ELT(out, 1, sexp_var);
  SET_VECTOR_ELT(out, 2, sexp_pval);
  SET_VECTOR_ELT(out, 3, ScalarInteger(n_nonpos));
  SET_VECTOR_ELT(out, 4, ScalarInteger(n_draws_nonpos));
  SET_VECTOR_ELT(out, 5, sexp_w2);
  SET_VECTOR_ELT(out, 6, sexp_nex);
  SET_VECTOR_ELT(out, 7, sexp_all);
  SET_VECTOR_ELT(out, 8, sexp_exc);
  SET_VECTOR_ELT(out, 9, sexp_null);
  SET_VECTOR_ELT(out, 10, sexp_stop);
  for (int k = 0; k < 11; k++) SET_STRING_ELT(names, k, mkChar(nm[k]));
  setAttrib(out, R_NamesSymbol, names);

  free_all();
  UNPROTECT(11);
  return out;
}
