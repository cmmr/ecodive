// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Interpolated and extrapolated Hill numbers of orders q = 0, 1, 2, with
// sample coverage, for hill_extrapolation() and hill_standardize().
// Written from the papers (clean room; see plans/hill_extrapolation.md):
//
//   * Rarefaction (m < n) and extrapolation (m > n) of each order, and of
//     sample coverage: Chao et al. (2014), Table 1 (Eqs. 8, 11b, 11c); for
//     Shannon diversity above n, Chao et al. (2023), Table 2;
//     the reference sample's coverage, Eq. 12 and Appendix G.
//   * Asymptotes: Chao1 (bias-corrected, Appendix G), Chao, Wang & Jost's
//     (2013) entropy (their Eq. 10b), and n (n - 1) / sum X (X - 1).
//   * The bootstrap of Chao et al. (2014), Appendix G: a multinomial sample
//     of size n from the estimated assemblage, every estimate recomputed
//     from it. Each replicate b of sample s has its own PCG stream,
//     (seed, s * n_boot + b), so results are the same for any number of
//     threads; standard deviations are combined over fixed blocks of
//     replicates in a fixed order.
//
// Settled from iNEXT 3.0.2's outputs (the plan's black-box protocol):
// between two whole depths, every estimate (and coverage) is the straight
// line between its values at those depths; a coverage below that of one
// observation is taken at depth 1; the extrapolated coverage of f2 = 0
// uses Appendix G's (n - 1)(f1 - 1) / ((n - 1)(f1 - 1) + 2), as the
// reference coverage does.
//
// The C entry point C_hill_ext() takes, for each sample, its "size rows"
// (depths) and its "coverage rows" (coverage levels, and optionally the
// depth of the point estimate), and returns estimates, coverages, the
// asymptotes, the frequency counts, and with n_boot > 0, bootstrap
// standard deviations of all of them.

#include "ecodive.h"
#include <float.h>

#define HX_BLOCK 8    // Replicates per bootstrap task.
#define HX_LOG_MAX 16777216.0  // Largest whole number in the log table.
#define HX_COV_TOL 1e-14       // Coverages this close are equal.

static double *log_tab;        // log_tab[j] = log(j), j <= log_max.
static double  log_max;

static inline double hx_log (double j) {
  return (j <= log_max) ? log_tab[(size_t) j] : log(j);
}


//======================================================
// One sample's summaries.
//======================================================

typedef struct {
  int           K;        // Distinct abundances.
  const double *k, *f;    // Abundances (increasing) and their counts.
  double N, S;            // Total and richness.
  double f1, f2;          // Singletons and doubletons.
  double f0;              // Chao1's undetected richness.
  double Chat;            // Coverage of the whole sample.
  double Bc;              // Extrapolated coverage's base (1 - Chat = f1/N Bc).
  double Hobs;            // Shannon entropy of the sample.
  double S2;              // sum X (X - 1) / (N (N - 1)).
  double Hhat;            // Chao, Wang & Jost's entropy (NAN until needed).
  double D1m;             // Shannon diversity rarefied to N - 1 (NAN until needed).
} hx_t;


static void hx_prepare (hx_t *h, const double *k, const double *f, int K,
                        double N, double S) {
  h->K = K; h->k = k; h->f = f; h->N = N; h->S = S;
  double f1 = 0, f2 = 0, hs = 0, s2 = 0;
  for (int i = 0; i < K; i++) {
    if (k[i] == 1) f1 = f[i];
    if (k[i] == 2) f2 = f[i];
    double p = k[i] / N;
    hs -= f[i] * p * log(p);
    s2 += f[i] * k[i] * (k[i] - 1);
  }
  h->f1   = f1;
  h->f2   = f2;
  h->Hobs = hs;
  h->S2   = (N > 1) ? s2 / (N * (N - 1)) : 0;
  h->f0   = (f2 > 0) ? (N - 1) / N * f1 * f1 / (2 * f2)
                     : (N - 1) / N * f1 * (f1 - 1) / 2;
  if (f1 == 0) {
    h->Bc   = 0;
    h->Chat = 1;
  }
  else {
    h->Bc   = (f2 > 0) ? (N - 1) * f1 / ((N - 1) * f1 + 2 * f2)
                       : (N - 1) * (f1 - 1) / ((N - 1) * (f1 - 1) + 2);
    h->Chat = 1 - f1 / N * h->Bc;
  }
  h->Hhat = NAN;
  h->D1m  = NAN;
}


// Chao, Wang & Jost's (2013) entropy, Eq. 10b. The first term,
// sum_{k=1}^{N-1} 1/k sum_{X <= N-k} X/N C(N-X, k) / C(N-1, k), equals
// sum_{X <= N-1} X/N sum_{j=X}^{N-1} 1/j (the same sum, by abundance), so
// one pass of harmonic sums from N - 1 down serves every abundance. The
// second, f1/N (1 - A)^(1-N) [-log A - sum_{r=1}^{N-1} (1 - A)^r / r],
// is f1/N sum_{t >= 1} x^t / (N - 1 + t) with x = 1 - A: summed directly
// when that converges quickly, else (N A <= 5, where (1 - A)^(1-N) is at
// most e^5) from the closed form.
static double hx_hhat (hx_t *h) {
  if (!isnan(h->Hhat)) return h->Hhat;
  double N = h->N;
  double H = 0;
  if (N > 1) {
    double hs = 0;  // sum_{j=X}^{N-1} 1/j
    double j  = N - 1;
    for (int i = h->K - 1; i >= 0; i--) {
      double X = h->k[i];
      if (X > N - 1) continue;
      for (; j >= X; j -= 1) hs += 1 / j;
      H += h->f[i] * X / N * hs;
    }
  }
  double f1 = h->f1, f2 = h->f2;
  if (f1 > 0) {
    double A = (f2 > 0) ? 2 * f2 / ((N - 1) * f1 + 2 * f2)
                        : 2 / ((N - 1) * (f1 - 1) + 2);
    if (A < 1) {
      double x = 1 - A, a = N - 1, T = 0;
      if (N * A <= 5) {
        double s = 0, xr = 1;
        for (double r = 1; r <= a; r += 1) { xr *= x; s += xr / r; }
        T = exp(-a * log1p(-A)) * (-log(A) - s);
      }
      else {
        double xt = 1;
        for (double t = 1; ; t += 1) {
          xt *= x;
          double term = xt / (a + t);
          T += term;
          if (term < 1e-17 * T || xt == 0) break;
        }
      }
      H += f1 / N * T;
    }
  }
  h->Hhat = H;
  return H;
}


// The expectation of -(K/d) log(K/d) for K hypergeometric (a feature with
// X of the sample's N observations, in a subsample of d), from the pmf
// built outward from its mode by ratios and normalized by its own sum.
// Each tail stops once its terms fall below 1e-17 of the running sum (the
// pmf is log-concave, so its terms only fall away from the mode).
static double hx_hyper_entropy (double N, double X, double d) {
  double lo = d - (N - X); if (lo < 0) lo = 0;
  double hi = (X < d) ? X : d;
  double j0 = floor((d + 1) * (X + 1) / (N + 2));
  if (j0 < lo) j0 = lo;
  if (j0 > hi) j0 = hi;
  double ld = hx_log(d), id = 1 / d;
  #define HX_G(j) (((j) > 0) ? ((j) * id) * (ld - hx_log(j)) : 0)
  double s0 = 1, s1 = HX_G(j0);
  double p  = 1;
  for (double j = j0; j < hi; j += 1) {   // up: P(j+1) / P(j)
    p *= (X - j) * (d - j) / ((j + 1) * (N - X - d + j + 1));
    s0 += p; s1 += p * HX_G(j + 1);
    if (p < 1e-17 * s0) break;
  }
  p = 1;
  for (double j = j0; j > lo; j -= 1) {   // down: P(j-1) / P(j)
    p *= j * (N - X - d + j) / ((X - j + 1) * (d - j + 1));
    s0 += p; s1 += p * HX_G(j - 1);
    if (p < 1e-17 * s0) break;
  }
  #undef HX_G
  return s1 / s0;
}


// The absence probabilities C(N - U[i], d) / C(N, d), as rc_absence()
// (rarefaction_curve.c) computes them, by a running product over j = 0 ..
// U[i] - 1 of (N - d - j) / (N - j), but with the numerators and the
// denominators each multiplied in chunks of 16 and divided once per chunk
// (about as many roundings, a sixteenth of the divisions).
static void hx_absence (double N, double d, const double *U, int u, double *A) {
  double p = 1, k = 0;  // p = A(k)
  int    i = 0;
  for (; i < u; i++) {
    double target = U[i];
    if (target > N - d) break;
    while (k < target && p > 0) {
      double stop = target - k;
      if (stop > 16) stop = 16;
      double num = 1, den = 1;
      for (double j = 0; j < stop; j += 1) {
        num *= N - d - k - j;
        den *= N - k - j;
      }
      p *= num / den;
      k += stop;
    }
    if (!(p > 0)) break;
    A[i] = p;
  }
  for (; i < u; i++) A[i] = 0;
}


// Estimates at a whole depth d (1 <= d <= N): out[0..2] the Hill numbers of
// orders 0, 1, 2 (where need[q]), out[3] the coverage. A: scratch for K.
static void hx_at_int (hx_t *h, double d, const int *need, double *A,
                       double *out) {
  double N = h->N;
  int    K = h->K;
  if (d >= N) {
    out[0] = h->S;
    out[1] = need[1] ? exp(h->Hobs) : 0;
    out[2] = 1 / (1 / N + (N - 1) / N * h->S2);
    out[3] = h->Chat;
    return;
  }
  hx_absence(N, d, h->k, K, A);
  if (need[0]) {
    double e = 0;
    for (int i = 0; i < K; i++) e += h->f[i] * A[i];
    out[0] = (d > 1) ? h->S - e : 1;
  }
  else out[0] = 0;
  if (need[1]) {
    double hs = 0;
    if (d > 1)
      for (int i = 0; i < K; i++) hs += h->f[i] * hx_hyper_entropy(N, h->k[i], d);
    out[1] = exp(hs);
  }
  else out[1] = 0;
  out[2] = 1 / (1 / d + (d - 1) / d * h->S2);
  out[3] = rc_coverage(N, d, K, h->k, h->f, A);
}


// Shannon diversity extrapolated to N + ms (Chao et al. 2023, Table 2, the
// unified form that replaces Chao et al. 2014's Eq. 10c):
//   D(N + m*) = D_obs + (D_asy - D_obs) [1 - (1 - beta)^m*],
// with beta = (D_obs - D(N - 1)) / (D_asy - D(N - 1)), so that the curve
// also passes through the rarefied value at N - 1. For q = 0 this beta is
// Chao et al. (2014)'s f1 / (N f0 + f1), as the expected loss of richness
// from dropping one observation is f1 / N. Settled against iNEXT 3.0.2.
static double hx_q1_extrap (hx_t *h, double ms) {
  double obs = exp(h->Hobs), asy = exp(hx_hhat(h));
  if (!(asy != obs) || !isfinite(asy)) return obs;
  double N = h->N;
  if (isnan(h->D1m)) {
    double hs = 0;
    if (N > 2)
      for (int i = 0; i < h->K; i++) hs += h->f[i] * hx_hyper_entropy(N, h->k[i], N - 1);
    h->D1m = (N >= 2) ? exp(hs) : obs;
  }
  double beta = (obs - h->D1m) / (asy - h->D1m);
  if (!(beta > 0)) return obs;
  if (beta >= 1)   return asy;
  return obs + (asy - obs) * (1 - exp(ms * log1p(-beta)));
}


// Estimates at any depth m >= 1, or Inf (the asymptotes, coverage 1).
static void hx_at (hx_t *h, double m, const int *need, double *A, double *out) {
  double N = h->N;
  if (m < N) {
    double lo = floor(m), fr = m - lo;
    hx_at_int(h, lo, need, A, out);
    if (fr > 0) {
      double hi[4];
      hx_at_int(h, lo + 1, need, A, hi);
      for (int j = 0; j < 4; j++) out[j] = (1 - fr) * out[j] + fr * hi[j];
    }
    return;
  }
  if (m == N) { hx_at_int(h, N, need, A, out); return; }
  if (isinf(m)) {
    out[0] = h->S + h->f0;
    out[1] = need[1] ? exp(hx_hhat(h)) : 0;
    out[2] = (h->S2 > 0) ? 1 / h->S2 : R_PosInf;
    out[3] = 1;
    return;
  }
  double ms = m - N;  // m*
  out[0] = (h->f1 > 0 && h->f0 > 0)
    ? h->S + h->f0 * (1 - exp(ms * log1p(-h->f1 / (N * h->f0 + h->f1))))
    : h->S;
  out[1] = need[1] ? hx_q1_extrap(h, ms) : 0;
  out[2] = 1 / (1 / m + (m - 1) / m * h->S2);
  out[3] = (h->f1 > 0 && h->Bc > 0) ? 1 - h->f1 / N * exp((ms + 1) * log(h->Bc)) : 1;
}


// Coverage at a whole depth 1 <= d <= N.
static double hx_cov_int (hx_t *h, double d, double *A) {
  if (d >= h->N) return h->Chat;
  hx_absence(h->N, d, h->k, h->K, A);
  return rc_coverage(h->N, d, h->K, h->k, h->f, A);
}


// The depth at which the expected coverage is L (0 < L <= 1): Inf for 1
// (beyond a complete sample); in closed form when extrapolating; else
// between the whole depths whose coverages straddle L, on the line between
// them. At most depth 1. [lo, hi] are whole depths known to bracket the
// answer (coverage(lo) < L <= coverage(hi)), or 0 when unknown.
static double hx_depth_for (hx_t *h, double L, double *A, double lo, double hi) {
  double N = h->N;
  if (L >= 1) return (h->Chat >= 1) ? N : R_PosInf;
  if (L > h->Chat) {
    double ms = log(N * (1 - L) / h->f1) / log(h->Bc) - 1;
    return (ms > 0) ? N + ms : N;
  }
  double c1 = hx_cov_int(h, 1, A);
  if (L <= c1 + HX_COV_TOL) return 1;
  if (lo < 1)                       lo = 1;
  if (hi <= lo || hi > N)           hi = N;
  // Smallest whole depth with coverage >= L, in (lo, hi]: by interpolation
  // on the chord between the bracket's ends, with a bisection after any
  // step that moves the same end twice running (coverage is smooth and
  // concave, so plain interpolation can creep from one side). Coverages
  // within HX_COV_TOL of L count as reaching it: where coverage is flat
  // over a range of depths, rounding would otherwise pick any of them.
  L -= HX_COV_TOL;
  double cl = hx_cov_int(h, lo, A), ch = hx_cov_int(h, hi, A);
  int    last = 0;  // which end moved last: -1 lo, +1 hi
  int    same = 0;
  while (hi - lo > 1) {
    double mid;
    if (same >= 2 || !(ch > cl)) {
      mid = floor((lo + hi) / 2);
    }
    else {
      mid = floor(lo + (L - cl) / (ch - cl) * (hi - lo));
      if (mid <= lo) mid = lo + 1;
      if (mid >= hi) mid = hi - 1;
    }
    double cm = hx_cov_int(h, mid, A);
    int    end;
    if (cm >= L) { hi = mid; ch = cm; end = +1; }
    else         { lo = mid; cl = cm; end = -1; }
    same = (end == last) ? same + 1 : 0;
    if (same >= 3) same = 0;
    last = end;
  }
  L += HX_COV_TOL;
  if (ch <= L + HX_COV_TOL) return hi;
  double fr = (L - cl) / (ch - cl);
  return lo + ((fr > 1) ? 1 : fr);
}




//======================================================
// The whole computation.
//======================================================

static int     n_samples, n_boot, n_blocks;
static int     task_off, task_end;   // This batch of bootstrap tasks.
static int     need[3];
static uint64_t hx_seed;

static int    *pos_vec;          // ecomatrix offsets
static double *val_vec;
static double *frq_k, *frq_f, *frq_N, *frq_S;
static int    *frq_K;

static int    *size_pos;         // Sample s's size rows: size_pos[s] .. size_pos[s+1]-1
static double *size_depth;
static int    *cov_pos;          // Its coverage rows.
static double *cov_level, *cov_given;

static double *o_size;           // R x 4: q0, q1, q2, coverage
static double *o_cdepth;         // C: depth of each coverage row
static double *o_cov;            // C x 4: q0, q1, q2, coverage reached
static double *o_asy;            // n x 3
static double *o_info;           // n x 13: N, S, Chat, f1 .. f10

// Bootstrap: one assemblage per sample.
static int    *asm_pos;          // Sample s's probabilities: asm_p[asm_pos[s] ..]
static double *asm_p;
static int     use_hook;         // Replicates from hook_ptr, not drawn.
static double **hook_ptr;        // Sample s's replicate counts (d x n_boot).

static int    *val_pos;          // Sample s's replicate values: val_pos[s] ..
static double *acc_mean, *acc_m2;// n_blocks x total values, by block

static int     scr_K, scr_D;     // Per-thread scratch sizes.
static double *scr_A, *scr_k, *scr_f, *scr_buf, *scr_y;
static int64_t *scr_cnt;
static int     n_rows_total, n_cov_total;


// Every estimate of one sample (or one replicate of it), into vals: size
// rows (4 each), coverage rows (q0..q2 each; also the depths and
// coverages reached when cdepth / creach aren't NULL), asymptotes (3).
static void hx_sample (hx_t *h, int s, int point, double *A, double *vals,
                       double *cdepth, double *creach) {
  int    r0 = size_pos[s], r1 = size_pos[s + 1];
  int    c0 = cov_pos[s],  c1 = cov_pos[s + 1];
  int    nr = r1 - r0;
  double out[4];

  for (int r = r0; r < r1; r++) {
    hx_at(h, size_depth[r], need, A, out);
    for (int j = 0; j < 4; j++) vals[4 * (r - r0) + j] = out[j];
  }

  for (int c = c0; c < c1; c++) {
    double d;
    if (point && !isnan(cov_given[c]) && !isnan(cov_level[c])) {
      d = cov_given[c];
    }
    else if (point && !isnan(cov_given[c]) && cov_given[c] >= h->N) {
      d = cov_given[c];   // the level is this depth's coverage, reached here first
    }
    else {
      // A level not given is the coverage at the given depth (below N),
      // which may be reached at a smaller depth already: solve for it.
      if (isnan(cov_level[c])) {
        hx_at_int(h, cov_given[c], (int[]){0, 0, 0}, A, out);
        cov_level[c] = out[3];
      }
      // Bracket with the size rows' coverages (whole depths below N).
      double L = cov_level[c], lo = 0, hi = 0;
      for (int r = r0; r < r1; r++) {
        double dr = size_depth[r], cr = vals[4 * (r - r0) + 3];
        if (dr != floor(dr) || dr >= h->N) continue;
        if (cr <  L - HX_COV_TOL && (lo == 0 || dr > lo)) lo = dr;
        if (cr >= L - HX_COV_TOL && (hi == 0 || dr < hi)) hi = dr;
      }
      if (hi == 0) hi = h->N;
      d = hx_depth_for(h, L, A, lo, hi);
    }
    hx_at(h, d, need, A, out);
    if (point && isnan(cov_level[c])) cov_level[c] = out[3];
    double *v = vals + 4 * nr + 3 * (c - c0);
    v[0] = out[0]; v[1] = out[1]; v[2] = out[2];
    if (cdepth) cdepth[c] = d;
    if (creach) creach[c] = out[3];
  }

  hx_at(h, R_PosInf, need, A, out);
  double *v = vals + 4 * nr + 3 * (c1 - c0);
  v[0] = out[0]; v[1] = out[1]; v[2] = out[2];
}


static int cmp_desc (const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x < y) - (x > y);
}

static void *point_worker (void *arg) {
  int     t = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n;
  double *A = scr_A + (size_t) t * scr_K;

  for (int s = t; s < n_samples; s += n) {

    double N = frq_N[s], S = frq_S[s];
    int    K = frq_K[s];
    double *k = frq_k + pos_vec[s], *f = frq_f + pos_vec[s];

    double *info = o_info + (size_t) s * 13;
    info[0] = N; info[1] = S;
    for (int j = 3; j < 13; j++) info[j] = 0;
    for (int i = 0; i < K; i++) if (k[i] <= 10) info[2 + (int) k[i]] = f[i];

    int nv = val_pos[s + 1] - val_pos[s];
    int nr = size_pos[s + 1] - size_pos[s];
    int nc = cov_pos[s + 1] - cov_pos[s];

    if (N < 1) {
      info[2] = NA_REAL;
      for (int r = 0; r < nr; r++) for (int j = 0; j < 4; j++) o_size[4 * (size_pos[s] + r) + j] = NA_REAL;
      for (int c = 0; c < nc; c++) {
        for (int j = 0; j < 4; j++) o_cov[4 * (cov_pos[s] + c) + j] = NA_REAL;
        o_cdepth[cov_pos[s] + c] = NA_REAL;
      }
      for (int j = 0; j < 3; j++) o_asy[3 * s + j] = NA_REAL;
      continue;
    }

    hx_t h;
    hx_prepare(&h, k, f, K, N, S);
    info[2] = h.Chat;

    double *vals = (double *) malloc((size_t) (nv > 0 ? nv : 1) * sizeof(double));
    if (!vals) continue;  // # nocov
    hx_sample(&h, s, 1, A, vals, o_cdepth, NULL);

    for (int r = 0; r < nr; r++)
      for (int j = 0; j < 4; j++)
        o_size[4 * (size_pos[s] + r) + j] = vals[4 * r + j];
    for (int c = 0; c < nc; c++) {
      for (int j = 0; j < 3; j++) o_cov[4 * (cov_pos[s] + c) + j] = vals[4 * nr + 3 * c + j];
      double d = o_cdepth[cov_pos[s] + c], out[4];
      hx_at(&h, d, (int[]){0, 0, 0}, A, out);
      o_cov[4 * (cov_pos[s] + c) + 3] = out[3];
    }
    for (int j = 0; j < 3; j++) o_asy[3 * s + j] = vals[4 * nr + 3 * nc + j];

    // The bootstrap assemblage (Appendix G), most abundant first.
    if (n_boot > 0 && !use_hook) {
      double *p   = asm_p + asm_pos[s];
      int     d   = asm_pos[s + 1] - asm_pos[s];
      double  lam = 0, den = 0;
      if (h.Chat < 1) {
        for (int i = 0; i < K; i++) den += f[i] * k[i] / N * pow(1 - k[i] / N, N);
        lam = (den > 0) ? (1 - h.Chat) / den : 0;
      }
      int j = 0;
      for (int i = K - 1; i >= 0; i--) {
        double pi = k[i] / N * (1 - lam * pow(1 - k[i] / N, N));
        if (pi < 0) pi = 0;
        for (int c = 0; c < (int) f[i]; c++) p[j++] = pi;
      }
      int f0u = d - j;
      for (; j < d; j++) p[j] = (1 - h.Chat) / f0u;
      qsort(p, (size_t) d, sizeof(double), cmp_desc);
    }

    free(vals);
  }
  return NULL;
}


static void *boot_worker (void *arg) {
  int      t   = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n;
  double  *A   = scr_A   + (size_t) t * scr_K;
  double  *kk  = scr_k   + (size_t) t * scr_D;
  double  *ff  = scr_f   + (size_t) t * scr_D;
  double  *buf = scr_buf + (size_t) t * scr_D;
  double  *yv  = scr_y   + (size_t) t * scr_D;
  int64_t *cnt = scr_cnt + (size_t) t * scr_D;
  int      total_v = val_pos[n_samples];

  for (int task = task_off + t; task < task_end; task += n) {

    int s   = task / n_blocks, blk = task % n_blocks;
    int nv  = val_pos[s + 1] - val_pos[s];
    double *mean = acc_mean + (size_t) blk * total_v + val_pos[s];
    double *m2   = acc_m2   + (size_t) blk * total_v + val_pos[s];
    for (int j = 0; j < nv; j++) { mean[j] = 0; m2[j] = 0; }
    if (frq_N[s] < 1 || nv == 0) continue;

    int b0 = blk * HX_BLOCK, b1 = b0 + HX_BLOCK;
    if (b1 > n_boot) b1 = n_boot;

    double *vals = (double *) malloc((size_t) nv * sizeof(double));
    if (!vals) continue;  // # nocov

    const double *p = asm_p + asm_pos[s];
    int           d = asm_pos[s + 1] - asm_pos[s];

    for (int b = b0; b < b1; b++) {

      if (use_hook) {
        const double *col = hook_ptr[s] + (size_t) b * d;
        for (int i = 0; i < d; i++) yv[i] = col[i];
      }
      else {
        rng_t rng;
        rng_seed(&rng, hx_seed, (uint64_t) s * (uint64_t) n_boot + (uint64_t) b);
        rng_multinom(&rng, (int64_t) frq_N[s], p, d, cnt);
        for (int i = 0; i < d; i++) yv[i] = (double) cnt[i];
      }

      double NN, SS;
      int KK = rc_freq(yv, d, buf, kk, ff, &NN, &SS);
      hx_t h;
      hx_prepare(&h, kk, ff, KK, NN, SS);
      hx_sample(&h, s, 0, A, vals, NULL, NULL);

      // Welford, in replicate order.
      double cnt_b = (double) (b - b0 + 1);
      for (int j = 0; j < nv; j++) {
        double delta = vals[j] - mean[j];
        mean[j] += delta / cnt_b;
        m2[j]   += delta * (vals[j] - mean[j]);
      }
    }
    free(vals);
  }
  return NULL;
}


// The frequency counts of every sample.
static void *freq_worker (void *arg) {
  int i = ((worker_t *) arg)->i, n = ((worker_t *) arg)->n;
  for (int s = i; s < n_samples; s += n) {
    int p = pos_vec[s], nnz = pos_vec[s + 1] - p;
    frq_K[s] = rc_freq(val_vec + p, nnz, frq_f + p, frq_k + p, frq_f + p, frq_N + s, frq_S + s);
  }
  return NULL;
}


// counts, margin: the table (whole numbers). size_pos, size_depth: n + 1
// offsets and each size row's depth (>= 1, or Inf). cov_pos, cov_level,
// cov_given: n + 1 offsets, each coverage row's level (0, 1], and the
// depth of its point estimate if known (else NaN, to solve for). A NaN
// level is the coverage at the given depth; below N, the depth is then
// solved for too: the smallest whole depth with that coverage (coverage can
// be flat over a range of depths). need:
// which of q = 0, 1, 2 to compute (q = 1 is the costly one). n_boot, seed:
// the bootstrap; boot: NULL, or a list with one numeric matrix per sample
// (rows: the assemblage's species; columns: replicates) to use instead of
// drawing (for testing).
SEXP C_hill_ext (SEXP sexp_counts, SEXP sexp_margin,
                 SEXP sexp_size_pos, SEXP sexp_size_depth,
                 SEXP sexp_cov_pos, SEXP sexp_cov_level, SEXP sexp_cov_given,
                 SEXP sexp_need, SEXP sexp_n_boot, SEXP sexp_seed,
                 SEXP sexp_boot, SEXP sexp_n_threads) {

  init_n_ptrs(30);
  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  n_samples = em->n_samples;
  pos_vec   = em->pos_vec;
  val_vec   = em->val_vec;
  int nnz   = em->nnz;

  frq_k = (double *) safe_malloc((size_t) (nnz + 1) * sizeof(double));
  frq_f = (double *) safe_malloc((size_t) (nnz + 1) * sizeof(double));
  frq_K = (int *)    safe_malloc((size_t) (n_samples + 1) * sizeof(int));
  frq_N = (double *) safe_malloc((size_t) (n_samples + 1) * sizeof(double));
  frq_S = (double *) safe_malloc((size_t) (n_samples + 1) * sizeof(double));
  run_parallel(freq_worker, n_threads, n_samples);

  size_pos   = INTEGER(sexp_size_pos);
  size_depth = REAL(sexp_size_depth);
  cov_pos    = INTEGER(sexp_cov_pos);
  n_cov_total = LENGTH(sexp_cov_level);
  cov_level  = (double *) safe_malloc((size_t) (n_cov_total + 1) * sizeof(double));
  for (int c = 0; c < n_cov_total; c++) cov_level[c] = REAL(sexp_cov_level)[c];
  cov_given  = REAL(sexp_cov_given);
  for (int j = 0; j < 3; j++) need[j] = INTEGER(sexp_need)[j];
  n_boot     = asInteger(sexp_n_boot);
  hx_seed    = (uint64_t) asInteger(sexp_seed);
  use_hook   = (sexp_boot != R_NilValue);
  hook_ptr   = (double **) safe_malloc((size_t) (n_samples + 1) * sizeof(double *));
  n_rows_total = LENGTH(sexp_size_depth);

  if (LENGTH(sexp_size_pos) != n_samples + 1 || size_pos[n_samples] != n_rows_total ||
      LENGTH(sexp_cov_pos)  != n_samples + 1 || cov_pos[n_samples]  != n_cov_total) {
    free_all();                               // # nocov
    error("Rows don't match the samples.");   // # nocov
  }

  // Values per sample, for the bootstrap.
  val_pos = (int *) safe_malloc((size_t) (n_samples + 1) * sizeof(int));
  val_pos[0] = 0;
  for (int s = 0; s < n_samples; s++)
    val_pos[s + 1] = val_pos[s] + 4 * (size_pos[s + 1] - size_pos[s]) +
                     3 * (cov_pos[s + 1] - cov_pos[s]) + 3;

  // The bootstrap assemblages: S + ceil(f0) species, when f0 > 0 and the
  // coverage is short of 1.
  asm_pos = (int *) safe_malloc((size_t) (n_samples + 1) * sizeof(int));
  asm_pos[0] = 0;
  scr_D = 1;
  for (int s = 0; s < n_samples; s++) {
    int d = 0;
    if (n_boot > 0 && frq_N[s] >= 1) {
      if (use_hook) {
        SEXP hm = VECTOR_ELT(sexp_boot, s);
        if (!isReal(hm) || !isMatrix(hm) || ncols(hm) < n_boot) {
          free_all();                                         // # nocov
          error("`boot` must hold numeric matrices.");        // # nocov
        }
        d = nrows(hm);
        hook_ptr[s] = REAL(hm);
      }
      else {
        hx_t h;
        hx_prepare(&h, frq_k + pos_vec[s], frq_f + pos_vec[s], frq_K[s], frq_N[s], frq_S[s]);
        double f0u = (h.f0 > 0 && h.Chat < 1) ? ceil(h.f0) : 0;
        if (h.S + f0u > (double) INT_MAX / 4) {
          free_all();                                                       // # nocov
          error("A sample's estimated undetected richness is too large.");  // # nocov
        }
        d = (int) (h.S + f0u);
      }
    }
    asm_pos[s + 1] = asm_pos[s] + d;
    if (d > scr_D) scr_D = d;
  }
  asm_p = (double *) safe_malloc((size_t) (asm_pos[n_samples] + 1) * sizeof(double));

  // log(j) for the whole numbers up to the largest sample (or 2^24).
  double max_N = 1;
  for (int s = 0; s < n_samples; s++) if (frq_N[s] > max_N) max_N = frq_N[s];
  log_max = (max_N < HX_LOG_MAX) ? max_N : HX_LOG_MAX;
  log_tab = (double *) safe_malloc(((size_t) log_max + 1) * sizeof(double));
  log_tab[0] = R_NegInf;
  for (size_t j = 1; j <= (size_t) log_max; j++) log_tab[j] = log((double) j);

  scr_K = 1;
  for (int s = 0; s < n_samples; s++) if (frq_K[s] > scr_K) scr_K = frq_K[s];
  if (scr_D > scr_K) scr_K = scr_D;
  scr_A = (double *) safe_malloc((size_t) n_threads * scr_K * sizeof(double));

  SEXP sexp_size   = PROTECT(allocMatrix(REALSXP, 4, n_rows_total));
  SEXP sexp_cdepth = PROTECT(allocVector(REALSXP, n_cov_total));
  SEXP sexp_cov    = PROTECT(allocMatrix(REALSXP, 4, n_cov_total));
  SEXP sexp_asy    = PROTECT(allocMatrix(REALSXP, 3, n_samples));
  SEXP sexp_info   = PROTECT(allocMatrix(REALSXP, 13, n_samples));
  o_size   = REAL(sexp_size);
  o_cdepth = REAL(sexp_cdepth);
  o_cov    = REAL(sexp_cov);
  o_asy    = REAL(sexp_asy);
  o_info   = REAL(sexp_info);

  // Without the bootstrap, each sample's curves are worth a thread when
  // they hold many depths, or the samples are many.
  run_parallel(point_worker, n_threads,
               (n_samples > 1 && (n_boot > 0 || n_rows_total + n_cov_total >= 100)) ? 100 : n_samples);

  SEXP sexp_sd = R_NilValue;
  if (n_boot > 0) {
    n_blocks = (n_boot + HX_BLOCK - 1) / HX_BLOCK;
    int total_v = val_pos[n_samples];
    acc_mean = (double *) safe_malloc((size_t) n_blocks * (total_v + 1) * sizeof(double));
    acc_m2   = (double *) safe_malloc((size_t) n_blocks * (total_v + 1) * sizeof(double));
    scr_k    = (double *) safe_malloc((size_t) n_threads * scr_D * sizeof(double));
    scr_f    = (double *) safe_malloc((size_t) n_threads * scr_D * sizeof(double));
    scr_buf  = (double *) safe_malloc((size_t) n_threads * scr_D * sizeof(double));
    scr_y    = (double *) safe_malloc((size_t) n_threads * scr_D * sizeof(double));
    scr_cnt  = (int64_t *) safe_malloc((size_t) n_threads * scr_D * sizeof(int64_t));

    // In batches, so that an interrupt is seen between them.
    int n_tasks = n_samples * n_blocks, batch = 16 * n_threads;
    for (task_off = 0; task_off < n_tasks; task_off += batch) {
      task_end = task_off + batch;
      if (task_end > n_tasks) task_end = n_tasks;
      // (Each task is a block of replicates: worth a thread of its own.)
      run_parallel(boot_worker, n_threads, (task_end - task_off > 1) ? 100 : 1);
      if (task_end < n_tasks && user_interrupted()) {
        free_all(); UNPROTECT(5); error("Interrupted."); // # nocov
      }
    }

    // Combine the blocks, in block order (Chan et al.'s pairwise update).
    sexp_sd = PROTECT(allocVector(REALSXP, total_v));
    double *sd = REAL(sexp_sd);
    for (int v = 0; v < total_v; v++) {
      double cn = 0, mu = 0, M2 = 0;
      for (int blk = 0; blk < n_blocks; blk++) {
        double nb = (blk == n_blocks - 1) ? n_boot - blk * HX_BLOCK : HX_BLOCK;
        double mb = acc_mean[(size_t) blk * total_v + v];
        double qb = acc_m2[(size_t) blk * total_v + v];
        if (cn == 0) { cn = nb; mu = mb; M2 = qb; continue; }
        double delta = mb - mu, tot = cn + nb;
        mu += delta * nb / tot;
        M2 += qb + delta * delta * cn * nb / tot;
        cn  = tot;
      }
      sd[v] = (n_boot > 1) ? sqrt(M2 / (n_boot - 1)) : NA_REAL;
    }
    // Empty samples have no standard deviations.
    for (int s = 0; s < n_samples; s++)
      if (frq_N[s] < 1)
        for (int v = val_pos[s]; v < val_pos[s + 1]; v++) sd[v] = NA_REAL;
  }
  else {
    sexp_sd = PROTECT(allocVector(REALSXP, 0));
  }

  SEXP sexp_vpos = PROTECT(allocVector(INTSXP, n_samples + 1));
  for (int s = 0; s <= n_samples; s++) INTEGER(sexp_vpos)[s] = val_pos[s];

  SEXP sexp_clevel = PROTECT(allocVector(REALSXP, n_cov_total));
  for (int c = 0; c < n_cov_total; c++) REAL(sexp_clevel)[c] = cov_level[c];

  const char *nms[] = {"size", "cov_depth", "cov", "asy", "info", "sd", "val_pos", "cov_level"};
  SEXP sexp_result = PROTECT(allocVector(VECSXP, 8));
  SEXP sexp_names  = PROTECT(allocVector(STRSXP, 8));
  SET_VECTOR_ELT(sexp_result, 0, sexp_size);
  SET_VECTOR_ELT(sexp_result, 1, sexp_cdepth);
  SET_VECTOR_ELT(sexp_result, 2, sexp_cov);
  SET_VECTOR_ELT(sexp_result, 3, sexp_asy);
  SET_VECTOR_ELT(sexp_result, 4, sexp_info);
  SET_VECTOR_ELT(sexp_result, 5, sexp_sd);
  SET_VECTOR_ELT(sexp_result, 6, sexp_vpos);
  SET_VECTOR_ELT(sexp_result, 7, sexp_clevel);
  for (int j = 0; j < 8; j++) SET_STRING_ELT(sexp_names, j, mkChar(nms[j]));
  setAttrib(sexp_result, R_NamesSymbol, sexp_names);

  free_all();
  UNPROTECT(10);
  return sexp_result;
}


// Binomial draws, for the sampler's tests: n_draws values of Bin(n, p).
SEXP C_rng_binom (SEXP sexp_n_draws, SEXP sexp_n, SEXP sexp_p, SEXP sexp_seed) {
  int  nd = asInteger(sexp_n_draws);
  SEXP out = PROTECT(allocVector(REALSXP, nd));
  rng_t rng;
  rng_seed(&rng, (uint64_t) asInteger(sexp_seed), 0);
  int64_t n = (int64_t) asReal(sexp_n);
  double  p = asReal(sexp_p);
  for (int i = 0; i < nd; i++) REAL(out)[i] = (double) rng_binom(&rng, n, p);
  UNPROTECT(1);
  return out;
}


// Multinomial draws, for the sampler's tests: n_draws columns of
// Multinom(n, p), each from its own stream (seed, draw).
SEXP C_rng_multinom (SEXP sexp_n_draws, SEXP sexp_n, SEXP sexp_p, SEXP sexp_seed) {
  int     nd = asInteger(sexp_n_draws), d = LENGTH(sexp_p);
  SEXP    out = PROTECT(allocMatrix(REALSXP, d, nd));
  int64_t *cnt = (int64_t *) R_alloc((size_t) (d > 0 ? d : 1), sizeof(int64_t));
  for (int b = 0; b < nd; b++) {
    rng_t rng;
    rng_seed(&rng, (uint64_t) asInteger(sexp_seed), (uint64_t) b);
    rng_multinom(&rng, (int64_t) asReal(sexp_n), REAL(sexp_p), d, cnt);
    for (int i = 0; i < d; i++) REAL(out)[(size_t) b * d + i] = (double) cnt[i];
  }
  UNPROTECT(1);
  return out;
}
