// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Dirichlet-multinomial mixtures (dmm()), written from Holmes, Harris &
// Quince (2012) (clean room). DirichletMultinomial (Bioconductor, LGPL) and
// mothur (GPL) were only run, to compare outputs; their source was not read.
// The prior's constants and the form of the reported criteria were pinned
// from DirichletMultinomial's outputs (plans/community_typing.md).
//
// Sample i has counts x_ij (total J_i) over S features. Type k has Dirichlet
// parameters alpha_kj = exp(lambda_kj), with A_k = sum_j alpha_kj, and weight
// pi_k. The log likelihood of sample i under type k is
//
//   l_ik = log(J_i! / prod_j x_ij!) + lgamma(A_k) - lgamma(A_k + J_i)
//          + sum_j [lgamma(alpha_kj + x_ij) - lgamma(alpha_kj)],
//
// where only the features with x_ij > 0 add to the sum, so the work is
// O(non-zeros). Each lambda_kj has the prior Gamma(eta, nu) on alpha_kj,
// carried to lambda: log p = eta log nu - lgamma(eta) + eta lambda - nu alpha.
//
// EM: z_ik = pi_k exp(l_ik) / sum_k pi_k exp(l_ik) (E-step); then for each
// k, lambda_k minimises
//
//   f(lambda) = -sum_i z_ik [l_ik without the constant] + nu sum_j alpha_j
//               - eta sum_j lambda_j,
//
// with gradient
//
//   df/dlambda_j = -alpha_j t_j + nu alpha_j - eta,
//   t_j = sum_i z_ik [psi(alpha_j + x_ij) - psi(alpha_j) - psi(A + J_i) + psi(A)],
//
// and Hessian
//
//   H = diag(d) + c alpha alpha',
//   d_j = -alpha_j t_j - alpha_j^2 u_j + nu alpha_j,
//   u_j = sum_i z_ik [psi1(alpha_j + x_ij) - psi1(alpha_j)],
//   c   = sum_i z_ik [psi1(A + J_i) - psi1(A)];
//
// and pi_k = sum_i z_ik / n. EM stops when the negative log posterior,
// NLE = -sum_i log sum_k pi_k exp(l_ik) - log prior, changes by less than
// tol times its value.
//
// The M-step: Holmes et al. used GSL's BFGS. A diagonal plus a rank-one
// matrix is solved in O(S) by the Sherman-Morrison formula, so ecodive
// takes Newton steps instead (safeguarded; see newton()), which reach the
// optimum to rounding in a few steps from EM's previous estimate.
//
// The Laplace approximation needs log|H| of the Hessian of the negative log
// posterior at the optimum: block-diagonal, one H as above per type, so
// log|H_k| = sum_j log d_j + log(1 + c sum_j alpha_j^2 / d_j) in O(S), with
// no S x S matrix. The weights add the (K - 1)-dimensional block of the
// multinomial in pi (with pi_K = 1 - the rest), whose log determinant at
// pi_k = sum_i z_ik / n is (K - 1) log n - sum_k log pi_k.
//
// Threads: the E-step is split over samples; each M-step evaluation over
// features (for the per-feature sums, from a by-feature copy of the table)
// and samples (for the per-sample terms), each written to its own slot and
// then summed in a fixed order on the calling thread. The k-means start is
// split over samples the same way. Results are identical for any number of
// threads. Random starts draw from PCG stream (seed, start), never R's RNG.
//
// Holmes, I., Harris, K., & Quince, C. (2012). Dirichlet multinomial
// mixtures: generative models for microbial metagenomics. PLoS ONE, 7,
// e30126. doi:10.1371/journal.pone.0030126

#include "ecodive.h"
#include <Rmath.h> // lgammafn: no global state for positive arguments


//======================================================
// digamma and trigamma for x > 0: the recurrences
// psi(x) = psi(x + 1) - 1/x and psi1(x) = psi1(x + 1) +
// 1/x^2 up to x >= 10, then the asymptotic series
// (Abramowitz & Stegun 6.3.18 and 6.4.12), to about
// 1e-16 relative. No state, so threads may call them.
//======================================================

static double dmm_digamma (double x) {
  double r = 0;
  while (x < 10) { r -= 1 / x; x += 1; }
  double f = 1 / (x * x);
  double t = f * (-1.0 / 12 + f * (1.0 / 120 + f * (-1.0 / 252 + f * (1.0 / 240 +
             f * (-1.0 / 132 + f * (691.0 / 32760 + f * (-1.0 / 12)))))));
  return r + log(x) - 0.5 / x + t;
}

static double dmm_trigamma (double x) {
  double r = 0;
  while (x < 10) { r += 1 / (x * x); x += 1; }
  double f = 1 / (x * x);
  double t = f * (1.0 / 6 + f * (-1.0 / 30 + f * (1.0 / 42 + f * (-1.0 / 30 +
             f * (5.0 / 66 + f * (-691.0 / 2730 + f * (7.0 / 6)))))));
  return r + (1 + t) / x + 0.5 * f;
}



//======================================================
// The table, by sample (CSR) and by feature (CSC).
//======================================================

static int     dm_n, dm_S;
static int    *dm_pos, *dm_otu;   // by sample
static double *dm_val;
static int    *dm_fpos, *dm_fsam; // by feature
static double *dm_fval;
static double *dm_J;              // each sample's total
static double *dm_mc;             // log(J! / prod x!)
static double  dm_eta, dm_nu;

// The current fit.
static int     dm_K;
static double *dm_lam;   // K x S, row k at dm_lam + k * S
static double *dm_alpha; // K x S
static double *dm_pi;    // K
static double *dm_z;     // n x K, column-major
static double *dm_lse;   // n: log sum_k pi_k exp(l_ik)
static double *dm_lga;   // K x S: lgamma(alpha)

// One M-step evaluation, for type dm_k.
static int     dm_k;
static double *dm_ev_a;  // S: alpha
static double  dm_ev_A;
static double *dm_fv, *dm_fg, *dm_fh; // S: per-feature sums
static double *dm_sv, *dm_sg, *dm_sh; // n: per-sample terms
static int     dm_hess;               // 1: also the Hessian's sums

// k-means.
static double *dm_cent;  // K x S centres
static double *dm_cc;    // K: |centre|^2
static double *dm_pp;    // n: |p_i|^2
static double *dm_dist;  // n: distance to nearest centre
static int    *dm_asg;   // n: nearest centre


//======================================================
// E-step: z, lse. Over samples.
//======================================================

static void *estep_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = dm_n, S = dm_S, K = dm_K;

  double *ll = (double *) malloc((size_t) (K + 1) * sizeof(double));
  if (ll == NULL) return NULL; // # nocov

  for (int i = thread_i; i < n; i += n_threads) {

    double mx = R_NegInf;

    for (int k = 0; k < K; k++) {
      if (!(dm_pi[k] > 0)) { ll[k] = R_NegInf; continue; }
      const double *a  = dm_alpha + (size_t) k * S;
      const double *lg = dm_lga   + (size_t) k * S;
      double A = 0;
      for (int j = 0; j < S; j++) A += a[j];
      double s = lgammafn(A) - lgammafn(A + dm_J[i]);
      for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++) {
        int j = dm_otu[p];
        s += lgammafn(a[j] + dm_val[p]) - lg[j];
      }
      ll[k] = s + log(dm_pi[k]);
      if (ll[k] > mx) mx = ll[k];
    }

    double tot = 0;
    for (int k = 0; k < K; k++) tot += exp(ll[k] - mx);
    double lse = mx + log(tot);
    dm_lse[i] = lse;
    for (int k = 0; k < K; k++) dm_z[i + (size_t) n * k] = exp(ll[k] - lse);
  }

  free(ll);
  return NULL;
}


// The sum of every alpha_k (K values), recomputed in a fixed order.
static double alpha_sum (int k) {
  const double *a = dm_alpha + (size_t) k * dm_S;
  double A = 0;
  for (int j = 0; j < dm_S; j++) A += a[j];
  return A;
}


// The log prior of the current fit.
static double log_prior (void) {
  double lp = 0, c = dm_eta * log(dm_nu) - lgammafn(dm_eta);
  for (int k = 0; k < dm_K; k++) {
    const double *l = dm_lam + (size_t) k * dm_S, *a = dm_alpha + (size_t) k * dm_S;
    for (int j = 0; j < dm_S; j++) lp += c + dm_eta * l[j] - dm_nu * a[j];
  }
  return lp;
}


// E-step, then the negative log posterior.
static double estep (int n_threads) {
  int K = dm_K, S = dm_S;
  for (size_t q = 0; q < (size_t) K * S; q++) dm_lga[q] = lgammafn(dm_alpha[q]);
  int nt = threads_for(n_threads, (double) (dm_pos[dm_n] + dm_n) * K * 40, 1e5);
  run_parallel(estep_worker, nt, dm_n < 100 ? 100 : dm_n);
  double ll = 0;
  for (int i = 0; i < dm_n; i++) ll += dm_lse[i] + dm_mc[i];
  return -ll - log_prior();
}



//======================================================
// M-step: the objective of type dm_k, its gradient and
// (with dm_hess) the pieces of its Hessian. Tasks
// 0 .. S - 1 are features; S .. S + n - 1 are samples.
//======================================================

static void *mstep_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = dm_n, S = dm_S, hess = dm_hess;
  const double *z = dm_z + (size_t) n * dm_k;
  const double *a = dm_ev_a;
  double A = dm_ev_A;

  for (int t = thread_i; t < S + n; t += n_threads) {

    if (t < S) {
      int j = t;
      double aj = a[j], v = 0, g = 0, h = 0;
      double lga = lgammafn(aj), dga = dmm_digamma(aj);
      double tga = hess ? dmm_trigamma(aj) : 0;
      for (int p = dm_fpos[j]; p < dm_fpos[j + 1]; p++) {
        double w = z[dm_fsam[p]];
        if (w == 0) continue;
        double c = aj + dm_fval[p];
        v += w * (lgammafn(c) - lga);
        g += w * (dmm_digamma(c) - dga);
        if (hess) h += w * (dmm_trigamma(c) - tga);
      }
      dm_fv[j] = v;
      dm_fg[j] = g;
      dm_fh[j] = h;
    }

    else {
      int i = t - S;
      double w = z[i];
      if (w == 0) { dm_sv[i] = dm_sg[i] = dm_sh[i] = 0; continue; }
      double C = A + dm_J[i];
      dm_sv[i] = w * (lgammafn(A) - lgammafn(C));
      dm_sg[i] = w * (dmm_digamma(C) - dmm_digamma(A));
      dm_sh[i] = hess ? w * (dmm_trigamma(C) - dmm_trigamma(A)) : 0;
    }
  }

  return NULL;
}


static int dm_n_threads;

// f(lambda) of type dm_k; its gradient in grad; with `hess`, the Hessian
// diag(d) + c alpha alpha' in d and *c. Sums run in a fixed order.
static double mstep_eval_full (const double *lam, double *grad, int hess, double *d, double *c) {

  int S = dm_S, n = dm_n;
  double *a = dm_alpha + (size_t) dm_k * S;
  double A = 0;
  for (int j = 0; j < S; j++) { a[j] = exp(lam[j]); A += a[j]; }
  memcpy(dm_ev_a, a, (size_t) S * sizeof(double));
  dm_ev_A = A;
  dm_hess = hess;

  int nt = threads_for(dm_n_threads, (double) (dm_pos[n] + n) * (hess ? 90 : 60), 1e5);
  run_parallel(mstep_worker, nt, S + n < 100 ? 100 : S + n);

  double ll = 0, gA = 0, hA = 0;
  for (int i = 0; i < n; i++) { ll += dm_sv[i]; gA += dm_sg[i]; hA += dm_sh[i]; }
  for (int j = 0; j < S; j++) ll += dm_fv[j];

  double f = -ll;
  for (int j = 0; j < S; j++) {
    f += dm_nu * a[j] - dm_eta * lam[j];
    double tj = dm_fg[j] - gA;
    grad[j] = -a[j] * tj + dm_nu * a[j] - dm_eta;
    if (hess) d[j] = -a[j] * tj - a[j] * a[j] * dm_fh[j] + dm_nu * a[j];
  }
  if (hess) *c = hA;
  return f;
}


// Newton's method on type dm_k's lambda, in place. The Hessian is
// diag(d) + c alpha alpha', so each step solves in O(S) by the
// Sherman-Morrison formula. Where d_j <= 0 it uses |d_j| (or 1), and drops
// the rank-one term if it would leave the matrix indefinite, so the step is
// always a descent direction; steps are halved until f decreases enough
// (Armijo), and capped at 5 in any lambda. Stops when the predicted decrease
// (the Newton decrement) is below DM_NEWTON_TOL times |f|, after taking
// that last step. Returns 0, or 1 if no step could decrease f.
#define DM_NEWTON_TOL 1e-13
#define DM_NEWTON_MAX 200

static double *dm_nw;  // 6 x S scratch

static int newton (double *lam) {

  int S = dm_S;
  double *g  = dm_nw,         *d  = g + S,      *p = d + S;
  double *l2 = p + S,         *g2 = l2 + S,     *d2 = g2 + S;
  double  c, c2, f = mstep_eval_full(lam, g, 1, d, &c);
  double *a  = dm_alpha + (size_t) dm_k * S;

  for (int it = 0; it < DM_NEWTON_MAX; it++) {

    // The (modified) Newton direction p = -H^-1 g.
    double q = 0, w = 0;
    for (int j = 0; j < S; j++) {
      double dj = d[j] > 0 ? d[j] : (d[j] < 0 ? -d[j] : 1);
      d[j] = dj;
      q += a[j] * a[j] / dj;
      w += a[j] * g[j] / dj;
    }
    double cc = (1 + c * q > 1e-8) ? c : 0;
    double coef = cc * w / (1 + cc * q);
    double gp = 0, big = 0;
    for (int j = 0; j < S; j++) {
      p[j] = -(g[j] - coef * a[j]) / d[j];
      gp  += g[j] * p[j];
      if (fabs(p[j]) > big) big = fabs(p[j]);
    }
    if (!(gp < 0)) return 0; // at a stationary point, to rounding
    if (big > 5) { double r = 5 / big; for (int j = 0; j < S; j++) p[j] *= r; gp *= r; }

    double dec  = -gp / 2;
    int    last = dec <= DM_NEWTON_TOL * fabs(f);

    // Halve the step until f decreases enough.
    double step = 1, f2 = 0;
    int    ok = 0;
    for (int h = 0; h < 40; h++) {
      for (int j = 0; j < S; j++) l2[j] = lam[j] + step * p[j];
      f2 = mstep_eval_full(l2, g2, 1, d2, &c2);
      if (f2 <= f + 1e-4 * step * gp || (last && step == 1)) { ok = 1; break; }
      step /= 2;
    }
    if (!ok) {
      mstep_eval_full(lam, g, 0, NULL, NULL); // alpha back to lambda's
      return 1;
    }

    memcpy(lam, l2, (size_t) S * sizeof(double));
    memcpy(g,   g2, (size_t) S * sizeof(double));
    memcpy(d,   d2, (size_t) S * sizeof(double));
    c = c2;
    f = f2;
    if (last) break;
  }

  return 0;
}


static void mstep (void) {
  int S = dm_S;
  for (int k = 0; k < dm_K; k++) {
    dm_k = k;
    double *lam = dm_lam + (size_t) k * S;
    if (dm_pi[k] > 0) newton(lam);
    double *a = dm_alpha + (size_t) k * S;
    for (int j = 0; j < S; j++) a[j] = exp(lam[j]);
  }
}


// log|H| of type k's block (see the top of this file); NaN if H is not
// positive definite.
static double block_log_det (int k) {
  int S = dm_S;
  dm_k = k;
  double *g = dm_nw, *d = g + S, c;
  mstep_eval_full(dm_lam + (size_t) k * S, g, 1, d, &c);
  const double *a = dm_alpha + (size_t) k * S;
  double ld = 0, q = 0;
  for (int j = 0; j < S; j++) {
    if (!(d[j] > 0)) return R_NaN;
    ld += log(d[j]);
    q  += a[j] * a[j] / d[j];
  }
  double r = 1 + c * q;
  if (!(r > 0)) return R_NaN;
  return ld + log(r);
}



//======================================================
// The start: k-means++ seeding (Arthur & Vassilvitskii
// 2007) on the samples' proportions, then Lloyd's
// iterations; types are k-means clusters.
//======================================================

static void *assign_worker (void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int S = dm_S, K = dm_K;
  for (int i = thread_i; i < dm_n; i += n_threads) {
    int best = 0;
    double bd = R_PosInf;
    for (int k = 0; k < K; k++) {
      const double *c = dm_cent + (size_t) k * S;
      double dot = 0;
      for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++)
        dot += dm_val[p] / dm_J[i] * c[dm_otu[p]];
      double d = dm_pp[i] + dm_cc[k] - 2 * dot;
      if (d < 0) d = 0;
      if (d < bd) { bd = d; best = k; }
    }
    dm_asg[i]  = best;
    dm_dist[i] = bd;
  }
  return NULL;
}


static void set_centre (int k, int i) {
  double *c = dm_cent + (size_t) dm_S * k;
  for (int j = 0; j < dm_S; j++) c[j] = 0;
  for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++) c[dm_otu[p]] = dm_val[p] / dm_J[i];
  dm_cc[k] = dm_pp[i];
}


static void centre_norms (int K) {
  for (int k = 0; k < K; k++) {
    const double *c = dm_cent + (size_t) dm_S * k;
    double s = 0;
    for (int j = 0; j < dm_S; j++) s += c[j] * c[j];
    dm_cc[k] = s;
  }
}


// Distances of every sample to its nearest centre among the first K.
static void assign (int K, int n_threads) {
  int keep = dm_K;
  dm_K = K;
  int nt = threads_for(n_threads, (double) (dm_pos[dm_n] + dm_n) * K, 1e5);
  run_parallel(assign_worker, nt, dm_n < 100 ? 100 : dm_n);
  dm_K = keep;
}


// One k-means run: dm_asg gets the clusters. Returns the within-cluster
// sum of squares, or -1 if some cluster could not be filled (fewer
// distinct samples than types).
static double kmeans_once (rng_t *rng, int n_threads, int *size, int *prev) {

  int n = dm_n, S = dm_S, K = dm_K;

  set_centre(0, (int) (rng_unif(rng) * n));
  for (int k = 1; k < K; k++) {
    assign(k, n_threads);
    double tot = 0;
    for (int i = 0; i < n; i++) tot += dm_dist[i];
    if (!(tot > 0)) return -1;
    double u = rng_unif(rng) * tot, cum = 0;
    int pick = -1;
    for (int i = 0; i < n; i++) {
      if (dm_dist[i] <= 0) continue;
      pick = i;
      cum += dm_dist[i];
      if (u < cum) break;
    }
    set_centre(k, pick);
  }

  for (int i = 0; i < n; i++) prev[i] = -1;

  int ok = 1;
  for (int iter = 0; iter < 100; iter++) {

    assign(K, n_threads);

    // Empty clusters take the sample farthest from its centre (the first
    // on ties), from a cluster of more than one.
    for (int k = 0; k < K; k++) size[k] = 0;
    for (int i = 0; i < n; i++) size[dm_asg[i]]++;
    for (int k = 0; k < K; k++) {
      if (size[k] > 0) continue;
      int far = -1;
      for (int i = 0; i < n; i++)
        if (size[dm_asg[i]] > 1 && dm_dist[i] > 0 && (far < 0 || dm_dist[i] > dm_dist[far]))
          far = i;
      if (far < 0) { ok = 0; break; }
      size[dm_asg[far]]--;
      dm_asg[far]  = k;
      dm_dist[far] = 0;
      size[k] = 1;
    }
    if (!ok) break;

    int changed = 0;
    for (int i = 0; i < n; i++) if (dm_asg[i] != prev[i]) { changed = 1; prev[i] = dm_asg[i]; }
    if (!changed) break;

    // The clusters' mean proportions.
    for (size_t q = 0; q < (size_t) K * S; q++) dm_cent[q] = 0;
    for (int i = 0; i < n; i++) {
      double *c = dm_cent + (size_t) S * dm_asg[i];
      for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++) c[dm_otu[p]] += dm_val[p] / dm_J[i];
    }
    for (int k = 0; k < K; k++) {
      double *c = dm_cent + (size_t) S * k;
      for (int j = 0; j < S; j++) c[j] /= size[k];
    }
    centre_norms(K);
  }

  if (!ok) return -1;
  double ss = 0;
  for (int i = 0; i < n; i++) ss += dm_dist[i];
  return ss;
}


// Fills dm_z with the hard clustering of the best (smallest sum of
// squares) of KM_RUNS k-means runs, drawn in turn from `rng`. Returns 0 if
// no run could fill every cluster.
#define KM_RUNS 10

static int kmeans_start (rng_t *rng, int n_threads) {

  int n = dm_n, K = dm_K;
  int *size = (int *) safe_malloc((size_t) (K + 1) * sizeof(int));
  int *prev = (int *) safe_malloc((size_t) (n + 1) * sizeof(int));
  int *best = (int *) safe_malloc((size_t) (n + 1) * sizeof(int));
  double best_ss = -1;

  for (int r = 0; r < KM_RUNS; r++) {
    double ss = kmeans_once(rng, n_threads, size, prev);
    if (ss >= 0 && (best_ss < 0 || ss < best_ss)) {
      best_ss = ss;
      memcpy(best, dm_asg, (size_t) n * sizeof(int));
    }
  }

  if (best_ss >= 0)
    for (int k = 0; k < K; k++)
      for (int i = 0; i < n; i++)
        dm_z[i + (size_t) n * k] = best[i] == k ? 1 : 0;

  free_one(size);
  free_one(prev);
  free_one(best);
  return best_ss >= 0;
}


// Starting lambda for each type from moment estimates of the weighted
// samples' proportions, and pi from the weights.
static void moment_start (void) {

  int n = dm_n, S = dm_S, K = dm_K;
  double *m = (double *) safe_malloc((size_t) (2 * S + 1) * sizeof(double));
  double *v = m + S;

  for (int k = 0; k < K; k++) {

    const double *z = dm_z + (size_t) n * k;
    double W = 0;
    for (int j = 0; j < S; j++) m[j] = v[j] = 0;
    for (int i = 0; i < n; i++) {
      if (z[i] == 0) continue;
      W += z[i];
      for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++) {
        double q = dm_val[p] / dm_J[i];
        m[dm_otu[p]] += z[i] * q;
        v[dm_otu[p]] += z[i] * q * q;
      }
    }
    dm_pi[k] = W / n;

    // The precision: the geometric mean over features of
    // m (1 - m) / var - 1, where positive; else 1.
    double sl = 0;
    int    nl = 0;
    for (int j = 0; j < S; j++) {
      if (W > 0) { m[j] /= W; v[j] = v[j] / W - m[j] * m[j]; }
      if (m[j] > 0 && m[j] < 1 && v[j] > 0) {
        double r = m[j] * (1 - m[j]) / v[j] - 1;
        if (r > 0) { sl += log(r); nl++; }
      }
    }
    double A = nl ? exp(sl / nl) : 1;
    if (A < 1e-2) A = 1e-2;
    if (A > 1e6)  A = 1e6;

    double *l = dm_lam + (size_t) k * S, *a = dm_alpha + (size_t) k * S;
    for (int j = 0; j < S; j++) {
      double aj = A * m[j];
      if (aj < 1e-6) aj = 1e-6;
      l[j] = log(aj);
      a[j] = aj;
    }
  }

  free_one(m);
}



//======================================================
// One fit of K types from the current dm_z.
// Returns the number of EM iterations; *nle the final
// negative log posterior; *conv 1 if converged.
//======================================================

static int em_fit (int max_iters, double tol, int n_threads, double *nle, int *conv) {
  moment_start();
  double prev = R_PosInf, cur = R_PosInf;
  int it = 0;
  *conv = 0;
  while (it < max_iters) {
    it++;
    mstep();
    for (int k = 0; k < dm_K; k++) {
      double W = 0;
      for (int i = 0; i < dm_n; i++) W += dm_z[i + (size_t) dm_n * k];
      dm_pi[k] = W / dm_n;
    }
    cur = estep(n_threads);
    if (fabs(prev - cur) <= tol * fabs(cur)) { *conv = 1; break; }
    prev = cur;
    if (user_interrupted()) { *conv = -1; break; }
  }
  *nle = cur;
  return it;
}



//======================================================
// R interface.
//
// sexp_types:  the numbers of types to fit (integers).
// sexp_init:   NULL, or an n x K matrix of starting
//              weights z (for a single K; replaces the
//              k-means starts).
//
// Returns, for each number of types, a list(status,
// nle, log_det, n_iters, converged, start, alpha (K x S),
// pi, z (n x K)); types in decreasing pi. status 0: no
// start could fill K clusters.
//======================================================

SEXP C_dmm (SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_types,
            SEXP sexp_n_starts, SEXP sexp_max_iters, SEXP sexp_tol,
            SEXP sexp_seed, SEXP sexp_eta, SEXP sexp_nu, SEXP sexp_init,
            SEXP sexp_n_threads) {

  init_n_ptrs(40);

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int    n_threads = asInteger(sexp_n_threads);
  int    n_starts  = asInteger(sexp_n_starts);
  int    max_iters = asInteger(sexp_max_iters);
  double tol       = asReal(sexp_tol);
  uint64_t seed    = (uint64_t) (uint32_t) asInteger(sexp_seed);
  int    n_types   = length(sexp_types);
  int   *types     = INTEGER(sexp_types);
  if (n_threads < 1) n_threads = 1;
  dm_n_threads = n_threads;

  int n = em->n_samples, S = em->n_otus, nnz = em->nnz;
  dm_n   = n;
  dm_S   = S;
  dm_pos = em->pos_vec;
  dm_otu = em->otu_vec;
  dm_val = em->val_vec;
  dm_eta = asReal(sexp_eta);
  dm_nu  = asReal(sexp_nu);

  // Totals and multinomial constants.
  dm_J  = (double *) safe_malloc((size_t) (n + 1) * sizeof(double));
  dm_mc = (double *) safe_malloc((size_t) (n + 1) * sizeof(double));
  dm_pp = (double *) safe_malloc((size_t) (n + 1) * sizeof(double));
  for (int i = 0; i < n; i++) {
    double J = 0, lf = 0;
    for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++) { J += dm_val[p]; lf += lgammafn(dm_val[p] + 1); }
    dm_J[i]  = J;
    dm_mc[i] = lgammafn(J + 1) - lf;
    double pp = 0;
    for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++) { double q = dm_val[p] / J; pp += q * q; }
    dm_pp[i] = pp;
  }

  // By feature (stable: samples in order within each feature).
  dm_fpos = (int *)    safe_malloc((size_t) (S + 1) * sizeof(int));
  dm_fsam = (int *)    safe_malloc((size_t) (nnz + 1) * sizeof(int));
  dm_fval = (double *) safe_malloc((size_t) (nnz + 1) * sizeof(double));
  for (int j = 0; j <= S; j++) dm_fpos[j] = 0;
  for (int p = 0; p < nnz; p++) dm_fpos[dm_otu[p] + 1]++;
  for (int j = 0; j < S; j++) dm_fpos[j + 1] += dm_fpos[j];
  {
    int *fill = (int *) safe_malloc((size_t) (S + 1) * sizeof(int));
    for (int j = 0; j < S; j++) fill[j] = dm_fpos[j];
    for (int i = 0; i < n; i++)
      for (int p = dm_pos[i]; p < dm_pos[i + 1]; p++) {
        int q = fill[dm_otu[p]]++;
        dm_fsam[q] = i;
        dm_fval[q] = dm_val[p];
      }
    free_one(fill);
  }

  int Kmax = 1;
  for (int t = 0; t < n_types; t++) if (types[t] > Kmax) Kmax = types[t];

  size_t KS = (size_t) Kmax * S;
  dm_lam   = (double *) safe_malloc((KS + 1) * sizeof(double));
  dm_alpha = (double *) safe_malloc((KS + 1) * sizeof(double));
  dm_lga   = (double *) safe_malloc((KS + 1) * sizeof(double));
  dm_cent  = (double *) safe_malloc((KS + 1) * sizeof(double));
  dm_pi    = (double *) safe_malloc((size_t) (Kmax + 1) * sizeof(double));
  dm_cc    = (double *) safe_malloc((size_t) (Kmax + 1) * sizeof(double));
  dm_z     = (double *) safe_malloc(((size_t) n * Kmax + 1) * sizeof(double));
  dm_lse   = (double *) safe_malloc((size_t) (n + 1) * sizeof(double));
  dm_dist  = (double *) safe_malloc((size_t) (n + 1) * sizeof(double));
  dm_asg   = (int *)    safe_malloc((size_t) (n + 1) * sizeof(int));
  dm_ev_a  = (double *) safe_malloc((size_t) (S + 1) * sizeof(double));
  dm_fv    = (double *) safe_malloc((size_t) (3 * S + 1) * sizeof(double));
  dm_fg    = dm_fv + S;
  dm_fh    = dm_fg + S;
  dm_sv    = (double *) safe_malloc((size_t) (3 * n + 1) * sizeof(double));
  dm_sg    = dm_sv + n;
  dm_sh    = dm_sg + n;
  dm_nw    = (double *) safe_malloc((size_t) (6 * S + 1) * sizeof(double));

  // The best start's fit.
  double *b_lam   = (double *) safe_malloc((KS + 1) * sizeof(double));
  double *b_alpha = (double *) safe_malloc((KS + 1) * sizeof(double));
  double *b_pi    = (double *) safe_malloc((size_t) (Kmax + 1) * sizeof(double));
  double *b_z     = (double *) safe_malloc(((size_t) n * Kmax + 1) * sizeof(double));
  int    *ord     = (int *)    safe_malloc((size_t) (Kmax + 1) * sizeof(int));

  SEXP sexp_result = PROTECT(allocVector(VECSXP, n_types));
  int  interrupted = 0;

  for (int t = 0; t < n_types && !interrupted; t++) {

    int K = types[t];
    dm_K  = K;

    double best = R_PosInf;
    int    b_iters = 0, b_conv = 0, b_start = 0, any = 0;
    int    starts  = isNull(sexp_init) ? n_starts : 1;

    for (int s = 0; s < starts; s++) {

      if (!isNull(sexp_init)) {
        const double *z0 = REAL(sexp_init);
        for (size_t q = 0; q < (size_t) n * K; q++) dm_z[q] = z0[q];
      }
      else if (K == 1) {
        for (int i = 0; i < n; i++) dm_z[i] = 1;
      }
      else {
        rng_t rng;
        rng_seed(&rng, seed, (uint64_t) s);
        if (!kmeans_start(&rng, n_threads)) continue;
      }

      double nle;
      int    conv;
      int    iters = em_fit(max_iters, tol, n_threads, &nle, &conv);
      if (conv < 0) { interrupted = 1; break; }

      if (!any || nle < best) {
        any     = 1;
        best    = nle;
        b_iters = iters;
        b_conv  = conv;
        b_start = s + 1;
        memcpy(b_lam,   dm_lam,   (size_t) K * S * sizeof(double));
        memcpy(b_alpha, dm_alpha, (size_t) K * S * sizeof(double));
        memcpy(b_pi,    dm_pi,    (size_t) K * sizeof(double));
        memcpy(b_z,     dm_z,     (size_t) n * K * sizeof(double));
      }
    }
    if (interrupted) break;

    SEXP sexp_fit = PROTECT(allocVector(VECSXP, 9));
    SET_VECTOR_ELT(sexp_result, t, sexp_fit);
    UNPROTECT(1);
    SET_VECTOR_ELT(sexp_fit, 0, ScalarInteger(any));
    if (!any) continue;

    // Types in decreasing pi (stable).
    for (int k = 0; k < K; k++) ord[k] = k;
    for (int a = 1; a < K; a++) {
      int o = ord[a], b = a - 1;
      while (b >= 0 && b_pi[ord[b]] < b_pi[o]) { ord[b + 1] = ord[b]; b--; }
      ord[b + 1] = o;
    }
    for (int k = 0; k < K; k++) {
      int o = ord[k];
      memcpy(dm_lam   + (size_t) k * S, b_lam   + (size_t) o * S, (size_t) S * sizeof(double));
      memcpy(dm_alpha + (size_t) k * S, b_alpha + (size_t) o * S, (size_t) S * sizeof(double));
      memcpy(dm_z     + (size_t) k * n, b_z     + (size_t) o * n, (size_t) n * sizeof(double));
      dm_pi[k] = b_pi[o];
    }

    // log|H|.
    double ld = 0;
    for (int k = 0; k < K; k++) ld += block_log_det(k);
    ld += (K - 1) * log((double) n);
    for (int k = 0; k < K; k++) ld -= log(dm_pi[k]);

    SEXP sexp_alpha = PROTECT(allocMatrix(REALSXP, K, S));
    SEXP sexp_pi    = PROTECT(allocVector(REALSXP, K));
    SEXP sexp_z     = PROTECT(allocMatrix(REALSXP, n, K));
    double *ra = REAL(sexp_alpha);
    for (int k = 0; k < K; k++)
      for (int j = 0; j < S; j++)
        ra[k + (size_t) K * j] = dm_alpha[(size_t) k * S + j];
    memcpy(REAL(sexp_pi), dm_pi, (size_t) K * sizeof(double));
    memcpy(REAL(sexp_z),  dm_z,  (size_t) n * K * sizeof(double));

    SET_VECTOR_ELT(sexp_fit, 1, ScalarReal(best));
    SET_VECTOR_ELT(sexp_fit, 2, ScalarReal(ld));
    SET_VECTOR_ELT(sexp_fit, 3, ScalarInteger(b_iters));
    SET_VECTOR_ELT(sexp_fit, 4, ScalarLogical(b_conv));
    SET_VECTOR_ELT(sexp_fit, 5, ScalarInteger(b_start));
    SET_VECTOR_ELT(sexp_fit, 6, sexp_alpha);
    SET_VECTOR_ELT(sexp_fit, 7, sexp_pi);
    SET_VECTOR_ELT(sexp_fit, 8, sexp_z);
    UNPROTECT(3);
  }

  free_all();

  if (interrupted) {
    UNPROTECT(1);
    error("Interrupted.");
  }

  UNPROTECT(1);
  return sexp_result;
}
