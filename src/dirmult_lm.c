// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// dirmult_lm() is a C port of the linear-model path of aldex() from ALDEx3
// 1.4.0 (R/aldex.R, aldex-internals.R: aldex.sampler() and aldex.pvals(),
// fflm.R: fflm(), scale.R: clr.sm(), tss.sm() and sample.sm(),
// dirichlet.R: rLogDirichlet()), restructured for ecodive (see
// plans/dirmult_lm.md), under this license:
//
//   Copyright (c) 2026 ALDEx3 authors
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//
//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// The rest of this file is Copyright 2026 ecodive authors, under the MIT
// license above. The Gamma and normal draws are rng_dist.c's (on PCG, not
// R's rgamma() and rnorm()); least squares use a QR of the design
// (feature_lm.c) where fflm() solves the normal equations; random
// intercepts use feature_lmm.c (REML, Satterthwaite) where ALDEx3's
// method = "lme4" calls lmerTest.
//
// In ALDEx3, each Monte Carlo draw d gives every sample s a log2
// composition, log2 comp_sf = (log G_sf - log sum_g G_sg) / log 2 with
// G_sf ~ Gamma(x_sf + 0.5, 1), and a log2 scale from the scale model:
//
//   clr.sm():    -mean_f log2 comp_sf + x_s' L,  L ~ N(0, gamma^2 I_P)
//   tss.sm():    x_s' L
//   sample.sm(): N(log2 load_s, load_sd_s^2), independently per sample
//
// with one L (a coefficient per design column) per draw. The log2
// abundances Y = log2 comp + log2 scale are regressed on the design X,
// feature by feature: estimate b = (X'X)^-1 X'y, residuals e, and standard
// errors classical (s^2 (X'X)^-1), HC0 (sum_s A_ks^2 e_s^2, A = (X'X)^-1
// X') or HC3 (the same with e_s^2 / (1 - h_s)^2, h the leverages). Then,
// with t = b / se and N - P degrees of freedom, p_upper = pt(t, N - P,
// lower.tail = FALSE) and p_lower = 1 - p_upper (lmerTest's: p_lower =
// pt(t, df) with Satterthwaite's df, p_upper = 1 - p_lower). aldex.pvals()
// takes l = min(1, 2 p_lower) and u = min(1, 2 p_upper) per draw, and
// reports min(mean l, mean u) and min(mean adj(l), mean adj(u)), adj being
// p.adjust() across the features, per design column and draw.
// summary.aldex() reports the means of b and se over the draws.
//
// Here, for each draw d = 0 .. n_draws - 1:
//
//   Phase A (threads, one task per sample): the sample's row of log G,
//   from its own PCG stream (seed, d * (n + 1) + s), features in order (the
//   streams and order of dirmult_ttest(), so both see the same draws); then
//   its row of Y. The scale model's normals come from stream (seed,
//   d * (n + 1) + n), drawn on this thread before phase A: P of them in
//   design-column order (clr, tss), or n in sample order (load).
//
//   Phase B (threads, one task per block of features): each feature's
//   coefficients, standard errors and residual standard deviation, and its
//   mean clr value; or its random-intercept fit (lmm_fit()).
//
//   Then, on this thread (R's pt() can warn, which a worker thread must
//   not): t, the two tails, and running sums in draw order; and, with a
//   p_adjust method, R's p.adjust() of each design column's tails.
//
// Every task depends only on its own sample or features and the sums run in
// draw order on one thread, so results are identical for any number of
// threads. Memory: one n x D draw of doubles, plus P x D running sums.

#include "ecodive.h"
#include <Rmath.h>

#define LM_BLOCK  64    // features per phase B task (fixed effects)
#define LMM_BLOCK 4     // features per phase B task (random intercepts)

static int           dl_D, dl_P, dl_mode, dl_se_type, dl_mixed;
static uint64_t      dl_n;           // samples
static const int    *dl_pos, *dl_otu;
static const double *dl_val;
static double        dl_pc;
static uint64_t      dl_seed;
static int           dl_draw;
static const double *dl_given;       // n x D x n_draws given G, or NULL
static const double *dl_X;           // n x P design
static const double *dl_A;           // P x n: (X'X)^-1 X'
static const double *dl_h;           // n leverages
static const double *dl_xtxi;        // P: diagonal of (X'X)^-1
static double       *dl_Y;           // n x D, row by row
static double       *dl_sc;          // n: this draw's log2 scale (clr: the x'L part)
static double       *dl_cc;          // n: Y - clr, per sample
static double       *dl_lsc;         // n: Y - log2 comp, per sample (details)
static double       *dl_coef, *dl_se, *dl_df;   // P x D, this draw
static double       *dl_rsd, *dl_mean;          // D, this draw
static const lmm_design_t *dl_lmm;
static double       *dl_lmem;        // per thread: lmm workspace
static size_t        dl_lmem_sz;


// Phase A: each sample's row of Y.
static void *dl_draw_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D = dl_D, n = (int) dl_n;
  double ln2 = log(2.0);

  for (int s = thread_i; s < n; s += n_threads) {

    double *row = dl_Y + (size_t) s * D;

    if (dl_given) {
      const double *g = dl_given + (size_t) dl_draw * dl_n * D + s;
      for (int f = 0; f < D; f++) row[f] = log(g[(size_t) f * dl_n]);
    }
    else {
      // As dirmult_ttest(): the sample's values densely, then the draws in
      // feature order, so that sparse and dense tables agree.
      for (int f = 0; f < D; f++) row[f] = 0;
      for (int e = dl_pos[s]; e < dl_pos[s + 1]; e++) row[dl_otu[e]] = dl_val[e];

      rng_t rng;
      rng_seed(&rng, dl_seed, (uint64_t) dl_draw * (dl_n + 1) + (uint64_t) s);
      for (int f = 0; f < D; f++) row[f] = rng_log_gamma(&rng, row[f] + dl_pc);
    }

    if (dl_mode == 0) {
      // clr.sm(): log2 comp - mean(log2 comp) + x'L; the closure cancels.
      double sum = 0;
      for (int f = 0; f < D; f++) sum += row[f];
      double mean = sum / D;
      double add  = dl_sc[s];
      if (dl_lsc) {
        // Only for the details: the log2 scale, Y - log2 comp.
        double mx = R_NegInf;
        for (int f = 0; f < D; f++) if (row[f] > mx) mx = row[f];
        double se = 0;
        for (int f = 0; f < D; f++) se += exp(row[f] - mx);
        dl_lsc[s] = (mx + log(se) - mean) / ln2 + add;
      }
      for (int f = 0; f < D; f++) row[f] = (row[f] - mean) / ln2 + add;
      dl_cc[s] = add;
    }
    else {
      // log2 composition (rLogDirichlet(): log-sum-exp), plus the scale.
      double mx = R_NegInf;
      for (int f = 0; f < D; f++) if (row[f] > mx) mx = row[f];
      double sum = 0;
      for (int f = 0; f < D; f++) sum += exp(row[f] - mx);
      double lse = mx + log(sum);
      double add = dl_sc[s], tot = 0;
      for (int f = 0; f < D; f++) {
        double lc = (row[f] - lse) / ln2;
        tot   += lc;
        row[f] = lc + add;
      }
      dl_cc[s] = add - tot / D;   // Y - clr
      if (dl_lsc) dl_lsc[s] = add;
    }
  }
  return NULL;
}


// Phase B, fixed effects: coefficients A y, residuals, standard errors,
// for a block of features at a time, summing over samples in order.
static void *dl_lm_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D = dl_D, P = dl_P, n = (int) dl_n;

  double c[LM_BLOCK * 16], v[LM_BLOCK * 16];    // P <= 16 here (else see below)
  double *cb = c, *vb = v;
  double *heap = NULL;
  if (P > 16) {
    heap = (double *) malloc(sizeof(double) * 2 * (size_t) LM_BLOCK * P);
    if (!heap) return NULL;   // # nocov
    cb = heap; vb = heap + (size_t) LM_BLOCK * P;
  }
  double rss[LM_BLOCK], mn[LM_BLOCK];

  int n_blocks = (D + LM_BLOCK - 1) / LM_BLOCK;
  for (int b = thread_i; b < n_blocks; b += n_threads) {

    int f0 = b * LM_BLOCK;
    int nf = D - f0 < LM_BLOCK ? D - f0 : LM_BLOCK;

    for (int k = 0; k < P * nf; k++) cb[k] = vb[k] = 0;
    for (int j = 0; j < nf; j++) rss[j] = mn[j] = 0;

    // b = A y
    for (int s = 0; s < n; s++) {
      const double *row = dl_Y + (size_t) s * D + f0;
      for (int k = 0; k < P; k++) {
        double a = dl_A[k + (size_t) P * s];
        double *ck = cb + (size_t) k * nf;
        for (int j = 0; j < nf; j++) ck[j] += a * row[j];
      }
    }

    // Residuals, their squares (HC: weighted), and the mean clr.
    for (int s = 0; s < n; s++) {
      const double *row = dl_Y + (size_t) s * D + f0;
      double hw = 1;
      if (dl_se_type == 2) {
        double d = 1 - dl_h[s];
        hw = d > 1e-10 ? 1 / (d * d) : R_NaN;   // leverage 1: undefined
      }
      double cc = dl_cc[s];
      for (int j = 0; j < nf; j++) {
        double fit = 0;
        for (int k = 0; k < P; k++) fit += dl_X[s + (size_t) n * k] * cb[(size_t) k * nf + j];
        double e = row[j] - fit;
        double e2 = e * e;
        rss[j] += e2;
        mn[j]  += row[j] - cc;
        if (dl_se_type > 0) {
          double w = e2 * hw;
          for (int k = 0; k < P; k++) {
            double a = dl_A[k + (size_t) P * s];
            vb[(size_t) k * nf + j] += a * a * w;
          }
        }
      }
    }

    for (int j = 0; j < nf; j++) {
      int f = f0 + j;
      double s2 = rss[j] / (n - P);
      dl_rsd[f]  = sqrt(s2);
      dl_mean[f] = mn[j] / n;
      for (int k = 0; k < P; k++) {
        size_t o = (size_t) k + (size_t) P * f;
        dl_coef[o] = cb[(size_t) k * nf + j];
        dl_se[o]   = dl_se_type == 0 ? sqrt(dl_xtxi[k] * s2) : sqrt(vb[(size_t) k * nf + j]);
      }
    }
  }

  if (heap) free(heap);
  return NULL;
}


// Phase B, random intercepts: each feature's REML fit.
static void *dl_lmm_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D = dl_D, P = dl_P, n = (int) dl_n;

  double *mem  = dl_lmem + (size_t) thread_i * dl_lmem_sz;
  double *y    = mem + lmm_doubles(n, P, dl_lmm->J);
  double *coef = y + n;
  double *df   = coef + P;
  double *vc   = df + P;

  int n_blocks = (D + LMM_BLOCK - 1) / LMM_BLOCK;
  for (int b = thread_i; b < n_blocks; b += n_threads) {
    int f1 = (b + 1) * LMM_BLOCK < D ? (b + 1) * LMM_BLOCK : D;
    for (int f = b * LMM_BLOCK; f < f1; f++) {

      double mn = 0;
      for (int s = 0; s < n; s++) {
        y[s] = dl_Y[(size_t) s * D + f];
        mn  += y[s] - dl_cc[s];
      }
      dl_mean[f] = mn / n;

      double sigma, dev;
      lmm_fit(dl_lmm, y, -1, mem, coef, vc, df, &sigma, &dev);
      dl_rsd[f] = sigma;
      for (int k = 0; k < P; k++) {
        size_t o = (size_t) k + (size_t) P * f;
        dl_coef[o] = coef[k];
        dl_se[o]   = sqrt(vc[k + (size_t) P * k]);
        dl_df[o]   = df[k];
      }
    }
  }
  return NULL;
}


static inline double nan_min (double a, double b) {
  return (ISNAN(a) || ISNAN(b)) ? R_NaN : (b < a ? b : a);
}
static inline double cap1 (double x) { return x > 1 ? 1 : x; }   // keeps NaN



//======================================================
// R interface for dirmult_lm().
//
// sexp_X:       n x P design (full rank).
// sexp_se_type: 0 classical, 1 HC0, 2 HC3.
// sexp_mode:    0 clr, 1 tss, 2 load, 3 given scale.
// sexp_sd:      gamma (clr, tss), in log2 units.
// sexp_mu:      log2 load per sample (load), or NULL.
// sexp_lsd:     load_sd per sample (load), or NULL.
// sexp_groups:  1-based random-intercept groups, or NULL.
// sexp_padj:    function(p) of the adjusted p-values, or
//               NULL for none.
// sexp_given:   NULL, or n x D x n_draws Gamma draws.
// sexp_gscale:  NULL, or n x n_draws log2 scales (mode 3).
// sexp_details: TRUE also returns each draw's estimate,
//               se, df, p_lower, p_upper (P x D x n_draws),
//               scale (n x n_draws, log2: Y - log2 comp)
//               and shift (P or n x n_draws: L, or the
//               load's normal deviates times load_sd).
//
// Returns list(estimate, se, df (P x D), p_value,
// p_adjusted (P x D), resid_sd, mean (D)[, details]).
//======================================================
SEXP C_dirmult_lm(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_pseudocount, SEXP sexp_n_draws,
    SEXP sexp_X, SEXP sexp_se_type, SEXP sexp_mode, SEXP sexp_sd, SEXP sexp_mu,
    SEXP sexp_lsd, SEXP sexp_groups, SEXP sexp_padj, SEXP sexp_seed,
    SEXP sexp_given, SEXP sexp_gscale, SEXP sexp_details, SEXP sexp_n_threads ) {

  double pc      = asReal(sexp_pseudocount);
  int    n_draws = asInteger(sexp_n_draws);
  int    se_type = asInteger(sexp_se_type);
  int    mode    = asInteger(sexp_mode);
  double gam     = asReal(sexp_sd);
  int    details = asLogical(sexp_details);
  int    n_thr   = asInteger(sexp_n_threads);
  int    mixed   = !isNull(sexp_groups);
  if (n_thr < 1) n_thr = 1;

  // The table, copied out of the ecomatrix into R's memory.
  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int n   = em->n_samples;
  int D   = em->n_otus;
  int nnz = em->nnz;

  int    *pos = (int *)    R_alloc((size_t) n + 1,   sizeof(int));
  int    *otu = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  double *val = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  memcpy(pos, em->pos_vec, ((size_t) n + 1) * sizeof(int));
  if (nnz) {
    memcpy(otu, em->otu_vec, (size_t) nnz * sizeof(int));
    memcpy(val, em->val_vec, (size_t) nnz * sizeof(double));
  }
  free_all();

  int P = ncols(sexp_X);
  if (nrows(sexp_X) != n) error("The design must have one row per sample.");  // # nocov
  const double *X = REAL(sexp_X);


  //--------------------------------------------------------
  // The design: A = (X'X)^-1 X' = R^-1 Q1', the leverages,
  // and the diagonal of (X'X)^-1, from one QR.
  //--------------------------------------------------------
  double *A    = (double *) R_alloc((size_t) P * n + 1, sizeof(double));
  double *xtxi = (double *) R_alloc((size_t) P + 1, sizeof(double));
  double *hv   = (double *) R_alloc((size_t) n + 1, sizeof(double));
  {
    double *fmem  = (double *) R_alloc(flm_doubles(n, P, 1), sizeof(double));
    int    *fimem = (int *)    R_alloc((size_t) P + 1, sizeof(int));
    flm_t q;
    flm_init(&q, fmem, fimem, n, P, 1);
    flm_decompose(&q, X, n, NULL, n, P, 1e-7);
    if (q.rank < P) error("The design is rank deficient.");   // # nocov (checked in R)
    for (int k = 0; k < P; k++) {
      int pk = q.pivot[k];
      long double d = 0;
      for (int m = k; m < P; m++) { double r = q.rinv[k + (size_t) P * m]; d += (long double) r * r; }
      xtxi[pk] = (double) d;
      for (int s = 0; s < n; s++) {
        long double a = 0;
        for (int m = k; m < P; m++)
          a += (long double) q.rinv[k + (size_t) P * m] * q.q1[s + (size_t) n * m];
        A[pk + (size_t) P * s] = (double) a;
      }
    }
    memcpy(hv, q.h, (size_t) n * sizeof(double));
  }

  lmm_design_t lmm;
  int nt_b_max = n_thr;
  if (mixed) {
    int *grp = (int *) R_alloc((size_t) n, sizeof(int));
    int J = 0;
    for (int s = 0; s < n; s++) { grp[s] = INTEGER(sexp_groups)[s] - 1; if (grp[s] + 1 > J) J = grp[s] + 1; }
    lmm_design_init(&lmm, X, grp, n, P, J);
    dl_lmm     = &lmm;
    dl_lmem_sz = lmm_doubles(n, P, J) + (size_t) n + 2 * (size_t) P + (size_t) P * P + 8;
    dl_lmem    = (double *) R_alloc((size_t) nt_b_max * dl_lmem_sz, sizeof(double));
  }

  double *Y    = (double *) R_alloc((size_t) n * D + 1, sizeof(double));
  double *sc   = (double *) R_alloc((size_t) n + 1, sizeof(double));
  double *cc   = (double *) R_alloc((size_t) n + 1, sizeof(double));
  double *coef = (double *) R_alloc((size_t) P * D + 1, sizeof(double));
  double *se   = (double *) R_alloc((size_t) P * D + 1, sizeof(double));
  double *dfv  = (double *) R_alloc((size_t) P * D + 1, sizeof(double));
  double *rsd  = (double *) R_alloc((size_t) D + 1, sizeof(double));
  double *mn   = (double *) R_alloc((size_t) D + 1, sizeof(double));
  double *lam  = (double *) R_alloc((size_t) (P > n ? P : n) + 1, sizeof(double));
  double *pl   = (double *) R_alloc((size_t) P * D + 1, sizeof(double));
  double *pu   = (double *) R_alloc((size_t) P * D + 1, sizeof(double));

  dl_D       = D;
  dl_P       = P;
  dl_n       = (uint64_t) n;
  dl_mode    = mode;
  dl_se_type = se_type;
  dl_mixed   = mixed;
  dl_pos     = pos;
  dl_otu     = otu;
  dl_val     = val;
  dl_pc      = pc;
  dl_seed    = (uint64_t) asInteger(sexp_seed);
  dl_given   = isNull(sexp_given) ? NULL : REAL(sexp_given);
  dl_X       = X;
  dl_A       = A;
  dl_h       = hv;
  dl_xtxi    = xtxi;
  dl_Y       = Y;
  dl_sc      = sc;
  dl_cc      = cc;
  dl_coef    = coef;
  dl_se      = se;
  dl_df      = dfv;
  dl_rsd     = rsd;
  dl_mean    = mn;
  dl_lsc     = details ? (double *) R_alloc((size_t) n + 1, sizeof(double)) : NULL;

  const double *mu  = isNull(sexp_mu)     ? NULL : REAL(sexp_mu);
  const double *lsd = isNull(sexp_lsd)    ? NULL : REAL(sexp_lsd);
  const double *gsc = isNull(sexp_gscale) ? NULL : REAL(sexp_gscale);

  int nt_a = threads_for(n_thr, (double) n * D * 40, 1e6);
  int nt_b = mixed
    ? threads_for(n_thr, (double) D * n * 2000, 1e6)
    : threads_for(n_thr, (double) n * D * (3 * P + 4), 1e6);
  int nb   = mixed ? (D + LMM_BLOCK - 1) / LMM_BLOCK : (D + LM_BLOCK - 1) / LM_BLOCK;
  if (nt_b > nb) nt_b = nb < 1 ? 1 : nb;


  //--------------------------------------------------------
  // Results and running sums.
  //--------------------------------------------------------
  const char *nms[] = {
    "estimate", "se", "df", "p_value", "p_adjusted", "resid_sd", "mean", "details", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_est  = allocMatrix(REALSXP, P, D); SET_VECTOR_ELT(result, 0, s_est);
  SEXP s_se   = allocMatrix(REALSXP, P, D); SET_VECTOR_ELT(result, 1, s_se);
  SEXP s_df   = allocMatrix(REALSXP, P, D); SET_VECTOR_ELT(result, 2, s_df);
  SEXP s_p    = allocMatrix(REALSXP, P, D); SET_VECTOR_ELT(result, 3, s_p);
  SEXP s_padj = allocMatrix(REALSXP, P, D); SET_VECTOR_ELT(result, 4, s_padj);
  SEXP s_rsd  = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 5, s_rsd);
  SEXP s_mean = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 6, s_mean);
  double *est_sum = REAL(s_est), *se_sum = REAL(s_se), *df_sum = REAL(s_df);
  double *rsd_sum = REAL(s_rsd), *mean_sum = REAL(s_mean);
  size_t PD = (size_t) P * D;

  double *l_sum  = (double *) R_alloc(PD + 1, sizeof(double));
  double *u_sum  = (double *) R_alloc(PD + 1, sizeof(double));
  double *la_sum = (double *) R_alloc(PD + 1, sizeof(double));
  double *ua_sum = (double *) R_alloc(PD + 1, sizeof(double));
  for (size_t i = 0; i < PD; i++)
    est_sum[i] = se_sum[i] = df_sum[i] = l_sum[i] = u_sum[i] = la_sum[i] = ua_sum[i] = 0;
  for (int f = 0; f < D; f++) rsd_sum[f] = mean_sum[f] = 0;

  double *d_est = NULL, *d_se = NULL, *d_df = NULL, *d_pl = NULL, *d_pu = NULL,
         *d_sc = NULL, *d_sh = NULL;
  int n_shift = mode == 2 ? n : P;
  if (details) {
    const char *dn[] = { "estimate", "se", "df", "p_lower", "p_upper", "scale", "shift", "" };
    SEXP det = mkNamed(VECSXP, dn);
    SET_VECTOR_ELT(result, 7, det);
    SEXP dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dim)[0] = P; INTEGER(dim)[1] = D; INTEGER(dim)[2] = n_draws;
    double **ptr[] = { &d_est, &d_se, &d_df, &d_pl, &d_pu };
    for (int k = 0; k < 5; k++) {
      SEXP a = allocVector(REALSXP, (R_xlen_t) PD * n_draws);
      SET_VECTOR_ELT(det, k, a);
      setAttrib(a, R_DimSymbol, dim);
      *ptr[k] = REAL(a);
    }
    SEXP a = allocMatrix(REALSXP, n, n_draws);       SET_VECTOR_ELT(det, 5, a); d_sc = REAL(a);
    SEXP b = allocMatrix(REALSXP, n_shift, n_draws); SET_VECTOR_ELT(det, 6, b); d_sh = REAL(b);
    UNPROTECT(1);
  }

  // p.adjust() of one design column's tails, through R.
  SEXP padj_in = R_NilValue, padj_call = R_NilValue;
  int  n_prot  = 1;
  if (!isNull(sexp_padj)) {
    padj_in   = PROTECT(allocVector(REALSXP, D));
    padj_call = PROTECT(lang2(sexp_padj, padj_in));
    n_prot   += 2;
  }

  double dof_fixed = n - P;


  //--------------------------------------------------------
  // The draws.
  //--------------------------------------------------------
  for (int d = 0; d < n_draws; d++) {

    dl_draw = d;

    // The scale model's draws, on this thread.
    for (int s = 0; s < n; s++) sc[s] = 0;
    if (mode == 3) {
      for (int s = 0; s < n; s++) sc[s] = gsc[s + (size_t) n * d];
    }
    else if (mode == 2) {
      rng_t rng;
      rng_seed(&rng, dl_seed, (uint64_t) d * (dl_n + 1) + dl_n);
      for (int s = 0; s < n; s++) {
        lam[s] = lsd[s] * rng_norm(&rng);
        sc[s]  = mu[s] + lam[s];
      }
    }
    else {
      for (int k = 0; k < P; k++) lam[k] = 0;
      if (gam != 0) {
        rng_t rng;
        rng_seed(&rng, dl_seed, (uint64_t) d * (dl_n + 1) + dl_n);
        for (int k = 0; k < P; k++) lam[k] = gam * rng_norm(&rng);
        for (int s = 0; s < n; s++) {
          double v = 0;
          for (int k = 0; k < P; k++) v += X[s + (size_t) n * k] * lam[k];
          sc[s] = v;
        }
      }
    }

    run_parallel(dl_draw_worker, nt_a, n < 100 ? 100 : n);
    run_parallel(mixed ? dl_lmm_worker : dl_lm_worker, nt_b, nb < 100 ? 100 : nb);

    if (details) {
      // The log2 scale: Y - log2 comp (for clr: x'L - mean log2 comp).
      for (int s = 0; s < n; s++) d_sc[s + (size_t) n * d] = dl_lsc[s];
      for (int k = 0; k < n_shift; k++) d_sh[k + (size_t) n_shift * d] = mode == 3 ? R_NaN : lam[k];
    }

    // t, the tails, and the running sums.
    for (int f = 0; f < D; f++) {
      rsd_sum[f]  += rsd[f];
      mean_sum[f] += mn[f];
      for (int k = 0; k < P; k++) {
        size_t o  = (size_t) k + (size_t) P * f;
        double t  = coef[o] / se[o];
        double df = mixed ? dfv[o] : dof_fixed;
        double lo = R_NaN, up = R_NaN;
        if (!ISNAN(t) && df > 0) {
          if (mixed) { lo = pt(t, df, 1, 0); up = 1 - lo; }
          else       { up = pt(t, df, 0, 0); lo = 1 - up; }
        }
        double l = cap1(2 * lo), u = cap1(2 * up);
        pl[o] = l; pu[o] = u;
        est_sum[o] += coef[o];
        se_sum[o]  += se[o];
        df_sum[o]  += df;
        l_sum[o]   += l;
        u_sum[o]   += u;
        if (details) {
          size_t q = o + PD * d;
          d_est[q] = coef[o]; d_se[q] = se[o]; d_df[q] = df;
          d_pl[q] = lo; d_pu[q] = up;
        }
      }
    }

    if (isNull(sexp_padj)) {
      for (size_t i = 0; i < PD; i++) { la_sum[i] += pl[i]; ua_sum[i] += pu[i]; }
    }
    else {
      double *in = REAL(padj_in);
      for (int side = 0; side < 2; side++) {
        double *src = side ? pu : pl, *dst = side ? ua_sum : la_sum;
        for (int k = 0; k < P; k++) {
          for (int f = 0; f < D; f++) in[f] = src[(size_t) k + (size_t) P * f];
          SEXP r = eval(padj_call, R_GlobalEnv);
          if (TYPEOF(r) != REALSXP || XLENGTH(r) != D) error("p.adjust() failed.");  // # nocov
          const double *rv = REAL(r);
          for (int f = 0; f < D; f++) dst[(size_t) k + (size_t) P * f] += rv[f];
        }
      }
    }

    if (d + 1 < n_draws && user_interrupted()) {
      UNPROTECT(n_prot);
      error("Interrupted.");
    }
  }


  //--------------------------------------------------------
  // Means over the draws.
  //--------------------------------------------------------
  double *pv = REAL(s_p), *pa = REAL(s_padj);
  for (size_t i = 0; i < PD; i++) {
    est_sum[i] /= n_draws;
    se_sum[i]  /= n_draws;
    df_sum[i]  /= n_draws;
    pv[i] = nan_min(l_sum[i] / n_draws, u_sum[i] / n_draws);
    pa[i] = nan_min(la_sum[i] / n_draws, ua_sum[i] / n_draws);
  }
  for (int f = 0; f < D; f++) {
    rsd_sum[f]  /= n_draws;
    mean_sum[f] /= n_draws;
  }

  UNPROTECT(n_prot);
  return result;
}
