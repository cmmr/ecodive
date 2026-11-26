// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// nmit() and volatility()'s bootstrap intervals.
//
// nmit() follows q2-longitudinal 2026.10.0's _nmit(), _temporal_corr() and
// _temporal_distance() (q2_longitudinal/_utilities.py); the correlations
// are computed as the pandas and SciPy functions they call compute them
// (pandas 2.3.3's libalgos.nancorr() and nancorr_spearman(), and
// DataFrame.corr()'s Kendall loop over scipy.stats.kendalltau()), and
// volatility()'s intervals follow Vega's bootstrapCI() (vega-statistics),
// the `ci0` and `ci1` aggregates of q2-longitudinal's volatility chart.
// Those sources are under these licenses:
//
//   q2-longitudinal: Copyright (c) 2017-2026, QIIME 2 development team.
//   pandas: Copyright (c) 2008-2011, AQR Capital Management, LLC, Lambda
//     Foundry, Inc. and PyData Development Team; Copyright (c) 2011-2025,
//     Open source contributors.
//   SciPy: Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
//   Vega: Copyright (c) 2015-2023, University of Washington Interactive
//     Data Lab.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//
//   3. Neither the name of the copyright holder nor the names of its
//      contributors may be used to endorse or promote products derived
//      from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The rest of this file is Copyright 2026 ecodive authors, under the MIT
// license above.
//
// nmit() runs in two phases, each split into tasks that don't depend on the
// number of threads, so results are identical for any `cpus`:
//
//   A. One task per subject: the D x D correlation matrix of the features
//      over the subject's samples (in their input order), stored as a
//      packed lower triangle with its diagonal, D (D + 1) / 2 values.
//      NaN (a feature constant within the subject) is 0, as q2's
//      fillna(0).
//   B. One task per column of the subjects' distance triangle: the matrix
//      norm of each difference of two subjects' correlation matrices.
//      Frobenius: sqrt(2 sum_{a>b} d_ab^2 + sum_a d_aa^2), summed in the
//      packed order. Nuclear: the sum of the absolute eigenvalues of the
//      symmetric difference (lz_sym_eigen()), O(D^3) per pair of subjects.
//
// The correlations, as pandas and SciPy compute them (a > b; x is column
// a and y column b for pandas' kernels, the other way round for Kendall,
// as DataFrame.corr() calls kendalltau(column b, column a)):
//
//   Pearson   Welford's updates over the samples, then cov / sqrt(ssx ssy):
//             1 on the diagonal (exactly), NaN for a constant feature.
//   Spearman  Average ranks (1-based; ties share their mean rank), centred
//             on (n + 1) / 2, then the same ratio of sums.
//   Kendall   Tau-b from counts of the sample pairs: with s_x and s_y the
//             signs of each pair's differences, con - dis = sum s_x s_y,
//             and tau = (con - dis) / sqrt(tot - xtie) / sqrt(tot - ytie),
//             clamped to [-1, 1]; NaN when either feature is all ties or
//             there is one sample. The diagonal is 1, whatever the feature.

#include "ecodive.h"

#define NMIT_PEARSON  1
#define NMIT_SPEARMAN 2
#define NMIT_KENDALL  3

#define NMIT_FROBENIUS 1
#define NMIT_NUCLEAR   2


static ecomatrix_t *nm_em;
static int          nm_D, nm_S, nm_nmax, nm_method, nm_mnorm, nm_zero_clr;
static const int   *nm_spos;   // S + 1 offsets into nm_ssam
static const int   *nm_ssam;   // samples, grouped by subject, in input order
static double      *nm_corr;   // S x D (D + 1) / 2 packed correlations
static size_t       nm_plen;   // D (D + 1) / 2
static double      *nm_work;   // per-share workspace
static size_t       nm_wlen;   // its length per share, in doubles
static int         *nm_iwork;  // per-share integer workspace
static size_t       nm_iwlen;
static double      *nm_out;    // S (S - 1) / 2 distances
static int         *nm_failed; // per share: an eigendecomposition failed


// Column f of a subject's dense n x D block X (column-major).
#define COL(X, f, n) ((X) + (size_t) (f) * (n))


// Average ranks (1-based) of v[0 .. n - 1] into r, as pandas' rank_1d().
static void nm_rank (const double *v, double *r, int *idx, int n) {

  // Insertion sort of the indices by value (n is small): stable.
  for (int i = 0; i < n; i++) {
    int j = i;
    while (j > 0 && v[idx[j - 1]] > v[i]) { idx[j] = idx[j - 1]; j--; }
    idx[j] = i;
  }

  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && v[idx[j + 1]] == v[idx[i]]) j++;
    double avg = ((double) (i + 1) + (double) (j + 1)) / 2;
    for (int k = i; k <= j; k++) r[idx[k]] = avg;
    i = j + 1;
  }
}


static void *nm_corr_worker (void *arg) {

  int share = ((worker_t *) arg)->i;
  int n_sh  = ((worker_t *) arg)->n;
  int D     = nm_D;

  double *X  = nm_work + (size_t) share * nm_wlen;                // nmax x D
  double *Rk = X + (size_t) nm_nmax * D;                           // nmax x D ranks
  int    *ix = nm_iwork + (size_t) share * nm_iwlen;               // nmax
  int    *xt = ix + nm_nmax;                                       // D ties
  signed char *sg = (signed char *) (xt + D);                      // P x D signs

  for (int s = share; s < nm_S; s += n_sh) {

    int    n  = nm_spos[s + 1] - nm_spos[s];
    double *C = nm_corr + (size_t) s * nm_plen;

    // The subject's dense block.
    for (int i = 0; i < n; i++) {
      int    sam  = nm_ssam[nm_spos[s] + i];
      double zero = nm_zero_clr ? nm_em->clr_vec[sam] : 0;
      for (int f = 0; f < D; f++) COL(X, f, n)[i] = zero;
      for (int e = nm_em->pos_vec[sam]; e < nm_em->pos_vec[sam + 1]; e++)
        COL(X, nm_em->otu_vec[e], n)[i] = nm_em->val_vec[e];
    }

    size_t idx = 0;

    if (nm_method == NMIT_PEARSON) {

      for (int a = 0; a < D; a++) {
        const double *x = COL(X, a, n);
        for (int b = 0; b <= a; b++, idx++) {
          const double *y = COL(X, b, n);
          double meanx = 0, meany = 0, ssqdmx = 0, ssqdmy = 0, covxy = 0;
          for (int i = 0; i < n; i++) {
            double vx = x[i], vy = y[i];
            double nobs = i + 1;
            double dx = vx - meanx, dy = vy - meany;
            meanx += 1. / nobs * dx;
            meany += 1. / nobs * dy;
            ssqdmx += (vx - meanx) * dx;
            ssqdmy += (vy - meany) * dy;
            covxy  += (vx - meanx) * dy;
          }
          double div = sqrt(ssqdmx * ssqdmy);
          double r   = div != 0 ? covxy / div : NAN;
          C[idx] = ISNAN(r) ? 0 : r;
        }
      }
    }

    else if (nm_method == NMIT_SPEARMAN) {

      for (int f = 0; f < D; f++) nm_rank(COL(X, f, n), COL(Rk, f, n), ix, n);
      double mean = (n + 1) / 2.;

      for (int a = 0; a < D; a++) {
        const double *x = COL(Rk, a, n);
        for (int b = 0; b <= a; b++, idx++) {
          const double *y = COL(Rk, b, n);
          double sxy = 0, sxx = 0, syy = 0;
          for (int i = 0; i < n; i++) {
            double vx = x[i] - mean, vy = y[i] - mean;
            sxy += vx * vy;
            sxx += vx * vx;
            syy += vy * vy;
          }
          double div = sqrt(sxx * syy);
          double r   = div != 0 ? sxy / div : NAN;
          C[idx] = ISNAN(r) ? 0 : r;
        }
      }
    }

    else {

      // Each feature's sign of x_j - x_i for the sample pairs i < j.
      int64_t tot = (int64_t) n * (n - 1) / 2;
      for (int f = 0; f < D; f++) {
        const double *x = COL(X, f, n);
        signed char  *g = sg + (size_t) f * tot;
        int ties = 0;
        for (int i = 0, p = 0; i < n; i++)
          for (int j = i + 1; j < n; j++, p++) {
            g[p] = (signed char) ((x[j] > x[i]) - (x[j] < x[i]));
            ties += g[p] == 0;
          }
        xt[f] = ties;
      }

      for (int a = 0; a < D; a++) {
        const signed char *ga = sg + (size_t) a * tot;
        for (int b = 0; b <= a; b++, idx++) {
          if (a == b) { C[idx] = 1; continue; }
          if (n < 2 || xt[a] == tot || xt[b] == tot) { C[idx] = 0; continue; }
          const signed char *gb = sg + (size_t) b * tot;
          int64_t cmd = 0;
          for (int64_t p = 0; p < tot; p++) cmd += ga[p] * gb[p];
          double tau = (double) cmd / sqrt((double) (tot - xt[b])) / sqrt((double) (tot - xt[a]));
          if (tau > 1)  tau = 1;
          if (tau < -1) tau = -1;
          C[idx] = tau;
        }
      }
    }
  }

  return NULL;
}


static void *nm_dist_worker (void *arg) {

  int share = ((worker_t *) arg)->i;
  int n_sh  = ((worker_t *) arg)->n;
  int D     = nm_D;
  int S     = nm_S;

  double *A    = nm_work + (size_t) share * nm_wlen;   // D x D
  double *vecs = A + (size_t) D * D;                   // D x D
  double *vals = vecs + (size_t) D * D;                // D
  double *tmp  = vals + D;                             // 2 D

  for (int r = share; r < S - 1; r += n_sh) {

    size_t        out = (size_t) r * S - (size_t) r * (r + 1) / 2;
    const double *Ci  = nm_corr + (size_t) r * nm_plen;

    for (int c = r + 1; c < S; c++, out++) {

      const double *Cj = nm_corr + (size_t) c * nm_plen;

      if (nm_mnorm == NMIT_FROBENIUS) {
        double off = 0, diag = 0;
        size_t idx = 0;
        for (int a = 0; a < D; a++) {
          for (int b = 0; b < a; b++, idx++) {
            double d = Ci[idx] - Cj[idx];
            off += d * d;
          }
          double d = Ci[idx] - Cj[idx];
          diag += d * d;
          idx++;
        }
        nm_out[out] = sqrt(2 * off + diag);
      }

      else {
        size_t idx = 0;
        for (int a = 0; a < D; a++)
          for (int b = 0; b <= a; b++, idx++)
            A[(size_t) b * D + a] = Ci[idx] - Cj[idx];   // lower triangle
        if (lz_sym_eigen(A, D, vals, vecs, tmp)) nm_failed[share] = 1;
        double sum = 0;
        for (int k = 0; k < D; k++) sum += fabs(vals[k]);
        nm_out[out] = sum;
      }
    }
  }

  return NULL;
}


//======================================================
// R interface for nmit().
//
// sexp_subj:    each sample's subject, 0 .. S - 1 (S =
//               sexp_n_subj), or -1 to leave it out.
// sexp_method:  1 Pearson, 2 Spearman, 3 Kendall.
// sexp_mnorm:   1 Frobenius, 2 nuclear.
// Returns the packed distances, S (S - 1) / 2 values in
// the order of a `dist`.
//======================================================
SEXP C_nmit (
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_pseudocount,
    SEXP sexp_subj, SEXP sexp_n_subj, SEXP sexp_method, SEXP sexp_mnorm,
    SEXP sexp_n_threads ) {

  int norm  = asInteger(sexp_norm);
  int S     = asInteger(sexp_n_subj);
  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  init_n_ptrs(20);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  if (norm) normalize(em, norm, n_thr, asReal(sexp_pseudocount));

  int n = em->n_samples;
  int D = em->n_otus;
  const int *subj = INTEGER(sexp_subj);

  // Samples grouped by subject, each in input order.
  int *spos = (int *) R_alloc((size_t) S + 1, sizeof(int));
  int *ssam = (int *) R_alloc((size_t) n + 1, sizeof(int));
  int *fill = (int *) R_alloc((size_t) S + 1, sizeof(int));
  for (int s = 0; s <= S; s++) spos[s] = 0;
  for (int i = 0; i < n; i++) if (subj[i] >= 0) spos[subj[i] + 1]++;
  for (int s = 0; s < S; s++) spos[s + 1] += spos[s];
  for (int s = 0; s < S; s++) fill[s] = spos[s];
  for (int i = 0; i < n; i++) if (subj[i] >= 0) ssam[fill[subj[i]]++] = i;

  int nmax = 1;
  for (int s = 0; s < S; s++) if (spos[s + 1] - spos[s] > nmax) nmax = spos[s + 1] - spos[s];

  nm_em       = em;
  nm_D        = D;
  nm_S        = S;
  nm_nmax     = nmax;
  nm_method   = asInteger(sexp_method);
  nm_mnorm    = asInteger(sexp_mnorm);
  nm_zero_clr = (norm == 2 || norm == 7);   // clr, log: zeros have values
  nm_spos     = spos;
  nm_ssam     = ssam;
  nm_plen     = (size_t) D * (D + 1) / 2;

  nm_corr = (double *) safe_malloc((size_t) S * nm_plen * sizeof(double));

  // Phase A: each subject's correlations.
  size_t P   = (size_t) nmax * (nmax - 1) / 2;
  double wA  = (double) S * D * D * (nmax < 1 ? 1 : nmax);
  int    thA = threads_for(n_thr, wA, 1e6);
  if (thA > S) thA = S < 1 ? 1 : S;
  nm_wlen  = (size_t) 2 * nmax * D;
  nm_iwlen = (size_t) nmax + D + (P * D + sizeof(int) - 1) / sizeof(int) + 1;
  nm_work  = (double *) safe_malloc((size_t) thA * nm_wlen * sizeof(double));
  nm_iwork = (int *)    safe_malloc((size_t) thA * nm_iwlen * sizeof(int));
  run_parallel(nm_corr_worker, thA, S < 100 ? 100 : S);
  free_one(nm_work);
  free_one(nm_iwork);

  // Phase B: the distances.
  size_t n_dist = (size_t) S * (S - 1) / 2;
  SEXP   result = PROTECT(allocVector(REALSXP, (R_xlen_t) n_dist));
  nm_out = REAL(result);

  double wB  = (double) n_dist * (nm_mnorm == NMIT_NUCLEAR ? (double) D * D * D : (double) nm_plen);
  int    thB = threads_for(n_thr, wB, 1e6);
  if (thB > S - 1) thB = S - 1 < 1 ? 1 : S - 1;
  nm_wlen   = nm_mnorm == NMIT_NUCLEAR ? (size_t) 2 * D * D + 3 * (size_t) D : 1;
  nm_work   = (double *) safe_malloc((size_t) thB * nm_wlen * sizeof(double));
  nm_failed = (int *)    safe_malloc((size_t) thB * sizeof(int));
  for (int t = 0; t < thB; t++) nm_failed[t] = 0;
  if (S > 1) run_parallel(nm_dist_worker, thB, S - 1 < 100 ? 100 : S - 1);

  int failed = 0;
  for (int t = 0; t < thB; t++) failed |= nm_failed[t];

  free_all();

  if (failed) { // # nocov start
    UNPROTECT(1);
    error("An eigendecomposition for the nuclear norm did not converge.");
  } // # nocov end

  UNPROTECT(1);
  return result;
}



//======================================================
// R interface for volatility(): Vega's bootstrapCI() for
// each cell (a group at a state) of values.
//
// sexp_vals:  the values, cell after cell.
// sexp_pos:   n_cells + 1 offsets into sexp_vals.
// sexp_args:  n_resamples, seed.
// sexp_alpha: 1 - conf_level.
//
// For each cell (stream = its 0-based index): n_resamples
// means of n values drawn with replacement (index
// floor(u n), u from rng_unif()), each summed in draw
// order; sorted; then d3's (type 7) quantiles at alpha / 2
// and 1 - alpha / 2. A cell of no values gives NA.
// Returns a 2 x n_cells matrix: lower, upper.
//======================================================
static int cmp_double (const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

static double d3_quantile (const double *v, int n, double p) {
  if (n < 1) return NA_REAL;
  if (p <= 0 || n < 2) return v[0];
  if (p >= 1) return v[n - 1];
  double i  = (n - 1) * p;
  int    i0 = (int) floor(i);
  double v0 = v[i0], v1 = v[i0 + 1];
  return v0 + (v1 - v0) * (i - i0);
}

SEXP C_volatility_ci (SEXP sexp_vals, SEXP sexp_pos, SEXP sexp_args, SEXP sexp_alpha) {

  const double *vals  = REAL(sexp_vals);
  const int    *pos   = INTEGER(sexp_pos);
  int           cells = length(sexp_pos) - 1;
  int           m     = INTEGER(sexp_args)[0];
  int           seed  = INTEGER(sexp_args)[1];
  double        alpha = asReal(sexp_alpha);

  // Cell c draws from PCG stream base + c; base is args[2] * 2^32 when
  // given (feature_volatility(): one cell per feature, group and state),
  // else 0 (volatility()).
  uint64_t base = length(sexp_args) > 2 ? ((uint64_t) INTEGER(sexp_args)[2]) << 32 : 0;

  SEXP    result = PROTECT(allocMatrix(REALSXP, 2, cells));
  double *out    = REAL(result);
  double *mu     = (double *) R_alloc((size_t) (m < 1 ? 1 : m), sizeof(double));

  for (int c = 0; c < cells; c++) {

    int           n = pos[c + 1] - pos[c];
    const double *v = vals + pos[c];

    if (n < 1 || m < 1) {
      out[2 * c] = out[2 * c + 1] = NA_REAL;
      continue;
    }

    rng_t r;
    rng_seed(&r, (uint64_t) seed, base + (uint64_t) c);

    for (int j = 0; j < m; j++) {
      double a = 0;
      for (int i = 0; i < n; i++) {
        int k = (int) floor(rng_unif(&r) * n);
        if (k >= n) k = n - 1;   // # nocov
        a += v[k];
      }
      mu[j] = a / n;
    }

    qsort(mu, (size_t) m, sizeof(double), cmp_double);
    out[2 * c]     = d3_quantile(mu, m, alpha / 2);
    out[2 * c + 1] = d3_quantile(mu, m, 1 - alpha / 2);
  }

  UNPROTECT(1);
  return result;
}
