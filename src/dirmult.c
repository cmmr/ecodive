// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// dirmult_ttest() follows dirmult_ttest(), _dirmult_draw() and
// _welch_draw_stats() from scikit-bio 0.7.4
// (skbio/stats/composition/_dirmult.py), under this license:
//
//   Copyright (c) 2013--, scikit-bio development team.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   * Redistributions of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//
//   * Neither the names scikit-bio, skbio, or biocore nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
//   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
//   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Its scale model (scale_sd) follows clr.sm() from ALDEx3 1.4.0
// (R/scale.R), under this license:
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
// clr.sm() gives each sample the log2 scale -mean(log2 composition) + X' L,
// with one N(0, gamma^2) coefficient per row of the design X per Monte
// Carlo draw: samples with the same design share the shift. For the two
// sides of a contrast (X = intercept + treatment indicator), the
// intercept's shift cancels from every difference, and the treatment's,
// L ~ N(0, scale_sd^2) in log2 units, is added to each draw's difference
// of clr means (scaled by log 2 here, which works in natural logs); the
// variances, being within sides, don't change.
//
// The rest of this file is Copyright 2026 ecodive authors, under the MIT
// license above. The Gamma draws are rng_dist.c's (on PCG, not NumPy's
// PCG64 sampler, so the draws themselves differ from scikit-bio's); the
// t distribution is R's pt() and qt() (scikit-bio: SciPy's stdtr() and
// stdtrit()).
//
// For each draw d = 0 .. n_draws - 1:
//
//   Phase A (threads, one task per sample of the contrast): the sample's
//   row of log G, G_f ~ Gamma(x_f + pseudocount, 1), from its own PCG
//   stream (seed, d * (n + 1) + s), features in order; then the clr, by
//   subtracting the row's mean. (Or, for the tests, the log of given G.)
//
//   Phase B (threads, one task per block of features): Welch's difference
//   (treatment minus reference), standard error and degrees of freedom,
//   in scikit-bio's order of operations.
//
//   Then, on this thread (R's pt() and qt() are not thread-safe when they
//   warn): the scale model's shift of every difference (one normal draw
//   from stream (seed, d * (n + 1) + n)), t, the one-sided tail probability, the interval's bounds, and
//   the running sums, minima and maxima across draws, in draw order, as
//   scikit-bio's loop.
//
// Every task depends only on its own sample or features, and the sums run
// in draw order on one thread, so results are identical for any number of
// threads. Memory: one (samples in the contrast) x D buffer of doubles.

#include "ecodive.h"
#include <Rmath.h>

#define DM_BLOCK 64     // features per phase B task

static int           dm_D;
static int           dm_m;          // samples in the contrast
static const int    *dm_sam;        // their indices in the table
static const int    *dm_side;       // 1 treatment, 2 reference, per contrast sample
static int           dm_n1, dm_n2;
static const int    *dm_pos;        // the table, by sample (ecomatrix)
static const int    *dm_otu;
static const double *dm_val;
static double        dm_pc;
static uint64_t      dm_seed;
static uint64_t      dm_n;          // samples in the table
static int           dm_draw;
static const double *dm_given;      // n x D x n_draws given G, or NULL
static double       *dm_clr;        // m x D, row by row
static double       *dm_diff, *dm_se, *dm_dof;


// Phase A: each sample's row of clr values.
static void *draw_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D = dm_D;

  for (int i = thread_i; i < dm_m; i += n_threads) {

    int     s   = dm_sam[i];
    double *row = dm_clr + (size_t) i * D;

    if (dm_given) {
      const double *g = dm_given + (size_t) dm_draw * dm_n * D + s;
      for (int f = 0; f < D; f++) row[f] = log(g[(size_t) f * dm_n]);
    }
    else {
      // The sample's values, densely, then the draws in feature order, so
      // that a sparse table gives the same numbers as a dense one.
      for (int f = 0; f < D; f++) row[f] = 0;
      for (int e = dm_pos[s]; e < dm_pos[s + 1]; e++) row[dm_otu[e]] = dm_val[e];

      rng_t rng;
      rng_seed(&rng, dm_seed, (uint64_t) dm_draw * (dm_n + 1) + (uint64_t) s);
      for (int f = 0; f < D; f++) row[f] = rng_log_gamma(&rng, row[f] + dm_pc);
    }

    double sum = 0;
    for (int f = 0; f < D; f++) sum += row[f];
    double mean = sum / D;
    for (int f = 0; f < D; f++) row[f] -= mean;
  }
  return NULL;
}


// Phase B: Welch's difference, standard error and degrees of freedom, per
// feature, as _welch_draw_stats(): means, then variances (ddof = 1) as
// sums of squared deviations from them, summed in sample order.
static void *welch_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D = dm_D, m = dm_m;
  double n1 = dm_n1, n2 = dm_n2;

  double s1[DM_BLOCK], s2[DM_BLOCK], q1[DM_BLOCK], q2[DM_BLOCK];
  int n_blocks = (D + DM_BLOCK - 1) / DM_BLOCK;

  for (int b = thread_i; b < n_blocks; b += n_threads) {

    int f0 = b * DM_BLOCK;
    int nf = D - f0 < DM_BLOCK ? D - f0 : DM_BLOCK;

    for (int k = 0; k < nf; k++) s1[k] = s2[k] = q1[k] = q2[k] = 0;

    for (int i = 0; i < m; i++) {
      const double *row = dm_clr + (size_t) i * D + f0;
      double       *acc = dm_side[i] == 1 ? s1 : s2;
      for (int k = 0; k < nf; k++) acc[k] += row[k];
    }
    for (int k = 0; k < nf; k++) { s1[k] /= n1; s2[k] /= n2; }   // means

    for (int i = 0; i < m; i++) {
      const double *row  = dm_clr + (size_t) i * D + f0;
      const double *mean = dm_side[i] == 1 ? s1 : s2;
      double       *acc  = dm_side[i] == 1 ? q1 : q2;
      for (int k = 0; k < nf; k++) { double r = row[k] - mean[k]; acc[k] += r * r; }
    }

    for (int k = 0; k < nf; k++) {
      double vn1 = (q1[k] / (n1 - 1)) / n1;
      double vn2 = (q2[k] / (n2 - 1)) / n2;
      double var_sum = vn1 + vn2;
      double den = (vn1 * vn1) / (n1 - 1) + (vn2 * vn2) / (n2 - 1);
      dm_diff[f0 + k] = s1[k] - s2[k];
      dm_dof[f0 + k]  = (var_sum * var_sum) / den;
      dm_se[f0 + k]   = sqrt(var_sum);
    }
  }
  return NULL;
}


// NaN-propagating minimum and maximum, as numpy.minimum() and maximum().
static inline double nan_min (double a, double b) {
  return (ISNAN(a) || ISNAN(b)) ? R_NaN : (b < a ? b : a);
}
static inline double nan_max (double a, double b) {
  return (ISNAN(a) || ISNAN(b)) ? R_NaN : (b > a ? b : a);
}


//======================================================
// R interface for dirmult_ttest().
//
// sexp_side:    per sample of the table: 1 treatment,
//               2 reference, 0 not in the contrast.
// sexp_prob:    the t quantile of the interval (0.975 for
//               a 95% interval).
// sexp_scale:   the scale model's sd, in log2 units (0 for
//               none).
// sexp_given:   NULL, or n x D x n_draws Gamma draws to use
//               instead of sampling (tests).
// sexp_details: TRUE also returns each draw's difference,
//               se, dof, t, tail probability and bounds
//               (n_draws x D matrices, natural logs), and
//               each draw's scale shift.
//
// Returns list(t, log2_fc, ci_lower, ci_upper, p_value[,
// draws, shift]).
//======================================================
SEXP C_dirmult_ttest(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_side, SEXP sexp_pseudocount,
    SEXP sexp_n_draws, SEXP sexp_prob, SEXP sexp_scale, SEXP sexp_seed,
    SEXP sexp_given, SEXP sexp_details, SEXP sexp_n_threads ) {

  double pc      = asReal(sexp_pseudocount);
  int    n_draws = asInteger(sexp_n_draws);
  double prob    = asReal(sexp_prob);
  double sc_sd   = asReal(sexp_scale);
  int    details = asLogical(sexp_details);
  int    n_thr   = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  // The table, copied out of the ecomatrix into R's memory, so that an
  // interrupt between draws leaks nothing.
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

  // The samples in the contrast, in table order.
  const int *side_all = INTEGER(sexp_side);
  int m = 0, n1 = 0, n2 = 0;
  for (int s = 0; s < n; s++) if (side_all[s]) m++;
  int *sam  = (int *) R_alloc((size_t) m + 1, sizeof(int));
  int *side = (int *) R_alloc((size_t) m + 1, sizeof(int));
  m = 0;
  for (int s = 0; s < n; s++)
    if (side_all[s]) {
      sam[m] = s; side[m++] = side_all[s];
      if (side_all[s] == 1) n1++; else n2++;
    }

  double *clr  = (double *) R_alloc((size_t) m * D + 1, sizeof(double));
  double *diff = (double *) R_alloc((size_t) D, sizeof(double));
  double *se   = (double *) R_alloc((size_t) D, sizeof(double));
  double *dof  = (double *) R_alloc((size_t) D, sizeof(double));

  dm_D     = D;
  dm_m     = m;
  dm_sam   = sam;
  dm_side  = side;
  dm_n1    = n1;
  dm_n2    = n2;
  dm_pos   = pos;
  dm_otu   = otu;
  dm_val   = val;
  dm_pc    = pc;
  dm_seed  = (uint64_t) asInteger(sexp_seed);
  dm_n     = (uint64_t) n;
  dm_given = isNull(sexp_given) ? NULL : REAL(sexp_given);
  dm_clr   = clr;
  dm_diff  = diff;
  dm_se    = se;
  dm_dof   = dof;

  // Threads only for enough work: about 30 ns per Gamma draw, 3 ns per
  // value in phase B.
  int nt_a = threads_for(n_thr, (double) m * D * 30, 1e6);
  int nt_b = threads_for(n_thr, (double) m * D * 3,  1e6);
  int n_blocks = (D + DM_BLOCK - 1) / DM_BLOCK;

  // Running sums, minima and maxima, as scikit-bio.
  SEXP sexp_t     = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_delta = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_lower = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_upper = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_pval  = PROTECT(allocVector(REALSXP, D));
  double *tstat_sum = REAL(sexp_t),     *delta_sum = REAL(sexp_delta);
  double *lower     = REAL(sexp_lower), *upper     = REAL(sexp_upper);
  double *pval_sum  = REAL(sexp_pval);
  for (int f = 0; f < D; f++) {
    tstat_sum[f] = delta_sum[f] = pval_sum[f] = 0;
    lower[f] = R_PosInf;
    upper[f] = R_NegInf;
  }

  SEXP sexp_draws = R_NilValue;
  double *d_diff = NULL, *d_se = NULL, *d_dof = NULL, *d_t = NULL,
         *d_p = NULL, *d_lo = NULL, *d_hi = NULL;
  if (details) {
    sexp_draws = allocVector(VECSXP, 7);
    PROTECT(sexp_draws);
    const char *nms[] = { "diff", "se", "dof", "t", "p_tail", "lower", "upper" };
    double    **ptr[] = { &d_diff, &d_se, &d_dof, &d_t, &d_p, &d_lo, &d_hi };
    SEXP names = PROTECT(allocVector(STRSXP, 7));
    for (int k = 0; k < 7; k++) {
      SEXP mat = allocMatrix(REALSXP, n_draws, D);
      SET_VECTOR_ELT(sexp_draws, k, mat);
      *ptr[k] = REAL(mat);
      SET_STRING_ELT(names, k, mkChar(nms[k]));
    }
    setAttrib(sexp_draws, R_NamesSymbol, names);
    UNPROTECT(1);
  }
  else {
    PROTECT(sexp_draws);
  }

  SEXP sexp_shift = PROTECT(allocVector(REALSXP, n_draws));
  double *shift_d = REAL(sexp_shift);
  double  log2_   = log(2.0);

  for (int d = 0; d < n_draws; d++) {

    dm_draw = d;
    run_parallel(draw_worker,  nt_a, m < 100 ? 100 : m);
    run_parallel(welch_worker, nt_b, n_blocks < 100 ? 100 : n_blocks);

    // The scale model: the same shift of every difference in this draw.
    double shift = 0;
    if (sc_sd > 0) {
      rng_t rng;
      rng_seed(&rng, dm_seed, (uint64_t) d * (dm_n + 1) + dm_n);
      shift = rng_norm(&rng) * sc_sd * log2_;
      for (int f = 0; f < D; f++) diff[f] += shift;
    }
    shift_d[d] = shift;

    for (int f = 0; f < D; f++) {

      double t = diff[f] / se[f];
      tstat_sum[f] += t;

      // The one-sided tail; doubled once, after the loop. No pt() or qt()
      // for degrees of freedom that aren't positive (NaN from a side of
      // one sample, or a feature constant on both sides).
      double p = R_NaN, q = R_NaN;
      if (dof[f] > 0) {
        p = pt(-fabs(t), dof[f], 1, 0);
        q = qt(prob, dof[f], 1, 0);
      }
      pval_sum[f] += p;

      double margin = q * se[f];
      double lo = diff[f] - margin;
      double hi = diff[f] + margin;
      lower[f] = nan_min(lower[f], lo);
      upper[f] = nan_max(upper[f], hi);

      delta_sum[f] += diff[f];

      if (details) {
        size_t k = (size_t) f * n_draws + d;
        d_diff[k] = diff[f]; d_se[k] = se[f]; d_dof[k] = dof[f];
        d_t[k] = t; d_p[k] = p; d_lo[k] = lo; d_hi[k] = hi;
      }
    }

    if (d + 1 < n_draws && user_interrupted()) {
      UNPROTECT(7);
      error("Interrupted.");
    }
  }

  // Averages across draws, then base 2, in scikit-bio's order.
  for (int f = 0; f < D; f++) {
    delta_sum[f] /= n_draws;
    tstat_sum[f] /= n_draws;
    pval_sum[f]  *= 2.0 / n_draws;
    delta_sum[f] /= log2_;
    upper[f]     /= log2_;
    lower[f]     /= log2_;
  }

  int  n_out  = details ? 7 : 5;
  SEXP result = PROTECT(allocVector(VECSXP, n_out));
  SEXP names  = PROTECT(allocVector(STRSXP, n_out));
  SET_VECTOR_ELT(result, 0, sexp_t);     SET_STRING_ELT(names, 0, mkChar("t"));
  SET_VECTOR_ELT(result, 1, sexp_delta); SET_STRING_ELT(names, 1, mkChar("log2_fc"));
  SET_VECTOR_ELT(result, 2, sexp_lower); SET_STRING_ELT(names, 2, mkChar("ci_lower"));
  SET_VECTOR_ELT(result, 3, sexp_upper); SET_STRING_ELT(names, 3, mkChar("ci_upper"));
  SET_VECTOR_ELT(result, 4, sexp_pval);  SET_STRING_ELT(names, 4, mkChar("p_value"));
  if (details) {
    SET_VECTOR_ELT(result, 5, sexp_draws);
    SET_STRING_ELT(names, 5, mkChar("draws"));
    SET_VECTOR_ELT(result, 6, sexp_shift);
    SET_STRING_ELT(names, 6, mkChar("shift"));
  }
  setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(9);
  return result;
}
