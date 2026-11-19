// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The variance of log ratios, vlr(), follows vlr() and pairwise_vlr() from
// scikit-bio 0.7.4 (skbio/stats/composition/_base.py), under this license:
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
// The rest of this file is Copyright 2026 ecodive authors, under the MIT
// license above. It computes each pair's variance from the differences of
// centred logs, not from scikit-bio's covariance matrix (c_ii + c_jj -
// 2 c_ij), which cancels badly for nearly proportional pairs; rho_p() and
// phi_s() (Erb & Notredame 2016; Quinn et al. 2017) are formulas.
//
// With L the logs of (x + pseudocount), each feature centred over the
// samples (a_sf = L_sf - mean_s L_sf), and z = L - g the log-ratio
// transform against a per-sample reference g (the mean log of all features
// for the clr, or of the `denominator` features), centred too:
//
//   SS_ij = sum_s (a_si - a_sj)^2       vlr_ij   = SS_ij / (n - 1)
//   V_i   = sum_s (a_si - g_s)^2        rho_p_ij = 1 - SS_ij / (V_i + V_j)
//                                       phi_s_ij = SS_ij / (2 V_i + 2 V_j - SS_ij)
//
// The table is kept sparse: features are transposed into lists of their
// non-zero samples (centred logs), and every zero of feature f has the same
// centred log, log(pseudocount) - mean_f. A pair walks its two lists in
// sample order and adds the samples where both are zero at the end, as
// their count times one squared difference.
//
// Each task is one row i of the packed lower triangle (the pairs i, j > i),
// rows dealt to threads in turn. A pair's sum is always taken in the same
// order and pairs share nothing, so results are identical for any number
// of threads.

#include "ecodive.h"

#define PROP_VLR   1
#define PROP_RHO_P 2
#define PROP_PHI_S 3

static int           pr_method;
static int           pr_n;        // samples
static int           pr_K;        // features kept in the result
static const int    *pr_keep;     // 0-based feature of each kept one
static const int    *pr_fpos;     // per feature: start of its list (D + 1)
static const int    *pr_fsam;     // sample of each listed value
static const double *pr_fval;     // its centred log
static const double *pr_fzero;    // per feature: the centred log of a zero
static const double *pr_V;        // per feature: V_f (rho_p, phi_s)
static double       *pr_out;


// SS_ij: the sum of squared differences of two features' centred logs.
static double pair_ss(int fi, int fj) {

  const int    *i = pr_fsam + pr_fpos[fi], *i_end = pr_fsam + pr_fpos[fi + 1];
  const int    *j = pr_fsam + pr_fpos[fj], *j_end = pr_fsam + pr_fpos[fj + 1];
  const double *a = pr_fval + pr_fpos[fi];
  const double *b = pr_fval + pr_fpos[fj];
  double zi = pr_fzero[fi], zj = pr_fzero[fj];
  double ss = 0, d;
  int    n_union = 0;

  while (i != i_end || j != j_end) {
    if (j == j_end || (i != i_end && *i < *j)) { d = *a - zj; i++; a++; }
    else if (i == i_end || *j < *i)            { d = zi - *b; j++; b++; }
    else                                       { d = *a - *b; i++; a++; j++; b++; }
    ss += d * d;
    n_union++;
  }

  if (n_union < pr_n) {
    d   = zi - zj;
    ss += (double) (pr_n - n_union) * (d * d);
  }
  return ss;
}


static void *pair_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int K         = pr_K;
  double n1     = pr_n - 1;

  for (int r = thread_i; r < K - 1; r += n_threads) {

    // Row r's pairs start after those of rows 0 .. r - 1.
    size_t idx = (size_t) r * K - (size_t) r * (r + 1) / 2;
    int    fi  = pr_keep[r];

    for (int c = r + 1; c < K; c++, idx++) {
      int    fj = pr_keep[c];
      double ss = pair_ss(fi, fj);
      double v;
      switch (pr_method) {
        case PROP_RHO_P: v = 1 - ss / (pr_V[fi] + pr_V[fj]);               break;
        case PROP_PHI_S: {
          // var(z_i + z_j); not above 0 only for reciprocal pairs, by rounding.
          double den = 2 * pr_V[fi] + 2 * pr_V[fj] - ss;
          v = den > 0 ? ss / den : ss > 0 ? R_PosInf : R_NaN;
          break;
        }
        default:         v = ss / n1;
      }
      pr_out[idx] = v;
    }
  }
  return NULL;
}


//======================================================
// R interface for vlr(), rho_p() and phi_s(): the packed
// lower triangle over the kept features, as a `dist`.
//
// sexp_method: 1 vlr, 2 rho_p, 3 phi_s.
// sexp_ref:    0-based features whose mean log is the
//              reference g, or NULL for all (the clr).
//              Unused by vlr.
// sexp_keep:   0-based features in the result, increasing.
// sexp_labels: their names, or NULL.
//======================================================
SEXP C_proportionality(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_method, SEXP sexp_ref,
    SEXP sexp_keep, SEXP sexp_labels, SEXP sexp_pseudocount,
    SEXP sexp_n_threads ) {

  int    method = asInteger(sexp_method);
  double pc     = asReal(sexp_pseudocount);
  int    n_thr  = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n   = em->n_samples;
  int D   = em->n_otus;
  int nnz = em->nnz;

  if (pc == 0 && (double) nnz < (double) n * D) {
    free_all();
    errorcall(R_NilValue, "`counts` has zeros, and `pseudocount = 0` leaves their logs undefined. Use `pseudocount = NULL` or a positive pseudocount, or replace_zeros() first.");
  }

  double lz = pc > 0 ? log(pc) : 0;   // log of a zero (none when pc is 0)

  // Transpose into per-feature lists, samples increasing.
  int    *fpos  = (int *)    R_alloc((size_t) D + 1, sizeof(int));
  int    *fsam  = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  double *fval  = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  double *fzero = (double *) R_alloc((size_t) D, sizeof(double));
  int    *fill  = (int *)    R_alloc((size_t) D, sizeof(int));

  for (int f = 0; f <= D; f++) fpos[f] = 0;
  for (int e = 0; e < nnz; e++) fpos[em->otu_vec[e] + 1]++;
  for (int f = 0; f < D; f++)   fpos[f + 1] += fpos[f];
  for (int f = 0; f < D; f++)   fill[f] = fpos[f];

  for (int s = 0; s < n; s++)
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++) {
      int f = em->otu_vec[e];
      fsam[fill[f]]   = s;
      fval[fill[f]++] = log(em->val_vec[e] + pc);
    }

  // Centre each feature over the samples, summing in sample order.
  for (int f = 0; f < D; f++) {
    double sum = 0;
    int    e   = fpos[f], e_end = fpos[f + 1];
    for (int s = 0; s < n; s++) {
      if (e < e_end && fsam[e] == s) sum += fval[e++];
      else                           sum += lz;
    }
    double mean = sum / n;
    for (e = fpos[f]; e < e_end; e++) fval[e] -= mean;
    fzero[f] = lz - mean;
  }

  // V_f for rho_p and phi_s.
  double *V = NULL;
  if (method != PROP_VLR) {

    // The reference g: each sample's mean log over the reference features,
    // summed in feature order, then centred.
    int        n_ref = isNull(sexp_ref) ? D : length(sexp_ref);
    const int *ref   = isNull(sexp_ref) ? NULL : INTEGER(sexp_ref);
    double    *row   = (double *) R_alloc((size_t) D, sizeof(double));
    double    *g     = (double *) R_alloc((size_t) n, sizeof(double));
    double     gsum  = 0;

    for (int s = 0; s < n; s++) {
      for (int f = 0; f < D; f++) row[f] = lz;
      for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++)
        row[em->otu_vec[e]] = log(em->val_vec[e] + pc);
      double sum = 0;
      for (int k = 0; k < n_ref; k++) sum += row[ref ? ref[k] : k];
      g[s]  = sum / n_ref;
      gsum += g[s];
    }
    double gmean = gsum / n;
    for (int s = 0; s < n; s++) g[s] -= gmean;

    V = (double *) R_alloc((size_t) D, sizeof(double));
    for (int f = 0; f < D; f++) {
      double ss = 0, d;
      int    e  = fpos[f], e_end = fpos[f + 1];
      for (int s = 0; s < n; s++) {
        if (e < e_end && fsam[e] == s) d = fval[e++] - g[s];
        else                           d = fzero[f]  - g[s];
        ss += d * d;
      }
      V[f] = ss;
    }
  }

  int K = length(sexp_keep);
  size_t n_dist = (size_t) K * (K - 1) / 2;

  SEXP result = PROTECT(allocVector(REALSXP, (R_xlen_t) n_dist));

  pr_method = method;
  pr_n      = n;
  pr_K      = K;
  pr_keep   = INTEGER(sexp_keep);
  pr_fpos   = fpos;
  pr_fsam   = fsam;
  pr_fval   = fval;
  pr_fzero  = fzero;
  pr_V      = V;
  pr_out    = REAL(result);

  // Threads only for enough work (about a millisecond each).
  double work = (double) n_dist * (n < 1 ? 1 : n);
  run_parallel(pair_worker, threads_for(n_thr, work, 1e6), K - 1 < 100 ? 100 : K - 1);

  // Each value protected while install() may allocate (rchk).
  SEXP sexp_size  = PROTECT(ScalarInteger(K));
  SEXP sexp_diag  = PROTECT(ScalarLogical(0));
  SEXP sexp_upper = PROTECT(ScalarLogical(0));
  SEXP sexp_class = PROTECT(mkString("dist"));

  setAttrib(result, install("Size"),   sexp_size);
  setAttrib(result, install("Labels"), sexp_labels);
  setAttrib(result, install("Diag"),   sexp_diag);
  setAttrib(result, install("Upper"),  sexp_upper);
  setAttrib(result, R_ClassSymbol,     sexp_class);

  free_all();
  UNPROTECT(5);
  return result;
}
