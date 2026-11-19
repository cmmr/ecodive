// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// balance() gives one coordinate of ilr() with a sequential binary
// partition, computed as ilr() computes it (see log_ratio_transforms.c),
// which follows ilr() and sbp_basis() from scikit-bio 0.7.4
// (skbio/stats/composition/_base.py), under this license:
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
// The rest of this file, including log_ratio() (the log of the ratio of
// summed abundances, as Qurro's sample plots show it), is Copyright 2026
// ecodive authors, under the MIT license above.
//
// One value per sample, single-threaded: each sample's stored values are
// spread into a dense scratch vector, and only the entries written are
// cleared afterwards, so a sparse table costs its non-zeros plus the
// selected features.
//
//   log_ratio: log(sum_N (x + c) / sum_D (x + c)), the sums taken in feature
//              order; NA when either sum is 0.
//   balance:   the inner product of clr(x + c) (all D features, as
//              normalize(norm = 'clr')) with one SBP row converted by
//              sbp_basis(), its entries in feature order, so it equals the
//              ilr() coordinate to the last bit. NA when a selected feature
//              is 0 (pseudocount 0). With zeros only outside the selection,
//              the clr is undefined but the balance is not: it is then the
//              same inner product with the uncentred logs (the basis row
//              sums to 0, so the centring cancels).

#include "ecodive.h"


// Sample s's values (0 where not stored) into the dense vector v; the
// features written are listed in `touched` and their count returned.
static int spread(const ecomatrix_t *em, int s, double *v, int *touched) {
  int k = 0;
  for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++) {
    v[em->otu_vec[e]] = em->val_vec[e];
    touched[k++]      = em->otu_vec[e];
  }
  return k;
}

static void clear(double *v, const int *touched, int k) {
  for (int t = 0; t < k; t++) v[touched[t]] = 0;
}


//======================================================
// R interface for log_ratio(). sexp_num and sexp_den:
// the 0-based features of each side, increasing.
//======================================================
SEXP C_log_ratio(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_num,
    SEXP sexp_den, SEXP sexp_pseudocount ) {

  double pc  = asReal(sexp_pseudocount);
  int    nn  = length(sexp_num), nd = length(sexp_den);
  const int *num = INTEGER(sexp_num), *den = INTEGER(sexp_den);

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int ns = em->n_samples, D = em->n_otus;

  double *v       = (double *) R_alloc((size_t) D, sizeof(double));
  int    *touched = (int *)    R_alloc((size_t) D, sizeof(int));
  memset(v, 0, sizeof(double) * (size_t) D);

  SEXP result = PROTECT(allocVector(REALSXP, ns));
  double *out = REAL(result);

  for (int s = 0; s < ns; s++) {
    int k = spread(em, s, v, touched);
    double sn = 0, sd = 0;
    for (int j = 0; j < nn; j++) sn += v[num[j]] + pc;
    for (int j = 0; j < nd; j++) sd += v[den[j]] + pc;
    out[s] = (sn > 0 && sd > 0) ? log(sn / sd) : NA_REAL;
    clear(v, touched, k);
  }

  free_all();
  UNPROTECT(1);
  return result;
}


//======================================================
// R interface for balance(). sexp_feat and sexp_coef:
// the basis row's non-zero entries (0-based features,
// increasing, and their values), as ilr_basis()'s CSR.
//======================================================
SEXP C_balance(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_feat,
    SEXP sexp_coef, SEXP sexp_pseudocount ) {

  double pc = asReal(sexp_pseudocount);
  int    nf = length(sexp_feat);
  const int    *feat = INTEGER(sexp_feat);
  const double *coef = REAL(sexp_coef);

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int ns = em->n_samples, D = em->n_otus;

  double *c = (double *) R_alloc((size_t) D, sizeof(double));

  SEXP result = PROTECT(allocVector(REALSXP, ns));
  double *out = REAL(result);

  double lz = log(pc);  // -Inf for pseudocount 0

  for (int s = 0; s < ns; s++) {

    int e0 = em->pos_vec[s], e1 = em->pos_vec[s + 1];
    if (e1 - e0 < D)
      for (int o = 0; o < D; o++) c[o] = lz;
    for (int e = e0; e < e1; e++)
      c[em->otu_vec[e]] = log(em->val_vec[e] + pc);

    int sel_zero = 0;
    for (int j = 0; j < nf; j++)
      if (!R_FINITE(c[feat[j]])) { sel_zero = 1; break; }
    if (sel_zero) { out[s] = NA_REAL; continue; }

    double logsum = 0;
    for (int o = 0; o < D; o++) logsum += c[o];

    if (R_FINITE(logsum)) {
      double log_mean = logsum / D;
      for (int o = 0; o < D; o++) c[o] -= log_mean;
    }

    double b = 0;
    for (int j = 0; j < nf; j++) b += coef[j] * c[feat[j]];
    out[s] = b;
  }

  free_all();
  UNPROTECT(1);
  return result;
}
