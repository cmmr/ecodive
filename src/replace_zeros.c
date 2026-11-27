// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Multiplicative zero replacement, replace_zeros_mult(), follows
// multi_replace() and closure() from scikit-bio 0.7.4
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
// np_pairwise_sum() follows pairwise_sum_DOUBLE() from numpy 2.x
// (numpy/_core/src/umath/loops_utils.h.src), under this license:
//
//   Copyright (c) 2005-2025, NumPy Developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided
//     with the distribution.
//
//   * Neither the name of the NumPy Developers nor the names of any
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
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
// The zero replacements read any table type through the ecomatrix and
// return a dense matrix in the input's orientation. Each sample is expanded
// into a dense scratch vector, so its total is summed in numpy's pairwise
// order (np.sum over a contiguous row), and the arithmetic is
// multi_replace()'s, step for step: results are bit-identical to
// scikit-bio's.
//
// The count zero multiplicative and Bayesian-multiplicative replacements
// (C_replace_zeros_cmult) are a clean-room implementation: their formulas
// were established from zCompositions' documentation and from the outputs
// of its cmultRepl() (GPL; its source was not read), as recorded in
// ecodive's plans/zero_replacement.md.

#include "ecodive.h"


// numpy's pairwise summation (pairwise_sum_DOUBLE, see the license above),
// which np.sum() uses along a
// contiguous axis: sequential below 8 values, eight running sums up to 128,
// and halves (at a multiple of 8) above.
static double np_pairwise_sum(const double *a, int n) {

  if (n < 8) {
    double res = 0.;
    for (int i = 0; i < n; i++) res += a[i];
    return res;
  }

  if (n <= 128) {
    double r[8];
    int i;
    for (int j = 0; j < 8; j++) r[j] = a[j];
    for (i = 8; i < n - (n % 8); i += 8)
      for (int j = 0; j < 8; j++) r[j] += a[i + j];
    double res = ((r[0] + r[1]) + (r[2] + r[3])) +
                 ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; i++) res += a[i];
    return res;
  }

  int n2 = n / 2;
  n2 -= n2 % 8;
  return np_pairwise_sum(a, n2) + np_pairwise_sum(a + n2, n - n2);
}


// Sample s of `em` as a dense vector of all its features.
static void dense_sample(const ecomatrix_t *em, int s, double *x) {
  memset(x, 0, sizeof(double) * (size_t) em->n_otus);
  for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++)
    x[em->otu_vec[e]] += em->val_vec[e];
}


//======================================================
// R interface: the zeros of a table, for the checks that
// every zero replacement shares. A list of
//   sample_zeros:     the number of zeros in each sample;
//   feature_nonzeros: the number of samples in which each
//                     feature is non-zero.
//======================================================
SEXP C_zero_counts(SEXP sexp_counts, SEXP sexp_margin) {

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int ns = em->n_samples;
  int D  = em->n_otus;

  SEXP sexp_sz = PROTECT(allocVector(INTSXP, ns));
  SEXP sexp_fn = PROTECT(allocVector(INTSXP, D));
  int *sz = INTEGER(sexp_sz);
  int *fn = INTEGER(sexp_fn);
  memset(fn, 0, sizeof(int) * (size_t) D);

  double *x = (double *) R_alloc((size_t) D, sizeof(double));

  for (int s = 0; s < ns; s++) {
    dense_sample(em, s, x);
    int nz = 0;
    for (int o = 0; o < D; o++) {
      if (x[o] == 0) nz++;
      else           fn[o]++;
    }
    sz[s] = nz;
  }

  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SEXP names  = PROTECT(allocVector(STRSXP, 2));
  SET_VECTOR_ELT(result, 0, sexp_sz);
  SET_VECTOR_ELT(result, 1, sexp_fn);
  SET_STRING_ELT(names, 0, mkChar("sample_zeros"));
  SET_STRING_ELT(names, 1, mkChar("feature_nonzeros"));
  setAttrib(result, R_NamesSymbol, names);

  free_all();
  UNPROTECT(4);
  return result;
}


//======================================================
// R interface: multi_replace()'s default delta,
// (1.0 / D) ** 2. Python's ** calls the C library's
// pow(), which is not always x * x (R's ^2), so it is
// called here; `two` is volatile so that the compiler
// can't turn pow(x, 2) into x * x.
//======================================================
SEXP C_mult_default_delta(SEXP sexp_D) {
  volatile double two = 2.0;
  return ScalarReal(pow(1.0 / (double) asInteger(sexp_D), two));
}


//======================================================
// R interface for replace_zeros_mult(): scikit-bio's
// multi_replace() of every sample,
//
//   x = counts / sum(counts)          (closure)
//   r = delta                         where x == 0
//   r = (1 - n_zeros * delta) * x     elsewhere
//
// sexp_delta: one value, or one per sample. R has
// checked that 1 - n_zeros * delta > 0 and that no
// sample is all zeros.
//
// Returns a dense matrix in the input's orientation:
// samples x features for margin 1, features x samples
// for margin 2.
//======================================================
SEXP C_replace_zeros_mult(SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_delta) {

  int margin = asInteger(sexp_margin);

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int     ns      = em->n_samples;
  int     D       = em->n_otus;
  double *delta   = REAL(sexp_delta);
  int     n_delta = length(sexp_delta);

  SEXP result = PROTECT(margin == 1 ? allocMatrix(REALSXP, ns, D) : allocMatrix(REALSXP, D, ns));
  double *out = REAL(result);

  double *x = (double *) R_alloc((size_t) D, sizeof(double));

  for (int s = 0; s < ns; s++) {

    dense_sample(em, s, x);

    // closure(), then its zeros (as skbio, of the proportions).
    double tot = np_pairwise_sum(x, D);
    int    nz  = 0;
    for (int o = 0; o < D; o++) {
      x[o] = x[o] / tot;
      if (x[o] == 0) nz++;
    }

    double d     = delta[n_delta == 1 ? 0 : s];
    double zcnts = 1 - (double) nz * d;

    for (int o = 0; o < D; o++) {
      double v = (x[o] == 0) ? d : zcnts * x[o];
      if (margin == 1) out[s + (size_t) ns * o] = v;
      else             out[o + (size_t) D * s]  = v;
    }
  }

  free_all();
  UNPROTECT(1);
  return result;
}


//======================================================
// R interface for replace_zeros_czm(), _bl(), _sq() and
// _gbm(). With c the counts, n_i = sum_j c_ij the depth
// of sample i and x_ij = c_ij / n_i:
//
//   T_ij = sum_{k != i} c_kj / sum_{k != i} n_k
//          (the other samples' pooled proportion)
//
//   code 0, CZM:  r_ij = frac * threshold / n_i
//   code 1, BL:   r_ij = s_i T_ij / (n_i + s_i), s_i = D
//   code 2, SQ:   ... s_i = sqrt(n_i)
//   code 3, GBM:  ... s_i = 1 / exp(mean_j log T_ij)
//
// for every zero; with `adjust`, an r_ij above
// m_j = min_{k: c_kj > 0} x_kj (strictly) becomes
// frac * m_j. The non-zero parts are x_ij * (1 - the
// sample's sum of r_ij).
//
// R has checked that no sample is all zeros and, for
// the Bayesian codes, that there are at least two
// samples and no feature is zero everywhere (and, for
// GBM, that every feature is non-zero in at least two
// samples, so T_ij > 0).
//
// Every sample is expanded into a dense vector before
// it is summed, so dense and sparse input give the same
// bits. Returns a list of
//   result:     the dense matrix in the input's
//               orientation (samples x features for
//               margin 1, features x samples for 2);
//   kept:       each sample's 1 - sum of r_ij, the
//               factor of its non-zero parts;
//   n_adjusted: the number of imputed values adjusted.
//======================================================
SEXP C_replace_zeros_cmult(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_code,
    SEXP sexp_frac,   SEXP sexp_threshold, SEXP sexp_adjust ) {

  int    margin    = asInteger(sexp_margin);
  int    code      = asInteger(sexp_code);
  double frac      = asReal(sexp_frac);
  double threshold = asReal(sexp_threshold);
  int    adjust    = asLogical(sexp_adjust);

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int ns = em->n_samples;
  int D  = em->n_otus;

  double *x     = (double *) R_alloc((size_t) D,  sizeof(double));
  double *depth = (double *) R_alloc((size_t) ns, sizeof(double));
  double *csum  = (double *) R_alloc((size_t) D,  sizeof(double));
  double *cmin  = (double *) R_alloc((size_t) D,  sizeof(double));
  double  total = 0;

  for (int o = 0; o < D; o++) { csum[o] = 0; cmin[o] = R_PosInf; }

  // Pass 1: depths, feature totals, and each feature's smallest observed
  // proportion. Summed in sample order, then feature order.
  for (int s = 0; s < ns; s++) {
    dense_sample(em, s, x);
    double n = 0;
    for (int o = 0; o < D; o++) n += x[o];
    depth[s] = n;
    total   += n;
    for (int o = 0; o < D; o++) {
      if (x[o] == 0) continue;
      csum[o] += x[o];
      double p = x[o] / n;
      if (p < cmin[o]) cmin[o] = p;
    }
  }

  SEXP sexp_res  = PROTECT(margin == 1 ? allocMatrix(REALSXP, ns, D) : allocMatrix(REALSXP, D, ns));
  SEXP sexp_kept = PROTECT(allocVector(REALSXP, ns));
  double *out  = REAL(sexp_res);
  double *kept = REAL(sexp_kept);
  double  n_adjusted = 0;

  double *r = (double *) R_alloc((size_t) D, sizeof(double));

  // Pass 2: each sample.
  for (int s = 0; s < ns; s++) {

    dense_sample(em, s, x);
    double n    = depth[s];
    double rest = total - n;  // the other samples' depths

    double str = 0;  // s_i
    if      (code == 1) str = (double) D;
    else if (code == 2) str = sqrt(n);
    else if (code == 3) {
      double lsum = 0;
      for (int o = 0; o < D; o++) lsum += log((csum[o] - x[o]) / rest);
      str = 1 / exp(lsum / D);
    }

    double imputed = 0;
    for (int o = 0; o < D; o++) {
      if (x[o] != 0) continue;
      double v = (code == 0)
        ? frac * threshold / n
        : str * ((csum[o] - x[o]) / rest) / (n + str);
      if (adjust && v > cmin[o]) { v = frac * cmin[o]; n_adjusted++; }
      r[o]     = v;
      imputed += v;
    }

    double k = 1 - imputed;
    kept[s]  = k;

    for (int o = 0; o < D; o++) {
      double v = (x[o] == 0) ? r[o] : x[o] / n * k;
      if (margin == 1) out[s + (size_t) ns * o] = v;
      else             out[o + (size_t) D * s]  = v;
    }
  }

  SEXP result = PROTECT(allocVector(VECSXP, 3));
  SEXP names  = PROTECT(allocVector(STRSXP, 3));
  SET_VECTOR_ELT(result, 0, sexp_res);
  SET_VECTOR_ELT(result, 1, sexp_kept);
  SET_VECTOR_ELT(result, 2, ScalarReal(n_adjusted));
  SET_STRING_ELT(names, 0, mkChar("result"));
  SET_STRING_ELT(names, 1, mkChar("kept"));
  SET_STRING_ELT(names, 2, mkChar("n_adjusted"));
  setAttrib(result, R_NamesSymbol, names);

  free_all();
  UNPROTECT(4);
  return result;
}
