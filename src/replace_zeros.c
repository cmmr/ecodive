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
