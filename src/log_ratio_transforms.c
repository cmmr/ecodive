// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The additive and isometric log-ratio transforms, alr() and ilr(), follow
// alr(), ilr() and _gram_schmidt_basis() from scikit-bio 0.7.4
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
// The rest of this file is Copyright 2026 ecodive authors, under the MIT
// license above.
//
// Each sample is one task: its logs log(x + pseudocount) are expanded into a
// dense scratch vector of the thread's own, then
//
//   alr: log(x_j) - log(x_den) for every j but the denominator, in order;
//   ilr: the clr (logs minus their mean), then multiplied by the basis:
//        * the default (Helmert / Gram-Schmidt) basis by its O(D)
//          recurrence: coordinate j (1-based) is a_j * (c_1 + ... + c_j) +
//          b_j * c_{j+1}, with a_j = (1 / j) * sqrt(j / (j + 1)) and
//          b_j = -sqrt(j / (j + 1)), skbio's entries to the last bit;
//        * any other basis as its non-zero entries by row (CSR), so an SBP
//          basis costs only its non-zeros.
//
// Samples don't share anything, so results are identical for any number of
// threads.

#include "ecodive.h"

static ecomatrix_t  *lr_em;
static int           lr_margin;
static int           lr_n_coords;
static double        lr_pc;
static int           lr_den;      // alr: 0-based denominator; ilr: -1
static const int    *lr_bp;       // ilr: CSR row pointers (n_coords + 1), or NULL = Helmert
static const int    *lr_bi;       //      0-based feature of each entry
static const double *lr_bx;       //      value of each entry
static const double *lr_ha;       // Helmert a_j and b_j (n_coords each)
static const double *lr_hb;
static double       *lr_scratch;  // D doubles per thread
static double       *lr_out;


// Sample s's logs of (x + pseudocount), every feature.
static void sample_logs(int s, double *L) {

  const ecomatrix_t *em = lr_em;
  int    D  = em->n_otus;
  int    e0 = em->pos_vec[s], e1 = em->pos_vec[s + 1];

  // Unstored zeros: log(pseudocount). With pseudocount 0, R has checked
  // that there are none.
  if (e1 - e0 < D) {
    double lz = log(lr_pc);
    for (int o = 0; o < D; o++) L[o] = lz;
  }
  for (int e = e0; e < e1; e++)
    L[em->otu_vec[e]] = log(em->val_vec[e] + lr_pc);
}


static inline void store(int s, int k, double v) {
  if (lr_margin == 1) lr_out[s + (size_t) lr_em->n_samples * k] = v;
  else                lr_out[k + (size_t) lr_n_coords * s]      = v;
}


static void *alr_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D         = lr_em->n_otus;
  double *L     = lr_scratch + (size_t) thread_i * D;

  for (int s = thread_i; s < lr_em->n_samples; s += n_threads) {
    sample_logs(s, L);
    double ld = L[lr_den];
    int k = 0;
    for (int o = 0; o < D; o++)
      if (o != lr_den) store(s, k++, L[o] - ld);
  }
  return NULL;
}


static void *ilr_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D         = lr_em->n_otus;
  double *c     = lr_scratch + (size_t) thread_i * D;

  for (int s = thread_i; s < lr_em->n_samples; s += n_threads) {

    // clr, as normalize(norm = 'clr').
    sample_logs(s, c);
    double logsum = 0;
    for (int o = 0; o < D; o++) logsum += c[o];
    double log_mean = logsum / D;
    for (int o = 0; o < D; o++) c[o] -= log_mean;

    if (lr_bp == NULL) {
      double run = 0;
      for (int j = 0; j < lr_n_coords; j++) {
        run += c[j];
        store(s, j, lr_ha[j] * run + lr_hb[j] * c[j + 1]);
      }
    }
    else {
      for (int k = 0; k < lr_n_coords; k++) {
        double v = 0;
        for (int e = lr_bp[k]; e < lr_bp[k + 1]; e++)
          v += lr_bx[e] * c[lr_bi[e]];
        store(s, k, v);
      }
    }
  }
  return NULL;
}


// The Helmert basis's two distinct entries in row j (0-based), in
// scikit-bio's arithmetic: e = [1/i] * i + [-1] + [0] * ..., times
// sqrt(i / (i + 1)), with i = j + 1.
static void helmert_entries(int n_coords, double *a, double *b) {
  for (int j = 0; j < n_coords; j++) {
    double i = j + 1;
    double r = sqrt(i / (i + 1));
    a[j] = (1 / i) * r;
    b[j] = -1 * r;
  }
}


//======================================================
// R interface: the default ILR basis, (D - 1) x D, as
// scikit-bio's _gram_schmidt_basis(D).
//======================================================
SEXP C_helmert_basis(SEXP sexp_D) {

  int D  = asInteger(sexp_D);
  int nc = D - 1;

  SEXP result = PROTECT(allocMatrix(REALSXP, nc, D));
  double *out = REAL(result);
  memset(out, 0, sizeof(double) * (size_t) nc * D);

  double *a = (double *) R_alloc((size_t) nc + 1, sizeof(double));
  double *b = (double *) R_alloc((size_t) nc + 1, sizeof(double));
  helmert_entries(nc, a, b);

  for (int j = 0; j < nc; j++) {
    for (int o = 0; o <= j; o++) out[j + (size_t) nc * o] = a[j];
    out[j + (size_t) nc * (j + 1)] = b[j];
  }

  UNPROTECT(1);
  return result;
}


//======================================================
// R interface for alr() and ilr(). The coordinates as a
// dense matrix in the input's orientation: samples x
// coordinates for margin 1, coordinates x samples for
// margin 2.
//
// sexp_den:   alr: the 0-based denominator; ilr: -1.
// sexp_basis: ilr: list(p, i, x), the basis's non-zero
//             entries by row (0-based), or NULL for the
//             Helmert basis.
//======================================================
SEXP C_log_ratio_coords(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_den,
    SEXP sexp_basis, SEXP sexp_pseudocount, SEXP sexp_n_threads ) {

  int margin = asInteger(sexp_margin);
  int n_thr  = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int ns = em->n_samples;
  int D  = em->n_otus;
  int nc = D - 1;

  lr_em     = em;
  lr_margin = margin;
  lr_pc     = asReal(sexp_pseudocount);
  lr_den    = asInteger(sexp_den);
  lr_bp     = NULL;
  lr_bi     = NULL;
  lr_bx     = NULL;
  lr_ha     = NULL;
  lr_hb     = NULL;

  if (lr_den < 0 && !isNull(sexp_basis)) {
    lr_bp = INTEGER(VECTOR_ELT(sexp_basis, 0));
    lr_bi = INTEGER(VECTOR_ELT(sexp_basis, 1));
    lr_bx = REAL(VECTOR_ELT(sexp_basis, 2));
    nc    = length(VECTOR_ELT(sexp_basis, 0)) - 1;
  }
  else if (lr_den < 0) {
    double *ab = (double *) R_alloc(2 * (size_t) nc + 2, sizeof(double));
    helmert_entries(nc, ab, ab + nc + 1);
    lr_ha = ab;
    lr_hb = ab + nc + 1;
  }
  lr_n_coords = nc;

  // run_parallel() hands out at most n_thr shares.
  lr_scratch = (double *) R_alloc((size_t) n_thr * D, sizeof(double));

  SEXP result = PROTECT(margin == 1 ? allocMatrix(REALSXP, ns, nc) : allocMatrix(REALSXP, nc, ns));
  lr_out = REAL(result);

  run_parallel(lr_den >= 0 ? alr_worker : ilr_worker, n_thr, ns);

  free_all();
  UNPROTECT(1);
  return result;
}
