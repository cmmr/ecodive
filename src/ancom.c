// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// ANCOM, ancom(), follows ancom() and _log_compare() from scikit-bio 0.7.4
// (skbio/stats/composition/_ancom.py), and its per-feature p-value
// adjustment follows _adjust_pvalues() (skbio/stats/composition/_utils.py),
// under this license:
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
// license above. The tests' statistics follow SciPy's f_oneway() and
// kruskal() (their definitions and degenerate cases, reproduced from their
// documented behaviour and outputs), and the adjustments R's p.adjust().
//
// For every pair of features i < j, the log ratio d_s = L_si - L_sj of each
// sample (L = log(x + pseudocount)) is tested across the groups:
//
//   anova    One-way ANOVA F, from the group means and two-pass sums of
//            squares of d itself (not from a Gram matrix, which cancels
//            badly for nearly proportional pairs). As SciPy: F = Inf when d
//            is constant within every group, NaN when d is the same in
//            every sample.
//   kruskal  Kruskal-Wallis H with average ranks and the tie correction,
//            ties by exact equality, in SciPy's order of operations (so H is
//            the same to the bit). NaN when every d is tied.
//
// The table is kept sparse: features are transposed into lists of their
// non-zero samples' logs, and every zero has the log of the pseudocount.
//
// Phase 1 (threads) stores each pair's statistic in the packed triangle;
// the main thread turns them into p-values with pf() / pchisq() (R's API
// is not thread-safe); phase 2 (threads) adjusts each feature's row of
// p-values, with a 0 for the diagonal in the family as scikit-bio does,
// and counts W. Each task is one row of the triangle (phase 1) or one
// feature (phase 2), dealt to threads in turn, and depends on nothing
// else, so results are identical for any number of threads.

#include "ecodive.h"
#include <Rmath.h>

#define ANCOM_ANOVA   1
#define ANCOM_KRUSKAL 2

#define ADJ_NONE  0
#define ADJ_BONF  1
#define ADJ_HOLM  2
#define ADJ_BH    3
#define ADJ_BY    4

static int           an_test;
static int           an_n;        // samples
static int           an_D;        // features
static int           an_k;        // groups
static const int    *an_grp;      // 0-based group of each sample
static const int    *an_ng;       // samples in each group
static const int    *an_fpos;     // per feature: start of its list (D + 1)
static const int    *an_fsam;     // sample of each listed value
static const double *an_fval;     // its log
static double        an_lz;       // the log of a zero
static double       *an_stat;     // packed triangle, row by row
static double       *an_buf;      // per share: n (d), k (sums)
static void         *an_vbuf;     // per share: n (vg_t)

static int           an_adj;
static double        an_alpha;
static const double *an_pval;     // packed triangle of p-values
static double       *an_padj;     // D x D adjusted p-values, or NULL
static int          *an_W;
static double       *an_rbuf;     // per share: D values
static int          *an_rord;     // per share: D order


// The index in the packed triangle (row by row, i < j) of pair (i, j).
static inline size_t tri_idx (int i, int j, int D) {
  if (i > j) { int t = i; i = j; j = t; }
  return (size_t) i * D - (size_t) i * (i + 1) / 2 + (size_t) (j - i - 1);
}


// d_s = L_si - L_sj for every sample.
static void pair_diff (int fi, int fj, double *d) {

  int n = an_n;
  for (int s = 0; s < n; s++) d[s] = an_lz - an_lz;   // 0, as log(pc) - log(pc)

  int e = an_fpos[fi], e_end = an_fpos[fi + 1];
  int f = an_fpos[fj], f_end = an_fpos[fj + 1];

  // Walk both lists in sample order.
  while (e < e_end || f < f_end) {
    int si = e < e_end ? an_fsam[e] : n;
    int sj = f < f_end ? an_fsam[f] : n;
    if      (si < sj) { d[si] = an_fval[e] - an_lz;         e++; }
    else if (sj < si) { d[sj] = an_lz      - an_fval[f];    f++; }
    else              { d[si] = an_fval[e] - an_fval[f];    e++; f++; }
  }
}


static double anova_f (const double *d, double *gsum) {

  int n = an_n, k = an_k;
  const int *grp = an_grp;

  // SciPy's degenerate cases: every group constant (Inf), or every sample
  // the same (NaN). Exact equality, as SciPy's np.diff(...) == 0.
  int all_same = 1, all_const = 1;
  for (int g = 0; g < k; g++) gsum[g] = R_NaN;      // first value of group g
  for (int s = 0; s < n; s++) {
    if (d[s] != d[0]) all_same = 0;
    int g = grp[s];
    if (ISNAN(gsum[g])) gsum[g] = d[s];
    else if (d[s] != gsum[g]) all_const = 0;
  }
  if (all_same)  return R_NaN;
  if (all_const) return R_PosInf;

  // The grand mean m, then each group's mean of e = d - m, so that a
  // small difference between the group means isn't lost to rounding in
  // two large means; ebar, the mean of e (0 but for rounding), corrects
  // the between-group sum (the corrected two-pass algorithm). Sums in
  // sample order.
  double tot = 0;
  for (int s = 0; s < n; s++) tot += d[s];
  double m = tot / n;

  double etot = 0;
  for (int g = 0; g < k; g++) gsum[g] = 0;
  for (int s = 0; s < n; s++) { double e = d[s] - m; gsum[grp[s]] += e; etot += e; }
  double ebar = etot / n;
  for (int g = 0; g < k; g++) gsum[g] /= an_ng[g];

  double ssb = 0, ssw = 0, r;
  for (int g = 0; g < k; g++) { r = gsum[g] - ebar; ssb += an_ng[g] * (r * r); }
  for (int s = 0; s < n; s++) { r = (d[s] - m) - gsum[grp[s]]; ssw += r * r; }

  return (ssb / (k - 1)) / (ssw / (n - k));
}


// A merge sort of idx[0 .. n - 1] by d, ties by index (deterministic):
// qsort()'s comparator can't be given d without a static, which threads
// would share. tmp has room for n.
static void sort_order (int *idx, int *tmp, int n, const double *d) {
  for (int w = 1; w < n; w *= 2) {
    for (int lo = 0; lo < n; lo += 2 * w) {
      int mid = lo + w < n ? lo + w : n;
      int hi  = lo + 2 * w < n ? lo + 2 * w : n;
      int a = lo, b = mid, o = lo;
      while (a < mid && b < hi) {
        int ia = idx[a], ib = idx[b];
        if (d[ib] < d[ia] || (d[ib] == d[ia] && ib < ia)) tmp[o++] = idx[b++];
        else                                                tmp[o++] = idx[a++];
      }
      while (a < mid) tmp[o++] = idx[a++];
      while (b < hi)  tmp[o++] = idx[b++];
    }
    for (int i = 0; i < n; i++) idx[i] = tmp[i];
  }
}


static double kruskal_h (const double *d, vg_t *vg, double *rsum) {
  for (int s = 0; s < an_n; s++) { vg[s].v = d[s]; vg[s].g = an_grp[s]; }
  return rank_kruskal_h(vg, an_n, an_ng, an_k, rsum);
}


static void *stat_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D = an_D, n = an_n, k = an_k;

  double *d    = an_buf + (size_t) thread_i * (n + k);
  double *gbuf = d + n;
  vg_t   *vg   = (vg_t *) an_vbuf + (size_t) thread_i * n;

  for (int r = thread_i; r < D - 1; r += n_threads) {
    size_t t = tri_idx(r, r + 1, D);
    for (int c = r + 1; c < D; c++, t++) {
      pair_diff(r, c, d);
      an_stat[t] = an_test == ANCOM_KRUSKAL
                   ? kruskal_h(d, vg, gbuf)
                   : anova_f(d, gbuf);
    }
  }
  return NULL;
}


// The values of v[0 .. m - 1] in increasing order of v, ties by position.
static void sort_row (const double *v, int *ord, int *tmp, int m) {
  for (int i = 0; i < m; i++) ord[i] = i;
  sort_order(ord, tmp, m, v);
}


// R's p.adjust() of c(0, the row without the diagonal and NaNs), for each
// feature, then W = the number of other features with adjusted p < alpha.
static void *adjust_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D = an_D;

  double *v   = an_rbuf + (size_t) thread_i * 2 * D;
  double *adj = v + D;
  int    *ord = an_rord + (size_t) thread_i * 3 * D;
  int    *pos = ord + D;          // the other feature of each family value
  int    *tmp = pos + D;

  for (int i = thread_i; i < D; i += n_threads) {

    // The family: the diagonal's 0 first, then the non-NaN p-values.
    int m = 0;
    v[m] = 0; pos[m++] = -1;
    for (int j = 0; j < D; j++) {
      if (j == i) continue;
      double p = an_pval[tri_idx(i, j, D)];
      if (ISNAN(p)) continue;
      v[m] = p; pos[m++] = j;
    }

    double nd = m;
    switch (an_adj) {

      case ADJ_NONE:
        for (int a = 0; a < m; a++) adj[a] = v[a];
        break;

      case ADJ_BONF:
        for (int a = 0; a < m; a++) { double q = nd * v[a]; adj[a] = q < 1 ? q : 1; }
        break;

      case ADJ_HOLM: {
        // pmin(1, cummax((n + 1 - i) * p[o]))
        sort_row(v, ord, tmp, m);
        double cm = R_NegInf;
        for (int a = 0; a < m; a++) {
          double q = (double) (m - a) * v[ord[a]];
          if (q > cm) cm = q;
          adj[ord[a]] = cm < 1 ? cm : 1;
        }
        break;
      }

      case ADJ_BH:
      case ADJ_BY: {
        // pmin(1, cummin(q * n / i * p[o])), o decreasing, i = n:1
        double qn = nd;
        if (an_adj == ADJ_BY) {
          long double h = 0;   // as R's sum()
          for (int a = 1; a <= m; a++) h += 1 / (double) a;
          qn = (double) h * nd;
        }
        sort_row(v, ord, tmp, m);
        double cm = R_PosInf;
        for (int a = m - 1; a >= 0; a--) {
          double q = qn / (double) (a + 1) * v[ord[a]];
          if (q < cm) cm = q;
          adj[ord[a]] = cm < 1 ? cm : 1;
        }
        break;
      }
    }

    int W = 0;
    for (int a = 1; a < m; a++) if (adj[a] < an_alpha) W++;
    an_W[i] = W;

    if (an_padj) {
      double *col = an_padj + (size_t) i * D;    // row i of a symmetric matrix
      for (int j = 0; j < D; j++) col[j] = R_NaN;
      col[i] = 1;
      for (int a = 1; a < m; a++) col[pos[a]] = adj[a];
    }
  }
  return NULL;
}


//======================================================
// R interface for ancom().
//
// sexp_groups:   0-based group of each sample.
// sexp_n_groups: the number of groups (each non-empty).
// sexp_test:     1 anova, 2 kruskal.
// sexp_adjust:   0 none, 1 bonferroni, 2 holm, 3 BH, 4 BY.
// sexp_details:  TRUE also returns the statistics and the
//                raw and adjusted p-value matrices (tests).
//======================================================
SEXP C_ancom(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_groups, SEXP sexp_n_groups,
    SEXP sexp_test, SEXP sexp_adjust, SEXP sexp_alpha, SEXP sexp_pseudocount,
    SEXP sexp_details, SEXP sexp_n_threads ) {

  int    k       = asInteger(sexp_n_groups);
  int    test    = asInteger(sexp_test);
  int    adjust  = asInteger(sexp_adjust);
  double alpha   = asReal(sexp_alpha);
  double pc      = asReal(sexp_pseudocount);
  int    details = asLogical(sexp_details);
  int    n_thr   = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n   = em->n_samples;
  int D   = em->n_otus;
  int nnz = em->nnz;

  if (pc == 0 && (double) nnz < (double) n * D) {
    free_all();
    errorcall(R_NilValue, "`counts` has zeros, and `pseudocount = 0` leaves their logs undefined.");
  }

  // Transpose into per-feature lists of logs, samples increasing.
  int    *fpos = (int *)    R_alloc((size_t) D + 1, sizeof(int));
  int    *fsam = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  double *fval = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  int    *fill = (int *)    R_alloc((size_t) D, sizeof(int));

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

  // Group sizes.
  const int *grp = INTEGER(sexp_groups);
  int *ng = (int *) R_alloc((size_t) k, sizeof(int));
  for (int g = 0; g < k; g++) ng[g] = 0;
  for (int s = 0; s < n; s++) ng[grp[s]]++;

  size_t n_pairs = (size_t) D * (D - 1) / 2;
  double *stat = (double *) R_alloc(n_pairs + 1, sizeof(double));

  // Phase 1: the statistics. Threads only for enough work.
  double work = (double) n_pairs * n * (test == ANCOM_KRUSKAL ? 4 : 1);
  int    nt1  = threads_for(n_thr, work, 1e6);

  an_test = test;
  an_n    = n;
  an_D    = D;
  an_k    = k;
  an_grp  = grp;
  an_ng   = ng;
  an_fpos = fpos;
  an_fsam = fsam;
  an_fval = fval;
  an_lz   = pc > 0 ? log(pc) : 0;
  an_stat = stat;
  an_buf  = (double *) R_alloc((size_t) nt1 * (n + k), sizeof(double));
  an_vbuf = R_alloc((size_t) nt1 * n, sizeof(vg_t));

  run_parallel(stat_worker, nt1, D - 1 < 100 ? 100 : D - 1);

  // Keep the statistics for `details`, then turn them into p-values in
  // place, on this thread.
  SEXP sexp_stat = R_NilValue;
  if (details) {
    sexp_stat = allocVector(REALSXP, (R_xlen_t) n_pairs);
    for (size_t t = 0; t < n_pairs; t++) REAL(sexp_stat)[t] = stat[t];
  }
  PROTECT(sexp_stat);

  double *pval = stat;
  for (size_t t = 0; t < n_pairs; t++) {
    double s = stat[t];
    if      (ISNAN(s))                pval[t] = R_NaN;
    else if (test == ANCOM_KRUSKAL)   pval[t] = pchisq(s, k - 1, 0, 0);
    else                              pval[t] = pf(s, k - 1, n - k, 0, 0);
  }

  // Phase 2: adjust each feature's row, and count W.
  SEXP sexp_W    = PROTECT(allocVector(INTSXP, D));
  SEXP sexp_padj = R_NilValue;
  if (details) sexp_padj = allocMatrix(REALSXP, D, D);
  PROTECT(sexp_padj);

  int nt2 = threads_for(n_thr, (double) D * D * 4, 1e6);

  an_adj   = adjust;
  an_alpha = alpha;
  an_pval  = pval;
  an_padj  = details ? REAL(sexp_padj) : NULL;
  an_W     = INTEGER(sexp_W);
  an_rbuf  = (double *) R_alloc((size_t) nt2 * 2 * D, sizeof(double));
  an_rord  = (int *)    R_alloc((size_t) nt2 * 3 * D, sizeof(int));

  run_parallel(adjust_worker, nt2, D < 100 ? 100 : D);

  SEXP result;
  if (!details) {
    result = sexp_W;
  }
  else {
    SEXP sexp_pval = PROTECT(allocVector(REALSXP, (R_xlen_t) n_pairs));
    for (size_t t = 0; t < n_pairs; t++) REAL(sexp_pval)[t] = pval[t];
    result = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, sexp_W);
    SET_VECTOR_ELT(result, 1, sexp_stat);
    SET_VECTOR_ELT(result, 2, sexp_pval);
    SET_VECTOR_ELT(result, 3, sexp_padj);
    SEXP nms = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, mkChar("W"));
    SET_STRING_ELT(nms, 1, mkChar("statistic"));
    SET_STRING_ELT(nms, 2, mkChar("p_value"));
    SET_STRING_ELT(nms, 3, mkChar("p_adjusted"));
    setAttrib(result, R_NamesSymbol, nms);
    UNPROTECT(3);
  }

  free_all();
  UNPROTECT(3);
  return result;
}
