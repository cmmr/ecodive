// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The SIMPER breakdown follows ecopy's simper() (ecopy 0.1.2.2,
// ecopy/matrix_comp/simper.py), distributed under the MIT License:
//
//   Copyright (c) 2017 Nathan Lemoine
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

/*
 * SIMPER (Clarke 1993): each taxon's contribution to the
 * mean Bray-Curtis dissimilarity between two groups.
 *
 * For samples i (group a) and j (group b), with totals T_i
 * and T_j, taxon k contributes
 *
 *     c_ijk = |x_ik - x_jk| / (T_i + T_j),
 *
 * and sum_k c_ijk is their Bray-Curtis dissimilarity. As in
 * ecopy, each taxon's average is the mean of c_ijk over the
 * n_a n_b between-group pairs, and its sd their standard
 * deviation (n_a n_b - 1 denominator).
 *
 * Sparsity: c_ijk is 0 unless taxon k is in sample i or j.
 * So taxon k's sum runs over the samples that have it,
 *
 *     S_k = sum_{i in a, x_ik > 0} sum_{j in b} |x_ik - x_jk| w_ij
 *         + sum_{j in b, x_jk > 0} sum_{i in a, x_ik = 0} x_jk w_ij,
 *
 * with w_ij = 1 / (T_i + T_j): about (nnz of k) x n / 2
 * operations, and nothing for taxa absent from both groups.
 *
 * Permutation p-values (Clarke 1993): the pair's group
 * labels are shuffled among its samples (permuted sample a
 * takes the label lab[p[a]], as vegan's grouping[p]), and
 * each taxon's S_k recomputed. A taxon's p-value counts the
 * permutations whose S_k is at least the observed one, within
 * a relative 1e-7 (so rounding doesn't hide a tie): Clarke's
 * (1993) rule. vegan::simper() counts only larger ones, which
 * makes p-values too small when groups are small.
 *
 * Determinism: the observed sums are split over taxa, each
 * summed by one thread in a fixed order. Each permutation is
 * computed by one thread, in a fixed order, and only integer
 * counts are combined across threads. So results are
 * bit-identical for any number of threads.
 */

#include "ecodive.h"

#define SP_TOL 1e-7       // Relative tolerance for ties.
#define SP_PB  16         // Permutations per task.


//------------------------------------------------------
// One pair of groups' samples, taxon by taxon (CSC):
// taxon k's samples are row[cp[k] .. cp[k + 1] - 1]
// (positions 0 .. m - 1 in the pair), with values val.
//------------------------------------------------------
typedef struct {
  int           m;      // Samples in the pair.
  int           p;      // Taxa.
  const int    *cp;     // p + 1 column pointers.
  const int    *row;    // Sample positions.
  const double *val;    // Values (all > 0).
  const double *tot;    // Each sample's total (m).
} sp_pair_t;


//------------------------------------------------------
// S_k for taxa k0 .. k1 - 1, given each sample's label L
// (0: group a, 1: group b), and A / B: the samples of each
// group, in increasing order. xs is m zeros on entry and
// on exit. If ss is not NULL, also the sum of squared
// deviations from mean[k] over all n_a n_b pairs, given
// the observed averages `mean`.
//------------------------------------------------------
static void sp_sums(
    const sp_pair_t *P, int k0, int k1, const int *L,
    const int *A, int na, const int *B, int nb, double *xs,
    double *S, const double *mean, double *ss ) {

  const double *tot = P->tot;

  for (int k = k0; k < k1; k++) {

    int c0 = P->cp[k], c1 = P->cp[k + 1];
    if (c0 == c1) { S[k] = 0; if (ss) ss[k] = 0; continue; }

    for (int e = c0; e < c1; e++) xs[P->row[e]] = P->val[e];

    double s = 0;
    for (int e = c0; e < c1; e++) {
      int    r = P->row[e];
      double v = P->val[e];
      if (L[r] == 0) {
        for (int t = 0; t < nb; t++) {
          int j = B[t];
          s += fabs(v - xs[j]) / (tot[r] + tot[j]);
        }
      }
      else {
        for (int t = 0; t < na; t++) {
          int i = A[t];
          if (xs[i] == 0) s += v / (tot[i] + tot[r]);
        }
      }
    }
    S[k] = s;

    if (ss) {
      // Pairs with taxon k in either sample, and their
      // squared deviations; the other pairs have c = 0.
      double mu = mean[k], q = 0, present = 0;
      for (int e = c0; e < c1; e++) {
        int    r = P->row[e];
        double v = P->val[e];
        if (L[r] == 0) {
          present += nb;
          for (int t = 0; t < nb; t++) {
            int    j = B[t];
            double d = fabs(v - xs[j]) / (tot[r] + tot[j]) - mu;
            q += d * d;
          }
        }
        else {
          for (int t = 0; t < na; t++) {
            int i = A[t];
            if (xs[i] == 0) {
              double d = v / (tot[i] + tot[r]) - mu;
              q += d * d;
              present += 1;
            }
          }
        }
      }
      ss[k] = q + ((double) na * nb - present) * mu * mu;
    }

    for (int e = c0; e < c1; e++) xs[P->row[e]] = 0;
  }
}


//------------------------------------------------------
// Shared with the worker threads.
//------------------------------------------------------
static const sp_pair_t *sp_P;
static const int       *sp_lab;     // Observed labels (m).
static const int       *sp_A, *sp_B;
static int              sp_na, sp_nb;
static double          *sp_S;       // Observed sums (p).
static double          *sp_mean;    // Observed averages (p).
static double          *sp_ss;      // Squared deviations (p).
static int              sp_ntask;   // Tasks this round.
static int              sp_task0;   // First task of the round.
static int              sp_ktask;   // Taxa per task (observed).
static perm_source_t   *sp_ps;
static int              sp_total;   // Permutations + 1.
static double          *sp_work;    // Per-thread doubles.
static int             *sp_iwork;   // Per-thread ints.
static int             *sp_cnt;     // Per-thread counts (p each).


// Observed sums (mode 0) or squared deviations (mode 1),
// over blocks of sp_ktask taxa.
static int sp_mode;

static void *sp_obs_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int m = sp_P->m, p = sp_P->p;

  double *xs = sp_work + (size_t) thread_i * (m + p);
  double *S  = xs + m;   // Scratch: S[k] for the sd pass.

  for (int t = thread_i; t < sp_ntask; t += n_threads) {
    int k0 = t * sp_ktask, k1 = k0 + sp_ktask;
    if (k1 > p) k1 = p;
    if (sp_mode == 0)
      sp_sums(sp_P, k0, k1, sp_lab, sp_A, sp_na, sp_B, sp_nb, xs, sp_S, NULL, NULL);
    else
      sp_sums(sp_P, k0, k1, sp_lab, sp_A, sp_na, sp_B, sp_nb, xs, S, sp_mean, sp_ss);
  }

  return NULL;
}


// Permutations 1 .. sp_total - 1, in tasks of SP_PB.
static void *sp_perm_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int m = sp_P->m, p = sp_P->p;

  double *xs  = sp_work  + (size_t) thread_i * (m + p);
  double *S   = xs + m;
  int    *pp  = sp_iwork + (size_t) thread_i * 4 * m;
  int    *L   = pp + m;
  int    *A   = L + m;
  int    *B   = A + m;
  int    *cnt = sp_cnt + (size_t) thread_i * p;

  for (int t = sp_task0 + thread_i; t < sp_task0 + sp_ntask; t += n_threads) {

    int b0 = 1 + t * SP_PB, b1 = b0 + SP_PB;
    if (b1 > sp_total) b1 = sp_total;

    for (int b = b0; b < b1; b++) {

      perm_source_get(sp_ps, b, pp);

      int na = 0, nb = 0;
      for (int a = 0; a < m; a++) {
        L[a] = sp_lab[pp[a]];
        if (L[a] == 0) A[na++] = a; else B[nb++] = a;
      }

      sp_sums(sp_P, 0, p, L, A, na, B, nb, xs, S, NULL, NULL);

      for (int k = 0; k < p; k++)
        if (S[k] >= sp_S[k] - SP_TOL * sp_S[k]) cnt[k]++;
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// counts, margin: the table (any type ecomatrix reads).
// norm:   normalization code (none, percent, chord,
//         hellinger or binary), applied first.
// grp:    each sample's group, 0 .. g - 1 (integer).
// n_grp:  g.
// pairs:  a list with one element per pair of groups, each
//         list(idx, lab, perms, strata, seed, n_perm):
//           idx:    the pair's samples, 0-based, increasing;
//           lab:    each one's label, 0 (group a) or 1 (b);
//           perms:  m x n_perm 0-based permutations, or NULL;
//           strata: 0-based stratum codes (m), or NULL;
//           seed, n_perm: integers.
// n_threads: integer.
//
// Returns list(means, totals, pairs): means is the p x g
// matrix of group means, totals each sample's total, and
// pairs a list of list(sum, ss, count) (empty if any
// sample's total is 0): each taxon's sum of
// contributions over the pair's n_a n_b sample pairs, sum
// of squared deviations from their mean, and number of
// permutations exceeding the observed sum.
//======================================================
SEXP C_simper(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_grp,
    SEXP sexp_n_grp, SEXP sexp_pairs, SEXP sexp_n_threads ) {

  int norm      = asInteger(sexp_norm);
  int g         = asInteger(sexp_n_grp);
  int n_threads = asInteger(sexp_n_threads);

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n = em->n_samples, p = em->n_otus, nnz = em->nnz;

  for (int e = 0; e < nnz; e++) {
    double v = em->val_vec[e];
    if (ISNAN(v) || v < 0 || !R_FINITE(v)) {
      free_all();
      error("`counts` must not contain NAs, negative or infinite values.");
    }
  }

  if (norm) normalize(em, norm, n_threads, 0);

  if (LENGTH(sexp_grp) != n || g < 1) {
    free_all();                                  // # nocov
    error("Invalid arguments to C_simper().");   // # nocov
  }
  const int *grp = INTEGER(sexp_grp);


  //--------------------------------------------------
  // Copy the table (R_alloc'd), its totals, and the
  // group means, then release ecomatrix's memory.
  //--------------------------------------------------
  int    *pos = (int *)    R_alloc((size_t) n + 1, sizeof(int));
  int    *otu = (int *)    R_alloc(nnz > 0 ? nnz : 1, sizeof(int));
  double *val = (double *) R_alloc(nnz > 0 ? nnz : 1, sizeof(double));
  memcpy(pos, em->pos_vec, ((size_t) n + 1) * sizeof(int));
  if (nnz) {
    memcpy(otu, em->otu_vec, (size_t) nnz * sizeof(int));
    memcpy(val, em->val_vec, (size_t) nnz * sizeof(double));
  }
  free_all();

  SEXP sexp_tot   = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_means = PROTECT(allocMatrix(REALSXP, p, g));
  double *tot   = REAL(sexp_tot);
  double *means = REAL(sexp_means);
  int    *gn    = (int *) R_alloc(g, sizeof(int));

  memset(means, 0, (size_t) p * g * sizeof(double));
  for (int c = 0; c < g; c++) gn[c] = 0;

  for (int s = 0; s < n; s++) {
    if (grp[s] < 0 || grp[s] >= g) error("Invalid arguments to C_simper()."); // # nocov
    double t = 0;
    for (int e = pos[s]; e < pos[s + 1]; e++) {
      t += val[e];
      means[(size_t) grp[s] * p + otu[e]] += val[e];
    }
    tot[s] = t;
    gn[grp[s]]++;
  }
  for (int c = 0; c < g; c++)
    for (int k = 0; k < p; k++)
      means[(size_t) c * p + k] /= gn[c];


  //--------------------------------------------------
  // Each pair of groups.
  //--------------------------------------------------
  // An empty sample makes Bray-Curtis 0 / 0; R reports it.
  int  n_pairs    = LENGTH(sexp_pairs);
  for (int s = 0; s < n; s++) if (tot[s] == 0) n_pairs = 0;
  SEXP sexp_pres  = PROTECT(allocVector(VECSXP, n_pairs));

  for (int q = 0; q < n_pairs; q++) {

    SEXP pr   = VECTOR_ELT(sexp_pairs, q);
    SEXP sidx = VECTOR_ELT(pr, 0);
    int  m    = LENGTH(sidx);
    const int *idx = INTEGER(sidx);
    const int *lab = INTEGER(VECTOR_ELT(pr, 1));
    int  seed   = asInteger(VECTOR_ELT(pr, 4));
    int  n_perm = asInteger(VECTOR_ELT(pr, 5));

    if (LENGTH(VECTOR_ELT(pr, 1)) != m || m < 2 || n_perm < 0)
      error("Invalid arguments to C_simper()."); // # nocov

    // The pair's table, taxon by taxon.
    int    *cp  = (int *)    R_alloc((size_t) p + 1, sizeof(int));
    double *pt  = (double *) R_alloc(m, sizeof(double));
    for (int k = 0; k <= p; k++) cp[k] = 0;
    for (int a = 0; a < m; a++) {
      int s = idx[a];
      if (s < 0 || s >= n) error("Invalid arguments to C_simper()."); // # nocov
      pt[a] = tot[s];
      for (int e = pos[s]; e < pos[s + 1]; e++)
        if (val[e] > 0) cp[otu[e] + 1]++;
    }
    for (int k = 0; k < p; k++) cp[k + 1] += cp[k];
    int     pn   = cp[p];
    int    *row  = (int *)    R_alloc(pn > 0 ? pn : 1, sizeof(int));
    double *pv   = (double *) R_alloc(pn > 0 ? pn : 1, sizeof(double));
    int    *fill = (int *)    R_alloc((size_t) p + 1, sizeof(int));
    memcpy(fill, cp, ((size_t) p + 1) * sizeof(int));
    for (int a = 0; a < m; a++) {
      int s = idx[a];
      for (int e = pos[s]; e < pos[s + 1]; e++)
        if (val[e] > 0) {
          int f = fill[otu[e]]++;
          row[f] = a;
          pv[f]  = val[e];
        }
    }

    sp_pair_t P = { m, p, cp, row, pv, pt };

    int *A = (int *) R_alloc(m, sizeof(int));
    int *B = (int *) R_alloc(m, sizeof(int));
    int  na = 0, nb = 0;
    for (int a = 0; a < m; a++) {
      if (lab[a] != 0 && lab[a] != 1) error("Invalid arguments to C_simper()."); // # nocov
      if (lab[a] == 0) A[na++] = a; else B[nb++] = a;
    }

    SEXP sexp_S   = PROTECT(allocVector(REALSXP, p));
    SEXP sexp_ss  = PROTECT(allocVector(REALSXP, p));
    SEXP sexp_cnt = PROTECT(allocVector(INTSXP,  p));

    sp_P    = &P;
    sp_lab  = lab;
    sp_A    = A;  sp_na = na;
    sp_B    = B;  sp_nb = nb;
    sp_S    = REAL(sexp_S);
    sp_ss   = REAL(sexp_ss);
    sp_mean = (double *) R_alloc(p, sizeof(double));


    //------------------------------------------------
    // Observed sums, then squared deviations: tasks of
    // a fixed number of taxa.
    //------------------------------------------------
    double obs_work = (double) pn * m + p;
    sp_ktask = 256;
    sp_ntask = (p + sp_ktask - 1) / sp_ktask;
    int nt   = threads_for(n_threads, obs_work, 1e6);
    if (nt > sp_ntask) nt = sp_ntask;

    sp_work = (double *) R_alloc(((size_t) m + p) * (nt > 0 ? nt : 1), sizeof(double));
    memset(sp_work, 0, ((size_t) m + p) * (nt > 0 ? nt : 1) * sizeof(double));

    sp_mode = 0;
    run_parallel(sp_obs_worker, nt, nt > 1 ? 100 : 1);
    for (int k = 0; k < p; k++) sp_mean[k] = sp_S[k] / ((double) na * nb);
    sp_mode = 1;
    run_parallel(sp_obs_worker, nt, nt > 1 ? 100 : 1);


    //------------------------------------------------
    // Permutations, in rounds of tasks, checking for
    // interrupts between rounds.
    //------------------------------------------------
    int *cnt = INTEGER(sexp_cnt);
    for (int k = 0; k < p; k++) cnt[k] = 0;

    if (n_perm > 0) {

      perm_source_t ps;
      perm_source_init(&ps, m, VECTOR_ELT(pr, 2), VECTOR_ELT(pr, 3), seed);
      sp_ps    = &ps;
      sp_total = n_perm + 1;

      int    n_tasks  = (n_perm + SP_PB - 1) / SP_PB;
      double per_task = obs_work * SP_PB;
      int    pt_n     = threads_for(n_threads, per_task * n_tasks, 1e6);
      if (pt_n > n_tasks) pt_n = n_tasks;

      sp_work  = (double *) R_alloc(((size_t) m + p) * pt_n, sizeof(double));
      sp_iwork = (int *)    R_alloc((size_t) 4 * m * pt_n, sizeof(int));
      sp_cnt   = (int *)    R_alloc((size_t) p * pt_n, sizeof(int));
      memset(sp_work, 0, ((size_t) m + p) * pt_n * sizeof(double));
      memset(sp_cnt,  0, (size_t) p * pt_n * sizeof(int));

      double rounds    = 1e8 * pt_n / (per_task + 1);
      int    per_round = (rounds < n_tasks) ? (int) rounds : n_tasks;
      if (per_round < pt_n) per_round = pt_n;

      for (sp_task0 = 0; sp_task0 < n_tasks; sp_task0 += sp_ntask) {
        sp_ntask = (n_tasks - sp_task0 > per_round) ? per_round : n_tasks - sp_task0;
        int r_nt = (sp_ntask < pt_n) ? sp_ntask : pt_n;
        run_parallel(sp_perm_worker, r_nt, r_nt > 1 ? 100 : 1);
        if (sp_task0 + sp_ntask < n_tasks) R_CheckUserInterrupt();
      }

      for (int t = 0; t < pt_n; t++)
        for (int k = 0; k < p; k++)
          cnt[k] += sp_cnt[(size_t) t * p + k];
    }

    const char *names[] = { "sum", "ss", "count", "" };
    SEXP sexp_r = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(sexp_r, 0, sexp_S);
    SET_VECTOR_ELT(sexp_r, 1, sexp_ss);
    SET_VECTOR_ELT(sexp_r, 2, sexp_cnt);
    SET_VECTOR_ELT(sexp_pres, q, sexp_r);
    UNPROTECT(4);
  }

  const char *names[] = { "means", "totals", "pairs", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_means);
  SET_VECTOR_ELT(sexp_result, 1, sexp_tot);
  SET_VECTOR_ELT(sexp_result, 2, sexp_pres);

  UNPROTECT(4);
  return sexp_result;
}
