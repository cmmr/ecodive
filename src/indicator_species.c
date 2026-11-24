// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Indicator species analysis (Dufrene & Legendre 1997; De Caceres &
 * Legendre 2009; De Caceres, Legendre & Moretti 2010): the association of
 * each feature with each combination of groups, the best combination, and
 * a permutation test of it. Written from the papers and from the outputs
 * of indicspecies::multipatt() (see plans/indicator_species.md); no code
 * from indicspecies or labdsv.
 *
 * For a feature with values a_i, groups k = 1 .. K of n_k samples, and a
 * combination C of groups (n_C samples, c groups), from each group's sum
 * S_k, number of samples with the feature P_k, and mean m_k = S_k / n_k:
 *
 *   IndVal    sqrt(A B),  A = sum_C S_k / sum_k S_k,  B = sum_C P_k / n_C
 *   IndVal.g  sqrt(A B),  A = sum_C m_k / sum_k m_k
 *   r         the Pearson correlation of a with the indicator of C:
 *             (S_C - n_C mu) / sqrt(SS n_C (n - n_C) / n), with mu the
 *             feature's mean and SS = sum_i (a_i - mu)^2
 *   r.g       the same correlation with sample weights 1 / n_k (each group
 *             weighs the same): (sum_C m_k - c M / K) /
 *             sqrt(V c (K - c) / K), where M = sum_k m_k and
 *             V = sum_k Q_k / n_k - K (mu - M / K)^2, Q_k being group k's
 *             sum of (a_i - mu)^2.
 *
 * Only the samples that have the feature are visited: absent samples add
 * nothing to S_k and P_k, and (n_k - P_k) mu^2 to Q_k.
 *
 * The best combination has the largest value (largest |value| with
 * `negative`); values within a relative 1e-12 of the largest count as
 * ties, which go to the first combination in the list (the fewest groups,
 * then the groups' order). The permutation test shuffles the samples'
 * groups (sample a takes group g[p[a]], as multipatt()'s cluster[p]) and
 * takes the largest value over all the combinations considered (the
 * largest |value| with `negative`): a feature's p-value counts the
 * permutations whose largest value is at least the observed best one,
 * compared exactly, as multipatt() does. A permutation that leaves a
 * combination's samples unchanged gives bit-identical group sums (they
 * are summed in the same order), and so the same value.
 *
 * Determinism: features are independent. Permutations are computed in
 * tasks of a fixed size, each by one thread, and only integer counts are
 * combined across threads. So results are identical for any number of
 * threads.
 */

#include "ecodive.h"

#define IS_TOL 1e-12     // Relative tolerance for ties between combinations.
#define IS_PB  16        // Permutations per task.

#define IDX_INDVAL   0
#define IDX_INDVAL_G 1
#define IDX_R        2
#define IDX_R_G      3


//------------------------------------------------------
// The shared, read-only problem.
//------------------------------------------------------
typedef struct {
  int           n, D, K, nc;  // Samples, features, groups, combinations.
  const int    *fpos;         // Feature f's values: fval[fpos[f] ..
  const int    *fsam;         //   fpos[f + 1] - 1], of samples fsam[...].
  const double *fval;
  const double *mu;           // Each feature's mean.
  const double *ss;           // Each feature's sum of (a_i - mu)^2.
  const int    *mask;         // Each combination's groups (bit k: group k).
  const int    *nk;           // Group sizes (K).
  const double *nC;           // Samples in each combination (nc).
  int           index;
  int           negative;
  int           all;          // The mask of all K groups.
} is_prob_t;


//------------------------------------------------------
// One feature's values for every combination, given each
// sample's group L. w is 3K doubles of scratch. Writes
// the values to str (nc), and returns the largest score
// (value, or |value| with `negative`; -Inf if none is
// defined). If best is not NULL, also the first
// combination within the tolerance of it (-1 if none),
// and its A and B (IndVal indices only).
//------------------------------------------------------
static double is_feature(
    const is_prob_t *P, int f, const int *L, double *w, double *str,
    int *best, double *A_best, double *B_best ) {

  int     K   = P->K;
  double *S   = w;
  double *Pk  = w + K;
  double *Q   = w + 2 * K;
  double  mu  = P->mu[f];
  int     idx = P->index;

  for (int k = 0; k < K; k++) { S[k] = 0; Pk[k] = 0; Q[k] = 0; }

  for (int e = P->fpos[f]; e < P->fpos[f + 1]; e++) {
    int    k = L[P->fsam[e]];
    double v = P->fval[e];
    if (v <= 0) continue;   // A stored zero.
    Pk[k] += 1;
    S[k]  += v;
    if (idx == IDX_R_G) Q[k] += (v - mu) * (v - mu);
  }

  // Per-group quantities that every combination uses.
  double tot = 0, M = 0, V = 0;
  if (idx == IDX_INDVAL) {
    for (int k = 0; k < K; k++) tot += S[k];
  }
  else if (idx == IDX_INDVAL_G || idx == IDX_R_G) {
    for (int k = 0; k < K; k++) { S[k] /= P->nk[k]; M += S[k]; }   // S: means
    if (idx == IDX_R_G) {
      for (int k = 0; k < K; k++)
        V += (Q[k] + (P->nk[k] - Pk[k]) * mu * mu) / P->nk[k];
      double d = mu - M / K;
      V -= K * d * d;
    }
  }

  double smax = R_NegInf;

  for (int c = 0; c < P->nc; c++) {

    int    m  = P->mask[c];
    double sC = 0, pC = 0, nC = P->nC[c], val = R_NaN;
    int    cg = 0;
    for (int k = 0; k < K; k++)
      if (m & (1 << k)) { sC += S[k]; pC += Pk[k]; cg++; }

    switch (idx) {
      case IDX_INDVAL:
        if (tot > 0) val = sqrt(sC / tot * (pC / nC));
        break;
      case IDX_INDVAL_G:
        if (M > 0) val = sqrt(sC / M * (pC / nC));
        break;
      case IDX_R:
        if (P->ss[f] > 0 && nC < P->n)
          val = (sC - nC * mu) / sqrt(P->ss[f] * nC * (P->n - nC) / P->n);
        break;
      case IDX_R_G:
        if (V > 0 && cg < K)
          val = (sC - (double) cg / K * M) / sqrt(V * cg * (K - cg) / K);
        break;
    }

    str[c] = val;
    if (!ISNAN(val)) {
      double s = P->negative ? fabs(val) : val;
      if (s > smax) smax = s;
    }
  }

  if (best) {
    *best = -1; *A_best = R_NaN; *B_best = R_NaN;
    if (smax > R_NegInf) {
      double lim = smax - IS_TOL * fabs(smax);
      for (int c = 0; c < P->nc; c++) {
        double val = str[c];
        if (ISNAN(val)) continue;
        if ((P->negative ? fabs(val) : val) >= lim) { *best = c; break; }
      }
      if (idx == IDX_INDVAL || idx == IDX_INDVAL_G) {
        int    m  = P->mask[*best];
        double sC = 0, pC = 0;
        for (int k = 0; k < K; k++)
          if (m & (1 << k)) { sC += S[k]; pC += Pk[k]; }
        *A_best = sC / (idx == IDX_INDVAL ? tot : M);
        *B_best = pC / P->nC[*best];
      }
    }
  }

  return smax;
}


//------------------------------------------------------
// Shared with the worker threads.
//------------------------------------------------------
static const is_prob_t *is_P;
static const int       *is_grp;     // Observed groups (n).
static perm_source_t   *is_ps;
static int              is_total;   // Permutations + 1.
static int              is_task0;   // First task of the round.
static int              is_ntask;   // Tasks this round.
static int              is_ftask;   // Features per task (observed).
static const double    *is_obs;     // Observed best scores (D).
static const int       *is_test;    // Whether to test each feature (D).
static double          *is_work;    // Per-thread doubles.
static int             *is_iwork;   // Per-thread ints.
static int             *is_cnt;     // Per-thread counts (D each).

// Observed results.
static double *is_str;              // D x nc, or NULL.
static double *is_stat, *is_A, *is_B, *is_smax;
static int    *is_best;


// The observed values, in tasks of is_ftask features.
static void *is_obs_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  const is_prob_t *P = is_P;
  int D = P->D, nc = P->nc;

  double *w   = is_work + (size_t) thread_i * (3 * P->K + nc);
  double *str = w + 3 * P->K;

  for (int t = thread_i; t < is_ntask; t += n_threads) {
    int f0 = t * is_ftask, f1 = f0 + is_ftask;
    if (f1 > D) f1 = D;
    for (int f = f0; f < f1; f++) {
      int    b;
      double A, B;
      is_smax[f] = is_feature(P, f, is_grp, w, str, &b, &A, &B);
      is_best[f] = b;
      is_stat[f] = (b >= 0) ? str[b] : R_NaN;
      is_A[f]    = A;
      is_B[f]    = B;
      if (is_str)
        for (int c = 0; c < nc; c++) is_str[f + (size_t) D * c] = str[c];
    }
  }

  return NULL;
}


// Permutations 1 .. is_total - 1, in tasks of IS_PB.
static void *is_perm_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  const is_prob_t *P = is_P;
  int n = P->n, D = P->D;

  double *w   = is_work  + (size_t) thread_i * (3 * P->K + P->nc);
  double *str = w + 3 * P->K;
  int    *pp  = is_iwork + (size_t) thread_i * 2 * n;
  int    *L   = pp + n;
  int    *cnt = is_cnt + (size_t) thread_i * D;

  for (int t = is_task0 + thread_i; t < is_task0 + is_ntask; t += n_threads) {

    int b0 = 1 + t * IS_PB, b1 = b0 + IS_PB;
    if (b1 > is_total) b1 = is_total;

    for (int b = b0; b < b1; b++) {

      perm_source_get(is_ps, b, pp);
      for (int a = 0; a < n; a++) L[a] = is_grp[pp[a]];

      for (int f = 0; f < D; f++) {
        if (!is_test[f]) continue;
        double s   = is_feature(P, f, L, w, str, NULL, NULL, NULL);
        if (s >= is_obs[f]) cnt[f]++;
      }
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// counts, margin: the table (any type ecomatrix reads).
// grp:      each sample's group, 0 .. K - 1 (integer).
// n_grp:    K (at most 30).
// masks:    the combinations considered, as bit masks
//           (integer; bit k: group k).
// index:    0 IndVal, 1 IndVal.g, 2 r, 3 r.g.
// negative: whether negative values count (r, r.g).
// perms:    n x n_perm 0-based permutations, or NULL.
// strata:   0-based stratum codes (n), or NULL.
// seed, n_perm: integers.
// keep_str: whether to return every combination's value.
// n_threads: integer.
//
// Returns list(best, stat, a, b, count, means, str):
// each feature's best combination (1-based position in
// masks, NA if none), its value, A and B (NA for r and
// r.g), the number of permutations whose largest score is
// at least the observed one (NA: not tested, as when the
// best combination is all K groups), the groups' means
// (D x K), and the D x nc values (NULL unless keep_str).
//======================================================
SEXP C_indicator_species(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_grp, SEXP sexp_n_grp,
    SEXP sexp_masks, SEXP sexp_index, SEXP sexp_negative, SEXP sexp_perms,
    SEXP sexp_strata, SEXP sexp_seed, SEXP sexp_n_perm, SEXP sexp_keep_str,
    SEXP sexp_n_threads ) {

  int K         = asInteger(sexp_n_grp);
  int nc        = LENGTH(sexp_masks);
  int index     = asInteger(sexp_index);
  int negative  = asLogical(sexp_negative) == TRUE;
  int n_perm    = asInteger(sexp_n_perm);
  int keep_str  = asLogical(sexp_keep_str) == TRUE;
  int n_threads = asInteger(sexp_n_threads);

  if (K < 1 || K > 30 || nc < 1 || index < 0 || index > 3 || n_perm < 0)
    error("Invalid arguments to C_indicator_species()."); // # nocov

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, 0, 0, n_threads, &tab);
  int n = tab.n, D = tab.D;

  int nnz = tab.fpos[D];
  for (int e = 0; e < nnz; e++) {
    double v = tab.fval[e];
    if (ISNAN(v) || v < 0 || !R_FINITE(v))
      error("`counts` must not contain NAs, negative or infinite values.");
  }

  if (LENGTH(sexp_grp) != n) error("Invalid arguments to C_indicator_species()."); // # nocov
  const int *grp  = INTEGER(sexp_grp);
  const int *mask = INTEGER(sexp_masks);

  int *nk = (int *) R_alloc(K, sizeof(int));
  for (int k = 0; k < K; k++) nk[k] = 0;
  for (int s = 0; s < n; s++) {
    if (grp[s] < 0 || grp[s] >= K) error("Invalid arguments to C_indicator_species()."); // # nocov
    nk[grp[s]]++;
  }
  for (int k = 0; k < K; k++)
    if (nk[k] == 0) error("Invalid arguments to C_indicator_species()."); // # nocov

  double *nC = (double *) R_alloc(nc, sizeof(double));
  for (int c = 0; c < nc; c++) {
    if (mask[c] <= 0 || mask[c] >= (1 << K)) error("Invalid arguments to C_indicator_species()."); // # nocov
    nC[c] = 0;
    for (int k = 0; k < K; k++) if (mask[c] & (1 << k)) nC[c] += nk[k];
  }


  //--------------------------------------------------
  // Each feature's mean, sum of squared deviations, and
  // the groups' means.
  //--------------------------------------------------
  double *mu = (double *) R_alloc((size_t) D + 1, sizeof(double));
  double *ss = (double *) R_alloc((size_t) D + 1, sizeof(double));
  SEXP sexp_means = PROTECT(allocMatrix(REALSXP, D, K));
  double *means = REAL(sexp_means);
  memset(means, 0, (size_t) D * K * sizeof(double));

  for (int f = 0; f < D; f++) {
    double s = 0;
    for (int e = tab.fpos[f]; e < tab.fpos[f + 1]; e++) {
      s += tab.fval[e];
      means[f + (size_t) D * grp[tab.fsam[e]]] += tab.fval[e];
    }
    double m = s / n, q = 0;
    for (int e = tab.fpos[f]; e < tab.fpos[f + 1]; e++)
      q += (tab.fval[e] - m) * (tab.fval[e] - m);
    q += (double) (n - (tab.fpos[f + 1] - tab.fpos[f])) * m * m;
    mu[f] = m;
    ss[f] = q;
  }
  for (int k = 0; k < K; k++)
    for (int f = 0; f < D; f++) means[f + (size_t) D * k] /= nk[k];

  is_prob_t P = {
    n, D, K, nc, tab.fpos, tab.fsam, tab.fval, mu, ss, mask, nk, nC,
    index, negative, (int) ((1u << K) - 1u) };
  is_P   = &P;
  is_grp = grp;


  //--------------------------------------------------
  // Observed values, in tasks of a fixed number of
  // features.
  //--------------------------------------------------
  SEXP sexp_best = PROTECT(allocVector(INTSXP,  D));
  SEXP sexp_stat = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_A    = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_B    = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_cnt  = PROTECT(allocVector(INTSXP,  D));
  SEXP sexp_str  = PROTECT(keep_str ? allocMatrix(REALSXP, D, nc) : R_NilValue);

  is_best = INTEGER(sexp_best);
  is_stat = REAL(sexp_stat);
  is_A    = REAL(sexp_A);
  is_B    = REAL(sexp_B);
  is_str  = keep_str ? REAL(sexp_str) : NULL;
  is_smax = (double *) R_alloc((size_t) D + 1, sizeof(double));

  double per_feat = (double) nc * K + (double) nnz / (D > 0 ? D : 1);
  is_ftask = 64;
  is_ntask = (D + is_ftask - 1) / is_ftask;
  int nt = threads_for(n_threads, per_feat * D, 1e6);
  if (nt > is_ntask) nt = is_ntask;
  if (nt < 1) nt = 1;

  size_t wlen = (size_t) 3 * K + nc;
  is_work = (double *) R_alloc(wlen * nt, sizeof(double));
  if (is_ntask > 0) run_parallel(is_obs_worker, nt, nt > 1 ? 100 : 1);


  //--------------------------------------------------
  // Which features to test: a best combination that is
  // not all the groups.
  //--------------------------------------------------
  int *test = (int *) R_alloc((size_t) D + 1, sizeof(int));
  int *cnt  = INTEGER(sexp_cnt);
  int  n_test = 0;
  for (int f = 0; f < D; f++) {
    test[f] = is_best[f] >= 0 && mask[is_best[f]] != P.all;
    cnt[f]  = test[f] ? 0 : NA_INTEGER;
    n_test += test[f];
    is_best[f] = (is_best[f] >= 0) ? is_best[f] + 1 : NA_INTEGER;
  }


  //--------------------------------------------------
  // Permutations, in rounds of tasks, checking for
  // interrupts between rounds.
  //--------------------------------------------------
  if (n_perm > 0 && n_test > 0) {

    perm_source_t ps;
    perm_source_init(&ps, n, sexp_perms, sexp_strata, asInteger(sexp_seed));
    is_ps    = &ps;
    is_total = n_perm + 1;
    is_obs   = is_smax;
    is_test  = test;

    int    n_tasks  = (n_perm + IS_PB - 1) / IS_PB;
    double per_task = (per_feat * n_test + n) * IS_PB;
    int    pt_n     = threads_for(n_threads, per_task * n_tasks, 1e6);
    if (pt_n > n_tasks) pt_n = n_tasks;
    if (pt_n < 1) pt_n = 1;

    is_work  = (double *) R_alloc(wlen * pt_n, sizeof(double));
    is_iwork = (int *)    R_alloc((size_t) 2 * n * pt_n, sizeof(int));
    is_cnt   = (int *)    R_alloc((size_t) D * pt_n, sizeof(int));
    memset(is_cnt, 0, (size_t) D * pt_n * sizeof(int));

    double rounds    = 1e8 * pt_n / (per_task + 1);
    int    per_round = (rounds < n_tasks) ? (int) rounds : n_tasks;
    if (per_round < pt_n) per_round = pt_n;

    for (is_task0 = 0; is_task0 < n_tasks; is_task0 += is_ntask) {
      is_ntask = (n_tasks - is_task0 > per_round) ? per_round : n_tasks - is_task0;
      int r_nt = (is_ntask < pt_n) ? is_ntask : pt_n;
      run_parallel(is_perm_worker, r_nt, r_nt > 1 ? 100 : 1);
      if (is_task0 + is_ntask < n_tasks) R_CheckUserInterrupt();
    }

    for (int t = 0; t < pt_n; t++)
      for (int f = 0; f < D; f++)
        if (test[f]) cnt[f] += is_cnt[(size_t) t * D + f];
  }

  const char *names[] = { "best", "stat", "a", "b", "count", "means", "str", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_best);
  SET_VECTOR_ELT(sexp_result, 1, sexp_stat);
  SET_VECTOR_ELT(sexp_result, 2, sexp_A);
  SET_VECTOR_ELT(sexp_result, 3, sexp_B);
  SET_VECTOR_ELT(sexp_result, 4, sexp_cnt);
  SET_VECTOR_ELT(sexp_result, 5, sexp_means);
  SET_VECTOR_ELT(sexp_result, 6, sexp_str);

  UNPROTECT(8);
  return sexp_result;
}
