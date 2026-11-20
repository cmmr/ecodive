// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Partitioning around medoids (PAM), silhouettes and the
 * Calinski-Harabasz index, on a packed lower triangle of
 * distances, for enterotype().
 *
 * Written from the published algorithms (clean room):
 * Kaufman & Rousseeuw (1990), Finding Groups in Data, ch. 2;
 * Rousseeuw (1987), Silhouettes; Calinski & Harabasz (1974).
 * cluster::pam() and clusterSim::index.G1() are GPL, and were
 * only run, to compare outputs.
 *
 * PAM
 * ---
 * BUILD: the first medoid minimises sum_j d_ij. Each later
 * one maximises sum_j max(D_j - d_ij, 0), with D_j the
 * distance from j to its nearest medoid so far, summed over
 * every non-medoid j, the candidate i included (it gains
 * D_i, since d_ii = 0): the reduction in the total distance
 * to the nearest medoid. (cluster::pam()'s outputs settle
 * that i is included.)
 *
 * SWAP: for every medoid i and non-medoid h, the change in
 * the total distance to the nearest medoid when h replaces
 * i is T_ih = sum_j C_jih, where, with D_j and E_j the
 * distances from j to its nearest and second nearest
 * medoids,
 *
 *     C_jih = min(d_jh, E_j) - D_j       if i is j's nearest,
 *     C_jih = min(d_jh - D_j, 0)         otherwise.
 *
 * (When j is as close to another medoid as to i, E_j = D_j
 * and the two forms agree, so ties in "nearest" don't
 * matter.) The best swap is made if T_ih < 0; repeat until
 * none improves.
 *
 * Ties, as cluster::pam()'s outputs show them: in BUILD, the
 * highest candidate index; in SWAP, the lowest h, then the
 * lowest medoid i; and a sample equally near two medoids
 * goes to the lower-indexed one. (cluster::pam() breaks a
 * few ties in its first BUILD step otherwise, probably by
 * rounding in its own arithmetic; exact ties aren't
 * required to match.)
 *
 * Determinism: every T_ih, and every BUILD gain, is a sum
 * over j = 0 .. n - 1 in order, computed by one thread.
 * Threads take candidates h (or i) by index, and the best is
 * chosen afterwards in index order. So results are
 * bit-identical for any number of threads.
 */

#include "ecodive.h"
#include <float.h>


//------------------------------------------------------
// d_ij from the packed lower triangle (i != j).
//------------------------------------------------------
static inline double pam_d(const double *d, int i, int j, int n) {
  if (i == j) return 0;
  if (i < j) { int t = i; i = j; j = t; }
  return d[dist_col_start(j, n) + (size_t) (i - j - 1)];
}


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// each candidate's own slot of the output.
//------------------------------------------------------
static int           pm_n;       // Samples.
static int           pm_k;       // Medoids.
static const double *pm_d;       // Packed lower triangle.
static int           pm_step;    // BUILD step (0 = first medoid).
static const int    *pm_ismed;   // 1 for medoids.
static const int    *pm_med;     // Medoids, by index (k).
static const int    *pm_near;    // Slot of each sample's nearest medoid.
static const double *pm_dnear;   // D_j.
static const double *pm_dsec;    // E_j.
static double       *pm_score;   // Per candidate: gain, or best T.
static int          *pm_slot;    // Per candidate: best medoid slot.
static double       *pm_work;    // Per-thread scratch (k each).


//------------------------------------------------------
// BUILD: the score of each candidate i.
//   step 0: sum_j d_ij (smaller is better).
//   later:  sum_j max(D_j - d_ij, 0) over non-medoids j,
//           i included (larger is better).
//------------------------------------------------------
static void *pam_build_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pm_n;

  for (int i = thread_i; i < n; i += n_threads) {

    if (pm_ismed[i]) continue;
    double s = 0;

    if (pm_step == 0) {
      for (int j = 0; j < n; j++)
        s += pam_d(pm_d, i, j, n);
    }
    else {
      for (int j = 0; j < n; j++) {
        if (pm_ismed[j]) continue;
        double g = pm_dnear[j] - pam_d(pm_d, i, j, n);
        if (g > 0) s += g;
      }
    }

    pm_score[i] = s;
  }

  return NULL;
}


//------------------------------------------------------
// SWAP: for each non-medoid h, T_ih for every medoid i,
// and the best of them (the lowest slot on ties).
//------------------------------------------------------
static void *pam_swap_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pm_n;
  int k         = pm_k;
  double *T     = pm_work + (size_t) thread_i * k;

  for (int h = thread_i; h < n; h += n_threads) {

    if (pm_ismed[h]) continue;

    for (int m = 0; m < k; m++) T[m] = 0;

    for (int j = 0; j < n; j++) {
      double djh = pam_d(pm_d, j, h, n);
      double Dj  = pm_dnear[j];
      double c   = djh - Dj;
      if (c > 0) c = 0;
      int    nj  = pm_near[j];
      for (int m = 0; m < k; m++) {
        if (m == nj) {
          double Ej = pm_dsec[j];
          T[m] += (djh < Ej ? djh : Ej) - Dj;
        }
        else {
          T[m] += c;
        }
      }
    }

    int    best = 0;
    for (int m = 1; m < k; m++)
      if (T[m] < T[best]) best = m;

    pm_score[h] = T[best];
    pm_slot[h]  = best;
  }

  return NULL;
}


//------------------------------------------------------
// Nearest and second nearest medoid of every sample.
// med[] is sorted, so ties go to the lowest index (as in
// cluster::pam()'s outputs); a medoid is always its own
// nearest.
// Returns the total distance to the nearest medoid.
//------------------------------------------------------
static double pam_assign(
    const double *d, int n, int k, const int *med,
    int *near, double *dnear, double *dsec ) {

  double total = 0;

  for (int j = 0; j < n; j++) {
    double best = R_PosInf, second = R_PosInf;
    int    nb   = 0;
    for (int m = 0; m < k; m++) {
      // A medoid is nearest to itself, even if another medoid is a
      // duplicate of it (at distance 0).
      double v = (med[m] == j) ? -1 : pam_d(d, j, med[m], n);
      if (v < best) { second = best; best = v; nb = m; }
      else if (v < second) { second = v; }
    }
    if (best < 0) best = 0;
    near[j]  = nb;
    dnear[j] = best;
    dsec[j]  = (k > 1) ? second : best;
    total   += best;
  }

  return total;
}


static void sort_int(int *x, int len) {
  for (int a = 1; a < len; a++) {
    int v = x[a], b = a - 1;
    while (b >= 0 && x[b] > v) { x[b + 1] = x[b]; b--; }
    x[b + 1] = v;
  }
}


//======================================================
// R interface: PAM.
//
// dm:   packed lower triangle (double), checked in R.
// n:    number of samples.
// k:    number of medoids, 1 .. n.
// cpus: threads.
//
// Returns list(medoids, clustering, objective, n_swaps):
// medoids (1-based, ascending); each sample's cluster,
// 1 .. k, numbered in the order of the medoids; the mean
// distance to the nearest medoid after BUILD and after
// SWAP; and the number of swaps made.
//======================================================
SEXP C_pam(SEXP sexp_dm, SEXP sexp_n, SEXP sexp_k, SEXP sexp_cpus) {

  int n         = asInteger(sexp_n);
  int k         = asInteger(sexp_k);
  int n_threads = asInteger(sexp_cpus);

  if (n < 1 || k < 1 || k > n || n_threads < 1 ||
      XLENGTH(sexp_dm) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_pam()."); // # nocov

  const double *d = REAL(sexp_dm);

  int    *ismed = (int *)    R_alloc(n, sizeof(int));
  int    *med   = (int *)    R_alloc(k, sizeof(int));
  int    *near  = (int *)    R_alloc(n, sizeof(int));
  double *dnear = (double *) R_alloc(n, sizeof(double));
  double *dsec  = (double *) R_alloc(n, sizeof(double));
  double *score = (double *) R_alloc(n, sizeof(double));
  int    *slot  = (int *)    R_alloc(n, sizeof(int));
  double *work  = (double *) R_alloc((size_t) n_threads * k, sizeof(double));

  for (int j = 0; j < n; j++) { ismed[j] = 0; dnear[j] = R_PosInf; }

  pm_n     = n;
  pm_k     = k;
  pm_d     = d;
  pm_ismed = ismed;
  pm_med   = med;
  pm_near  = near;
  pm_dnear = dnear;
  pm_dsec  = dsec;
  pm_score = score;
  pm_slot  = slot;
  pm_work  = work;

  // Threads are worth it once each task is O(n) or more.
  int n_tasks = (n >= 100) ? n : 1;


  //--------------------------------------------------
  // BUILD.
  //--------------------------------------------------
  for (int step = 0; step < k; step++) {

    pm_step = step;
    run_parallel(pam_build_worker, n_threads, n_tasks);

    int best = -1;
    for (int i = 0; i < n; i++) {
      if (ismed[i]) continue;
      if (best < 0) { best = i; continue; }
      // Ties go to the highest index, as in cluster::pam()'s
      // outputs.
      if (step == 0 ? (score[i] <= score[best]) : (score[i] >= score[best]))
        best = i;
    }

    ismed[best] = 1;
    med[step]   = best;

    for (int j = 0; j < n; j++) {
      double v = pam_d(d, j, best, n);
      if (v < dnear[j]) dnear[j] = v;
    }
  }

  sort_int(med, k);
  double build = pam_assign(d, n, k, med, near, dnear, dsec);


  //--------------------------------------------------
  // SWAP. A swap must lower the total by more than
  // rounding error, so that rounding can't make it cycle.
  //--------------------------------------------------
  double total   = build;
  int    n_swaps = 0;

  if (k < n) {
    for (;;) {

      run_parallel(pam_swap_worker, n_threads, n_tasks);

      int bh = -1;
      for (int h = 0; h < n; h++) {
        if (ismed[h]) continue;
        if (bh < 0 || score[h] < score[bh]) bh = h;
      }

      if (!(score[bh] < -16 * DBL_EPSILON * total)) break;

      int out = med[slot[bh]];
      ismed[out] = 0;
      ismed[bh]  = 1;
      med[slot[bh]] = bh;
      sort_int(med, k);
      total = pam_assign(d, n, k, med, near, dnear, dsec);
      n_swaps++;

      if (n_swaps > 100 * n) break; // # nocov
    }
  }


  //--------------------------------------------------
  // Results.
  //--------------------------------------------------
  SEXP sexp_med = PROTECT(allocVector(INTSXP,  k));
  SEXP sexp_cl  = PROTECT(allocVector(INTSXP,  n));
  SEXP sexp_obj = PROTECT(allocVector(REALSXP, 2));
  SEXP sexp_res = PROTECT(allocVector(VECSXP,  4));

  for (int m = 0; m < k; m++) INTEGER(sexp_med)[m] = med[m] + 1;
  for (int j = 0; j < n; j++) INTEGER(sexp_cl)[j]  = near[j] + 1;
  REAL(sexp_obj)[0] = build / n;
  REAL(sexp_obj)[1] = total / n;

  SET_VECTOR_ELT(sexp_res, 0, sexp_med);
  SET_VECTOR_ELT(sexp_res, 1, sexp_cl);
  SET_VECTOR_ELT(sexp_res, 2, sexp_obj);
  SET_VECTOR_ELT(sexp_res, 3, ScalarInteger(n_swaps));

  UNPROTECT(4);
  return sexp_res;
}



//======================================================
// Silhouettes (Rousseeuw 1987).
//
// s_i = (b_i - a_i) / max(a_i, b_i), where a_i is the mean
// distance from i to the other members of its cluster, and
// b_i the smallest mean distance from i to another cluster.
// 0 for a sample alone in its cluster, and when a_i = b_i.
//======================================================
static const int *sl_cl;      // Each sample's cluster, 0 .. k - 1.
static const int *sl_size;    // Cluster sizes.
static double    *sl_out;     // Silhouettes.
static int       *sl_nbr;     // Each sample's neighbouring cluster.

static void *silhouette_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = pm_n;
  int k         = pm_k;
  double *S     = pm_work + (size_t) thread_i * k;

  for (int i = thread_i; i < n; i += n_threads) {

    for (int c = 0; c < k; c++) S[c] = 0;
    for (int j = 0; j < n; j++)
      if (j != i) S[sl_cl[j]] += pam_d(pm_d, i, j, n);

    int    own = sl_cl[i];
    double b   = R_PosInf;
    int    nb  = -1;
    for (int c = 0; c < k; c++) {
      if (c == own || sl_size[c] == 0) continue;
      double v = S[c] / sl_size[c];
      if (v < b) { b = v; nb = c; }
    }
    sl_nbr[i] = nb + 1;

    if (sl_size[own] < 2) { sl_out[i] = 0; continue; }

    double a  = S[own] / (sl_size[own] - 1);
    double mx = (a > b) ? a : b;
    sl_out[i] = (mx > 0) ? (b - a) / mx : 0;
  }

  return NULL;
}


//======================================================
// R interface: silhouettes.
//
// cl: each sample's cluster, 1 .. k (integer).
// Returns list(silhouette, neighbor).
//======================================================
SEXP C_silhouette(SEXP sexp_dm, SEXP sexp_n, SEXP sexp_cl, SEXP sexp_k, SEXP sexp_cpus) {

  int n         = asInteger(sexp_n);
  int k         = asInteger(sexp_k);
  int n_threads = asInteger(sexp_cpus);

  if (n < 1 || k < 1 || n_threads < 1 || LENGTH(sexp_cl) != n ||
      XLENGTH(sexp_dm) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_silhouette()."); // # nocov

  int *cl   = (int *) R_alloc(n, sizeof(int));
  int *size = (int *) R_alloc(k, sizeof(int));
  for (int c = 0; c < k; c++) size[c] = 0;
  for (int j = 0; j < n; j++) {
    cl[j] = INTEGER(sexp_cl)[j] - 1;
    if (cl[j] < 0 || cl[j] >= k) error("Invalid arguments to C_silhouette()."); // # nocov
    size[cl[j]]++;
  }

  SEXP sexp_s   = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_nb  = PROTECT(allocVector(INTSXP,  n));
  SEXP sexp_res = PROTECT(allocVector(VECSXP,  2));

  pm_n    = n;
  pm_k    = k;
  pm_d    = REAL(sexp_dm);
  pm_work = (double *) R_alloc((size_t) n_threads * k, sizeof(double));
  sl_cl   = cl;
  sl_size = size;
  sl_out  = REAL(sexp_s);
  sl_nbr  = INTEGER(sexp_nb);

  run_parallel(silhouette_worker, n_threads, (n >= 100) ? n : 1);

  SET_VECTOR_ELT(sexp_res, 0, sexp_s);
  SET_VECTOR_ELT(sexp_res, 1, sexp_nb);
  UNPROTECT(3);
  return sexp_res;
}



//======================================================
// R interface: the Calinski-Harabasz index.
//
// cl:   each sample's cluster, 1 .. k (integer).
// form: 0 for the centroid form, 1 for the medoid form.
// x:    for the medoid form, an n x p matrix (double) of the
//       samples' coordinates, or NULL to use the distances.
//
// Centroid form (the pseudo-F of distances; for Euclidean
// distances, exactly Calinski & Harabasz's index):
//
//     T = (1/n) sum_{i<j} d_ij^2
//     W = sum_c (1/n_c) sum_{i<j in c} d_ij^2
//     CH = ((T - W) / (k - 1)) / (W / (n - k))
//
// Medoid form, as clusterSim::index.G1(x, cl, d,
// centrotypes = "medoids") computes it (pinned from its
// outputs): m_c, the medoid of cluster c, and M, that of all
// samples, are the members with the smallest total distance
// d to the others (the lowest index on ties). Then, with
// e(i, j) the Euclidean distance between rows i and j of x
// (or d(i, j) when x is NULL),
//
//     T = sum_i e(i, M)^2
//     W = sum_i e(i, m_c(i))^2
//     CH = ((T - W) / (k - 1)) / (W / (n - k))
//
// Sums run in index order. Returns CH (Inf when W = 0).
//======================================================
static double ch_e2(const double *d, const double *x, int n, int p, int i, int j) {
  if (x == NULL) { double v = pam_d(d, i, j, n); return v * v; }
  double s = 0;
  for (int f = 0; f < p; f++) {
    double v = x[(size_t) f * n + i] - x[(size_t) f * n + j];
    s += v * v;
  }
  return s;
}

SEXP C_calinski(SEXP sexp_dm, SEXP sexp_n, SEXP sexp_cl, SEXP sexp_k, SEXP sexp_form, SEXP sexp_x) {

  int n    = asInteger(sexp_n);
  int k    = asInteger(sexp_k);
  int form = asInteger(sexp_form);

  if (n < 2 || k < 2 || k >= n || LENGTH(sexp_cl) != n ||
      XLENGTH(sexp_dm) != (R_xlen_t) n * (n - 1) / 2 ||
      (!isNull(sexp_x) && (!isReal(sexp_x) || !isMatrix(sexp_x) || nrows(sexp_x) != n)))
    error("Invalid arguments to C_calinski()."); // # nocov

  const double *d  = REAL(sexp_dm);
  const int    *cl = INTEGER(sexp_cl);
  for (int j = 0; j < n; j++)
    if (cl[j] < 1 || cl[j] > k) error("Invalid arguments to C_calinski()."); // # nocov

  double T, W;

  if (form == 0) {

    double *wc   = (double *) R_alloc(k, sizeof(double));
    int    *size = (int *)    R_alloc(k, sizeof(int));
    for (int c = 0; c < k; c++) { wc[c] = 0; size[c] = 0; }
    for (int j = 0; j < n; j++) size[cl[j] - 1]++;

    T = 0;
    for (int j = 0; j < n - 1; j++) {
      const double *dj = d + dist_col_start(j, n);
      int cj = cl[j];
      for (int i = j + 1; i < n; i++) {
        double v = dj[i - j - 1];
        v *= v;
        T += v;
        if (cl[i] == cj) wc[cj - 1] += v;
      }
    }

    T /= n;
    W  = 0;
    for (int c = 0; c < k; c++)
      if (size[c] > 0) W += wc[c] / size[c];
  }

  else {

    const double *x = isNull(sexp_x) ? NULL : REAL(sexp_x);
    int           p = isNull(sexp_x) ? 0    : ncols(sexp_x);

    // Each sample's total distance to its cluster, and to all.
    double *sc  = (double *) R_alloc(n, sizeof(double));
    double *sa  = (double *) R_alloc(n, sizeof(double));
    int    *med = (int *)    R_alloc(k, sizeof(int));
    for (int j = 0; j < n; j++) { sc[j] = 0; sa[j] = 0; }
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        double v = pam_d(d, i, j, n);
        sa[j] += v;
        if (cl[i] == cl[j]) sc[j] += v;
      }
    }

    int M = 0;
    for (int j = 1; j < n; j++) if (sa[j] < sa[M]) M = j;
    for (int c = 0; c < k; c++) med[c] = -1;
    for (int j = 0; j < n; j++) {
      int c = cl[j] - 1;
      if (med[c] < 0 || sc[j] < sc[med[c]]) med[c] = j;
    }

    T = 0;
    W = 0;
    for (int j = 0; j < n; j++) {
      T += ch_e2(d, x, n, p, j, M);
      W += ch_e2(d, x, n, p, j, med[cl[j] - 1]);
    }
  }

  return ScalarReal(((T - W) / (k - 1)) / (W / (n - k)));
}
