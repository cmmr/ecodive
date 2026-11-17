// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * PERMDISP (Anderson 2006): each sample's distance to its
 * group's centre, and the permutation test of the ANOVA on
 * those distances, as vegan::betadisper() and
 * vegan::permutest.betadisper().
 *
 * The samples live in the principal coordinates of D. For
 * non-Euclidean D some of them are imaginary: their
 * eigenvalues are negative. As vegan does, a sample's
 * squared distance to its centre c is
 *
 *     z^2 = |x+ - c+|^2 - |x- - c-|^2,
 *
 * over the real (+) and imaginary (-) axes separately.
 *
 * Centroids. With w = e_i - 1_g / n_g, z_i^2 = w' G w over
 * all axes together, and since w sums to 0,
 *
 *     z_i^2 = (1/n_g) sum_{j in g} d_ij^2
 *             - 1/(2 n_g^2) sum_{j,k in g} d_jk^2,
 *
 * which needs only the distances within groups: no
 * eigendecomposition, and O(sum n_g^2) work.
 *
 * Spatial medians. The median minimizes the sum of
 * (Euclidean) distances to the group's points, separately
 * on the real and the imaginary axes. That needs the
 * coordinates: a full eigendecomposition of Gower's
 * centred matrix G = -1/2 J A J, keeping vegan's axes
 * (|lambda| > max(TOL, TOL * lambda_1), TOL = sqrt(eps)).
 * The median is then found with Weiszfeld's algorithm,
 * with Vardi & Zhang's (2000) step at data points, to full
 * precision. vegan instead stops optim()'s BFGS early, so
 * its medians are only accurate to ~1e-5.
 *
 * Permutations: as vegan, the residuals of the ANOVA (each
 * distance minus its group's mean) are permuted, with the
 * groups fixed, and F is recomputed. Pairwise t statistics
 * permute the distances themselves. Permutation k comes
 * from its own PCG stream, (seed, k), or from a supplied
 * matrix; each is computed by one thread in a fixed order,
 * so results are bit-identical for any number of threads.
 * Each pool share (i of n) owns fixed samples, tasks or
 * permutations, and the scratch space indexed by i.
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON

// Weiszfeld: stop when a step moves the median by less than
// this times the mean distance to it.
#define PD_TOL   1e-14
#define PD_MAXIT 1000


// d(i, j), i != j, from the packed lower triangle.
static inline double pd_dist(const double *d, int i, int j, int n) {
  if (i > j) { int t = i; i = j; j = t; }
  return d[dist_col_start(i, n) + (size_t) (j - i - 1)];
}



//======================================================
// Distances to the group centres.
//======================================================

static int           pd_n;
static const double *pd_d;      // Packed lower triangle of D.
static int           pd_ng;     // Groups.
static const int    *pd_beg;    // Group g is pd_mem[pd_beg[g]] ..
static const int    *pd_mem;    // ... pd_mem[pd_beg[g + 1] - 1], ascending.
static double       *pd_rs;     // Centroids: sum of d^2 within the group.

static const double *pd_vec;    // Medians: n x n eigenvectors (column-major),
static const double *pd_scl;    // ... sqrt(|lambda|) for each kept axis,
static const int    *pd_axes;   // ... the kept axes: positive ones first,
static int           pd_npos;   // ... pd_npos of them,
static int           pd_nneg;   // ... then pd_nneg negative ones.
static double       *pd_dsq;    // n x 2: squared distances on + / - axes.
static int          *pd_iters;  // Iterations, per task (-1: not converged).


//------------------------------------------------------
// Centroids: each sample's sum of squared distances to
// the other members of its group. Sample i is done by
// thread i % n_threads, summing in a fixed order.
//------------------------------------------------------
static void *centroid_worker(void *arg) {

  worker_t *w = (worker_t *) arg;
  int n = pd_n;

  for (int g = 0; g < pd_ng; g++) {
    const int *mem = pd_mem + pd_beg[g];
    int        len = pd_beg[g + 1] - pd_beg[g];
    for (int a = 0; a < len; a++) {
      int i = mem[a];
      if (i % w->n != w->i) continue;
      double s = 0;
      for (int b = 0; b < len; b++) {
        if (b == a) continue;
        double x = pd_dist(pd_d, i, mem[b], n);
        s += x * x;
      }
      pd_rs[i] = s;
    }
  }

  return NULL;
}


//------------------------------------------------------
// The spatial median of the `m` rows of Y (m x r,
// row-major), by Weiszfeld's algorithm with Vardi &
// Zhang's modification. Writes each row's squared
// distance to the median into dsq. Returns iterations,
// or -1 if not converged. ws needs 2 * r + m doubles.
//------------------------------------------------------
static int spatial_median(const double *Y, int m, int r, double *dsq, double *ws) {

  double *med = ws;
  double *nxt = ws + r;
  double *dd  = ws + 2 * r;
  int     its = 0;

  // One axis: the ordinary median, as vegan's start
  // (the midpoint of the middle two, for even m), where
  // its optimizer stops. Weiszfeld could stop anywhere
  // between the middle two.
  if (r == 1) {
    for (int j = 0; j < m; j++) dd[j] = Y[j];
    // Insertion sort; m is small, or this is rare.
    for (int j = 1; j < m; j++) {
      double v = dd[j];
      int    k = j - 1;
      while (k >= 0 && dd[k] > v) { dd[k + 1] = dd[k]; k--; }
      dd[k + 1] = v;
    }
    med[0] = (m % 2) ? dd[m / 2] : (dd[m / 2 - 1] + dd[m / 2]) / 2;
  }
  else {

    // Start at the centroid.
    for (int k = 0; k < r; k++) med[k] = 0;
    for (int j = 0; j < m; j++)
      for (int k = 0; k < r; k++) med[k] += Y[(size_t) j * r + k];
    for (int k = 0; k < r; k++) med[k] /= m;

    while (m > 2) {

      // Distances to the current median.
      double sum_dd = 0, sum_w = 0;
      int    n_zero = 0;
      for (int j = 0; j < m; j++) {
        const double *y = Y + (size_t) j * r;
        double s = 0;
        for (int k = 0; k < r; k++) { double t = y[k] - med[k]; s += t * t; }
        dd[j]   = sqrt(s);
        sum_dd += dd[j];
        if (dd[j] > 0) sum_w += 1 / dd[j];
        else           n_zero++;
      }
      if (sum_dd == 0) break; // All points at the median.

      // Weighted mean of the other points.
      for (int k = 0; k < r; k++) nxt[k] = 0;
      for (int j = 0; j < m; j++) {
        if (dd[j] == 0) continue;
        const double *y = Y + (size_t) j * r;
        double        w = 1 / dd[j];
        for (int k = 0; k < r; k++) nxt[k] += w * y[k];
      }
      for (int k = 0; k < r; k++) nxt[k] /= sum_w;

      // At a data point: move only if the pull of the
      // others outweighs it (Vardi & Zhang 2000).
      if (n_zero > 0) {
        double rn = 0;
        for (int k = 0; k < r; k++) {
          double t = (nxt[k] - med[k]) * sum_w;
          rn += t * t;
        }
        rn = sqrt(rn);
        double gam = (rn > 0) ? n_zero / rn : 1;
        if (gam >= 1) break; // The median is this point.
        for (int k = 0; k < r; k++)
          nxt[k] = (1 - gam) * nxt[k] + gam * med[k];
      }

      double step = 0;
      for (int k = 0; k < r; k++) {
        double t = nxt[k] - med[k];
        step += t * t;
        med[k] = nxt[k];
      }
      its++;

      if (sqrt(step) <= PD_TOL * sum_dd / m) break;
      if (its >= PD_MAXIT) { its = -1; break; }
    }
  }

  for (int j = 0; j < m; j++) {
    const double *y = Y + (size_t) j * r;
    double s = 0;
    for (int k = 0; k < r; k++) { double t = y[k] - med[k]; s += t * t; }
    dsq[j] = s;
  }

  return its;
}


//------------------------------------------------------
// Medians: task t = 2 * g + s (group g, s = 0 for the
// positive axes, 1 for the negative ones) is done by
// thread t % n_threads.
//------------------------------------------------------
static int pd_failed;

// Set by any worker that runs out of memory (all store 1; read after
// they finish). Atomic, so it is not a data race.
#if defined(__GNUC__) || defined(__clang__)
  #define SET_FAILED() __atomic_store_n(&pd_failed, 1, __ATOMIC_RELAXED)
#else
  #define SET_FAILED() (pd_failed = 1)
#endif

static void *median_worker(void *arg) {

  worker_t *w = (worker_t *) arg;
  int n = pd_n;

  for (int t = w->i; t < 2 * pd_ng; t += w->n) {

    int g = t / 2, s = t % 2;
    int r = s ? pd_nneg : pd_npos;
    const int *axes = pd_axes + (s ? pd_npos : 0);
    const int *mem  = pd_mem + pd_beg[g];
    int        m    = pd_beg[g + 1] - pd_beg[g];

    if (r == 0 || m == 0) { pd_iters[t] = 0; continue; }

    double *Y   = (double *) malloc(((size_t) m * r + 2 * (size_t) r + 2 * (size_t) m) * sizeof(double));
    if (!Y) { SET_FAILED(); continue; } // # nocov
    double *dsq = Y + (size_t) m * r;
    double *ws  = dsq + m;

    for (int j = 0; j < m; j++)
      for (int k = 0; k < r; k++)
        Y[(size_t) j * r + k] = pd_vec[(size_t) axes[k] * n + mem[j]] * pd_scl[axes[k]];

    pd_iters[t] = spatial_median(Y, m, r, dsq, ws);

    for (int j = 0; j < m; j++)
      pd_dsq[(size_t) s * n + mem[j]] = dsq[j];

    free(Y);
  }

  return NULL;
}


//======================================================
// R interface.
//
// dm:   packed lower triangle of D (double).
// grp:  each sample's group, 0 .. n_groups - 1.
// opts: c(n_groups, median (0/1), n_threads).
//
// Returns list(z2, n_pos, n_neg, converged): z2 is each
// sample's signed squared distance to its group's centre
// (negative when the imaginary part is larger).
//======================================================
SEXP C_permdisp(SEXP sexp_dm, SEXP sexp_grp, SEXP sexp_opts) {

  int n         = length(sexp_grp);
  int ng        = INTEGER(sexp_opts)[0];
  int median    = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];
  const int *grp = INTEGER(sexp_grp);

  if (n < 1 || ng < 1 || n_threads < 1 ||
      XLENGTH(sexp_dm) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_permdisp()."); // # nocov

  for (int a = 0; a < n; a++)
    if (grp[a] < 0 || grp[a] >= ng)
      error("Invalid groups in C_permdisp()."); // # nocov

  pd_n  = n;
  pd_d  = REAL(sexp_dm);
  pd_ng = ng;


  // Members of each group, in ascending order.
  int *beg = (int *) R_alloc(ng + 1, sizeof(int));
  int *mem = (int *) R_alloc(n, sizeof(int));
  int *fil = (int *) R_alloc(ng, sizeof(int));
  for (int g = 0; g <= ng; g++) beg[g] = 0;
  for (int a = 0; a < n; a++) beg[grp[a] + 1]++;
  for (int g = 0; g < ng; g++) beg[g + 1] += beg[g];
  memcpy(fil, beg, (size_t) ng * sizeof(int));
  for (int a = 0; a < n; a++) mem[fil[grp[a]]++] = a;
  pd_beg = beg;
  pd_mem = mem;

  SEXP sexp_z2 = PROTECT(allocVector(REALSXP, n));
  double *z2   = REAL(sexp_z2);
  int n_pos = 0, n_neg = 0, converged = 1;


  if (!median) {

    //--------------------------------------------------
    // Centroids, from the distances within groups.
    //--------------------------------------------------
    pd_rs = (double *) R_alloc(n, sizeof(double));

    double pairs = 0;
    for (int g = 0; g < ng; g++) {
      double len = beg[g + 1] - beg[g];
      pairs += len * len;
    }
    run_parallel(centroid_worker, threads_for(n_threads, pairs, 1e6), 100);

    for (int g = 0; g < ng; g++) {
      double len = beg[g + 1] - beg[g], tot = 0;
      for (int a = beg[g]; a < beg[g + 1]; a++) tot += pd_rs[mem[a]];
      for (int a = beg[g]; a < beg[g + 1]; a++)
        z2[mem[a]] = pd_rs[mem[a]] / len - tot / (2 * len * len);
    }
  }

  else {

    //--------------------------------------------------
    // Gower's centred matrix G (lower triangle), as
    // vegan: A = D^2, then column and row means removed.
    //--------------------------------------------------
    // R_alloc() everything first, so that an R_alloc()
    // error cannot leak the malloc()-ed n x n buffers.
    double *rmean = (double *) R_alloc(n, sizeof(double));
    double *vals  = (double *) R_alloc(n, sizeof(double));
    double *tmp   = (double *) R_alloc(2 * (size_t) n, sizeof(double));
    int    *axes  = (int *)    R_alloc(n, sizeof(int));
    double *scl   = (double *) R_alloc(n, sizeof(double));
    pd_dsq   = (double *) R_alloc(2 * (size_t) n, sizeof(double));
    pd_iters = (int *)    R_alloc(2 * (size_t) ng, sizeof(int));

    double *G = (double *) malloc((size_t) n * n * sizeof(double));
    if (!G) error("Not enough memory for a %d x %d matrix.", n, n); // # nocov

    for (int i = 0; i < n; i++) rmean[i] = 0;
    for (int j = 0; j < n; j++) {
      G[(size_t) j * n + j] = 0;
      for (int i = j + 1; i < n; i++) {
        double x = pd_d[dist_col_start(j, n) + (size_t) (i - j - 1)];
        x *= x;
        G[(size_t) j * n + i] = x;
        rmean[i] += x;
        rmean[j] += x;
      }
    }
    double gmean = 0;
    for (int i = 0; i < n; i++) { rmean[i] /= n; gmean += rmean[i]; }
    gmean /= n;
    for (int j = 0; j < n; j++)
      for (int i = j; i < n; i++)
        G[(size_t) j * n + i] = -0.5 * (G[(size_t) j * n + i] - rmean[i] - rmean[j] + gmean);

    double *vecs = (double *) malloc((size_t) n * n * sizeof(double));
    if (!vecs) { free(G); error("Not enough memory for a %d x %d matrix.", n, n); } // # nocov

    int ql_failed = lz_sym_eigen(G, n, vals, vecs, tmp);
    free(G);
    if (ql_failed) { // # nocov start
      free(vecs);
      error("The eigendecomposition for the spatial medians did not converge.");
    } // # nocov end


    //--------------------------------------------------
    // vegan's axes: |lambda| > max(TOL, TOL * lambda_1).
    //--------------------------------------------------
    double tol = sqrt(DBL_EPSILON);
    double cut = (tol * vals[0] > tol) ? tol * vals[0] : tol;

    for (int k = 0; k < n; k++) {
      scl[k] = sqrt(fabs(vals[k]));
      if (fabs(vals[k]) > cut && vals[k] > 0) axes[n_pos++] = k;
    }
    for (int k = 0; k < n; k++)
      if (fabs(vals[k]) > cut && vals[k] < 0) axes[n_pos + n_neg++] = k;

    pd_vec   = vecs;
    pd_scl   = scl;
    pd_axes  = axes;
    pd_npos  = n_pos;
    pd_nneg  = n_neg;
    for (size_t i = 0; i < 2 * (size_t) n; i++) pd_dsq[i] = 0;
    pd_failed = 0;

    // About 30 Weiszfeld iterations per median.
    int nt = threads_for(n_threads, 30.0 * n * (n_pos + n_neg), 1e6);
    if (nt > 2 * ng) nt = 2 * ng;
    run_parallel(median_worker, nt, 100);

    free(vecs);
    if (pd_failed) error("Out of memory in C_permdisp()."); // # nocov

    for (int t = 0; t < 2 * ng; t++)
      if (pd_iters[t] < 0) converged = 0;

    for (int i = 0; i < n; i++)
      z2[i] = pd_dsq[i] - pd_dsq[(size_t) n + i];
  }

  const char *names[] = { "z2", "n_pos", "n_neg", "converged", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_z2);
  SET_VECTOR_ELT(sexp_result, 1, ScalarInteger(n_pos));
  SET_VECTOR_ELT(sexp_result, 2, ScalarInteger(n_neg));
  SET_VECTOR_ELT(sexp_result, 3, ScalarLogical(converged));

  UNPROTECT(2);
  return sexp_result;
}



//======================================================
// Permutation statistics.
//======================================================

static int           pp_n;
static int           pp_ng;
static const double *pp_res;    // ANOVA residuals of the distances.
static const double *pp_z;      // The distances.
static const int    *pp_grp;
static const double *pp_cnt;    // Group sizes.
static int           pp_pair;   // Pairwise t statistics?
static int           pp_npair;
static perm_source_t pp_ps;     // Supplied or random permutations.
static int           pp_k0;     // Permutations for this round:
static int           pp_k1;     // ... pp_k0 to pp_k1 - 1 (1-based).
static double       *pp_f;      // n_perm.
static double       *pp_t;      // n_pair x n_perm.
static int          *pp_iwork;  // n ints per thread.
static double       *pp_work;   // 3 * n_groups doubles per thread.


static void *perm_worker(void *arg) {

  worker_t *w  = (worker_t *) arg;
  int       n  = pp_n, ng = pp_ng;
  int      *p  = pp_iwork + (size_t) w->i * n;
  double   *mu = pp_work  + (size_t) w->i * 3 * ng;
  double   *ss = mu + ng;
  double   *mz = ss + ng;

  for (int k = pp_k0 + w->i; k < pp_k1; k += w->n) {

    // Permuted row a is original row p[a], as vegan's
    // resids[take].
    perm_source_get(&pp_ps, k, p);

    //--------------------------------------------------
    // F of the permuted residuals.
    //--------------------------------------------------
    double tot = 0;
    for (int g = 0; g < ng; g++) mu[g] = 0;
    for (int a = 0; a < n; a++) {
      double x = pp_res[p[a]];
      mu[pp_grp[a]] += x;
      tot += x;
    }
    for (int g = 0; g < ng; g++) mu[g] /= pp_cnt[g];
    tot /= n;

    // Summed over samples, not groups, so the level order can't change it.
    double mss = 0, rss = 0;
    for (int a = 0; a < n; a++) {
      double t = mu[pp_grp[a]] - tot;
      mss += t * t;
    }
    for (int a = 0; a < n; a++) {
      double t = pp_res[p[a]] - mu[pp_grp[a]];
      rss += t * t;
    }
    pp_f[k - 1] = (mss / (ng - 1)) / (rss / (n - ng));


    //--------------------------------------------------
    // t statistics of the permuted distances, for each
    // pair of groups, in combn() order.
    //--------------------------------------------------
    if (pp_pair) {

      for (int g = 0; g < ng; g++) { mz[g] = 0; ss[g] = 0; }
      for (int a = 0; a < n; a++) mz[pp_grp[a]] += pp_z[p[a]];
      for (int g = 0; g < ng; g++) mz[g] /= pp_cnt[g];
      for (int a = 0; a < n; a++) {
        double t = pp_z[p[a]] - mz[pp_grp[a]];
        ss[pp_grp[a]] += t * t;
      }

      double *out = pp_t + (size_t) (k - 1) * pp_npair;
      int     c   = 0;
      for (int u = 0; u < ng - 1; u++) {
        for (int v = u + 1; v < ng; v++) {
          double nu = pp_cnt[u], nv = pp_cnt[v];
          double pooled = sqrt((ss[u] + ss[v]) / (nu + nv - 2));
          out[c++] = (mz[u] - mz[v]) / (pooled * sqrt(1 / nu + 1 / nv));
        }
      }
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// res:    the residuals (each distance minus its group's
//         mean).
// z:      the distances.
// grp:    each sample's group, 0 .. n_groups - 1; every
//         group non-empty.
// perms:  n x n_perm integer matrix of 0-based
//         permutations, or NULL for random ones.
// strata: each sample's stratum, 0 .. n_strata - 1, or
//         NULL.
// opts:   c(n_groups, n_perm, seed, n_threads, pairwise).
//
// Returns list(f = n_perm F statistics, t = n_pairs x
// n_perm t statistics, or NULL).
//======================================================
SEXP C_permdisp_perm(
    SEXP sexp_res, SEXP sexp_z, SEXP sexp_grp, SEXP sexp_perms,
    SEXP sexp_strata, SEXP sexp_opts ) {

  int n         = length(sexp_grp);
  int ng        = INTEGER(sexp_opts)[0];
  int n_perm    = INTEGER(sexp_opts)[1];
  int seed      = INTEGER(sexp_opts)[2];
  int n_threads = INTEGER(sexp_opts)[3];
  int pairwise  = INTEGER(sexp_opts)[4];

  if (n < 1 || ng < 1 || n_perm < 0 || n_threads < 1 ||
      length(sexp_res) != n || length(sexp_z) != n)
    error("Invalid arguments to C_permdisp_perm()."); // # nocov

  pp_n     = n;
  pp_ng    = ng;
  pp_res   = REAL(sexp_res);
  pp_z     = REAL(sexp_z);
  pp_grp   = INTEGER(sexp_grp);
  pp_pair  = pairwise;
  pp_npair = ng * (ng - 1) / 2;

  double *cnt = (double *) R_alloc(ng, sizeof(double));
  for (int g = 0; g < ng; g++) cnt[g] = 0;
  for (int a = 0; a < n; a++) {
    if (pp_grp[a] < 0 || pp_grp[a] >= ng)
      error("Invalid groups in C_permdisp_perm()."); // # nocov
    cnt[pp_grp[a]] += 1;
  }
  pp_cnt = cnt;


  perm_source_init(&pp_ps, n, sexp_perms, sexp_strata, seed);


  SEXP sexp_f = PROTECT(allocVector(REALSXP, n_perm));
  SEXP sexp_t = PROTECT(pairwise ? allocMatrix(REALSXP, pp_npair, n_perm) : R_NilValue);
  pp_f = REAL(sexp_f);
  pp_t = pairwise ? REAL(sexp_t) : NULL;

  if (n_threads > n_perm) n_threads = (n_perm > 0) ? n_perm : 1;
  pp_iwork = (int *)    R_alloc((size_t) n * n_threads, sizeof(int));
  pp_work  = (double *) R_alloc((size_t) 3 * ng * n_threads, sizeof(double));


  // Rounds of permutations, checking for interrupts
  // between them: ~1e7 operations per thread each.
  double per_perm  = (double) n * (pairwise ? 4 : 2) + 1;
  double rounds    = 1e7 * n_threads / per_perm;
  int    per_round = (rounds < n_perm) ? (int) rounds : n_perm; // no overflow
  if (per_round < n_threads) per_round = n_threads;

  for (pp_k0 = 1; pp_k0 <= n_perm; pp_k0 = pp_k1) {

    pp_k1 = (n_perm + 1 - pp_k0 > per_round) ? pp_k0 + per_round : n_perm + 1;

    int cnt_k = pp_k1 - pp_k0;
    int nt    = threads_for(n_threads, per_perm * cnt_k, 1e6);
    if (nt > cnt_k) nt = cnt_k;
    run_parallel(perm_worker, nt, 100);

    if (pp_k1 <= n_perm) R_CheckUserInterrupt();
  }

  const char *names[] = { "f", "t", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_f);
  SET_VECTOR_ELT(sexp_result, 1, sexp_t);

  UNPROTECT(3);
  return sexp_result;
}
