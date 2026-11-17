// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * envfit(): how well variables (vectors) and factors fit an
 * ordination, for the observed data and for each permutation
 * of the samples.
 *
 * Vectors. With X the n x k centred ordination coordinates and
 * Q an orthonormal basis of its columns (from R's qr()), the
 * fitted values of a centred variable v_c are Q Q' v_c, and
 *
 *     r^2 = cor(QQ'v_c, v_c)^2 = ||Q' v_c||^2 / ||v_c||^2,
 *
 * as in vegan's vectorfit(). Permuting v's rows (v_p[i] =
 * v[p[i]]) leaves ||v_c|| unchanged, and
 *
 *     Q' v_p = sum_a Q[pinv[a], ] v[a],
 *
 * so each permutation scatters Q's rows instead (QP[a, ] =
 * Q[pinv[a], ]) and sums over v's non-zero values only. A
 * sparse table of taxa is never made dense. Q's columns sum
 * to 0 (s = Q'1 is ~1e-16), so for a variable with mean mu,
 * Q' v_c = Q' v - mu * s; columns with no zeros are centred
 * first instead, which avoids cancellation.
 *
 * CLR: a zero in sample a stands for z0[a], so each value is
 * stored as v - z0[a], and the permuted sum gains
 * sum_a QP[a, ] (z0[a] - mean(z0)) once per permutation.
 *
 * Factors. With X_c the centred coordinates, S_l the sum of
 * X_c's rows in level l, and n_l its size, the sum of squares
 * between levels is sum_l ||S_l||^2 / n_l, and
 *
 *     r^2 = 1 - SS_within / SS_total = SS_between / SS_total,
 *
 * as in vegan's factorfit(). Permuted labels: g_p[i] = g[p[i]].
 *
 * Determinism: permutation k (k >= 1) comes from its own PCG
 * stream, seeded with (seed, k), or from a supplied matrix,
 * through permute.c's perm_source_get(), as in permanova()
 * and mrpp(), so the same seed gives the same permutations. Permutation 0 is the identity. Each
 * block of permutations is done by one thread, and every sum
 * runs in a fixed order, so results are bit-identical for any
 * number of threads.
 */

#include "ecodive.h"

// Lanes per sample in the scattered basis: as many
// permutations as fit, times the basis's columns. 16
// lanes keep the sums in SSE2's 16 registers.
#define EF_W 16


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the parts of the output and workspace that each block
// or thread owns.
//------------------------------------------------------
static int           ef_n;       // Samples (rows of the ordination).
static int           ef_r;       // Columns of Q.
static const double *ef_q;       // n x r orthonormal basis, column-major.
static const double *ef_s;       // r column sums of Q.

static int           ef_m;       // Vectors.
static const int    *ef_cp;      // CSC column starts (m + 1).
static const int    *ef_ci;      // ... each value's sample.
static double       *ef_cv;      // ... values (transformed in place).
static double       *ef_madj;    // Each vector's multiple of s to subtract.
static double       *ef_iss;     // 1 / ||v_c||^2, or 0 if constant.
static const double *ef_z0c;     // CLR: centred zero values, or NULL.
static double       *ef_qtv;     // r x m: observed Q' v_c.

static int           ef_f;       // Factors.
static int           ef_k;       // Columns of X_c.
static const double *ef_xc;      // n x k centred coordinates.
static const int    *ef_grp;     // n x f level codes, 0-based.
static const int    *ef_lofs;    // Factor f's levels start at ef_lofs[f].
static const double *ef_inl;     // 1 / level size, by ef_lofs.
static double        ef_itot;    // 1 / SS_total of X_c.

static perm_source_t ef_ps;     // Permutations (permute.c).

static int           ef_total;   // Permutations, with the identity.
static int           ef_bsize;   // Permutations per block.
static int           ef_wid;     // Lanes per sample (>= ef_bsize * ef_r).
static int           ef_blk0;    // This round's blocks:
static int           ef_blk1;    // ... ef_blk0 .. ef_blk1 - 1.
static double       *ef_r2v;     // ef_total x m.
static double       *ef_r2f;     // ef_total x f.

static double       *ef_work;    // Per-thread doubles ...
static size_t        ef_wlen;    // ... this many each.
static int          *ef_iwork;   // Per-thread ints ...
static size_t        ef_iwlen;   // ... this many each.


//------------------------------------------------------
// acc[0 .. EF_W - 1] += sum_e qp[ci[e] * EF_W + .] * cv[e]
// for one vector's values. GCC and clang get explicit
// two-lane vectors: compilers don't vectorize this at R's
// default -O2. Each lane does the same operations in the
// same order as the scalar fallback, so results don't
// depend on which is used.
//------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)

typedef double ef_v2d __attribute__((vector_size(16)));

static inline ef_v2d ef_load2(const double *p) { ef_v2d x; memcpy(&x, p, sizeof(x)); return x; }

#define EF_REP(F) F(0) F(1) F(2) F(3) F(4) F(5) F(6) F(7)

static inline void scatter_sum_w(
    const double *qp, const int *ci, const double *cv, int len, double *acc) {

  #define S_INIT(h) ef_v2d a##h = { 0, 0 };
  #define S_STEP(h) a##h += ef_load2(row + 2 * h) * v;
  #define S_DONE(h) acc[2 * h] += a##h[0]; acc[2 * h + 1] += a##h[1];

  EF_REP(S_INIT)
  for (int e = 0; e < len; e++) {
    const double *row = qp + (size_t) ci[e] * EF_W;
    ef_v2d        v   = { cv[e], cv[e] };
    EF_REP(S_STEP)
  }
  EF_REP(S_DONE)
}

#else

static inline void scatter_sum_w(
    const double *qp, const int *ci, const double *cv, int len, double *acc) {

  double a[EF_W];
  for (int c = 0; c < EF_W; c++) a[c] = 0;
  for (int e = 0; e < len; e++) {
    const double *row = qp + (size_t) ci[e] * EF_W;
    double        v   = cv[e];
    for (int c = 0; c < EF_W; c++) a[c] += row[c] * v;
  }
  for (int c = 0; c < EF_W; c++) acc[c] += a[c];
}

#endif

// Any width (more than EF_W columns).
static inline void scatter_sum(
    const double *qp, int wid, const int *ci, const double *cv, int len, double *acc) {

  for (int e = 0; e < len; e++) {
    const double *row = qp + (size_t) ci[e] * wid;
    double        v   = cv[e];
    for (int c = 0; c < wid; c++) acc[c] += row[c] * v;
  }
}


//------------------------------------------------------
// The statistics for each permutation in blocks
// w->i, w->i + w->n, ... of this round.
//------------------------------------------------------
static void *perm_worker(void *arg) {

  worker_t *w = (worker_t *) arg;

  int n   = ef_n;
  int r   = ef_r;
  int wid = ef_wid;
  int B   = ef_bsize;

  double *qp   = ef_work  + ef_wlen  * (size_t) w->i;  // n x wid
  double *acc  = qp + (size_t) n * wid;                // wid
  double *cb   = acc + wid;                            // wid
  double *sl   = cb + wid;                             // levels x k
  int    *p    = ef_iwork + ef_iwlen * (size_t) w->i;  // n
  int    *pinv = p + n;                                // n

  for (int blk = ef_blk0 + w->i; blk < ef_blk1; blk += w->n) {

    int k0 = blk * B;
    int kb = ef_total - k0;
    if (kb > B) kb = B;

    if (ef_m > 0)
      memset(qp, 0, (size_t) n * wid * sizeof(double));

    for (int b = 0; b < kb; b++) {

      int k = k0 + b;
      perm_source_get(&ef_ps, k, p);


      //--------------------------------------------
      // Vectors: scatter Q's rows, QP[a] = Q[pinv[a]].
      //--------------------------------------------
      if (ef_m > 0) {
        for (int i = 0; i < n; i++) pinv[p[i]] = i;
        for (int a = 0; a < n; a++) {
          double *row = qp + (size_t) a * wid + (size_t) b * r;
          int     src = pinv[a];
          for (int c = 0; c < r; c++) row[c] = ef_q[src + (size_t) c * n];
        }
      }


      //--------------------------------------------
      // Factors: level sums of the permuted labels.
      //--------------------------------------------
      for (int f = 0; f < ef_f; f++) {

        const int    *g   = ef_grp + (size_t) f * n;
        int           l0  = ef_lofs[f];
        int           nl  = ef_lofs[f + 1] - l0;
        int           kx  = ef_k;

        memset(sl, 0, (size_t) nl * kx * sizeof(double));
        for (int i = 0; i < n; i++) {
          double *s = sl + (size_t) g[p[i]] * kx;
          for (int c = 0; c < kx; c++) s[c] += ef_xc[i + (size_t) c * n];
        }

        double between = 0;
        for (int l = 0; l < nl; l++) {
          double ss = 0;
          for (int c = 0; c < kx; c++) ss += sl[l * kx + c] * sl[l * kx + c];
          between += ss * ef_inl[l0 + l];
        }

        ef_r2f[(size_t) f * ef_total + k] = between * ef_itot;
      }
    }

    if (ef_m == 0) continue;


    //----------------------------------------------
    // CLR: the zeros' share, once per permutation.
    //----------------------------------------------
    memset(cb, 0, (size_t) wid * sizeof(double));
    if (ef_z0c) {
      for (int a = 0; a < n; a++) {
        const double *row = qp + (size_t) a * wid;
        double        z   = ef_z0c[a];
        for (int c = 0; c < wid; c++) cb[c] += row[c] * z;
      }
    }


    //----------------------------------------------
    // Each vector: Q' v_p for the block's permutations.
    //----------------------------------------------
    for (int t = 0; t < ef_m; t++) {

      int    beg  = ef_cp[t];
      int    len  = ef_cp[t + 1] - beg;
      double madj = ef_madj[t];
      double iss  = ef_iss[t];

      memcpy(acc, cb, (size_t) wid * sizeof(double));

      if (wid == EF_W) scatter_sum_w(qp,      ef_ci + beg, ef_cv + beg, len, acc);
      else             scatter_sum(qp, wid, ef_ci + beg, ef_cv + beg, len, acc);

      for (int b = 0; b < kb; b++) {
        double ss = 0;
        for (int c = 0; c < r; c++) {
          double u = acc[b * r + c] - madj * ef_s[c];
          ss += u * u;
          if (k0 + b == 0) ef_qtv[(size_t) t * r + c] = u;
        }
        ef_r2v[(size_t) t * ef_total + k0 + b] = ss * iss;
      }
    }
  }

  return NULL;
}


//------------------------------------------------------
// Each vector's mean, ||v_c||^2 and constancy, for the
// columns w->i, w->i + w->n, ... Writes the transformed
// values (see the top of this file) to ef_cv, never to
// its input, so that a rerun gives the same results.
//------------------------------------------------------
static const double *st_cv;   // Values as read (ef_cv gets them transformed).
static const double *st_z0;   // CLR zero values, by sample, or NULL.
static double        st_z0m;  // ... their mean.
static int          *st_const;
static double       *st_mean;
static double       *st_ss;

static void *stats_worker(void *arg) {

  worker_t *w = (worker_t *) arg;
  int       n = ef_n;
  double   *d = ef_work + ef_wlen * (size_t) w->i; // n, for CLR

  for (int t = w->i; t < ef_m; t += w->n) {

    int     beg = ef_cp[t];
    int     len = ef_cp[t + 1] - beg;
    const int    *ci = ef_ci + beg;
    const double *cv = st_cv + beg;
    double       *ov = ef_cv + beg;
    double  mean = 0, ss = 0, lo, hi;

    if (st_z0) {
      // Dense copy: zeros stand for z0.
      for (int a = 0; a < n; a++) d[a] = st_z0[a];
      for (int e = 0; e < len; e++) d[ci[e]] = cv[e];
      for (int a = 0; a < n; a++) mean += d[a];
      mean /= n;
      lo = hi = d[0];
      for (int a = 0; a < n; a++) {
        double u = d[a] - mean;
        ss += u * u;
        if (d[a] < lo) lo = d[a];
        if (d[a] > hi) hi = d[a];
      }
      for (int e = 0; e < len; e++) ov[e] = cv[e] - st_z0[ci[e]];
      ef_madj[t] = mean - st_z0m;
    }
    else {
      for (int e = 0; e < len; e++) mean += cv[e];
      mean /= n;
      lo = hi = (len < n) ? 0 : cv[0];
      for (int e = 0; e < len; e++) {
        double u = cv[e] - mean;
        ss += u * u;
        if (cv[e] < lo) lo = cv[e];
        if (cv[e] > hi) hi = cv[e];
      }
      ss += (double) (n - len) * mean * mean;

      if (len == n) {
        // No zeros: centre the values themselves.
        for (int e = 0; e < len; e++) ov[e] = cv[e] - mean;
        ef_madj[t] = 0;
      }
      else {
        memcpy(ov, cv, (size_t) len * sizeof(double));
        ef_madj[t] = mean;
      }
    }

    st_const[t] = (lo == hi);
    st_mean[t]  = mean;
    st_ss[t]    = ss;
    ef_iss[t]   = (lo == hi || !(ss > 0)) ? 0 : 1 / ss;
  }

  return NULL;
}


//======================================================
// R interface.
//
// q:       n x r orthonormal basis of the centred
//          coordinates, or NULL if there are no vectors.
// counts:  The vectors: a table in any form ecodive's
//          counts take (dense, Matrix, slam), or NULL.
// margin, norm, pseudocount: as for counts elsewhere.
// map:     For each sample of counts, its row in the
//          ordination (0-based), or -1 to leave it out.
// xc:      n x k centred coordinates, for factors.
// grp:     n x f integer matrix of 0-based level codes, or
//          NULL if there are no factors.
// nlev:    Levels of each factor.
// perms:   n x n_perm 0-based permutations, or NULL.
// strata:  Each sample's stratum, 0-based, or NULL.
// opts:    c(n_perm, seed, n_threads).
//
// Returns list(r2v, r2f, qtv, mean, ss, const):
//   r2v    (n_perm + 1) x m r^2 of each vector, row 1
//          observed.
//   r2f    (n_perm + 1) x f r^2 of each factor.
//   qtv    r x m: Q' v_c, observed.
//   mean, ss, const: each vector's mean, ||v_c||^2, and
//          whether it is constant.
//======================================================
SEXP C_envfit(
    SEXP sexp_q, SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm,
    SEXP sexp_pseudocount, SEXP sexp_map, SEXP sexp_xc, SEXP sexp_grp,
    SEXP sexp_nlev, SEXP sexp_perms, SEXP sexp_strata, SEXP sexp_opts ) {

  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];
  int n         = nrows(sexp_xc);
  int kx        = ncols(sexp_xc);

  if (n < 2 || n_perm < 0 || n_threads < 1 || kx < 1)
    error("Invalid arguments to C_envfit()."); // # nocov

  ef_n     = n;
  ef_total = n_perm + 1;
  ef_xc    = REAL(sexp_xc);
  ef_k     = kx;


  //--------------------------------------------------
  // Vectors: a CSC copy with rows = ordination rows.
  //--------------------------------------------------
  ef_m = 0; ef_r = 0; ef_z0c = NULL;
  double *z0 = NULL;
  double  z0m = 0;

  if (!isNull(sexp_counts)) {

    ef_q = REAL(sexp_q);
    ef_r = ncols(sexp_q);
    if (nrows(sexp_q) != n || ef_r < 1)
      error("Invalid arguments to C_envfit()."); // # nocov

    int    norm        = asInteger(sexp_norm);
    double pseudocount = asReal(sexp_pseudocount);

    init_n_ptrs(10);
    ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

    int        ns   = em->n_samples;
    int        m    = em->n_otus;
    const int *map  = INTEGER(sexp_map);

    if (XLENGTH(sexp_map) != ns) {
      free_all();                                 // # nocov
      error("Invalid arguments to C_envfit()."); // # nocov
    }

    // Before normalizing, which could hide them (as 'binary'
    // did); and again after, for what `norm` makes.
    for (int pass = 0; pass < 2; pass++) {
      if (pass == 1 && norm) normalize(em, norm, n_threads, pseudocount);
      for (int j = 0; j < ns; j++) {
        if (map[j] < 0) continue;
        for (int e = em->pos_vec[j]; e < em->pos_vec[j + 1]; e++)
          if (!R_FINITE(em->val_vec[e])) {
            free_all();
            error(pass == 0
              ? "`counts` must not contain NA, NaN, or infinite values."
              : "`counts` must not contain NA, NaN, or infinite values (after `norm`).");
          }
      }
    }

    int *cp = (int *) R_alloc((size_t) m + 1, sizeof(int));
    memset(cp, 0, ((size_t) m + 1) * sizeof(int));

    for (int j = 0; j < ns; j++) {
      if (map[j] < 0) continue;
      for (int e = em->pos_vec[j]; e < em->pos_vec[j + 1]; e++)
        cp[em->otu_vec[e] + 1]++;
    }
    for (int t = 0; t < m; t++) cp[t + 1] += cp[t];

    size_t  nnz  = (size_t) cp[m];
    int    *ci   = (int *)    R_alloc(nnz + 1, sizeof(int));
    double *cv   = (double *) R_alloc(nnz + 1, sizeof(double));
    int    *fill = (int *)    R_alloc((size_t) m + 1, sizeof(int));
    memcpy(fill, cp, (size_t) m * sizeof(int));

    for (int j = 0; j < ns; j++) {
      int a = map[j];
      if (a < 0) continue;
      for (int e = em->pos_vec[j]; e < em->pos_vec[j + 1]; e++) {
        int t = em->otu_vec[e];
        ci[fill[t]] = a;
        cv[fill[t]] = em->val_vec[e];
        fill[t]++;
      }
    }

    if (em->clr_vec) {
      z0 = (double *) R_alloc((size_t) n, sizeof(double));
      for (int j = 0; j < ns; j++)
        if (map[j] >= 0) z0[map[j]] = em->clr_vec[j];
      for (int a = 0; a < n; a++) z0m += z0[a];
      z0m /= n;
      double *z0c = (double *) R_alloc((size_t) n, sizeof(double));
      for (int a = 0; a < n; a++) z0c[a] = z0[a] - z0m;
      ef_z0c = z0c;
    }

    free_all();

    ef_m  = m;
    ef_cp = cp;
    ef_ci = ci;
    st_cv = cv;
    ef_cv = (double *) R_alloc(nnz + 1, sizeof(double));

    double *s = (double *) R_alloc((size_t) ef_r, sizeof(double));
    for (int c = 0; c < ef_r; c++) {
      s[c] = 0;
      for (int a = 0; a < n; a++) s[c] += ef_q[a + (size_t) c * n];
    }
    ef_s = s;
  }

  SEXP sexp_mean  = PROTECT(allocVector(REALSXP, ef_m));
  SEXP sexp_ss    = PROTECT(allocVector(REALSXP, ef_m));
  SEXP sexp_const = PROTECT(allocVector(LGLSXP,  ef_m));
  SEXP sexp_qtv   = PROTECT(allocMatrix(REALSXP, ef_r, ef_m));
  ef_madj = (double *) R_alloc((size_t) ef_m + 1, sizeof(double));
  ef_iss  = (double *) R_alloc((size_t) ef_m + 1, sizeof(double));
  ef_qtv  = REAL(sexp_qtv);


  //--------------------------------------------------
  // Factors.
  //--------------------------------------------------
  ef_f = 0;
  int n_lev = 0;

  if (!isNull(sexp_grp)) {

    ef_f   = ncols(sexp_grp);
    ef_grp = INTEGER(sexp_grp);
    const int *nlev = INTEGER(sexp_nlev);

    int *lofs = (int *) R_alloc((size_t) ef_f + 1, sizeof(int));
    lofs[0] = 0;
    for (int f = 0; f < ef_f; f++) lofs[f + 1] = lofs[f] + nlev[f];
    n_lev = lofs[ef_f];

    double *inl = (double *) R_alloc((size_t) n_lev + 1, sizeof(double));
    for (int l = 0; l < n_lev; l++) inl[l] = 0;
    for (int f = 0; f < ef_f; f++)
      for (int i = 0; i < n; i++) {
        int g = ef_grp[(size_t) f * n + i];
        if (g < 0 || g >= nlev[f]) error("Invalid arguments to C_envfit()."); // # nocov
        inl[lofs[f] + g] += 1;
      }
    for (int l = 0; l < n_lev; l++) inl[l] = (inl[l] > 0) ? 1 / inl[l] : 0;

    ef_lofs = lofs;
    ef_inl  = inl;

    double tot = 0;
    for (size_t e = 0; e < (size_t) n * kx; e++) tot += ef_xc[e] * ef_xc[e];
    ef_itot = (tot > 0) ? 1 / tot : 0;
  }


  //--------------------------------------------------
  // Permutations: supplied, or shuffles within strata.
  // Permuted row a is original row p[a], as vegan's
  // P[perm[k, ], ]; the workers scatter Q's rows (see
  // the top of this file).
  //--------------------------------------------------
  perm_source_init(&ef_ps, n, sexp_perms, sexp_strata, seed);


  //--------------------------------------------------
  // Blocks and per-thread workspace.
  //--------------------------------------------------
  if (ef_r == 0 || ef_r > EF_W) {
    ef_bsize = (ef_r == 0) ? EF_W : 1;
    ef_wid   = (ef_r + 3) / 4 * 4;
  }
  else {
    ef_bsize = EF_W / ef_r;
    ef_wid   = EF_W;
  }

  int n_blocks = (ef_total + ef_bsize - 1) / ef_bsize;
  int nt_max   = n_threads;
  if (nt_max < 1) nt_max = 1;

  int maxlev = 0;
  for (int f = 0; f < ef_f; f++)
    if (ef_lofs[f + 1] - ef_lofs[f] > maxlev) maxlev = ef_lofs[f + 1] - ef_lofs[f];

  ef_wlen  = (size_t) n * ef_wid + 2 * (size_t) ef_wid + (size_t) maxlev * kx + (size_t) n + 1;
  ef_iwlen = 2 * (size_t) n;
  ef_work  = (double *) R_alloc(ef_wlen  * nt_max, sizeof(double));
  ef_iwork = (int *)    R_alloc(ef_iwlen * nt_max, sizeof(int));


  //--------------------------------------------------
  // Each vector's mean and sum of squares.
  //--------------------------------------------------
  st_z0    = z0;
  st_z0m   = z0m;
  st_const = LOGICAL(sexp_const);
  st_mean  = REAL(sexp_mean);
  st_ss    = REAL(sexp_ss);
  if (ef_m > 0)
    run_parallel(stats_worker, n_threads,
                 ((double) ef_m * n > 1e5) ? 100 : 1);


  //--------------------------------------------------
  // Rounds of blocks, checking for interrupts between
  // them. Each round has ~1e8 multiply-adds per thread.
  //--------------------------------------------------
  SEXP sexp_r2v = PROTECT(allocMatrix(REALSXP, ef_total, ef_m));
  SEXP sexp_r2f = PROTECT(allocMatrix(REALSXP, ef_total, ef_f));
  ef_r2v = REAL(sexp_r2v);
  ef_r2f = REAL(sexp_r2f);

  double nnz      = (ef_m > 0) ? (double) ef_cp[ef_m] : 0;
  double blk_work = (double) ef_wid * (nnz + 2 * n) + (double) ef_bsize * n * (ef_f * kx + 2) + 1;
  int    per_round = (int) (1e8 * n_threads / blk_work);
  if (per_round < n_threads) per_round = n_threads;

  int tasks_hint = (blk_work > 1e4) ? 100 : 1;

  for (ef_blk0 = 0; ef_blk0 < n_blocks; ef_blk0 = ef_blk1) {

    ef_blk1 = ef_blk0 + per_round;
    if (ef_blk1 > n_blocks) ef_blk1 = n_blocks;

    int nt = (ef_blk1 - ef_blk0 < n_threads) ? ef_blk1 - ef_blk0 : n_threads;
    run_parallel(perm_worker, nt, (ef_blk1 - ef_blk0 > 1) ? tasks_hint : 1);

    if (ef_blk1 < n_blocks) R_CheckUserInterrupt();
  }

  const char *names[] = { "r2v", "r2f", "qtv", "mean", "ss", "const", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_r2v);
  SET_VECTOR_ELT(sexp_result, 1, sexp_r2f);
  SET_VECTOR_ELT(sexp_result, 2, sexp_qtv);
  SET_VECTOR_ELT(sexp_result, 3, sexp_mean);
  SET_VECTOR_ELT(sexp_result, 4, sexp_ss);
  SET_VECTOR_ELT(sexp_result, 5, sexp_const);

  UNPROTECT(7);
  return sexp_result;
}
