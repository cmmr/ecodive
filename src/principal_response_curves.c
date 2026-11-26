// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Permutation tests of principal response curves (a partial
 * RDA, Y ~ treatment:time + Condition(time)): the first
 * constrained eigenvalue and the constrained sum of squares,
 * for the observed design and for each permutation.
 *
 * With A the centred (and normalized) table of pca.c, Qz an
 * orthonormal basis of the conditions, and Q an orthonormal
 * basis of the constraints after removing Qz, the model's
 * constrained axes are the eigenpairs of the m x m matrix
 *
 *     C = (Q' A) (Q' A)'
 *
 * so trace(C) is the constrained sum of squares and C's
 * largest eigenvalue the first axis's. Both come from one
 * product L = A' Q (p x m) per permutation, which costs
 * O(m * nnz), and an m x m eigenproblem; no n x n matrix is
 * formed.
 *
 * Permutations, as permanova.c's (and vegan's anova.cca()
 * under the reduced model): the design's rows are scattered,
 * v[p[a]] = x[a] (vegan's Y[p, ]), then made orthogonal to
 * Qz and orthonormalized within the design (modified
 * Gram-Schmidt, twice). A column with less than PRC_TOL of
 * its norm left is aliased and contributes 0.
 *
 * Determinism: each permutation's design and its m x m
 * matrix are computed by one thread, in a fixed order, and
 * A' Q by pca_cross(), whose sums run in a fixed order for
 * each column, whatever the batch of columns it is part of.
 * Results are bit-identical for any number of threads.
 */

#include "ecodive.h"

#define PRC_TOL 1e-7

// Columns of A' Q per batch (several permutations' designs).
#define PRC_COLS 256


//------------------------------------------------------
// Shared with the worker threads.
//------------------------------------------------------
static int           prc_n, prc_p, prc_m, prc_pz;
static const double *prc_x;     // n x m design columns.
static const double *prc_qz;    // n x pz, or NULL.
static perm_source_t prc_ps;
static int           prc_k0;    // First permutation of the batch.
static int           prc_nb;    // Permutations in the batch.
static double       *prc_u;     // n x (nb * m): the designs.
static double       *prc_l;     // p x (nb * m): A' Q.
static double       *prc_work;  // Per-thread workspace ...
static size_t        prc_wlen;  // ... of this many doubles.
static int          *prc_iwork; // Per-thread permutation.
static double       *prc_out;   // 2 x total: lambda1, trace.


//------------------------------------------------------
// Make v (length n) orthogonal to Qz and to the nb columns
// of basis, then normalize it. Returns 0 (and zeroes v) if
// less than PRC_TOL of its norm is left.
//------------------------------------------------------
static int prc_orthonormalize(double *v, const double *basis, int nb, int n) {

  double norm0 = 0;
  for (int i = 0; i < n; i++) norm0 += v[i] * v[i];
  norm0 = sqrt(norm0);

  for (int pass = 0; pass < 2; pass++) {
    for (int c = 0; c < prc_pz + nb; c++) {
      const double *q = (c < prc_pz)
        ? prc_qz + (size_t) c * n
        : basis + (size_t) (c - prc_pz) * n;
      double dot = 0;
      for (int i = 0; i < n; i++) dot += q[i] * v[i];
      for (int i = 0; i < n; i++) v[i] -= dot * q[i];
    }
  }

  double norm = 0;
  for (int i = 0; i < n; i++) norm += v[i] * v[i];
  norm = sqrt(norm);

  if (!(norm > PRC_TOL * norm0)) {
    memset(v, 0, (size_t) n * sizeof(double));
    return 0;
  }

  for (int i = 0; i < n; i++) v[i] /= norm;
  return 1;
}


//------------------------------------------------------
// Each permutation of the batch: its orthonormal design.
//------------------------------------------------------
static void *design_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = prc_n, m = prc_m;
  int *p = prc_iwork + (size_t) thread_i * n;

  for (int b = thread_i; b < prc_nb; b += n_threads) {

    perm_source_get(&prc_ps, prc_k0 + b, p);

    double *U = prc_u + (size_t) b * m * n;
    for (int c = 0; c < m; c++) {
      double       *v = U + (size_t) c * n;
      const double *x = prc_x + (size_t) c * n;
      for (int a = 0; a < n; a++) v[p[a]] = x[a];
      prc_orthonormalize(v, U, c, n);
    }
  }

  return NULL;
}


//------------------------------------------------------
// Each permutation of the batch: C = L_b' L_b, its trace
// and its largest eigenvalue.
//------------------------------------------------------
static void *eigen_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int p = prc_p, m = prc_m;

  double *C    = prc_work + (size_t) thread_i * prc_wlen;  // m x m
  double *vals = C + (size_t) m * m;                       // m
  double *vecs = vals + m;                                 // m x m
  double *tmp  = vecs + (size_t) m * m;                    // 2 m

  for (int b = thread_i; b < prc_nb; b += n_threads) {

    const double *L = prc_l + (size_t) b * m * p;

    double tr = 0;
    for (int c = 0; c < m; c++) {
      const double *lc = L + (size_t) c * p;
      for (int d = 0; d <= c; d++) {
        const double *ld = L + (size_t) d * p;
        double s = 0;
        for (int j = 0; j < p; j++) s += lc[j] * ld[j];
        C[(size_t) d * m + c] = s;
        C[(size_t) c * m + d] = s;
      }
      tr += C[(size_t) c * m + c];
    }

    double lam = 0;
    if (m == 1) {
      lam = C[0];
    }
    else if (lz_sym_eigen(C, m, vals, vecs, tmp)) {
      lam = NA_REAL; // # nocov
    }
    else {
      lam = vals[0];
    }

    double *out = prc_out + 2 * (size_t) (prc_k0 + b);
    out[0] = lam;
    out[1] = tr;
  }

  return NULL;
}


//======================================================
// R interface.
//
// counts, margin, norm, pseudocount: as for C_rda().
// x:      n x m design columns (centred; without
//         conditions, orthonormal), as design_basis()
//         gives them.
// qz:     n x pz orthonormal basis of the conditions, or
//         NULL.
// perms:  n x n_perm 0-based permutations, or NULL.
// strata: 0-based codes, or NULL.
// ints:   c(n_perm, seed, n_threads).
//
// Returns list(stats = 2 x (n_perm + 1) matrix of lambda1
// and trace(C), sums of squares (not divided by n - 1),
// column 1 the observed; total = trace(A' A); ss_z = the
// conditions' sum of squares).
//======================================================
SEXP C_prc_perm(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_pseudocount,
    SEXP sexp_x, SEXP sexp_qz, SEXP sexp_perms, SEXP sexp_strata, SEXP sexp_ints ) {

  int    norm        = asInteger(sexp_norm);
  double pseudocount = asReal(sexp_pseudocount);
  int    n_perm      = INTEGER(sexp_ints)[0];
  int    seed        = INTEGER(sexp_ints)[1];
  int    n_threads   = INTEGER(sexp_ints)[2];
  int    m           = ncols(sexp_x);

  if (n_perm < 0 || n_threads < 1 || m < 1)
    error("Invalid arguments to C_prc_perm()."); // # nocov

  init_n_ptrs(10);

  pca_table_t tbl;
  double total = pca_table_init(
    sexp_counts, sexp_margin, norm, pseudocount, 0, n_threads, &tbl );

  int n = tbl.n, p = tbl.p;

  if (nrows(sexp_x) != n || (!isNull(sexp_qz) && nrows(sexp_qz) != n)) {
    free_all();                                  // # nocov
    error("Invalid arguments to C_prc_perm()."); // # nocov
  }

  prc_n  = n;
  prc_p  = p;
  prc_m  = m;
  prc_x  = REAL(sexp_x);
  prc_qz = isNull(sexp_qz) ? NULL : REAL(sexp_qz);
  prc_pz = isNull(sexp_qz) ? 0 : ncols(sexp_qz);
  perm_source_init(&prc_ps, n, sexp_perms, sexp_strata, seed);


  //--------------------------------------------------
  // The conditions' sum of squares: || A' Qz ||^2.
  //--------------------------------------------------
  double ss_z = 0;
  if (prc_pz > 0) {
    double *Lz = (double *) R_alloc((size_t) p * prc_pz, sizeof(double));
    pca_cross(prc_qz, prc_pz, Lz, n_threads);
    for (size_t e = 0; e < (size_t) p * prc_pz; e++) ss_z += Lz[e] * Lz[e];
  }


  //--------------------------------------------------
  // Batches of permutations (the first holds the
  // identity, permutation 0).
  //--------------------------------------------------
  int total_perms = n_perm + 1;
  int bsize       = PRC_COLS / m;
  if (bsize < 1) bsize = 1;
  // Keep A' Q under about 2^24 doubles (128 MB).
  double cap = 16777216.0 / ((double) m * p);
  if (bsize > cap) bsize = (cap >= 1) ? (int) cap : 1;
  if (bsize > total_perms) bsize = total_perms;

  prc_u     = (double *) R_alloc((size_t) n * m * bsize, sizeof(double));
  prc_l     = (double *) R_alloc((size_t) p * m * bsize, sizeof(double));
  prc_wlen  = 2 * (size_t) m * m + 3 * (size_t) m;
  prc_work  = (double *) R_alloc(prc_wlen * n_threads, sizeof(double));
  prc_iwork = (int *)    R_alloc((size_t) n * n_threads, sizeof(int));

  SEXP sexp_stats = PROTECT(allocMatrix(REALSXP, 2, total_perms));
  prc_out = REAL(sexp_stats);

  double work = (double) bsize * m * ((double) m * p + (double) n * (prc_pz + m));
  int    nt   = threads_for(n_threads, work, 1e6);

  for (prc_k0 = 0; prc_k0 < total_perms; prc_k0 += bsize) {

    prc_nb = (total_perms - prc_k0 < bsize) ? total_perms - prc_k0 : bsize;
    int ntb = (nt < prc_nb) ? nt : prc_nb;

    run_parallel(design_worker, ntb, (ntb > 1) ? 100 : 1);
    pca_cross(prc_u, prc_nb * m, prc_l, n_threads);
    run_parallel(eigen_worker, ntb, (ntb > 1) ? 100 : 1);

    if (prc_k0 + prc_nb < total_perms && user_interrupted()) { // # nocov start
      free_all();
      UNPROTECT(1);
      error("Interrupted.");
    } // # nocov end
  }

  free_all();

  const char *names[] = { "stats", "total", "ss_z", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_stats);
  SET_VECTOR_ELT(sexp_result, 1, ScalarReal(total));
  SET_VECTOR_ELT(sexp_result, 2, ScalarReal(ss_z));

  UNPROTECT(2);
  return sexp_result;
}
