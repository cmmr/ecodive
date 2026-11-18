// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#ifndef ECODIVE_H
#define ECODIVE_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h> // R_registerRoutines, R_useDynamicSymbols

#include <inttypes.h> // uint32_t, uint64_t
#include <math.h>     // exp, fabs, floor, log, lgamma, pow, round, sqrt
#include <string.h>   // memcpy, memset, strlen, strncpy
#include <stdlib.h>   // calloc, free, malloc, NULL, qsort, strtod
#include <stddef.h>   // size_t
#include "pcg_basic.h"

// A uniform integer in 0 .. bound - 1, without modulo bias: pcg-cpp's
// bounded draw, which uwot uses for negative samples.
static inline uint32_t pcg32_bounded(pcg32_random_t *rng, uint32_t bound) {
  uint32_t threshold = (0u - bound) % bound;
  for (;;) {
    uint32_t r = pcg32_random_r(rng);
    if (r >= threshold) return r % bound;
  }
}


// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

typedef void *(*pthread_func_t)(void *);


// ecomatrix data structure
typedef struct {
  int     n_samples;
  int     n_otus;
  int     nnz;
  int    *sam_vec;
  int    *pos_vec;
  int    *otu_vec;
  double *val_vec;
  double *clr_vec;
  SEXP    sexp_sample_names;
} ecomatrix_t;

// ecotree data structures
typedef struct {
  int    edge;
  int    parent;
  double length;
} node_t;

typedef struct {
  int     n_edges;
  double *edge_lengths;
  node_t *node_vec;
} ecotree_t;

typedef struct {
  int i;
  int n;
} worker_t;


/* --- ecomatrix.c --- */
ecomatrix_t* new_ecomatrix(SEXP sexp_matrix, SEXP sexp_margin);
double* rw_val_vec(ecomatrix_t *em);
double* rw_clr_vec(ecomatrix_t *em);

/* --- ecotree.c --- */
ecotree_t* new_ecotree(SEXP sexp_phylo_tree);

/* --- get.c --- */
SEXP get(SEXP, const char *);
void set(SEXP, const char *, SEXP);

/* --- lanczos.c --- */
// W = B * V for `m` column-major columns of length n.
typedef void (*lanczos_op_t)(const double *V, double *W, int m, void *ctx);

typedef struct {
  int          n;         // In: vector length (B is n x n).
  int          k;         // In: eigenpairs wanted, 1 <= k <= dims.
  int          dims;      // In: dimension of B's range, at most.
  int          deflate;   // In: keep vectors orthogonal to ones.
  const double *basis;    // In: n x n_basis orthonormal columns to
  int          n_basis;   //     keep vectors orthogonal to, or NULL/0.
  double       scale;     // In: trace(B), for convergence tests.
  lanczos_op_t op;        // In: product function.
  void        *ctx;       // In: passed on to op.
  double      *vecs;      // Out: n x k unit eigenvectors.
  double      *vals;      // Out: k eigenvalues, decreasing.
  double      *resid;     // Out: k relative residual norms.
  int          steps;     // Out: products with blocks.
  int          converged; // Out: 1 if converged.

  // Optional. Left at 0 / NULL, each keeps the default,
  // and lanczos() runs exactly as without these fields.
  int          max_basis; // In: basis size before a thick restart (default LZ_BASIS, 60).
  int          rr_every;  // In: Rayleigh-Ritz every this many steps, and
                          //     whenever a restart or the end is due (default 1).
  double       tol;       // In: residual tolerance, replacing both LZ_TOL and
                          //     LZ_TOL_STALL (default 0: those).
  int          early_stop;// In: if set, also stop (as converged) once the largest
  double       stop_above;//     Ritz value is above stop_above with its residual
  double       stop_tol;  //     at most stop_tol (INFINITY: at once). Ritz values
                          //     never exceed B's largest eigenvalue, so this
                          //     proves that eigenvalue is above stop_above.
  uint64_t     stream;    // In: PCG stream of the starting block (default 54).
  double      *ritz;      // Out: if not NULL, up to n_ritz of the final Ritz
  int          n_ritz;    //      vectors (n x n_ritz, largest first; signs not
  int          n_ritz_out;//      fixed), and how many were written.
  int          max_steps; // In: give up (not converged) after this many
                          //     products with blocks (default LZ_MAX_STEPS).
} lanczos_t;

void   lanczos(lanczos_t *lz);
double lz_dot(const double *x, const double *y, int n);
void   lz_center(double *x, int n);
int    lz_sym_eigen(const double *A, int s, double *vals, double *vecs, double *tmp);
int    lz_gen_eigenvalues(double *A, int s, double *wr, double *wi, double *ort);
int    lz_sign_index(const double *x, int n);
void   lz_combine(const double *V, int n, int s, const double *Y, int m, double *X);
void   lz_project_out(const double *F, int f, double *x, int n);
void   lz_fill_random(double *x, size_t n, pcg32_random_t *rng);
int    lz_orthonormalize(const double *V, int s, double *R, int m, int n,
                         int deflate, const double *F, int f, pcg32_random_t *rng);
void   lz_extend_h(const double *V, const double *W, double *H, int M, int s_w, int s_v, int n);
void   lz_restart_h(const double *V, const double *W, double *H, int M, int keep, int n);

/* --- dist objects --- */
// Offset of element (j + 1, j) in a dist object's packed
// lower triangle for n samples. Column j's elements, rows
// j + 1 to n - 1, follow it in order.
static inline size_t dist_col_start(int j, int n) {
  return (size_t) j * (size_t) n - (size_t) j * (size_t) (j + 1) / 2;
}

/* --- permute.c --- */
// Where permutation tests get their permutations: a
// supplied matrix, or random shuffles within strata, one
// PCG stream (seed, k) per permutation k. See permute.c.
typedef struct {
  int           n;      // Samples.
  const int    *perms;  // n x n_perm 0-based permutations, or NULL.
  uint64_t      seed;
  const int    *order;  // Samples, ordered by stratum.
  const int    *sbeg;   // Stratum s is order[sbeg[s]] ..
  int           ns;     // ... order[sbeg[s + 1] - 1].
} perm_source_t;

void perm_source_init(perm_source_t *ps, int n, SEXP sexp_perms, SEXP sexp_strata, int seed);
void perm_source_get(const perm_source_t *ps, int k, int *p);

/* --- gower.c --- */
// Products with Gower's centered matrix G = -1/2 J D^2 J.
typedef struct {
  int           n;          // Samples.
  const double *dist;       // Packed lower triangle of D.
  const double *F;          // n x f orthonormal columns to project
  int           f;          //   out of each product, or NULL/0.
  int           max_m;      // Vectors per pass over `dist`.
  double       *xr, *yr;    // n x max_m scratch.
  int           n_threads;
  double        shift;      // Distances are read as d + shift (default 0;
                            //   set with gower_set_shift()).
  double        max_d;      // Out: the largest distance as given.
  double       *yr2;        // n x max_m scratch for gower_mult2(), or NULL.
} gower_t;

double gower_init(gower_t *g, const double *dist, int n, int max_m, int n_threads);
void   gower_mult(const double *V, double *W, int m, void *ctx);
double gower_set_shift(gower_t *g, double shift);
void   gower_mult2(const double *V, double *W, double *W2, int m, gower_t *g);

/* --- knn.c --- */
void knn_from_dist(const double *dist, int n, int k, int include_self,
                   int *idx, double *dd, int n_threads);

/* --- memory.c --- */
void  init_n_ptrs(int n);
void* safe_malloc(size_t bytes);
int   is_safe_ptr(void *ptr);
void  free_all(void);
void* free_one(void *ptr);
void* maybe_free_one(void *ptr);

/* --- normalize.c --- */
void normalize(ecomatrix_t *em, int norm, int n_threads, double pseudocount_);

/* --- pca.c: the centred (scaled) table, shared with rda() and plsda() --- */
typedef struct {
  int           n, p, nnz;  // Samples, features, non-zeros.
  const double *w_feat;     // 1 / sd per feature (scale), or NULL.
  void         *ctx;        // For pca_gram_mult().
} pca_table_t;

double pca_table_init(SEXP counts, SEXP margin, int norm, double pseudocount,
                      int scale, int n_threads, pca_table_t *t);
void   pca_gram_mult(const double *V, double *W, int m, void *ctx);
void   pca_cross(const double *U, int m, double *L, int n_threads);
void   pca_mult(const double *V, int m, double *S, int n_threads);

/* --- plsda.c: PLS on pca.c's table, for plsda() (and oplsda()) --- */
// The deflated table X_a = A - T P', A being pca_table_init()'s.
typedef struct {
  int     n, p;       // Samples, features.
  int     k, max_k;   // Deflations so far; room in T and P.
  double *T;          // n x max_k scores (column-major).
  double *P;          // p x max_k loadings.
  int     n_threads;
} pls_x_t;

void   pls_x_mult(const pls_x_t *x, const double *v, double *out);
void   pls_xt_mult(const pls_x_t *x, const double *U, int m, double *L);
void   pls_x_deflate(pls_x_t *x, const double *t, const double *p);
double pls_weights(const pls_x_t *x, const double *Y, int g, double *M, double *w);

/* --- cca.c: the chi-square table of cca() --- */
double        cca_table_init(SEXP counts, SEXP margin, int n_threads, int *n,
                             int *p, int **empty);
const double *cca_sqrt_r(void);
const double *cca_sqrt_c(void);
void          cca_gram_mult(const double *V, double *W, int m, void *ctx);
void          cca_cross(const double *U, int m, double *L);
void          cca_perm_inertia(const double *x, int px, perm_source_t *ps,
                               int n_perm, double *out);
void          cca_perm_inertia_fixed(const double *Q, int rq, int rz,
                                     perm_source_t *ps, int n_perm,
                                     double *out_all, double *out_z);

/* --- parallel.c --- */
void run_parallel(pthread_func_t func, int n_threads, int n_tasks);
int  user_interrupted(void);
void stop_thread_pool(void);
SEXP C_stop_thread_pool(void);

// How many of n_threads are worth starting for `work` units, when each
// thread should get at least `min_work` (about a millisecond's worth):
// starting and joining threads costs more than smaller shares save, and
// far more when the cores are busy with other processes. Results are the
// same for any number of threads.
static inline int threads_for (int n_threads, double work, double min_work) {
  double t = work / min_work;
  if (!(t >= 2)) return 1;
  return t < n_threads ? (int) t : n_threads;
}


#endif
