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
} lanczos_t;

void   lanczos(lanczos_t *lz);
double lz_dot(const double *x, const double *y, int n);
void   lz_center(double *x, int n);
int    lz_sym_eigen(const double *A, int s, double *vals, double *vecs, double *tmp);
int    lz_sign_index(const double *x, int n);
void   lz_combine(const double *V, int n, int s, const double *Y, int m, double *X);
void   lz_project_out(const double *F, int f, double *x, int n);

/* --- dist objects --- */
// Offset of element (j + 1, j) in a dist object's packed
// lower triangle for n samples. Column j's elements, rows
// j + 1 to n - 1, follow it in order.
static inline size_t dist_col_start(int j, int n) {
  return (size_t) j * (size_t) n - (size_t) j * (size_t) (j + 1) / 2;
}

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
} gower_t;

double gower_init(gower_t *g, const double *dist, int n, int max_m, int n_threads);
void   gower_mult(const double *V, double *W, int m, void *ctx);

/* --- memory.c --- */
void  init_n_ptrs(int n);
void* safe_malloc(size_t bytes);
int   is_safe_ptr(void *ptr);
void  free_all(void);
void* free_one(void *ptr);
void* maybe_free_one(void *ptr);

/* --- normalize.c --- */
void normalize(ecomatrix_t *em, int norm, int n_threads, double pseudocount_);

/* --- parallel.c --- */
void run_parallel(pthread_func_t func, int n_threads, int n_tasks);
int  user_interrupted(void);
void stop_thread_pool(void);

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
