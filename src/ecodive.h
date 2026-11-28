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
  int          n_threads; // In: the threads the caller may use (default 0:
                          //     one). Each step on blocks of vectors takes
                          //     as many as its work is worth (threads_for());
                          //     results are the same for any number.
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
  int           n_threads;  // Threads for a pass of GOWER_W vectors.
  int           n_req;      // Threads asked for.
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

/* --- mantel.c: Kendall's tau-b pieces, shared with bioenv.c --- */
int64_t kendall_sort_count(int *a, int *buf, size_t len);
double  kendall_tied_sorted(const int *a, size_t len);

/* --- rng_dist.c: continuous distributions on PCG --- */
typedef struct {
  pcg32_random_t pcg;
  int            has_spare;   // rng_norm()'s second value of a pair
  double         spare;
} rng_t;

void   rng_seed(rng_t *r, uint64_t seed, uint64_t stream);
double rng_unif(rng_t *r);
double rng_norm(rng_t *r);
double rng_gamma(rng_t *r, double shape);
double rng_log_gamma(rng_t *r, double shape);

/* --- phylo_structure.c: patristic distances among a sample's tips --- */
// A binary tree over tips 0 .. T - 1, by preorder number q = 0 .. N - 1.
typedef struct {
  int     T, N;     // tips; nodes (2T - 1)
  double *dep;      // q: distance from the root
  int    *size;     // q: nodes in its subtree (preorder span)
  int    *par;      // q: parent's preorder number (-1: root)
  double *len;      // q: length of the edge to its parent
  int    *tip;      // q: tip id, or -1
  int    *tpos;     // tip id: rank among the tips in preorder
  int    *tpre;     // rank: preorder number
  int    *rmq;      // sparse table: LCAs of tips adjacent in preorder
  int     levels;
  int    *lg;       // 0 .. T: floor(log2(x))
} ptree_t;

// Scratch for one sample's virtual tree (per thread).
typedef struct {
  int     cap;      // most tips
  int    *key, *node, *par, *c1, *stk, *leaf, *tq, *tp, *tpo;
  double *kw, *len, *w, *b1, *b2, *tw;
} ptree_ws_t;

int    ptree_init(ptree_t *pt, int T, const int *epar, const int *echd,
                  const double *elen);
void   ptree_ws_init(ptree_ws_t *ws, int cap);
int    ptree_vtree(const ptree_t *pt, ptree_ws_t *ws, int S, const int *tips,
                   const double *wt);
double ptree_pair_sum(ptree_ws_t *ws, int m);

/* --- rng_binomial.c: binomial and multinomial draws (NumPy's sampler) --- */
int64_t rng_binom(rng_t *rng, int64_t n, double p);
void    rng_multinom(rng_t *rng, int64_t n, const double *p, int d, int64_t *out);

/* --- rank_tests.c: rank tests shared by ancom.c and lefse.c --- */
typedef struct { double v; int g; } vg_t;   // a value and its group

void   rank_sort_vg(vg_t *a, int n);
double rank_kruskal_h(vg_t *vg, int n, const int *ng, int k, double *rsum);
double rank_sum_ties(vg_t *vg, int n, double *ties);

/* --- lda.c: the one-sided Jacobi SVD, shared with lefse.c --- */
int lda_jacobi_svd(double *A, int m, int c, double *V, double *s, double *tmp);

/* --- rarefaction_curve.c: closed forms of rarefaction, shared with hill_extrapolation --- */
int    rc_freq(const double *vals, int nnz, double *buf, double *k, double *f,
               double *N, double *S);
void   rc_absence(double N, double d, const double *U, int u, double *A);
double rc_coverage(double N, double d, int K, const double *k,
                   const double *f, const double *a);
double rc_coverage_full(double N, int K, const double *k, const double *f);

/* --- memory.c --- */
void  init_n_ptrs(int n);
void* safe_malloc(size_t bytes);
int   is_safe_ptr(void *ptr);
void  free_all(void);
void* free_one(void *ptr);
void* maybe_free_one(void *ptr);

/* --- normalize.c --- */
void normalize(ecomatrix_t *em, int norm, int n_threads, double pseudocount_);

/* --- da_tests.c: a normalized table by feature, for the DA tests --- */
typedef struct {
  int     n, D;     // Samples, features.
  int    *fpos;     // Feature f's stored values are fval[fpos[f] ..
  int    *fsam;     //   fpos[f + 1] - 1], of samples fsam[...] (in order).
  double *fval;
  double *z0;       // Each sample's value for a zero (NaN: missing).
} da_table_t;

void da_table_init(SEXP sexp_counts, SEXP sexp_margin, int norm, double pc,
                   int n_threads, da_table_t *t);
void da_table_column(const da_table_t *t, int f, double *y);

/* --- feature_lm.c: least squares for one design and many responses --- */
// A pivoted Householder QR of an n x p design, as LINPACK's dqrdc2 (which
// lm() uses): columns whose remaining norm falls below tol times their own
// are moved to the end and left out (aliased, coefficient NaN).
typedef struct {
  int     n, p, rank;
  double *a;        // n x p: R above the diagonal, Householder vectors below.
  double *qraux;    // p: each reflection's first element (0: none).
  int    *pivot;    // p: the design column in each position.
  double *rinv;     // p x p: the inverse of R (rank x rank, upper).
  double *q1;       // n x p: Q's first rank columns, or NULL.
  double *h;        // n: leverages, or NULL (with q1).
  double *work;     // n + p scratch.
} flm_t;

size_t flm_doubles(int n, int p, int want_q1);
void   flm_init(flm_t *q, double *mem, int *imem, int n, int p, int want_q1);
void   flm_decompose(flm_t *q, const double *X, int ldx, const int *rows,
                     int n, int p, double tol);
void   flm_fit(const flm_t *q, const double *y, double *qty, double *coef,
               double *rss);
void   flm_se(const flm_t *q, const double *y, const double *qty, double rss,
              int se_type, double *se, double *resid);

/* --- feature_lmm.c: REML random-intercept models for many responses --- */
typedef struct {
  int           n, p, J;
  const int    *grp;    // n: each sample's group, 0-based
  const double *nj;     // J: group sizes
  const double *xbar;   // J x p: group means of the design
  const double *xc;     // n x p: the design, centred within groups
  const double *wxx;    // p x p: xc' xc
} lmm_design_t;

void   lmm_design_init(lmm_design_t *d, const double *X, const int *grp,
                       int n, int p, int J);
size_t lmm_design_doubles(int n, int p, int J);
void   lmm_design_init_mem(lmm_design_t *d, const double *X, const int *grp,
                           int n, int p, int J, double *mem, long double *acc);
size_t lmm_doubles(int n, int p, int J);
double lmm_fit(const lmm_design_t *d, const double *y, double theta,
               double *mem, double *coef, double *vcov, double *df,
               double *sigma, double *dev);

/* --- rpca.c: shared with rpca_phylo.c --- */
void rpca_filter(int N0, int P0, const int *tpos, const int *tidx, const double *tval,
                 double msc, double mfc, double mff, int drop_single,
                 int *s_new, int *f_new, int *single, int *n_single, int *n, int *m);
SEXP rpca_core(int n, int m, int *rpos, int *ridx, double *rval,
               int r, int max_iter, double tol, int want_dm, int debug,
               double step, int n_threads, const char *no_var);
SEXP rpca_list_append(SEXP x, const char **names, SEXP *vals, int k);

/* --- pca.c: the centred (scaled) table, shared with rda(), plsda() and oplsda() --- */
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

/* --- plsda.c: PLS on pca.c's table, for plsda() and oplsda() --- */
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
