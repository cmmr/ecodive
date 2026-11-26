// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The decision trees and forests in this file are ported from scikit-learn
// 1.9.1: the depth-first tree builder and Tree (sklearn/tree/_tree.pyx),
// the best splitter (sklearn/tree/_splitter.pyx), the dense and sparse
// partitioners (sklearn/tree/_partitioner.pyx), the sort they use
// (sklearn/utils/_sorting.pyx), the Gini, entropy and squared-error
// criteria (sklearn/tree/_criterion.pyx), our_rand_r() and rand_int()
// (sklearn/utils/_random.pxd, sklearn/tree/_utils.pyx), and the forest's
// bootstrap, out-of-bag predictions and importances
// (sklearn/ensemble/_forest.py, _bootstrap.py, _base.py,
// sklearn/utils/class_weight.py). Those files are covered by this license:
//
//   BSD 3-Clause License
//
//   Copyright (c) 2007-2026 The scikit-learn developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The parity mode's random draws (used only by the tests, through the
// internal random_forest_details()) follow numpy's legacy RandomState
// (numpy/random/mtrand.pyx, numpy/random/src/mt19937/, and the masked
// bounded integers of numpy/random/src/distributions/), and np_sum()
// follows numpy's pairwise summation (pairwise_sum_DOUBLE(),
// numpy/_core/src/umath/loops_utils.h.src). Copyright (c) 2005-2025, NumPy
// Developers; BSD 3-Clause License, as above (see inst/COPYRIGHTS for the
// full text and for MT19937's own notice: Copyright (C) 1997 - 2002, Makoto
// Matsumoto and Takuji Nishimura; BSD 3-Clause).


/*
 * random_forest(): forests of CART trees, as scikit-learn grows them.
 *
 * Each tree is grown on a bootstrap sample of the training rows: the
 * draws' counts become the rows' weights, and rows drawn zero times are
 * left out of the tree. The tree is grown depth-first (left subtree
 * first), and each node is split at the best threshold over up to
 * `max_features` features, drawn without replacement (Fisher-Yates, also
 * keeping track of features known to be constant in the node). Values are
 * float32, as scikit-learn's.
 *
 * Randomness: every tree has its own stream, so a forest is the same for
 * any number of threads:
 *   - PCG (the default): stream `stream_base + t` of `seed`, for the
 *     bootstrap draws, then the tree's feature draws.
 *   - numpy (parity with scikit-learn, for the tests): tree t's seed s_t
 *     comes from MT19937(random_state).randint(2^31 - 1); the bootstrap is
 *     MT19937(s_t).randint(0, n, n) (or choice(p = )), and the feature
 *     draws are our_rand_r() seeded by MT19937(s_t).randint(0, 2^31 - 1).
 *
 * Entry points (all .Call):
 *   C_rf_grow        grow a forest on some rows of a table
 *   C_rf_predict     average the trees' leaf values for rows of a table,
 *                    over all trees or only those whose bootstrap missed
 *                    each row (out of bag)
 *   C_rf_importance  impurity (mean decrease in impurity) importances
 *   C_rf_permutation out-of-bag permutation importances
 *   C_rf_shuffle     a PCG Fisher-Yates shuffle (for the CV folds)
 *   C_rf_rng_test    numpy's and scikit-learn's draws, for the tests
 */

#include "ecodive.h"
#include <float.h>  // DBL_EPSILON
#include <limits.h> // INT_MAX

#define RF_LEAF   (-1)
#define RF_UNDEF  (-2)

static const float  RF_FEATURE_THRESHOLD = 1e-7f;  // float32, as sklearn's
static const float  RF_EXTRACT_NNZ_SWITCH = 0.1f;
#define RF_EPSILON DBL_EPSILON

#define CRIT_GINI    0
#define CRIT_ENTROPY 1
#define CRIT_MSE     2


//======================================================
// numpy's pairwise summation (np.sum of a 1-d double
// array): 0 + pairwise_sum(a, n).
//======================================================
static double np_pairwise (const double *a, int n) {
  if (n < 8) {
    double res = 0.;
    for (int i = 0; i < n; i++) res += a[i];
    return res;
  }
  else if (n <= 128) {
    double r[8];
    int i;
    for (int j = 0; j < 8; j++) r[j] = a[j];
    for (i = 8; i < n - (n % 8); i += 8)
      for (int j = 0; j < 8; j++) r[j] += a[i + j];
    double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; i++) res += a[i];
    return res;
  }
  int n2 = n / 2;
  n2 -= n2 % 8;
  return np_pairwise(a, n2) + np_pairwise(a + n2, n - n2);
}

static double np_sum (const double *a, int n) {
  return 0. + np_pairwise(a, n);
}


//======================================================
// MT19937 as numpy's legacy RandomState uses it.
//======================================================
typedef struct { uint32_t key[624]; int pos; } mt_t;

static void mt_seed (mt_t *m, uint32_t seed) {
  for (int pos = 0; pos < 624; pos++) {
    m->key[pos] = seed;
    seed = (uint32_t) (1812433253UL * (seed ^ (seed >> 30)) + (uint32_t) pos + 1);
  }
  m->pos = 624;
}

static void mt_gen (mt_t *m) {
  uint32_t y;
  int i;
  uint32_t *k = m->key;
  for (i = 0; i < 624 - 397; i++) {
    y = (k[i] & 0x80000000UL) | (k[i + 1] & 0x7fffffffUL);
    k[i] = k[i + 397] ^ (y >> 1) ^ (-(y & 1) & 0x9908b0dfUL);
  }
  for (; i < 623; i++) {
    y = (k[i] & 0x80000000UL) | (k[i + 1] & 0x7fffffffUL);
    k[i] = k[i + (397 - 624)] ^ (y >> 1) ^ (-(y & 1) & 0x9908b0dfUL);
  }
  y = (k[623] & 0x80000000UL) | (k[0] & 0x7fffffffUL);
  k[623] = k[396] ^ (y >> 1) ^ (-(y & 1) & 0x9908b0dfUL);
  m->pos = 0;
}

static uint32_t mt_next (mt_t *m) {
  if (m->pos == 624) mt_gen(m);
  uint32_t y = m->key[m->pos++];
  y ^= (y >> 11);
  y ^= (y << 7)  & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  return y;
}

// randint(0, n): the masked rejection method, for n <= 2^32.
static uint32_t mt_randint (mt_t *m, uint64_t n) {
  uint64_t rng = n - 1;
  if (rng == 0) return 0;
  if (rng == 0xFFFFFFFFULL) return mt_next(m);
  uint32_t mask = (uint32_t) rng;
  mask |= mask >> 1; mask |= mask >> 2; mask |= mask >> 4;
  mask |= mask >> 8; mask |= mask >> 16;
  uint32_t v;
  while ((v = (mt_next(m) & mask)) > (uint32_t) rng);
  return v;
}

// random_sample(): 53 bits in [0, 1).
static double mt_double (mt_t *m) {
  int32_t a = (int32_t) (mt_next(m) >> 5), b = (int32_t) (mt_next(m) >> 6);
  return (a * 67108864.0 + b) / 9007199254740992.0;
}


//======================================================
// scikit-learn's our_rand_r(): a 32-bit xorshift.
//======================================================
static inline uint32_t our_rand_r (uint32_t *seed) {
  if (seed[0] == 0) seed[0] = 1;
  seed[0] ^= (uint32_t) (seed[0] << 13);
  seed[0] ^= (uint32_t) (seed[0] >> 17);
  seed[0] ^= (uint32_t) (seed[0] << 5);
  return seed[0] % ((uint32_t) 2147483647 + 1);
}


//======================================================
// A tree's random draws: PCG, or numpy's.
//======================================================
typedef struct {
  int      numpy;
  uint32_t xs;    // numpy mode: our_rand_r()'s state
  mt_t     mt;    // numpy mode: the bootstrap's MT19937
  rng_t    pcg;   // PCG mode
} tree_rng_t;

static inline int rand_int (tree_rng_t *r, int low, int high) {
  if (r->numpy)
    return low + (int) (our_rand_r(&r->xs) % (uint32_t) (high - low));
  return low + (int) pcg32_bounded(&r->pcg.pcg, (uint32_t) (high - low));
}

static inline uint32_t boot_int (tree_rng_t *r, int n) {
  if (r->numpy) return mt_randint(&r->mt, (uint64_t) n);
  return pcg32_bounded(&r->pcg.pcg, (uint32_t) n);
}

static inline double boot_unif (tree_rng_t *r) {
  if (r->numpy) return mt_double(&r->mt);
  return rng_unif(&r->pcg);
}


//======================================================
// The sort (sklearn/utils/_sorting.pyx): introsort with
// median-of-3 pivots and 3-way partitioning; not stable.
//======================================================
static inline void sswap (float *v, int *ix, int i, int j) {
  float t = v[i]; v[i] = v[j]; v[j] = t;
  int   u = ix[i]; ix[i] = ix[j]; ix[j] = u;
}

static inline float median3 (const float *v, int n) {
  float a = v[0], b = v[n / 2], c = v[n - 1];
  if (a < b) {
    if (b < c)      return b;
    else if (a < c) return c;
    else            return a;
  }
  else if (b < c) {
    if (a < c) return a;
    else       return c;
  }
  else return b;
}

static void sift_down (float *v, int *ix, int start, int end) {
  int child, maxind, root = start;
  for (;;) {
    child  = root * 2 + 1;
    maxind = root;
    if (child < end && v[maxind] < v[child])         maxind = child;
    if (child + 1 < end && v[maxind] < v[child + 1]) maxind = child + 1;
    if (maxind == root) break;
    sswap(v, ix, root, maxind);
    root = maxind;
  }
}

static void heapsort_f (float *v, int *ix, int n) {
  int start = (n - 2) / 2, end = n;
  for (;;) {
    sift_down(v, ix, start, end);
    if (start == 0) break;
    start--;
  }
  end = n - 1;
  while (end > 0) {
    sswap(v, ix, 0, end);
    sift_down(v, ix, 0, end);
    end--;
  }
}

static void insertion_sort (float *v, int *ix, int n) {
  for (int i = 1; i < n; i++) {
    float tv = v[i];
    int   ti = ix[i];
    int   j  = i;
    while (j > 0 && v[j - 1] > tv) {
      v[j]  = v[j - 1];
      ix[j] = ix[j - 1];
      j--;
    }
    v[j]  = tv;
    ix[j] = ti;
  }
}

static void introsort_3way (float *v, int *ix, int n, int maxd) {
  while (n > 15) {
    if (maxd <= 0) { heapsort_f(v, ix, n); return; }
    maxd--;
    float pivot = median3(v, n);
    int i = 0, l = 0, r = n;
    while (i < r) {
      if (v[i] < pivot)      { sswap(v, ix, i, l); i++; l++; }
      else if (v[i] > pivot) { r--; sswap(v, ix, i, r); }
      else                   { i++; }
    }
    introsort_3way(v, ix, l, maxd);
    v  += r;
    ix += r;
    n  -= r;
  }
  insertion_sort(v, ix, n);
}

static void simultaneous_sort (float *v, int *ix, int n) {
  if (n == 0) return;
  int maxd = 2 * (int) log2((double) n);
  introsort_3way(v, ix, n, maxd);
}


//======================================================
// The table, as float32: dense (feature-major) or CSC.
//======================================================
typedef struct {
  int    n, p, sparse;
  float *X;                 // dense: X[f * n + i]
  int   *indptr, *indices;  // CSC
  float *data;
} rf_x_t;


//======================================================
// A grown tree.
//======================================================
typedef struct {
  int     n_nodes, cap, k;
  int    *left, *right, *feature, *n_samples;
  double *threshold, *impurity, *wn, *value;
} rf_tree_t;

static void tree_free (rf_tree_t *t) {
  free(t->left); free(t->right); free(t->feature); free(t->n_samples);
  free(t->threshold); free(t->impurity); free(t->wn); free(t->value);
  memset(t, 0, sizeof(rf_tree_t));
}

static int tree_grow_cap (rf_tree_t *t) {
  int cap = t->cap ? 2 * t->cap : 63;
  void *a;
#define RF_REALLOC(field, type, width) \
  a = realloc(t->field, (size_t) cap * (width) * sizeof(type)); \
  if (!a) return -1; \
  t->field = (type *) a;
  RF_REALLOC(left,      int,    1)
  RF_REALLOC(right,     int,    1)
  RF_REALLOC(feature,   int,    1)
  RF_REALLOC(n_samples, int,    1)
  RF_REALLOC(threshold, double, 1)
  RF_REALLOC(impurity,  double, 1)
  RF_REALLOC(wn,        double, 1)
  RF_REALLOC(value,     double, t->k)
#undef RF_REALLOC
  t->cap = cap;
  return 0;
}

static int tree_add_node (
    rf_tree_t *t, int parent, int is_left, int is_leaf, int feature,
    double threshold, double impurity, int n_node_samples, double wn ) {

  int id = t->n_nodes;
  if (id >= t->cap && tree_grow_cap(t)) return -1;
  t->impurity[id]  = impurity;
  t->n_samples[id] = n_node_samples;
  t->wn[id]        = wn;
  if (parent != RF_UNDEF) {
    if (is_left) t->left[parent]  = id;
    else         t->right[parent] = id;
  }
  if (is_leaf) {
    t->left[id]      = RF_LEAF;
    t->right[id]     = RF_LEAF;
    t->feature[id]   = RF_UNDEF;
    t->threshold[id] = RF_UNDEF;
  }
  else {
    t->feature[id]   = feature;
    t->threshold[id] = threshold;
  }
  t->n_nodes++;
  return id;
}


//======================================================
// The criterion (sklearn/tree/_criterion.pyx).
//======================================================
typedef struct {
  int           kind, k;
  const double *y;   // class code (as double) or value, per row
  const double *w;   // weight per row
  const int    *samples;
  int           start, end, pos;
  double        wn_samples, wn_node, wn_left, wn_right, sq_sum_total;
  double       *sum_total, *sum_left, *sum_right;
} crit_t;

static void crit_reset (crit_t *c) {
  c->pos = c->start;
  memset(c->sum_left, 0, (size_t) c->k * sizeof(double));
  memcpy(c->sum_right, c->sum_total, (size_t) c->k * sizeof(double));
  c->wn_left  = 0.0;
  c->wn_right = c->wn_node;
}

static void crit_reverse_reset (crit_t *c) {
  c->pos = c->end;
  memset(c->sum_right, 0, (size_t) c->k * sizeof(double));
  memcpy(c->sum_left, c->sum_total, (size_t) c->k * sizeof(double));
  c->wn_right = 0.0;
  c->wn_left  = c->wn_node;
}

static void crit_init (crit_t *c, int start, int end) {
  c->start   = start;
  c->end     = end;
  c->wn_node = 0.;
  memset(c->sum_total, 0, (size_t) c->k * sizeof(double));
  if (c->kind == CRIT_MSE) {
    c->sq_sum_total = 0.0;
    for (int p = start; p < end; p++) {
      int    i    = c->samples[p];
      double w    = c->w[i];
      double y_ik = c->y[i];
      double wy   = w * y_ik;
      c->sum_total[0] += wy;
      c->sq_sum_total += wy * y_ik;
      c->wn_node      += w;
    }
  }
  else {
    for (int p = start; p < end; p++) {
      int    i = c->samples[p];
      double w = c->w[i];
      c->sum_total[(int) c->y[i]] += w;
      c->wn_node += w;
    }
  }
  crit_reset(c);
}

static void crit_update (crit_t *c, int new_pos) {
  int pos = c->pos;
  int cls = c->kind != CRIT_MSE;
  if ((new_pos - pos) <= (c->end - new_pos)) {
    for (int p = pos; p < new_pos; p++) {
      int    i = c->samples[p];
      double w = c->w[i];
      if (cls) c->sum_left[(int) c->y[i]] += w;
      else     c->sum_left[0] += w * c->y[i];
      c->wn_left += w;
    }
  }
  else {
    crit_reverse_reset(c);
    for (int p = c->end - 1; p > new_pos - 1; p--) {
      int    i = c->samples[p];
      double w = c->w[i];
      if (cls) c->sum_left[(int) c->y[i]] -= w;
      else     c->sum_left[0] -= w * c->y[i];
      c->wn_left -= w;
    }
  }
  c->wn_right = c->wn_node - c->wn_left;
  for (int j = 0; j < c->k; j++)
    c->sum_right[j] = c->sum_total[j] - c->sum_left[j];
  c->pos = new_pos;
}

static inline double ln2_log (double x) { return log(x) / log(2.0); }

static double crit_node_impurity (const crit_t *c) {
  const double *s  = c->sum_total;
  double        wn = c->wn_node;
  if (c->kind == CRIT_GINI) {
    double sq = 0.0;
    for (int j = 0; j < c->k; j++) sq += s[j] * s[j];
    double gini = 0.0;
    gini += 1.0 - sq / (wn * wn);
    return gini / 1;
  }
  if (c->kind == CRIT_ENTROPY) {
    double e = 0.0;
    for (int j = 0; j < c->k; j++) {
      double cnt = s[j];
      if (cnt > 0.0) { cnt /= wn; e -= cnt * ln2_log(cnt); }
    }
    return e / 1;
  }
  double imp = c->sq_sum_total / wn;
  double m   = s[0] / wn;
  imp -= m * m;
  return imp / 1;
}

static void crit_children_impurity (const crit_t *c, double *il, double *ir) {
  double wl = c->wn_left, wr = c->wn_right;
  if (c->kind == CRIT_GINI) {
    double sl = 0.0, sr = 0.0;
    for (int j = 0; j < c->k; j++) {
      double cnt = c->sum_left[j];
      sl += cnt * cnt;
      cnt = c->sum_right[j];
      sr += cnt * cnt;
    }
    double gl = 0.0, gr = 0.0;
    gl += 1.0 - sl / (wl * wl);
    gr += 1.0 - sr / (wr * wr);
    *il = gl / 1;
    *ir = gr / 1;
    return;
  }
  if (c->kind == CRIT_ENTROPY) {
    double el = 0.0, er = 0.0;
    for (int j = 0; j < c->k; j++) {
      double cnt = c->sum_left[j];
      if (cnt > 0.0) { cnt /= wl; el -= cnt * ln2_log(cnt); }
      cnt = c->sum_right[j];
      if (cnt > 0.0) { cnt /= wr; er -= cnt * ln2_log(cnt); }
    }
    *il = el / 1;
    *ir = er / 1;
    return;
  }
  double sq_left = 0.0;
  for (int p = c->start; p < c->pos; p++) {
    int    i    = c->samples[p];
    double w    = c->w[i];
    double y_ik = c->y[i];
    sq_left += w * y_ik * y_ik;
  }
  double sq_right = c->sq_sum_total - sq_left;
  double ml = c->sum_left[0] / wl, mr = c->sum_right[0] / wr;
  double l = sq_left / wl, r = sq_right / wr;
  l -= ml * ml;
  r -= mr * mr;
  *il = l / 1;
  *ir = r / 1;
}

static double crit_proxy (const crit_t *c) {
  if (c->kind == CRIT_MSE) {
    double pl = 0.0, pr = 0.0;
    pl += c->sum_left[0] * c->sum_left[0];
    pr += c->sum_right[0] * c->sum_right[0];
    return pl / c->wn_left + pr / c->wn_right;
  }
  double il, ir;
  crit_children_impurity(c, &il, &ir);
  return (- c->wn_right * ir - c->wn_left * il);
}

static double crit_improvement (const crit_t *c, double parent, double il, double ir) {
  return ((c->wn_node / c->wn_samples) *
          (parent - (c->wn_right / c->wn_node * ir)
                  - (c->wn_left  / c->wn_node * il)));
}

static void crit_node_value (const crit_t *c, double *dest) {
  for (int j = 0; j < c->k; j++) dest[j] = c->sum_total[j] / c->wn_node;
}


//======================================================
// The splitter's state for one tree.
//======================================================
typedef struct {
  // settings
  int max_features, min_samples_leaf, min_samples_split, max_depth;
  // data
  const rf_x_t *x;
  // splitter
  int    n_samples;            // rows with non-zero weight
  int   *samples;
  float *fv;                   // feature values
  int   *features, *constant_features;
  int    start, end;
  // sparse partitioner
  int   *index_to_samples, *sorted_samples;
  int    is_samples_sorted, end_negative, start_positive;
  // criterion
  crit_t crit;
  // the builder's stack
  struct stack_rec { int start, end, depth, parent, is_left, n_const; double impurity; } *stack;
  int    stack_cap;
} split_t;

typedef struct {
  double impurity_left, impurity_right, threshold, improvement;
  int    pos, feature;
} split_rec_t;


// ---- dense partitioner ----

static void dense_sort (split_t *s, int f) {
  const float *col = s->x->X + (size_t) f * s->x->n;
  for (int i = s->start; i < s->end; i++) s->fv[i] = col[s->samples[i]];
  simultaneous_sort(s->fv + s->start, s->samples + s->start, s->end - s->start);
}

static void dense_next_p (split_t *s, int *p_prev, int *p) {
  int end = s->end;
  (*p)++;
  while (*p < end && s->fv[*p] <= (float) (s->fv[*p - 1] + RF_FEATURE_THRESHOLD))
    (*p)++;
  *p_prev = *p - 1;
}

static void dense_partition_final (split_t *s, int feature, double threshold) {
  const float *col = s->x->X + (size_t) feature * s->x->n;
  int ps = s->start, pe = s->end;
  int *samples = s->samples;
  while (ps < pe) {
    if (col[samples[ps]] <= threshold) ps++;
    else {
      pe--;
      int t = samples[ps]; samples[ps] = samples[pe]; samples[pe] = t;
    }
  }
}


// ---- sparse partitioner ----

static inline void sparse_swap (int *its, int *samples, int a, int b) {
  int t = samples[a]; samples[a] = samples[b]; samples[b] = t;
  its[samples[a]] = a;
  its[samples[b]] = b;
}

static int cmp_int (const void *a, const void *b) {
  return (*(const int *) a) - (*(const int *) b);
}

static inline void binary_search (
    const int *arr, int start, int end, int value, int *index, int *new_start ) {
  *index = -1;
  while (start < end) {
    int pivot = start + (end - start) / 2;
    if (arr[pivot] == value) { *index = pivot; start = pivot + 1; break; }
    if (arr[pivot] < value) start = pivot + 1;
    else                    end   = pivot;
  }
  *new_start = start;
}

static void extract_nnz (split_t *s, int feature) {
  const rf_x_t *x = s->x;
  int  ip_start   = x->indptr[feature], ip_end = x->indptr[feature + 1];
  int  n_indices  = ip_end - ip_start;
  int  n_samples  = s->end - s->start;
  int  start = s->start, end = s->end;
  int *its = s->index_to_samples, *samples = s->samples;
  const int   *Xi = x->indices;
  const float *Xd = x->data;
  int  en = start, sp = end;

  double lhs = (1 - s->is_samples_sorted) * (double) n_samples * log2((double) n_samples) +
               (double) n_samples * log2((double) n_indices);
  double rhs = (double) (float) (RF_EXTRACT_NNZ_SWITCH * (float) n_indices);

  if (lhs < rhs) {
    // binary search
    int *ss = s->sorted_samples;
    if (!s->is_samples_sorted) {
      memcpy(ss + start, samples + start, (size_t) n_samples * sizeof(int));
      qsort(ss + start, (size_t) n_samples, sizeof(int), cmp_int);
      s->is_samples_sorted = 1;
    }
    while (ip_start < ip_end && ss[start] > Xi[ip_start]) ip_start++;
    while (ip_start < ip_end && ss[end - 1] < Xi[ip_end - 1]) ip_end--;
    int p = start, k;
    while (p < end && ip_start < ip_end) {
      binary_search(Xi, ip_start, ip_end, ss[p], &k, &ip_start);
      if (k != -1) {
        if (Xd[k] > 0) {
          sp--;
          s->fv[sp] = Xd[k];
          sparse_swap(its, samples, its[Xi[k]], sp);
        }
        else if (Xd[k] < 0) {
          s->fv[en] = Xd[k];
          sparse_swap(its, samples, its[Xi[k]], en);
          en++;
        }
      }
      p++;
    }
  }
  else {
    for (int k = ip_start; k < ip_end; k++) {
      int idx = its[Xi[k]];
      if (start <= idx && idx < end) {
        if (Xd[k] > 0) {
          sp--;
          s->fv[sp] = Xd[k];
          sparse_swap(its, samples, its[Xi[k]], sp);
        }
        else if (Xd[k] < 0) {
          s->fv[en] = Xd[k];
          sparse_swap(its, samples, its[Xi[k]], en);
          en++;
        }
      }
    }
  }
  s->end_negative   = en;
  s->start_positive = sp;
}

static void sparse_sort (split_t *s, int feature) {
  int *its = s->index_to_samples, *samples = s->samples;
  extract_nnz(s, feature);
  simultaneous_sort(s->fv + s->start, samples + s->start, s->end_negative - s->start);
  if (s->start_positive < s->end)
    simultaneous_sort(s->fv + s->start_positive, samples + s->start_positive,
                      s->end - s->start_positive);
  for (int p = s->start; p < s->end_negative; p++) its[samples[p]] = p;
  for (int p = s->start_positive; p < s->end; p++) its[samples[p]] = p;
  if (s->end_negative < s->start_positive) {
    s->start_positive--;
    s->fv[s->start_positive] = 0.f;
    if (s->end_negative != s->start_positive) {
      s->fv[s->end_negative] = 0.f;
      s->end_negative++;
    }
  }
}

static void sparse_next_p (split_t *s, int *p_prev, int *p) {
  int p_next = (*p + 1 != s->end_negative) ? *p + 1 : s->start_positive;
  while (p_next < s->end &&
         s->fv[p_next] <= (float) (s->fv[*p] + RF_FEATURE_THRESHOLD)) {
    *p = p_next;
    p_next = (*p + 1 != s->end_negative) ? *p + 1 : s->start_positive;
  }
  *p_prev = *p;
  *p      = p_next;
}

static void sparse_partition_final (split_t *s, int feature, double threshold) {
  extract_nnz(s, feature);
  int p, pe;
  int *its = s->index_to_samples, *samples = s->samples;
  float *fv = s->fv;
  if (threshold < 0.)      { p = s->start;          pe = s->end_negative; }
  else if (threshold > 0.) { p = s->start_positive; pe = s->end; }
  else return;
  while (p < pe) {
    if (fv[p] <= threshold) p++;
    else {
      pe--;
      float t = fv[p]; fv[p] = fv[pe]; fv[pe] = t;
      sparse_swap(its, samples, p, pe);
    }
  }
}


//======================================================
// node_split_best() (sklearn/tree/_splitter.pyx).
//======================================================
static void node_split_best (
    split_t *s, tree_rng_t *rng, double impurity, int *n_constant_features,
    split_rec_t *best ) {

  int start = s->start, end = s->end;
  int *features = s->features, *constant_features = s->constant_features;
  int n_features = s->x->p;
  float *fv = s->fv;
  int sparse = s->x->sparse;
  crit_t *crit = &s->crit;

  double best_proxy = -INFINITY, current_proxy;
  int f_i = n_features, f_j, p, p_prev;
  int n_visited = 0, n_found_constants = 0, n_drawn_constants = 0;
  int n_known_constants = *n_constant_features;
  int n_total_constants = n_known_constants;
  split_rec_t cur;

  best->impurity_left  = INFINITY;
  best->impurity_right = INFINITY;
  best->pos            = end;
  best->feature        = 0;
  best->threshold      = 0.;
  best->improvement    = -INFINITY;
  cur = *best;

  if (sparse) s->is_samples_sorted = 0;

  while (f_i > n_total_constants &&
         (n_visited < s->max_features ||
          n_visited <= n_found_constants + n_drawn_constants)) {

    n_visited++;
    f_j = rand_int(rng, n_drawn_constants, f_i - n_found_constants);

    if (f_j < n_known_constants) {
      int t = features[n_drawn_constants];
      features[n_drawn_constants] = features[f_j];
      features[f_j] = t;
      n_drawn_constants++;
      continue;
    }

    f_j += n_found_constants;
    cur.feature = features[f_j];
    if (sparse) sparse_sort(s, cur.feature);
    else        dense_sort(s, cur.feature);

    if (end == start ||
        fv[end - 1] <= (float) (fv[start] + RF_FEATURE_THRESHOLD)) {
      int t = features[f_j];
      features[f_j] = features[n_total_constants];
      features[n_total_constants] = t;
      n_found_constants++;
      n_total_constants++;
      continue;
    }

    f_i--;
    { int t = features[f_i]; features[f_i] = features[f_j]; features[f_j] = t; }

    crit_reset(crit);
    p = start;
    while (p < end) {
      if (sparse) sparse_next_p(s, &p_prev, &p);
      else        dense_next_p(s, &p_prev, &p);
      if (p == end) continue;
      int n_left = p - start, n_right = end - p;
      if (n_left < s->min_samples_leaf || n_right < s->min_samples_leaf) continue;
      cur.pos = p;
      crit_update(crit, cur.pos);
      if (crit->wn_left < 0 || crit->wn_right < 0) continue;  // min_weight_leaf = 0
      current_proxy = crit_proxy(crit);
      if (current_proxy > best_proxy) {
        best_proxy = current_proxy;
        cur.threshold = fv[p_prev] / 2.0 + fv[p] / 2.0;
        *best = cur;
      }
    }
  }

  if (best->pos < end) {
    if (sparse) sparse_partition_final(s, best->feature, best->threshold);
    else        dense_partition_final(s, best->feature, best->threshold);
    crit_reset(crit);
    crit_update(crit, best->pos);
    crit_children_impurity(crit, &best->impurity_left, &best->impurity_right);
    best->improvement = crit_improvement(
      crit, impurity, best->impurity_left, best->impurity_right );
  }

  memcpy(features, constant_features, (size_t) n_known_constants * sizeof(int));
  memcpy(constant_features + n_known_constants, features + n_known_constants,
         (size_t) n_found_constants * sizeof(int));
  *n_constant_features = n_total_constants;
}


//======================================================
// DepthFirstTreeBuilder.build() (sklearn/tree/_tree.pyx).
//======================================================
static int build_tree (split_t *s, tree_rng_t *rng, rf_tree_t *tree) {

  crit_t *crit = &s->crit;
  int    n_stack = 0;
  int    first = 1;
  double parent_impurity = INFINITY;
  int    n_const = 0;
  split_rec_t split;
  memset(&split, 0, sizeof(split));

#define PUSH(a, b, c, d, e, f, g) do { \
    if (n_stack >= s->stack_cap) { \
      int nc = s->stack_cap ? 2 * s->stack_cap : 64; \
      void *q = realloc(s->stack, (size_t) nc * sizeof(*s->stack)); \
      if (!q) return -1; \
      s->stack = q; s->stack_cap = nc; \
    } \
    s->stack[n_stack].start = (a); s->stack[n_stack].end = (b); \
    s->stack[n_stack].depth = (c); s->stack[n_stack].parent = (d); \
    s->stack[n_stack].is_left = (e); s->stack[n_stack].impurity = (f); \
    s->stack[n_stack].n_const = (g); n_stack++; } while (0)

  PUSH(0, s->n_samples, 0, RF_UNDEF, 0, INFINITY, 0);

  while (n_stack > 0) {
    struct stack_rec r = s->stack[--n_stack];
    int start = r.start, end = r.end, depth = r.depth;
    parent_impurity = r.impurity;
    n_const = r.n_const;
    int n_node_samples = end - start;

    s->start = start;
    s->end   = end;
    crit_init(crit, start, end);
    double wn = crit->wn_node;

    int is_leaf = (depth >= s->max_depth ||
                   n_node_samples < s->min_samples_split ||
                   n_node_samples < 2 * s->min_samples_leaf ||
                   wn < 0);   // 2 * min_weight_leaf, which is 0
    if (first) {
      parent_impurity = crit_node_impurity(crit);
      first = 0;
    }
    is_leaf = is_leaf || parent_impurity <= RF_EPSILON;

    if (!is_leaf) {
      node_split_best(s, rng, parent_impurity, &n_const, &split);
      is_leaf = is_leaf || split.pos >= end ||
                (split.improvement + RF_EPSILON < 0.0);
    }

    int node_id = tree_add_node(
      tree, r.parent, r.is_left, is_leaf, split.feature, split.threshold,
      parent_impurity, n_node_samples, wn );
    if (node_id < 0) return -1;

    // The node's value; crit_node_value() reads sum_total, which
    // node_split_best() leaves as the node's.
    crit_node_value(crit, tree->value + (size_t) node_id * tree->k);

    if (!is_leaf) {
      PUSH(split.pos, end,       depth + 1, node_id, 0, split.impurity_right, n_const);
      PUSH(start,     split.pos, depth + 1, node_id, 1, split.impurity_left,  n_const);
    }
  }
#undef PUSH
  return 0;
}


//======================================================
// Reading the table into float32, for some rows.
//======================================================

// sexp_x: a dense double matrix (n_total x p), or list(i, p, x, n) with a
// CSC table (0-based row indices i, column pointers p, values x, n rows).
// rows: 0-based rows to keep, in the order wanted. Returns 0, or -1 if out
// of memory.
static int read_x (SEXP sexp_x, const int *rows, int n_rows, rf_x_t *x) {

  memset(x, 0, sizeof(rf_x_t));
  x->n = n_rows;

  if (isMatrix(sexp_x)) {
    int     n_total = nrows(sexp_x);
    int     p       = ncols(sexp_x);
    const double *v = REAL(sexp_x);
    x->p = p;
    x->X = (float *) malloc((size_t) n_rows * (size_t) (p > 0 ? p : 1) * sizeof(float));
    if (!x->X) return -1;
    for (int f = 0; f < p; f++) {
      const double *col = v + (size_t) f * n_total;
      float        *out = x->X + (size_t) f * n_rows;
      for (int i = 0; i < n_rows; i++) out[i] = (float) col[rows[i]];
    }
    return 0;
  }

  x->sparse = 1;
  const int    *ci  = INTEGER(VECTOR_ELT(sexp_x, 0));
  const int    *cp  = INTEGER(VECTOR_ELT(sexp_x, 1));
  const double *cx  = REAL(VECTOR_ELT(sexp_x, 2));
  int           n_total = asInteger(VECTOR_ELT(sexp_x, 3));
  int           p   = length(VECTOR_ELT(sexp_x, 1)) - 1;
  x->p = p;

  // Each row's new index, or -1. Rows must be in increasing order for
  // the column's indices to stay sorted; R passes them so.
  int *map = (int *) malloc((size_t) (n_total > 0 ? n_total : 1) * sizeof(int));
  x->indptr = (int *) malloc((size_t) (p + 1) * sizeof(int));
  if (!map || !x->indptr) { free(map); return -1; }
  for (int i = 0; i < n_total; i++) map[i] = -1;
  for (int i = 0; i < n_rows; i++) map[rows[i]] = i;

  int nnz = 0;
  for (int f = 0; f < p; f++)
    for (int e = cp[f]; e < cp[f + 1]; e++)
      if (map[ci[e]] >= 0 && cx[e] != 0) nnz++;

  x->indices = (int *)   malloc((size_t) (nnz > 0 ? nnz : 1) * sizeof(int));
  x->data    = (float *) malloc((size_t) (nnz > 0 ? nnz : 1) * sizeof(float));
  if (!x->indices || !x->data) { free(map); return -1; }

  nnz = 0;
  for (int f = 0; f < p; f++) {
    x->indptr[f] = nnz;
    for (int e = cp[f]; e < cp[f + 1]; e++)
      if (map[ci[e]] >= 0 && cx[e] != 0) {
        x->indices[nnz] = map[ci[e]];
        x->data[nnz]    = (float) cx[e];
        nnz++;
      }
  }
  x->indptr[p] = nnz;
  free(map);
  return 0;
}

static void free_x (rf_x_t *x) {
  free(x->X); free(x->indptr); free(x->indices); free(x->data);
  memset(x, 0, sizeof(rf_x_t));
}


//======================================================
// Growing a forest: one task per tree.
//======================================================
static struct {
  const rf_x_t *x;
  const double *y;          // per training row: class code or value
  const double *cw;         // per training row: class weight, or NULL
  int    k, kind, n_trees, n_boot, balanced_sub;
  int    max_features, min_samples_leaf, min_samples_split, max_depth;
  uint64_t seed, stream_base;
  const int *numpy_seeds;   // per tree, or NULL
  rf_tree_t *trees;
  int    *inbag;            // n x n_trees counts
  int     next_tree;
  int     failed;
} gf;

#if defined(__GNUC__) || defined(__clang__)
  #define RF_CLAIM(ctr) __atomic_fetch_add(&(ctr), 1, __ATOMIC_RELAXED)
  #define RF_FAIL()     __atomic_store_n(&gf.failed, 1, __ATOMIC_RELAXED)
  #define RF_PFAIL()    __atomic_store_n(&gp.failed, 1, __ATOMIC_RELAXED)
#else
  #define RF_CLAIM(ctr) (-1)
  #define RF_FAIL()     (gf.failed = 1)
  #define RF_PFAIL()    (gp.failed = 1)
#endif

static void grow_one (int t, split_t *s, double *w, double *cdf) {

  int n = gf.x->n;
  tree_rng_t rng;
  memset(&rng, 0, sizeof(rng));

  if (gf.numpy_seeds) {
    rng.numpy = 1;
    uint32_t ts = (uint32_t) gf.numpy_seeds[t];
    mt_t m2;
    mt_seed(&m2, ts);
    rng.xs = mt_randint(&m2, 2147483647ULL);   // randint(0, RAND_R_MAX)
    mt_seed(&rng.mt, ts);
  }
  else {
    rng_seed(&rng.pcg, gf.seed, gf.stream_base + (uint64_t) t);
  }

  // The bootstrap: counts of draws, as the rows' weights.
  int *cnt = gf.inbag + (size_t) t * n;
  memset(cnt, 0, (size_t) n * sizeof(int));
  if (gf.cw) {
    // choice(n, n_boot, p = cw / sum(cw)): cumsum, normalised, searchsorted
    // (side = 'right') of uniform draws.
    double tot = np_sum(gf.cw, n);
    double acc = 0;
    for (int i = 0; i < n; i++) { acc += gf.cw[i] / tot; cdf[i] = acc; }
    double last = cdf[n - 1];
    for (int i = 0; i < n; i++) cdf[i] /= last;
    for (int b = 0; b < gf.n_boot; b++) {
      double u = boot_unif(&rng);
      int lo = 0, hi = n;   // first index with cdf > u
      while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (cdf[mid] > u) hi = mid; else lo = mid + 1;
      }
      if (lo >= n) lo = n - 1;  // # nocov (u < 1 = cdf[n - 1])
      cnt[lo]++;
    }
  }
  else {
    for (int b = 0; b < gf.n_boot; b++) cnt[boot_int(&rng, n)]++;
  }

  for (int i = 0; i < n; i++) w[i] = (double) cnt[i];

  if (gf.balanced_sub) {
    // compute_sample_weight('balanced', y, indices = draws).
    double *cc = cdf;   // k class counts
    memset(cc, 0, (size_t) gf.k * sizeof(double));
    for (int i = 0; i < n; i++) cc[(int) gf.y[i]] += cnt[i];
    int k_sub = 0;
    for (int c = 0; c < gf.k; c++) if (cc[c] > 0) k_sub++;
    for (int i = 0; i < n; i++) {
      double c = cc[(int) gf.y[i]];
      double cw = c > 0 ? (double) gf.n_boot / ((double) k_sub * c) : 0.0;
      w[i] = (double) cnt[i] * cw;
    }
  }

  // Splitter.init(): the rows with non-zero weight, in order.
  double wsum = 0.0;
  int j = 0;
  for (int i = 0; i < n; i++) {
    if (w[i] != 0.0) s->samples[j++] = i;
    wsum += w[i];
  }
  s->n_samples = j;
  s->crit.wn_samples = wsum;
  s->crit.w = w;
  for (int f = 0; f < gf.x->p; f++) s->features[f] = f;
  if (gf.x->sparse) {
    for (int i = 0; i < n; i++) s->index_to_samples[i] = -1;
    for (int q = 0; q < j; q++) s->index_to_samples[s->samples[q]] = q;
  }

  rf_tree_t *tree = gf.trees + t;
  tree->k = gf.k;
  if (build_tree(s, &rng, tree)) RF_FAIL();
}

static void *grow_worker (void *arg) {
  int i  = ((worker_t *) arg)->i;
  int nt = ((worker_t *) arg)->n;
  int n  = gf.x->n, p = gf.x->p, k = gf.k;

  split_t s;
  memset(&s, 0, sizeof(s));
  s.x                 = gf.x;
  s.max_features      = gf.max_features;
  s.min_samples_leaf  = gf.min_samples_leaf;
  s.min_samples_split = gf.min_samples_split;
  s.max_depth         = gf.max_depth;
  s.crit.kind         = gf.kind;
  s.crit.k            = k;
  s.crit.y            = gf.y;

  size_t nn = (size_t) (n > 0 ? n : 1), pp = (size_t) (p > 0 ? p : 1);
  s.samples           = (int *)   malloc(nn * sizeof(int));
  s.fv                = (float *) malloc(nn * sizeof(float));
  s.features          = (int *)   malloc(pp * sizeof(int));
  s.constant_features = (int *)   malloc(pp * sizeof(int));
  s.index_to_samples  = (int *)   malloc(nn * sizeof(int));
  s.sorted_samples    = (int *)   malloc(nn * sizeof(int));
  s.crit.sum_total    = (double *) malloc((size_t) k * sizeof(double));
  s.crit.sum_left     = (double *) malloc((size_t) k * sizeof(double));
  s.crit.sum_right    = (double *) malloc((size_t) k * sizeof(double));
  s.crit.samples      = s.samples;
  double *w   = (double *) malloc(nn * sizeof(double));
  double *cdf = (double *) malloc((nn > (size_t) k ? nn : (size_t) k) * sizeof(double));

  if (!s.samples || !s.fv || !s.features || !s.constant_features ||
      !s.index_to_samples || !s.sorted_samples || !s.crit.sum_total ||
      !s.crit.sum_left || !s.crit.sum_right || !w || !cdf) {
    RF_FAIL();
  }
  else {
#if defined(__GNUC__) || defined(__clang__)
    (void) nt; (void) i;
    for (;;) {
      int t = RF_CLAIM(gf.next_tree);
      if (t >= gf.n_trees || __atomic_load_n(&gf.failed, __ATOMIC_RELAXED)) break;
      grow_one(t, &s, w, cdf);
    }
#else
    for (int t = i; t < gf.n_trees; t += nt) grow_one(t, &s, w, cdf);
#endif
  }

  free(s.samples); free(s.fv); free(s.features); free(s.constant_features);
  free(s.index_to_samples); free(s.sorted_samples);
  free(s.crit.sum_total); free(s.crit.sum_left); free(s.crit.sum_right);
  free(s.stack); free(w); free(cdf);
  return NULL;
}


// The forest as an R list.
static SEXP forest_sexp (rf_tree_t *trees, int n_trees, int k, SEXP inbag) {

  int total = 0;
  for (int t = 0; t < n_trees; t++) total += trees[t].n_nodes;

  const char *names[] = {
    "offsets", "left", "right", "feature", "threshold", "value",
    "impurity", "n_node_samples", "weighted_n_node_samples", "inbag", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  SEXP off = allocVector(INTSXP,  n_trees + 1);        SET_VECTOR_ELT(res, 0, off);
  SEXP lft = allocVector(INTSXP,  total);              SET_VECTOR_ELT(res, 1, lft);
  SEXP rgt = allocVector(INTSXP,  total);              SET_VECTOR_ELT(res, 2, rgt);
  SEXP fea = allocVector(INTSXP,  total);              SET_VECTOR_ELT(res, 3, fea);
  SEXP thr = allocVector(REALSXP, total);              SET_VECTOR_ELT(res, 4, thr);
  SEXP val = allocMatrix(REALSXP, k, total);           SET_VECTOR_ELT(res, 5, val);
  SEXP imp = allocVector(REALSXP, total);              SET_VECTOR_ELT(res, 6, imp);
  SEXP nns = allocVector(INTSXP,  total);              SET_VECTOR_ELT(res, 7, nns);
  SEXP wnn = allocVector(REALSXP, total);              SET_VECTOR_ELT(res, 8, wnn);
  SET_VECTOR_ELT(res, 9, inbag);

  int o = 0;
  for (int t = 0; t < n_trees; t++) {
    rf_tree_t *tr = trees + t;
    int m = tr->n_nodes;
    INTEGER(off)[t] = o;
    if (m > 0) {
      memcpy(INTEGER(lft) + o, tr->left,      (size_t) m * sizeof(int));
      memcpy(INTEGER(rgt) + o, tr->right,     (size_t) m * sizeof(int));
      memcpy(INTEGER(fea) + o, tr->feature,   (size_t) m * sizeof(int));
      memcpy(INTEGER(nns) + o, tr->n_samples, (size_t) m * sizeof(int));
      memcpy(REAL(thr) + o,    tr->threshold, (size_t) m * sizeof(double));
      memcpy(REAL(imp) + o,    tr->impurity,  (size_t) m * sizeof(double));
      memcpy(REAL(wnn) + o,    tr->wn,        (size_t) m * sizeof(double));
      memcpy(REAL(val) + (size_t) o * k, tr->value, (size_t) m * k * sizeof(double));
    }
    o += m;
  }
  INTEGER(off)[n_trees] = o;

  UNPROTECT(1);
  return res;
}


//======================================================
// C_rf_grow(x, rows, y, k, params, class_w, seed,
//           stream_base, numpy_seeds, n_threads)
//
//   x            dense double matrix, or CSC list (see read_x())
//   rows         0-based training rows, increasing
//   y            per training row: class code 0 .. k-1, or value
//   k            the number of classes; 0 for regression
//   params       integer: n_trees, max_features, min_samples_leaf,
//                min_samples_split, max_depth (NA: none), criterion
//                (0 gini, 1 entropy, 2 mse), balanced_subsample (0/1)
//   max_samples  NULL, or c(type, value): type 1 a count, 2 a fraction
//   class_w      NULL, or per training row: class weight (drives a
//                weighted bootstrap)
//   seed         integer
//   stream_base  double: PCG stream of tree 0
//   numpy_seeds  NULL, or integer per tree (parity mode)
//======================================================
SEXP C_rf_grow (
    SEXP sexp_x, SEXP sexp_rows, SEXP sexp_y, SEXP sexp_k, SEXP sexp_params,
    SEXP sexp_max_samples, SEXP sexp_class_w, SEXP sexp_seed,
    SEXP sexp_stream_base, SEXP sexp_numpy_seeds, SEXP sexp_n_threads ) {

  int  n_rows = length(sexp_rows);
  int *par    = INTEGER(sexp_params);
  int  k      = asInteger(sexp_k);

  rf_x_t x;
  if (read_x(sexp_x, INTEGER(sexp_rows), n_rows, &x)) {
    free_x(&x);
    error("Unable to allocate memory for the forest's table.");
  }

  memset(&gf, 0, sizeof(gf));
  gf.x                 = &x;
  gf.y                 = REAL(sexp_y);
  gf.k                 = k > 0 ? k : 1;
  gf.n_trees           = par[0];
  gf.max_features      = par[1];
  gf.min_samples_leaf  = par[2];
  gf.min_samples_split = par[3] > 2 * par[2] ? par[3] : 2 * par[2];
  gf.max_depth         = par[4] == NA_INTEGER ? INT_MAX : par[4];
  gf.kind              = par[5];
  gf.balanced_sub      = par[6];
  gf.cw                = isNull(sexp_class_w) ? NULL : REAL(sexp_class_w);
  gf.seed              = (uint64_t) asInteger(sexp_seed);
  gf.stream_base       = (uint64_t) asReal(sexp_stream_base);
  gf.numpy_seeds       = isNull(sexp_numpy_seeds) ? NULL : INTEGER(sexp_numpy_seeds);

  // _get_n_samples_bootstrap()
  gf.n_boot = n_rows;
  if (!isNull(sexp_max_samples)) {
    double type = REAL(sexp_max_samples)[0], v = REAL(sexp_max_samples)[1];
    if (type == 1) gf.n_boot = (int) v;
    else {
      double wsum = gf.cw ? np_sum(gf.cw, n_rows) : (double) n_rows;
      double nb   = v * wsum;
      gf.n_boot = nb >= 1 ? (int) nb : 1;
    }
  }

  SEXP inbag = PROTECT(allocMatrix(INTSXP, n_rows, gf.n_trees));
  gf.inbag = INTEGER(inbag);
  gf.trees = (rf_tree_t *) calloc((size_t) (gf.n_trees > 0 ? gf.n_trees : 1), sizeof(rf_tree_t));
  if (!gf.trees) { free_x(&x); error("Unable to allocate memory for the forest."); }

  int n_threads = asInteger(sexp_n_threads);
  if (n_threads > gf.n_trees) n_threads = gf.n_trees;
  if (n_threads < 1) n_threads = 1;
  if (n_rows > 0 && gf.n_trees > 0)
    run_parallel(grow_worker, n_threads, n_threads > 1 ? 100 : 1);

  int failed = gf.failed;
  SEXP res = R_NilValue;
  if (!failed) res = PROTECT(forest_sexp(gf.trees, gf.n_trees, gf.k, inbag));

  for (int t = 0; t < gf.n_trees; t++) tree_free(gf.trees + t);
  free(gf.trees);
  free_x(&x);
  memset(&gf, 0, sizeof(gf));

  if (failed) error("Unable to allocate memory for the forest's trees.");
  UNPROTECT(2);
  return res;
}


//======================================================
// Row access for prediction: dense columns, or each row's
// sorted (feature, value) list.
//======================================================
typedef struct {
  int    n, p, sparse;
  const double *X; int n_total;  // dense: X[f * n_total + rows[i]]
  const int *rows;
  int   *rptr, *rcol;            // sparse: CSR of the rows
  float *rval;
} rf_rows_t;

static int read_rows (SEXP sexp_x, const int *rows, int n_rows, rf_rows_t *r) {

  memset(r, 0, sizeof(rf_rows_t));
  r->n    = n_rows;
  r->rows = rows;

  if (isMatrix(sexp_x)) {
    r->X       = REAL(sexp_x);
    r->n_total = nrows(sexp_x);
    r->p       = ncols(sexp_x);
    return 0;
  }

  r->sparse = 1;
  const int    *ci = INTEGER(VECTOR_ELT(sexp_x, 0));
  const int    *cp = INTEGER(VECTOR_ELT(sexp_x, 1));
  const double *cx = REAL(VECTOR_ELT(sexp_x, 2));
  int n_total = asInteger(VECTOR_ELT(sexp_x, 3));
  int p       = length(VECTOR_ELT(sexp_x, 1)) - 1;
  r->p = p;

  int *map = (int *) malloc((size_t) (n_total > 0 ? n_total : 1) * sizeof(int));
  r->rptr  = (int *) calloc((size_t) n_rows + 1, sizeof(int));
  if (!map || !r->rptr) { free(map); return -1; }
  for (int i = 0; i < n_total; i++) map[i] = -1;
  // A row may be asked for more than once (not by R; kept simple):
  // only its last position is filled, and the others share it below.
  for (int i = 0; i < n_rows; i++) map[rows[i]] = i;

  for (int f = 0; f < p; f++)
    for (int e = cp[f]; e < cp[f + 1]; e++)
      if (map[ci[e]] >= 0 && cx[e] != 0) r->rptr[map[ci[e]] + 1]++;
  for (int i = 0; i < n_rows; i++) r->rptr[i + 1] += r->rptr[i];
  int nnz = r->rptr[n_rows];

  r->rcol = (int *)   malloc((size_t) (nnz > 0 ? nnz : 1) * sizeof(int));
  r->rval = (float *) malloc((size_t) (nnz > 0 ? nnz : 1) * sizeof(float));
  int *fill = (int *) malloc((size_t) (n_rows > 0 ? n_rows : 1) * sizeof(int));
  if (!r->rcol || !r->rval || !fill) { free(map); free(fill); return -1; }
  memcpy(fill, r->rptr, (size_t) n_rows * sizeof(int));

  // Columns in order, so each row's features are sorted.
  for (int f = 0; f < p; f++)
    for (int e = cp[f]; e < cp[f + 1]; e++) {
      int i = map[ci[e]];
      if (i >= 0 && cx[e] != 0) {
        r->rcol[fill[i]] = f;
        r->rval[fill[i]] = (float) cx[e];
        fill[i]++;
      }
    }
  free(map);
  free(fill);
  return 0;
}

static void free_rows (rf_rows_t *r) {
  free(r->rptr); free(r->rcol); free(r->rval);
  memset(r, 0, sizeof(rf_rows_t));
}

static inline float row_value (const rf_rows_t *r, int i, int f) {
  if (!r->sparse) return (float) r->X[(size_t) f * r->n_total + r->rows[i]];
  int lo = r->rptr[i], hi = r->rptr[i + 1];
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    int c = r->rcol[mid];
    if (c == f) return r->rval[mid];
    if (c < f) lo = mid + 1; else hi = mid;
  }
  return 0.f;
}


// The forest from R.
typedef struct {
  int n_trees, k;
  const int *off, *left, *right, *feature;
  const double *threshold, *value;
  const int *inbag; int n_inbag;
} rf_forest_t;

static void read_forest (SEXP f, rf_forest_t *F) {
  SEXP off = VECTOR_ELT(f, 0);
  F->n_trees   = length(off) - 1;
  F->off       = INTEGER(off);
  F->left      = INTEGER(VECTOR_ELT(f, 1));
  F->right     = INTEGER(VECTOR_ELT(f, 2));
  F->feature   = INTEGER(VECTOR_ELT(f, 3));
  F->threshold = REAL(VECTOR_ELT(f, 4));
  F->value     = REAL(VECTOR_ELT(f, 5));
  F->k         = nrows(VECTOR_ELT(f, 5));
  SEXP ib      = VECTOR_ELT(f, 9);
  F->inbag     = isNull(ib) ? NULL : INTEGER(ib);
  F->n_inbag   = isNull(ib) ? 0 : nrows(ib);
}

// The leaf of row i in tree t (node index relative to the tree).
static inline int tree_leaf (const rf_forest_t *F, int t, const rf_rows_t *r, int i) {
  int o = F->off[t], node = 0;
  while (F->left[o + node] != RF_LEAF) {
    float v = row_value(r, i, F->feature[o + node]);
    node = (v <= F->threshold[o + node]) ? F->left[o + node] : F->right[o + node];
  }
  return node;
}


//======================================================
// C_rf_predict(forest, x, rows, oob, n_threads): the
// mean of the trees' leaf values (k per row), summed in
// tree order. With oob, row i is row i of the forest's
// training set, and only trees whose bootstrap missed it
// count; the second element gives how many did.
//======================================================
static struct {
  const rf_forest_t *F;
  const rf_rows_t   *R;
  int     oob;
  double *out;
  int    *n_used;
} pf;

static void *predict_worker (void *arg) {
  int i0 = ((worker_t *) arg)->i, nt = ((worker_t *) arg)->n;
  const rf_forest_t *F = pf.F;
  int k = F->k, n = pf.R->n;
  for (int i = i0; i < n; i += nt) {
    double *o = pf.out + (size_t) i * k;
    int used = 0;
    for (int j = 0; j < k; j++) o[j] = 0.0;
    for (int t = 0; t < F->n_trees; t++) {
      if (pf.oob && F->inbag[(size_t) t * F->n_inbag + i] != 0) continue;
      int leaf = tree_leaf(F, t, pf.R, i);
      const double *v = F->value + (size_t) (F->off[t] + leaf) * k;
      for (int j = 0; j < k; j++) o[j] += v[j];
      used++;
    }
    pf.n_used[i] = used;
    if (used > 0) for (int j = 0; j < k; j++) o[j] /= used;
  }
  return NULL;
}

SEXP C_rf_predict (
    SEXP sexp_forest, SEXP sexp_x, SEXP sexp_rows, SEXP sexp_oob,
    SEXP sexp_n_threads ) {

  rf_forest_t F;
  read_forest(sexp_forest, &F);
  int n_rows = length(sexp_rows);
  int oob    = asLogical(sexp_oob);
  if (oob && (F.inbag == NULL || F.n_inbag != n_rows))
    error("Out-of-bag predictions need the forest's training rows."); // # nocov

  rf_rows_t R;
  if (read_rows(sexp_x, INTEGER(sexp_rows), n_rows, &R)) {
    free_rows(&R);
    error("Unable to allocate memory for predictions.");
  }

  const char *names[] = { "value", "n_trees", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  SEXP val = allocMatrix(REALSXP, F.k, n_rows);  SET_VECTOR_ELT(res, 0, val);
  SEXP nu  = allocVector(INTSXP, n_rows);        SET_VECTOR_ELT(res, 1, nu);

  pf.F = &F; pf.R = &R; pf.oob = oob; pf.out = REAL(val); pf.n_used = INTEGER(nu);
  int n_threads = asInteger(sexp_n_threads);
  n_threads = threads_for(n_threads, (double) n_rows * F.n_trees, 2e4);
  run_parallel(predict_worker, n_threads, n_threads > 1 ? 100 : 1);

  free_rows(&R);
  UNPROTECT(1);
  return res;
}


//======================================================
// C_rf_importance(forest, p): impurity importances,
// as BaseForest.feature_importances_.
//======================================================
SEXP C_rf_importance (SEXP sexp_forest, SEXP sexp_p) {

  rf_forest_t F;
  read_forest(sexp_forest, &F);
  int p = asInteger(sexp_p);
  const double *imp = REAL(VECTOR_ELT(sexp_forest, 6));
  const double *wn  = REAL(VECTOR_ELT(sexp_forest, 8));

  SEXP res = PROTECT(allocVector(REALSXP, p));
  double *out = REAL(res);
  double *ti  = (double *) R_alloc((size_t) (p > 0 ? p : 1), sizeof(double));
  memset(out, 0, (size_t) p * sizeof(double));

  int n_used = 0;
  for (int t = 0; t < F.n_trees; t++) {
    int o = F.off[t], m = F.off[t + 1] - o;
    if (m <= 1) continue;
    memset(ti, 0, (size_t) p * sizeof(double));
    for (int a = 0; a < m; a++) {
      int l = F.left[o + a];
      if (l == RF_LEAF) continue;
      int r = F.right[o + a];
      ti[F.feature[o + a]] +=
        wn[o + a] * imp[o + a] - wn[o + l] * imp[o + l] - wn[o + r] * imp[o + r];
    }
    for (int j = 0; j < p; j++) ti[j] /= wn[o];
    double norm = np_sum(ti, p);
    if (norm > 0.0) for (int j = 0; j < p; j++) ti[j] /= norm;
    for (int j = 0; j < p; j++) out[j] += ti[j];
    n_used++;
  }
  if (n_used > 0) {
    for (int j = 0; j < p; j++) out[j] /= n_used;
    double tot = np_sum(out, p);
    for (int j = 0; j < p; j++) out[j] /= tot;
  }

  UNPROTECT(1);
  return res;
}


//======================================================
// C_rf_permutation(forest, x, rows, y, k, seed,
// n_threads): out-of-bag permutation importance. For each
// tree and each feature it splits on, the feature's
// values are permuted among the tree's out-of-bag rows
// (PCG stream 2^61 + t * p + j), and the increase in the
// tree's out-of-bag error (misclassification rate, or mean
// squared error) is recorded. Features a tree doesn't use
// get 0. Returns the mean and SD over trees, per feature.
//======================================================
static struct {
  const rf_forest_t *F;
  const rf_rows_t   *R;
  const double *y;
  int     k, cls, p;
  uint64_t seed;
  double **d;      // per tree: increase per used feature
  int    **u;      // per tree: the used features (sorted), -1 ended
  int      next_tree, failed;
} gp;

static inline double leaf_loss (const double *v, int k, int cls, double y) {
  if (cls) {
    int best = 0;
    for (int j = 1; j < k; j++) if (v[j] > v[best]) best = j;
    return best != (int) y;
  }
  double e = v[0] - y;
  return e * e;
}

static void perm_one (int t, int *pos, int *col_of, int *idx) {

  const rf_forest_t *F = gp.F;
  int o = F->off[t], m = F->off[t + 1] - o, n = gp.R->n, k = gp.k;

  // The features this tree splits on.
  int nu = 0;
  for (int a = 0; a < m; a++) {
    int f = F->feature[o + a];
    if (F->left[o + a] != RF_LEAF && col_of[f] < 0) col_of[f] = nu++;
  }
  int *used = (int *) malloc((size_t) (nu + 1) * sizeof(int));
  double *d = (double *) calloc((size_t) (nu + 1), sizeof(double));
  if (!used || !d) { free(used); free(d); RF_PFAIL(); return; }
  nu = 0;
  for (int f = 0; f < gp.p; f++) if (col_of[f] >= 0) used[nu++] = f;
  for (int a = 0; a < nu; a++) col_of[used[a]] = a;    // sorted columns
  used[nu] = -1;

  // Out-of-bag rows.
  int n_oob = 0;
  for (int i = 0; i < n; i++)
    if (F->inbag[(size_t) t * F->n_inbag + i] == 0) pos[n_oob++] = i;

  float *Zt = NULL;   // n_oob x nu, row-major: the used features' values
  if (n_oob > 0 && nu > 0)
    Zt = (float *) malloc((size_t) n_oob * (size_t) nu * sizeof(float));
  if (n_oob > 0 && nu > 0 && !Zt) RF_PFAIL();

  if (Zt) {

    for (int r = 0; r < n_oob; r++)
      for (int a = 0; a < nu; a++)
        Zt[(size_t) r * nu + a] = row_value(gp.R, pos[r], used[a]);

    // Base error.
    double base = 0;
    for (int r = 0; r < n_oob; r++) {
      const float *z = Zt + (size_t) r * nu;
      int node = 0;
      while (F->left[o + node] != RF_LEAF) {
        float v = z[col_of[F->feature[o + node]]];
        node = (v <= F->threshold[o + node]) ? F->left[o + node] : F->right[o + node];
      }
      base += leaf_loss(F->value + (size_t) (o + node) * k, k, gp.cls, gp.y[pos[r]]);
    }
    base /= n_oob;

    for (int a = 0; a < nu; a++) {
      rng_t rng;
      rng_seed(&rng, gp.seed,
               ((uint64_t) 1 << 61) + (uint64_t) t * (uint64_t) gp.p + (uint64_t) used[a]);
      for (int r = 0; r < n_oob; r++) idx[r] = r;
      for (int r = n_oob - 1; r > 0; r--) {
        int b = (int) pcg32_bounded(&rng.pcg, (uint32_t) (r + 1));
        int tmp = idx[r]; idx[r] = idx[b]; idx[b] = tmp;
      }
      double err = 0;
      for (int r = 0; r < n_oob; r++) {
        const float *z = Zt + (size_t) r * nu;
        float zp = Zt[(size_t) idx[r] * nu + a];
        int node = 0;
        while (F->left[o + node] != RF_LEAF) {
          int   c = col_of[F->feature[o + node]];
          float v = c == a ? zp : z[c];
          node = (v <= F->threshold[o + node]) ? F->left[o + node] : F->right[o + node];
        }
        err += leaf_loss(F->value + (size_t) (o + node) * k, k, gp.cls, gp.y[pos[r]]);
      }
      d[a] = err / n_oob - base;
    }
  }

  free(Zt);
  for (int a = 0; a < nu; a++) col_of[used[a]] = -1;
  gp.d[t] = d;
  gp.u[t] = used;
}

static void *perm_worker (void *arg) {
  int i0 = ((worker_t *) arg)->i, nt = ((worker_t *) arg)->n;
  const rf_forest_t *F = gp.F;
  int n = gp.R->n;

  int *pos    = (int *) malloc((size_t) (n > 0 ? n : 1) * sizeof(int));
  int *idx    = (int *) malloc((size_t) (n > 0 ? n : 1) * sizeof(int));
  int *col_of = (int *) malloc((size_t) (gp.p > 0 ? gp.p : 1) * sizeof(int));
  if (!pos || !idx || !col_of) RF_PFAIL();
  else {
    for (int f = 0; f < gp.p; f++) col_of[f] = -1;
#if defined(__GNUC__) || defined(__clang__)
    (void) i0; (void) nt;
    for (;;) {
      int t = RF_CLAIM(gp.next_tree);
      if (t >= F->n_trees || __atomic_load_n(&gp.failed, __ATOMIC_RELAXED)) break;
      perm_one(t, pos, col_of, idx);
    }
#else
    for (int t = i0; t < F->n_trees; t += nt) perm_one(t, pos, col_of, idx);
#endif
  }
  free(pos); free(idx); free(col_of);
  return NULL;
}

SEXP C_rf_permutation (
    SEXP sexp_forest, SEXP sexp_x, SEXP sexp_rows, SEXP sexp_y, SEXP sexp_k,
    SEXP sexp_seed, SEXP sexp_n_threads ) {

  rf_forest_t F;
  read_forest(sexp_forest, &F);
  int n_rows = length(sexp_rows);
  if (F.inbag == NULL || F.n_inbag != n_rows)
    error("Permutation importance needs the forest's training rows."); // # nocov

  rf_rows_t R;
  if (read_rows(sexp_x, INTEGER(sexp_rows), n_rows, &R)) {
    free_rows(&R);
    error("Unable to allocate memory for permutation importance.");
  }

  memset(&gp, 0, sizeof(gp));
  gp.F    = &F;
  gp.R    = &R;
  gp.y    = REAL(sexp_y);
  gp.k    = F.k;
  gp.cls  = asInteger(sexp_k) > 0;
  gp.p    = R.p;
  gp.seed = (uint64_t) asInteger(sexp_seed);
  gp.d    = (double **) calloc((size_t) (F.n_trees > 0 ? F.n_trees : 1), sizeof(double *));
  gp.u    = (int **)    calloc((size_t) (F.n_trees > 0 ? F.n_trees : 1), sizeof(int *));
  if (!gp.d || !gp.u) {
    free(gp.d); free(gp.u); free_rows(&R);
    error("Unable to allocate memory for permutation importance.");
  }

  int n_threads = asInteger(sexp_n_threads);
  if (n_threads > F.n_trees) n_threads = F.n_trees;
  if (n_threads < 1) n_threads = 1;
  if (F.n_trees > 0) run_parallel(perm_worker, n_threads, n_threads > 1 ? 100 : 1);

  int p = R.p, T = F.n_trees;
  const char *names[] = { "mean", "sd", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  SEXP mn = allocVector(REALSXP, p); SET_VECTOR_ELT(res, 0, mn);
  SEXP sd = allocVector(REALSXP, p); SET_VECTOR_ELT(res, 1, sd);
  double *m = REAL(mn), *s = REAL(sd);
  int    *n_nz = (int *) R_alloc((size_t) (p > 0 ? p : 1), sizeof(int));
  memset(m, 0, (size_t) p * sizeof(double));
  memset(s, 0, (size_t) p * sizeof(double));
  memset(n_nz, 0, (size_t) p * sizeof(int));

  int failed = gp.failed;
  if (!failed) {
    // Means and SDs over all trees, in tree order (trees not splitting on
    // a feature add 0 to its sum, and mean^2 to its squares).
    for (int t = 0; t < T; t++)
      for (int a = 0; gp.u[t][a] >= 0; a++) m[gp.u[t][a]] += gp.d[t][a];
    for (int j = 0; j < p; j++) m[j] /= (T > 0 ? T : 1);
    for (int t = 0; t < T; t++)
      for (int a = 0; gp.u[t][a] >= 0; a++) {
        int j = gp.u[t][a];
        double e = gp.d[t][a] - m[j];
        s[j] += e * e;
        n_nz[j]++;
      }
    for (int j = 0; j < p; j++) {
      s[j] += (double) (T - n_nz[j]) * m[j] * m[j];
      s[j] = T > 1 ? sqrt(s[j] / (T - 1)) : NA_REAL;
    }
  }

  for (int t = 0; t < T; t++) { free(gp.d[t]); free(gp.u[t]); }
  free(gp.d); free(gp.u);
  free_rows(&R);
  memset(&gp, 0, sizeof(gp));

  if (failed) error("Unable to allocate memory for permutation importance.");
  UNPROTECT(1);
  return res;
}


//======================================================
// C_rf_shuffle(x, seed, stream): Fisher-Yates on a PCG
// stream (the CV folds' shuffles).
//======================================================
SEXP C_rf_shuffle (SEXP sexp_x, SEXP sexp_seed, SEXP sexp_stream) {
  SEXP res = PROTECT(duplicate(sexp_x));
  int *v = INTEGER(res);
  int  n = length(res);
  rng_t rng;
  rng_seed(&rng, (uint64_t) asInteger(sexp_seed), (uint64_t) asReal(sexp_stream));
  for (int i = n - 1; i > 0; i--) {
    int j = (int) pcg32_bounded(&rng.pcg, (uint32_t) (i + 1));
    int t = v[i]; v[i] = v[j]; v[j] = t;
  }
  UNPROTECT(1);
  return res;
}


//======================================================
// C_rf_rng_test(kind, seed, n, count): for the tests.
//   kind 0  RandomState(seed).random_sample(count)
//   kind 1  RandomState(seed).randint(0, n, count)
//   kind 2  our_rand_r() count times from state `seed`
//   kind 3  the forest's tree seeds: RandomState(seed)
//           .randint(2^31 - 1), count times
//   kind 4  np.sum() of the doubles in `n` (pairwise)
//   kind 5  RandomState(seed).permutation(n), 0-based
//           (count must be n): the legacy shuffle,
//           j = random_interval(i) for i = n - 1, ..., 1,
//           which train_test_split() draws from
//======================================================
SEXP C_rf_rng_test (SEXP sexp_kind, SEXP sexp_seed, SEXP sexp_n, SEXP sexp_count) {
  int kind  = asInteger(sexp_kind);
  int count = asInteger(sexp_count);
  if (kind == 4) return ScalarReal(np_sum(REAL(sexp_n), length(sexp_n)));
  uint32_t seed = (uint32_t) asReal(sexp_seed);
  double   n    = asReal(sexp_n);
  SEXP res = PROTECT(allocVector(REALSXP, count));
  double *out = REAL(res);
  mt_t m;
  mt_seed(&m, seed);
  uint32_t xs = seed;
  if (kind == 5) {
    for (int i = 0; i < count; i++) out[i] = i;
    for (int i = count - 1; i > 0; i--) {
      int j = (int) mt_randint(&m, (uint64_t) i + 1);   // random_interval(i)
      double t = out[i]; out[i] = out[j]; out[j] = t;
    }
    UNPROTECT(1);
    return res;
  }
  for (int i = 0; i < count; i++) {
    switch (kind) {
      case 0:  out[i] = mt_double(&m); break;
      case 1:  out[i] = (double) mt_randint(&m, (uint64_t) n); break;
      case 2:  out[i] = (double) our_rand_r(&xs); break;
      default: out[i] = (double) mt_randint(&m, 2147483647ULL); break;
    }
  }
  UNPROTECT(1);
  return res;
}
