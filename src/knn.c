// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Exact k nearest neighbours of every point, straight from the
 * packed lower triangle of a dist object, for umap() and
 * tsne().
 *
 * Neighbours are ordered by distance, ties going to the lower
 * index, as R's order() does in uwot's dist_nn() and as bhtsne
 * does. Each task finds the neighbours of KNN_ROWS points.
 * Distances to earlier points come from other columns of the
 * packed triangle; they are read as contiguous runs, one per
 * column, for all of the task's points at once, which is far
 * faster than gathering each point's row with a stride.
 *
 * Determinism: each point's neighbours are found by one
 * thread, offering candidates in a fixed order.
 */

#include "ecodive.h"

// Points per parallel task.
#define KNN_ROWS 64


static int           kn_n;     // Points.
static int           kn_k;     // Neighbours per point.
static int           kn_self;  // Count each point among its own?
static const double *kn_d;     // Packed lower triangle.
static int          *kn_idx;   // n x k, row-major, 0-based.
static double       *kn_dist;  // n x k, row-major.


// Is (d1, j1) closer than (d2, j2)?
static inline int nbr_less(double d1, int j1, double d2, int j2) {
  return d1 < d2 || (d1 == d2 && j1 < j2);
}


// Offer (d, j) to a max-heap of the k closest so far.
static inline void heap_offer(double *hd, int *hj, int *size, int k, double d, int j) {

  int i;

  if (*size < k) {
    i = (*size)++;
    while (i > 0) {
      int p = (i - 1) / 2;
      if (!nbr_less(hd[p], hj[p], d, j)) break;
      hd[i] = hd[p]; hj[i] = hj[p];
      i = p;
    }
    hd[i] = d; hj[i] = j;
    return;
  }

  if (!nbr_less(d, j, hd[0], hj[0])) return;

  // Replace the root and sift down.
  i = 0;
  for (;;) {
    int c = 2 * i + 1;
    if (c >= k) break;
    if (c + 1 < k && nbr_less(hd[c], hj[c], hd[c + 1], hj[c + 1])) c++;
    if (!nbr_less(d, j, hd[c], hj[c])) break;
    hd[i] = hd[c]; hj[i] = hj[c];
    i = c;
  }
  hd[i] = d; hj[i] = j;
}


// Sort a full max-heap into ascending order, in place.
static void heap_sort(double *hd, int *hj, int k) {

  for (int end = k - 1; end > 0; end--) {

    double d = hd[end]; int j = hj[end];
    hd[end] = hd[0];    hj[end] = hj[0];

    int i = 0;
    for (;;) {
      int c = 2 * i + 1;
      if (c >= end) break;
      if (c + 1 < end && nbr_less(hd[c], hj[c], hd[c + 1], hj[c + 1])) c++;
      if (!nbr_less(d, j, hd[c], hj[c])) break;
      hd[i] = hd[c]; hj[i] = hj[c];
      i = c;
    }
    hd[i] = d; hj[i] = j;
  }
}


static void *knn_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = kn_n;
  int k         = kn_k;
  int n_tasks   = (n + KNN_ROWS - 1) / KNN_ROWS;

  int sizes[KNN_ROWS];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int i0 = task * KNN_ROWS;
    int i1 = (i0 + KNN_ROWS < n) ? i0 + KNN_ROWS : n;
    int nr = i1 - i0;

    // Heaps live in the output rows.
    for (int r = 0; r < nr; r++) sizes[r] = 0;

    #define HEAP(i) kn_dist + (size_t) (i) * k, kn_idx + (size_t) (i) * k, sizes + ((i) - i0), k

    // Columns left of the block.
    for (int j = 0; j < i0; j++) {
      const double *d = kn_d + dist_col_start(j, n) + (i0 - j - 1);
      for (int r = 0; r < nr; r++)
        heap_offer(HEAP(i0 + r), d[r], j);
    }

    for (int i = i0; i < i1; i++) {

      // Earlier points within the block.
      for (int j = i0; j < i; j++)
        heap_offer(HEAP(i), kn_d[dist_col_start(j, n) + (i - j - 1)], j);

      if (kn_self) heap_offer(HEAP(i), 0.0, i);

      // Later points: column i, contiguous.
      const double *d = kn_d + dist_col_start(i, n);
      for (int j = i + 1; j < n; j++)
        heap_offer(HEAP(i), d[j - i - 1], j);

      heap_sort(kn_dist + (size_t) i * k, kn_idx + (size_t) i * k, k);
    }

    #undef HEAP
  }

  return NULL;
}


//------------------------------------------------------
// The k nearest neighbours of each of n points: idx and
// dist, n x k row-major, nearest first. With
// `include_self`, each point counts among its own (k <= n);
// otherwise not (k <= n - 1).
//------------------------------------------------------
void knn_from_dist (
    const double *dist, int n, int k, int include_self,
    int *idx, double *dd, int n_threads ) {

  kn_n    = n;
  kn_k    = k;
  kn_self = include_self;
  kn_d    = dist;
  kn_idx  = idx;
  kn_dist = dd;

  run_parallel(knn_worker, n_threads, n);
}
