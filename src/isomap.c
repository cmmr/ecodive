// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The neighbor graph, its joining of disconnected components, and
// the geodesic distances follow scikit-learn's Isomap
// (sklearn/manifold/_isomap.py, sklearn/neighbors/_base.py and
// sklearn/utils/graph.py, version 1.9.1), and are distributed under
// its license as well:
//
//   BSD 3-Clause License
//
//   Copyright (c) 2007-2026 The scikit-learn developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in
//     the documentation and/or other materials provided with the
//     distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//   POSSIBILITY OF SUCH DAMAGE.

/*
 * Isomap: geodesic distances through a neighbor graph.
 *
 * This file does everything up to the geodesic distances; the
 * embedding is then classical scaling (PCoA) of them, by
 * C_pcoa() in pcoa.c, called from R. The distances have already
 * been checked by R's dist_vector(): finite, >= 0, <= 1e150.
 *
 * 1. Neighbor graph. Either each sample's `K` nearest other
 *    samples (scikit-learn's kneighbors_graph(), made undirected
 *    as its shortest_path(directed = FALSE) does), found by
 *    knn_from_dist() in knn.c; or every other sample within
 *    distance `r` (radius_neighbors_graph()). Ties in distance
 *    are broken by sample index, so the graph is defined by a
 *    total order. Edges of weight 0 (duplicate samples) are
 *    edges, as in SciPy's sparse graphs.
 *
 *    Radius neighbors are found straight from the packed lower
 *    triangle, in fixed blocks of ISO_ROWS rows. Each block reads
 *    the runs of the columns to its left that hold its rows, then
 *    its own columns.
 *
 * 2. Components. If the graph has more than one connected
 *    component, every pair of components is joined by one edge
 *    between their closest pair of samples, as scikit-learn's
 *    _fix_connected_components() does. Components are numbered in
 *    order of their first sample (as SciPy's connected_components()
 *    labels them), and for components a > b, ties are broken by
 *    the sample in a, then the sample in b (the first minimum in
 *    scikit-learn's row-major argmin). Unlike scikit-learn, joining
 *    is refused (an R error) when its c(c - 1) / 2 edges would be
 *    more than `ratio` times the graph's own edges and more than
 *    `floor` (10 and 65,536 from R), which bounds its extra memory
 *    and search time relative to a connected graph.
 *
 * 3. Geodesic distances. One Dijkstra search (4-ary heap) per
 *    source sample j, which stops once every sample after j is
 *    settled, and writes column j of the packed lower triangle:
 *    d(i, j) for i > j. So each pair is computed once, as the
 *    path length from the lower-numbered sample, summed edge by
 *    edge along the path from it.
 *
 * The radius search, the joining scan and the searches each run
 * in chunks, with interrupts checked between them.
 *
 * Determinism: results are bit-identical for any number of
 * threads. Each row's neighbors, each pair of components' joining
 * edge, and each source's search are each computed by one thread
 * with a fixed order of operations, and all comparisons use total
 * orders (distance, then index). The minima for joining are
 * reduced over per-thread tables in a fixed order (and are
 * order-independent in any case). Which thread does which piece of
 * work never affects any value.
 */

#include "ecodive.h"

// Rows per radius-search task. Fixed, so results don't depend on
// the number of threads (they don't anyway; see above).
#define ISO_ROWS 64

// Work is split into chunks, one run_parallel() call each, with
// interrupts checked between them: about ISO_CHUNK_READS distances
// read per chunk of the radius search or the joining scan, and
// about ISO_CHUNK_EDGES edges relaxed (at most ISO_CHUNK sources) per
// chunk of the shortest-path searches. Chunks never change results.
#define ISO_CHUNK_READS ((size_t) 1 << 27)
#define ISO_CHUNK_EDGES ((size_t) 1 << 28)
#define ISO_CHUNK       256


typedef struct {
  double w; // Edge weight (distance).
  int    j; // Neighbor.
} edge_t;

typedef struct {
  double w; // Distance between u and v.
  int    u; // Sample in the higher-numbered component.
  int    v; // Sample in the lower-numbered component.
} join_t;

typedef struct {
  double d; // Tentative distance ...
  int    v; // ... of this sample.
} hent_t;

typedef struct {
  double *dist;
  hent_t *heap;
  int    *pos;  // Position in heap; -1 = not reached; -2 = settled.
} work_t;


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the output slots that each task owns.
//------------------------------------------------------
static int           iso_n;      // Samples.
static const double *iso_d;      // Packed lower triangle of distances.
static double        iso_r;      // Radius.
static size_t       *adj_off;    // Undirected graph, CSR: n + 1 offsets ...
static edge_t       *adj;        // ... and edges, sorted by neighbor.
static int          *adj_deg;    // Row lengths after removing duplicates.
static size_t       *jadj_off;   // Joining edges, CSR: n + 1 offsets ...
static edge_t       *jadj;       // ... and edges.
static const int    *iso_mem;    // Each sample's component.
static join_t       *jt;         // Per-thread tables of closest pairs.
static size_t        jt_size;    // Entries per table.
static double       *iso_geo;    // Packed lower triangle of geodesics.
static double       *iso_max;    // Per source: largest geodesic.
static int          *iso_reach;  // Per source: 1 if all targets reached.
static work_t       *iso_work;   // Per thread.
static int           iso_lo;     // This chunk's first task ...
static int           iso_hi;     // ... and the one after its last.
static int           iso_next;   // Next source to search ...
static int           iso_end;    // ... and the end of this chunk.

#ifdef HAVE_PTHREAD
  static pthread_mutex_t iso_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif


//------------------------------------------------------
// Radius neighbors: every other sample within iso_r. In
// the first pass (fill == 0), counts each row's neighbors
// into adj_off[i + 1]; in the second, writes them into
// adj at the rows' offsets. Rows come out sorted by
// neighbor index.
//------------------------------------------------------
static void radius_pass(int thread_i, int n_threads, int fill) {

  int    n     = iso_n;
  double r_max = iso_r;

  size_t cnt[ISO_ROWS];

  for (int task = iso_lo + thread_i; task < iso_hi; task += n_threads) {

    int i0 = task * ISO_ROWS;
    int i1 = (i0 + ISO_ROWS < n) ? i0 + ISO_ROWS : n;
    int nr = i1 - i0;

    for (int r = 0; r < nr; r++) cnt[r] = 0;

    #define ADD(ROW, NBR, X) {                                     \
      if (fill) adj[adj_off[ROW] + cnt[(ROW) - i0]] = (edge_t){ X, NBR }; \
      cnt[(ROW) - i0]++;                                           \
    }

    for (int j = 0; j < i0; j++) {
      const double *d = iso_d + dist_col_start(j, n) + (i0 - j - 1);
      for (int r = 0; r < nr; r++) {
        double x = d[r];
        if (x <= r_max) ADD(i0 + r, j, x)
      }
    }

    for (int j = i0; j < i1; j++) {
      const double *d = iso_d + dist_col_start(j, n) - (j + 1);
      for (int i = j + 1; i < n; i++) {
        double x = d[i];
        if (x <= r_max) {
          ADD(j, i, x)
          if (i < i1) ADD(i, j, x)
        }
      }
    }

    #undef ADD

    if (!fill)
      for (int r = 0; r < nr; r++) adj_off[i0 + r + 1] = cnt[r];
  }
}

static void *radius_count(void *arg) {
  radius_pass(((worker_t *) arg)->i, ((worker_t *) arg)->n, 0);
  return NULL;
}

static void *radius_fill(void *arg) {
  radius_pass(((worker_t *) arg)->i, ((worker_t *) arg)->n, 1);
  return NULL;
}


//------------------------------------------------------
// Sort each row of adj by neighbor, and drop duplicates
// (an edge listed from both of its ends), recording the
// new lengths in adj_deg.
//------------------------------------------------------
static int edge_cmp(const void *a, const void *b) {
  int ja = ((const edge_t *) a)->j;
  int jb = ((const edge_t *) b)->j;
  return (ja > jb) - (ja < jb);
}

static void *sort_rows(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  for (int i = thread_i; i < iso_n; i += n_threads) {
    edge_t *e = adj + adj_off[i];
    size_t  m = adj_off[i + 1] - adj_off[i];
    qsort(e, m, sizeof(edge_t), edge_cmp);
    size_t u = 0;
    for (size_t t = 0; t < m; t++)
      if (u == 0 || e[t].j != e[u - 1].j) e[u++] = e[t];
    adj_deg[i] = (int) u;
  }

  return NULL;
}


//------------------------------------------------------
// The closest pair of samples between each pair of
// components, ordered by distance, then the sample in the
// higher-numbered component, then the other. Each thread
// scans a fixed set of columns into its own table.
//------------------------------------------------------
static inline int join_lt(double w, int u, int v, const join_t *b) {
  return w < b->w || (w == b->w && (u < b->u || (u == b->u && v < b->v)));
}

static void *join_worker(void *arg) {

  int     thread_i  = ((worker_t *) arg)->i;
  int     n_threads = ((worker_t *) arg)->n;
  int     n         = iso_n;
  join_t *tab       = jt + (size_t) thread_i * jt_size;

  for (int j = iso_lo + thread_i; j < iso_hi; j += n_threads) {
    const double *d  = iso_d + dist_col_start(j, n) - (j + 1);
    int           cj = iso_mem[j];
    for (int i = j + 1; i < n; i++) {
      int ci = iso_mem[i];
      if (ci == cj) continue;
      int a, b, u, v;
      if (ci > cj) { a = ci; b = cj; u = i; v = j; }
      else         { a = cj; b = ci; u = j; v = i; }
      join_t *t = tab + (size_t) a * (a - 1) / 2 + b;
      if (join_lt(d[i], u, v, t)) *t = (join_t){ d[i], u, v };
    }
  }

  return NULL;
}


//------------------------------------------------------
// Indexed 4-ary min-heap of samples, ordered by their
// tentative distance, then index. Entries carry their
// distance, so that comparisons stay within the heap.
// The order is total, so the heap's shape never affects
// which sample comes out next.
//------------------------------------------------------
static inline int hent_lt(hent_t a, hent_t b) {
  return a.d < b.d || (a.d == b.d && a.v < b.v);
}

static void node_sift_up(work_t *w, int p) {
  hent_t *heap = w->heap;
  int    *pos  = w->pos;
  hent_t  x    = heap[p];
  while (p > 0) {
    int q = (p - 1) / 4;
    if (!hent_lt(x, heap[q])) break;
    heap[p]        = heap[q];
    pos[heap[p].v] = p;
    p = q;
  }
  heap[p] = x;
  pos[x.v] = p;
}

static void node_sift_down(work_t *w, int p, int m) {
  hent_t *heap = w->heap;
  int    *pos  = w->pos;
  hent_t  x    = heap[p];
  while (1) {
    int c = 4 * p + 1;
    if (c >= m) break;
    int e = (c + 4 < m) ? c + 4 : m;
    for (int t = c + 1; t < e; t++)
      if (hent_lt(heap[t], heap[c])) c = t;
    if (!hent_lt(heap[c], x)) break;
    heap[p]        = heap[c];
    pos[heap[p].v] = p;
    p = c;
  }
  heap[p] = x;
  pos[x.v] = p;
}


//------------------------------------------------------
// Dijkstra's search from source s, until every sample
// after s is settled. Writes column s of the geodesics.
//
// relax() needs no test for settled samples: a settled
// sample v has dist[v] <= du <= du + w (weights are >= 0,
// and rounding is monotonic), so it is never updated.
//------------------------------------------------------
static void relax(work_t *w, int *hn, double du, const edge_t *e, size_t m) {
  double *dist = w->dist;
  int    *pos  = w->pos;
  for (size_t t = 0; t < m; t++) {
    int    v  = e[t].j;
    double nd = du + e[t].w;
    if (nd < dist[v]) {
      dist[v] = nd;
      if (pos[v] < 0) {
        pos[v] = *hn;
        (*hn)++;
      }
      w->heap[pos[v]] = (hent_t){ nd, v };
      node_sift_up(w, pos[v]);
    }
  }
}

static void dijkstra(int s, work_t *w) {

  int     n    = iso_n;
  double *dist = w->dist;
  hent_t *heap = w->heap;
  int    *pos  = w->pos;

  for (int v = 0; v < n; v++) {
    dist[v] = INFINITY;
    pos[v]  = -1;
  }

  dist[s] = 0;
  heap[0] = (hent_t){ 0, s };
  pos[s]  = 0;

  int hn        = 1;
  int remaining = n - 1 - s;

  while (hn > 0) {

    int u = heap[0].v;
    hn--;
    if (hn > 0) {
      heap[0] = heap[hn];
      node_sift_down(w, 0, hn);
    }
    pos[u] = -2;

    if (u > s && --remaining == 0) break;

    relax(w, &hn, dist[u], adj + adj_off[u], (size_t) adj_deg[u]);
    if (jadj) relax(w, &hn, dist[u], jadj + jadj_off[u], jadj_off[u + 1] - jadj_off[u]);
  }

  double *g  = iso_geo + dist_col_start(s, n) - (s + 1);
  double  mx = 0;
  for (int i = s + 1; i < n; i++) {
    g[i] = dist[i];
    if (dist[i] > mx) mx = dist[i];
  }

  iso_max[s]   = mx;
  iso_reach[s] = (remaining == 0);
}

// Each share takes the next unsearched source until the
// chunk is done. Which share searches which source never
// affects the result: each search writes its own column.
static void *dijkstra_worker(void *arg) {

  work_t *w = iso_work + ((worker_t *) arg)->i;

  while (1) {

    int s;

    #ifdef HAVE_PTHREAD
      pthread_mutex_lock(&iso_mutex);
    #endif
    s = (iso_next < iso_end) ? iso_next++ : -1;
    #ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&iso_mutex);
    #endif

    if (s < 0) break;

    dijkstra(s, w);
  }

  return NULL;
}


//------------------------------------------------------
// Run `func` over tasks 0 to total - 1, `chunk` at a time
// (one run_parallel() call each; workers take tasks
// iso_lo + i, iso_lo + i + n_threads, ... below iso_hi),
// checking for interrupts between chunks. Everything is
// allocated with R_alloc(), so an interrupt leaks nothing.
//------------------------------------------------------
static void run_chunked(pthread_func_t func, int n_threads, int total, size_t chunk) {
  for (int lo = 0; lo < total; lo += (int) chunk) {
    iso_lo = lo;
    iso_hi = ((size_t) (total - lo) > chunk) ? lo + (int) chunk : total;
    run_parallel(func, n_threads, iso_n);
    if (iso_hi < total) R_CheckUserInterrupt();
  }
}

// About ISO_CHUNK_READS / `per` tasks (`per` reads each),
// but at least `lo` >= 1.
static size_t chunk_of(size_t per, size_t lo) {
  size_t c = ISO_CHUNK_READS / per;
  return (c > lo) ? c : lo;
}


//======================================================
// R interface. `sexp_dist` is a double vector holding a
// dist object's packed lower triangle for `sexp_n`
// samples (n >= 2), checked by dist_vector(). Neighbors
// are the `sexp_k` nearest (1 <= k <= n - 1), or, if
// sexp_k is 0, all within `sexp_radius`. `sexp_join` is
// c(bytes, ratio, floor): the per-thread tables of
// closest pairs between components may take `bytes` in
// all, with fewer threads scanning past that (but at
// least one, whose table is always made); and joining c
// components, which adds c(c - 1) / 2 edges, is refused
// with an R error when that is more than `ratio` times
// the neighbor graph's own edges and more than `floor`.
// If `sexp_graph` is TRUE, also returns the final graph's
// edges. Returns list(geo, membership, components, max,
// edges, weights): the geodesics as a dist vector; each
// sample's component (1-based) and the number of
// components before joining; the largest geodesic; and,
// if asked for, each undirected edge (i < j, 1-based)
// with its weight, joining edges included.
//======================================================
SEXP C_isomap(
    SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_radius,
    SEXP sexp_n_threads, SEXP sexp_join, SEXP sexp_graph ) {

  int    n          = asInteger(sexp_n);
  int    K          = asInteger(sexp_k);
  double radius     = asReal(sexp_radius);
  int    n_threads  = asInteger(sexp_n_threads);
  int    join_ok    = TYPEOF(sexp_join) == REALSXP && XLENGTH(sexp_join) == 3;
  double join_bytes = join_ok ? REAL(sexp_join)[0] : -1;
  double join_ratio = join_ok ? REAL(sexp_join)[1] : -1;
  double join_floor = join_ok ? REAL(sexp_join)[2] : -1;
  int    want_graph = asLogical(sexp_graph) == 1;

  if (n < 2 || K < 0 || K > n - 1 || n_threads < 1 ||
      (K == 0 && !(radius >= 0)) || TYPEOF(sexp_dist) != REALSXP ||
      !(join_bytes >= 0) || !(join_ratio >= 0) || !(join_floor >= 0) ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_isomap()."); // # nocov

  // More threads than tasks would sit idle.
  int n_tasks     = (n + ISO_ROWS - 1) / ISO_ROWS;
  int row_threads = (n_threads < n_tasks) ? n_threads : n_tasks;
  int col_threads = (n_threads < n - 1)   ? n_threads : n - 1;

  iso_n    = n;
  iso_d    = REAL(sexp_dist);
  iso_r    = radius;
  jadj     = NULL;
  jadj_off = NULL;

  adj_off = (size_t *) R_alloc((size_t) n + 1, sizeof(size_t));
  adj_deg = (int *)    R_alloc(n, sizeof(int));


  //--------------------------------------------------
  // Neighbor graph, undirected, as CSR rows sorted by
  // neighbor.
  //--------------------------------------------------
  if (K > 0) {

    // Each row's K nearest, by distance, then index.
    int    *nb_idx = (int *)    R_alloc((size_t) n * K, sizeof(int));
    double *nb_w   = (double *) R_alloc((size_t) n * K, sizeof(double));

    knn_from_dist(iso_d, n, K, 0, nb_idx, nb_w, n_threads);
    R_CheckUserInterrupt();

    // Each directed edge i -> j is listed in rows i and j.
    memset(adj_off, 0, ((size_t) n + 1) * sizeof(size_t));
    for (size_t e = 0; e < (size_t) n * K; e++)
      adj_off[nb_idx[e] + 1]++;
    for (int i = 0; i < n; i++)
      adj_off[i + 1] += adj_off[i] + (size_t) K;

    adj = (edge_t *) R_alloc(adj_off[n], sizeof(edge_t));

    size_t *fill = (size_t *) R_alloc(n, sizeof(size_t));
    memcpy(fill, adj_off, (size_t) n * sizeof(size_t));

    for (int i = 0; i < n; i++)
      for (int t = 0; t < K; t++) {
        size_t e = (size_t) i * K + t;
        int    j = nb_idx[e];
        adj[fill[i]++] = (edge_t){ nb_w[e], j };
        adj[fill[j]++] = (edge_t){ nb_w[e], i };
      }
  }
  else {

    memset(adj_off, 0, ((size_t) n + 1) * sizeof(size_t));
    size_t chunk = chunk_of((size_t) 2 * ISO_ROWS * n, 2 * (size_t) row_threads);
    run_chunked(radius_count, row_threads, n_tasks, chunk);

    for (int i = 0; i < n; i++) adj_off[i + 1] += adj_off[i];

    adj = (edge_t *) R_alloc(adj_off[n] ? adj_off[n] : 1, sizeof(edge_t));
    run_chunked(radius_fill, row_threads, n_tasks, chunk);
  }

  run_parallel(sort_rows, col_threads, n);


  //--------------------------------------------------
  // Connected components, numbered in order of their
  // first sample.
  //--------------------------------------------------
  SEXP sexp_mem = PROTECT(allocVector(INTSXP, n));
  int *mem      = INTEGER(sexp_mem);
  int *queue    = (int *) R_alloc(n, sizeof(int));
  int  ncomp    = 0;

  for (int i = 0; i < n; i++) mem[i] = -1;

  for (int i = 0; i < n; i++) {
    if (mem[i] >= 0) continue;
    int qh = 0, qt = 0;
    queue[qt++] = i;
    mem[i]      = ncomp;
    while (qh < qt) {
      int     u = queue[qh++];
      edge_t *e = adj + adj_off[u];
      for (int t = 0; t < adj_deg[u]; t++)
        if (mem[e[t].j] < 0) {
          mem[e[t].j] = ncomp;
          queue[qt++] = e[t].j;
        }
    }
    ncomp++;
  }

  iso_mem = mem;


  //--------------------------------------------------
  // Join every pair of components at their closest
  // samples.
  //--------------------------------------------------
  if (ncomp > 1) {

    jt_size = (size_t) ncomp * (ncomp - 1) / 2;

    // Refuse joins that would swamp the graph's own edges.
    double n_edges = 0;
    for (int i = 0; i < n; i++) n_edges += adj_deg[i];
    n_edges /= 2;
    if ((double) jt_size > join_ratio * n_edges && (double) jt_size > join_floor)
      error(
        "The neighbor graph has %d connected components. Joining them would "
        "add %.0f edges, more than %g times the graph's own %.0f edges. "
        "Increase `%s`.",
        ncomp, (double) jt_size, join_ratio, n_edges,
        (K > 0) ? "n_neighbors" : "radius" );

    int jt_threads = col_threads;
    size_t per     = jt_size * sizeof(join_t);
    if ((double) jt_threads * (double) per > join_bytes) {
      double fit = floor(join_bytes / (double) per);
      jt_threads = (fit < 1) ? 1 : (int) fit;
    }

    jt = (join_t *) R_alloc((size_t) jt_threads * jt_size, sizeof(join_t));
    for (size_t t = 0; t < (size_t) jt_threads * jt_size; t++)
      jt[t] = (join_t){ INFINITY, n, n };

    run_chunked(join_worker, jt_threads, n - 1, chunk_of((size_t) n, 2 * (size_t) jt_threads));

    for (int t = 1; t < jt_threads; t++)
      for (size_t p = 0; p < jt_size; p++) {
        join_t *b = jt + (size_t) t * jt_size + p;
        if (join_lt(b->w, b->u, b->v, jt + p)) jt[p] = *b;
      }

    // Joining edges, as CSR.
    jadj_off = (size_t *) R_alloc((size_t) n + 1, sizeof(size_t));
    memset(jadj_off, 0, ((size_t) n + 1) * sizeof(size_t));
    for (size_t p = 0; p < jt_size; p++) {
      jadj_off[jt[p].u + 1]++;
      jadj_off[jt[p].v + 1]++;
    }
    for (int i = 0; i < n; i++) jadj_off[i + 1] += jadj_off[i];

    jadj = (edge_t *) R_alloc(jadj_off[n], sizeof(edge_t));
    size_t *fill = (size_t *) R_alloc(n, sizeof(size_t));
    memcpy(fill, jadj_off, (size_t) n * sizeof(size_t));
    for (size_t p = 0; p < jt_size; p++) {
      jadj[fill[jt[p].u]++] = (edge_t){ jt[p].w, jt[p].v };
      jadj[fill[jt[p].v]++] = (edge_t){ jt[p].w, jt[p].u };
    }
  }


  //--------------------------------------------------
  // Geodesic distances: one search per source, `cpus`
  // at a time, in chunks between interrupt checks.
  //--------------------------------------------------
  SEXP sexp_geo = PROTECT(allocVector(REALSXP, XLENGTH(sexp_dist)));
  iso_geo   = REAL(sexp_geo);
  iso_max   = (double *) R_alloc(n, sizeof(double));
  iso_reach = (int *)    R_alloc(n, sizeof(int));

  int dj_threads = col_threads;
  iso_work = (work_t *) R_alloc(dj_threads, sizeof(work_t));
  for (int t = 0; t < dj_threads; t++) {
    iso_work[t].dist = (double *) R_alloc(n, sizeof(double));
    iso_work[t].heap = (hent_t *) R_alloc(n, sizeof(hent_t));
    iso_work[t].pos  = (int *)    R_alloc(n, sizeof(int));
  }

  // Sources per chunk: about ISO_CHUNK_EDGES edges relaxed
  // (each search relaxes up to n + all edges), at most
  // ISO_CHUNK, and at least one per thread.
  size_t per_source = (size_t) n + adj_off[n] + (jadj ? jadj_off[n] : 0);
  size_t dj_chunk   = ISO_CHUNK_EDGES / per_source;
  if (dj_chunk > ISO_CHUNK)            dj_chunk = ISO_CHUNK;
  if (dj_chunk < (size_t) dj_threads)  dj_chunk = (size_t) dj_threads;

  for (int s0 = 0; s0 < n - 1; s0 += (int) dj_chunk) {
    iso_next = s0;
    iso_end  = ((size_t) (n - 1 - s0) > dj_chunk) ? s0 + (int) dj_chunk : n - 1;
    run_parallel(dijkstra_worker, dj_threads, n);
    if (iso_end < n - 1) R_CheckUserInterrupt();
  }

  double mx = 0;
  for (int s = 0; s < n - 1; s++) {
    // The graph is connected, so a target is only out of
    // reach if its path length overflowed, which distances
    // of at most 1e150 (from dist_vector()) can't do.
    if (!iso_reach[s]) // # nocov start
      error("`dm` is too large: its geodesic distances overflow. Rescale `dm`.");
    // # nocov end
    if (iso_max[s] > mx) mx = iso_max[s];
  }


  //--------------------------------------------------
  // The final graph's edges, if asked for (by tests).
  //--------------------------------------------------
  SEXP sexp_edges = R_NilValue, sexp_wts = R_NilValue;

  if (want_graph) {

    R_xlen_t m = 0;
    for (int i = 0; i < n; i++) {
      for (int t = 0; t < adj_deg[i]; t++)
        if (adj[adj_off[i] + t].j > i) m++;
      if (jadj)
        for (size_t t = jadj_off[i]; t < jadj_off[i + 1]; t++)
          if (jadj[t].j > i) m++;
    }

    sexp_edges = PROTECT(allocMatrix(INTSXP, (int) m, 2));
    sexp_wts   = PROTECT(allocVector(REALSXP, m));
    int    *ed = INTEGER(sexp_edges);
    double *wt = REAL(sexp_wts);

    R_xlen_t p = 0;
    for (int i = 0; i < n; i++) {
      for (int t = 0; t < adj_deg[i]; t++) {
        edge_t e = adj[adj_off[i] + t];
        if (e.j > i) { ed[p] = i + 1; ed[p + m] = e.j + 1; wt[p] = e.w; p++; }
      }
      if (jadj)
        for (size_t t = jadj_off[i]; t < jadj_off[i + 1]; t++) {
          edge_t e = jadj[t];
          if (e.j > i) { ed[p] = i + 1; ed[p + m] = e.j + 1; wt[p] = e.w; p++; }
        }
    }
  }
  else {
    PROTECT(sexp_edges);
    PROTECT(sexp_wts);
  }

  for (int i = 0; i < n; i++) mem[i]++;

  const char *names[] = { "geo", "membership", "components", "max", "edges", "weights", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_geo);
  SET_VECTOR_ELT(sexp_result, 1, sexp_mem);
  SET_VECTOR_ELT(sexp_result, 2, ScalarInteger(ncomp));
  SET_VECTOR_ELT(sexp_result, 3, ScalarReal(mx));
  SET_VECTOR_ELT(sexp_result, 4, sexp_edges);
  SET_VECTOR_ELT(sexp_result, 5, sexp_wts);

  UNPROTECT(5);
  return sexp_result;
}
