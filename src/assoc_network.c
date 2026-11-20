// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The graph computations of assoc_network(), written from the papers cited
// below and NetCoMi's documentation (Peschel et al. 2021). NetCoMi and
// igraph are GPL and were used for their outputs only; none of their
// source was read. See plans/assoc_network.md.
//
// Graphs are undirected and weighted, given as edge lists (0-based `from`,
// `to`) over n nodes, and turned into a compressed adjacency list here.
//
//   C_network_paths()   Shortest paths from every node (Dijkstra, on edge
//                       lengths), with Brandes' (2001) accumulation for
//                       betweenness, and each node's sum of distances and
//                       of inverse distances to the nodes it reaches.
//   C_network_greedy()  Clauset, Newman & Moore's (2004) greedy modularity
//                       merges, on edge weights.
//   C_network_barrat()  Barrat et al.'s (2004) weighted local clustering
//                       coefficient.
//   C_network_eigen()   The leading eigenvector of the weighted adjacency
//                       matrix (block Lanczos, lanczos.c).
//
// Threads: C_network_paths() splits the sources into fixed blocks of
// NW_BLOCK. A block's betweenness contributions are summed, in source
// order, into the block's own row; the rows are then summed in block
// order. Per-source values go to their own slots. So results are
// identical for any number of threads.
//
// Brandes, U. (2001). A faster algorithm for betweenness centrality. J Math
//   Sociol, 25, 163-177.
// Clauset, A., Newman, M. E. J., & Moore, C. (2004). Finding community
//   structure in very large networks. Phys Rev E, 70, 066111.
// Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A. (2004).
//   The architecture of complex weighted networks. PNAS, 101, 3747-3752.

#include "ecodive.h"

#include <limits.h> // INT_MAX
#include <math.h>

#define NW_BLOCK 32       // Sources per task of C_network_paths().
#define NW_RTOL  1e-12    // Relative tolerance for equal path lengths.


//------------------------------------------------------
// Compressed adjacency: node v's neighbours are
// adj[ptr[v]] .. adj[ptr[v + 1] - 1], each with its edge
// value val[] (a length or a weight). Neighbours are in
// increasing order. R_alloc()'d.
//------------------------------------------------------
typedef struct {
  int     n;
  int     m;      // Edges (each stored twice).
  int    *ptr;
  int    *adj;
  double *val;
} nw_graph_t;

static void nw_build(nw_graph_t *g, int n, SEXP sexp_from, SEXP sexp_to, SEXP sexp_val) {

  int m = LENGTH(sexp_from);
  if (LENGTH(sexp_to) != m || LENGTH(sexp_val) != m)
    error("Invalid arguments to the network functions."); // # nocov

  const int    *from = INTEGER(sexp_from);
  const int    *to   = INTEGER(sexp_to);
  const double *val  = REAL(sexp_val);

  g->n   = n;
  g->m   = m;
  g->ptr = (int *)    R_alloc((size_t) n + 1, sizeof(int));
  g->adj = (int *)    R_alloc((size_t) 2 * m + 1, sizeof(int));
  g->val = (double *) R_alloc((size_t) 2 * m + 1, sizeof(double));

  for (int v = 0; v <= n; v++) g->ptr[v] = 0;
  for (int e = 0; e < m; e++) {
    if (from[e] < 0 || from[e] >= n || to[e] < 0 || to[e] >= n || from[e] == to[e])
      error("Invalid arguments to the network functions."); // # nocov
    g->ptr[from[e] + 1]++;
    g->ptr[to[e]   + 1]++;
  }
  for (int v = 0; v < n; v++) g->ptr[v + 1] += g->ptr[v];

  // Fill, then sort each list by neighbour (insertion sort: lists are
  // short, and often nearly sorted).
  int *fill = (int *) R_alloc((size_t) n + 1, sizeof(int));
  for (int v = 0; v < n; v++) fill[v] = g->ptr[v];
  for (int e = 0; e < m; e++) {
    int a = from[e], b = to[e];
    g->adj[fill[a]] = b; g->val[fill[a]++] = val[e];
    g->adj[fill[b]] = a; g->val[fill[b]++] = val[e];
  }
  for (int v = 0; v < n; v++) {
    for (int x = g->ptr[v] + 1; x < g->ptr[v + 1]; x++) {
      int    a = g->adj[x];
      double w = g->val[x];
      int    y = x - 1;
      while (y >= g->ptr[v] && g->adj[y] > a) {
        g->adj[y + 1] = g->adj[y];
        g->val[y + 1] = g->val[y];
        y--;
      }
      g->adj[y + 1] = a;
      g->val[y + 1] = w;
    }
  }
}



//======================================================
// Shortest paths.
//======================================================

static nw_graph_t  np_g;
static int         np_nblk;     // Blocks of sources.
static double     *np_bc;       // np_nblk x n: each block's betweenness.
static double     *np_harm;     // n: sum of 1 / dist over the nodes reached.
static double     *np_dsum;     // n: sum of dist.
static int        *np_reach;    // n: nodes reached (not counting itself).
static char       *np_work;     // Per-thread workspace ...
static size_t      np_wlen;     // ... of this many bytes.


// Binary min-heap of (key, node) with lazy deletion.
typedef struct { double key; int node; } nw_heap_t;

static void heap_push(nw_heap_t *h, int *size, double key, int node) {
  int i = (*size)++;
  while (i > 0) {
    int p = (i - 1) / 2;
    if (h[p].key < key || (h[p].key == key && h[p].node <= node)) break;
    h[i] = h[p];
    i = p;
  }
  h[i].key  = key;
  h[i].node = node;
}

static nw_heap_t heap_pop(nw_heap_t *h, int *size) {
  nw_heap_t top  = h[0];
  nw_heap_t last = h[--(*size)];
  int i = 0, n = *size;
  for (;;) {
    int c = 2 * i + 1;
    if (c >= n) break;
    if (c + 1 < n && (h[c + 1].key < h[c].key ||
        (h[c + 1].key == h[c].key && h[c + 1].node < h[c].node))) c++;
    if (last.key < h[c].key || (last.key == h[c].key && last.node <= h[c].node)) break;
    h[i] = h[c];
    i = c;
  }
  if (n > 0) h[i] = last;
  return top;
}


static void *np_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  const nw_graph_t *g = &np_g;
  int n  = g->n;
  int m2 = 2 * g->m;

  char      *w     = np_work + (size_t) thread_i * np_wlen;
  double    *dist  = (double *) w;                 w += (size_t) n * sizeof(double);
  double    *sigma = (double *) w;                 w += (size_t) n * sizeof(double);
  double    *delta = (double *) w;                 w += (size_t) n * sizeof(double);
  nw_heap_t *heap  = (nw_heap_t *) w;              w += ((size_t) m2 + 1) * sizeof(nw_heap_t);
  int       *npred = (int *) w;                    w += (size_t) n * sizeof(int);
  int       *pred  = (int *) w;                    w += ((size_t) m2 + 1) * sizeof(int);
  int       *stack = (int *) w;                    w += (size_t) n * sizeof(int);
  char      *done  = w;

  for (int blk = thread_i; blk < np_nblk; blk += n_threads) {

    double *bc = np_bc + (size_t) blk * n;
    for (int v = 0; v < n; v++) bc[v] = 0;

    int s1 = (blk + 1) * NW_BLOCK;
    if (s1 > n) s1 = n;

    for (int s = blk * NW_BLOCK; s < s1; s++) {

      for (int v = 0; v < n; v++) {
        dist[v]  = INFINITY;
        sigma[v] = 0;
        delta[v] = 0;
        npred[v] = 0;
        done[v]  = 0;
      }
      dist[s]  = 0;
      sigma[s] = 1;

      int hs = 0, ns = 0;
      heap_push(heap, &hs, 0, s);

      while (hs > 0) {
        nw_heap_t top = heap_pop(heap, &hs);
        int v = top.node;
        if (done[v] || top.key != dist[v]) continue;
        done[v]     = 1;
        stack[ns++] = v;

        for (int x = g->ptr[v]; x < g->ptr[v + 1]; x++) {
          int u = g->adj[x];
          if (done[u]) continue;
          double alt = dist[v] + g->val[x];
          double cur = dist[u];
          if (isfinite(cur) && fabs(alt - cur) <= NW_RTOL * (alt > cur ? alt : cur)) {
            // An equally short path (to a relative 1e-12).
            sigma[u] += sigma[v];
            pred[g->ptr[u] + npred[u]++] = v;
          }
          else if (alt < cur) {
            dist[u]  = alt;
            sigma[u] = sigma[v];
            npred[u] = 1;
            pred[g->ptr[u]] = v;
            heap_push(heap, &hs, alt, u);
          }
        }
      }

      // Distances, then dependencies in order of decreasing distance.
      double harm = 0, dsum = 0;
      for (int t = 1; t < ns; t++) {
        int u = stack[t];
        harm += 1 / dist[u];
        dsum += dist[u];
      }
      np_harm[s]  = harm;
      np_dsum[s]  = dsum;
      np_reach[s] = ns - 1;

      for (int t = ns - 1; t > 0; t--) {
        int    u    = stack[t];
        double coef = (1 + delta[u]) / sigma[u];
        for (int x = 0; x < npred[u]; x++) {
          int v = pred[g->ptr[u] + x];
          delta[v] += sigma[v] * coef;
        }
        bc[u] += delta[u];
      }
    }
  }

  return NULL;
}


//======================================================
// R interface.
//
// n:     nodes.
// from, to: 0-based ends of each edge (integer).
// len:   each edge's length (double, positive).
// cpus:  threads.
//
// Returns list(betweenness, harmonic, dist_sum, reached).
// betweenness is Brandes' sum over ordered pairs, halved
// (each unordered pair once), not normalised.
//======================================================
SEXP C_network_paths(SEXP sexp_n, SEXP sexp_from, SEXP sexp_to, SEXP sexp_len, SEXP sexp_cpus) {

  int n    = asInteger(sexp_n);
  int cpus = asInteger(sexp_cpus);
  if (n < 1 || cpus < 1) error("Invalid arguments to C_network_paths()."); // # nocov

  nw_build(&np_g, n, sexp_from, sexp_to, sexp_len);
  for (int x = 0; x < 2 * np_g.m; x++)
    if (!(np_g.val[x] > 0) || !isfinite(np_g.val[x]))
      error("Invalid arguments to C_network_paths()."); // # nocov

  SEXP res = PROTECT(allocVector(VECSXP, 4));
  SEXP sb  = PROTECT(allocVector(REALSXP, n)); SET_VECTOR_ELT(res, 0, sb);
  SEXP sh  = PROTECT(allocVector(REALSXP, n)); SET_VECTOR_ELT(res, 1, sh);
  SEXP sd  = PROTECT(allocVector(REALSXP, n)); SET_VECTOR_ELT(res, 2, sd);
  SEXP sr  = PROTECT(allocVector(INTSXP,  n)); SET_VECTOR_ELT(res, 3, sr);

  np_harm  = REAL(sh);
  np_dsum  = REAL(sd);
  np_reach = INTEGER(sr);
  np_nblk  = (n + NW_BLOCK - 1) / NW_BLOCK;
  np_bc    = (double *) R_alloc((size_t) np_nblk * n, sizeof(double));

  int m2 = 2 * np_g.m;
  np_wlen = 3 * (size_t) n * sizeof(double) + ((size_t) m2 + 1) * sizeof(nw_heap_t) +
            (2 * (size_t) n + (size_t) m2 + 1) * sizeof(int) + (size_t) n;
  np_wlen = (np_wlen + 15) / 16 * 16;

  // Work: about (m + n log n) per source.
  int nt = threads_for(cpus, (double) n * ((double) m2 + n * 4.0), 2e5);
  if (nt > np_nblk) nt = np_nblk;
  np_work = (char *) R_alloc((size_t) nt * np_wlen, 1);

  run_parallel(np_worker, nt, nt > 1 ? 100 : 1);

  double *bc = REAL(sb);
  for (int v = 0; v < n; v++) bc[v] = 0;
  for (int b = 0; b < np_nblk; b++)
    for (int v = 0; v < n; v++)
      bc[v] += np_bc[(size_t) b * n + v];
  for (int v = 0; v < n; v++) bc[v] /= 2;

  UNPROTECT(5);
  return res;
}



//======================================================
// Greedy modularity (Clauset, Newman & Moore 2004).
//
// Every node starts in its own community. At each step,
// the two connected communities whose merge most
// increases the modularity Q are merged:
//
//   dQ_ij = 2 (e_ij - a_i a_j),
//
// with e_ij the fraction of the total edge weight (each
// edge counted in both directions) between i and j, and
// a_i that of the edge ends in i. Ties go to the pair
// (i, j), i < j, that comes first; a merged community
// keeps the smaller id. Merging goes on while any two
// communities are connected, and the partition kept is
// the one with the largest Q (the first, on ties), as in
// the dendrogram cut that Clauset et al. describe; but
// never one community of all n nodes (igraph's rule:
// then the partition one merge before it).
//
// Each community's links are a sorted list of (community,
// e) pairs. dQ is recomputed from e and a when compared,
// rather than updated, so every pair's value is computed
// the same way at every step.
//
// R interface: n, from, to (0-based), w (positive).
// Returns list(membership (0-based, by smallest member),
// modularity of each step, from 0 merges).
//======================================================
typedef struct { int c; double e; } nw_link_t;

SEXP C_network_greedy(SEXP sexp_n, SEXP sexp_from, SEXP sexp_to, SEXP sexp_w) {

  int n = asInteger(sexp_n);
  if (n < 1) error("Invalid arguments to C_network_greedy()."); // # nocov

  nw_graph_t g;
  nw_build(&g, n, sexp_from, sexp_to, sexp_w);

  double tot = 0;
  for (int x = 0; x < 2 * g.m; x++) {
    if (!(g.val[x] > 0) || !isfinite(g.val[x]))
      error("Invalid arguments to C_network_greedy()."); // # nocov
    tot += g.val[x];
  }

  SEXP res  = PROTECT(allocVector(VECSXP, 2));
  SEXP smem = PROTECT(allocVector(INTSXP, n));
  SET_VECTOR_ELT(res, 0, smem);
  int *mem = INTEGER(smem);
  for (int v = 0; v < n; v++) mem[v] = v;

  if (g.m == 0) {
    SET_VECTOR_ELT(res, 1, ScalarReal(0));
    UNPROTECT(2);
    return res;
  }

  // Links and a_i. Lists are malloc'd so they can grow.
  nw_link_t **lnk = (nw_link_t **) R_alloc((size_t) n, sizeof(nw_link_t *));
  int        *len = (int *)        R_alloc((size_t) n, sizeof(int));
  int        *cap = (int *)        R_alloc((size_t) n, sizeof(int));
  double     *a   = (double *)     R_alloc((size_t) n, sizeof(double));
  int        *mi  = (int *)        R_alloc((size_t) n, sizeof(int));  // merges: i
  int        *mj  = (int *)        R_alloc((size_t) n, sizeof(int));  // merges: j
  double     *q   = (double *)     R_alloc((size_t) n, sizeof(double));
  nw_link_t  *tmp = NULL;
  int         ok  = 1;

  double q0 = 0;
  for (int v = 0; v < n; v++) {
    len[v] = g.ptr[v + 1] - g.ptr[v];
    cap[v] = len[v] > 0 ? len[v] : 1;
    lnk[v] = (nw_link_t *) malloc((size_t) cap[v] * sizeof(nw_link_t));
    if (!lnk[v]) ok = 0;
    a[v] = 0;
    for (int x = 0; ok && x < len[v]; x++) {
      lnk[v][x].c = g.adj[g.ptr[v] + x];
      lnk[v][x].e = g.val[g.ptr[v] + x] / tot;
      a[v] += lnk[v][x].e;
    }
    q0 -= a[v] * a[v];
  }
  tmp = (nw_link_t *) malloc((size_t) 2 * n * sizeof(nw_link_t));
  if (!tmp) ok = 0;

  int steps = 0;
  q[0] = q0;

  while (ok) {

    // The best pair: largest dQ, ties to the smallest (i, j).
    int    bi = -1, bj = -1;
    double best = 0;
    for (int i = 0; i < n; i++) {
      for (int x = 0; x < len[i]; x++) {
        int j = lnk[i][x].c;
        if (j <= i) continue;
        double dq = 2 * (lnk[i][x].e - a[i] * a[j]);
        if (bi < 0 || dq > best) { best = dq; bi = i; bj = j; }
      }
    }
    if (bi < 0) break;

    // Merge bj into bi: bi's links become the union.
    int     li = len[bi], lj = len[bj], k = 0, x = 0, y = 0;
    nw_link_t *Li = lnk[bi], *Lj = lnk[bj];
    while (x < li || y < lj) {
      int ci = x < li ? Li[x].c : INT_MAX;
      int cj = y < lj ? Lj[y].c : INT_MAX;
      if (ci == cj) { tmp[k].c = ci; tmp[k].e = Li[x].e + Lj[y].e; x++; y++; }
      else if (ci < cj) { tmp[k++] = Li[x++]; continue; }
      else { tmp[k++] = Lj[y++]; continue; }
      k++;
    }
    // Drop the links between bi and bj themselves.
    int k2 = 0;
    for (int t = 0; t < k; t++)
      if (tmp[t].c != bi && tmp[t].c != bj) tmp[k2++] = tmp[t];
    if (k2 > cap[bi]) {
      nw_link_t *p = (nw_link_t *) realloc(lnk[bi], (size_t) k2 * sizeof(nw_link_t));
      if (!p) { ok = 0; break; }
      lnk[bi] = p;
      cap[bi] = k2;
    }
    for (int t = 0; t < k2; t++) lnk[bi][t] = tmp[t];
    len[bi] = k2;
    len[bj] = 0;

    // In every neighbour's list, bj becomes bi.
    for (int t = 0; t < k2; t++) {
      int c = lnk[bi][t].c;
      nw_link_t *L = lnk[c];
      int lc = len[c], pi = -1, pj = -1;
      for (int u = 0; u < lc; u++) {
        if (L[u].c == bi) pi = u;
        else if (L[u].c == bj) pj = u;
      }
      if (pj >= 0 && pi >= 0) {
        L[pi].e += L[pj].e;
        memmove(L + pj, L + pj + 1, (size_t) (lc - pj - 1) * sizeof(nw_link_t));
        len[c] = lc - 1;
      }
      else if (pj >= 0) {
        // Re-insert as bi, keeping the list sorted.
        double e = L[pj].e;
        memmove(L + pj, L + pj + 1, (size_t) (lc - pj - 1) * sizeof(nw_link_t));
        lc--;
        int u = lc;
        while (u > 0 && L[u - 1].c > bi) { L[u] = L[u - 1]; u--; }
        L[u].c = bi;
        L[u].e = e;
        len[c] = lc + 1;
      }
    }

    a[bi] += a[bj];
    a[bj]  = 0;

    mi[steps] = bi;
    mj[steps] = bj;
    q[steps + 1] = q[steps] + best;
    steps++;
  }

  for (int v = 0; v < n; v++) free(lnk[v]);
  free(tmp);
  if (!ok) error("Out of memory in C_network_greedy()."); // # nocov

  // The step with the largest Q (the first, on ties). As igraph's
  // cluster_fast_greedy() (and so NetCoMi), a single community of every
  // node is never returned: the step before it is.
  int kbest = 0;
  for (int s = 1; s <= steps; s++)
    if (q[s] > q[kbest]) kbest = s;
  if (kbest == n - 1 && kbest > 0) kbest--;

  // Replay the merges up to it: union by the surviving id.
  int *parent = (int *) R_alloc((size_t) n, sizeof(int));
  for (int v = 0; v < n; v++) parent[v] = v;
  for (int s = 0; s < kbest; s++) parent[mj[s]] = mi[s];
  for (int v = 0; v < n; v++) {
    int r = v;
    while (parent[r] != r) r = parent[r];
    mem[v] = r;
  }

  SEXP sq = PROTECT(allocVector(REALSXP, steps + 1));
  for (int t = 0; t <= steps; t++) REAL(sq)[t] = q[t];
  SET_VECTOR_ELT(res, 1, sq);

  UNPROTECT(3);
  return res;
}



//======================================================
// Barrat's weighted local clustering coefficient:
//
//   C_i = 1 / (s_i (k_i - 1)) sum_{j, h} (w_ij + w_ih) / 2 a_ij a_ih a_jh
//
// over ordered pairs of i's neighbours, with s_i its
// strength and k_i its degree. NaN for k_i < 2.
//
// R interface: n, from, to (0-based), w. Returns C (n).
//======================================================
SEXP C_network_barrat(SEXP sexp_n, SEXP sexp_from, SEXP sexp_to, SEXP sexp_w) {

  int n = asInteger(sexp_n);
  if (n < 1) error("Invalid arguments to C_network_barrat()."); // # nocov

  nw_graph_t g;
  nw_build(&g, n, sexp_from, sexp_to, sexp_w);

  SEXP    res  = PROTECT(allocVector(REALSXP, n));
  double *cc   = REAL(res);
  double *mark = (double *) R_alloc((size_t) n, sizeof(double));
  char   *is   = (char *)   R_alloc((size_t) n, 1);
  for (int v = 0; v < n; v++) is[v] = 0;

  for (int i = 0; i < n; i++) {

    int k = g.ptr[i + 1] - g.ptr[i];
    if (k < 2) { cc[i] = R_NaN; continue; }

    double s = 0;
    for (int x = g.ptr[i]; x < g.ptr[i + 1]; x++) {
      is[g.adj[x]]   = 1;
      mark[g.adj[x]] = g.val[x];
      s += g.val[x];
    }

    double sum = 0;
    for (int x = g.ptr[i]; x < g.ptr[i + 1]; x++) {
      int    j   = g.adj[x];
      double wij = g.val[x];
      for (int y = g.ptr[j]; y < g.ptr[j + 1]; y++) {
        int h = g.adj[y];
        if (is[h]) sum += (wij + mark[h]) / 2;
      }
    }

    for (int x = g.ptr[i]; x < g.ptr[i + 1]; x++) is[g.adj[x]] = 0;

    cc[i] = sum / (s * (k - 1));
  }

  UNPROTECT(1);
  return res;
}



//======================================================
// The leading eigenvector of the weighted adjacency
// matrix A of a connected graph, by block Lanczos.
// A is non-negative and irreducible, so its largest
// eigenvalue has a positive eigenvector (Perron); it is
// returned scaled to a largest value of 1, with its
// eigenvalue as attribute "value". n = 1 gives 1.
//
// R interface: n, from, to (0-based), w, cpus.
//======================================================
static void ne_mult(const double *V, double *W, int m, void *ctx) {
  const nw_graph_t *g = (const nw_graph_t *) ctx;
  int n = g->n;
  for (int c = 0; c < m; c++) {
    const double *v = V + (size_t) c * n;
    double       *w = W + (size_t) c * n;
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int x = g->ptr[i]; x < g->ptr[i + 1]; x++)
        s += g->val[x] * v[g->adj[x]];
      w[i] = s;
    }
  }
}

SEXP C_network_eigen(SEXP sexp_n, SEXP sexp_from, SEXP sexp_to, SEXP sexp_w, SEXP sexp_cpus) {

  int n    = asInteger(sexp_n);
  int cpus = asInteger(sexp_cpus);
  if (n < 1 || cpus < 1) error("Invalid arguments to C_network_eigen()."); // # nocov

  nw_graph_t g;
  nw_build(&g, n, sexp_from, sexp_to, sexp_w);

  SEXP    res = PROTECT(allocVector(REALSXP, n));
  double *ev  = REAL(res);
  double  val = 0;

  if (n == 1 || g.m == 0) {
    for (int i = 0; i < n; i++) ev[i] = 1;
  }
  else {

    // Scale for the convergence test: the largest strength,
    // an upper bound of the largest eigenvalue.
    double smax = 0;
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int x = g.ptr[i]; x < g.ptr[i + 1]; x++) s += g.val[x];
      if (s > smax) smax = s;
    }

    double vals[1], resid[1];
    lanczos_t lz = {
      .n = n, .k = 1, .dims = n, .deflate = 0, .scale = smax,
      .op = ne_mult, .ctx = &g, .n_threads = 1,
      .vecs = ev, .vals = vals, .resid = resid };
    lanczos(&lz);
    val = vals[0];

    if (!lz.converged)
      warning("The eigenvector centralities did not converge."); // # nocov

    // Perron: one sign. Take the sign of the largest
    // magnitude, then scale to a maximum of 1.
    double big = 0;
    for (int i = 0; i < n; i++) if (fabs(ev[i]) > fabs(big)) big = ev[i];
    for (int i = 0; i < n; i++) ev[i] /= big;
    for (int i = 0; i < n; i++) if (ev[i] < 0) ev[i] = 0;  // rounding only
  }

  SEXP sval = PROTECT(ScalarReal(val));
  setAttrib(res, install("value"), sval);
  UNPROTECT(2);
  return res;
}
