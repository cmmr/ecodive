// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Phylogenetic community structure: MPD and MNTD (Webb 2000; Webb et al.
// 2002), observed and against null models (NRI / NTI), for mpd(), mntd(),
// weighted_mpd(), weighted_mntd(), mpd_ses() and mntd_ses().
//
// The swap moves (swap_attempt()) are adapted from IndependentSwap() and
// TrialSwap() in phylocom 4.2's src/combase.c, restructured for ecodive: one
// presence map per thread, PCG draws, and each null matrix a chain of its own
// from the observed table. phylocom's notice:
//
//   phylocom Version 4.2 (c) 2009 Campbell O. Webb, David D. Ackerly,
//   Steven W. Kembel. All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   - Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   - Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   - The names of the authors may not be used to endorse or promote
//     products derived from this software without specific prior written
//     permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
// Distances on a tree, without a tips x tips matrix
// --------------------------------------------------
// ptree_init() numbers the nodes in preorder and keeps each node's depth
// (distance from the root) and subtree size. With the tips in preorder,
// u_0 .. u_{T-1}, the lowest common ancestor (LCA) of u_k and u_{k+1} is the
// parent of the node that follows u_k in preorder; the LCA of u_i and u_j
// (i < j) is the one of LCA(u_k, u_{k+1}), i <= k < j, nearest the root,
// i.e. with the smallest preorder number: a range minimum, answered in O(1)
// from a sparse table (Bender & Farach-Colton 2000).
//
// For one sample, ptree_vtree() sorts its S tips into preorder and adds the
// LCAs of neighbours: the "virtual tree" of at most 2S - 1 nodes, whose
// edges are paths of the tree with lengths = depth differences. On it:
//
//   MPD:  sum over edges of L_e w_e (W - w_e), with w_e the weight below
//         the edge and W the total, is the sum over unordered pairs of
//         w_i w_j d_ij. Unweighted (w = 1): divided by S(S - 1)/2.
//         Weighted: 2 x that / W^2, i.e. over all ordered pairs, i = j
//         included (phylocom's and picante's denominator).
//   MNTD: each tip's distance to the nearest other tip, from the nearest
//         below each node (best two, from different children) and then the
//         nearest outside each node, top-down.
//
// Each costs O(S log S), for the sort. With `taxa_dm` instead of a tree,
// the S x S distances are read from the packed triangle: O(S^2).
//
//
// Null models (null_model; picante's ses.mpd() names)
// ---------------------------------------------------
// Each draw r = 0 .. n_draws - 1 uses PCG streams of its own, so draws can
// run on any thread and give the same values:
//
//   taxa_labels       stream r: the features go to D distinct tips drawn
//                     at random from all T tips (the labels shuffled).
//   richness          stream r n + i for sample i: its S values go to S
//                     distinct features drawn from all D, in random order.
//   frequency         stream r: each feature's values go to that many
//                     distinct samples, drawn from all n.
//   sample_pool       as richness, from the P features present in a sample.
//   phylogeny_pool    as richness, from all T tips.
//   independent_swap  stream r: n_swaps successful checkerboard swaps.
//   trial_swap        stream r: n_swaps attempted swaps.
//
// Draws of "ordered k of N without replacement" are partial Fisher-Yates
// shuffles of a per-thread identity array, undone afterwards (O(k)).
// Swaps start each draw from the observed table (independent chains). The
// values stay with their feature, as in phylocom and picante.
//
// Per sample, the draws' values are summed in draw order (Welford) after
// each batch of draws, so results are the same for any number of threads.
//
//
// Exact moments (null_model = 'exact', MPD only)
// -----------------------------------------------
// For S tips drawn at random without replacement from all T (the
// phylogeny_pool and, per sample, taxa_labels distribution), with the
// sample's abundances x_1 .. x_S following its taxa: MPD = Z / c with
//   Z = sum_{i != j} a_ij d(t_i, t_j),  a_ij = 1 (c = S(S - 1)) or
//   a_ij = x_i x_j (c = W^2).
// With distances centred on their mean mu over pairs of distinct tips,
//   P2 = sum_{a != b} (d_ab - mu)^2,  r_a = sum_b (d_ab - mu) (b != a),
//   R2 = sum_a r_a^2,
// and the falling factorials (T)_k:
//   E[Z]   = mu sum a_ij
//   Var[Z] = A2 P2 / (T)_2 + A1 (R2 - P2) / (T)_3 + A0 (2 P2 - 4 R2) / (T)_4
// where A2 = 2 sum a_ij^2, A1 = 4 (sum_i rho_i^2 - sum a_ij^2) with
// rho_i = sum_j a_ij, and A0 = (sum a_ij)^2 - A2 - A1: the pairs of
// ordered pairs of sample positions sharing two, one or no positions. On a
// tree, sum_{a<b} d_ab, every r_a and sum_{a<b} d_ab^2 take one or two
// passes over the edges.
//
// Bender, M. A., & Farach-Colton, M. (2000). The LCA problem revisited.
// LATIN 2000, LNCS 1776, 88-94.

#include "ecodive.h"

#define STOP(...) do { free_all(); error(__VA_ARGS__); } while (0)

#define PS_NONE        0
#define PS_TAXA_LABELS 1
#define PS_RICHNESS    2
#define PS_FREQUENCY   3
#define PS_SAMPLE_POOL 4
#define PS_PHYLO_POOL  5
#define PS_IND_SWAP    6
#define PS_TRIAL_SWAP  7
#define PS_EXACT       8
#define PS_SUPPLIED    9


//======================================================
// The tree.
//======================================================

// Builds `pt` from the edges of a binary tree whose tips are 0 .. T - 1
// and internal nodes T .. 2T - 2 (C_philr_tree()'s numbering; any root).
// Memory is R_alloc'd. Returns 0, or 1 if the edges aren't such a tree.
int ptree_init (ptree_t *pt, int T, const int *epar, const int *echd,
                const double *elen) {

  int N = 2 * T - 1, E = N - 1;
  pt->T = T;
  pt->N = N;

  int    *par   = (int *)    R_alloc(N, sizeof(int));
  double *blen  = (double *) R_alloc(N, sizeof(double));
  int    *ch1   = (int *)    R_alloc(N, sizeof(int));
  int    *ch2   = (int *)    R_alloc(N, sizeof(int));
  for (int v = 0; v < N; v++) { par[v] = -1; ch1[v] = ch2[v] = -1; blen[v] = 0; }

  for (int e = 0; e < E; e++) {
    int p = epar[e], c = echd[e];
    if (p < T || p >= N || c < 0 || c >= N || c == p || par[c] >= 0) return 1;
    par[c]  = p;
    blen[c] = elen[e];
    if      (ch1[p] < 0) ch1[p] = c;
    else if (ch2[p] < 0) ch2[p] = c;
    else return 1;
  }

  int root = -1;
  for (int v = 0; v < N; v++)
    if (par[v] < 0) { if (root >= 0) return 1; root = v; }
  if (root < 0) return 1;

  // Preorder, first child first.
  int    *pre    = (int *)    R_alloc(N, sizeof(int));
  int    *at     = (int *)    R_alloc(N, sizeof(int));
  int    *stack  = (int *)    R_alloc(N, sizeof(int));
  int     top    = 0, k = 0;
  stack[top++] = root;
  while (top > 0) {
    int v = stack[--top];
    if (k >= N) return 1;
    pre[v]  = k;
    at[k++] = v;
    if (v >= T) {
      if (ch1[v] < 0 || ch2[v] < 0) return 1;
      stack[top++] = ch2[v];
      stack[top++] = ch1[v];
    }
  }
  if (k != N) return 1;

  double *dep  = (double *) R_alloc(N, sizeof(double)); // by preorder number
  int    *size = (int *)    R_alloc(N, sizeof(int));    // by preorder number
  int    *ppar = (int *)    R_alloc(N, sizeof(int));    // by preorder number
  double *plen = (double *) R_alloc(N, sizeof(double)); // by preorder number
  int    *ptip = (int *)    R_alloc(N, sizeof(int));    // tip id, or -1
  int    *tpos = (int *)    R_alloc(T, sizeof(int));    // tip id -> rank among tips
  int    *tpre = (int *)    R_alloc(T, sizeof(int));    // rank -> preorder number

  for (int q = 0; q < N; q++) {
    int v   = at[q];
    ppar[q] = par[v] < 0 ? -1 : pre[par[v]];
    plen[q] = blen[v];
    dep[q]  = ppar[q] < 0 ? 0 : dep[ppar[q]] + plen[q];
    size[q] = 1;
    ptip[q] = v < T ? v : -1;
  }
  for (int q = N - 1; q > 0; q--) size[ppar[q]] += size[q];

  int t = 0;
  for (int q = 0; q < N; q++)
    if (ptip[q] >= 0) { tpos[ptip[q]] = t; tpre[t++] = q; }

  // LCA of neighbouring tips, and the sparse table of range minima.
  int M = T - 1, levels = 1;
  while ((1 << levels) <= M) levels++;
  int *rmq = (int *) R_alloc((size_t) levels * (M > 0 ? M : 1), sizeof(int));
  for (int j = 0; j < M; j++) rmq[j] = ppar[tpre[j] + 1];
  for (int l = 1; l < levels; l++) {
    int h = 1 << (l - 1);
    int *a = rmq + (size_t) (l - 1) * M, *b = rmq + (size_t) l * M;
    for (int j = 0; j + 2 * h <= M; j++) b[j] = a[j] < a[j + h] ? a[j] : a[j + h];
  }

  pt->dep = dep; pt->size = size; pt->par = ppar; pt->len = plen;
  pt->tip = ptip; pt->tpos = tpos; pt->tpre = tpre;
  int *lg = (int *) R_alloc((size_t) T + 1, sizeof(int));
  lg[0] = lg[1] = 0;
  for (int j = 2; j <= T; j++) lg[j] = lg[j / 2] + 1;

  pt->rmq = rmq; pt->levels = levels; pt->lg = lg;
  return 0;
}


// Preorder number of the LCA of the tips ranked a < b among the tips.
static inline int ptree_lca (const ptree_t *pt, int a, int b) {
  int M = pt->T - 1, l = pt->lg[b - a];
  const int *r = pt->rmq + (size_t) l * M;
  int x = r[a], y = r[b - (1 << l)];
  return x < y ? x : y;
}


//======================================================
// Per-thread workspace for one sample's virtual tree.
//======================================================

void ptree_ws_init (ptree_ws_t *ws, int cap) {
  int m = 2 * cap + 1;
  ws->cap  = cap;
  ws->key  = (int *)    R_alloc(m, sizeof(int));
  ws->kw   = (double *) R_alloc(m, sizeof(double));
  ws->node = (int *)    R_alloc(m, sizeof(int));
  ws->par  = (int *)    R_alloc(m, sizeof(int));
  ws->len  = (double *) R_alloc(m, sizeof(double));
  ws->w    = (double *) R_alloc(m, sizeof(double));
  ws->b1   = (double *) R_alloc(m, sizeof(double));
  ws->b2   = (double *) R_alloc(m, sizeof(double));
  ws->c1   = (int *)    R_alloc(m, sizeof(int));
  ws->stk  = (int *)    R_alloc(m, sizeof(int));
  ws->leaf = (int *)    R_alloc(m, sizeof(int));
  ws->tq   = (int *)    R_alloc(m, sizeof(int));
  ws->tp   = (int *)    R_alloc(m, sizeof(int));
  ws->tpo  = (int *)    R_alloc(m, sizeof(int));
  ws->tw   = (double *) R_alloc(m, sizeof(double));
}


// Sorts (key, weight) pairs by key, keys in 0 .. 2^31 - 1: insertion sort
// for a few, else LSD radix sort on 8-bit digits (as many as the largest
// key needs), through the buffers tk and tw.
static void sort_kw (int *k, double *w, int n, int *tk, double *tw) {
  if (n < 48) {
    for (int i = 1; i < n; i++) {
      int kk = k[i]; double ww = w[i]; int j = i - 1;
      while (j >= 0 && k[j] > kk) { k[j + 1] = k[j]; w[j + 1] = w[j]; j--; }
      k[j + 1] = kk; w[j + 1] = ww;
    }
    return;
  }
  int mx = 0;
  for (int i = 0; i < n; i++) if (k[i] > mx) mx = k[i];
  int *sk = k, *dk = tk; double *sw = w, *dw = tw;
  for (int shift = 0; shift < 32 && (mx >> shift) > 0; shift += 8) {
    int cnt[257];
    memset(cnt, 0, sizeof(cnt));
    for (int i = 0; i < n; i++) cnt[((sk[i] >> shift) & 255) + 1]++;
    for (int d = 0; d < 256; d++) cnt[d + 1] += cnt[d];
    for (int i = 0; i < n; i++) {
      int q = cnt[(sk[i] >> shift) & 255]++;
      dk[q] = sk[i]; dw[q] = sw[i];
    }
    int *t = sk; sk = dk; dk = t;
    double *u = sw; sw = dw; dw = u;
  }
  if (sk != k) {
    memcpy(k, sk, (size_t) n * sizeof(int));
    memcpy(w, sw, (size_t) n * sizeof(double));
  }
}


// The virtual tree of the S tips `tips` (distinct tip ids) with weights
// `wt` (NULL: 1 each). On return, its m nodes are numbered 0 .. m - 1 with
// every parent before its children (node 0 is the root): ws->node (their
// preorder numbers in the tree), ws->par (parent's number; -1 for the
// root), ws->len (edge length to the parent), ws->w (the tips' weights; 0
// at the other nodes), and ws->leaf[0 .. S - 1] (the tips' numbers).
// Returns m (at most 2S - 1).
//
// Built with a stack holding the path from the root to the latest tip
// (tips in preorder): each new tip's LCA with the previous one is where
// the path branches. Nodes leave the stack after all their descendants,
// so the order they leave in is a postorder, reversed for the numbering.
int ptree_vtree (const ptree_t *pt, ptree_ws_t *ws, int S,
                 const int *tips, const double *wt) {

  int    *key = ws->key, *stk = ws->stk;
  int    *tq  = ws->tq, *tp = ws->tp, *tpo = ws->tpo;
  double *kw  = ws->kw, *tw = ws->tw;

  for (int i = 0; i < S; i++) { key[i] = pt->tpos[tips[i]]; kw[i] = wt ? wt[i] : 1; }
  sort_kw(key, kw, S, tp, tw);

  // Temporary ids in order of creation: tq (preorder number), tw
  // (weight), tp (parent's id); tpo lists ids as they leave the stack.
  int m = 0, top = 0, np = 0;
  #define LEAVE(x, y) do { tp[x] = (y); tpo[np++] = (x); } while (0)

  tq[m] = pt->tpre[key[0]]; tw[m] = kw[0]; stk[top++] = m++;
  for (int k = 1; k < S; k++) {
    int lq = ptree_lca(pt, key[k - 1], key[k]);
    while (top >= 2 && tq[stk[top - 2]] >= lq) {
      LEAVE(stk[top - 1], stk[top - 2]);
      top--;
    }
    if (tq[stk[top - 1]] != lq) {
      tq[m] = lq; tw[m] = 0;
      LEAVE(stk[top - 1], m);
      stk[top - 1] = m++;
    }
    tq[m] = pt->tpre[key[k]]; tw[m] = kw[k]; stk[top++] = m++;
  }
  while (top >= 2) { LEAVE(stk[top - 1], stk[top - 2]); top--; }
  tp[stk[0]] = -1;
  tpo[np++]  = stk[0];
  #undef LEAVE

  // Renumber: the j-th to leave is m - 1 - j.
  int *nw = stk; // free now; at least m long
  for (int j = 0; j < m; j++) nw[tpo[j]] = m - 1 - j;
  int t = 0;
  for (int j = m - 1; j >= 0; j--) {
    int x = tpo[j], i = m - 1 - j, q = tq[x];
    ws->node[i] = q;
    ws->par[i]  = tp[x] < 0 ? -1 : nw[tp[x]];
    ws->len[i]  = tp[x] < 0 ? 0 : pt->dep[q] - pt->dep[tq[tp[x]]];
    ws->w[i]    = tw[x];
    if (pt->tip[q] >= 0) ws->leaf[t++] = i;
  }
  return m;
}


// Sum over unordered pairs of tips of w_i w_j d_ij, on the virtual tree
// (whose ws->w is overwritten with the weights below each node).
double ptree_pair_sum (ptree_ws_t *ws, int m) {
  double W = 0;
  for (int i = 0; i < m; i++) W += ws->w[i];
  for (int i = m - 1; i > 0; i--) ws->w[ws->par[i]] += ws->w[i];
  double s = 0;
  for (int i = 1; i < m; i++) s += ws->len[i] * ws->w[i] * (W - ws->w[i]);
  return s;
}


//======================================================
// Shared state for the worker threads.
//======================================================
typedef struct {
  // Table (observed), samples x features, CSR.
  int     n, D, T, nnz;
  const int    *rpos, *ridx;
  const double *rval;
  // Table by feature (CSC), for `frequency`.
  const int    *cpos, *cidx;
  const double *cval;
  // Supplied null tables: R tables of n rows each, CSR.
  const int    *spos, *sidx;
  const double *sval;
  int     smax;
  // Distances.
  const ptree_t *pt;       // tree, or NULL
  const double  *dm;       // packed T x T triangle, or NULL
  // Settings.
  int     stat, weighted, model, n_swaps, cap;
  uint64_t seed;
  const int *pool; int P;  // sample_pool's features
  int     nested;          // swaps: no checkerboard exists
  // Batch.
  int     r0, nb;          // draws r0 .. r0 + nb - 1
  double *vals;            // nb x n
  // Per-thread scratch.
  ptree_ws_t *ws;
  int    **arr;            // identity arrays (length max(T, n))
  int    **rec;            // undo records / drawn indices
  int    **tips;  double **tw;
  int    **npos;  int **nidx; double **nval;   // frequency's null table
  int    **smap;  int **sfeat; double **sval2; // swaps
} ps_ctx_t;

static ps_ctx_t ps;


//======================================================
// MPD or MNTD of one community.
//======================================================
static double score (int thread, int S, const int *tips, const double *w) {

  if (S < 2) return NA_REAL;
  int weighted = ps.weighted;

  if (ps.pt) {
    ptree_ws_t *ws = ps.ws + thread;
    int m = ptree_vtree(ps.pt, ws, S, tips, weighted ? w : NULL);

    if (ps.stat == 0) {
      double W = 0;
      for (int i = 0; i < m; i++) W += ws->w[i];
      double s = ptree_pair_sum(ws, m);
      return weighted ? 2 * s / (W * W) : s / ((double) S * (S - 1) / 2);
    }

    // MNTD: nearest tip below each node (best two, from different
    // children), then nearest outside each node, top-down.
    double *b1 = ws->b1, *b2 = ws->b2, *len = ws->len;
    int    *c1 = ws->c1, *par = ws->par;
    for (int i = 0; i < m; i++) { b1[i] = b2[i] = R_PosInf; c1[i] = -1; }
    for (int t = 0; t < S; t++) b1[ws->leaf[t]] = 0;
    for (int i = m - 1; i > 0; i--) {
      double c = b1[i] + len[i];
      int    p = par[i];
      if      (c < b1[p]) { b2[p] = b1[p]; b1[p] = c; c1[p] = i; }
      else if (c < b2[p])   b2[p] = c;
    }
    // Nearest outside each node, into ws->kw (free after
    // ptree_vtree(), and at least m long).
    double *up = ws->kw;
    up[0] = R_PosInf;
    for (int i = 1; i < m; i++) {
      int    p   = par[i];
      double sib = c1[p] == i ? b2[p] : b1[p];
      double u   = up[p] < sib ? up[p] : sib;
      up[i] = len[i] + u;
    }
    double s = 0, W = 0;
    for (int t = 0; t < S; t++) {
      int    i = ws->leaf[t];
      double x = weighted ? ws->w[i] : 1;
      s += x * up[i];
      W += x;
    }
    return s / W;
  }

  // taxa_dm: the packed lower triangle over T labels.
  const double *dm = ps.dm;
  size_t T = (size_t) ps.T;
  #define DM(a, b) ((a) > (b) ? \
    dm[(size_t) (b) * T - (size_t) (b) * ((b) + 1) / 2 + (size_t) ((a) - (b) - 1)] : \
    dm[(size_t) (a) * T - (size_t) (a) * ((a) + 1) / 2 + (size_t) ((b) - (a) - 1)])

  if (ps.stat == 0) {
    double s = 0, W = 0;
    for (int i = 0; i < S; i++) {
      double xi = weighted ? w[i] : 1;
      W += xi;
      for (int j = 0; j < i; j++) {
        double xj = weighted ? w[j] : 1;
        s += xi * xj * DM(tips[i], tips[j]);
      }
    }
    return weighted ? 2 * s / (W * W) : s / ((double) S * (S - 1) / 2);
  }

  double s = 0, W = 0;
  for (int i = 0; i < S; i++) {
    double mn = R_PosInf;
    for (int j = 0; j < S; j++) {
      if (j == i) continue;
      double d = DM(tips[i], tips[j]);
      if (d < mn) mn = d;
    }
    double xi = weighted ? w[i] : 1;
    s += xi * mn;
    W += xi;
  }
  return s / W;
  #undef DM
}


//======================================================
// Ordered draw of k of 0 .. N - 1 without replacement,
// into out[], with arr (identity, length >= N) restored.
//======================================================
static void draw_k (pcg32_random_t *rng, int *arr, int *rec, int N, int k, int *out) {
  for (int t = 0; t < k; t++) {
    int j = t + (int) pcg32_bounded(rng, (uint32_t) (N - t));
    int x = arr[j]; arr[j] = arr[t]; arr[t] = x;
    rec[t] = j;
    out[t] = x;
  }
  for (int t = k - 1; t >= 0; t--) {
    int j = rec[t];
    int x = arr[j]; arr[j] = arr[t]; arr[t] = x;
  }
}


//======================================================
// One swap attempt (adapted from phylocom; see the
// header). Returns 1 if a checkerboard was swapped.
//======================================================
static int swap_attempt (pcg32_random_t *rng, int *map, int *feat, double *val) {

  int n = ps.n, D = ps.D;
  int i1 = (int) pcg32_bounded(rng, (uint32_t) n);
  int i2 = (int) pcg32_bounded(rng, (uint32_t) (n - 1)); if (i2 >= i1) i2++;
  int f1 = (int) pcg32_bounded(rng, (uint32_t) D);
  int f2 = (int) pcg32_bounded(rng, (uint32_t) (D - 1)); if (f2 >= f1) f2++;

  size_t a = (size_t) i1 * D + f1, b = (size_t) i1 * D + f2;
  size_t c = (size_t) i2 * D + f1, d = (size_t) i2 * D + f2;

  // (i1, f1) and (i2, f2) present, the others absent, or the reverse.
  // A feature's values stay with the feature: each moves to the other
  // sample.
  if (map[a] >= 0 && map[d] >= 0 && map[b] < 0 && map[c] < 0) {
    int ea = map[a], ed = map[d];
    feat[ea] = f2; feat[ed] = f1;
    double t = val[ea]; val[ea] = val[ed]; val[ed] = t;
    map[a] = -1; map[d] = -1; map[b] = ea; map[c] = ed;
    return 1;
  }
  if (map[b] >= 0 && map[c] >= 0 && map[a] < 0 && map[d] < 0) {
    int eb = map[b], ec = map[c];
    feat[eb] = f1; feat[ec] = f2;
    double t = val[eb]; val[eb] = val[ec]; val[ec] = t;
    map[b] = -1; map[c] = -1; map[a] = eb; map[d] = ec;
    return 1;
  }
  return 0;
}


//======================================================
// Worker: draws r0 .. r0 + nb - 1, statically shared.
//======================================================
static void *null_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = ps.n, D = ps.D, T = ps.T;
  int    *arr  = ps.arr[thread_i];
  int    *rec  = ps.rec[thread_i];
  int    *tips = ps.tips[thread_i];
  double *tw   = ps.tw[thread_i];

  for (int b = thread_i; b < ps.nb; b += n_threads) {

    int     r   = ps.r0 + b;
    double *out = ps.vals + (size_t) b * n;
    pcg32_random_t rng;

    switch (ps.model) {

      case PS_TAXA_LABELS: {
        pcg32_srandom_r(&rng, ps.seed, (uint64_t) r);
        int *map = ps.npos[thread_i]; // D: each feature's tip
        draw_k(&rng, arr, rec, T, D, map);
        for (int i = 0; i < n; i++) {
          int S = 0;
          for (int e = ps.rpos[i]; e < ps.rpos[i + 1]; e++) {
            tips[S] = map[ps.ridx[e]]; tw[S++] = ps.rval[e];
          }
          out[i] = score(thread_i, S, tips, tw);
        }
        break;
      }

      case PS_RICHNESS:
      case PS_SAMPLE_POOL:
      case PS_PHYLO_POOL: {
        int N = ps.model == PS_RICHNESS ? D : ps.model == PS_SAMPLE_POOL ? ps.P : T;
        for (int i = 0; i < n; i++) {
          pcg32_srandom_r(&rng, ps.seed, (uint64_t) r * n + i);
          int S = ps.rpos[i + 1] - ps.rpos[i];
          draw_k(&rng, arr, rec, N, S, tips);
          for (int t = 0; t < S; t++) {
            if (ps.model == PS_SAMPLE_POOL) tips[t] = ps.pool[tips[t]];
            tw[t] = ps.rval[ps.rpos[i] + t];
          }
          out[i] = score(thread_i, S, tips, tw);
        }
        break;
      }

      case PS_FREQUENCY: {
        pcg32_srandom_r(&rng, ps.seed, (uint64_t) r);
        int *cnt = ps.npos[thread_i], *idx = ps.nidx[thread_i];
        double *nv = ps.nval[thread_i];
        int *who = ps.smap[thread_i]; // nnz: each entry's new sample
        memset(cnt, 0, ((size_t) n + 1) * sizeof(int));
        for (int f = 0; f < D; f++) {
          int k = ps.cpos[f + 1] - ps.cpos[f];
          if (k == 0) continue;
          draw_k(&rng, arr, rec, n, k, who + ps.cpos[f]);
          for (int e = ps.cpos[f]; e < ps.cpos[f + 1]; e++) cnt[who[e] + 1]++;
        }
        for (int i = 0; i < n; i++) cnt[i + 1] += cnt[i];
        // Fill by feature, so each sample's entries are in feature order.
        int *fill = ps.sfeat[thread_i];
        memcpy(fill, cnt, (size_t) n * sizeof(int));
        for (int f = 0; f < D; f++)
          for (int e = ps.cpos[f]; e < ps.cpos[f + 1]; e++) {
            int q = fill[who[e]]++;
            idx[q] = f; nv[q] = ps.cval[e];
          }
        for (int i = 0; i < n; i++)
          out[i] = score(thread_i, cnt[i + 1] - cnt[i], idx + cnt[i], nv + cnt[i]);
        break;
      }

      case PS_IND_SWAP:
      case PS_TRIAL_SWAP: {
        pcg32_srandom_r(&rng, ps.seed, (uint64_t) r);
        int    *map  = ps.smap[thread_i];
        int    *feat = ps.sfeat[thread_i];
        double *val  = ps.sval2[thread_i];
        memcpy(feat, ps.ridx, (size_t) ps.nnz * sizeof(int));
        memcpy(val,  ps.rval, (size_t) ps.nnz * sizeof(double));
        if (!ps.nested) {
          for (int i = 0; i < n; i++)
            for (int e = ps.rpos[i]; e < ps.rpos[i + 1]; e++)
              map[(size_t) i * D + feat[e]] = e;
          if (ps.model == PS_IND_SWAP) {
            for (int s = 0; s < ps.n_swaps; )
              s += swap_attempt(&rng, map, feat, val);
          }
          else {
            for (int s = 0; s < ps.n_swaps; s++)
              swap_attempt(&rng, map, feat, val);
          }
          for (int i = 0; i < n; i++)
            for (int e = ps.rpos[i]; e < ps.rpos[i + 1]; e++)
              map[(size_t) i * D + feat[e]] = -1;
        }
        for (int i = 0; i < n; i++) {
          int e0 = ps.rpos[i];
          out[i] = score(thread_i, ps.rpos[i + 1] - e0, feat + e0, val + e0);
        }
        break;
      }

      case PS_SUPPLIED: {
        for (int i = 0; i < n; i++) {
          size_t row = (size_t) r * n + i;
          int e0 = ps.spos[row];
          out[i] = score(thread_i, ps.spos[row + 1] - e0, ps.sidx + e0, ps.sval + e0);
        }
        break;
      }

      default: break;
    }
  }

  return NULL;
}


// Observed values, samples statically shared.
static void *obs_worker (void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  for (int i = thread_i; i < ps.n; i += n_threads) {
    int e0 = ps.rpos[i];
    ps.vals[i] = score(thread_i, ps.rpos[i + 1] - e0, ps.ridx + e0, ps.rval + e0);
  }
  return NULL;
}


//======================================================
// Exact moments of MPD over random tips (see header).
// Sets mu (mean distance between distinct tips), P2 and
// R2 (centred), as long doubles.
//======================================================
static void exact_sums (long double *mu_, long double *P2_, long double *R2_) {

  int T = ps.T;
  long double Tl = T, sum = 0, sq = 0;
  long double *r = (long double *) R_alloc(T, sizeof(long double));

  if (ps.pt) {
    const ptree_t *pt = ps.pt;
    int N = pt->N;
    long double *cnt = (long double *) R_alloc(N, sizeof(long double));
    long double *s1  = (long double *) R_alloc(N, sizeof(long double));
    long double *s2  = (long double *) R_alloc(N, sizeof(long double));
    long double *R   = (long double *) R_alloc(N, sizeof(long double));
    for (int q = 0; q < N; q++) { cnt[q] = pt->tip[q] >= 0; s1[q] = s2[q] = 0; }

    // Postorder: merge children into parents, counting cross pairs.
    for (int q = N - 1; q > 0; q--) {
      int p = pt->par[q];
      long double L = pt->len[q];
      long double c = cnt[q], a1 = s1[q] + c * L, a2 = s2[q] + 2 * L * s1[q] + c * L * L;
      sum += L * c * (Tl - c);
      // Cross pairs with the children of p merged so far.
      sq  += cnt[p] * a2 + 2 * s1[p] * a1 + c * s2[p];
      cnt[p] += c; s1[p] += a1; s2[p] += a2;
    }
    // Distances from every node to all tips: down, then up.
    R[0] = s1[0];
    for (int q = 1; q < N; q++) {
      long double c = cnt[q]; // tips below q
      R[q] = R[pt->par[q]] + (long double) pt->len[q] * (Tl - 2 * c);
    }
    for (int q = 0; q < N; q++) if (pt->tip[q] >= 0) r[pt->tip[q]] = R[q];
  }
  else {
    for (int a = 0; a < T; a++) r[a] = 0;
    size_t k = 0;
    for (int a = 0; a < T; a++)
      for (int b = a + 1; b < T; b++, k++) {
        long double d = ps.dm[k];
        sum += d; sq += d * d;
        r[a] += d; r[b] += d;
      }
  }

  long double pairs = Tl * (Tl - 1);         // ordered pairs
  long double mu    = 2 * sum / pairs;
  long double P2    = 2 * sq - pairs * mu * mu;
  long double R2    = 0;
  for (int a = 0; a < T; a++) {
    long double x = r[a] - (Tl - 1) * mu;
    R2 += x * x;
  }
  *mu_ = mu; *P2_ = P2; *R2_ = R2;
}


static void exact_moments (int S, const double *x, long double mu, long double P2,
                           long double R2, double *mean, double *sd) {

  long double Tl = ps.T, Sl = S;
  long double suma, sumsq, sumrho2, c;

  if (S < 2) { *mean = *sd = NA_REAL; return; }

  if (ps.weighted) {
    long double W = 0, x2 = 0, x4 = 0;
    for (int i = 0; i < S; i++) { long double v = x[i]; W += v; x2 += v * v; x4 += v * v * v * v; }
    suma  = W * W - x2;
    sumsq = x2 * x2 - x4;
    sumrho2 = 0;
    for (int i = 0; i < S; i++) { long double rho = x[i] * (W - x[i]); sumrho2 += rho * rho; }
    c = W * W;
  }
  else {
    suma = sumsq = Sl * (Sl - 1);
    sumrho2 = Sl * (Sl - 1) * (Sl - 1);
    c = Sl * (Sl - 1);
  }

  long double A2 = 2 * sumsq;
  long double A1 = 4 * (sumrho2 - sumsq);
  long double A0 = suma * suma - A2 - A1;
  long double T2 = Tl * (Tl - 1), T3 = T2 * (Tl - 2), T4 = T3 * (Tl - 3);
  long double v  = A2 * P2 / T2;
  if (S >= 3) v += A1 * (R2 - P2) / T3;
  if (S >= 4) v += A0 * (2 * P2 - 4 * R2) / T4;
  if (v < 0) v = 0;

  *mean = (double) (mu * suma / c);
  *sd   = (double) (sqrtl(v) / c);
}


//======================================================
// R interface.
//
// counts     samples x D features (features are tips
//            0 .. D - 1 of the tree / rows of dm).
// epar, echd, elen: a binary tree over T tips (ids
//            0 .. T - 1), or NULL with dm (packed, T).
// stat       0 MPD, 1 MNTD.
// model      PS_* code; n_draws, n_swaps, seed.
// supplied   with PS_SUPPLIED: a table of n_draws x n
//            rows (draw-major) and T columns.
//
// Returns list(n_taxa, obs, mean, sd, rank, n_used,
//              nested).
//======================================================
SEXP C_phylo_structure (
    SEXP sexp_counts, SEXP sexp_margin,
    SEXP sexp_epar, SEXP sexp_echd, SEXP sexp_elen, SEXP sexp_dm, SEXP sexp_T,
    SEXP sexp_stat, SEXP sexp_weighted, SEXP sexp_model, SEXP sexp_n_draws,
    SEXP sexp_n_swaps, SEXP sexp_seed, SEXP sexp_cpus, SEXP sexp_supplied ) {

  int T         = asInteger(sexp_T);
  int stat      = asInteger(sexp_stat);
  int weighted  = asLogical(sexp_weighted);
  int model     = asInteger(sexp_model);
  int n_draws   = asInteger(sexp_n_draws);
  int n_swaps   = asInteger(sexp_n_swaps);
  int n_threads = asInteger(sexp_cpus);
  double seed_d = asReal(sexp_seed);
  int tree      = !isNull(sexp_epar);

  if (T < 2 || stat < 0 || stat > 1 || weighted == NA_LOGICAL || model < 0 ||
      model > PS_SUPPLIED || n_draws < 0 || n_swaps < 0 || n_threads < 1 ||
      !R_FINITE(seed_d))
    error("Invalid arguments to C_phylo_structure()."); // # nocov
  if (tree && (XLENGTH(sexp_epar) != 2 * (R_xlen_t) T - 2 ||
               XLENGTH(sexp_echd) != XLENGTH(sexp_epar) ||
               XLENGTH(sexp_elen) != XLENGTH(sexp_epar)))
    error("Invalid arguments to C_phylo_structure()."); // # nocov
  if (!tree && XLENGTH(sexp_dm) != (R_xlen_t) T * (T - 1) / 2)
    error("Invalid arguments to C_phylo_structure()."); // # nocov
  if (model == PS_NONE || model == PS_EXACT) n_draws = 0;

  // The observed table, CSR, non-zero values only.
  init_n_ptrs(20);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int n = em->n_samples, D = em->n_otus;
  if (D > T) STOP("Invalid arguments to C_phylo_structure()."); // # nocov

  int    *rpos = (int *) R_alloc((size_t) n + 1, sizeof(int));
  int     nnz  = 0;
  for (int e = 0; e < em->nnz; e++) {
    double x = em->val_vec[e];
    if (ISNAN(x) || !R_FINITE(x)) STOP("`counts` must not contain NA, NaN, or infinite values.");
    if (x < 0)                    STOP("`counts` must not contain negative values.");
    if (x > 0) nnz++;
  }
  int    *ridx = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  double *rval = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  {
    int k = 0;
    rpos[0] = 0;
    for (int i = 0; i < n; i++) {
      for (int e = em->pos_vec[i]; e < em->pos_vec[i + 1]; e++)
        if (em->val_vec[e] > 0) { ridx[k] = em->otu_vec[e]; rval[k++] = em->val_vec[e]; }
      rpos[i + 1] = k;
    }
  }
  free_all();

  // Supplied null tables.
  int *spos = NULL, *sidx = NULL, smax = 0;
  double *sval = NULL;
  if (model == PS_SUPPLIED) {
    SEXP sexp_one = PROTECT(ScalarInteger(1));
    init_n_ptrs(20);
    ecomatrix_t *sm = new_ecomatrix(sexp_supplied, sexp_one);
    if (sm->n_otus != T || (R_xlen_t) sm->n_samples != (R_xlen_t) n * n_draws)
      STOP("Invalid arguments to C_phylo_structure()."); // # nocov
    size_t rows = (size_t) sm->n_samples;
    spos = (int *) R_alloc(rows + 1, sizeof(int));
    int k = 0;
    for (int e = 0; e < sm->nnz; e++) {
      double x = sm->val_vec[e];
      if (!R_FINITE(x) || x < 0) STOP("Supplied null tables must be finite and non-negative.");
      if (x > 0) k++;
    }
    sidx = (int *)    R_alloc((size_t) k + 1, sizeof(int));
    sval = (double *) R_alloc((size_t) k + 1, sizeof(double));
    k = 0; spos[0] = 0;
    for (size_t i = 0; i < rows; i++) {
      for (int e = sm->pos_vec[i]; e < sm->pos_vec[i + 1]; e++)
        if (sm->val_vec[e] > 0) { sidx[k] = sm->otu_vec[e]; sval[k++] = sm->val_vec[e]; }
      spos[i + 1] = k;
      if (k - spos[i] > smax) smax = k - spos[i];
    }
    free_all();
    UNPROTECT(1);
  }

  // Distances.
  ptree_t pt;
  if (tree) {
    if (ptree_init(&pt, T, INTEGER(sexp_epar), INTEGER(sexp_echd), REAL(sexp_elen)))
      error("Invalid arguments to C_phylo_structure()."); // # nocov
  }

  memset(&ps, 0, sizeof(ps));
  ps.n = n; ps.D = D; ps.T = T; ps.nnz = nnz;
  ps.rpos = rpos; ps.ridx = ridx; ps.rval = rval;
  ps.spos = spos; ps.sidx = sidx; ps.sval = sval; ps.smax = smax;
  ps.pt = tree ? &pt : NULL;
  ps.dm = tree ? NULL : REAL(sexp_dm);
  ps.stat = stat; ps.weighted = weighted; ps.model = model; ps.n_swaps = n_swaps;
  ps.seed = (uint64_t) (int64_t) seed_d;

  // Richness: observed, and the most any community can have.
  int *n_taxa = (int *) R_alloc(n > 0 ? n : 1, sizeof(int));
  int  smax_obs = 0;
  for (int i = 0; i < n; i++) {
    n_taxa[i] = rpos[i + 1] - rpos[i];
    if (n_taxa[i] > smax_obs) smax_obs = n_taxa[i];
  }
  int cap = smax_obs;
  if (model == PS_FREQUENCY) cap = D;
  if (model == PS_SUPPLIED && smax > cap) cap = smax;
  if (cap < 1) cap = 1;
  ps.cap = cap;

  // sample_pool: the features present in any sample.
  if (model == PS_SAMPLE_POOL) {
    int *seen = (int *) R_alloc(D > 0 ? D : 1, sizeof(int));
    memset(seen, 0, (size_t) D * sizeof(int));
    for (int e = 0; e < nnz; e++) seen[ridx[e]] = 1;
    int *pool = (int *) R_alloc(D > 0 ? D : 1, sizeof(int));
    int P = 0;
    for (int f = 0; f < D; f++) if (seen[f]) pool[P++] = f;
    ps.pool = pool; ps.P = P;
  }

  // frequency: the table by feature.
  if (model == PS_FREQUENCY) {
    int    *cpos = (int *)    R_alloc((size_t) D + 1, sizeof(int));
    int    *cidx = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
    double *cval = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
    memset(cpos, 0, ((size_t) D + 1) * sizeof(int));
    for (int e = 0; e < nnz; e++) cpos[ridx[e] + 1]++;
    for (int f = 0; f < D; f++) cpos[f + 1] += cpos[f];
    int *fill = (int *) R_alloc((size_t) D + 1, sizeof(int));
    memcpy(fill, cpos, (size_t) D * sizeof(int));
    for (int i = 0; i < n; i++)
      for (int e = rpos[i]; e < rpos[i + 1]; e++) {
        int q = fill[ridx[e]]++;
        cidx[q] = i; cval[q] = rval[e];
      }
    ps.cpos = cpos; ps.cidx = cidx; ps.cval = cval;
  }

  // Swaps: none possible when every pair of samples is nested (then the
  // observed table is the only one with its margins). Samples ranked by
  // richness, decreasing: nested iff each feature is in exactly the
  // first k of them, k its frequency.
  if (model == PS_IND_SWAP || model == PS_TRIAL_SWAP) {
    int nested = (n < 2 || D < 2);
    if (!nested) {
      int *ord  = (int *) R_alloc(n, sizeof(int));
      int *rnk  = (int *) R_alloc(n, sizeof(int));
      int *freq = (int *) R_alloc(D, sizeof(int));
      int *mx   = (int *) R_alloc(D, sizeof(int));
      // Counting sort by richness, decreasing.
      int *bucket = (int *) R_alloc((size_t) D + 2, sizeof(int));
      memset(bucket, 0, ((size_t) D + 2) * sizeof(int));
      for (int i = 0; i < n; i++) bucket[D - n_taxa[i] + 1]++;
      for (int k = 0; k <= D; k++) bucket[k + 1] += bucket[k];
      for (int i = 0; i < n; i++) ord[bucket[D - n_taxa[i]]++] = i;
      for (int k = 0; k < n; k++) rnk[ord[k]] = k;
      memset(freq, 0, (size_t) D * sizeof(int));
      for (int f = 0; f < D; f++) mx[f] = -1;
      for (int i = 0; i < n; i++)
        for (int e = rpos[i]; e < rpos[i + 1]; e++) {
          int f = ridx[e];
          freq[f]++;
          if (rnk[i] > mx[f]) mx[f] = rnk[i];
        }
      nested = 1;
      for (int f = 0; f < D && nested; f++) if (freq[f] && mx[f] != freq[f] - 1) nested = 0;
    }
    ps.nested = nested;
  }


  //--------------------------------------------------
  // Threads and their scratch.
  //--------------------------------------------------
  int nt = n_threads;
  if (n_draws == 0) nt = threads_for(n_threads, (double) nnz * 20 + n, 2e5);
  else              nt = threads_for(n_threads, (double) n_draws * ((double) nnz * 20 + n), 2e5);
  if (nt < 1) nt = 1;

  // The swaps' samples x features maps: fewer threads rather than more
  // than 1 GB of them (results are the same for any number).
  if (n_draws > 0 && (model == PS_IND_SWAP || model == PS_TRIAL_SWAP) && !ps.nested) {
    double bytes = (double) n * D * sizeof(int);
    while (nt > 1 && bytes * nt > 1073741824.0) nt--;
  }

  ps.ws   = (ptree_ws_t *) R_alloc(nt, sizeof(ptree_ws_t));
  ps.arr  = (int **)    R_alloc(nt, sizeof(int *));
  ps.rec  = (int **)    R_alloc(nt, sizeof(int *));
  ps.tips = (int **)    R_alloc(nt, sizeof(int *));
  ps.tw   = (double **) R_alloc(nt, sizeof(double *));
  ps.npos = (int **)    R_alloc(nt, sizeof(int *));
  ps.nidx = (int **)    R_alloc(nt, sizeof(int *));
  ps.nval = (double **) R_alloc(nt, sizeof(double *));
  ps.smap = (int **)    R_alloc(nt, sizeof(int *));
  ps.sfeat= (int **)    R_alloc(nt, sizeof(int *));
  ps.sval2= (double **) R_alloc(nt, sizeof(double *));

  int A = T > n ? T : n;
  for (int t = 0; t < nt; t++) {
    if (tree) ptree_ws_init(ps.ws + t, cap);
    ps.tips[t] = (int *)    R_alloc(cap, sizeof(int));
    ps.tw[t]   = (double *) R_alloc(cap, sizeof(double));
    ps.arr[t] = ps.rec[t] = ps.npos[t] = ps.nidx[t] = ps.smap[t] = ps.sfeat[t] = NULL;
    ps.nval[t] = ps.sval2[t] = NULL;
    if (n_draws > 0 && model != PS_SUPPLIED) {
      ps.arr[t] = (int *) R_alloc(A, sizeof(int));
      ps.rec[t] = (int *) R_alloc(A, sizeof(int));
      for (int k = 0; k < A; k++) ps.arr[t][k] = k;
    }
    if (n_draws > 0 && model == PS_TAXA_LABELS)
      ps.npos[t] = (int *) R_alloc(D > 0 ? D : 1, sizeof(int));
    if (n_draws > 0 && model == PS_FREQUENCY) {
      ps.npos[t]  = (int *)    R_alloc((size_t) n + 1, sizeof(int));
      ps.nidx[t]  = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
      ps.nval[t]  = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
      ps.smap[t]  = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
      ps.sfeat[t] = (int *)    R_alloc((size_t) n + 1, sizeof(int));
    }
    if (n_draws > 0 && (model == PS_IND_SWAP || model == PS_TRIAL_SWAP)) {
      ps.sfeat[t] = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
      ps.sval2[t] = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
      if (!ps.nested) {
        size_t cells = (size_t) n * D;
        ps.smap[t] = (int *) R_alloc(cells, sizeof(int));
        for (size_t k = 0; k < cells; k++) ps.smap[t][k] = -1;
      }
    }
  }


  //--------------------------------------------------
  // Observed values.
  //--------------------------------------------------
  const char *names[] = { "n_taxa", "obs", "mean", "sd", "rank", "n_used", "nested", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, names));
  SEXP sexp_ntaxa = PROTECT(allocVector(INTSXP, n));
  SEXP sexp_obs   = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_mean  = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_sd    = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_rank  = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_nused = PROTECT(allocVector(INTSXP, n));
  double *obs = REAL(sexp_obs), *mean = REAL(sexp_mean), *sd = REAL(sexp_sd);
  double *rank = REAL(sexp_rank);
  int    *nused = INTEGER(sexp_nused);
  for (int i = 0; i < n; i++) {
    INTEGER(sexp_ntaxa)[i] = n_taxa[i];
    mean[i] = sd[i] = rank[i] = NA_REAL;
    nused[i] = 0;
  }

  ps.vals = obs;
  run_parallel(obs_worker, threads_for(nt, (double) nnz * 20 + n, 2e5), n);

  //--------------------------------------------------
  // Exact moments.
  //--------------------------------------------------
  if (model == PS_EXACT) {
    long double mu, P2, R2;
    exact_sums(&mu, &P2, &R2);
    for (int i = 0; i < n; i++)
      exact_moments(n_taxa[i], rval + rpos[i], mu, P2, R2, mean + i, sd + i);
  }

  //--------------------------------------------------
  // Draws, in batches; summed per sample in draw order.
  //--------------------------------------------------
  if (n_draws > 0) {

    double *wm = (double *) R_alloc(n, sizeof(double)); // Welford mean
    double *m2 = (double *) R_alloc(n, sizeof(double));
    double *lt = (double *) R_alloc(n, sizeof(double)); // below obs (ties 1/2)
    for (int i = 0; i < n; i++) { wm[i] = m2[i] = lt[i] = 0; }

    size_t per = (size_t) (n > 0 ? n : 1);
    int    B   = (int) ((size_t) 4194304 / per);
    if (B < 8 * nt) B = 8 * nt;
    if (B < 1)      B = 1;
    if (B > n_draws) B = n_draws;
    ps.vals = (double *) R_alloc((size_t) B * per, sizeof(double));

    for (int r0 = 0; r0 < n_draws; r0 += B) {

      ps.r0 = r0;
      ps.nb = (r0 + B <= n_draws) ? B : n_draws - r0;
      run_parallel(null_worker, nt, nt > 1 ? 1000 : 1);

      for (int b = 0; b < ps.nb; b++) {
        const double *v = ps.vals + (size_t) b * n;
        for (int i = 0; i < n; i++) {
          double x = v[i];
          if (ISNAN(x)) continue;
          int    k = ++nused[i];
          double d = x - wm[i];
          wm[i] += d / k;
          m2[i] += d * (x - wm[i]);
          double o = obs[i];
          if (!ISNAN(o)) {
            double tol = 1e-12 * (fabs(o) > fabs(x) ? fabs(o) : fabs(x));
            if      (x < o - tol) lt[i] += 1;
            else if (x <= o + tol) lt[i] += 0.5;
          }
        }
      }

      if (user_interrupted()) {
        UNPROTECT(7);
        error("Interrupted.");
      }
    }

    for (int i = 0; i < n; i++) {
      mean[i] = nused[i] > 0 ? wm[i] : NA_REAL;
      sd[i]   = nused[i] > 1 ? sqrt(m2[i] / (nused[i] - 1)) : NA_REAL;
      rank[i] = ISNAN(obs[i]) ? NA_REAL : 1 + lt[i];
    }
  }

  SET_VECTOR_ELT(result, 0, sexp_ntaxa);
  SET_VECTOR_ELT(result, 1, sexp_obs);
  SET_VECTOR_ELT(result, 2, sexp_mean);
  SET_VECTOR_ELT(result, 3, sexp_sd);
  SET_VECTOR_ELT(result, 4, sexp_rank);
  SET_VECTOR_ELT(result, 5, sexp_nused);
  SET_VECTOR_ELT(result, 6, ScalarLogical(ps.nested));
  memset(&ps, 0, sizeof(ps));
  UNPROTECT(7);
  return result;
}
