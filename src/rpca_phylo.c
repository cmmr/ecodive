// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// C port of the tree handling in gemelli 0.0.13's phylogenetic_rpca()
// (gemelli/rpca.py, and bp_read_phylogeny(), tree_topology_filter() and
// fast_unifrac() in gemelli/preprocessing.py), which carries this notice:
//
//   BSD 3-Clause License
//
//   Copyright (c) 2019, Cameron Martino
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
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The tree operations that gemelli calls follow scikit-bio 0.7.4
// (TreeNode.prune(), TreeNode.assign_ids(), TreeNode.to_array() and
// vectorize_counts_and_tree(), in skbio/tree/_tree.py and
// skbio/diversity/_util.py and _phylogenetic.pyx), which carries this
// notice:
//
//   Copyright (c) 2013--, scikit-bio development team.
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
//   * Neither the names scikit-bio, skbio, or biocore nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// and iow 1.0.8 (the `bp` module: BP.shear() and BP.collapse(), in
// bp/_bp.pyx), which carries this notice:
//
//   Copyright (c) 2016, Daniel McDonald
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
//   * Neither the name of improved-octo-waddle nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// This is a modified port: the operations are rewritten on index arrays
// rather than on tree objects, and the counts are expanded sparsely and in
// parallel. See plans/rpca_phylo.md and rpca_phylo()'s "Differences from
// gemelli" for how the results differ.

/*
 * Phylogenetic robust Aitchison PCA: Martino et al. (2022), mSystems 7(3),
 * e00050-22, as implemented by gemelli's phylogenetic_rpca().
 *
 *   1. Filter the tip table as rpca() does (rpca_filter() in rpca.c), but
 *      keeping samples with a single non-zero tip.
 *   2. Shear the tree to the root, the present tips (non-zero after
 *      filtering) and their ancestors, then collapse single-child nodes
 *      into their child, adding their length to it (iow's BP.shear() and
 *      BP.collapse()). BP.collapse() keeps the root, but adds its length
 *      to a single child all the same.
 *   3. iow's to_skbio_treenode() then drops the root's length, so the
 *      root's is 0. scikit-bio's prune() merges a root that has a single
 *      child with that child, taking its length (which includes the
 *      root's own, from step 2).
 *   4. min_depth (gemelli's tree_topology_filter()): in postorder, each
 *      internal node with at most min_depth tips below it is removed, and
 *      its tips are appended, in order, to its parent's children.
 *   5. Number the nodes as scikit-bio's assign_ids(): in postorder, each
 *      node's children get consecutive ids, and the root the last. These
 *      are the columns of the node table, and give internal nodes their
 *      names (n<id>). A child's id is always below its parent's.
 *   6. Each sample's count of a node is the sum of its tips' counts. An
 *      entry is observed when that sum and the node's branch length are
 *      both positive; its value is log(count) + log(length), centered by
 *      the sample's mean over its observed entries (gemelli's
 *      matrix_rclr() with branch lengths). A sample whose observed values
 *      are all equal (so all its centered values are 0, which gemelli
 *      treats as missing) has no observed entries.
 *   7. rpca_core() (rpca.c) on the n x K table, with K counting every
 *      node, observed or not.
 *
 * Threads: step 6 runs one sample at a time, each on one thread, which
 * writes only that sample's entries, so results don't depend on the number
 * of threads. Each node's count is summed over its children in order of
 * their ids. The tree steps are linear in the tree's size and run on the
 * main thread.
 */

#include "ecodive.h"

#include <limits.h> // INT_MAX


//======================================================
// Shared with the worker threads (step 6).
//======================================================
static int           x_n;          // Kept samples.
static int           x_K;          // Nodes (columns).
static const int    *x_tpos;       // The tip table (new_ecomatrix()).
static const int    *x_tidx;
static const double *x_tval;
static const int    *x_sidx;       // Kept samples' rows in the tip table.
static const int    *x_tipcol;     // Each tip's column, or -1.
static const int    *x_colpar;     // Each column's parent column, or -1.
static const double *x_blen;       // Each column's branch length.
static const double *x_logb;       // log(branch length), where positive.
static int           x_fill;       // 0: count entries; 1: fill them in.
static int          *x_n_touch;    // Per sample: nodes with a count.
static int          *x_n_obs;      // Per sample: observed entries.
static const int    *x_opos;       // Observed entries' CSR.
static int          *x_oidx;
static double       *x_oval;
static const int    *x_npos;       // Raw node counts' CSR (or NULL).
static int          *x_nidx;
static double       *x_nval;
static double       *x_acc;        // Per thread: K sums,
static int          *x_stamp;      //   K marks,
static int          *x_list;       //   K columns,
static double       *x_v;          //   K values.


static int cmp_int(const void *a, const void *b) {
  int x = *(const int *) a, y = *(const int *) b;
  return (x > y) - (x < y);
}


//------------------------------------------------------
// One sample's node counts. On return, list[0..len-1]
// holds the columns with a count (increasing), acc[]
// their counts, and v[0..n_obs-1] the observed values
// before centering (n_obs = 0 if they are all equal).
//------------------------------------------------------
static int sample_nodes(int k, double *acc, int *stamp, int *list, double *v, int *n_obs) {

  int i   = x_sidx[k];
  int tag = k + 1;
  int len = 0;

  for (int e = x_tpos[i]; e < x_tpos[i + 1]; e++) {
    int c = x_tipcol[x_tidx[e]];
    if (c < 0) continue;
    stamp[c] = tag; acc[c] = x_tval[e]; list[len++] = c;
    for (int p = x_colpar[c]; p >= 0 && stamp[p] != tag; p = x_colpar[p]) {
      stamp[p] = tag; acc[p] = 0; list[len++] = p;
    }
  }

  // Children before parents, in order of their ids.
  qsort(list, (size_t) len, sizeof(int), cmp_int);
  for (int a = 0; a < len; a++) {
    int p = x_colpar[list[a]];
    if (p >= 0) acc[p] += acc[list[a]];
  }

  int nv = 0, flat = 1;
  for (int a = 0; a < len; a++) {
    int c = list[a];
    if (!(acc[c] > 0 && x_blen[c] > 0)) continue;
    v[nv] = log(acc[c]) + x_logb[c];
    if (v[nv] != v[0]) flat = 0;
    nv++;
  }

  *n_obs = flat ? 0 : nv;
  return len;
}


static void *w_expand(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  size_t  K     = (size_t) x_K;
  double *acc   = x_acc   + K * thread_i;
  int    *stamp = x_stamp + K * thread_i;
  int    *list  = x_list  + K * thread_i;
  double *v     = x_v     + K * thread_i;

  for (int k = thread_i; k < x_n; k += n_threads) {

    int n_obs;
    int len = sample_nodes(k, acc, stamp, list, v, &n_obs);

    if (!x_fill) {
      x_n_touch[k] = len;
      x_n_obs[k]   = n_obs;
      continue;
    }

    if (n_obs > 0) {
      double mean = 0;
      for (int a = 0; a < n_obs; a++) mean += v[a];
      mean /= n_obs;
      int at = x_opos[k], a = 0;
      for (int b = 0; b < len; b++) {
        int c = list[b];
        if (!(acc[c] > 0 && x_blen[c] > 0)) continue;
        x_oidx[at] = c;
        x_oval[at] = v[a++] - mean;
        at++;
      }
    }

    if (x_npos) {
      int at = x_npos[k];
      for (int b = 0; b < len; b++) {
        x_nidx[at] = list[b];
        x_nval[at] = acc[list[b]];
        at++;
      }
    }
  }

  return NULL;
}


// Frees new_ecomatrix()'s allocations, then raises an R error.
#define STOP(...) do { free_all(); error(__VA_ARGS__); } while (0)


//======================================================
// Doubly linked lists of children, on node indices.
//======================================================
static int *ch_first, *ch_last, *ch_next, *ch_prev, *ch_par;

static void ch_append(int p, int c) {
  ch_par[c]  = p;
  ch_next[c] = -1;
  ch_prev[c] = ch_last[p];
  if (ch_last[p] >= 0) ch_next[ch_last[p]] = c; else ch_first[p] = c;
  ch_last[p] = c;
}

static void ch_unlink(int c) {
  int p = ch_par[c];
  if (ch_prev[c] >= 0) ch_next[ch_prev[c]] = ch_next[c]; else ch_first[p] = ch_next[c];
  if (ch_next[c] >= 0) ch_prev[ch_next[c]] = ch_prev[c]; else ch_last[p]  = ch_prev[c];
  ch_par[c] = ch_next[c] = ch_prev[c] = -1;
}

// Postorder of the subtree at `root` (children in list
// order) into `out`; returns its length. `stack` and
// `pos` are scratch of the tree's size.
static int postorder(int root, int *out, int *stack, int *pos) {
  int len = 0, top = 0;
  stack[top] = root; pos[top] = ch_first[root];
  while (top >= 0) {
    int c = pos[top];
    if (c >= 0) {
      pos[top] = ch_next[c];
      stack[++top] = c; pos[top] = ch_first[c];
    }
    else {
      out[len++] = stack[top--];
    }
  }
  return len;
}


//======================================================
// R interface.
//
// `sexp_counts`/`sexp_margin`: the tip table, its
// features in the order of the tree's tips (as
// validate_tree() leaves it). `sexp_edge`,
// `sexp_edge_len`, `sexp_nnode`, `sexp_root_len`: the
// phylo object's edge matrix, edge lengths (NA counts
// as 0, as gemelli's missing lengths), number of
// internal nodes, and root length (0 if it has none;
// it only counts when the root has a single child).
// `sexp_min_depth`: gemelli's min_depth.
// `sexp_nc`: also return the node counts. The rest are
// as C_rpca()'s.
//
// Returns rpca_core()'s list, followed by:
//   samples:   the kept samples (1-based);
//   node_tip:  for each column, its tip (1-based) or 0;
//   node_par:  for each column, its parent column
//              (1-based), or 0 for the root;
//   node_len:  each column's branch length;
//   phy_edge:  the final tree's edge matrix in ape's
//              numbering (tips in preorder, then
//              internal nodes in preorder, root first),
//              rows in preorder;
//   phy_col:   for each ape node number, its column;
//   nc_pos, nc_idx, nc_val: the node counts as CSR
//              (samples x columns; with `sexp_nc`);
//   rclr_pos, rclr_idx, rclr: the observed values
//              as CSR (with `sexp_debug`).
//======================================================
SEXP C_rpca_phylo(
    SEXP sexp_counts,    SEXP sexp_margin,   SEXP sexp_edge,
    SEXP sexp_edge_len,  SEXP sexp_nnode,    SEXP sexp_root_len,
    SEXP sexp_k,         SEXP sexp_iters,    SEXP sexp_tol,
    SEXP sexp_msc,       SEXP sexp_mfc,      SEXP sexp_mff,
    SEXP sexp_min_depth, SEXP sexp_dm,       SEXP sexp_nc,
    SEXP sexp_debug,     SEXP sexp_step,     SEXP sexp_n_threads ) {

  int    k_req     = asInteger(sexp_k);
  int    max_iter  = asInteger(sexp_iters);
  double tol       = asReal(sexp_tol);
  double msc       = asReal(sexp_msc);
  double mfc       = asReal(sexp_mfc);
  double mff       = asReal(sexp_mff);
  int    min_depth = asInteger(sexp_min_depth);
  int    want_dm   = asLogical(sexp_dm);
  int    want_nc   = asLogical(sexp_nc);
  int    debug     = asLogical(sexp_debug);
  int    n_threads = asInteger(sexp_n_threads);
  double step      = asReal(sexp_step);
  double root_len  = asReal(sexp_root_len);
  int    n_int     = asInteger(sexp_nnode);

  if (k_req < 1 || max_iter < 1 || n_threads < 1 || !(step > 0) || min_depth < 0)
    error("Invalid arguments to C_rpca_phylo()."); // # nocov

  if (!isInteger(sexp_edge) || !isMatrix(sexp_edge) || ncols(sexp_edge) != 2 ||
      !isReal(sexp_edge_len) || n_int < 1 || n_int == NA_INTEGER)
    error("`tree` is not a valid phylo object.");

  int        n_edges = nrows(sexp_edge);
  const int *edge    = INTEGER(sexp_edge);
  const double *elen = REAL(sexp_edge_len);

  if (XLENGTH(sexp_edge_len) != n_edges)
    error("`tree` is not a valid phylo object: it has %d edges but %d edge lengths.",
          n_edges, (int) XLENGTH(sexp_edge_len));

  for (int e = 0; e < n_edges; e++)
    if (!ISNA(elen[e]) && !(R_FINITE(elen[e]) && elen[e] >= 0))
      error("`tree`'s branch lengths must be finite and non-negative.");
  if (ISNA(root_len)) root_len = 0;
  if (!(R_FINITE(root_len) && root_len >= 0))
    error("`tree`'s root length must be finite and non-negative.");


  //--------------------------------------------------
  // The tree: N nodes, tips 0 .. n_tips-1, then the
  // internal nodes, as in ape (1-based there). Each
  // node's children in the order of the edge matrix.
  //--------------------------------------------------
  int N      = n_edges + 1;
  int n_tips = N - n_int;
  if (n_tips < 1)
    error("`tree` is not a valid phylo object."); // # nocov

  size_t Nz = (size_t) N;
  int    *par  = (int *)    R_alloc(Nz, sizeof(int));
  double *blen = (double *) R_alloc(Nz, sizeof(double));
  ch_first = (int *) R_alloc(Nz, sizeof(int));
  ch_last  = (int *) R_alloc(Nz, sizeof(int));
  ch_next  = (int *) R_alloc(Nz, sizeof(int));
  ch_prev  = (int *) R_alloc(Nz, sizeof(int));
  ch_par   = (int *) R_alloc(Nz, sizeof(int));

  for (int v = 0; v < N; v++) {
    par[v] = -1; blen[v] = 0;
    ch_first[v] = ch_last[v] = ch_next[v] = ch_prev[v] = ch_par[v] = -1;
  }

  for (int e = 0; e < n_edges; e++) {
    int p = edge[e] - 1, c = edge[n_edges + e] - 1;
    if (p < n_tips || p >= N || c < 0 || c >= N || par[c] >= 0 || c == p)
      error("`tree` is not a valid phylo object: malformed edge matrix.");
    par[c]  = p;
    blen[c] = ISNA(elen[e]) ? 0 : elen[e];
    ch_append(p, c);
  }

  int root = -1;
  for (int v = n_tips; v < N; v++)
    if (par[v] < 0) {
      if (root >= 0) error("`tree` is not a valid phylo object: it has more than one root.");
      root = v;
    }
  if (root < 0)
    error("`tree` is not a valid phylo object: it has no root."); // # nocov
  blen[root] = root_len;

  // Preorder (children in order), which also checks
  // that every node hangs from the root.
  int *pre   = (int *) R_alloc(Nz, sizeof(int));
  int *stack = (int *) R_alloc(Nz, sizeof(int));
  int *spos  = (int *) R_alloc(Nz, sizeof(int));
  int  n_pre = 0;
  {
    int top = 0;
    stack[0] = root;
    while (top >= 0) {
      int v = stack[top--];
      pre[n_pre++] = v;
      if (n_pre > N) break; // # nocov
      // Push children in reverse, so the first is visited first.
      for (int c = ch_last[v]; c >= 0; c = ch_prev[c]) {
        if (top + 1 >= N) { n_pre = N + 1; break; } // # nocov
        stack[++top] = c;
      }
    }
  }
  if (n_pre != N)
    error("`tree` is not a valid phylo object: not every node descends from the root.");
  for (int v = 0; v < n_tips; v++)
    if (ch_first[v] >= 0)
      error("`tree` is not a valid phylo object: a tip has children."); // # nocov
  for (int v = n_tips; v < N; v++)
    if (ch_first[v] < 0)
      error("`tree` is not a valid phylo object: an internal node has no children.");


  //--------------------------------------------------
  // Read and filter the tip table.
  //--------------------------------------------------
  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int           N0   = em->n_samples, P0 = em->n_otus, e1 = em->nnz;
  const int    *tpos = em->pos_vec, *tidx = em->otu_vec;
  const double *tval = em->val_vec;

  if (P0 != n_tips)
    STOP("`counts` has %d features but `tree` has %d tips.", P0, n_tips); // # nocov

  for (int e = 0; e < e1; e++) {
    if (!R_FINITE(tval[e])) STOP("`counts` must not contain NA, NaN, or infinite values.");
    if (tval[e] < 0)        STOP("`counts` must not contain negative values.");
  }

  int *f_new = (int *) R_alloc(P0 > 0 ? P0 : 1, sizeof(int));
  int *s_new = (int *) R_alloc(N0 > 0 ? N0 : 1, sizeof(int));
  int  n, m_tips, n_single;

  rpca_filter(N0, P0, tpos, tidx, tval, msc, mfc, mff, 0,
              s_new, f_new, NULL, &n_single, &n, &m_tips);

  if (n < 2 || m_tips < 1)
    STOP(
      "After filtering, %d sample%s and %d tip%s remain; at least 2 samples "
      "are needed, with counts on at least 2 tree nodes.",
      n, (n == 1) ? "" : "s", m_tips, (m_tips == 1) ? "" : "s" );


  //--------------------------------------------------
  // 2. Shear, then collapse single-child nodes (in
  // preorder, so lengths add down chains). The root is
  // kept, but still adds its length to a single child.
  //--------------------------------------------------
  int    *keep = (int *)    R_alloc(Nz, sizeof(int));
  int    *gone = (int *)    R_alloc(Nz, sizeof(int));
  double *len  = (double *) R_alloc(Nz, sizeof(double));

  for (int v = 0; v < N; v++) { keep[v] = 0; gone[v] = 0; len[v] = blen[v]; }
  keep[root] = 1;
  for (int t = 0; t < n_tips; t++)
    if (f_new[t] >= 0)
      for (int v = t; v >= 0 && !keep[v]; v = par[v]) keep[v] = 1;

  int *only = (int *) R_alloc(Nz, sizeof(int)); // A node's single kept child.
  for (int i = 0; i < N; i++) {
    int v = pre[i];
    if (!keep[v] || v < n_tips) continue;
    int kc = 0;
    for (int c = ch_first[v]; c >= 0; c = ch_next[c])
      if (keep[c]) { kc++; only[v] = c; }
    if (kc == 1) {
      len[only[v]] += len[v];
      if (v != root) gone[v] = 1;
    }
    else only[v] = -1;
  }


  //--------------------------------------------------
  // 3. iow's to_skbio_treenode() drops the root's
  // length (None, so 0 in the node table). A root with
  // a single child then takes the child's place and
  // length, which includes the root's from step 2
  // (scikit-bio's prune()). Then the remaining tree's children
  // lists: each kept node under its nearest ancestor
  // that is not gone, in preorder, so in order.
  //--------------------------------------------------
  if (only[root] < 0) len[root] = 0;
  else {
    int c = only[root];
    while (gone[c]) c = only[c];
    len[root] = len[c];
    gone[c] = 1;
    if (c < n_tips) // A single tip: nothing to ordinate.
      STOP("After filtering, only 1 tip has counts; at least 2 are needed.");
  }

  int *up = (int *) R_alloc(Nz, sizeof(int)); // Nearest node at or above that isn't gone.
  for (int v = 0; v < N; v++) ch_first[v] = ch_last[v] = ch_next[v] = ch_prev[v] = ch_par[v] = -1;
  for (int i = 0; i < N; i++) {
    int v = pre[i];
    if (!keep[v]) continue;
    up[v] = (v == root || !gone[v]) ? v : up[par[v]];
    if (v != root && !gone[v]) ch_append(up[par[v]], v);
  }


  //--------------------------------------------------
  // 4. min_depth: in postorder, remove each internal
  // node with <= min_depth tips, appending its children
  // (by then, all tips) to its parent's.
  //--------------------------------------------------
  int *post   = (int *) R_alloc(Nz, sizeof(int));
  int  n_post = postorder(root, post, stack, spos);

  if (min_depth > 0) {

    int *ntip = (int *) R_alloc(Nz, sizeof(int));
    for (int i = 0; i < n_post; i++) {
      int v = post[i];
      ntip[v] = (v < n_tips) ? 1 : 0;
      for (int c = ch_first[v]; c >= 0; c = ch_next[c]) ntip[v] += ntip[c];
    }

    if (ntip[root] <= min_depth)
      STOP(
        "`min_depth` (%d) must be less than the number of tips with counts (%d).",
        min_depth, ntip[root] );

    for (int i = 0; i < n_post; i++) {
      int x = post[i];
      if (x < n_tips || ntip[x] > min_depth) continue;
      int p = ch_par[x];
      for (int c = ch_first[x]; c >= 0; ) {
        int next = ch_next[c];
        ch_unlink(c);
        ch_append(p, c);
        c = next;
      }
      ch_unlink(x);
    }

    n_post = postorder(root, post, stack, spos);
  }


  //--------------------------------------------------
  // 5. Column ids, as scikit-bio's assign_ids().
  //--------------------------------------------------
  int  K   = n_post;
  int *id  = (int *) R_alloc(Nz, sizeof(int));
  int  cur = 0;
  for (int i = 0; i < n_post; i++)
    for (int c = ch_first[post[i]]; c >= 0; c = ch_next[c])
      id[c] = cur++;
  id[root] = cur++;
  if (cur != K)
    STOP("Internal error numbering tree nodes."); // # nocov

  if (K < 2)
    STOP("After filtering, the tree has 1 node; at least 2 are needed."); // # nocov

  size_t  Kz     = (size_t) K;
  int    *colnod = (int *)    R_alloc(Kz, sizeof(int));
  int    *colpar = (int *)    R_alloc(Kz, sizeof(int));
  double *collen = (double *) R_alloc(Kz, sizeof(double));
  double *logb   = (double *) R_alloc(Kz, sizeof(double));
  int    *tipcol = (int *)    R_alloc(n_tips, sizeof(int));

  for (int t = 0; t < n_tips; t++) tipcol[t] = -1;
  for (int i = 0; i < n_post; i++) {
    int v = post[i], c = id[v];
    colnod[c] = v;
    colpar[c] = (v == root) ? -1 : id[ch_par[v]];
    collen[c] = len[v];
    logb[c]   = (len[v] > 0) ? log(len[v]) : 0;
    if (v < n_tips) tipcol[v] = c;
  }


  //--------------------------------------------------
  // 6. Node counts and their robust clr: two passes,
  // counting then filling.
  //--------------------------------------------------
  int *s_idx = (int *) R_alloc(n, sizeof(int));
  for (int i = 0; i < N0; i++) if (s_new[i] >= 0) s_idx[s_new[i]] = i;

  // Rough work: each non-zero tip count, times log K.
  int nt = threads_for(n_threads, (double) e1 * (1 + log((double) K)) * 8, 65536);
  if (nt > n) nt = n;

  x_n = n; x_K = K;
  x_tpos = tpos; x_tidx = tidx; x_tval = tval; x_sidx = s_idx;
  x_tipcol = tipcol; x_colpar = colpar; x_blen = collen; x_logb = logb;
  x_n_touch = (int *) R_alloc(n, sizeof(int));
  x_n_obs   = (int *) R_alloc(n, sizeof(int));
  x_acc     = (double *) R_alloc(Kz * nt, sizeof(double));
  x_stamp   = (int *)    R_alloc(Kz * nt, sizeof(int));
  x_list    = (int *)    R_alloc(Kz * nt, sizeof(int));
  x_v       = (double *) R_alloc(Kz * nt, sizeof(double));
  for (size_t i = 0; i < Kz * nt; i++) x_stamp[i] = 0;

  x_fill = 0;
  run_parallel(w_expand, nt, e1);

  int    *opos = (int *) R_alloc((size_t) n + 1, sizeof(int));
  int    *npos = want_nc ? (int *) R_alloc((size_t) n + 1, sizeof(int)) : NULL;
  double  o_tot = 0, n_tot = 0;
  opos[0] = 0;
  if (npos) npos[0] = 0;
  for (int k = 0; k < n; k++) {
    o_tot += x_n_obs[k];
    n_tot += x_n_touch[k];
  }
  if (o_tot > INT_MAX || (want_nc && n_tot > INT_MAX))
    STOP("The node table has too many non-zero entries (over 2^31 - 1)."); // # nocov
  for (int k = 0; k < n; k++) {
    opos[k + 1] = opos[k] + x_n_obs[k];
    if (npos) npos[k + 1] = npos[k] + x_n_touch[k];
  }

  int nnz = opos[n];
  x_opos = opos;
  x_oidx = (int *)    R_alloc(nnz > 0 ? nnz : 1, sizeof(int));
  x_oval = (double *) R_alloc(nnz > 0 ? nnz : 1, sizeof(double));
  x_npos = npos;
  x_nidx = NULL; x_nval = NULL;
  if (npos) {
    x_nidx = (int *)    R_alloc(npos[n] > 0 ? npos[n] : 1, sizeof(int));
    x_nval = (double *) R_alloc(npos[n] > 0 ? npos[n] : 1, sizeof(double));
  }

  for (size_t i = 0; i < Kz * nt; i++) x_stamp[i] = 0;
  x_fill = 1;
  run_parallel(w_expand, nt, e1);

  free_all();


  //--------------------------------------------------
  // The final tree in ape's numbering.
  //--------------------------------------------------
  int nprot = 0;

  int n_tip2 = 0;
  for (int c = 0; c < K; c++) if (colnod[c] < n_tips) n_tip2++;

  SEXP sexp_phy_edge = PROTECT(allocMatrix(INTSXP, K - 1, 2)); nprot++;
  SEXP sexp_phy_col  = PROTECT(allocVector(INTSXP, K));        nprot++;
  {
    int *ape = (int *) R_alloc(Kz, sizeof(int)); // By column.
    int  nt_ = 0, ni = 0, row = 0, top = 0;
    int *pe  = INTEGER(sexp_phy_edge);
    stack[0] = root;
    while (top >= 0) {
      int v = stack[top--], c = id[v];
      ape[c] = (v < n_tips) ? ++nt_ : n_tip2 + ++ni;
      INTEGER(sexp_phy_col)[ape[c] - 1] = c + 1;
      if (v != root) {
        pe[row]         = ape[id[ch_par[v]]];
        pe[K - 1 + row] = ape[c];
        row++;
      }
      for (int ch = ch_last[v]; ch >= 0; ch = ch_prev[ch]) stack[++top] = ch;
    }
  }


  //--------------------------------------------------
  // 7. OptSpace and PCA.
  //--------------------------------------------------
  int r = k_req;
  if (r > n - 1) r = n - 1;
  if (r > K - 1) r = K - 1;
  if (r > 10)    r = 10;

  SEXP sexp_rclr_pos = R_NilValue, sexp_rclr_idx = R_NilValue, sexp_rclr = R_NilValue;
  if (debug) {
    sexp_rclr_pos = PROTECT(allocVector(INTSXP,  n + 1)); nprot++;
    sexp_rclr_idx = PROTECT(allocVector(INTSXP,  nnz));   nprot++;
    sexp_rclr     = PROTECT(allocVector(REALSXP, nnz));   nprot++;
    memcpy(INTEGER(sexp_rclr_pos), opos, sizeof(int) * ((size_t) n + 1));
    for (int e = 0; e < nnz; e++) {
      INTEGER(sexp_rclr_idx)[e] = x_oidx[e];
      REAL(sexp_rclr)[e]        = x_oval[e];
    }
  }

  SEXP core = PROTECT(rpca_core(
    n, K, opos, x_oidx, x_oval, r, max_iter, tol, want_dm, debug, step, n_threads,
    "There is no variation to ordinate: within each sample, every node's "
    "count times branch length is the same." ));
  nprot++;


  //--------------------------------------------------
  // Results.
  //--------------------------------------------------
  SEXP sexp_sam = PROTECT(allocVector(INTSXP,  n)); nprot++;
  SEXP sexp_tip = PROTECT(allocVector(INTSXP,  K)); nprot++;
  SEXP sexp_par = PROTECT(allocVector(INTSXP,  K)); nprot++;
  SEXP sexp_len = PROTECT(allocVector(REALSXP, K)); nprot++;
  for (int k = 0; k < n; k++) INTEGER(sexp_sam)[k] = s_idx[k] + 1;
  for (int c = 0; c < K; c++) {
    INTEGER(sexp_tip)[c] = (colnod[c] < n_tips) ? colnod[c] + 1 : 0;
    INTEGER(sexp_par)[c] = colpar[c] + 1;
    REAL(sexp_len)[c]    = collen[c];
  }

  SEXP sexp_nc_pos = R_NilValue, sexp_nc_idx = R_NilValue, sexp_nc_val = R_NilValue;
  if (want_nc) {
    int nn = npos[n];
    sexp_nc_pos = PROTECT(allocVector(INTSXP,  n + 1)); nprot++;
    sexp_nc_idx = PROTECT(allocVector(INTSXP,  nn));    nprot++;
    sexp_nc_val = PROTECT(allocVector(REALSXP, nn));    nprot++;
    memcpy(INTEGER(sexp_nc_pos), npos, sizeof(int) * ((size_t) n + 1));
    if (nn > 0) {
      memcpy(INTEGER(sexp_nc_idx), x_nidx, sizeof(int)    * (size_t) nn);
      memcpy(REAL(sexp_nc_val),    x_nval, sizeof(double) * (size_t) nn);
    }
  }

  const char *names[] = {
    "samples", "node_tip", "node_par", "node_len", "phy_edge", "phy_col",
    "nc_pos", "nc_idx", "nc_val", "rclr_pos", "rclr_idx", "rclr" };
  SEXP vals[] = {
    sexp_sam, sexp_tip, sexp_par, sexp_len, sexp_phy_edge, sexp_phy_col,
    sexp_nc_pos, sexp_nc_idx, sexp_nc_val, sexp_rclr_pos, sexp_rclr_idx, sexp_rclr };
  SEXP res = rpca_list_append(core, names, vals, 12);

  UNPROTECT(nprot);
  return res;
}
