// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * The phylogenetic isometric log-ratio transform (PhILR) of Silverman,
 * Washburne, Mukherjee & David (2017), eLife 6, e21887: written from the
 * paper (CC-BY) and checked against the outputs of the philr package,
 * whose code was not read (see plans/philr.md).
 *
 * C_philr_tree():   the tree as philr() uses it: checked, pruned to the
 *                   features of the count table, single-child nodes
 *                   collapsed, polytomies resolved, internal nodes numbered
 *                   in preorder; with the branch-length weights.
 * C_philr_stats():  per-feature sums over samples, for the part weights.
 * C_philr_coords(): the coordinates, O(D) per sample, one sample per task.
 *
 * For a sample with logs L_j = log(x_j + pseudocount) and part weights p_j,
 * let w_j = p_j (L_j - log p_j). For an internal node with "+" clade A (its
 * first child) and "-" clade B (its second), N_A = sum of p over A and
 * S_A = sum of w over A (and likewise for B); then
 *
 *   coordinate = f sqrt(N_A N_B / (N_A + N_B)) (S_A / N_A - S_B / N_B),
 *
 * with f the node's ILR weight. The sums are accumulated child before
 * parent (reverse preorder), so no basis matrix is formed. Samples share
 * nothing, so results are identical for any number of threads.
 */

#include "ecodive.h"


// Frees new_ecomatrix()'s allocations, then raises an R error.
#define STOP(...) do { free_all(); error(__VA_ARGS__); } while (0)


//======================================================
// R interface: the tree that philr() uses.
//
// sexp_edge:     the phylo edge matrix (integer, E x 2).
// sexp_len:      its edge lengths (double, E), or NULL.
// sexp_ntip:     the number of tips (length(tip.label)).
// sexp_tip_feat: for each tip, its feature's 0-based
//                index in the count table, or -1 if
//                the table doesn't have it.
//
// Every feature has one tip (R checks), so the D kept
// tips give a tree with D - 1 internal nodes after:
//
//   1. dropping the tips that aren't features, and the
//      internal nodes left without kept tips;
//   2. collapsing each node left with one child into
//      that child, adding their branch lengths (a root
//      with one child is dropped, with its length);
//   3. resolving each node with k > 2 children c_1 ..
//      c_k into a ladder, as ape::multi2di(random =
//      FALSE) does: the node gets c_1 and a new node
//      (branch length 0), which gets c_2 and the next,
//      ..., the last getting c_{k-1} and c_k.
//
// Internal nodes are numbered k = 0 .. D - 2 in
// preorder (first child first; the root is 0); node k
// has id D + k, and a tip has its feature's index.
//
// Returns a list:
//   c1, c2:     each node's first ("+") and second
//               ("-") child (ids);
//   orig:       each node's ape number in the input
//               tree, or 0 for a node added in step 3;
//   ladder:     for an added node, the ape number of
//               the node it resolves, else 0;
//   rung:       for an added node, its place in the
//               ladder (1, 2, ...), else 0;
//   edge_par, edge_child, edge_len: the edges, in
//               preorder (the cladewise order);
//   blw:        sum of the two child branch lengths,
//               zero-length tip branches given the
//               smallest positive branch length;
//   md:         blw plus each child's mean distance to
//               its tips (with the same lengths);
//   n_poly:     the number of nodes resolved in step 3;
//   modified:   whether any of steps 1-3 changed the
//               tree.
// blw and md are NA when a branch length is missing,
// or when no branch length is positive.
//======================================================
SEXP C_philr_tree(SEXP sexp_edge, SEXP sexp_len, SEXP sexp_ntip, SEXP sexp_tip_feat) {

  if (!isInteger(sexp_edge) || !isMatrix(sexp_edge) || ncols(sexp_edge) != 2)
    error("`tree` is not a valid phylo object: malformed edge matrix.");

  int        E      = nrows(sexp_edge);
  const int *edge   = INTEGER(sexp_edge);
  int        n_tips = asInteger(sexp_ntip);
  const int *tfeat  = INTEGER(sexp_tip_feat);
  int        N      = E + 1;

  const double *elen = isNull(sexp_len) ? NULL : REAL(sexp_len);
  if (elen && XLENGTH(sexp_len) != E)
    error("`tree` is not a valid phylo object: it has %d edges but %d edge lengths.",
          E, (int) XLENGTH(sexp_len));

  if (n_tips < 2 || n_tips > N || XLENGTH(sexp_tip_feat) != n_tips)
    error("`tree` is not a valid phylo object."); // # nocov

  size_t Nz = (size_t) N;
  int    *par   = (int *)    R_alloc(Nz, sizeof(int));
  double *blen  = (double *) R_alloc(Nz, sizeof(double));
  int    *first = (int *)    R_alloc(Nz, sizeof(int));
  int    *last  = (int *)    R_alloc(Nz, sizeof(int));
  int    *next  = (int *)    R_alloc(Nz, sizeof(int));

  for (int v = 0; v < N; v++) {
    par[v] = -1; blen[v] = 0; first[v] = last[v] = next[v] = -1;
  }

  // Nodes are 0-based ape numbers: tips 0 .. n_tips - 1.
  for (int e = 0; e < E; e++) {
    int p = edge[e] - 1, c = edge[E + e] - 1;
    if (edge[e] == NA_INTEGER || edge[E + e] == NA_INTEGER ||
        p < n_tips || p >= N || c < 0 || c >= N || par[c] >= 0 || c == p)
      error("`tree` is not a valid phylo object: malformed edge matrix.");
    par[c]  = p;
    blen[c] = elen ? elen[e] : NA_REAL;
    if (last[p] >= 0) next[last[p]] = c; else first[p] = c;
    last[p] = c;
  }

  int root = -1;
  for (int v = n_tips; v < N; v++) {
    if (first[v] < 0)
      error("`tree` is not a valid phylo object: internal node %d has no children.", v + 1);
    if (par[v] < 0) {
      if (root >= 0) error("`tree` is not a valid phylo object: it has more than one root.");
      root = v;
    }
  }
  if (root < 0)
    error("`tree` is not a valid phylo object: it has no root."); // # nocov
  blen[root] = 0;

  for (int e = 0; e < E; e++) {
    double l = blen[edge[E + e] - 1];
    if (!ISNA(l) && !(R_FINITE(l) && l >= 0))
      error("`tree`'s branch lengths must be finite and non-negative.");
  }


  // Preorder, which also checks that every node hangs
  // from the root.
  int *pre   = (int *) R_alloc(Nz, sizeof(int));
  int *stack = (int *) R_alloc(Nz, sizeof(int));
  int  n_pre = 0;
  {
    int top = 0;
    stack[0] = root;
    while (top >= 0) {
      int v = stack[top--];
      pre[n_pre++] = v;
      // Push children in reverse: count them, then fill.
      int nc = 0;
      for (int c = first[v]; c >= 0; c = next[c]) nc++;
      if (top + nc >= N) break; // # nocov
      int at = top + nc;
      for (int c = first[v]; c >= 0; c = next[c]) stack[at--] = c;
      top += nc;
    }
  }
  if (n_pre != N)
    error("`tree` is not a valid phylo object: not every node descends from the root.");


  // Kept tips under each node, and each node's kept
  // children (those with kept tips), in edge order.
  int *n_kept = (int *) R_alloc(Nz, sizeof(int));
  int *n_kch  = (int *) R_alloc(Nz, sizeof(int));
  int  D      = 0;
  for (int v = 0; v < N; v++) n_kept[v] = n_kch[v] = 0;
  for (int t = 0; t < n_tips; t++)
    if (tfeat[t] >= 0) { n_kept[t] = 1; D++; }
  for (int a = N - 1; a > 0; a--) {
    int v = pre[a];
    if (n_kept[v] > 0) { n_kept[par[v]] += n_kept[v]; n_kch[par[v]]++; }
  }
  if (D < 2 || n_kept[root] != D)
    error("`tree` must have at least two of the features as tips."); // # nocov

  int *kpos = (int *) R_alloc(Nz + 1, sizeof(int));
  kpos[0] = 0;
  for (int v = 0; v < N; v++) kpos[v + 1] = kpos[v] + n_kch[v];
  int *kch = (int *) R_alloc((size_t) kpos[N] + 1, sizeof(int));
  for (int v = n_tips; v < N; v++) {
    int at = kpos[v];
    for (int c = first[v]; c >= 0; c = next[c])
      if (n_kept[c] > 0) kch[at++] = c;
  }

  int modified = (D < n_tips);
  int n_poly   = 0;


  // The root of the pruned tree: below any chain of
  // single-child nodes.
  while (root >= n_tips && n_kch[root] == 1) {
    root = kch[kpos[root]];
    modified = 1;
  }


  // Emit the final tree in preorder. A stack item is an
  // input node v from its kept child s on (s > 0: an
  // added ladder node), with its parent's final id (-1
  // for the root), which child it is (0 or 1), and its
  // branch length.
  int ni = D - 1;
  SEXP res_c1    = PROTECT(allocVector(INTSXP,  ni));
  SEXP res_c2    = PROTECT(allocVector(INTSXP,  ni));
  SEXP res_orig  = PROTECT(allocVector(INTSXP,  ni));
  SEXP res_lad   = PROTECT(allocVector(INTSXP,  ni));
  SEXP res_rung  = PROTECT(allocVector(INTSXP,  ni));
  SEXP res_epar  = PROTECT(allocVector(INTSXP,  2 * ni));
  SEXP res_echd  = PROTECT(allocVector(INTSXP,  2 * ni));
  SEXP res_elen  = PROTECT(allocVector(REALSXP, 2 * ni));
  SEXP res_blw   = PROTECT(allocVector(REALSXP, ni));
  SEXP res_md    = PROTECT(allocVector(REALSXP, ni));

  int    *c1   = INTEGER(res_c1),   *c2   = INTEGER(res_c2);
  int    *orig = INTEGER(res_orig), *lad  = INTEGER(res_lad), *rung = INTEGER(res_rung);
  int    *epar = INTEGER(res_epar), *echd = INTEGER(res_echd);
  double *elo  = REAL(res_elen);

  size_t  Sz   = Nz + (size_t) D + 2;
  int    *st_v = (int *)    R_alloc(Sz, sizeof(int));
  int    *st_s = (int *)    R_alloc(Sz, sizeof(int));
  int    *st_p = (int *)    R_alloc(Sz, sizeof(int));
  int    *st_w = (int *)    R_alloc(Sz, sizeof(int));
  double *st_l = (double *) R_alloc(Sz, sizeof(double));

  // Final id's branch length (tips 0 .. D-1, then nodes).
  double *flen = (double *) R_alloc((size_t) 2 * D, sizeof(double));

  int top = 0, k = 0, n_e = 0;
  st_v[0] = root; st_s[0] = 0; st_p[0] = -1; st_w[0] = 0; st_l[0] = 0;

  while (top >= 0) {

    int    v = st_v[top], s = st_s[top], p = st_p[top], w = st_w[top];
    double l = st_l[top];
    top--;

    int id;
    if (v < n_tips) {
      id = tfeat[v];
    }
    else {
      if (k >= ni) error("Internal error in C_philr_tree()."); // # nocov
      id      = D + k;
      orig[k] = s ? 0 : v + 1;
      lad[k]  = s ? v + 1 : 0;
      rung[k] = s;
      k++;
    }
    flen[id] = l;

    if (p >= 0) {
      if (w == 0) c1[p - D] = id; else c2[p - D] = id;
      epar[n_e] = p; echd[n_e] = id; elo[n_e] = l; n_e++;
    }

    if (v < n_tips) continue;

    int m = n_kch[v] - s;
    if (s == 0 && m > 2) n_poly++;

    // Push the second child, then the first, so that the
    // first is emitted first. A child is followed down
    // its chain of single-child nodes.
    for (int which = 1; which >= 0; which--) {
      int    c;
      double cl;
      if (which == 1 && m > 2) {
        c = v; cl = 0;          // The next ladder node.
      }
      else {
        c  = kch[kpos[v] + s + which];
        cl = blen[c];
        while (c >= n_tips && n_kch[c] == 1) {
          c   = kch[kpos[c]];
          cl += blen[c];
          modified = 1;
        }
      }
      top++;
      st_v[top] = c;
      st_s[top] = (which == 1 && m > 2) ? s + 1 : 0;
      st_p[top] = id;
      st_w[top] = which;
      st_l[top] = cl;
    }
  }

  if (k != ni || n_e != 2 * ni)
    error("Internal error in C_philr_tree()."); // # nocov
  if (n_poly) modified = 1;


  // Branch-length weights. Zero-length tip branches
  // get the smallest positive branch length.
  double *blw = REAL(res_blw), *md = REAL(res_md);
  double  min_pos = R_PosInf;
  int     has_na  = 0;
  for (int e = 0; e < n_e; e++) {
    if (ISNAN(elo[e])) has_na = 1;
    else if (elo[e] > 0 && elo[e] < min_pos) min_pos = elo[e];
  }

  if (has_na || !R_FINITE(min_pos)) {
    for (int j = 0; j < ni; j++) blw[j] = md[j] = NA_REAL;
  }
  else {
    for (int j = 0; j < D; j++)
      if (flen[j] == 0) flen[j] = min_pos;

    // Tips below each node, and their summed distances
    // to it; children have larger ids than parents.
    double *n_below = (double *) R_alloc((size_t) 2 * D, sizeof(double));
    double *d_below = (double *) R_alloc((size_t) 2 * D, sizeof(double));
    for (int j = 0; j < D; j++) { n_below[j] = 1; d_below[j] = 0; }

    for (int j = ni - 1; j >= 0; j--) {
      int a = c1[j], b = c2[j];
      blw[j]         = flen[a] + flen[b];
      n_below[D + j] = n_below[a] + n_below[b];
      d_below[D + j] = d_below[a] + n_below[a] * flen[a] + d_below[b] + n_below[b] * flen[b];
      md[j]          = blw[j] + d_below[a] / n_below[a] + d_below[b] / n_below[b];
    }
  }


  const char *nms[] = {
    "c1", "c2", "orig", "ladder", "rung", "edge_par", "edge_child", "edge_len",
    "blw", "md", "n_poly", "modified", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, nms));
  SET_VECTOR_ELT(res,  0, res_c1);
  SET_VECTOR_ELT(res,  1, res_c2);
  SET_VECTOR_ELT(res,  2, res_orig);
  SET_VECTOR_ELT(res,  3, res_lad);
  SET_VECTOR_ELT(res,  4, res_rung);
  SET_VECTOR_ELT(res,  5, res_epar);
  SET_VECTOR_ELT(res,  6, res_echd);
  SET_VECTOR_ELT(res,  7, res_elen);
  SET_VECTOR_ELT(res,  8, res_blw);
  SET_VECTOR_ELT(res,  9, res_md);
  SET_VECTOR_ELT(res, 10, ScalarInteger(n_poly));
  SET_VECTOR_ELT(res, 11, ScalarLogical(modified));

  UNPROTECT(11);
  return res;
}



//======================================================
// R interface: per-feature sums over the samples of
// x = value + pseudocount, for the part weights. A
// D x 4 matrix with columns: sum of x, sum of x^2,
// mean of log(x), and sum of (log(x) - that mean)^2.
// Stored values are summed in sample order, then the
// unstored zeros' share.
//======================================================
SEXP C_philr_stats(SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_pseudocount) {

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int    n  = em->n_samples;
  int    D  = em->n_otus;
  double pc = asReal(sexp_pseudocount);

  SEXP result = PROTECT(allocMatrix(REALSXP, D, 4));
  double *sx  = REAL(result);
  double *sxx = sx  + D;
  double *ml  = sxx + D;
  double *ssl = ml  + D;
  int    *nnz = (int *) R_alloc((size_t) D, sizeof(int));

  for (int j = 0; j < D; j++) { sx[j] = sxx[j] = ml[j] = ssl[j] = 0; nnz[j] = 0; }

  int nz = em->pos_vec[n];
  for (int e = 0; e < nz; e++) {
    int    j = em->otu_vec[e];
    double x = em->val_vec[e] + pc;
    sx[j]  += x;
    sxx[j] += x * x;
    ml[j]  += log(x);
    nnz[j]++;
  }

  double lz = (pc > 0) ? log(pc) : R_NegInf;
  for (int j = 0; j < D; j++) {
    int z = n - nnz[j];
    if (z > 0) {
      sx[j]  += z * pc;
      sxx[j] += z * pc * pc;
      ml[j]  += z * lz;
    }
    ml[j] /= n;
  }

  for (int e = 0; e < nz; e++) {
    int    j = em->otu_vec[e];
    double d = log(em->val_vec[e] + pc) - ml[j];
    ssl[j] += d * d;
  }
  for (int j = 0; j < D; j++) {
    int z = n - nnz[j];
    if (z > 0) { double d = lz - ml[j]; ssl[j] += z * d * d; }
  }

  free_all();
  UNPROTECT(1);
  return result;
}



//======================================================
// The coordinates.
//======================================================
static ecomatrix_t  *ph_em;
static int           ph_margin;
static int           ph_ni;      // Internal nodes (coordinates).
static double        ph_pc;
static const int    *ph_c1;      // Children's ids (see C_philr_tree()).
static const int    *ph_c2;
static const int    *ph_col;     // Each node's column.
static const double *ph_p;       // Part weights, per feature.
static const double *ph_logp;    // Their logs.
static const double *ph_a;       // f sqrt(N_A N_B / (N_A + N_B)), per node.
static const double *ph_ia;      // 1 / N_A,
static const double *ph_ib;      // 1 / N_B.
static double       *ph_scratch; // 2D doubles per thread.
static double       *ph_out;


static void *philr_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D         = ph_em->n_otus;
  int ni        = ph_ni;

  // S[0 .. D-1]: each feature's w; S[D + k]: node k's sum.
  double *S = ph_scratch + (size_t) thread_i * 2 * D;

  for (int s = thread_i; s < ph_em->n_samples; s += n_threads) {

    int e0 = ph_em->pos_vec[s], e1 = ph_em->pos_vec[s + 1];

    // Unstored zeros: log(pseudocount). With pseudocount 0, R has
    // checked that there are none.
    if (e1 - e0 < D) {
      double lz = log(ph_pc);
      for (int o = 0; o < D; o++) S[o] = lz;
    }
    for (int e = e0; e < e1; e++)
      S[ph_em->otu_vec[e]] = log(ph_em->val_vec[e] + ph_pc);
    for (int o = 0; o < D; o++)
      S[o] = ph_p[o] * (S[o] - ph_logp[o]);

    for (int k = ni - 1; k >= 0; k--) {
      double sa = S[ph_c1[k]], sb = S[ph_c2[k]];
      S[D + k]  = sa + sb;
      double v  = ph_a[k] * (sa * ph_ia[k] - sb * ph_ib[k]);
      int    c  = ph_col[k];
      if (ph_margin == 1) ph_out[s + (size_t) ph_em->n_samples * c] = v;
      else                ph_out[c + (size_t) ni * s]               = v;
    }
  }
  return NULL;
}


//======================================================
// R interface: the coordinates as a dense matrix in
// the input's orientation (samples x coordinates for
// margin 1). sexp_c1, sexp_c2: C_philr_tree()'s;
// sexp_col: each node's 0-based column; sexp_p: the
// part weights, per feature; sexp_f: the ILR weights,
// per node (in C_philr_tree()'s order).
//======================================================
SEXP C_philr_coords(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_c1, SEXP sexp_c2,
    SEXP sexp_col, SEXP sexp_p, SEXP sexp_f, SEXP sexp_pseudocount,
    SEXP sexp_n_threads ) {

  int margin = asInteger(sexp_margin);
  int n_thr  = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int ns = em->n_samples;
  int D  = em->n_otus;
  int ni = D - 1;

  if (LENGTH(sexp_c1) != ni || LENGTH(sexp_c2) != ni || LENGTH(sexp_col) != ni ||
      LENGTH(sexp_p) != D || LENGTH(sexp_f) != ni)
    STOP("Invalid arguments to C_philr_coords()."); // # nocov

  const int    *c1 = INTEGER(sexp_c1);
  const int    *c2 = INTEGER(sexp_c2);
  const double *p  = REAL(sexp_p);
  const double *f  = REAL(sexp_f);

  // Each clade's summed part weight, child before parent.
  double *N    = (double *) R_alloc((size_t) 2 * D, sizeof(double));
  double *logp = (double *) R_alloc((size_t) D,     sizeof(double));
  double *a    = (double *) R_alloc((size_t) ni + 1, sizeof(double));
  double *ia   = (double *) R_alloc((size_t) ni + 1, sizeof(double));
  double *ib   = (double *) R_alloc((size_t) ni + 1, sizeof(double));

  for (int o = 0; o < D; o++) { N[o] = p[o]; logp[o] = log(p[o]); }
  for (int k = ni - 1; k >= 0; k--) {
    double na = N[c1[k]], nb = N[c2[k]];
    N[D + k] = na + nb;
    a[k]     = f[k] * sqrt(na * nb / (na + nb));
    ia[k]    = 1 / na;
    ib[k]    = 1 / nb;
  }

  ph_em      = em;
  ph_margin  = margin;
  ph_ni      = ni;
  ph_pc      = asReal(sexp_pseudocount);
  ph_c1      = c1;
  ph_c2      = c2;
  ph_col     = INTEGER(sexp_col);
  ph_p       = p;
  ph_logp    = logp;
  ph_a       = a;
  ph_ia      = ia;
  ph_ib      = ib;
  ph_scratch = (double *) R_alloc((size_t) n_thr * 2 * D, sizeof(double));

  SEXP result = PROTECT(margin == 1 ? allocMatrix(REALSXP, ns, ni) : allocMatrix(REALSXP, ni, ns));
  ph_out = REAL(result);

  run_parallel(philr_worker, n_thr, ns);

  free_all();
  UNPROTECT(1);
  return result;
}
