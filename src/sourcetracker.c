// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The Gibbs sampler in this file is ported from SourceTracker2 2.0.1-dev
// (github.com/biota/sourcetracker2, master): `gibbs_sampler()`,
// `ConditionalProbability` and `collapse_source_data()` in
// sourcetracker/_sourcetracker.py, and the pipeline in `gibbs_helper()`
// in sourcetracker/_gibbs.py. It is restructured for ecodive (one random
// stream per sink and restart, threads, no Python), as described in
// plans/source_tracking.md. Those files are covered by this license:
//
//   Copyright (c) 2016, Biota Technology.
//   All rights reserved.
//
//   www.biota.com
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
//   * Neither the name of SourceTracker2 nor the names of its
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


/*
 * sourcetracker(): Bayesian microbial source tracking (Knights et al.
 * 2011) by SourceTracker2's Gibbs sampler.
 *
 * E environments (collapsed sources), F features, V = E + 1 states (the
 * environments, then the unknown source; under leave-one-out V = E, the
 * held-out environment being the sink).
 *
 * Random numbers. Every independent piece of work has its own PCG stream
 * of `seed`, so that results are the same for any number of threads, and
 * a sink's result doesn't depend on the other sinks:
 *
 *   rarefying sink row r of `counts`   stream r, exactly as rarefy() does
 *   rarefying environment e            ST_SRC  | e
 *   the Gibbs chain of sink row r,     ST_SINK | r << 24 | restart
 *     restart `restart`
 *   ... of held-out environment e      ST_LOO  | e << 24 | restart
 *
 * A chain draws, in this order: each read's initial state (a bounded
 * integer), then, per pass, a Fisher-Yates shuffle of the read order and
 * one uniform (rng_unif()) per read.
 */

#include "ecodive.h"

#define ST_SRC  ((uint64_t) 1 << 56)
#define ST_SINK ((uint64_t) 2 << 56)
#define ST_LOO  ((uint64_t) 3 << 56)


//======================================================
// Rarefaction of one sample, given as its non-zero
// counts in feature order: Knuth's selection sampling,
// the same draws as rarefy() makes for a sample of a
// table (src/rarefy.c), so the same result for the
// same stream. Rarefies in place; `cnt` holds whole
// numbers. Does nothing when the sample has no more
// than `target` observations.
//======================================================
static void st_rarefy (double *cnt, int k, uint32_t target, uint64_t seed, uint64_t stream) {

  double total = 0;
  for (int i = 0; i < k; i++) total += cnt[i];
  if (!(total > target)) return;

  uint32_t depth = (uint32_t) total;
  pcg32_random_t rng;
  pcg32_srandom_r(&rng, seed, stream);

  uint32_t tried = 0, kept = 0;
  for (int i = 0; i < k; i++) {
    double res = 0;
    for (uint32_t s = 0; s < (uint32_t) cnt[i] && kept < target; s++) {
      uint32_t not_tried  = depth - tried;
      uint32_t still_need = target - kept;
      if (pcg32_random_r(&rng) % not_tried < still_need) { res++; kept++; }
      tried++;
    }
    cnt[i] = res;
  }
}



//======================================================
// The fixed part of the conditional probability,
// SourceTracker2's `known_p_tv` (without its common
// denominator n - 1 + beta V):
//   q[f * E + e] = (m_ef + alpha1) / (m_e + F alpha1)
// from the E x F (column-major) table `src`.
//======================================================
static void st_fixed (const double *src, int E, int F, double alpha1, double *q) {
  for (int e = 0; e < E; e++) {
    double m = 0;
    for (int f = 0; f < F; f++) m += src[e + (size_t) E * f];
    double den = m + F * alpha1;
    for (int f = 0; f < F; f++)
      q[(size_t) f * E + e] = (src[e + (size_t) E * f] + alpha1) / den;
  }
}



//======================================================
// The unnormalized weights of the V states for a read
// of feature f, removed from its state: known state v
// (environment env[v] of q) q_f,env[v] (c_v + beta), the
// unknown (u_f + alpha2 n) / (U + alpha2 n F) (c_U + beta).
// Writes their running sums into cs; returns the total.
//======================================================
static inline double st_weights (
    const double *qf, const int *env, int V, const int *c, double beta,
    double uf, double U, double a2n, double a2nF, double *cs ) {

  double sum = 0;
  for (int v = 0; v < V - 1; v++) {
    sum  += qf[env[v]] * (c[v] + beta);
    cs[v] = sum;
  }
  sum += (uf + a2n) * (c[V - 1] + beta) / (U + a2nF);
  cs[V - 1] = sum;
  return sum;
}



//======================================================
// Shared state of a run.
//======================================================
static int           st_E, st_F, st_V, st_n_sinks, st_n_restarts, st_n_passes;
static int           st_n_draws;   // draws per restart
static int           st_loo, st_features;
static double        st_beta, st_alpha2;
static uint64_t      st_seed;
static const double *st_q;         // F x E, see st_fixed()
static const int    *st_pass_rec;  // n_passes: 1 if a draw follows the pass
static const int    *st_sink_id;   // n_sinks: row of `counts` or environment
static const int    *st_spos;      // n_sinks + 1: CSR offsets of the sinks
static const int    *st_sfeat;     // their features (0-based)
static const int    *st_scnt;      // and counts
static const int    *st_sn;        // n_sinks: reads per sink
static int          *st_draws;     // V x draws x n_sinks
static double       *st_facc;      // per unit: V x (sink's features)
static const size_t *st_fpos;      // per unit: offset into st_facc
static int           st_unit0, st_unit1;  // the block of units to run

// Per-thread scratch, sized for the largest sink.
static int    *st_e_buf, *st_t_buf, *st_o_buf, *st_u_buf, *st_c_buf, *st_env_buf;
static double *st_cs_buf;
static int     st_max_n, st_max_k;


static void st_chain (int unit, int thread) {

  int V = st_V, E = st_E;
  int j = unit / st_n_restarts, restart = unit % st_n_restarts;
  int n = st_sn[j];
  int k = st_spos[j + 1] - st_spos[j];
  const int *feat = st_sfeat + st_spos[j];
  const int *cnt  = st_scnt  + st_spos[j];

  int    *e   = st_e_buf   + (size_t) thread * st_max_n;  // each read's state
  int    *t   = st_t_buf   + (size_t) thread * st_max_n;  // its feature (local)
  int    *ord = st_o_buf   + (size_t) thread * st_max_n;  // the read order
  int    *u   = st_u_buf   + (size_t) thread * st_max_k;  // unknown, by feature
  int    *c   = st_c_buf   + (size_t) thread * V;         // reads per state
  int    *env = st_env_buf + (size_t) thread * V;         // state -> env of q
  double *cs  = st_cs_buf  + (size_t) thread * V;

  // The known states' environments: all of them, or all but the one held
  // out.
  int held = st_loo ? st_sink_id[j] : -1;
  for (int v = 0, x = 0; x < E; x++) if (x != held) env[v++] = x;

  // The reads, in feature order (numpy.repeat()).
  for (int i = 0, s = 0; i < k; i++)
    for (int r = 0; r < cnt[i]; r++) t[s++] = i;

  rng_t rng;
  uint64_t stream = (st_loo ? ST_LOO : ST_SINK) |
                    ((uint64_t) st_sink_id[j] << 24) | (uint64_t) restart;
  rng_seed(&rng, st_seed, stream);

  // A uniform random start.
  for (int v = 0; v < V; v++) c[v] = 0;
  for (int i = 0; i < k; i++) u[i] = 0;
  int U = 0;
  for (int s = 0; s < n; s++) {
    int v = (int) pcg32_bounded(&rng.pcg, (uint32_t) V);
    e[s] = v;
    c[v]++;
    if (v == V - 1) { u[t[s]]++; U++; }
    ord[s] = s;
  }

  double a2n  = st_alpha2 * n;
  double a2nF = a2n * st_F;
  double beta = st_beta;
  const double *q = st_q;
  int  *draw  = st_draws + ((size_t) j * st_n_restarts * st_n_draws +
                            (size_t) restart * st_n_draws) * V;
  double *facc = st_features ? st_facc + st_fpos[unit] : NULL;

  for (int pass = 0; pass < st_n_passes; pass++) {

    for (int i = n - 1; i > 0; i--) {
      int r   = (int) pcg32_bounded(&rng.pcg, (uint32_t) i + 1);
      int tmp = ord[i]; ord[i] = ord[r]; ord[r] = tmp;
    }

    for (int i = 0; i < n; i++) {
      int s  = ord[i];
      int v  = e[s];
      int ti = t[s];

      c[v]--;
      if (v == V - 1) { u[ti]--; U--; }

      double tot = st_weights(
        q + (size_t) feat[ti] * E, env, V, c, beta, u[ti], U, a2n, a2nF, cs );

      // numpy.searchsorted(cs, w), side 'left': the first v with cs[v] >= w.
      double w = rng_unif(&rng) * tot;
      v = 0;
      while (v < V - 1 && cs[v] < w) v++;

      e[s] = v;
      c[v]++;
      if (v == V - 1) { u[ti]++; U++; }
    }

    if (st_pass_rec[pass]) {
      memcpy(draw, c, (size_t) V * sizeof(int));
      draw += V;
      if (facc)
        for (int s = 0; s < n; s++) facc[(size_t) e[s] * k + t[s]] += 1;
    }
  }
}


static void *st_worker (void *arg) {
  int i = ((worker_t *) arg)->i;
  int n = ((worker_t *) arg)->n;
  for (int unit = st_unit0 + i; unit < st_unit1; unit += n)
    st_chain(unit, i);
  return NULL;
}



//======================================================
// R interface: collapse the source samples by
// environment, with the mean, truncated to a whole
// number (SourceTracker2's collapse_source_data(...,
// 'mean'), then astype(int32)).
//
// sexp_env: each sample's environment, 0 .. E - 1, or
//           -1 (not a source).
// Returns the E x F matrix.
//======================================================
SEXP C_st_collapse (SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_env, SEXP sexp_n_envs) {

  int E = asInteger(sexp_n_envs);
  const int *env = INTEGER(sexp_env);

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int n = em->n_samples, F = em->n_otus;

  SEXP sexp_res = PROTECT(allocMatrix(REALSXP, E, F));
  double *res = REAL(sexp_res);
  memset(res, 0, (size_t) E * F * sizeof(double));

  int *n_env = (int *) R_alloc((size_t) E + 1, sizeof(int));
  memset(n_env, 0, ((size_t) E + 1) * sizeof(int));

  for (int s = 0; s < n; s++) {
    int x = env[s];
    if (x < 0) continue;
    n_env[x]++;
    for (int p = em->pos_vec[s]; p < em->pos_vec[s + 1]; p++)
      res[x + (size_t) E * em->otu_vec[p]] += em->val_vec[p];
  }
  free_all();

  for (int f = 0; f < F; f++)
    for (int x = 0; x < E; x++) {
      double *r = res + x + (size_t) E * f;
      *r = n_env[x] ? trunc(*r / n_env[x]) : 0;
    }

  UNPROTECT(1);
  return sexp_res;
}



//======================================================
// R interface: the sampler.
//
// sexp_src:        The E x F collapsed sources.
// sexp_sink_rows:  The sinks' rows of `counts` (0-based);
//                  ignored with `loo`.
// sexp_depth,      Rarefaction depths of the sources and
// sexp_sink_depth: sinks; 0 for none.
// sexp_pass_rec:   One per pass: 1 if a draw is taken
//                  after it.
//
// Returns list(draws, features, sources, sinks):
//   draws     V x (n_restarts x draws) x n_sinks counts of
//             reads per state, restart by restart.
//   features  NULL, or per sink the V x F reads per state
//             and feature, summed over the draws.
//   sources   The rarefied E x F sources.
//   sinks     NULL, or (with `details`) the rarefied sinks,
//             n_sinks x F.
//======================================================
SEXP C_sourcetracker (
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_sink_rows, SEXP sexp_src,
    SEXP sexp_depth, SEXP sexp_sink_depth, SEXP sexp_alpha1, SEXP sexp_alpha2,
    SEXP sexp_beta, SEXP sexp_n_restarts, SEXP sexp_pass_rec, SEXP sexp_loo,
    SEXP sexp_features, SEXP sexp_seed, SEXP sexp_details, SEXP sexp_cpus ) {

  int      E          = nrows(sexp_src);
  int      F          = ncols(sexp_src);
  int      depth      = asInteger(sexp_depth);
  int      sink_depth = asInteger(sexp_sink_depth);
  double   alpha1     = asReal(sexp_alpha1);
  int      n_restarts = asInteger(sexp_n_restarts);
  int      n_passes   = length(sexp_pass_rec);
  int      loo        = asLogical(sexp_loo);
  int      features   = asLogical(sexp_features);
  int      details    = asLogical(sexp_details);
  int      n_thr      = asInteger(sexp_cpus);
  uint64_t seed       = (uint64_t) asInteger(sexp_seed);
  if (n_thr < 1) n_thr = 1;

  const int *pass_rec = INTEGER(sexp_pass_rec);
  int n_draws = 0;
  for (int p = 0; p < n_passes; p++) n_draws += pass_rec[p] != 0;

  int V = loo ? E : E + 1;


  //--------------------------------------------------------
  // Rarefied sources and the fixed part.
  //--------------------------------------------------------
  SEXP sexp_sources = PROTECT(allocMatrix(REALSXP, E, F));
  double *src = REAL(sexp_sources);
  memcpy(src, REAL(sexp_src), (size_t) E * F * sizeof(double));

  if (depth > 0) {
    double *row = (double *) R_alloc((size_t) F + 1, sizeof(double));
    int    *idx = (int *)    R_alloc((size_t) F + 1, sizeof(int));
    for (int x = 0; x < E; x++) {
      int k = 0;
      for (int f = 0; f < F; f++) {
        double v = src[x + (size_t) E * f];
        if (v > 0) { idx[k] = f; row[k++] = v; }
      }
      st_rarefy(row, k, (uint32_t) depth, seed, ST_SRC | (uint64_t) x);
      for (int i = 0; i < k; i++) src[x + (size_t) E * idx[i]] = row[i];
    }
  }

  double *q = (double *) R_alloc((size_t) F * E + 1, sizeof(double));
  st_fixed(src, E, F, alpha1, q);


  //--------------------------------------------------------
  // The sinks, as rarefied non-zero counts in feature
  // order (CSR).
  //--------------------------------------------------------
  int  n_sinks;
  int *sink_id, *spos, *sfeat, *scnt, *sn;

  if (loo) {
    n_sinks = E;
    sink_id = (int *) R_alloc((size_t) E + 1, sizeof(int));
    spos    = (int *) R_alloc((size_t) E + 1, sizeof(int));
    sn      = (int *) R_alloc((size_t) E + 1, sizeof(int));
    size_t nnz = 0;
    for (int x = 0; x < E; x++)
      for (int f = 0; f < F; f++) nnz += src[x + (size_t) E * f] > 0;
    sfeat = (int *) R_alloc(nnz + 1, sizeof(int));
    scnt  = (int *) R_alloc(nnz + 1, sizeof(int));
    int p = 0;
    for (int x = 0; x < E; x++) {
      sink_id[x] = x;
      spos[x]    = p;
      sn[x]      = 0;
      for (int f = 0; f < F; f++) {
        double v = src[x + (size_t) E * f];
        if (v > 0) { sfeat[p] = f; scnt[p++] = (int) v; sn[x] += (int) v; }
      }
    }
    spos[E] = p;
  }

  else {
    n_sinks = length(sexp_sink_rows);
    sink_id = INTEGER(sexp_sink_rows);
    spos    = (int *) R_alloc((size_t) n_sinks + 1, sizeof(int));
    sn      = (int *) R_alloc((size_t) n_sinks + 1, sizeof(int));

    init_n_ptrs(10);
    ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

    size_t nnz = 0;
    for (int j = 0; j < n_sinks; j++)
      nnz += (size_t) (em->pos_vec[sink_id[j] + 1] - em->pos_vec[sink_id[j]]);
    sfeat = (int *)    R_alloc(nnz + 1, sizeof(int));
    scnt  = (int *)    R_alloc(nnz + 1, sizeof(int));
    double *row = (double *) R_alloc(nnz + 1, sizeof(double));

    int p = 0;
    for (int j = 0; j < n_sinks; j++) {
      int s = sink_id[j], k = 0;
      double *r = row + p;
      int    *g = sfeat + p;
      // Entries come in feature order; a triplet matrix's repeated entries
      // are adjacent, and are summed.
      for (int a = em->pos_vec[s]; a < em->pos_vec[s + 1]; a++) {
        if (k && g[k - 1] == em->otu_vec[a]) r[k - 1] += em->val_vec[a];
        else { g[k] = em->otu_vec[a]; r[k++] = em->val_vec[a]; }
      }
      if (sink_depth > 0)
        st_rarefy(r, k, (uint32_t) sink_depth, seed, (uint64_t) s);
      spos[j] = p;
      sn[j]   = 0;
      for (int i = 0; i < k; i++)
        if (r[i] > 0) {
          sfeat[p] = g[i];
          scnt[p]  = (int) r[i];
          sn[j]   += (int) r[i];
          p++;
        }
    }
    spos[n_sinks] = p;
    free_all();
  }


  //--------------------------------------------------------
  // Outputs and scratch.
  //--------------------------------------------------------
  int n_units  = n_sinks * n_restarts;
  int D        = n_restarts * n_draws;

  SEXP sexp_draws = PROTECT(alloc3DArray(INTSXP, V, D, n_sinks));
  memset(INTEGER(sexp_draws), 0, (size_t) V * D * n_sinks * sizeof(int));

  int max_n = 1, max_k = 1;
  for (int j = 0; j < n_sinks; j++) {
    if (sn[j] > max_n) max_n = sn[j];
    if (spos[j + 1] - spos[j] > max_k) max_k = spos[j + 1] - spos[j];
  }

  size_t *fpos = NULL;
  double *facc = NULL;
  if (features) {
    fpos = (size_t *) R_alloc((size_t) n_units + 1, sizeof(size_t));
    size_t tot = 0;
    for (int u = 0; u < n_units; u++) {
      int j = u / n_restarts;
      fpos[u] = tot;
      tot += (size_t) V * (spos[j + 1] - spos[j]);
    }
    facc = (double *) R_alloc(tot + 1, sizeof(double));
    memset(facc, 0, (tot + 1) * sizeof(double));
  }

  // Threads only for enough work: about 10 ns per read per state per pass.
  double work = 0;
  for (int j = 0; j < n_sinks; j++) work += (double) sn[j];
  work *= (double) n_restarts * n_passes * V * 10;
  int nt = threads_for(n_thr, work, 1e6);
  if (nt > n_units) nt = n_units < 1 ? 1 : n_units;

  st_e_buf   = (int *)    R_alloc((size_t) nt * max_n, sizeof(int));
  st_t_buf   = (int *)    R_alloc((size_t) nt * max_n, sizeof(int));
  st_o_buf   = (int *)    R_alloc((size_t) nt * max_n, sizeof(int));
  st_u_buf   = (int *)    R_alloc((size_t) nt * max_k, sizeof(int));
  st_c_buf   = (int *)    R_alloc((size_t) nt * V,     sizeof(int));
  st_env_buf = (int *)    R_alloc((size_t) nt * V,     sizeof(int));
  st_cs_buf  = (double *) R_alloc((size_t) nt * V,     sizeof(double));
  st_max_n   = max_n;
  st_max_k   = max_k;

  st_E          = E;
  st_F          = F;
  st_V          = V;
  st_n_sinks    = n_sinks;
  st_n_restarts = n_restarts;
  st_n_passes   = n_passes;
  st_n_draws    = n_draws;
  st_loo        = loo;
  st_features   = features;
  st_beta       = asReal(sexp_beta);
  st_alpha2     = asReal(sexp_alpha2);
  st_seed       = seed;
  st_q          = q;
  st_pass_rec   = pass_rec;
  st_sink_id    = sink_id;
  st_spos       = spos;
  st_sfeat      = sfeat;
  st_scnt       = scnt;
  st_sn         = sn;
  st_draws      = INTEGER(sexp_draws);
  st_facc       = facc;
  st_fpos       = fpos;


  //--------------------------------------------------------
  // Run the chains, a block at a time so that the user can
  // interrupt. Each chain writes only its own outputs.
  //--------------------------------------------------------
  double unit_work = n_units ? work / n_units : 0;
  int    block     = unit_work > 0 ? (int) (1e9 / unit_work) : n_units;
  if (block < nt)  block = nt;
  if (block < 1)   block = 1;

  for (int u0 = 0; u0 < n_units; u0 += block) {
    st_unit0 = u0;
    st_unit1 = u0 + block < n_units ? u0 + block : n_units;
    run_parallel(st_worker, nt, nt > 1 ? 100 : 1);
    if (st_unit1 < n_units && user_interrupted()) error("Interrupted."); // # nocov
  }


  //--------------------------------------------------------
  // Per-feature contributions: each sink's V x F reads,
  // summed over its restarts in order.
  //--------------------------------------------------------
  SEXP sexp_feat = R_NilValue;
  if (features) {
    sexp_feat = PROTECT(allocVector(VECSXP, n_sinks));
    for (int j = 0; j < n_sinks; j++) {
      SEXP m = allocMatrix(REALSXP, V, F);
      SET_VECTOR_ELT(sexp_feat, j, m);
      double *out = REAL(m);
      memset(out, 0, (size_t) V * F * sizeof(double));
      int k = spos[j + 1] - spos[j];
      const int *feat = sfeat + spos[j];
      for (int r = 0; r < n_restarts; r++) {
        const double *a = facc + fpos[j * n_restarts + r];
        for (int v = 0; v < V; v++)
          for (int i = 0; i < k; i++)
            out[v + (size_t) V * feat[i]] += a[(size_t) v * k + i];
      }
    }
  }
  else PROTECT(sexp_feat);

  SEXP sexp_sinks = R_NilValue;
  if (details) {
    sexp_sinks = PROTECT(allocMatrix(REALSXP, n_sinks, F));
    double *m = REAL(sexp_sinks);
    memset(m, 0, (size_t) n_sinks * F * sizeof(double));
    for (int j = 0; j < n_sinks; j++)
      for (int p = spos[j]; p < spos[j + 1]; p++)
        m[j + (size_t) n_sinks * sfeat[p]] = scnt[p];
  }
  else PROTECT(sexp_sinks);

  SEXP sexp_res = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(sexp_res, 0, sexp_draws);
  SET_VECTOR_ELT(sexp_res, 1, sexp_feat);
  SET_VECTOR_ELT(sexp_res, 2, sexp_sources);
  SET_VECTOR_ELT(sexp_res, 3, sexp_sinks);
  SEXP nms = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, mkChar("draws"));
  SET_STRING_ELT(nms, 1, mkChar("features"));
  SET_STRING_ELT(nms, 2, mkChar("sources"));
  SET_STRING_ELT(nms, 3, mkChar("sinks"));
  setAttrib(sexp_res, R_NamesSymbol, nms);

  UNPROTECT(6);
  return sexp_res;
}



//======================================================
// R interface, for the tests: SourceTracker2's
// ConditionalProbability(alpha1, alpha2, beta, src),
// set_n(n), precalculate(), then
// calculate_cp_slice(xi, m_xiV, m_V, n_vnoti), with its
// denominator n - 1 + beta V. `xi` is 0-based.
//======================================================
SEXP C_st_cond_prob (
    SEXP sexp_src, SEXP sexp_alpha1, SEXP sexp_alpha2, SEXP sexp_beta,
    SEXP sexp_n, SEXP sexp_xi, SEXP sexp_m_xiV, SEXP sexp_m_V, SEXP sexp_c ) {

  int    E = nrows(sexp_src), F = ncols(sexp_src), V = E + 1;
  double n = asReal(sexp_n), beta = asReal(sexp_beta), a2 = asReal(sexp_alpha2);

  double *q   = (double *) R_alloc((size_t) F * E + 1, sizeof(double));
  int    *env = (int *)    R_alloc((size_t) V, sizeof(int));
  double *cs  = (double *) R_alloc((size_t) V, sizeof(double));
  st_fixed(REAL(sexp_src), E, F, asReal(sexp_alpha1), q);
  for (int v = 0; v < E; v++) env[v] = v;

  st_weights(
    q + (size_t) asInteger(sexp_xi) * E, env, V, INTEGER(sexp_c), beta,
    asReal(sexp_m_xiV), asReal(sexp_m_V), a2 * n, a2 * n * F, cs );

  double den = n - 1 + beta * V;
  SEXP sexp_res = PROTECT(allocVector(REALSXP, V));
  double *res = REAL(sexp_res);
  for (int v = 0; v < V; v++) res[v] = (cs[v] - (v ? cs[v - 1] : 0)) / den;
  UNPROTECT(1);
  return sexp_res;
}
