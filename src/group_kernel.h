// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Weighted within-group sums over a dist object's packed lower
 * triangle, for many permutations of the groups at once: the
 * inner loops of permanova()'s one-factor path (with d_ij^2)
 * and of mrpp() (with d_ij).
 *
 * The permutations go in chunks of GK_W lanes. A chunk holds
 * each sample's group (as a double, so that comparisons
 * vectorize alongside the sums) under each of its GK_W
 * permutations, row-major (n x GK_W). For column j of the
 * triangle (nd = n - j - 1 rows below it), with li pointing at
 * row j + 1's lanes and lj at row j's:
 *
 *   group_chunk:  s[k] += coef[lj[k]] *
 *                         sum_t d[t] * (li[t, k] == lj[k])
 *
 * group_chunk2 does columns j and j + 1 at once: each row's
 * groups are loaded once for both. Every sum runs in the same
 * order as with one column at a time.
 *
 * GCC and clang get explicit two-lane vectors: compilers don't
 * vectorize these loops at R's default -O2. Each lane does the
 * same operations in the same order as the scalar fallback, so
 * results don't depend on which is used.
 */

#ifndef ECODIVE_GROUP_KERNEL_H
#define ECODIVE_GROUP_KERNEL_H

#include "ecodive.h"

#define GK_W 8  // Lanes (permutations) per chunk.


#if defined(__GNUC__) || defined(__clang__)

typedef double    v2d __attribute__((vector_size(16)));
typedef long long v2l __attribute__((vector_size(16)));

static inline v2d load2(const double *p) { v2d x; memcpy(&x, p, sizeof(x)); return x; }

// Repeat F(h) for each of GK_W / 2 lane pairs.
#define V2_REP(F) F(0) F(1) F(2) F(3)

static inline void group_chunk(
    const double *d, const double *li, const double *lj, int nd, double *s,
    const double *coef) {

  #define G_INIT(h) v2d acc##h = { 0, 0 }, g##h = load2(lj + 2 * h);
  #define G_STEP(h) acc##h += (v2d) ((v2l) a & (v2l) (load2(li + 2 * h) == g##h));
  #define G_DONE(h) s[2 * h]     += acc##h[0] * coef[(int) lj[2 * h]]; \
                    s[2 * h + 1] += acc##h[1] * coef[(int) lj[2 * h + 1]];

  V2_REP(G_INIT)
  for (int t = 0; t < nd; t++, li += GK_W) {
    v2d a = { d[t], d[t] };
    V2_REP(G_STEP)
  }
  V2_REP(G_DONE)

  #undef G_INIT
  #undef G_STEP
  #undef G_DONE
}

static inline void group_chunk2(
    const double *d0, const double *d1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s,
    const double *coef) {

  #define G2_INIT(h) v2d acc##h = { 0, 0 }, bcc##h = { 0, 0 }, \
                     g##h = load2(lj0 + 2 * h), k##h = load2(lj1 + 2 * h);
  #define G2_FIRST(h) acc##h += (v2d) ((v2l) a & (v2l) (load2(li + 2 * h) == g##h));
  #define G2_STEP(h) { v2d x = load2(li + 2 * h); \
                       acc##h += (v2d) ((v2l) a & (v2l) (x == g##h)); \
                       bcc##h += (v2d) ((v2l) b & (v2l) (x == k##h)); }
  #define G2_DONE(h) \
    s[2 * h]     += acc##h[0] * coef[(int) lj0[2 * h]];     \
    s[2 * h + 1] += acc##h[1] * coef[(int) lj0[2 * h + 1]]; \
    s[2 * h]     += bcc##h[0] * coef[(int) lj1[2 * h]];     \
    s[2 * h + 1] += bcc##h[1] * coef[(int) lj1[2 * h + 1]];

  V2_REP(G2_INIT)
  { v2d a = { d0[0], d0[0] }; V2_REP(G2_FIRST) }
  li += GK_W;
  for (int t = 1; t < nd; t++, li += GK_W) {
    v2d a = { d0[t],     d0[t]     };
    v2d b = { d1[t - 1], d1[t - 1] };
    V2_REP(G2_STEP)
  }
  V2_REP(G2_DONE)

  #undef G2_INIT
  #undef G2_FIRST
  #undef G2_STEP
  #undef G2_DONE
}

#else

#define GK_REP(F) F(0) F(1) F(2) F(3) F(4) F(5) F(6) F(7)

static inline void group_chunk(
    const double *d, const double *li, const double *lj, int nd, double *s,
    const double *coef) {

  #define G_INIT(k) double acc##k = 0, g##k = lj[k];
  #define G_STEP(k) acc##k += (li[k] == g##k) ? a : 0.0;
  #define G_DONE(k) s[k] += acc##k * coef[(int) g##k];

  GK_REP(G_INIT)
  for (int t = 0; t < nd; t++, li += GK_W) {
    double a = d[t];
    GK_REP(G_STEP)
  }
  GK_REP(G_DONE)

  #undef G_INIT
  #undef G_STEP
  #undef G_DONE
}

static inline void group_chunk2(
    const double *d0, const double *d1, const double *li,
    const double *lj0, const double *lj1, int nd, double *s,
    const double *coef) {
  group_chunk(d0, li, lj0, nd, s, coef);
  group_chunk(d1, li + GK_W, lj1, nd - 1, s, coef);
}

#endif

#endif
