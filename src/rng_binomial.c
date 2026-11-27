// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Binomial and multinomial draws on PCG (rng_unif() in rng_dist.c), for
// hill_extrapolation()'s bootstrap. Ported from NumPy 2.3.5:
// random_binomial(), random_binomial_btpe(), random_binomial_inversion()
// and random_multinomial() in numpy/random/src/distributions/distributions.c.
// Changes for ecodive: rng_unif() replaces next_double(); the setup that
// NumPy caches between calls with the same n and p (binomial_t) is
// recomputed on each call, as each multinomial category has its own p; and
// a ratio that rounding pushes past 1 in rng_multinom() is clamped. Those
// files are covered by this license:
//
//   Copyright (c) 2005-2025, NumPy Developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//       * Redistributions of source code must retain the above copyright
//          notice, this list of conditions and the following disclaimer.
//
//       * Redistributions in binary form must reproduce the above
//          copyright notice, this list of conditions and the following
//          disclaimer in the documentation and/or other materials provided
//          with the distribution.
//
//       * Neither the name of the NumPy Developers nor the names of any
//          contributors may be used to endorse or promote products derived
//          from this software without specific prior written permission.
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
//   rng_binom      Binomial(n, p): inversion for n min(p, 1 - p) <= 30,
//                  else the BTPE algorithm of Kachitvichyanukul & Schmeiser
//                  (1988).
//   rng_multinom   Multinomial(n, p[0 .. d-1]) by conditional binomials:
//                  O(d) draws rather than O(n).
//
// Kachitvichyanukul, V., & Schmeiser, B. W. (1988). Binomial random variate
// generation. Communications of the ACM, 31, 216-222.

#include "ecodive.h"


// BTPE, for 0 < p <= 0.5 and n p > 30.
static int64_t binom_btpe (rng_t *rng, int64_t n, double p) {

  double r, q, fm, p1, xm, xl, xr, c, laml, lamr, p2, p3, p4;
  double a, u, v, s, F, rho, t, A, nrq, x1, x2, f1, f2, z, z2, w, w2, x;
  int64_t m, y, k, i;

  r    = (p < 1.0 - p) ? p : 1.0 - p;
  q    = 1.0 - r;
  fm   = n * r + r;
  m    = (int64_t) floor(fm);
  p1   = floor(2.195 * sqrt(n * r * q) - 4.6 * q) + 0.5;
  xm   = m + 0.5;
  xl   = xm - p1;
  xr   = xm + p1;
  c    = 0.134 + 20.5 / (15.3 + m);
  a    = (fm - xl) / (fm - xl * r);
  laml = a * (1.0 + a / 2.0);
  a    = (xr - fm) / (xr * q);
  lamr = a * (1.0 + a / 2.0);
  p2   = p1 * (1.0 + 2.0 * c);
  p3   = p2 + c / laml;
  p4   = p3 + c / lamr;

Step10:
  nrq = n * r * q;
  u = rng_unif(rng) * p4;
  v = rng_unif(rng);
  if (u > p1)
    goto Step20;
  y = (int64_t) floor(xm - p1 * v + u);
  goto Step60;

Step20:
  if (u > p2)
    goto Step30;
  x = xl + (u - p1) / c;
  v = v * c + 1.0 - fabs(m - x + 0.5) / p1;
  if (v > 1.0)
    goto Step10;
  y = (int64_t) floor(x);
  goto Step50;

Step30:
  if (u > p3)
    goto Step40;
  y = (int64_t) floor(xl + log(v) / laml);
  if ((y < 0) || (v == 0.0))
    goto Step10;
  v = v * (u - p2) * laml;
  goto Step50;

Step40:
  y = (int64_t) floor(xr - log(v) / lamr);
  if ((y > n) || (v == 0.0))
    goto Step10;
  v = v * (u - p3) * lamr;

Step50:
  k = (y > m) ? y - m : m - y;
  if ((k > 20) && (k < ((nrq) / 2.0 - 1)))
    goto Step52;

  s = r / q;
  a = s * (n + 1);
  F = 1.0;
  if (m < y) {
    for (i = m + 1; i <= y; i++) F *= (a / i - s);
  }
  else if (m > y) {
    for (i = y + 1; i <= m; i++) F /= (a / i - s);
  }
  if (v > F)
    goto Step10;
  goto Step60;

Step52:
  rho = (k / (nrq)) * ((k * (k / 3.0 + 0.625) + 0.16666666666666666) / nrq + 0.5);
  t   = -k * k / (2 * nrq);
  A   = log(v);
  if (A < (t - rho))
    goto Step60;
  if (A > (t + rho))
    goto Step10;

  x1 = y + 1;
  f1 = m + 1;
  z  = n + 1 - m;
  w  = n - y + 1;
  x2 = x1 * x1;
  f2 = f1 * f1;
  z2 = z * z;
  w2 = w * w;
  if (A > (xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) +
           (y - m) * log(w * r / (x1 * q)) +
           (13680. - (462. - (132. - (99. - 140. / f2) / f2) / f2) / f2) / f1 / 166320. +
           (13680. - (462. - (132. - (99. - 140. / z2) / z2) / z2) / z2) / z  / 166320. +
           (13680. - (462. - (132. - (99. - 140. / x2) / x2) / x2) / x2) / x1 / 166320. +
           (13680. - (462. - (132. - (99. - 140. / w2) / w2) / w2) / w2) / w  / 166320.)) {
    goto Step10;
  }

Step60:
  if (p > 0.5) y = n - y;
  return y;
}


// Inversion, for 0 < p <= 0.5 and n p <= 30.
static int64_t binom_inversion (rng_t *rng, int64_t n, double p) {
  double  q     = 1.0 - p;
  double  qn    = exp(n * log(q));
  double  np    = n * p;
  double  b     = np + 10.0 * sqrt(np * q + 1);
  int64_t bound = (int64_t) ((double) n < b ? (double) n : b);
  int64_t X     = 0;
  double  px    = qn;
  double  U     = rng_unif(rng);
  while (U > px) {
    X++;
    if (X > bound) {
      X  = 0;
      px = qn;
      U  = rng_unif(rng);
    }
    else {
      U -= px;
      px = ((n - X + 1) * p * px) / (X * q);
    }
  }
  return X;
}


int64_t rng_binom (rng_t *rng, int64_t n, double p) {
  if (n <= 0 || !(p > 0)) return 0;
  if (p >= 1) return n;
  if (p <= 0.5) {
    return (p * n <= 30.0) ? binom_inversion(rng, n, p) : binom_btpe(rng, n, p);
  }
  double q = 1.0 - p;
  return n - ((q * n <= 30.0) ? binom_inversion(rng, n, q) : binom_btpe(rng, n, q));
}


// n draws over the d categories with probabilities p (summing to 1), into
// out[0 .. d-1].
void rng_multinom (rng_t *rng, int64_t n, const double *p, int d, int64_t *out) {
  double  remaining_p = 1.0;
  int64_t dn = n;
  for (int j = 0; j < d; j++) out[j] = 0;
  for (int j = 0; j < d - 1; j++) {
    double pj = (remaining_p > 0) ? p[j] / remaining_p : 1.0;
    if (pj > 1) pj = 1;
    out[j] = rng_binom(rng, dn, pj);
    dn -= out[j];
    if (dn <= 0) break;
    remaining_p -= p[j];
  }
  if (dn > 0) out[d - 1] = dn;
}
