// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Random draws from continuous distributions, on PCG (pcg_basic.c), for
// methods that sample: dirmult_ttest() so far. Each stream is seeded with
// rng_seed(r, seed, stream), so that a method can give every independent
// unit of work (e.g. one sample of one Monte Carlo draw) its own stream,
// and get the same numbers whichever thread does the work.
//
//   rng_unif       Uniform on (0, 1): 52 random bits, (k + 1/2) / 2^52, so
//                  never exactly 0 or 1. Two PCG outputs per value.
//   rng_norm       Standard normal: Marsaglia's polar method (Marsaglia &
//                  Bray 1964). Each accepted pair gives two values; the
//                  second is kept in the stream's state for the next call.
//   rng_log_gamma  The natural log of a Gamma(shape, 1) draw: Marsaglia &
//                  Tsang (2000) for shape >= 1; for shape < 1, their boost
//                  log G(shape + 1) + log(U) / shape, kept on the log scale
//                  so that tiny shapes never underflow to log(0).
//   rng_gamma      exp() of that for shape < 1; d * v directly otherwise.
//
// Marsaglia, G., & Bray, T. A. (1964). A convenient method for generating
// normal variables. SIAM Review, 6, 260-264.
//
// Marsaglia, G., & Tsang, W. W. (2000). A simple method for generating
// gamma variables. ACM Transactions on Mathematical Software, 26, 363-372.

#include "ecodive.h"


void rng_seed (rng_t *r, uint64_t seed, uint64_t stream) {
  pcg32_srandom_r(&r->pcg, seed, stream);
  r->has_spare = 0;
  r->spare     = 0;
}


double rng_unif (rng_t *r) {
  uint64_t a = pcg32_random_r(&r->pcg) >> 6;    // 26 bits
  uint64_t b = pcg32_random_r(&r->pcg) >> 6;    // 26 bits
  double   k = (double) ((a << 26) | b);        // 0 .. 2^52 - 1, exact
  return (k + 0.5) * 0x1p-52;                   // exact: 53 bits
}


double rng_norm (rng_t *r) {

  if (r->has_spare) {
    r->has_spare = 0;
    return r->spare;
  }

  double u, v, s;
  do {
    u = 2 * rng_unif(r) - 1;
    v = 2 * rng_unif(r) - 1;
    s = u * u + v * v;
  } while (s >= 1 || s == 0);

  double m = sqrt(-2 * log(s) / s);
  r->spare     = v * m;
  r->has_spare = 1;
  return u * m;
}


// Marsaglia and Tsang's G(a), a >= 1, as d * v; returns log(d * v) when
// `log_scale`.
static double mt_gamma (rng_t *r, double a, int log_scale) {

  double d = a - 1.0 / 3.0;
  double c = 1.0 / sqrt(9.0 * d);

  for (;;) {
    double x, v;
    do {
      x = rng_norm(r);
      v = 1 + c * x;
    } while (v <= 0);
    v = v * v * v;

    double u  = rng_unif(r);
    double x2 = x * x;
    if (u < 1 - 0.0331 * (x2 * x2) ||
        log(u) < 0.5 * x2 + d * (1 - v + log(v)))
      return log_scale ? log(d * v) : d * v;
  }
}


double rng_log_gamma (rng_t *r, double shape) {
  if (shape >= 1) return mt_gamma(r, shape, 1);
  double lg = mt_gamma(r, shape + 1, 1);
  return lg + log(rng_unif(r)) / shape;
}


double rng_gamma (rng_t *r, double shape) {
  if (shape >= 1) return mt_gamma(r, shape, 0);
  return exp(rng_log_gamma(r, shape));
}



//======================================================
// R interface, for the tests: n draws from one stream.
//
// sexp_dist:  1 uniform, 2 normal, 3 gamma, 4 log gamma.
// sexp_param: the gamma's shape (ignored otherwise).
//======================================================
SEXP C_rng_dist (
    SEXP sexp_n, SEXP sexp_dist, SEXP sexp_param, SEXP sexp_seed,
    SEXP sexp_stream ) {

  int    n     = asInteger(sexp_n);
  int    dist  = asInteger(sexp_dist);
  double param = asReal(sexp_param);

  rng_t r;
  rng_seed(&r, (uint64_t) asInteger(sexp_seed), (uint64_t) asReal(sexp_stream));

  SEXP    result = PROTECT(allocVector(REALSXP, n));
  double *x      = REAL(result);

  for (int i = 0; i < n; i++)
    switch (dist) {
      case 1:  x[i] = rng_unif(&r);             break;
      case 2:  x[i] = rng_norm(&r);             break;
      case 3:  x[i] = rng_gamma(&r, param);     break;
      default: x[i] = rng_log_gamma(&r, param); break;
    }

  UNPROTECT(1);
  return result;
}
