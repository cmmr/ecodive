// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Non-metric multidimensional scaling (NMDS).
 *
 * A port of Peter Minchin's monoMDS (vegan's monoMDS.f, global
 * model) that follows its arithmetic step for step: the same
 * stress formulas, tie treatments, monotone regression, step
 * sizes, back-ups and stopping rules, evaluated in the same
 * order, so that each start gives the same result as
 * vegan::monoMDS() from the same initial configuration.
 *
 * Memory. monoMDS keeps, per start, the dissimilarities and
 * their i and j indices sorted by dissimilarity, plus the
 * configuration's distances, the fitted distances, and a work
 * array. Here the sorted order is one shared, read-only array
 * of 4-byte pair codes (i << 16 | j), made once by
 * C_nmds_prep(). Configuration distances are recomputed from
 * the coordinates whenever they're needed, which is cheap for a
 * few dimensions. Each running start keeps only its fitted
 * distances (8 bytes per pair), its regression blocks' sizes
 * (4 bytes per pair) and, with weak ties, its own order of the
 * pairs within groups of tied dissimilarities (4 bytes per tied
 * pair).
 *
 * Threads. Starts run in parallel, one per thread, and each
 * start is computed entirely by one thread, in the same order
 * as a serial run. So results don't depend on the number of
 * threads. Threads take the next unstarted start from a
 * mutex-protected counter, which only affects which thread runs
 * which start, not the result.
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the output slots that each start owns.
//------------------------------------------------------
static int             nm_n;       // Samples.
static int             nm_k;       // Dimensions.
static size_t          nm_ndis;    // Pairs: n * (n - 1) / 2.
static const uint32_t *nm_codes;   // Pairs, sorted by dissimilarity.
static int             nm_ngrp;    // Groups of tied dissimilarities.
static const int      *nm_gstart;  // Each group's first position.
static const int      *nm_glen;    // Each group's length (>= 2).
static size_t          nm_ntied;   // Total of nm_glen.

static int    nm_weak;     // 1 = weak (primary) ties, 0 = strong.
static int    nm_isform;   // Kruskal's stress formula: 1 or 2.
static int    nm_maxit;
static double nm_smin, nm_sfgrmin, nm_sratmax;

static int           nm_m;       // Starts.
static const double *nm_init;    // n x k x m initial configurations, or
static uint64_t      nm_seed;    // ... NULL for random ones from this
static int           nm_first;   // ... seed, numbered from nm_first.
static double       *nm_points;  // n x k x m final configurations.
static double       *nm_stress;  // m
static int          *nm_iters;   // m
static int          *nm_icause;  // m

// Next start to run.
static int nm_next;
#ifdef HAVE_PTHREAD
  static pthread_mutex_t nm_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif


// Per-thread work space.
typedef struct {
  double   *x;       // n x k, row-major.
  double   *grad;    // n x k, row-major.
  double   *grlast;  // n x k, row-major.
  double   *dhat;    // ndis
  int      *iwork;   // ndis
  uint32_t *tcodes;  // ntied (weak ties only)
  int      *gtoff;   // ngrp (weak ties only)
} work_t;

static work_t *nm_work;


//------------------------------------------------------
// Shell sort of x[1..n] ascending, permuting c[1..n]
// alongside. A port of monoMDS's ASORT4, which sorts in
// place; the order it leaves tied values in depends on
// their order before, and is kept for the next call.
//------------------------------------------------------
static void asort4(double *x, uint32_t *c, size_t n) {

  if (n < 2) return;

  x--; c--; // 1-based, as in Fortran.

  // FN=REAL(N); NLOOPS=MAX(NINT(LOG(FN)/LOG(2.)),1)
  float  fn     = (float) n;
  long   nloops = lroundf(logf(fn) / logf(2.0f));
  if (nloops < 1) nloops = 1;
  size_t m      = (size_t) 1 << (nloops - 1);

  for (long ii = 1; ii <= nloops; ii++) {

    float  fm   = (float) m;
    size_t last = (n > m + 1) ? n - m : 1;

    for (size_t i = 1; i <= last; i++) {
      if (x[i] > x[i + m]) {

        double   xtemp = x[i + m];
        uint32_t ctemp = c[i + m];
        x[i + m] = x[i];
        c[i + m] = c[i];

        if (i <= m) {
          x[i] = xtemp;
          c[i] = ctemp;
        }
        else {
          ptrdiff_t j;
          for (j = (ptrdiff_t) (i - m); j >= 1; j -= (ptrdiff_t) m) {
            if (xtemp >= x[j]) break;
            x[j + m] = x[j];
            c[j + m] = c[j];
          }
          x[j + m] = xtemp;
          c[j + m] = ctemp;
        }
      }
    }

    m = (size_t) (int) (fm / 2.0f);
  }
}


//------------------------------------------------------
// Stable sort of x[0..n-1] ascending, permuting c alongside:
// an LSD radix sort on the bits of x, in four 16-bit
// passes. For non-negative doubles, the order of their bits
// is the order of their values (-0.0 counts as 0). Tied
// values keep their order. O(n), where ASORT4 grows much
// faster: 31 s for 4,000 samples (8 million pairs).
//------------------------------------------------------
static void radix_sort(double *x, uint32_t *c, size_t n) {

  uint64_t *k1    = (uint64_t *) R_alloc(n, sizeof(uint64_t));
  uint64_t *k2    = (uint64_t *) R_alloc(n, sizeof(uint64_t));
  uint32_t *c1    = c;
  uint32_t *c2    = (uint32_t *) R_alloc(n, sizeof(uint32_t));
  size_t   *count = (size_t *)   R_alloc(65537, sizeof(size_t));

  for (size_t i = 0; i < n; i++) {
    uint64_t b;
    memcpy(&b, x + i, sizeof b);
    k1[i] = (b >> 63) ? 0 : b;
  }

  for (int shift = 0; shift < 64; shift += 16) {
    memset(count, 0, 65537 * sizeof(size_t));
    for (size_t i = 0; i < n; i++) count[((k1[i] >> shift) & 0xFFFF) + 1]++;
    for (int d = 0; d < 65536; d++) count[d + 1] += count[d];
    for (size_t i = 0; i < n; i++) {
      size_t dst = count[(k1[i] >> shift) & 0xFFFF]++;
      k2[dst] = k1[i];
      c2[dst] = c1[i];
    }
    uint64_t *kt = k1; k1 = k2; k2 = kt;
    uint32_t *ct = c1; c1 = c2; c2 = ct;
  }

  // After four passes, the sorted codes are back in c.
  for (size_t i = 0; i < n; i++) memcpy(x + i, k1 + i, sizeof(double));
}


//------------------------------------------------------
// Distance between the pair coded `c` in configuration x
// (row-major, k columns), summed as monoMDS's CLCDIS.
//------------------------------------------------------
static inline double pair_dist(const double *x, int k, uint32_t c) {
  const double *xi = x + (size_t) (c >> 16)     * k;
  const double *xj = x + (size_t) (c & 0xFFFFu) * k;
  double s = 0;
  for (int d = 0; d < k; d++) {
    double t = xi[d] - xj[d];
    s += t * t;
  }
  return sqrt(s);
}


//------------------------------------------------------
// Sort a group of tied dissimilarities' pairs by their
// distances x[0..n-1] (computed from configuration `cfg`
// with k columns), permuting their codes c alongside,
// with the same result as asort4().
//
// From one iteration to the next, distances change only
// a little, so the group is nearly sorted already, and
// insertion sort is much faster than asort4(). When the
// distances are all different, any correct sort gives
// the same order. When some are exactly equal, though,
// asort4() leaves them in an order of its own; and
// insertion sort is slow when values move far. So in
// either case, the codes are restored from `bak` (n
// codes of scratch space), and asort4() is run instead.
//------------------------------------------------------
#define NM_MOVES 16 // Most moves per value before giving up.

static void sort_group(
    double *x, uint32_t *c, uint32_t *bak, size_t n,
    const double *cfg, int k ) {

  memcpy(bak, c, n * sizeof(uint32_t));

  size_t budget = NM_MOVES * n, moves = 0;
  int    ok     = 1;

  for (size_t i = 1; i < n && ok; i++) {

    double   v  = x[i];
    uint32_t cv = c[i];
    size_t   j  = i;

    while (j > 0 && x[j - 1] > v) {
      x[j] = x[j - 1];
      c[j] = c[j - 1];
      j--;
    }
    x[j] = v;
    c[j] = cv;

    moves += i - j;
    if (moves > budget) ok = 0;
  }

  for (size_t i = 1; i < n && ok; i++)
    if (x[i] == x[i - 1]) ok = 0;

  if (ok) return;

  memcpy(c, bak, n * sizeof(uint32_t));
  for (size_t i = 0; i < n; i++)
    x[i] = pair_dist(cfg, k, c[i]);
  asort4(x, c, n);
}


//------------------------------------------------------
// Visit every pair in sorted order: position `S`, code
// `C`. Pairs within a tied group come from that group's
// own order (`TB`, from offset `TO[g]`), the rest from
// the shared order.
//------------------------------------------------------
#define FOR_PAIRS(TB, TO, ...) {                                   \
  size_t S = 0;                                                   \
  for (int g_ = 0; g_ <= nm_ngrp; g_++) {                         \
    size_t gs_ = (g_ < nm_ngrp) ? (size_t) nm_gstart[g_] : nd;    \
    for (; S < gs_; S++) { uint32_t C = codes[S]; __VA_ARGS__ }          \
    if (g_ < nm_ngrp) {                                           \
      const uint32_t *t_ = (TB) + (TO)[g_] - gs_;                 \
      size_t ge_ = gs_ + (size_t) nm_glen[g_];                    \
      for (; S < ge_; S++) { uint32_t C = t_[S]; __VA_ARGS__ }           \
    }                                                             \
  }                                                               \
}


//------------------------------------------------------
// One start: monoMDS's iteration loop.
//------------------------------------------------------
static void mono_run(int start, work_t *w) {

  int    n  = nm_n;
  int    k  = nm_k;
  size_t nd = nm_ndis;
  size_t nk = (size_t) n * k;

  const uint32_t *codes = nm_codes;

  double   *x      = w->x;
  double   *grad   = w->grad;
  double   *grlast = w->grlast;
  double   *dhat   = w->dhat   - 1; // 1-based
  int      *iwork  = w->iwork  - 1; // 1-based

  // Where each tied group's pairs are listed.
  const uint32_t *tbase = codes;
  const int      *toff  = nm_gstart;
  if (nm_weak && nm_ngrp > 0) {
    for (int g = 0; g < nm_ngrp; g++)
      memcpy(
        w->tcodes + w->gtoff[g], codes + nm_gstart[g],
        (size_t) nm_glen[g] * sizeof(uint32_t) );
    tbase = w->tcodes;
    toff  = w->gtoff;
  }

  // Initial configuration, to row-major. Random ones are
  // uniform on [0, 1), as from vegan's initMDS(), drawn
  // column by column from PCG stream `nm_first + start`.
  if (nm_init) {
    const double *init = nm_init + (size_t) start * nk;
    for (int i = 0; i < n; i++)
      for (int d = 0; d < k; d++)
        x[(size_t) i * k + d] = init[(size_t) d * n + i];
  }
  else {
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, nm_seed, (uint64_t) (nm_first + start));
    for (int d = 0; d < k; d++)
      for (int i = 0; i < n; i++)
        x[(size_t) i * k + d] = pcg32_random_r(&rng) / 4294967296.0;
  }


  //--------------------------------------------------
  // Constants, with monoMDS's single-precision ones.
  //--------------------------------------------------
  double sqrtn  = (double) sqrtf((float) nd);
  double sratf1 = 0.5 * (1.0 + nm_sratmax);
  double sratf2 = 1.0 - sratf1;
  double fngrp  = 1.0;
  double strinc = (double) 1.1f;
  double cosmin = (double) -0.95f;
  double stress = 1.0;
  double cosav  = 0.0;
  double acosav = 0.0;
  double sratav = (double) 0.8f;
  double step   = 0.0;
  double fndim  = (double) k;
  double sfgr   = sqrtn;
  double fn     = (double) n;

  double strlst, sratio, sfact, tfact, dmean, cagrgl, sfglst;
  double ssfact, ssfacb = 1.0, strs;
  int    iters = 0, icause = 0;

  double g0 = sqrt(1.0 / fndim);
  for (size_t i = 0; i < nk; i++) grad[i] = g0;


  for (int it = 0; it <= nm_maxit; it++) {

    int nback = 0;

    again:
    strlst = stress;
    stress = 0.0;


    //------------------------------------------------
    // NRMCON: center, and scale to unit RMS distance
    // from the centroid.
    //------------------------------------------------
    ssfact = 0.0;
    for (int d = 0; d < k; d++) {
      double fmean = 0.0;
      for (int i = 0; i < n; i++) fmean += x[(size_t) i * k + d];
      fmean = fmean / fn;
      for (int i = 0; i < n; i++) {
        double *xi = x + (size_t) i * k + d;
        *xi    = *xi - fmean;
        ssfact = ssfact + (*xi) * (*xi);
      }
    }
    ssfact = sqrt(fn / ssfact);
    for (size_t i = 0; i < nk; i++) x[i] = x[i] * ssfact;

    for (size_t i = 0; i < nk; i++) grlast[i] = grad[i] * ssfact;
    sfglst = sfgr * ssfact;
    memset(grad, 0, nk * sizeof(double));


    //------------------------------------------------
    // MONREG, initial partition: fitted values start
    // as the distances; each untied pair is a block.
    // Weak ties: a tied group's pairs are sorted by
    // distance, each its own block. Strong ties: a
    // tied group is one block.
    //------------------------------------------------
    {
      size_t s = 0;
      for (int g = 0; g <= nm_ngrp; g++) {

        size_t gs = (g < nm_ngrp) ? (size_t) nm_gstart[g] : nd;

        for (; s < gs; s++) {
          dhat[s + 1]  = pair_dist(x, k, codes[s]);
          iwork[s + 1] = 1;
        }
        if (g == nm_ngrp) break;

        size_t len = (size_t) nm_glen[g];

        if (nm_weak) {
          uint32_t *tc = w->tcodes + w->gtoff[g];
          for (size_t t = 0; t < len; t++)
            dhat[s + t + 1] = pair_dist(x, k, tc[t]);
          sort_group(dhat + s + 1, tc, (uint32_t *) (iwork + s + 1), len, x, k);
          for (size_t t = 0; t < len; t++)
            iwork[s + t + 1] = 1;
        }
        else {
          const uint32_t *tc = codes + s;
          double sum;
          if (len == 2) {
            sum = pair_dist(x, k, tc[0]) + pair_dist(x, k, tc[1]);
          } else {
            sum = 0.0;
            for (size_t t = 0; t < len; t++) sum = sum + pair_dist(x, k, tc[t]);
          }
          dhat[s + 1]    = sum;
          dhat[s + len]  = sum;
          iwork[s + 1]   = (int) len;
          iwork[s + len] = (int) len;
        }

        s += len;
      }
    }


    //------------------------------------------------
    // MONREG, Kruskal's up-and-down blocks. Each block
    // keeps its size and sum at both of its ends.
    //------------------------------------------------
    {
      ptrdiff_t N     = (ptrdiff_t) nd;
      ptrdiff_t icurr = 1;

      while (icurr <= N) {

        int    iactiv = 1, nsatis = 0;
        double dhatav = dhat[icurr] / iwork[icurr];

        while (nsatis < 2) {
          if (iactiv == 0) {
            if (icurr == 1) {
              nsatis++;
            }
            else if (dhatav > dhat[icurr - 1] / iwork[icurr - 1]) {
              nsatis++;
            }
            else {
              ptrdiff_t icnew = icurr - iwork[icurr - 1];
              ptrdiff_t jcnew = icurr + iwork[icurr] - 1;
              iwork[icnew] = (int) (jcnew - icnew + 1);
              iwork[jcnew] = iwork[icnew];
              dhat[icnew]  = dhat[icurr - 1] + dhat[icurr];
              dhat[jcnew]  = dhat[icnew];
              icurr  = icnew;
              dhatav = dhat[icurr] / iwork[icurr];
              nsatis = 0;
            }
          }
          else {
            ptrdiff_t inext = icurr + iwork[icurr];
            if (inext > N) {
              nsatis++;
            }
            else if (dhatav < dhat[inext] / iwork[inext]) {
              nsatis++;
            }
            else {
              iwork[icurr] = iwork[icurr] + iwork[inext];
              ptrdiff_t jcurr = icurr + iwork[icurr] - 1;
              iwork[jcurr] = iwork[icurr];
              dhat[icurr]  = dhat[icurr] + dhat[inext];
              dhat[jcurr]  = dhat[icurr];
              dhatav = dhat[icurr] / iwork[icurr];
              nsatis = 0;
            }
          }
          iactiv = 1 - iactiv;
        }

        icurr = icurr + iwork[icurr];
      }

      // Spread each block's average over its members.
      icurr = 1;
      do {
        int sz = iwork[icurr];
        if (sz > 2) {
          double av = dhat[icurr] / sz;
          for (ptrdiff_t j = icurr; j <= icurr + sz - 1; j++) dhat[j] = av;
        }
        else if (sz == 2) {
          dhat[icurr]     = 0.5 * dhat[icurr];
          dhat[icurr + 1] = dhat[icurr];
        }
        icurr = icurr + sz;
      } while (icurr < N);
    }

    dhat++; // Back to 0-based for FOR_PAIRS.


    //------------------------------------------------
    // CLCSTR: Kruskal's stress.
    //------------------------------------------------
    double snumer = 0.0, sdenom = 0.0;
    dmean = 0.0;

    if (nm_isform >= 2) {
      FOR_PAIRS(tbase, toff, {
        dmean = dmean + pair_dist(x, k, C);
      })
      dmean = dmean / (double) nd;
      FOR_PAIRS(tbase, toff, {
        double dist = pair_dist(x, k, C);
        double r = dist - dhat[S], t = dist - dmean;
        snumer = snumer + r * r;
        sdenom = sdenom + t * t;
      })
    }
    else {
      FOR_PAIRS(tbase, toff, {
        double dist = pair_dist(x, k, C);
        double r = dist - dhat[S];
        snumer = snumer + r * r;
        sdenom = sdenom + dist * dist;
      })
    }
    sfact  = snumer;
    tfact  = sdenom;
    strs   = sqrt(snumer / sdenom);
    stress = stress + sfact / tfact;


    //------------------------------------------------
    // CLCGRD: the negative gradient.
    //------------------------------------------------
    if (strs > 0.0) {

      double sotsq  = sfact / (tfact * tfact);
      double recipt = 1.0 / tfact;

      FOR_PAIRS(tbase, toff, {
        double dist = pair_dist(x, k, C);
        if (dist > 0.0) {
          double wgt = (nm_isform <= 1)
            ? (sotsq - recipt * (dist - dhat[S]) / dist)
            : (sotsq * (dist - dmean) / dist - recipt * (dist - dhat[S]) / dist);
          double *xi = x    + (size_t) (C >> 16)     * k;
          double *xj = x    + (size_t) (C & 0xFFFFu) * k;
          double *gi = grad + (size_t) (C >> 16)     * k;
          double *gj = grad + (size_t) (C & 0xFFFFu) * k;
          for (int d = 0; d < k; d++) {
            double delta = wgt * (xi[d] - xj[d]);
            gi[d] = gi[d] + delta;
            gj[d] = gj[d] - delta;
          }
        }
      })
    }

    dhat--; // 1-based again.

    stress = sqrt(stress / fngrp);
    sratio = (it == 0) ? 0.8 : stress / strlst;


    //------------------------------------------------
    // CLCSFA: the gradient's scale factor, and its
    // direction cosine with the last gradient.
    //------------------------------------------------
    sfgr   = 0.0;
    cagrgl = 0.0;
    for (int d = 0; d < k; d++)
      for (int i = 0; i < n; i++) {
        double g = grad[(size_t) i * k + d];
        sfgr   = sfgr   + g * g;
        cagrgl = cagrgl + g * grlast[(size_t) i * k + d];
      }
    sfgr = sqrt(sfgr / fn);
    {
      double factor = sfgr * sfglst * fn;
      if (factor != 0.0) cagrgl = cagrgl / factor;
    }


    //------------------------------------------------
    // BACKUP: if stress rose appreciably, or the
    // gradient nearly reversed, back up along the last
    // gradient with a smaller step (3 times at most).
    //------------------------------------------------
    if ((sratio > strinc || cagrgl < cosmin) && nback < 3) {
      double stpnew = step * 0.1;
      nback++;
      if (nback == 1) ssfacb = 1.0;
      else            ssfacb = ssfacb * ssfact;
      double factor = ssfacb * (step - stpnew) / sfglst;
      for (size_t i = 0; i < nk; i++) {
        x[i]    = x[i] - factor * grlast[i];
        grad[i] = grlast[i];
      }
      step   = stpnew;
      sfgr   = sfglst;
      stress = strlst;
      goto again;
    }

    sratav = pow(sratio, 0.33334) * pow(sratav, 0.66666);
    cosav  = cagrgl * 0.66666 + cosav * 0.33334;
    acosav = fabs(cagrgl) * 0.66666 + acosav * 0.33334;


    //------------------------------------------------
    // CLCSTP: Kruskal's step size.
    //------------------------------------------------
    if (it == 0) {
      step = 25.0 * stress * sfgr;
    }
    else {
      double sa     = (sratav < 1.0) ? sratav : 1.0;
      double sa2    = sa * sa;
      double factr1 = pow(4.0, cosav);
      double factr2 = (double) 1.6f / ((1.0 + sa2 * (sa2 * sa)) * (1.0 + acosav - fabs(cosav)));
      double factr3 = sqrt((sratio < 1.0) ? sratio : 1.0);
      step = step * factr1 * factr2 * factr3;
    }


    //------------------------------------------------
    // Stop?
    //------------------------------------------------
    iters = it;
    if (stress < nm_smin) {
      icause = 2;
      break;
    }
    else if (fabs(sratav - sratf1) <= sratf2 && fabs(sratio - sratf1) <= sratf2) {
      icause = 3;
      break;
    }
    else if (sfgr <= nm_sfgrmin) {
      icause = 4;
      break;
    }
    else if (it == nm_maxit) {
      icause = 1;
      break;
    }


    //------------------------------------------------
    // NEWCON: step along the negative gradient.
    //------------------------------------------------
    {
      double factor = step / sfgr;
      for (size_t i = 0; i < nk; i++)
        x[i] = x[i] + factor * grad[i];
    }
  }


  //--------------------------------------------------
  // Results, column-major.
  //--------------------------------------------------
  double *out = nm_points + (size_t) start * nk;
  for (int i = 0; i < n; i++)
    for (int d = 0; d < k; d++)
      out[(size_t) d * n + i] = x[(size_t) i * k + d];

  nm_stress[start] = stress;
  nm_iters[start]  = iters;
  nm_icause[start] = icause;
}


//------------------------------------------------------
// Worker: run starts until none are left.
//------------------------------------------------------
static void *mono_worker(void *arg) {

  work_t *w = nm_work + ((worker_t *) arg)->i;

  while (1) {

    int start;

    #ifdef HAVE_PTHREAD
      pthread_mutex_lock(&nm_mutex);
    #endif
    start = nm_next++;
    #ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&nm_mutex);
    #endif

    if (start >= nm_m) break;

    mono_run(start, w);
  }

  return NULL;
}


//------------------------------------------------------
// Groups of tied values in sorted x[0..n-1], as monoMDS
// finds them: runs whose successive differences are at
// most sqrt(DBL_EPSILON). (monoMDS fails to end the last
// run when the largest value is that small; here the last
// value always ends one.) Stores each group of two or
// more's first position and length, unless gstart is
// NULL, and returns how many there are.
//------------------------------------------------------
static int tied_groups(const double *x, size_t n, int *gstart, int *glen) {

  double toler = sqrt(DBL_EPSILON);
  int    ngrp  = 0;
  size_t j     = 0;

  for (size_t i = 1; i <= n; i++)
    if (i == n || fabs(x[i] - x[i - 1]) > toler) {
      if (i - j > 1) {
        if (gstart) {
          gstart[ngrp] = (int) j;
          glen[ngrp]   = (int) (i - j);
        }
        ngrp++;
      }
      j = i;
    }

  return ngrp;
}


//======================================================
// R interface. Sorts the pairs of a dist object's
// packed lower triangle (`sexp_dist`, `sexp_n` samples)
// by dissimilarity, as monoMDS does, and finds the
// groups of tied dissimilarities. Returns list(codes,
// gstart, glen): the sorted pairs as i << 16 | j
// (0-based, i > j), and each tied group's 0-based first
// position and length.
//======================================================
SEXP C_nmds_prep(SEXP sexp_dist, SEXP sexp_n, SEXP sexp_stable) {

  int    n  = asInteger(sexp_n);
  int    stable = asLogical(sexp_stable) == 1;
  size_t nd = (size_t) n * (n - 1) / 2;

  if (n < 2 || n > 65536 || (size_t) XLENGTH(sexp_dist) != nd)
    error("Invalid arguments to C_nmds_prep()."); // # nocov

  SEXP sexp_codes = PROTECT(allocVector(INTSXP, (R_xlen_t) nd));
  uint32_t *codes = (uint32_t *) INTEGER(sexp_codes);
  double   *diss  = (double *) R_alloc(nd, sizeof(double));

  memcpy(diss, REAL(sexp_dist), nd * sizeof(double));

  // dist order: column j, then rows i = j + 1, ..., n - 1.
  size_t p = 0;
  for (uint32_t j = 0; j < (uint32_t) n; j++)
    for (uint32_t i = j + 1; i < (uint32_t) n; i++)
      codes[p++] = (i << 16) | j;

  // With weak ties, each start re-sorts every group of tied
  // dissimilarities by its own fitted distances, so the
  // order within a group doesn't matter, and a stable
  // O(n log n) sort gives the same results as ASORT4. With
  // strong ties, ASORT4's order of tied pairs is used.
  if (stable) radix_sort(diss, codes, nd);
  else        asort4(diss, codes, nd);


  int  ngrp        = tied_groups(diss, nd, NULL, NULL);
  SEXP sexp_gstart = PROTECT(allocVector(INTSXP, ngrp));
  SEXP sexp_glen   = PROTECT(allocVector(INTSXP, ngrp));
  tied_groups(diss, nd, INTEGER(sexp_gstart), INTEGER(sexp_glen));

  const char *names[] = { "codes", "gstart", "glen", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_codes);
  SET_VECTOR_ELT(sexp_result, 1, sexp_gstart);
  SET_VECTOR_ELT(sexp_result, 2, sexp_glen);

  UNPROTECT(4);
  return sexp_result;
}


//======================================================
// R interface. Runs monoMDS from m starting
// configurations, given the sorted pairs from
// C_nmds_prep(). `sexp_init` is either an n x k x m array
// of them, or NULL for random ones. `sexp_opts` is
// c(maxit, weakties, isform, n_threads, seed, first, m),
// where seed, first and m are for random starts (start
// t uses PCG stream first + t), and `sexp_tols` is
// c(smin, sfgrmin, sratmax). Returns list(points, stress,
// iters, icause), with points n x k x m.
//======================================================
SEXP C_nmds_run(
    SEXP sexp_prep, SEXP sexp_n, SEXP sexp_k, SEXP sexp_init,
    SEXP sexp_opts, SEXP sexp_tols ) {

  int n         = asInteger(sexp_n);
  int k         = asInteger(sexp_k);
  int maxit     = INTEGER(sexp_opts)[0];
  int weak      = INTEGER(sexp_opts)[1];
  int isform    = INTEGER(sexp_opts)[2];
  int n_threads = INTEGER(sexp_opts)[3];
  int seed      = INTEGER(sexp_opts)[4];
  int first     = INTEGER(sexp_opts)[5];
  int m         = INTEGER(sexp_opts)[6];

  size_t nd = (size_t) n * (n - 1) / 2;
  size_t nk = (size_t) n * k;

  SEXP sexp_codes  = VECTOR_ELT(sexp_prep, 0);
  SEXP sexp_gstart = VECTOR_ELT(sexp_prep, 1);
  SEXP sexp_glen   = VECTOR_ELT(sexp_prep, 2);

  if (!isNull(sexp_init)) m = (int) (XLENGTH(sexp_init) / nk);

  if (n < 2 || n > 65536 || k < 1 || n_threads < 1 || m < 1 ||
      (size_t) XLENGTH(sexp_codes) != nd ||
      (!isNull(sexp_init) && (size_t) XLENGTH(sexp_init) != nk * m))
    error("Invalid arguments to C_nmds_run()."); // # nocov

  nm_n       = n;
  nm_k       = k;
  nm_ndis    = nd;
  nm_codes   = (const uint32_t *) INTEGER(sexp_codes);
  nm_ngrp    = (int) XLENGTH(sexp_gstart);
  nm_gstart  = INTEGER(sexp_gstart);
  nm_glen    = INTEGER(sexp_glen);
  nm_weak    = weak;
  nm_isform  = isform;
  nm_maxit   = maxit;
  nm_smin    = REAL(sexp_tols)[0];
  nm_sfgrmin = REAL(sexp_tols)[1];
  nm_sratmax = REAL(sexp_tols)[2];
  nm_m       = m;
  nm_init    = isNull(sexp_init) ? NULL : REAL(sexp_init);
  nm_seed    = (uint64_t) seed;
  nm_first   = first;
  nm_next    = 0;

  nm_ntied = 0;
  for (int g = 0; g < nm_ngrp; g++) nm_ntied += (size_t) nm_glen[g];

  SEXP sexp_points = PROTECT(allocVector(REALSXP, (R_xlen_t) (nk * m)));
  SEXP sexp_stress = PROTECT(allocVector(REALSXP, m));
  SEXP sexp_iters  = PROTECT(allocVector(INTSXP,  m));
  SEXP sexp_icause = PROTECT(allocVector(INTSXP,  m));
  nm_points = REAL(sexp_points);
  nm_stress = REAL(sexp_stress);
  nm_iters  = INTEGER(sexp_iters);
  nm_icause = INTEGER(sexp_icause);


  //--------------------------------------------------
  // Work space for each start that runs at once.
  //--------------------------------------------------
  if (n_threads > m) n_threads = m;

  nm_work = (work_t *) R_alloc(n_threads, sizeof(work_t));
  for (int t = 0; t < n_threads; t++) {
    work_t *w = nm_work + t;
    w->x      = (double *) R_alloc(nk, sizeof(double));
    w->grad   = (double *) R_alloc(nk, sizeof(double));
    w->grlast = (double *) R_alloc(nk, sizeof(double));
    w->dhat   = (double *) R_alloc(nd, sizeof(double));
    w->iwork  = (int *)    R_alloc(nd, sizeof(int));
    w->tcodes = NULL;
    w->gtoff  = NULL;
    if (weak && nm_ngrp > 0) {
      w->tcodes = (uint32_t *) R_alloc(nm_ntied, sizeof(uint32_t));
      w->gtoff  = (int *)      R_alloc(nm_ngrp,  sizeof(int));
      int off = 0;
      for (int g = 0; g < nm_ngrp; g++) {
        w->gtoff[g] = off;
        off += nm_glen[g];
      }
    }
  }

  // run_parallel() only starts threads when n_tasks >= 100,
  // which keeps tiny problems on one thread.
  run_parallel(mono_worker, n_threads, (nd >= 100) ? 100 : 1);

  const char *names[] = { "points", "stress", "iters", "icause", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, sexp_stress);
  SET_VECTOR_ELT(sexp_result, 2, sexp_iters);
  SET_VECTOR_ELT(sexp_result, 3, sexp_icause);

  UNPROTECT(5);
  return sexp_result;
}



//------------------------------------------------------
// Post-processing: sums over pairs for the final scaling,
// as vegan's postMDS(). Each task handles one column of
// the dist vector; the columns' results are combined in
// order on the main thread.
//------------------------------------------------------
static const double *sc_dist;   // Packed lower triangle.
static const double *sc_x;      // n x k, row-major.
static double        sc_thresh; // Regress pairs with dist below.
static double        sc_mx;     // Mean configuration distance ...
static double        sc_my;     // ... and dissimilarity, of those.
static double       *sc_out;    // n x 6, row-major.

// Column j's elements, rows j + 1 to n - 1, as d[j + 1] to d[n - 1].
#define COL_PTR(j, n) (sc_dist + (size_t) (j) * (n) - (size_t) (j) * ((j) + 1) / 2 - ((j) + 1))

static inline double sc_pair(int i, int j, int k) {
  const double *xi = sc_x + (size_t) i * k;
  const double *xj = sc_x + (size_t) j * k;
  double s = 0;
  for (int a = 0; a < k; a++) { double t = xi[a] - xj[a]; s += t * t; }
  return sqrt(s);
}

static void *scale_pass1(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = nm_n, k = nm_k;

  for (int j = thread_i; j < n - 1; j += n_threads) {

    const double *d = COL_PTR(j, n);
    double mx = 0, cnt = 0, sx = 0, sy = 0;

    for (int i = j + 1; i < n; i++) {
      double o = sc_pair(i, j, k);
      if (o > mx) mx = o;
      if (d[i] < sc_thresh) { cnt += 1; sx += o; sy += d[i]; }
    }

    double *out = sc_out + (size_t) j * 6;
    out[0] = mx; out[1] = cnt; out[2] = sx; out[3] = sy;
  }

  return NULL;
}

static void *scale_pass2(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = nm_n, k = nm_k;

  for (int j = thread_i; j < n - 1; j += n_threads) {

    const double *d = COL_PTR(j, n);
    double sxx = 0, sxy = 0;

    for (int i = j + 1; i < n; i++)
      if (d[i] < sc_thresh) {
        double ox = sc_pair(i, j, k) - sc_mx;
        sxx += ox * ox;
        sxy += ox * (d[i] - sc_my);
      }

    double *out = sc_out + (size_t) j * 6;
    out[4] = sxx; out[5] = sxy;
  }

  return NULL;
}


//======================================================
// R interface. For the configuration `sexp_points`
// (n x k) and dist vector `sexp_dist`, returns c(max,
// count, intercept, slope): the largest distance between
// points, and the least-squares line of dissimilarity on
// configuration distance over the `count` pairs with
// dissimilarity below `sexp_thresh` (NA if count < 2).
//======================================================
SEXP C_nmds_scale(SEXP sexp_dist, SEXP sexp_points, SEXP sexp_thresh, SEXP sexp_n_threads) {

  int n         = nrows(sexp_points);
  int k         = ncols(sexp_points);
  int n_threads = asInteger(sexp_n_threads);

  if (n < 2 || (size_t) XLENGTH(sexp_dist) != (size_t) n * (n - 1) / 2)
    error("Invalid arguments to C_nmds_scale()."); // # nocov

  double *x = (double *) R_alloc((size_t) n * k, sizeof(double));
  for (int i = 0; i < n; i++)
    for (int a = 0; a < k; a++)
      x[(size_t) i * k + a] = REAL(sexp_points)[(size_t) a * n + i];

  nm_n      = n;
  nm_k      = k;
  sc_dist   = REAL(sexp_dist);
  sc_x      = x;
  sc_thresh = asReal(sexp_thresh);
  sc_out    = (double *) R_alloc((size_t) n * 6, sizeof(double));

  run_parallel(scale_pass1, n_threads, n);

  double mx = 0, cnt = 0, sx = 0, sy = 0;
  for (int j = 0; j < n - 1; j++) {
    double *o = sc_out + (size_t) j * 6;
    if (o[0] > mx) mx = o[0];
    cnt += o[1]; sx += o[2]; sy += o[3];
  }

  double b0 = NA_REAL, b1 = NA_REAL;

  if (cnt >= 2) {

    sc_mx = sx / cnt;
    sc_my = sy / cnt;

    run_parallel(scale_pass2, n_threads, n);

    double sxx = 0, sxy = 0;
    for (int j = 0; j < n - 1; j++) {
      sxx += sc_out[(size_t) j * 6 + 4];
      sxy += sc_out[(size_t) j * 6 + 5];
    }

    b1 = sxy / sxx;
    b0 = sc_my - b1 * sc_mx;
  }

  SEXP sexp_result = PROTECT(allocVector(REALSXP, 4));
  REAL(sexp_result)[0] = mx;
  REAL(sexp_result)[1] = cnt;
  REAL(sexp_result)[2] = b0;
  REAL(sexp_result)[3] = b1;

  UNPROTECT(1);
  return sexp_result;
}
