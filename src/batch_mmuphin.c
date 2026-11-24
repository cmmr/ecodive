// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// batch_mmuphin() (R/batch_mmuphin.r): a C port of adjust_batch() and its
// helpers (construct_ind(), fit_stand_feature(), standardize_feature(),
// fit_EB(), aprior(), bprior(), fit_shrink(), it_sol(), postmean(),
// postvar(), relocate_scale(), add_back_covariates(), back_transform_abd(),
// normalize_features(), transform_features(), TSS(), LOG()) from MMUPHin
// 2.0.0 (R/adjust_batch.R, R/helpers_adjust_batch.R, R/helpers.R; these
// are identical in MMUPHin 1.16.0), restructured for ecodive, see
// plans/batch_correction.md.
//
// MMUPHin's notice:
//
//   MMUPHin: Copyright (c) 2019 Harvard School of Public Health.
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//
//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// The steps, with every sum in a fixed order (so results are the same for
// any number of threads):
//
//   1. Per feature (threads): its raw values, which are zero, and its
//      relative abundances' means per batch.
//   2. Per sample (threads): the total, and the total with the pseudocount
//      added to every feature (TSS of counts + pseudocount).
//   3. Per feature (threads): log2 of TSS; with zero inflation, which
//      batches the feature is adjusted in (its non-zero samples' batches,
//      if at least two, with a full-rank design that has more rows than
//      columns: construct_ind()); the least-squares fit on that design
//      (feature_lm.c's pivoted QR, where MMUPHin solves the normal
//      equations), the standardisation (standardize_feature()), and each
//      batch's mean and standard deviation of the standardised values
//      (fit_EB()).
//   4. Per batch (main thread): the priors over features, in feature order.
//   5. Per batch (threads): the shrinkage iteration (it_sol()), over the
//      features adjusted in the batch, in feature order.
//   6. Per feature (threads): the adjustment (relocate_scale(),
//      add_back_covariates()), 2^x, and zeros back to zero.
//   7. Per sample (threads): the adjusted total; per feature (threads):
//      re-closure to each sample's original total, rounding, and the
//      adjusted relative abundances' means per batch.

#include "ecodive.h"
#include <float.h>   // DBL_EPSILON


static int           bm_n, bm_D, bm_B, bm_q, bm_p, bm_zi, bm_round;
static double        bm_pc, bm_tol, bm_qr_tol;
static int           bm_max_iters;
static const da_table_t *bm_tab;
static const int    *bm_batch;    // n: each sample's batch, 0-based
static const double *bm_X;        // n x p: batch indicators, then covariates
static double       *bm_out;      // n x D: values, then the result
static char         *bm_nz;       // n x D: 1 where the raw value isn't 0
static double       *bm_tot;      // n: sample totals
static double       *bm_tot_pc;   // n: totals with the pseudocount
static double       *bm_adj_tot;  // n: adjusted totals (before re-closure)
static double       *bm_mean0;    // D x B: mean relative abundance, before
static double       *bm_mean1;    // D x B: ... and after
static double       *bm_overall;  // D: mean relative abundance, all samples
static int          *bm_ind;      // D x B: adjusted in batch b
static double       *bm_ghat, *bm_dhat, *bm_gstar, *bm_dstar;   // D x B
static double       *bm_gm;       // D: grand mean
static double       *bm_vp;       // D: pooled variance
static double       *bm_bcov;     // D x q: covariate coefficients
static int          *bm_n_fb;     // D x B: values of feature f in batch b
static double       *bm_gbar, *bm_t2, *bm_a, *bm_b;   // B: priors
static int          *bm_iters;    // B
static int          *bm_status;   // B: 0 ok, 1 too many iterations
static double       *bm_mem;      // per thread: flm workspace
static int          *bm_imem;     // per thread: p + n ints
static size_t        bm_mem_sz;
static double       *bm_buf;      // per thread: scratch
static size_t        bm_buf_sz;
static const flm_t  *bm_full;     // QR of the whole design (no zero inflation)
static const int    *bm_bsam;     // samples ordered by batch ...
static const int    *bm_bpos;     // ... batch b's are bsam[bpos[b] .. bpos[b + 1] - 1]


// R's mean(): the long double mean, refined by the mean of the residuals.
static double r_mean (const double *x, int n) {
  if (n < 1) return R_NaN;
  long double s = 0;
  for (int i = 0; i < n; i++) s += x[i];
  s /= n;
  if (R_FINITE((double) s)) {
    long double t = 0;
    for (int i = 0; i < n; i++) t += (x[i] - s);
    s += t / n;
  }
  return (double) s;
}


// R's var(), with its own (long double) mean.
static double r_var (const double *x, int n) {
  if (n < 2) return R_NaN;
  long double sum = 0;
  for (int i = 0; i < n; i++) sum += x[i];
  long double tmp = sum / n;
  if (R_FINITE((double) tmp)) {
    sum = 0;
    for (int i = 0; i < n; i++) sum += (x[i] - tmp);
    tmp = tmp + sum / n;
  }
  double xm = (double) tmp;
  sum = 0;
  for (int i = 0; i < n; i++) sum += (x[i] - xm) * (x[i] - xm);
  return (double) (sum / (n - 1));
}



//======================================================
// 1. Raw values and which are zero.
//======================================================
static void *bm_raw_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D, B = bm_B;
  double *ra = bm_buf + (size_t) thread_i * bm_buf_sz;

  for (int f = thread_i; f < D; f += n_threads) {
    double *y  = bm_out + (size_t) n * f;
    char   *nz = bm_nz  + (size_t) n * f;
    da_table_column(bm_tab, f, y);
    for (int s = 0; s < n; s++) nz[s] = (y[s] != 0);

    // Relative abundances (TSS: an all-zero sample stays zero).
    for (int s = 0; s < n; s++) ra[s] = bm_tot[s] > 0 ? y[s] / bm_tot[s] : 0;
    for (int b = 0; b < B; b++) {
      long double m = 0;
      int lo = bm_bpos[b], hi = bm_bpos[b + 1];
      for (int k = lo; k < hi; k++) m += ra[bm_bsam[k]];
      bm_mean0[f + (size_t) D * b] = (double) (m / (hi - lo));
    }
    long double m = 0;
    for (int s = 0; s < n; s++) m += ra[s];
    bm_overall[f] = (double) (m / n);
  }
  return NULL;
}



//======================================================
// 2. Sample totals, in feature order.
//======================================================
static void *bm_tot_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D;

  for (int s = thread_i; s < n; s += n_threads) {
    long double tp = 0;
    for (int f = 0; f < D; f++) tp += bm_out[s + (size_t) n * f] + bm_pc;
    bm_tot_pc[s] = (double) tp;
  }
  return NULL;
}


// Raw totals (apply(feature_abd, 2, sum)), straight from the table.
static void bm_raw_totals (void) {
  int n = bm_n, D = bm_D;
  long double *acc = (long double *) R_alloc((size_t) n + 1, sizeof(long double));
  for (int s = 0; s < n; s++) acc[s] = 0;
  for (int f = 0; f < D; f++)
    for (int e = bm_tab->fpos[f]; e < bm_tab->fpos[f + 1]; e++)
      acc[bm_tab->fsam[e]] += bm_tab->fval[e];
  for (int s = 0; s < n; s++) bm_tot[s] = (double) acc[s];
}



//======================================================
// A row's covariate part: sum_j X[i, k + j] beta[j],
// as R's matrix product sums it (in column order).
//======================================================
static inline double bm_cov_part (int i, const double *beta) {
  double c = 0;
  for (int j = 0; j < bm_q; j++) c += bm_X[i + (size_t) bm_n * (bm_B + j)] * beta[j];
  return c;
}



//======================================================
// 3. Per feature: which batches, the fit, the
//    standardised values, and the batch moments.
//======================================================
static void *bm_fit_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D, B = bm_B, q = bm_q, p = bm_p;

  double *buf   = bm_buf + (size_t) thread_i * bm_buf_sz;
  double *yo    = buf;                       // n
  double *qty   = yo   + n;                  // n
  double *fit   = qty  + n;                  // n
  double *res   = fit  + n;                  // n
  double *coef  = res  + n;                  // p
  double *Xf    = coef + p;                  // n x p
  int    *imem  = bm_imem + (size_t) thread_i * (2 * (size_t) p + (size_t) n + B + 2);
  int    *rows  = imem + p;                  // n
  int    *cols  = rows + n;                  // p
  int    *pres  = cols + p;                  // B

  flm_t mq;
  flm_init(&mq, bm_mem + (size_t) thread_i * bm_mem_sz, imem, n, p, 0);

  for (int f = thread_i; f < D; f += n_threads) {

    double     *y  = bm_out + (size_t) n * f;
    const char *nz = bm_nz  + (size_t) n * f;
    int        *ind = bm_ind + f;   // stride D

    // log2 of the TSS of counts + pseudocount.
    for (int s = 0; s < n; s++) y[s] = log2((y[s] + bm_pc) / bm_tot_pc[s]);

    // construct_ind(): the samples and batches that enter.
    int nr = 0;
    for (int b = 0; b < B; b++) pres[b] = 0;
    for (int s = 0; s < n; s++)
      if (!bm_zi || nz[s]) { rows[nr++] = s; pres[bm_batch[s]] = 1; }

    int k = 0;
    for (int b = 0; b < B; b++) if (pres[b]) k++;

    int ok = 1;
    const flm_t *qr = NULL;
    const double *Xp = Xf;   // the feature's design, nr x pf
    flm_t qq;
    int pf = k + q;

    if (bm_zi) {
      ok = (k > 1) && (nr > pf);
      if (ok) {
        int c = 0;
        for (int b = 0; b < B; b++) if (pres[b]) cols[c++] = b;
        for (int j = 0; j < q; j++) cols[c++] = B + j;
        for (int j = 0; j < pf; j++) {
          const double *x = bm_X + (size_t) n * cols[j];
          double *xf = Xf + (size_t) nr * j;
          for (int i = 0; i < nr; i++) xf[i] = x[rows[i]];
        }
        flm_decompose(&mq, Xf, nr, NULL, nr, pf, bm_qr_tol);
        ok = (mq.rank == pf);
        qq = mq;
        qr = &qq;
      }
    }
    else {
      Xp      = bm_X;
      qq      = *bm_full;
      qq.work = mq.work;
      qr      = &qq;
    }

    for (int b = 0; b < B; b++) {
      ind[(size_t) D * b] = ok && pres[b];
      bm_ghat[f + (size_t) D * b] = R_NaN;
      bm_dhat[f + (size_t) D * b] = R_NaN;
      bm_n_fb[f + (size_t) D * b] = 0;
    }
    for (int j = 0; j < q; j++) bm_bcov[f + (size_t) D * j] = R_NaN;
    bm_gm[f] = bm_vp[f] = R_NaN;
    if (!ok) continue;

    // standardize_feature(): least squares on the feature's design.
    for (int i = 0; i < nr; i++) yo[i] = y[rows[i]];
    double rss;
    flm_fit(qr, yo, qty, coef, &rss);

    // Fitted values, and the batch part's mean (the grand mean).
    for (int i = 0; i < nr; i++) {
      double fv = 0, bv = 0;
      for (int j = 0; j < pf; j++) {
        double t = Xp[i + (size_t) nr * j] * coef[j];
        fv += t;
        if (j < k) bv += t;
      }
      fit[i] = fv;
      res[i] = bv;
    }
    double gm = r_mean(res, nr);
    for (int i = 0; i < nr; i++) res[i] = yo[i] - fit[i];
    double vp = r_var(res, nr);
    // isTRUE(all.equal(var_pooled, 0)): at most sqrt(.Machine$double.eps).
    if (!(vp > sqrt(DBL_EPSILON))) vp = 1;

    double *bc = coef + k;
    for (int j = 0; j < q; j++) bm_bcov[f + (size_t) D * j] = bc[j];
    bm_gm[f] = gm;
    bm_vp[f] = vp;

    double sd = sqrt(vp);
    for (int i = 0; i < nr; i++) {
      int s = rows[i];
      double sm = gm;
      if (q > 0) sm = gm + bm_cov_part(s, bc);
      y[s] = (yo[i] - sm) / sd;
    }

    // fit_EB(): each batch's mean and standard deviation.
    for (int b = 0; b < B; b++) {
      if (!pres[b]) continue;
      int m = 0;
      for (int i = 0; i < nr; i++) if (bm_batch[rows[i]] == b) fit[m++] = y[rows[i]];
      double g = r_mean(fit, m);
      double d = sqrt(r_var(fit, m));
      if (ISNAN(d) || d == 0) d = 1;
      bm_ghat[f + (size_t) D * b] = g;
      bm_dhat[f + (size_t) D * b] = d;
      bm_n_fb[f + (size_t) D * b] = m;
    }
  }

  return NULL;
}



//======================================================
// 5. it_sol(), for one batch per task.
//======================================================
static void *bm_shrink_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D, B = bm_B;

  for (int b = thread_i; b < B; b += n_threads) {

    double *gs = bm_gstar + (size_t) D * b;
    double *ds = bm_dstar + (size_t) D * b;
    const double *gh = bm_ghat + (size_t) D * b;
    const double *dh = bm_dhat + (size_t) D * b;
    const int    *in = bm_ind  + (size_t) D * b;
    const int    *nb = bm_n_fb + (size_t) D * b;
    double gbar = bm_gbar[b], t2 = bm_t2[b], a = bm_a[b], bb = bm_b[b];
    int lo = bm_bpos[b], hi = bm_bpos[b + 1];

    for (int f = 0; f < D; f++) {
      gs[f] = in[f] ? gh[f] : R_NaN;
      ds[f] = in[f] ? dh[f] : R_NaN;
    }

    double change = 1;
    int    count  = 0;
    bm_status[b]  = 0;

    while (change > bm_tol) {

      change = R_NegInf;

      for (int f = 0; f < D; f++) {
        if (!in[f]) continue;
        double g_old = gs[f], d_old = ds[f];
        double nn    = nb[f];

        // postmean(), then postvar() of the sum of squares.
        double g_new = (t2 * nn * gh[f] + d_old * gbar) / (t2 * nn + d_old);

        const double *y  = bm_out + (size_t) n * f;
        const char   *nz = bm_nz  + (size_t) n * f;
        long double sum2 = 0;
        for (int kk = lo; kk < hi; kk++) {
          int s = bm_bsam[kk];
          if (bm_zi && !nz[s]) continue;
          double r = y[s] - g_new;
          sum2 += r * r;
        }
        double s2    = (double) sum2;
        double d_new = (s2 == 0) ? R_NaN : (0.5 * s2 + bb) / (nn / 2 + a - 1);

        // max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old,
        // na.rm = TRUE): the divisors are signed, as in MMUPHin.
        double r1 = fabs(g_new - g_old) / g_old;
        double r2 = fabs(d_new - d_old) / d_old;
        if (!ISNAN(r1) && r1 > change) change = r1;
        if (!ISNAN(r2) && r2 > change) change = r2;

        gs[f] = g_new;
        ds[f] = d_new;
      }

      count++;
      if (count > bm_max_iters) { bm_status[b] = 1; break; }
    }

    bm_iters[b] = count;
  }

  return NULL;
}



//======================================================
// 6. Adjust, and back to the abundance scale.
//======================================================
static void *bm_adjust_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D, B = bm_B;

  for (int f = thread_i; f < D; f += n_threads) {

    double     *y  = bm_out + (size_t) n * f;
    const char *nz = bm_nz  + (size_t) n * f;

    int adjusted = 0;
    for (int b = 0; b < B; b++) if (bm_ind[f + (size_t) D * b]) adjusted = 1;

    if (adjusted) {
      double sd = sqrt(bm_vp[f]);
      double gm = bm_gm[f];
      double *bcp = bm_buf + (size_t) thread_i * bm_buf_sz;   // q
      for (int j = 0; j < bm_q; j++) bcp[j] = bm_bcov[f + (size_t) D * j];

      for (int s = 0; s < n; s++) {
        if (bm_zi && !nz[s]) continue;
        int b = bm_batch[s];
        double v = y[s];
        if (bm_ind[f + (size_t) D * b])
          v = (v - bm_gstar[f + (size_t) D * b]) / sqrt(bm_dstar[f + (size_t) D * b]);
        double sm = gm;
        if (bm_q > 0) sm = gm + bm_cov_part(s, bcp);
        y[s] = v * sd + sm;
      }
    }

    // back_transform_abd(): 2^x, and the original zeros back to 0.
    for (int s = 0; s < n; s++) y[s] = nz[s] ? pow(2.0, y[s]) : 0;
  }

  return NULL;
}



//======================================================
// 7. Re-closure: totals, then each value.
//======================================================
static void *bm_adj_tot_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D;

  for (int s = thread_i; s < n; s += n_threads) {
    long double t = 0;
    for (int f = 0; f < D; f++) t += bm_out[s + (size_t) n * f];
    bm_adj_tot[s] = (double) t;
  }
  return NULL;
}


static void *bm_close_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D;

  for (int f = thread_i; f < D; f += n_threads) {
    double *y = bm_out + (size_t) n * f;
    for (int s = 0; s < n; s++) {
      double u = bm_adj_tot[s];
      double v = (u == 0) ? 0 : (y[s] / u) * bm_tot[s];
      if (bm_round) v = nearbyint(v);
      y[s] = v;
    }
  }
  return NULL;
}


static void *bm_mean1_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bm_n, D = bm_D, B = bm_B;

  for (int f = thread_i; f < D; f += n_threads) {
    const double *y = bm_out + (size_t) n * f;
    for (int b = 0; b < B; b++) {
      long double m = 0;
      int lo = bm_bpos[b], hi = bm_bpos[b + 1];
      for (int k = lo; k < hi; k++) {
        int s = bm_bsam[k];
        double u = bm_adj_tot[s];
        m += u > 0 ? y[s] / u : 0;
      }
      bm_mean1[f + (size_t) D * b] = (double) (m / (hi - lo));
    }
  }
  return NULL;
}



//======================================================
// R interface, for batch_mmuphin().
//
//   counts, margin   The table.
//   batch            n: each sample's batch, 0-based.
//   X                n x (B + q): the batch indicators, then the
//                    covariates (no intercept).
//   n_batches        B.
//   zero_inflation   0 or 1.
//   pseudocount      Added to every value before TSS.
//   round            Round the result (counts in, counts out).
//   tol, max_iters   it_sol()'s conv and maxit.
//
// Returns list(adjusted (n x D), ind (D x B logical), gamma_hat,
// gamma_star, delta_hat, delta_star (D x B), n_iters (B), status
// (0: ok; 1: every feature is single-batch; 2: some batches can't be
// adjusted, see `bad`; 3: a batch's iteration didn't converge, see
// `which`; 4: a shrunken scale is NA, see `which`), bad (B logical),
// which (batch, feature, 1-based), means_before, means_after (D x B),
// overall (D)).
//======================================================
SEXP C_batch_mmuphin(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_batch, SEXP sexp_X,
    SEXP sexp_n_batches, SEXP sexp_zero_inflation, SEXP sexp_pseudocount,
    SEXP sexp_round, SEXP sexp_tol, SEXP sexp_max_iters, SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, 0, 0, n_thr, &tab);
  int n = tab.n, D = tab.D;

  int B = asInteger(sexp_n_batches);
  SEXP dims = getAttrib(sexp_X, R_DimSymbol);
  if (INTEGER(dims)[0] != n) error("The design must have one row per sample.");  // # nocov
  int p = INTEGER(dims)[1];
  int q = p - B;
  if (LENGTH(sexp_batch) != n) error("`batch` must have one value per sample.");  // # nocov

  bm_n = n; bm_D = D; bm_B = B; bm_q = q; bm_p = p;
  bm_zi        = asLogical(sexp_zero_inflation);
  bm_round     = asLogical(sexp_round);
  bm_pc        = asReal(sexp_pseudocount);
  bm_tol       = asReal(sexp_tol);
  bm_max_iters = asInteger(sexp_max_iters);
  bm_qr_tol    = 1e-7;
  bm_tab       = &tab;
  bm_batch     = INTEGER(sexp_batch);
  bm_X         = REAL(sexp_X);

  const char *nms[] = {
    "adjusted", "ind", "gamma_hat", "gamma_star", "delta_hat", "delta_star",
    "n_iters", "status", "bad", "which", "means_before", "means_after",
    "overall", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_out  = allocMatrix(REALSXP, n, D); SET_VECTOR_ELT(result, 0, s_out);
  SEXP s_ind  = allocMatrix(LGLSXP, D, B);  SET_VECTOR_ELT(result, 1, s_ind);
  SEXP s_gh   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 2, s_gh);
  SEXP s_gs   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 3, s_gs);
  SEXP s_dh   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 4, s_dh);
  SEXP s_ds   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 5, s_ds);
  SEXP s_it   = allocVector(INTSXP, B);     SET_VECTOR_ELT(result, 6, s_it);
  SEXP s_bad  = allocVector(LGLSXP, B);     SET_VECTOR_ELT(result, 8, s_bad);
  SEXP s_wh   = allocVector(INTSXP, 2);     SET_VECTOR_ELT(result, 9, s_wh);
  SEXP s_m0   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 10, s_m0);
  SEXP s_m1   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 11, s_m1);
  SEXP s_ov   = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 12, s_ov);
  int status  = 0;

  bm_out   = REAL(s_out);
  bm_ind   = LOGICAL(s_ind);
  bm_ghat  = REAL(s_gh);  bm_gstar = REAL(s_gs);
  bm_dhat  = REAL(s_dh);  bm_dstar = REAL(s_ds);
  bm_iters = INTEGER(s_it);
  bm_mean0 = REAL(s_m0);  bm_mean1 = REAL(s_m1);
  bm_overall = REAL(s_ov);
  for (int b = 0; b < B; b++) { bm_iters[b] = 0; LOGICAL(s_bad)[b] = 0; }
  INTEGER(s_wh)[0] = INTEGER(s_wh)[1] = NA_INTEGER;
  for (size_t i = 0; i < (size_t) D * B; i++) bm_gstar[i] = bm_dstar[i] = bm_mean1[i] = R_NaN;

  // Samples by batch.
  int *bpos = (int *) R_alloc((size_t) B + 1, sizeof(int));
  int *bsam = (int *) R_alloc((size_t) n + 1, sizeof(int));
  memset(bpos, 0, ((size_t) B + 1) * sizeof(int));
  for (int s = 0; s < n; s++) bpos[bm_batch[s] + 1]++;
  for (int b = 0; b < B; b++) bpos[b + 1] += bpos[b];
  int *fill = (int *) R_alloc((size_t) B + 1, sizeof(int));
  memcpy(fill, bpos, ((size_t) B + 1) * sizeof(int));
  for (int s = 0; s < n; s++) bsam[fill[bm_batch[s]]++] = s;
  bm_bpos = bpos;
  bm_bsam = bsam;

  bm_nz      = (char *)   R_alloc((size_t) n * D + 1, sizeof(char));
  bm_tot     = (double *) R_alloc((size_t) n + 1, sizeof(double));
  bm_tot_pc  = (double *) R_alloc((size_t) n + 1, sizeof(double));
  bm_adj_tot = (double *) R_alloc((size_t) n + 1, sizeof(double));
  bm_gm      = (double *) R_alloc((size_t) D + 1, sizeof(double));
  bm_vp      = (double *) R_alloc((size_t) D + 1, sizeof(double));
  bm_bcov    = (double *) R_alloc((size_t) D * q + 1, sizeof(double));
  bm_n_fb    = (int *)    R_alloc((size_t) D * B + 1, sizeof(int));
  bm_gbar    = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bm_t2      = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bm_a       = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bm_b       = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bm_status  = (int *)    R_alloc((size_t) B + 1, sizeof(int));

  // Per-thread workspace.
  int nt = threads_for(n_thr, (double) D * n * (p + 4) * 8, 1e6);
  bm_mem_sz = flm_doubles(n, p, 0);
  bm_buf_sz = 4 * (size_t) n + (size_t) p + (size_t) n * p + (size_t) q + 1;
  bm_mem    = (double *) R_alloc((size_t) nt * bm_mem_sz, sizeof(double));
  bm_imem   = (int *)    R_alloc((size_t) nt * (2 * (size_t) p + (size_t) n + B + 2), sizeof(int));
  bm_buf    = (double *) R_alloc((size_t) nt * bm_buf_sz, sizeof(double));
  int tasks = D < 100 ? 100 : D;
  int stasks = n < 100 ? 100 : n;

  // 1-2. Raw values, totals.
  bm_raw_totals();
  run_parallel(bm_raw_worker, nt, tasks);
  run_parallel(bm_tot_worker, nt, stasks);

  // The whole design's QR, for every feature without zero inflation.
  flm_t full;
  if (!bm_zi) {
    double *fmem  = (double *) R_alloc(bm_mem_sz, sizeof(double));
    int    *fimem = (int *)    R_alloc((size_t) p + 1, sizeof(int));
    flm_init(&full, fmem, fimem, n, p, 0);
    flm_decompose(&full, bm_X, n, NULL, n, p, bm_qr_tol);
    bm_full = &full;
  }

  // 3. Fits.
  run_parallel(bm_fit_worker, nt, tasks);

  // Batches adjusted in fewer than two features are switched off.
  int any_feature = 0;
  for (int b = 0; b < B; b++) {
    int cnt = 0;
    for (int f = 0; f < D; f++) cnt += bm_ind[f + (size_t) D * b];
    if (cnt < 2) for (int f = 0; f < D; f++) bm_ind[f + (size_t) D * b] = 0;
  }
  for (int f = 0; f < D && !any_feature; f++)
    for (int b = 0; b < B; b++) if (bm_ind[f + (size_t) D * b]) { any_feature = 1; break; }
  if (!any_feature) status = 1;

  if (!status) {
    for (int b = 0; b < B; b++) {
      int any = 0;
      for (int f = 0; f < D; f++) if (bm_ind[f + (size_t) D * b]) { any = 1; break; }
      if (!any) { LOGICAL(s_bad)[b] = 1; status = 2; }
    }
  }

  if (!status) {

    // 4. Priors per batch, over the features adjusted in it.
    double *v = (double *) R_alloc((size_t) D + 1, sizeof(double));
    double *w = (double *) R_alloc((size_t) D + 1, sizeof(double));
    for (int b = 0; b < B; b++) {
      int m = 0;
      for (int f = 0; f < D; f++)
        if (bm_ind[f + (size_t) D * b]) {
          v[m]   = bm_ghat[f + (size_t) D * b];
          w[m++] = bm_dhat[f + (size_t) D * b];
        }
      bm_gbar[b] = r_mean(v, m);
      bm_t2[b]   = r_var(v, m);
      double mm = r_mean(w, m), s2 = r_var(w, m);
      if (s2 == 0) s2 = 1;
      bm_a[b] = (2 * s2 + mm * mm) / s2;
      bm_b[b] = (mm * s2 + pow(mm, 3.0)) / s2;
    }

    // 5. Shrinkage, one batch per task.
    int bt = threads_for(n_thr, (double) B * n * D * 50, 1e6);
    if (bt > B) bt = B;
    run_parallel(bm_shrink_worker, bt, B * 100);

    for (int b = 0; b < B && !status; b++)
      if (bm_status[b]) { status = 3; INTEGER(s_wh)[0] = b + 1; }

    for (int b = 0; b < B && !status; b++)
      for (int f = 0; f < D; f++)
        if (bm_ind[f + (size_t) D * b] &&
            (ISNAN(bm_gstar[f + (size_t) D * b]) || ISNAN(bm_dstar[f + (size_t) D * b]))) {
          status = 4;
          INTEGER(s_wh)[0] = b + 1;
          INTEGER(s_wh)[1] = f + 1;
          break;
        }
  }

  // Not shrunk: NA, as MMUPHin's matrices.
  for (size_t i = 0; i < (size_t) D * B; i++)
    if (!bm_ind[i]) bm_ghat[i] = bm_dhat[i] = bm_gstar[i] = bm_dstar[i] = R_NaN;

  if (!status) {
    // 6-7. Adjust, back-transform, re-close.
    run_parallel(bm_adjust_worker, nt, tasks);
    run_parallel(bm_adj_tot_worker, nt, stasks);
    run_parallel(bm_close_worker, nt, tasks);
    run_parallel(bm_adj_tot_worker, nt, stasks);
    run_parallel(bm_mean1_worker, nt, tasks);
  }

  SET_VECTOR_ELT(result, 7, ScalarInteger(status));
  UNPROTECT(1);
  return result;
}
