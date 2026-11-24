// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// batch_combat() (R/batch_combat.r): a C port of ComBat() and its helpers
// aprior(), bprior(), postmean(), postvar(), it.sol() and int.eprior() from
// sva 3.60.0 (R/ComBat.R, R/helper.R), restructured for ecodive, see
// plans/batch_correction.md.
//
// sva's notice:
//
//   sva: Copyright (c) Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker,
//   Elana J. Fertig, Andrew E. Jaffe, Yuqing Zhang, John D. Storey and
//   Leonardo Collado Torres (the authors named in sva's DESCRIPTION, which
//   has no other copyright line). Licensed under the Artistic License 2.0,
//   below.
//
//   This file is a modified version of that code (Artistic License 2.0,
//   section 4), distributed as part of the package ecodive, which is not
//   sva and doesn't prevent sva from being installed. How it differs: the
//   whole computation is one pass of threaded C per step (features in
//   blocks, the shrinkage per batch, the non-parametric integral per batch
//   and feature) over a table read in any of ecodive's formats, with every
//   sum in a fixed order; the least-squares fits use one pivoted QR of the
//   design (src/feature_lm.c) where ComBat() solves the normal equations;
//   the non-parametric integral weighs the other features by their
//   likelihoods on the log scale (relative to the largest), where
//   int.eprior() multiplies densities that can underflow to 0 and give NaN;
//   the shrinkage iteration takes its tolerance and an iteration cap from
//   the caller (it.sol() fixes conv = 1e-4 and has no cap) and stops with
//   an error on NaN; missing values are refused (Beta.NA() is not ported);
//   and the prior plots are not drawn (plot.ecodive_batch() draws other
//   diagnostics). The user-visible differences are listed in
//   ?batch_combat ("Differences from sva").
//
//   Artistic License 2.0
//
//   Copyright (c) 2000-2006, The Perl Foundation.
//
//   Everyone is permitted to copy and distribute verbatim copies of this
//   license document, but changing it is not allowed.
//
//   Preamble
//   ********
//
//   This license establishes the terms under which a given free software
//   Package may be copied, modified, distributed, and/or redistributed.  The
//   intent is that the Copyright Holder maintains some artistic control over
//   the development of that Package while still keeping the Package
//   available as open source and free software.
//
//   You are always permitted to make arrangements wholly outside of this
//   license directly with the Copyright Holder of a given Package.  If the
//   terms of this license do not permit the full use that you propose to
//   make of the Package, you should contact the Copyright Holder and seek a
//   different licensing arrangement.
//
//   Definitions
//   ***********
//
//   "Copyright Holder" means the individual(s) or organization(s) named in
//   the copyright notice for the entire Package.
//
//   "Contributor" means any party that has contributed code or other
//   material to the Package, in accordance with the Copyright Holder's
//   procedures.
//
//   "You" and "your" means any person who would like to copy, distribute, or
//   modify the Package.
//
//   "Package" means the collection of files distributed by the Copyright
//   Holder, and derivatives of that collection and/or of those files.  A
//   given Package may consist of either the Standard Version, or a Modified
//   Version.
//
//   "Distribute" means providing a copy of the Package or making it
//   accessible to anyone else, or in the case of a company or organization,
//   to others outside of your company or organization.
//
//   "Distributor Fee" means any fee that you charge for Distributing this
//   Package or providing support for this Package to another party.  It does
//   not mean licensing fees.
//
//   "Standard Version" refers to the Package if it has not been modified, or
//   has been modified only in ways explicitly requested by the Copyright
//   Holder.
//
//   "Modified Version" means the Package, if it has been changed, and such
//   changes were not explicitly requested by the Copyright Holder.
//
//   "Original License" means this Artistic License as Distributed with the
//   Standard Version of the Package, in its current version or as it may be
//   modified by The Perl Foundation in the future.
//
//   "Source" form means the source code, documentation source, and
//   configuration files for the Package.
//
//   "Compiled" form means the compiled bytecode, object code, binary, or any
//   other form resulting from mechanical transformation or translation of
//   the Source form.
//
//   Permission for Use and Modification Without Distribution
//   ********************************************************
//
//   (1) You are permitted to use the Standard Version and create and use
//   Modified Versions for any purpose without restriction, provided that you
//   do not Distribute the Modified Version.
//
//   Permissions for Redistribution of the Standard Version
//   ******************************************************
//
//   (2) You may Distribute verbatim copies of the Source form of the
//   Standard Version of this Package in any medium without restriction,
//   either gratis or for a Distributor Fee, provided that you duplicate all
//   of the original copyright notices and associated disclaimers.  At your
//   discretion, such verbatim copies may or may not include a Compiled form
//   of the Package.
//
//   (3) You may apply any bug fixes, portability changes, and other
//   modifications made available from the Copyright Holder.  The resulting
//   Package will still be considered the Standard Version, and as such will
//   be subject to the Original License.
//
//   Distribution of Modified Versions of the Package as Source
//   **********************************************************
//
//   (4) You may Distribute your Modified Version as Source (either gratis or
//   for a Distributor Fee, and with or without a Compiled form of the
//   Modified Version) provided that you clearly document how it differs from
//   the Standard Version, including, but not limited to, documenting any
//   non-standard features, executables, or modules, and provided that you do
//   at least ONE of the following:
//
//   (a) make the Modified Version available to the Copyright Holder of the
//   Standard Version, under the Original License, so that the Copyright
//   Holder may include your modifications in the Standard Version.
//
//   (b) ensure that installation of your Modified Version does not prevent
//   the user installing or running the Standard Version.  In addition, the
//   Modified Version must bear a name that is different from the name of the
//   Standard Version.
//
//   (c) allow anyone who receives a copy of the Modified Version to make the
//   Source form of the Modified Version available to others under
//
//   (i) the Original License or
//
//   (ii) a license that permits the licensee to freely copy, modify and
//   redistribute the Modified Version using the same licensing terms that
//   apply to the copy that the licensee received, and requires that the
//   Source form of the Modified Version, and of any works derived from it,
//   be made freely available in that license fees are prohibited but
//   Distributor Fees are allowed.
//
//   Distribution of Compiled Forms of the Standard Version or Modified
//   ******************************************************************
//   Versions without the Source
//   ***************************
//
//   (5) You may Distribute Compiled forms of the Standard Version without
//   the Source, provided that you include complete instructions on how to
//   get the Source of the Standard Version.  Such instructions must be valid
//   at the time of your distribution.  If these instructions, at any time
//   while you are carrying out such distribution, become invalid, you must
//   provide new instructions on demand or cease further distribution.  If
//   you provide valid instructions or cease distribution within thirty days
//   after you become aware that the instructions are invalid, then you do
//   not forfeit any of your rights under this license.
//
//   (6) You may Distribute a Modified Version in Compiled form without the
//   Source, provided that you comply with Section 4 with respect to the
//   Source of the Modified Version.
//
//   Aggregating or Linking the Package
//   **********************************
//
//   (7) You may aggregate the Package (either the Standard Version or
//   Modified Version) with other packages and Distribute the resulting
//   aggregation provided that you do not charge a licensing fee for the
//   Package.  Distributor Fees are permitted, and licensing fees for other
//   components in the aggregation are permitted.  The terms of this license
//   apply to the use and Distribution of the Standard or Modified Versions
//   as included in the aggregation.
//
//   (8) You are permitted to link Modified and Standard Versions with other
//   works, to embed the Package in a larger work of your own, or to build
//   stand-alone binary or bytecode versions of applications that include the
//   Package, and Distribute the result without restriction, provided the
//   result does not expose a direct interface to the Package.
//
//   Items That are Not Considered Part of a Modified Version
//   ********************************************************
//
//   (9) Works (including, but not limited to, modules and scripts) that
//   merely extend or make use of the Package, do not, by themselves, cause
//   the Package to be a Modified Version.  In addition, such works are not
//   considered parts of the Package itself, and are not subject to the terms
//   of this license.
//
//   General Provisions
//   ******************
//
//   (10) Any use, modification, and distribution of the Standard or Modified
//   Versions is governed by this Artistic License.  By using, modifying or
//   distributing the Package, you accept this license.  Do not use, modify,
//   or distribute the Package, if you do not accept this license.
//
//   (11) If your Modified Version has been derived from a Modified Version
//   made by someone other than you, you are nevertheless required to ensure
//   that your Modified Version complies with the requirements of this
//   license.
//
//   (12) This license does not grant you the right to use any trademark,
//   service mark, tradename, or logo of the Copyright Holder.
//
//   (13) This license includes the non-exclusive, worldwide, free-of-charge
//   patent license to make, have made, use, offer to sell, sell, import and
//   otherwise transfer the Package with respect to any patent claims
//   licensable by the Copyright Holder that are necessarily infringed by the
//   Package.  If you institute patent litigation (including a cross-claim or
//   counterclaim) against any party alleging that the Package constitutes
//   direct or contributory patent infringement, then this Artistic License
//   to you shall terminate on the date that such litigation is filed.
//
//   (14) Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT
//   HOLDER AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED
//   WARRANTIES.  THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//   PARTICULAR PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT
//   PERMITTED BY YOUR LOCAL LAW.  UNLESS REQUIRED BY LAW, NO COPYRIGHT
//   HOLDER OR CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT,
//   INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE
//   OF THE PACKAGE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The steps, with every sum in a fixed order (so results are the same for
// any number of threads):
//
//   1. Per feature (threads): its values (after `norm`), their means per
//      batch, and whether they are constant within any batch of two or
//      more samples (such features are returned unchanged, as in sva).
//   2. Per feature (threads): the least-squares fit on the whole design
//      (batch indicators, the reference batch's being all ones, then the
//      covariates), the grand mean (the batch coefficients weighted by
//      batch size, or the reference batch's), the pooled variance (divisor
//      n, or the reference batch's size), the standardised values, and each
//      batch's location (mean) and scale (variance, n - 1).
//   3. Per batch (main thread): the priors over features, in feature order.
//   4. Parametric: per batch (threads), it.sol()'s iteration over features
//      in order. Non-parametric: per batch and feature (threads), the
//      weighted means over every other feature, in feature order.
//   5. Per feature (threads): the adjustment, then each batch's means.

#include "ecodive.h"


static int           bc_n, bc_D, bc_B, bc_q, bc_p, bc_ref, bc_mean_only, bc_K;
static double        bc_tol;
static int           bc_max_iters;
static const da_table_t *bc_tab;
static const int    *bc_batch;    // n: each sample's batch, 0-based
static const double *bc_X;        // n x p: batch indicators, then covariates
static const flm_t  *bc_full;     // its QR
static double       *bc_out;      // n x D: values, then the result
static double       *bc_s;        // n x D: standardised values
static int          *bc_keep;     // D: adjusted (not constant within a batch)
static const int    *bc_kept;     // K: the kept features, in order
static double       *bc_ghat, *bc_dhat, *bc_gstar, *bc_dstar;   // D x B
static double       *bc_gm, *bc_vp, *bc_bcov;   // D, D, D x q
static double       *bc_mean0, *bc_mean1, *bc_overall;
static double       *bc_gbar, *bc_t2, *bc_a, *bc_b;   // B: priors
static int          *bc_iters, *bc_status;            // B
static int          *bc_np_bad;   // B x K: non-parametric weights all 0
static const int    *bc_bsam;     // samples ordered by batch ...
static const int    *bc_bpos;     // ... batch b's are bsam[bpos[b] .. bpos[b + 1] - 1]
static double       *bc_buf;      // per thread: scratch
static size_t        bc_buf_sz;


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


// genefilter's rowVars() of one row: rowSums((x - rowMeans(x))^2) / (n - 1),
// with rowMeans() and rowSums() summing in long double.
static double row_var (const double *x, int n) {
  if (n < 2) return R_NaN;
  long double m = 0;
  for (int i = 0; i < n; i++) m += x[i];
  double xm = (double) (m / n);
  long double s = 0;
  for (int i = 0; i < n; i++) { double d = x[i] - xm; s += d * d; }
  return (double) s / (n - 1);
}


// rowMeans() of one row.
static double row_mean (const double *x, int n) {
  long double m = 0;
  for (int i = 0; i < n; i++) m += x[i];
  return (double) (m / n);
}


// A row's covariate part, as R's product with the design whose batch
// columns are zeroed sums it (in column order).
static inline double bc_cov_part (int i, const double *beta) {
  double c = 0;
  for (int j = 0; j < bc_q; j++) c += bc_X[i + (size_t) bc_n * (bc_B + j)] * beta[j];
  return c;
}



//======================================================
// 1. Values, means, and constant features.
//======================================================
static void *bc_read_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bc_n, D = bc_D, B = bc_B;

  for (int f = thread_i; f < D; f += n_threads) {
    double *y = bc_out + (size_t) n * f;
    da_table_column(bc_tab, f, y);

    int keep = 1;
    for (int b = 0; b < B; b++) {
      int lo = bc_bpos[b], hi = bc_bpos[b + 1];
      long double m = 0;
      for (int k = lo; k < hi; k++) m += y[bc_bsam[k]];
      bc_mean0[f + (size_t) D * b] = (double) (m / (hi - lo));
      // var(x) == 0 in a batch of two or more samples.
      if (hi - lo > 1) {
        double v0 = y[bc_bsam[lo]];
        int same = 1;
        for (int k = lo + 1; k < hi && same; k++) if (y[bc_bsam[k]] != v0) same = 0;
        if (same) keep = 0;
      }
    }
    long double m = 0;
    for (int s = 0; s < n; s++) m += y[s];
    bc_overall[f] = (double) (m / n);
    bc_keep[f] = keep;
  }
  return NULL;
}



//======================================================
// 2. The fits, standardisation, and batch moments.
//======================================================
static void *bc_fit_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bc_n, D = bc_D, B = bc_B, q = bc_q, p = bc_p, K = bc_K;

  double *buf  = bc_buf + (size_t) thread_i * bc_buf_sz;
  double *qty  = buf;          // n
  double *fit  = qty + n;      // n
  double *coef = fit + n;      // p
  double *tmp  = coef + p;     // n
  double *wrk  = tmp + n;      // n + 3 p: flm_fit()'s scratch

  flm_t qq = *bc_full;
  qq.work  = wrk;

  for (int kk = thread_i; kk < K; kk += n_threads) {

    int f = bc_kept[kk];
    const double *y = bc_out + (size_t) n * f;
    double       *s = bc_s   + (size_t) n * f;

    double rss;
    flm_fit(&qq, y, qty, coef, &rss);

    for (int i = 0; i < n; i++) {
      double fv = 0;
      for (int j = 0; j < p; j++) fv += bc_X[i + (size_t) n * j] * coef[j];
      fit[i] = fv;
    }

    // The grand mean, and the pooled variance (divisor n, or n_ref).
    double gm = 0, vp = 0;
    if (bc_ref < 0) {
      for (int b = 0; b < B; b++)
        gm += ((double) (bc_bpos[b + 1] - bc_bpos[b]) / n) * coef[b];
      double w = 1.0 / n;
      for (int i = 0; i < n; i++) { double r = y[i] - fit[i]; vp += w * (r * r); }
    }
    else {
      gm = coef[bc_ref];
      int lo = bc_bpos[bc_ref], hi = bc_bpos[bc_ref + 1];
      double w = 1.0 / (hi - lo);
      for (int k = lo; k < hi; k++) {
        int i = bc_bsam[k];
        double r = y[i] - fit[i];
        vp += w * (r * r);
      }
    }

    const double *bc = coef + B;
    for (int j = 0; j < q; j++) bc_bcov[f + (size_t) D * j] = bc[j];
    bc_gm[f] = gm;
    bc_vp[f] = vp;

    double sd = sqrt(vp);
    for (int i = 0; i < n; i++) s[i] = (y[i] - (gm + bc_cov_part(i, bc))) / sd;

    // Each batch's location (by least squares on the batch design) and
    // scale (variance).
    double gref = 0;
    if (bc_ref >= 0) {
      int lo = bc_bpos[bc_ref], hi = bc_bpos[bc_ref + 1];
      double t = 0;
      for (int k = lo; k < hi; k++) t += s[bc_bsam[k]];
      gref = t / (hi - lo);
    }
    for (int b = 0; b < B; b++) {
      int lo = bc_bpos[b], hi = bc_bpos[b + 1], m = hi - lo;
      double t = 0;
      for (int k = lo; k < hi; k++) { tmp[k - lo] = s[bc_bsam[k]]; t += tmp[k - lo]; }
      double g = t / m;
      if (bc_ref >= 0 && b != bc_ref) g -= gref;
      if (b == bc_ref) g = gref;
      bc_ghat[f + (size_t) D * b] = g;
      bc_dhat[f + (size_t) D * b] = bc_mean_only ? 1 : row_var(tmp, m);
    }
  }

  return NULL;
}



//======================================================
// 4a. it.sol(), for one batch per task.
//======================================================
static void *bc_itsol_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bc_n, D = bc_D, B = bc_B, K = bc_K;

  for (int b = thread_i; b < B; b += n_threads) {

    double *gs = bc_gstar + (size_t) D * b;
    double *ds = bc_dstar + (size_t) D * b;
    const double *gh = bc_ghat + (size_t) D * b;
    const double *dh = bc_dhat + (size_t) D * b;
    double gbar = bc_gbar[b], t2 = bc_t2[b], a = bc_a[b], bb = bc_b[b];
    int lo = bc_bpos[b], hi = bc_bpos[b + 1];
    double nn = hi - lo;

    for (int kk = 0; kk < K; kk++) { int f = bc_kept[kk]; gs[f] = gh[f]; ds[f] = dh[f]; }

    int count = 0;
    bc_status[b] = 0;

    for (;;) {
      double change = R_NegInf;
      int    nan    = 0;

      for (int kk = 0; kk < K; kk++) {
        int f = bc_kept[kk];
        double g_old = gs[f], d_old = ds[f];
        double g_new = (t2 * nn * gh[f] + d_old * gbar) / (t2 * nn + d_old);

        const double *s = bc_s + (size_t) n * f;
        long double sum2 = 0;
        for (int k = lo; k < hi; k++) { double r = s[bc_bsam[k]] - g_new; sum2 += r * r; }
        double d_new = (0.5 * (double) sum2 + bb) / (nn / 2 + a - 1);

        // The divisors are signed, as in sva.
        double r1 = fabs(g_new - g_old) / g_old;
        double r2 = fabs(d_new - d_old) / d_old;
        if (ISNAN(r1) || ISNAN(r2)) nan = 1;
        else {
          if (r1 > change) change = r1;
          if (r2 > change) change = r2;
        }
        gs[f] = g_new;
        ds[f] = d_new;
      }

      count++;
      if (nan)                   { bc_status[b] = 2; break; }
      if (!(change > bc_tol))    break;
      if (count >= bc_max_iters) { bc_status[b] = 1; break; }
    }

    bc_iters[b] = count;
  }

  return NULL;
}



//======================================================
// 4b. int.eprior(), per batch and feature: the other
//     features' locations and scales, weighed by the
//     likelihood of this feature's values under each.
//======================================================
static void *bc_nonpar_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bc_n, D = bc_D, B = bc_B, K = bc_K;

  double *buf = bc_buf + (size_t) thread_i * bc_buf_sz;
  double *x   = buf;        // n
  double *ll  = x + n;      // K

  for (int t = thread_i; t < B * K; t += n_threads) {

    int b = t / K, kk = t % K, f = bc_kept[kk];
    int lo = bc_bpos[b], hi = bc_bpos[b + 1], m = hi - lo;
    const double *s  = bc_s + (size_t) n * f;
    const double *gh = bc_ghat + (size_t) D * b;
    const double *dh = bc_dhat + (size_t) D * b;

    for (int k = lo; k < hi; k++) x[k - lo] = s[bc_bsam[k]];

    // log LH = -(n/2) log(2 pi d) - sum2 / (2 d).
    double mx = R_NegInf;
    for (int jj = 0; jj < K; jj++) {
      ll[jj] = R_NaN;
      if (jj == kk) continue;
      int j = bc_kept[jj];
      double g = gh[j], d = bc_mean_only ? 1 : dh[j];
      double sum2 = 0;
      for (int i = 0; i < m; i++) { double r = x[i] - g; sum2 += r * r; }
      double v = -(m / 2.0) * log(2 * M_PI * d) - sum2 / (2 * d);
      if (!ISNAN(v) && d > 0) {
        ll[jj] = v;
        if (v > mx) mx = v;
      }
    }

    long double sw = 0, sg = 0, sd = 0;
    if (R_FINITE(mx))
      for (int jj = 0; jj < K; jj++) {
        if (ISNAN(ll[jj])) continue;
        int j = bc_kept[jj];
        double w = exp(ll[jj] - mx);
        double d = bc_mean_only ? 1 : dh[j];
        sw += w;
        sg += gh[j] * w;
        sd += d * w;
      }

    if (sw > 0) {
      bc_gstar[f + (size_t) D * b] = (double) (sg / sw);
      bc_dstar[f + (size_t) D * b] = (double) (sd / sw);
    }
    else {
      bc_gstar[f + (size_t) D * b] = bc_dstar[f + (size_t) D * b] = R_NaN;
      bc_np_bad[b + (size_t) B * kk] = 1;
    }
  }

  return NULL;
}



//======================================================
// 5. Adjust, then the batch means.
//======================================================
static void *bc_adjust_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = bc_n, D = bc_D, B = bc_B, q = bc_q;
  double *bcp = bc_buf + (size_t) thread_i * bc_buf_sz;   // q

  for (int f = thread_i; f < D; f += n_threads) {

    double *y = bc_out + (size_t) n * f;

    if (bc_keep[f]) {
      const double *s = bc_s + (size_t) n * f;
      double sd = sqrt(bc_vp[f]), gm = bc_gm[f];
      for (int j = 0; j < q; j++) bcp[j] = bc_bcov[f + (size_t) D * j];
      for (int i = 0; i < n; i++) {
        int b = bc_batch[i];
        if (b == bc_ref) continue;   // the reference batch is unchanged
        double v = (s[i] - bc_gstar[f + (size_t) D * b]) / sqrt(bc_dstar[f + (size_t) D * b]);
        y[i] = v * sd + (gm + bc_cov_part(i, bcp));
      }
    }

    for (int b = 0; b < B; b++) {
      int lo = bc_bpos[b], hi = bc_bpos[b + 1];
      long double m = 0;
      for (int k = lo; k < hi; k++) m += y[bc_bsam[k]];
      bc_mean1[f + (size_t) D * b] = (double) (m / (hi - lo));
    }
  }

  return NULL;
}



//======================================================
// R interface, for batch_combat().
//
//   counts, margin, norm, pseudocount   The table, as normalize().
//   batch          n: each sample's batch, 0-based.
//   X              n x (B + q): the batch indicators (the reference
//                  batch's column all ones), then the covariates.
//   n_batches      B.
//   ref            The reference batch, 0-based, or -1.
//   parametric     1: it.sol(); 0: int.eprior().
//   mean_only      1: no scale adjustment.
//   tol, max_iters it.sol()'s tolerance, and a cap on its iterations.
//
// Returns list(adjusted (n x D), keep (D logical), gamma_hat, gamma_star,
// delta_hat, delta_star (D x B), n_iters (B), status (0: ok; 1: fewer
// than two features vary within every batch; 2: a batch's iteration
// didn't converge, see `which`; 3: NaN in a batch's iteration, see
// `which`; 4: a non-parametric estimate had no weight, see `which`;
// 5: the design is rank deficient), which (batch, feature, 1-based),
// means_before, means_after (D x B), overall (D)).
//======================================================
SEXP C_batch_combat(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_pseudocount,
    SEXP sexp_batch, SEXP sexp_X, SEXP sexp_n_batches, SEXP sexp_ref,
    SEXP sexp_parametric, SEXP sexp_mean_only, SEXP sexp_tol,
    SEXP sexp_max_iters, SEXP sexp_n_threads ) {

  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, asInteger(sexp_norm), asReal(sexp_pseudocount), n_thr, &tab);
  int n = tab.n, D = tab.D;

  int B = asInteger(sexp_n_batches);
  SEXP dims = getAttrib(sexp_X, R_DimSymbol);
  if (INTEGER(dims)[0] != n) error("The design must have one row per sample.");  // # nocov
  int p = INTEGER(dims)[1];
  int q = p - B;
  if (LENGTH(sexp_batch) != n) error("`batch` must have one value per sample.");  // # nocov

  bc_n = n; bc_D = D; bc_B = B; bc_q = q; bc_p = p;
  bc_ref       = asInteger(sexp_ref);
  bc_mean_only = asLogical(sexp_mean_only);
  bc_tol       = asReal(sexp_tol);
  bc_max_iters = asInteger(sexp_max_iters);
  bc_tab       = &tab;
  bc_batch     = INTEGER(sexp_batch);
  bc_X         = REAL(sexp_X);
  int parametric = asLogical(sexp_parametric);

  const char *nms[] = {
    "adjusted", "keep", "gamma_hat", "gamma_star", "delta_hat", "delta_star",
    "n_iters", "status", "which", "means_before", "means_after", "overall", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, nms));
  SEXP s_out  = allocMatrix(REALSXP, n, D); SET_VECTOR_ELT(result, 0, s_out);
  SEXP s_keep = allocVector(LGLSXP, D);     SET_VECTOR_ELT(result, 1, s_keep);
  SEXP s_gh   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 2, s_gh);
  SEXP s_gs   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 3, s_gs);
  SEXP s_dh   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 4, s_dh);
  SEXP s_ds   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 5, s_ds);
  SEXP s_it   = allocVector(INTSXP, B);     SET_VECTOR_ELT(result, 6, s_it);
  SEXP s_wh   = allocVector(INTSXP, 2);     SET_VECTOR_ELT(result, 8, s_wh);
  SEXP s_m0   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 9, s_m0);
  SEXP s_m1   = allocMatrix(REALSXP, D, B); SET_VECTOR_ELT(result, 10, s_m1);
  SEXP s_ov   = allocVector(REALSXP, D);    SET_VECTOR_ELT(result, 11, s_ov);
  int status  = 0;

  bc_out   = REAL(s_out);
  bc_keep  = LOGICAL(s_keep);
  bc_ghat  = REAL(s_gh);  bc_gstar = REAL(s_gs);
  bc_dhat  = REAL(s_dh);  bc_dstar = REAL(s_ds);
  bc_iters = INTEGER(s_it);
  bc_mean0 = REAL(s_m0);  bc_mean1 = REAL(s_m1);
  bc_overall = REAL(s_ov);
  for (int b = 0; b < B; b++) bc_iters[b] = 0;
  INTEGER(s_wh)[0] = INTEGER(s_wh)[1] = NA_INTEGER;
  for (size_t i = 0; i < (size_t) D * B; i++)
    bc_ghat[i] = bc_dhat[i] = bc_gstar[i] = bc_dstar[i] = bc_mean1[i] = R_NaN;

  // Samples by batch.
  int *bpos = (int *) R_alloc((size_t) B + 1, sizeof(int));
  int *bsam = (int *) R_alloc((size_t) n + 1, sizeof(int));
  memset(bpos, 0, ((size_t) B + 1) * sizeof(int));
  for (int s = 0; s < n; s++) bpos[bc_batch[s] + 1]++;
  for (int b = 0; b < B; b++) bpos[b + 1] += bpos[b];
  int *fill = (int *) R_alloc((size_t) B + 1, sizeof(int));
  memcpy(fill, bpos, ((size_t) B + 1) * sizeof(int));
  for (int s = 0; s < n; s++) bsam[fill[bc_batch[s]]++] = s;
  bc_bpos = bpos;
  bc_bsam = bsam;

  bc_gm     = (double *) R_alloc((size_t) D + 1, sizeof(double));
  bc_vp     = (double *) R_alloc((size_t) D + 1, sizeof(double));
  bc_bcov   = (double *) R_alloc((size_t) D * q + 1, sizeof(double));
  bc_gbar   = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bc_t2     = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bc_a      = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bc_b      = (double *) R_alloc((size_t) B + 1, sizeof(double));
  bc_status = (int *)    R_alloc((size_t) B + 1, sizeof(int));

  int nt = threads_for(n_thr, (double) D * n * (p + 4) * 4, 1e6);
  bc_buf_sz = 4 * (size_t) n + 4 * (size_t) p + (size_t) D + (size_t) q + 1;
  bc_buf    = (double *) R_alloc((size_t) n_thr * bc_buf_sz, sizeof(double));
  int tasks = D < 100 ? 100 : D;

  // 1. Values; which features vary within every batch.
  run_parallel(bc_read_worker, nt, tasks);
  int K = 0;
  for (int f = 0; f < D; f++) K += bc_keep[f];
  int *kept = (int *) R_alloc((size_t) K + 1, sizeof(int));
  K = 0;
  for (int f = 0; f < D; f++) if (bc_keep[f]) kept[K++] = f;
  bc_K    = K;
  bc_kept = kept;
  if (K < 2) status = 1;

  // The design's QR.
  size_t  msz   = flm_doubles(n, p, 0);
  double *fmem  = (double *) R_alloc(msz, sizeof(double));
  int    *fimem = (int *)    R_alloc((size_t) p + 1, sizeof(int));
  flm_t   full;
  flm_init(&full, fmem, fimem, n, p, 0);
  flm_decompose(&full, bc_X, n, NULL, n, p, 1e-7);
  bc_full = &full;
  if (!status && full.rank < p) status = 5;

  if (!status) {

    // 2. Fits.
    bc_s = (double *) R_alloc((size_t) n * D + 1, sizeof(double));
    int kt = threads_for(n_thr, (double) K * n * (p + 4) * 4, 1e6);
    run_parallel(bc_fit_worker, kt, K < 100 ? 100 : K);

    // 3. Priors per batch, over the kept features.
    double *v = (double *) R_alloc((size_t) K + 1, sizeof(double));
    for (int b = 0; b < B; b++) {
      for (int kk = 0; kk < K; kk++) v[kk] = bc_ghat[kept[kk] + (size_t) D * b];
      bc_gbar[b] = row_mean(v, K);
      bc_t2[b]   = row_var(v, K);
      for (int kk = 0; kk < K; kk++) v[kk] = bc_dhat[kept[kk] + (size_t) D * b];
      double m = r_mean(v, K), s2 = r_var(v, K);
      bc_a[b] = (2 * s2 + m * m) / s2;
      bc_b[b] = (m * s2 + pow(m, 3.0)) / s2;
    }

    // 4. Shrinkage.
    if (parametric && bc_mean_only) {
      for (int b = 0; b < B; b++)
        for (int kk = 0; kk < K; kk++) {
          size_t i = kept[kk] + (size_t) D * b;
          bc_gstar[i] = (bc_t2[b] * 1 * bc_ghat[i] + 1 * bc_gbar[b]) / (bc_t2[b] * 1 + 1);
          bc_dstar[i] = 1;
        }
    }
    else if (parametric) {
      int bt = threads_for(n_thr, (double) K * n * 50, 1e6);
      if (bt > B) bt = B;
      run_parallel(bc_itsol_worker, bt, B * 100);
      for (int b = 0; b < B && !status; b++)
        if (bc_status[b]) { status = bc_status[b] == 1 ? 2 : 3; INTEGER(s_wh)[0] = b + 1; }
    }
    else {
      bc_np_bad = (int *) R_alloc((size_t) B * K + 1, sizeof(int));
      memset(bc_np_bad, 0, ((size_t) B * K + 1) * sizeof(int));
      int npt = threads_for(n_thr, (double) K * K * (n + B) * 3, 1e6);
      run_parallel(bc_nonpar_worker, npt, B * K < 100 ? 100 : B * K);
      for (int kk = 0; kk < K && !status; kk++)
        for (int b = 0; b < B; b++)
          if (bc_np_bad[b + (size_t) B * kk]) {
            status = 4;
            INTEGER(s_wh)[0] = b + 1;
            INTEGER(s_wh)[1] = kept[kk] + 1;
            break;
          }
    }

    // The reference batch: no location or scale change.
    if (bc_ref >= 0)
      for (int kk = 0; kk < K; kk++) {
        size_t i = kept[kk] + (size_t) D * bc_ref;
        bc_gstar[i] = 0;
        bc_dstar[i] = 1;
      }

    // 5. Adjust.
    if (!status) run_parallel(bc_adjust_worker, nt, tasks);
  }

  SET_VECTOR_ELT(result, 7, ScalarInteger(status));
  UNPROTECT(1);
  return result;
}
