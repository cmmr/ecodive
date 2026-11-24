// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The per-feature sums behind decontam() and decontam_not() (R/decontam.r):
// a C and R port of isContaminant(), isNotContaminant(),
// isContaminantFrequency(), isContaminantPrevalence() and fish.combine()
// from decontam 1.32.0 (R/decontam.R; identical in 1.22.0), restructured
// for ecodive, see plans/decontam.md. This file computes, for each feature
// and batch, what those functions' lm() fits and 2 x 2 tables are built
// from; R/decontam.r turns them into scores.
//
// decontam's notice:
//
//   decontam: Copyright (c) Benjamin Callahan and Nicole Marie Davis (the
//   authors named in decontam's DESCRIPTION, which has no other copyright
//   line). Licensed under the Artistic License 2.0, below.
//
//   This file is a modified version of that code (Artistic License 2.0,
//   section 4), distributed as part of the package ecodive, which is not
//   decontam and doesn't prevent decontam from being installed. How it
//   differs: isContaminantFrequency()'s two lm() fits (the intercept-only
//   model and the model with slope -1 on log concentration) are replaced by
//   their residual sums of squares, computed in two passes over each
//   feature's non-zero values (the means, then the sums of squared
//   deviations); isContaminantPrevalence()'s table is replaced by counts of
//   the controls and non-controls a feature is present in; every feature and
//   batch is done in one pass over a sparse table, threaded over features;
//   and samples whose `neg` is NA are left out of both tests (decontam
//   leaves them out of the prevalence table only). The user-visible
//   differences are listed in ?decontam ("Differences from decontam").
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

#include "ecodive.h"


//======================================================
// Shared state.
//======================================================
static const da_table_t *dc_tab;
static const int    *dc_batch;   // per sample: 0 .. B - 1, or -1 (left out)
static const int    *dc_neg;     // per sample: 1 control, 0 not, -1 unknown
static const double *dc_logc;    // per sample: log(conc), or NULL
static int           dc_B;
static double       *dc_buf;     // per thread: n + 2 * B doubles

// Outputs: freq_sum and prev (length D), and B x D matrices.
static double *dc_freq;
static int    *dc_prev;
static int    *dc_k;
static double *dc_ss0;
static double *dc_ss1;
static int    *dc_pneg;
static int    *dc_ppos;


static void *dc_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  const da_table_t *tab = dc_tab;
  int n = tab->n, D = tab->D, B = dc_B;

  double *y   = dc_buf + (size_t) thread_i * ((size_t) n + 2 * (size_t) B + 2);
  double *my  = y + n;        // per batch: mean of log f
  double *ml  = my + B;       // per batch: mean of log f + log c

  for (int f = thread_i; f < D; f += n_threads) {

    size_t  o    = (size_t) f * B;
    int    *k    = dc_k    + o;
    double *ss0  = dc_ss0  + o;
    double *ss1  = dc_ss1  + o;
    int    *pneg = dc_pneg + o;
    int    *ppos = dc_ppos + o;

    for (int b = 0; b < B; b++) {
      k[b] = 0; ss0[b] = 0; ss1[b] = 0; pneg[b] = 0; ppos[b] = 0;
      my[b] = 0; ml[b] = 0;
    }

    // Pass 1, over the feature's stored values in sample order: the mean
    // abundance's sum, the prevalence counts, and the sums for the means
    // of the frequency test's two models.
    long double fsum = 0;
    int         prev = 0;

    for (int e = tab->fpos[f]; e < tab->fpos[f + 1]; e++) {
      int    s = tab->fsam[e];
      int    b = dc_batch[s];
      double v = tab->fval[e];
      if (b < 0) continue;
      fsum += v;
      if (!(v > 0)) continue;
      prev++;
      if      (dc_neg[s] == 1) pneg[b]++;
      else if (dc_neg[s] == 0) {
        ppos[b]++;
        if (dc_logc) {
          double lf = log(v);
          y[s]   = lf;
          k[b]++;
          my[b] += lf;
          ml[b] += lf + dc_logc[s];
        }
      }
    }

    dc_freq[f] = (double) fsum;
    dc_prev[f] = prev;

    if (!dc_logc) continue;

    for (int b = 0; b < B; b++)
      if (k[b] > 0) { my[b] /= k[b]; ml[b] /= k[b]; }

    // Pass 2: the residual sums of squares of log f about its mean
    // (intercept only: frequency independent of concentration) and of
    // log f + log c about its mean (slope -1 on log c: a contaminant).
    for (int e = tab->fpos[f]; e < tab->fpos[f + 1]; e++) {
      int s = tab->fsam[e];
      int b = dc_batch[s];
      if (b < 0 || dc_neg[s] != 0 || !(tab->fval[e] > 0)) continue;
      double d0 = y[s] - my[b];
      double d1 = (y[s] + dc_logc[s]) - ml[b];
      ss0[b] += d0 * d0;
      ss1[b] += d1 * d1;
    }
  }

  return NULL;
}



//======================================================
// R interface. For a table of n samples and D features:
//
//   norm   0 (none) or 1 (percent; each sample over its
//          total)
//   batch  per sample: its batch, 0 .. B - 1, or -1 to
//          leave the sample out (an empty sample)
//   neg    per sample: 1 a negative control, 0 not, -1
//          unknown (left out of both tests)
//   logc   per sample: log of the DNA concentration, or
//          NULL for no frequency test
//
// Returns list(freq, prev, k, ss0, ss1, prev_neg,
// prev_pos): each feature's sum of normalized values and
// number of samples it is present in (over the samples not
// left out); and B x D matrices of, per batch and feature,
// the non-controls with the feature present (k) and the
// two residual sums of squares over them, and the controls
// and non-controls with the feature present.
//
// Each feature's sums run over its samples in order, so
// the results are the same for any number of threads.
//======================================================
SEXP C_decontam(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm, SEXP sexp_batch,
    SEXP sexp_n_batches, SEXP sexp_neg, SEXP sexp_logc, SEXP sexp_n_threads ) {

  int norm  = asInteger(sexp_norm);
  int B     = asInteger(sexp_n_batches);
  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;
  if (B < 1) error("No batches.");  // # nocov

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, norm, 0, n_thr, &tab);
  int n = tab.n, D = tab.D;

  if (LENGTH(sexp_batch) != n || LENGTH(sexp_neg) != n)
    error("`batch` and `neg` must have one value per sample.");  // # nocov
  if (!isNull(sexp_logc) && LENGTH(sexp_logc) != n)
    error("`logc` must have one value per sample.");  // # nocov

  SEXP sexp_freq = PROTECT(allocVector(REALSXP, D));
  SEXP sexp_prev = PROTECT(allocVector(INTSXP,  D));
  SEXP sexp_k    = PROTECT(allocMatrix(INTSXP,  B, D));
  SEXP sexp_ss0  = PROTECT(allocMatrix(REALSXP, B, D));
  SEXP sexp_ss1  = PROTECT(allocMatrix(REALSXP, B, D));
  SEXP sexp_pneg = PROTECT(allocMatrix(INTSXP,  B, D));
  SEXP sexp_ppos = PROTECT(allocMatrix(INTSXP,  B, D));

  // Threads only for enough work: about 30 ns per stored value.
  int nnz = tab.fpos[D];
  int nt  = threads_for(n_thr, ((double) nnz + D) * 30, 1e6);

  dc_tab   = &tab;
  dc_batch = INTEGER(sexp_batch);
  dc_neg   = INTEGER(sexp_neg);
  dc_logc  = isNull(sexp_logc) ? NULL : REAL(sexp_logc);
  dc_B     = B;
  dc_buf   = (double *) R_alloc((size_t) nt * ((size_t) n + 2 * (size_t) B + 2), sizeof(double));
  dc_freq  = REAL(sexp_freq);
  dc_prev  = INTEGER(sexp_prev);
  dc_k     = INTEGER(sexp_k);
  dc_ss0   = REAL(sexp_ss0);
  dc_ss1   = REAL(sexp_ss1);
  dc_pneg  = INTEGER(sexp_pneg);
  dc_ppos  = INTEGER(sexp_ppos);

  run_parallel(dc_worker, nt, D < 100 ? 100 : D);

  SEXP result = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(result, 0, sexp_freq);
  SET_VECTOR_ELT(result, 1, sexp_prev);
  SET_VECTOR_ELT(result, 2, sexp_k);
  SET_VECTOR_ELT(result, 3, sexp_ss0);
  SET_VECTOR_ELT(result, 4, sexp_ss1);
  SET_VECTOR_ELT(result, 5, sexp_pneg);
  SET_VECTOR_ELT(result, 6, sexp_ppos);

  SEXP names = PROTECT(allocVector(STRSXP, 7));
  const char *nm[7] = { "freq", "prev", "k", "ss0", "ss1", "prev_neg", "prev_pos" };
  for (int i = 0; i < 7; i++) SET_STRING_ELT(names, i, mkChar(nm[i]));
  setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(9);
  return result;
}
