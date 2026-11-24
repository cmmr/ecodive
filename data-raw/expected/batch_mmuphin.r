# Writes tests/tinytests/expected/batch_mmuphin.rds: inputs and reference
# results for test-batch_correct.r, so that the tests need no reference
# package.
#
# Run from the package root, with MMUPHin installed:
#
#     MMUPHIN_SRC=<unpacked MMUPHin 2.0.0 source> Rscript data-raw/expected/batch_mmuphin.r
#
# The reference is MMUPHin's adjust_batch() (MIT license). MMUPHin 2.0.0
# (Bioconductor 3.23) imports maaslin3, which needs R >= 4.4, so it can't be
# installed into R 4.3.3; the committed fixture comes from MMUPHin 1.16.0
# (Bioconductor 3.18), whose R/adjust_batch.R, R/helpers_adjust_batch.R and
# R/helpers.R compute the same (the files differ only in an author's e-mail
# address). With MMUPHIN_SRC set, the script also runs 2.0.0's own R files
# (sourced into an environment; MIT code) and checks that they give
# identical results.
#
# The internal estimates (gamma_hat, delta_hat, gamma_star, delta_star)
# come from MMUPHin's own helpers (construct_ind(), fit_stand_feature(),
# fit_EB(), fit_shrink(), adjust_EB(), back_transform_abd()), called in
# adjust_batch()'s order; the script checks that they reproduce
# adjust_batch()'s table exactly.
#
# Tables: a subset of MMUPHin's CRC_abd / CRC_meta (MIT; 5 studies, 15
# samples each, the 40 most prevalent species), and a simulated count
# table with planted batch shifts, features absent from some batches, and
# a covariate.

suppressPackageStartupMessages(library(MMUPHin))
cat('MMUPHin', format(packageVersion('MMUPHin')), '\n')

mm <- asNamespace('MMUPHin')

src <- Sys.getenv('MMUPHIN_SRC')
mm2 <- NULL
if (nzchar(src)) {
  mm2 <- new.env(parent = globalenv())
  for (f in c('helpers.R', 'helpers_adjust_batch.R', 'adjust_batch.R', 'controls.R'))
    sys.source(file.path(src, 'R', f), envir = mm2)
  cat('MMUPHin source', read.dcf(file.path(src, 'DESCRIPTION'), 'Version')[[1]], '\n')
}


# adjust_batch(), step by step, returning the table and the estimates.
# `feature_abd` is features x samples, as MMUPHin takes it.
run_mm <- function (ns, feature_abd, batch, covariates, data, zi, pc = NULL) {

  control <- list(
    zero_inflation = zi, pseudo_count = pc, verbose = FALSE,
    diagnostic_plot = NULL )
  ref <- ns$adjust_batch(feature_abd, batch, covariates, data, control)$feature_abd_adj

  type   <- ns$check_feature_abd(feature_abd)
  var_b  <- factor(data[[batch]])
  n_b    <- nlevels(var_b)
  bmod   <- ns$construct_design(data[, batch, drop = FALSE], with_intercept = FALSE)
  mod    <- ns$construct_design(if (is.null(covariates)) NULL else data[, covariates, drop = FALSE],
                                with_intercept = TRUE)[, -1, drop = FALSE]
  design <- cbind(bmod, mod)
  if (is.null(pc)) pc <- ns$set_pseudo(feature_abd)
  log_data <- ns$transform_features(
    ns$normalize_features(feature_abd, normalization = 'TSS', pseudo_count = pc),
    transform = 'LOG' )
  l_ind  <- ns$construct_ind(feature_abd, n_b, design, zi)
  stand  <- ns$fit_stand_feature(log_data, design, l_ind)
  fit    <- ns$fit_EB(stand$s_data, stand$l_stand_feature, bmod, n_b, l_ind)
  shrunk <- ns$fit_shrink(stand$s_data, fit, bmod, n_b, l_ind,
                          ns$match_control(ns$control_adjust_batch, control))
  adj    <- ns$adjust_EB(stand$s_data, shrunk, stand$l_stand_feature, bmod, n_b, l_ind)

  # Unrounded, then as adjust_batch() returns it.
  raw <- ns$back_transform_abd(adj, feature_abd, 'proportions')
  out <- ns$back_transform_abd(adj, feature_abd, type)
  stopifnot(identical(unname(out), unname(ref)))

  list(
    adjusted   = t(ref),          # samples x features
    unrounded  = if (type == 'counts') t(raw),
    gamma_hat  = unname(fit$gamma_hat),
    delta_hat  = unname(fit$delta_hat),
    gamma_star = unname(shrunk$gamma_star),
    delta_star = unname(shrunk$delta_star),
    ind_gamma  = unname(l_ind$ind_gamma),
    pseudocount = pc )
}


#--------------------------------------------------------
# Tables (samples x features).
#--------------------------------------------------------
data('CRC_abd', 'CRC_meta', package = 'MMUPHin', envir = environment())
studies <- sort(unique(CRC_meta$studyID))
keep_s  <- unlist(lapply(studies, function (s) head(which(CRC_meta$studyID == s), 15)))
prev    <- rowSums(CRC_abd[, keep_s] > 0)
keep_f  <- head(order(-prev, rownames(CRC_abd)), 40)
crc     <- t(CRC_abd[keep_f, keep_s])
crc_meta <- data.frame(
  studyID         = sub('\\..*', '', CRC_meta$studyID[keep_s]),
  study_condition = CRC_meta$study_condition[keep_s],
  row.names       = colnames(CRC_abd)[keep_s],
  stringsAsFactors = FALSE )
crc_meta$studyID <- factor(crc_meta$studyID)

set.seed(20261124)
n <- 48; D <- 30
sim_meta <- data.frame(
  batch = factor(rep(c('b1', 'b2', 'b3'), times = c(18, 14, 16))),
  cond  = rep(c('ctl', 'case'), length.out = n),
  age   = round(runif(n, 20, 70)),
  row.names = sprintf('S%02d', seq_len(n)),
  stringsAsFactors = FALSE )
mu <- matrix(rexp(D, 1 / 40), n, D, byrow = TRUE)
mu[sim_meta$batch == 'b2', 1:12]  <- mu[sim_meta$batch == 'b2', 1:12] * 3
mu[sim_meta$batch == 'b3', 8:16]  <- mu[sim_meta$batch == 'b3', 8:16] / 2
mu[sim_meta$cond == 'case', 18:24] <- mu[sim_meta$cond == 'case', 18:24] * 2
sim <- matrix(rnbinom(n * D, mu = mu, size = 1.5), n, D,
              dimnames = list(rownames(sim_meta), sprintf('F%02d', seq_len(D))))
sim[sim_meta$batch == 'b3', 26:27] <- 0     # absent from b3
sim[sim_meta$batch != 'b1', 28]    <- 0     # present in b1 only
sim[, 29] <- 0; sim[c(1, 20), 29] <- c(4, 7) # nearly absent


#--------------------------------------------------------
# Cases.
#--------------------------------------------------------
cases <- list(
  crc_zi        = list(table = 'crc', zi = TRUE,  cov = 'study_condition'),
  crc_nozi      = list(table = 'crc', zi = FALSE, cov = 'study_condition'),
  crc_zi_nocov  = list(table = 'crc', zi = TRUE,  cov = NULL),
  sim_zi        = list(table = 'sim', zi = TRUE,  cov = NULL),
  sim_zi_cov    = list(table = 'sim', zi = TRUE,  cov = c('cond', 'age')),
  sim_nozi      = list(table = 'sim', zi = FALSE, cov = NULL),
  sim_nozi_cov  = list(table = 'sim', zi = FALSE, cov = c('cond', 'age')),
  sim_zi_pc     = list(table = 'sim', zi = TRUE,  cov = 'cond', pc = 1) )

tables <- list(crc = crc, sim = sim)
metas  <- list(crc = crc_meta, sim = sim_meta)
batches <- c(crc = 'studyID', sim = 'batch')

results <- lapply(cases, function (cs) {
  x <- tables[[cs$table]]
  r <- run_mm(mm, t(x), batches[[cs$table]], cs$cov, metas[[cs$table]], cs$zi, cs$pc)
  if (!is.null(mm2)) {
    r2 <- run_mm(mm2, t(x), batches[[cs$table]], cs$cov, metas[[cs$table]], cs$zi, cs$pc)
    stopifnot(identical(r, r2))
  }
  r
})
if (!is.null(mm2)) cat('MMUPHin 2.0.0 source: identical results\n')

saveRDS(
  list(
    version = format(packageVersion('MMUPHin')),
    tables  = tables,
    metas   = metas,
    batches = batches,
    cases   = cases,
    results = results ),
  file.path('tests', 'tinytests', 'expected', 'batch_mmuphin.rds'),
  version = 2, compress = 'xz' )
