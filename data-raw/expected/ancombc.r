# Writes tests/tinytests/expected/ancombc.rds: inputs and reference results
# for test-ancombc.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     python data-raw/expected/ancombc.py   # scikit-bio's outputs
#     Rscript data-raw/expected/ancombc.r   # folds them in, and runs ANCOMBC
#
# scikit-bio 0.7.4 (BSD-3, derived from ANCOMBC): ancombc() and ancombc2()
# on each case of data-raw/expected/ancombc/cases.csv, with their biases,
# sampling fractions and fitted quantities. Every number in those files is a
# hexadecimal float, read back exactly.
#
# ANCOMBC (Bioconductor, Artistic-2.0): ancombc() and ancombc2() on the same
# tables, where their settings can match: ancombc() always adds 1, so only
# the cases with `added` = 1; ancombc2() with prv_cut = 0, lib_cut = 0 and
# no sensitivity analysis. ANCOMBC 2.14.0 (the newest release) needs
# R >= 4.5; this was run with 2.8.1, the newest that installs on R 4.3.3,
# whose estimation code (R/ancombc_bias_correct.R: .iter_mle(), .bias_em())
# is the same as 2.14.0's. Set ANCOMBC_LIB to the library holding it.

dir <- file.path('data-raw', 'expected', 'ancombc')

hex <- function (x) as.numeric(x)
csv <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

tabs    <- csv('tables.csv')
cases   <- csv('cases.csv')
results <- csv('results.csv')
bias    <- csv('bias.csv')
fit     <- csv('fit.csv')
vcv     <- csv('vcov.csv')
dofs    <- csv('dof.csv')
sfs     <- csv('sf.csv')
pglob   <- csv('ph_global.csv')
ppairs  <- csv('ph_pairs.csv')
ptrend  <- csv('ph_trend.csv')
pdraws  <- csv('ph_draws.csv')
pbig    <- csv('ph_trend_big.csv')
struc   <- csv('struc.csv')

read_table <- function (nm) {
  m <- as.matrix(read.csv(
    file.path(dir, sprintf('table-%s.csv', nm)),
    header = FALSE, colClasses = 'character' ))
  m <- matrix(hex(m), nrow(m), ncol(m))
  if (all(m == round(m))) storage.mode(m) <- 'integer'
  dimnames(m) <- list(sprintf('S%02d', seq_len(nrow(m))), sprintf('F%02d', seq_len(ncol(m))))
  m
}

read_meta <- function (nm) {
  tr <- tabs[tabs$name == nm, ]
  m  <- read.csv(file.path(dir, sprintf('meta-%s.csv', nm)), colClasses = 'character')
  ty <- strsplit(strsplit(tr$types, ' ', fixed = TRUE)[[1]], ':', fixed = TRUE)
  lv <- if (nzchar(tr$levels)) strsplit(strsplit(tr$levels, ' ', fixed = TRUE)[[1]], '=', fixed = TRUE) else list()
  for (t in ty) {
    v <- t[[1]]
    m[[v]] <- switch(
      t[[2]],
      str  = m[[v]],
      num  = hex(m[[v]]),
      bool = as.logical(m[[v]]),
      cat  = factor(m[[v]], levels = strsplit(
        Filter(function (l) l[[1]] == v, lv)[[1]][[2]], '|', fixed = TRUE)[[1]]) )
  }
  rownames(m) <- sprintf('S%02d', seq_len(nrow(m)))
  m
}

# patsy's column names, as model.matrix() names them.
r_name <- function (x) {
  x <- gsub('[T.True]',  'TRUE',  x, fixed = TRUE)
  x <- gsub('[T.False]', 'FALSE', x, fixed = TRUE)
  x <- gsub('([A-Za-z_][A-Za-z0-9_.]*)\\[T\\.([^]]*)\\]', '\\1\\2', x)
  x <- gsub('([A-Za-z_][A-Za-z0-9_.]*)\\[([^]]*)\\]',     '\\1\\2', x)
  x[x == 'Intercept'] <- '(Intercept)'
  x
}

# A long table (feature, covariate, value) as a features x covariates matrix,
# covariates in patsy's order.
as_mat <- function (f, cov, v, D) {
  lev <- unique(cov)
  m   <- matrix(NA_real_, D, length(lev), dimnames = list(NULL, r_name(lev)))
  m[cbind(f, match(cov, lev))] <- v
  m
}


#--------------------------------------------------------
# scikit-bio.
#--------------------------------------------------------
tables <- lapply(setNames(tabs$name, tabs$name), function (nm)
  list(counts = read_table(nm), data = read_meta(nm)) )

ref <- lapply(setNames(cases$case, cases$case), function (cs) {

  cr <- cases[cases$case == cs, ]
  x  <- tables[[cr$table]]$counts
  D  <- ncol(x)
  r  <- results[results$case == cs, ]
  f  <- as.integer(r$feature)

  out <- list(
    table        = cr$table,
    method       = cr$method,
    formula      = cr$formula,
    groups       = if (nzchar(cr$grouping)) cr$grouping,
    pseudocount  = if (cr$method == 'ancombc') hex(cr$added) else hex(cr$pseudocount),
    var_quantile = hex(cr$var_quantile),
    p_adjust     = cr$p_adjust,
    alpha        = hex(cr$alpha),
    max_iters    = as.integer(cr$max_iter),
    tol          = hex(cr$tol),
    log_fc       = as_mat(f, r$covariate, hex(r$lfc), D),
    se           = as_mat(f, r$covariate, hex(r$se),  D),
    p_value      = as_mat(f, r$covariate, hex(r$p),   D),
    p_adjusted   = as_mat(f, r$covariate, hex(r$q),   D),
    reject       = as_mat(f, r$covariate, r$signif == '1', D) == 1 )

  b <- bias[bias$case == cs, ]
  out$bias <- cbind(
    delta_em  = hex(b$delta_em),
    delta_wls = hex(b$delta_wls),
    var_delta = hex(b$var_delta) )
  rownames(out$bias) <- colnames(out$log_fc)
  out$n_iters <- setNames(as.integer(b$n_iter), colnames(out$log_fc))

  s <- sfs[sfs$case == cs, ]
  if (nrow(s)) out$sample_fractions <- hex(s$theta)[order(as.integer(s$sample))]

  if (!is.null(out$groups)) {

    # The fitted quantities (_beta_hat, _var_hat, _vcov_hat), for a few
    # cases (FIT_CASES in ancombc.py).
    fr <- fit[fit$case == cs, ]
    if (nrow(fr)) {
      P  <- ncol(out$log_fc)
      out$coef <- matrix(hex(fr$beta), D, P, byrow = TRUE, dimnames = list(NULL, colnames(out$log_fc)))
      out$var  <- matrix(hex(fr$var),  D, P, byrow = TRUE, dimnames = list(NULL, colnames(out$log_fc)))
      vr <- vcv[vcv$case == cs, ]
      k  <- max(as.integer(vr$i)) + 1L
      out$vcov <- array(NA_real_, c(D, k, k))
      out$vcov[cbind(as.integer(vr$feature) + 1L, as.integer(vr$i) + 1L, as.integer(vr$j) + 1L)] <- hex(vr$value)
    }

    dr <- dofs[dofs$case == cs, ]
    dr <- dr[order(as.integer(dr$feature)), ]
    out$dof       <- if (cr$method == 'ancombc2') hex(dr$dof)
    out$estimable <- dr$estimable == '1'

    # Post-hoc tests.
    g <- pglob[pglob$case == cs, ]
    g <- g[order(as.integer(g$feature)), ]
    out$global <- data.frame(w = hex(g$w), p = hex(g$p), q = hex(g$q), reject = g$signif == '1')

    # Pairwise: every column (w is log_fc / se). Dunnett's: its log_fc, se
    # and w are the fit's (checked through the primary table), so only the
    # bootstrapped p-values, for the cases with recorded draws.
    pp <- ppairs[ppairs$case == cs & ppairs$test == 'pairwise', ]
    out$pairwise <- data.frame(
      log_fc = hex(pp$lfc), se = hex(pp$se), p = hex(pp$p), q = hex(pp$q),
      reject = pp$signif == '1' )
    drawn <- any(pdraws$case == cs)
    if (drawn) {
      pp <- ppairs[ppairs$case == cs & ppairs$test == 'dunnett', ]
      out$dunnett <- data.frame(p = hex(pp$p), q = hex(pp$q), reject = pp$signif == '1')
    }

    # Trend: the (deterministic) statistic, and with draws the p-values.
    tt <- ptrend[ptrend$case == cs, ]
    tt <- tt[order(as.integer(tt$feature)), ]
    out$trend <- if (drawn) data.frame(w = hex(tt$w), p = hex(tt$p), q = hex(tt$q), reject = tt$signif == '1')
                 else       data.frame(w = hex(tt$w))

    # The null draws, as features x k x bootstraps arrays (0 where none).
    k <- nrow(pp) %/% D   # Dunnett's rows: k comparisons per feature
    for (test in c('dunnett', 'trend')) {
      dd <- pdraws[pdraws$case == cs & pdraws$test == test, ]
      if (!nrow(dd)) next
      B  <- max(as.integer(dd$b)) + 1L
      a  <- array(0, c(D, k, B))
      a[cbind(as.integer(dd$feature) + 1L, as.integer(dd$j) + 1L, as.integer(dd$b) + 1L)] <- hex(dd$value)
      out[[paste0(test, '_draws')]] <- a
    }

    bg <- pbig[pbig$case == cs, ]
    if (nrow(bg)) out$trend_big <- hex(bg$p)[order(as.integer(bg$feature))]
  }

  out
})


# struc_zero(): one logical matrix per table, grouping and neg_lb.
struc_ref <- lapply(split(struc, paste(struc$table, struc$neg_lb)), function (z) {
  z   <- z[order(as.integer(z$feature)), ]
  lev <- strsplit(z$groups[[1]], ' ', fixed = TRUE)[[1]]
  m   <- do.call(rbind, lapply(strsplit(z$zero, ' ', fixed = TRUE), function (v) v == '1'))
  colnames(m) <- lev
  list(table = z$table[[1]], groups = z$grouping[[1]], neg_lb = z$neg_lb[[1]] == '1', zero = m)
})


#--------------------------------------------------------
# ANCOMBC.
#--------------------------------------------------------
lib <- Sys.getenv('ANCOMBC_LIB')
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ANCOMBC))
ancombc_version <- as.character(packageVersion('ANCOMBC'))

ancombc_ref <- list()

for (cs in cases$case) {

  cr <- cases[cases$case == cs, ]
  if (as.integer(cr$max_iter) != 100L || hex(cr$tol) != 1e-5) next
  if (grepl('0 +', cr$formula, fixed = TRUE)) next   # ANCOMBC prepends `~`

  x    <- tables[[cr$table]]$counts
  meta <- tables[[cr$table]]$data
  for (v in names(meta)) if (is.logical(meta[[v]])) meta[[v]] <- factor(meta[[v]], levels = c(FALSE, TRUE))
  meta$sample_id <- rownames(meta)   # ANCOMBC drops a one-column data frame
  grp  <- if (nzchar(cr$grouping)) cr$grouping

  if (cr$method == 'ancombc') {
    if (hex(cr$added) != 1) next
    res <- suppressWarnings(suppressMessages(ANCOMBC::ancombc(
      data = x, taxa_are_rows = FALSE, meta_data = meta, formula = cr$formula,
      p_adj_method = cr$p_adjust, prv_cut = 0, lib_cut = 0, group = grp,
      struc_zero = FALSE, tol = 1e-5, max_iter = 100, conserve = FALSE,
      alpha = hex(cr$alpha), global = !is.null(grp), n_cl = 1, verbose = FALSE )))
    rr  <- res$res
    get <- function (what) {
      m <- as.matrix(rr[[what]][, -1, drop = FALSE])
      rownames(m) <- rr[[what]]$taxon
      m
    }
    out <- list(
      log_fc = get('lfc'), se = get('se'), w = get('W'), p_value = get('p_val'),
      p_adjusted = get('q_val'), reject = get('diff_abn'),
      bias = as.matrix(res$delta_em), sample_fractions = res$samp_frac,
      global = if (!is.null(res$res_global)) res$res_global[c('W', 'p_val', 'diff_abn')] )
  }
  else {
    # ANCOMBC's %dorng%: the Dunnett and trend bootstraps. Only the trend
    # test's (deterministic) statistic is compared, so it gets 2 bootstraps:
    # each is a CVXR solve per feature and pattern, too slow for 100.
    set.seed(1)
    post <- !is.null(grp)
    k    <- if (post) nlevels(factor(meta[[grp]])) - 1L else 0L
    inc  <- if (post) { m <- diag(k); m[cbind(2:k, 1:(k - 1L))] <- -1; m }
    run2 <- function (post) suppressWarnings(suppressMessages(ANCOMBC::ancombc2(
      data = x, taxa_are_rows = FALSE, meta_data = meta,
      fix_formula = cr$formula, p_adj_method = cr$p_adjust,
      pseudo = hex(cr$pseudocount), pseudo_sens = FALSE, prv_cut = 0,
      lib_cut = 0, s0_perc = hex(cr$var_quantile), group = grp,
      struc_zero = FALSE, alpha = hex(cr$alpha), n_cl = 1, verbose = FALSE,
      global = post, pairwise = post, dunnet = post, trend = post,
      em_control = list(tol = 1e-5, max_iter = 100),
      mdfdr_control = list(fwer_ctrl_method = cr$p_adjust, B = 100),
      trend_control = if (post) list(
        contrast = list(inc, -inc), node = list(k, k), solver = 'ECOS', B = 2) )))
    # ANCOMBC's mdFDR calls p.adjust(n = Inf) when no feature passes its
    # global screen (R = 0), which fails for BY; then run without post-hoc
    # tests.
    res <- tryCatch(run2(post), error = function (e) {
      message(cs, ': ANCOMBC post-hoc tests failed (', conditionMessage(e), ')')
      post <<- FALSE
      run2(FALSE) })
    rr  <- res$res
    get <- function (prefix) {
      cols <- grep(paste0('^', prefix, '_'), names(rr), value = TRUE)
      m    <- as.matrix(rr[cols])
      colnames(m) <- sub(paste0('^', prefix, '_'), '', cols)
      rownames(m) <- rr$taxon
      m
    }
    out <- list(
      log_fc = get('lfc'), se = get('se'), w = get('W'), p_value = get('p'),
      p_adjusted = get('q'), reject = get('diff'),
      sample_fractions = res$samp_frac )
    if (post) {
      # Only the deterministic parts the tests compare.
      out$global   <- res$res_global[c('W', 'p_val', 'diff_abn')]
      out$pairwise <- res$res_pair[grep('^(lfc|q|diff)_', names(res$res_pair))]
      out$trend    <- res$res_trend['W']
    }
  }

  ancombc_ref[[cs]] <- out
}


#--------------------------------------------------------
# Save. The tests compare ANCOMBC only where its primary log fold changes
# agree with scikit-bio's to 1e-5: where both bias EMs end at the same
# estimate. Elsewhere one of them runs to max_iters without converging,
# and the estimate it stops at depends on the path of nloptr's
# minimisations rather than SciPy's; or the case differs by design
# (var_quantile = 0, interactions with zeros as missing, a badly scaled
# covariate). See plans/ancombc.md.
#--------------------------------------------------------
agree <- function (cs) {
  r <- ancombc_ref[[cs]]
  max(abs(r$log_fc[, colnames(ref[[cs]]$log_fc)] - ref[[cs]]$log_fc), na.rm = TRUE) < 1e-5
}
for (cs in names(ancombc_ref)) {
  r <- ancombc_ref[[cs]]
  cat(sprintf(
    '%-15s ANCOMBC vs scikit-bio: max |log_fc diff| %.1e, max rel se diff %.1e, reject differs %d\n',
    cs, max(abs(r$log_fc[, colnames(ref[[cs]]$log_fc)] - ref[[cs]]$log_fc), na.rm = TRUE),
    max(abs(r$se[, colnames(ref[[cs]]$se)] / ref[[cs]]$se - 1), na.rm = TRUE),
    sum(r$reject[, colnames(ref[[cs]]$reject)] != ref[[cs]]$reject, na.rm = TRUE) ))
}
kept <- Filter(agree, names(ancombc_ref))
cat('Compared by the tests:', kept, '\n')
ancombc_ref <- ancombc_ref[kept]

saveRDS(
  list(
    tables  = tables,
    skbio   = ref,
    struc   = struc_ref,
    ancombc = ancombc_ref,
    versions = c(skbio = '0.7.4', ANCOMBC = ancombc_version) ),
  file.path('tests', 'tinytests', 'expected', 'ancombc.rds'),
  version = 2, compress = 'xz' )
