#test_that("coda_lasso", {

# Inputs and reference results, written by data-raw/expected/coda_lasso.r
# (zeroSum 2.0.7, coda4microbiome 0.2.4, glmnet).
ref <- readRDS(file.path('expected', 'coda_lasso.rds'))

max_cpus <- if (at_home()) 4L else 2L

# Largest absolute difference, relative to the reference's scale.
close_to <- function (a, b, tol) {
  a <- as.vector(a); b <- as.vector(b)
  length(a) == length(b) && max(abs(a - b)) <= tol * max(1, abs(b))
}

# zeroSum's cost at (b0, beta): -loglik / n + the elastic net, with the
# ridge term over sd(y) (gaussian) and penalty weights v.
zs_cost <- function (X, y, b0, b, lambda, alpha, family, v = rep(1, length(b))) {
  eta <- drop(b0 + X %*% b)
  if (family == 'gaussian') {
    sdy <- sqrt(mean((y - mean(y))^2))
    mean((y - eta)^2) / 2 + lambda * ((1 - alpha) / (2 * sdy) * sum((v * b)^2) + alpha * sum(v * abs(b)))
  }
  else
    -mean(y * eta - log1p(exp(eta))) + lambda * ((1 - alpha) / 2 * sum((v * b)^2) + alpha * sum(v * abs(b)))
}

zs_args <- function (cs, cpus = 2L) {
  d  <- ref[[cs$data]]
  yy <- if (cs$family == 'binomial') d$yb else d$y
  args <- list(counts = d$x, y = yy, family = cs$family, l1_ratio = cs$alpha,
               standardize = cs$standardize, lambda = cs$lambda, folds = d$foldid,
               cpus = cpus)
  if (cs$given) { args$data <- d$z; args$given <- ~ age + sex }
  if (cs$pseudocount > 0) args$pseudocount <- cs$pseudocount
  args
}


# --- zeroSum: coefficients, costs, CV errors and the chosen lambda -----------

if (at_home()) {
  for (nm in names(ref$zs)) {
    cs  <- ref$zs[[nm]]
    fit <- do.call(coda_lasso, zs_args(cs))
    d   <- ref[[cs$data]]
    p   <- ncol(d$x)
    cp  <- attr(fit, 'coef_path')
    path <- attr(fit, 'path')
    zb  <- cs$coef[2:(p + 1), , drop = FALSE]

    # The coefficients agree to zeroSum's convergence (precision 1e-8 on the
    # cost): the solution is flat in some directions.
    expect_true(close_to(cp, zb, if (cs$family == 'binomial') 5e-3 else 2e-3), info = nm)
    expect_true(max(abs(colSums(cp))) < 1e-10, info = nm)
    expect_equal(path$lambda, cs$lambda, tolerance = 1e-14, info = nm)
    expect_true(close_to(path$cv_error, cs$cv[, 3], 5e-3), info = nm)
    expect_true(close_to(path$cv_se, cs$cv[, 4], 5e-3), info = nm)
    if (nm != 'big') expect_equal(path$n_selected, colSums(zb != 0), info = nm)
    expect_equal(which(path$lambda == attr(fit, 'lambda')), cs$lambda_min_index, info = nm)

    # Both reach the same cost, to far tighter than the coefficients agree.
    if (!cs$given) {
      X  <- log(d$x + cs$pseudocount)
      yy <- if (cs$family == 'binomial') d$yb else d$y
      v  <- if (cs$standardize) sqrt(colMeans(sweep(X, 2, colMeans(X))^2)) else rep(1, p)
      for (l in unique(c(1L, cs$lambda_min_index, length(cs$lambda)))) {
        f1 <- do.call(coda_lasso, utils::modifyList(zs_args(cs), list(lambda = cs$lambda[l], folds = NULL, n_folds = 0)))
        ce <- zs_cost(X, yy, attr(f1, 'intercept'), f1$coefficient, cs$lambda[l], cs$alpha, cs$family, v)
        cz <- zs_cost(X, yy, cs$coef[1, l], zb[, l], cs$lambda[l], cs$alpha, cs$family, v)
        expect_true(abs(ce - cz) < 1e-6 * abs(cz), info = paste(nm, l))
      }
    }
  }

  # '1se': zeroSum's rule (the lambda no smaller than the best whose CV error
  # is closest to the best plus its SD).
  cs  <- ref$zs$gaussian
  fit <- do.call(coda_lasso, c(zs_args(cs), list(select = '1se')))
  expect_equal(which(attr(fit, 'path')$lambda == attr(fit, 'lambda')), cs$lambda_1se_index)
  cs  <- ref$zs$gaussian_given
  fit <- do.call(coda_lasso, c(zs_args(cs), list(select = '1se')))
  expect_equal(which(attr(fit, 'path')$lambda == attr(fit, 'lambda')), cs$lambda_1se_index)

  # The covariates' coefficients (unpenalized, outside the constraint).
  for (nm in c('gaussian_given', 'binomial_given')) {
    cs  <- ref$zs[[nm]]
    fit <- do.call(coda_lasso, zs_args(cs))
    s   <- cs$lambda_min_index
    expect_true(close_to(attr(fit, 'covariates'), cs$coef[32:33, s], 5e-3), info = nm)
    expect_equal(names(attr(fit, 'covariates')), c('age', 'sexm'), info = nm)
    expect_true(close_to(attr(fit, 'intercept'), cs$coef[1, s], 5e-3), info = nm)
  }

  # zeroSum's default lambda path (n >= p: down to 1e-4).
  fit <- coda_lasso(ref$small$x, ref$small$y, folds = ref$small$foldid, cpus = 2)
  expect_equal(attr(fit, 'path')$lambda, ref$zs_path, tolerance = 1e-12)

  # Without the constraint, the same solver is glmnet's lasso.
  L   <- log(ref$small$x)
  gn  <- ref$glmnet_nozs
  out <- ecodive:::cl_zerosum(L, ref$small$y, NULL, ref$small$foldid, 5L, FALSE, 1, gn$lambda,
                              100L, NULL, FALSE, 'min', 1e-13, 100L, 2L, zero_sum = FALSE)
  expect_true(close_to(out$coef_path, gn$coef[-1, ], 1e-5))
  expect_true(close_to(out$intercept, gn$coef[1, ], 1e-5))
}


# --- coda4microbiome and glmnet ---------------------------------------------

if (at_home()) {
  for (nm in names(ref$c4m)) {
    cs <- ref$c4m[[nm]]
    if (cs$family == 'binomial') { x <- ref$c4m_bin$x; y <- factor(ref$c4m_bin$yb) }
    else                         { x <- ref$small$x;   y <- ref$small$y }
    fit  <- coda_lasso_pairwise(x, y, given = cs$covar, folds = cs$foldid, cpus = 2)
    path <- attr(fit, 'path')
    r    <- cs$result

    # glmnet's path (including where it stops early), CV curve and choice.
    expect_equal(path$lambda, cs$lambda, tolerance = 1e-12, info = nm)
    expect_true(close_to(path$cv_error, cs$cvm, 2e-2), info = nm)
    expect_true(close_to(path$cv_se, cs$cvsd, 2e-2), info = nm)
    expect_equal(attr(fit, 'lambda'), cs$lambda_1se, tolerance = 1e-12, info = nm)
    expect_equal(attr(fit, 'cv_measure'),
                 if (cs$family == 'binomial') 'auc' else 'mse', info = nm)

    # The coefficients: glmnet at thresh = 1e-14, and at its default 1e-7.
    cp <- attr(fit, 'coef_path')
    expect_true(close_to(cp, cs$beta_tight, 1e-2), info = nm)
    expect_true(close_to(cp, cs$beta, 1e-2), info = nm)
    expect_true(max(abs(colSums(cp))) < 1e-12, info = nm)

    # coda_glmnet()'s summary: taxa, rescaled coefficients, predictions.
    b   <- fit$coefficient
    sel <- which(b != 0)
    expect_identical(sel, r$taxa.num, info = nm)
    expect_true(close_to(2 * b[sel] / sum(abs(b[sel])), r$`log-contrast coefficients`, 2e-3), info = nm)
    expect_true(close_to(attr(fit, 'fitted')$link, r$predictions, 1e-2), info = nm)
    perf <- attr(fit, 'performance')
    if (cs$family == 'binomial') {
      expect_equal(perf[['auc']],       r$`apparent AUC`, tolerance = 1e-6, info = nm)
      expect_equal(perf[['cv_auc']],    r$`mean cv-AUC`,  tolerance = 1e-6, info = nm)
      expect_equal(perf[['cv_auc_se']], r$`sd cv-AUC`,    tolerance = 1e-6, info = nm)
    }
    else {
      expect_equal(perf[['r_squared']], r$`apparent Rsq`, tolerance = 1e-3, info = nm)
      expect_equal(perf[['cv_mse']],    r$`mean cv-MSE`,  tolerance = 2e-2, info = nm)
      expect_equal(perf[['cv_mse_se']], r$`sd cv-MSE`,    tolerance = 2e-2, info = nm)
    }
  }

  # The covariates' offset model is coda4microbiome's lm() / glm().
  cs  <- ref$c4m$gaussian_covar
  fit <- coda_lasso_pairwise(ref$small$x, ref$small$y, given = cs$covar, folds = cs$foldid, cpus = 1)
  expect_equal(unname(attr(fit, 'covariates')),
               unname(stats::coef(stats::lm(ref$small$y ~ age + sex, cs$covar))), tolerance = 1e-10)
}


# --- The two functions agree in one case -------------------------------------

# A pure lasso without standardization: pairwise at 2 lambda = constrained
# at lambda (the cheapest theta for a beta has |theta|_1 = |beta|_1 / 2).
# (zeroSum's tolerance is on the cost; the coefficients agree less closely.)
lcost <- function (f, l) {
  b <- f$coefficient
  mean((ref$small$y - attr(f, 'intercept') - log(ref$small$x) %*% b)^2) / 2 + l * sum(abs(b))
}
for (l in c(0.3, 0.1, 0.03)) {
  a <- coda_lasso(ref$small$x, ref$small$y, lambda = l, n_folds = 0, cpus = 1)
  b <- coda_lasso_pairwise(ref$small$x, ref$small$y, lambda = 2 * l, n_folds = 0, l1_ratio = 1,
                           standardize = FALSE, tol = 1e-14, cpus = 1)
  expect_true(abs(lcost(a, l) - lcost(b, l)) < 1e-7 * lcost(b, l), info = l)
  expect_true(close_to(b$coefficient, a$coefficient, 5e-3), info = l)
  expect_identical(b$selected, a$selected, info = l)
}


# --- Invariants --------------------------------------------------------------

x <- ref$small$x
y <- ref$small$y
yb <- factor(ifelse(ref$small$yb == 1, 'case', 'control'), c('control', 'case'))

# Scale invariance: counts or proportions, the same fit.
f1 <- coda_lasso(x, y, n_folds = 5, cpus = 1)
f2 <- coda_lasso(x / rowSums(x), y, n_folds = 5, cpus = 1)
expect_true(close_to(f1$coefficient, f2$coefficient, 1e-6))
expect_true(close_to(attr(f1, 'path')$lambda, attr(f2, 'path')$lambda, 1e-10))
expect_warning(f1 <- coda_lasso_pairwise(x, yb, n_folds = 5, cpus = 1), 'Fewer than 10 samples per fold')
f2 <- suppressWarnings(coda_lasso_pairwise(x * 1e3, yb, n_folds = 5, cpus = 1))
expect_true(close_to(f1$coefficient, f2$coefficient, 1e-6))
expect_true(all(abs(colSums(attr(f1, 'coef_path'))) < 1e-12))

# Identical for any number of threads.
for (fn in list(coda_lasso, coda_lasso_pairwise)) {
  r1 <- suppressWarnings(fn(x, yb, n_folds = 5, cpus = 1))
  for (cpus in 2:max_cpus) expect_identical(suppressWarnings(fn(x, yb, n_folds = 5, cpus = cpus)), r1)
  r1 <- fn(x, y, n_folds = 5, l1_ratio = 0.5, cpus = 1)
  for (cpus in 2:max_cpus) expect_identical(fn(x, y, n_folds = 5, l1_ratio = 0.5, cpus = cpus), r1)
}

# Sparse and dense input agree (and margin = 2).
xz <- x
xz[sample(length(xz), 200)] <- 0
slam_x <- structure(
  list(i = which(xz != 0, arr.ind = TRUE)[, 1], j = which(xz != 0, arr.ind = TRUE)[, 2],
       v = xz[xz != 0], nrow = nrow(xz), ncol = ncol(xz), dimnames = dimnames(xz)),
  class = 'simple_triplet_matrix' )
fd <- coda_lasso(xz, y, pseudocount = 0.5, n_folds = 5, cpus = 1)
expect_identical(coda_lasso(slam_x, y, pseudocount = 0.5, n_folds = 5, cpus = 1), fd)
expect_identical(coda_lasso(t(xz), y, margin = 2, pseudocount = 0.5, n_folds = 5, cpus = 1), fd)
if (requireNamespace('Matrix', quietly = TRUE)) {
  mx <- Matrix::Matrix(xz, sparse = TRUE)
  expect_identical(coda_lasso_pairwise(mx, y, pseudocount = 0.5, n_folds = 5, cpus = 1),
                   coda_lasso_pairwise(xz, y, pseudocount = 0.5, n_folds = 5, cpus = 1))
}
expect_warning(coda_lasso(xz, y, n_folds = 5, cpus = 1), 'pseudocount')

# The same seed, the same folds; another seed, other folds.
expect_identical(attr(coda_lasso(x, yb, cpus = 1), 'folds'), attr(coda_lasso(x, yb, cpus = 1), 'folds'))
expect_false(identical(attr(coda_lasso(x, yb, seed = 1, cpus = 1), 'folds'),
                       attr(coda_lasso(x, yb, cpus = 1), 'folds')))
# Binomial folds are stratified by class.
fo <- attr(coda_lasso(x, yb, n_folds = 4, cpus = 1), 'folds')
expect_true(all(abs(table(fo, yb)[, 1] - 5) <= 1))


# --- Results -------------------------------------------------------------------

fit <- coda_lasso(x, yb, n_folds = 5, cpus = 1)
expect_identical(class(fit), c('ecodive_coda_lasso', 'data.frame'))
expect_identical(names(fit), c('feature', 'coefficient', 'selected'))
expect_identical(fit$feature, colnames(x))
expect_identical(fit$selected, fit$coefficient != 0)
expect_identical(attr(fit, 'method'), 'coda_lasso')
expect_identical(attr(fit, 'family'), 'binomial')
expect_identical(attr(fit, 'select'), 'min')
expect_identical(attr(fit, 'y'), yb)
expect_identical(names(attr(fit, 'performance')), c('auc', 'cv_auc', 'cv_auc_se'))
expect_identical(names(attr(fit, 'fitted')), c('sample', 'link', 'probability'))
expect_identical(attr(fit, 'fitted')$sample, rownames(x))
expect_identical(dim(attr(fit, 'coef_path')), c(ncol(x), nrow(attr(fit, 'path'))))

# predict() reproduces the fitted values, for both functions.
expect_equal(unname(predict(fit, x)), attr(fit, 'fitted')$link, tolerance = 1e-12)
expect_equal(unname(predict(fit, x, type = 'response')), attr(fit, 'fitted')$probability,
             tolerance = 1e-12)
expect_equal(unname(predict(fit, x[, rev(colnames(x))])), attr(fit, 'fitted')$link, tolerance = 1e-12)
expect_equal(unname(predict(fit, unname(x))), attr(fit, 'fitted')$link, tolerance = 1e-12)
fg <- coda_lasso_pairwise(x, y, data = ref$small$z, given = ~ age + sex, n_folds = 5, cpus = 1)
expect_equal(unname(predict(fg, x, data = ref$small$z)), attr(fg, 'fitted')$link, tolerance = 1e-12)
expect_true(all(c('numerator', 'denominator', 'coefficient') == names(attr(fg, 'pairs'))))
expect_true(all(attr(fg, 'pairs')$coefficient != 0))
fg <- coda_lasso(x, y, data = ref$small$z, given = c('age', 'sex'), n_folds = 5, cpus = 1)
expect_equal(unname(predict(fg, x, data = ref$small$z)), attr(fg, 'fitted')$link, tolerance = 1e-12)
fg <- coda_lasso(x, y, given = ref$small$z, n_folds = 5, cpus = 1)
expect_equal(unname(predict(fg, x, given = ref$small$z)), attr(fg, 'fitted')$link, tolerance = 1e-12)
expect_error(predict(fg, x), '`given`')
expect_error(predict(fit, x[, 1:3]), 'lacks')

# `y` as a column of `data`, NA in y left out.
md <- data.frame(row.names = rownames(x), out = y)
md$out[3] <- NA
expect_warning(f1 <- coda_lasso(x, 'out', md, n_folds = 5, cpus = 1), 'Left out 1 sample')
expect_identical(nrow(attr(f1, 'fitted')), 39L)
expect_false(rownames(x)[3] %in% attr(f1, 'fitted')$sample)

# print() and the family from y.
expect_stdout(print(fit), 'Compositional Lasso')
expect_stdout(print(fg), 'Numerator')
expect_identical(attr(coda_lasso(x, y > median(y), n_folds = 5, cpus = 1), 'family'), 'binomial')
expect_identical(attr(coda_lasso(x, ref$small$yb, family = 'binomial', n_folds = 5, cpus = 1), 'family'),
                 'binomial')


# --- Edge cases and errors -----------------------------------------------------

expect_warning(fc <- coda_lasso(x, rep(2, 40), n_folds = 5, cpus = 1), 'constant')
expect_true(all(attr(fc, 'coef_path') == 0) && attr(fc, 'intercept') == 2)
expect_warning(fc <- coda_lasso_pairwise(x, rep(2, 40), n_folds = 5, cpus = 1), 'constant')
expect_true(all(fc$coefficient == 0))

# p = 2: the single log ratio.
f2 <- coda_lasso(x[, 1:2], y, n_folds = 5, cpus = 1)
expect_true(nrow(f2) == 2L && abs(sum(f2$coefficient)) < 1e-12)
f2 <- coda_lasso_pairwise(x[, 1:2], y, n_folds = 5, cpus = 1)
expect_true(nrow(f2) == 2L && abs(sum(f2$coefficient)) < 1e-12)

# A class smaller than n_folds.
y3 <- factor(c(rep('a', 36), rep('b', 4)))
expect_warning(coda_lasso(x, y3, n_folds = 5, cpus = 1), 'smallest class')

# Collinear covariates are dropped.
z2 <- data.frame(row.names = rownames(x), a = ref$small$z$age, b = 2 * ref$small$z$age)
expect_warning(fz <- coda_lasso(x, y, given = z2, n_folds = 5, cpus = 1), 'linear combination')
expect_identical(names(attr(fz, 'covariates')), 'a')

expect_error(coda_lasso(x, factor(rep(c('a', 'b', 'c'), length.out = 40)), cpus = 1), 'three|3 classes')
expect_error(coda_lasso(x, y, folds = factor(rep(1:2, 20), levels = 1:3), cpus = 1), 'Every fold')
expect_error(coda_lasso(x, y, n_folds = 0, cpus = 1), 'cross-validation')
expect_error(coda_lasso(x, y, l1_ratio = 2, cpus = 1), 'l1_ratio')
expect_error(coda_lasso(x, y, n_folds = 41, cpus = 1), 'n_folds')
expect_error(coda_lasso(x, cpus = 1), '`y` is required')
expect_error(coda_lasso(x, y, family = 'binomial', cpus = 1), 'binomial')
expect_error(coda_lasso(x, yb, family = 'gaussian', cpus = 1), 'numeric')
expect_error(coda_lasso(x[, 1, drop = FALSE], y, cpus = 1), 'two features')
expect_error(coda_lasso(x, y, select = 'max', cpus = 1))


# --- Plots ---------------------------------------------------------------------

pdf(NULL)
seed <- if (exists('.Random.seed', globalenv())) get('.Random.seed', globalenv()) else NULL
mfrow <- par('mfrow')
fa <- coda_lasso(x, y, n_folds = 5, cpus = 1)
expect_identical(withVisible(plot(fa))$visible, FALSE)
expect_identical(plot(fa), fa)
expect_identical(plot(fit, n_features = 2), fit)
expect_identical(plot(fa, type = 'path'), fa)
expect_identical(plot(fit, type = 'cv'), fit)
expect_identical(plot(fg), fg)
fn <- coda_lasso(x, y, lambda = 100, n_folds = 0, cpus = 1)
expect_true(all(fn$coefficient == 0))
expect_identical(plot(fn), fn)
expect_error(plot(fn, type = 'cv'), 'no cross-validation')
expect_error(plot(fn, type = 'path'), 'single lambda')
expect_error(plot(fa, n_features = 0), 'n_features')
expect_identical(par('mfrow'), mfrow)
expect_identical(if (exists('.Random.seed', globalenv())) get('.Random.seed', globalenv()) else NULL, seed)
invisible(dev.off())

#})
