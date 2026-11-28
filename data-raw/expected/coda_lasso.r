# Writes tests/tinytests/expected/coda_lasso.rds: inputs and reference
# results for test-coda_lasso.r, so that the tests need no reference
# package.
#
# Run from the package root, with the references installed:
#
#     R_LIBS=<lib> Rscript data-raw/expected/coda_lasso.r
#
# References:
#   * zeroSum 2.0.7 (MIT; installed from its GitHub source,
#     cdn.jsdelivr.net/gh/rehbergT/zeroSum@2.0.7/zeroSum_2.0.7.tar.gz):
#     zeroSum() fits along a fixed lambda path with fixed folds
#     (cvStop = FALSE), and its default lambda path.
#   * coda4microbiome 0.2.4 (MIT): coda_glmnet(), with its cv.glmnet()
#     folds reproduced from the seed.
#   * glmnet (GPL-2; outputs only): cv.glmnet() on coda4microbiome's matrix
#     of pairwise log ratios, the same path with thresh = 1e-14, and
#     glmnet() on the log abundances (no zero-sum constraint).
#
# Fields of the .rds:
#
#   small     list(x (40 x 30 positive values, T01..T30), y (gaussian), yb
#             (0/1), z (data frame: age, sex), foldid (1..5))
#   big       list(x (200 x 500 counts, 48% zeros), y, foldid (1..5)); fitted
#             with pseudocount 0.5
#   zs        named list of zeroSum cases: list(data, family, alpha,
#             standardize, given, lambda, coef ((p+1) x L, intercept first),
#             cv (L x 5: lambda, training error, CV error, CV SD,
#             non-zero), lambda_min_index, lambda_1se_index)
#   zs_path   zeroSum's default lambda path for the small gaussian data
#   glmnet_nozs  list(lambda, coef ((p+1) x L)): glmnet() on log(x) with
#             standardize = FALSE, thresh = 1e-14 (zero-sum off)
#   c4m       named list of coda_glmnet() cases: list(data, family, covar,
#             foldid, result (coda_glmnet()'s list, without plots), lambda,
#             cvm, cvsd, lambda_min, lambda_1se, beta (p x L, glmnet's fit
#             as taxon coefficients), a0, beta_tight (the same with thresh =
#             1e-14), a0_tight)

suppressPackageStartupMessages({
  library(zeroSum)
  library(coda4microbiome)
  library(glmnet)
})
stopifnot(packageVersion('zeroSum') == '2.0.7',
          packageVersion('coda4microbiome') == '0.2.4')


# --- Inputs -----------------------------------------------------------------

set.seed(20261127)
n <- 40; p <- 30
x <- matrix(rgamma(n * p, 2), n, p) * rep(exp(rnorm(p, sd = 1.5)), each = n)
dimnames(x) <- list(sprintf('S%02d', 1:n), sprintf('T%02d', 1:p))
L <- log(x)
y  <- 1.5 * (L[, 'T03'] - L[, 'T07']) + 0.8 * (L[, 'T11'] - L[, 'T20']) + rnorm(n, sd = 0.7)
yb <- as.integer(y > stats::median(y))
z  <- data.frame(age = round(rnorm(n, 50, 10)), sex = sample(c('f', 'm'), n, TRUE),
                 row.names = rownames(x))
foldid <- sample(rep(1:5, length.out = n))
small <- list(x = x, y = y, yb = yb, z = z, foldid = foldid)

nb <- 200; pb <- 500
xb <- matrix(rnbinom(nb * pb, mu = 3, size = 0.3), nb, pb)
xb[, 1:4] <- xb[, 1:4] + rpois(nb * 4, 20)
dimnames(xb) <- list(sprintf('S%03d', 1:nb), sprintf('F%03d', 1:pb))
Lb <- log(xb + 0.5)
yB <- (Lb[, 1] - Lb[, 2]) + 0.5 * (Lb[, 3] - Lb[, 4]) + rnorm(nb, sd = 0.5)
big <- list(x = xb, y = yB, foldid = sample(rep(1:5, length.out = nb)))
cat('big zeros:', mean(xb == 0), '\n')


# --- zeroSum ----------------------------------------------------------------

# A path from the zero-sum lambda max (all pairs), log-spaced to 1%.
zs_lambdas <- function (L, yy, binom, alpha, std, Z = NULL, nl = 30) {
  n <- nrow(L)
  D <- if (is.null(Z)) matrix(1, n, 1) else cbind(1, Z)
  mu <- if (binom) stats::glm.fit(D, yy, family = binomial())$fitted.values else
    stats::lm.fit(D, yy)$fitted.values
  g <- drop(crossprod(L, (yy - mu) / n))
  v <- if (std) sqrt(colMeans(sweep(L, 2, colMeans(L))^2)) else rep(1, ncol(L))
  m <- 0
  for (j in 1:(ncol(L) - 1)) { k <- (j + 1):ncol(L); m <- max(m, abs(g[j] - g[k]) / (alpha * (v[j] + v[k]))) }
  exp(seq(log(m), log(m * 0.01), length.out = nl))
}

zs_case <- function (data, family, alpha, std, given = FALSE, nl = 30, pc = 0) {
  d  <- get(data)
  Lx <- log(d$x + pc)
  yy <- if (family == 'binomial') d$yb else d$y
  Z  <- NULL
  if (given) Z <- stats::model.matrix(~ age + sex, d$z)[, -1]
  lam <- zs_lambdas(Lx, yy, family == 'binomial', alpha, std, Z, nl)
  X <- if (given) cbind(Lx, Z) else Lx
  u <- c(rep(1, ncol(Lx)), rep(0, ncol(X) - ncol(Lx)))
  f <- zeroSum(X, yy, family = family, alpha = alpha, lambda = lam, foldid = d$foldid,
               standardize = std, cvStop = FALSE, threads = 1,
               zeroSum.weights = u, penalty.factor = u)
  coef <- sapply(f$coef, function (cc) as.numeric(cc))
  list(data = data, family = family, alpha = alpha, standardize = std, given = given,
       pseudocount = pc, lambda = lam, coef = coef,
       cv = unname(f$cv_stats[, 2:6]), lambda_min_index = f$lambdaMinIndex,
       lambda_1se_index = f$lambda1SEIndex)
}

zs <- list(
  gaussian       = zs_case('small', 'gaussian', 1,   FALSE),
  gaussian_enet  = zs_case('small', 'gaussian', 0.5, TRUE),
  binomial       = zs_case('small', 'binomial', 1,   FALSE),
  binomial_enet  = zs_case('small', 'binomial', 0.7, TRUE),
  gaussian_given = zs_case('small', 'gaussian', 1,   FALSE, given = TRUE),
  binomial_given = zs_case('small', 'binomial', 0.8, FALSE, given = TRUE) )
t0 <- Sys.time()
zs$big <- zs_case('big', 'gaussian', 1, FALSE, nl = 12, pc = 0.5)
cat('zeroSum big:', format(Sys.time() - t0), '\n')

# zeroSum's default path (its lambdaMax skips the first feature; check that
# the largest pair here doesn't involve it, so that the paths agree).
f0 <- zeroSum(L, y, foldid = foldid, cvStop = FALSE, threads = 1)
g0 <- drop(crossprod(L, (y - mean(y)) / n))
stopifnot(which.max(g0) != 1, which.min(g0) != 1)
zs_path <- f0$lambda


# --- glmnet without the constraint ------------------------------------------

lam_nozs <- exp(seq(log(max(abs(g0))), log(max(abs(g0)) * 0.01), length.out = 20))
gn <- glmnet(L, y, alpha = 1, standardize = FALSE, lambda = lam_nozs,
             control = list(thresh = 1e-14))
glmnet_nozs <- list(lambda = lam_nozs, coef = rbind(gn$a0, as.matrix(gn$beta)))


# --- coda4microbiome --------------------------------------------------------

taxon_beta <- function (B, idx, p) {
  apply(as.matrix(B), 2, function (th) {
    b <- numeric(p)
    for (c in seq_len(nrow(idx))) {
      b[idx[c, 1]] <- b[idx[c, 1]] + th[c]
      b[idx[c, 2]] <- b[idx[c, 2]] - th[c]
    }
    b })
}

c4m_case <- function (xx, yy, family, covar = NULL, seed) {
  lrm <- coda4microbiome:::logratios_matrix(xx)
  lrX <- lrm[[1]]; idx <- lrm[[2]]
  set.seed(seed)
  fid <- sample(rep(seq(10), length = nrow(xx)))
  set.seed(seed)
  r <- suppressWarnings(coda_glmnet(xx, yy, covar = covar, showPlots = FALSE))
  r$`predictions plot` <- NULL; r$`signature plot` <- NULL

  off <- NULL
  if (!is.null(covar)) {
    df0 <- data.frame(yy, covar)
    off <- if (family == 'binomial') predict(glm(yy ~ ., family = 'binomial', data = df0))
           else predict(lm(yy ~ ., data = df0))
  }
  tm <- if (family == 'binomial') 'auc' else 'deviance'
  cv <- suppressWarnings(cv.glmnet(lrX, yy, family = family, alpha = 0.9, type.measure = tm,
                                   foldid = fid, keep = TRUE, offset = off))
  ft <- glmnet(lrX, yy, family = family, alpha = 0.9, lambda = cv$lambda, offset = off,
               control = list(thresh = 1e-14, maxit = 1e7))
  list(family = family, covar = covar, foldid = fid, result = r,
       lambda = cv$lambda, cvm = cv$cvm, cvsd = cv$cvsd, measure = cv$name,
       lambda_min = cv$lambda.min, lambda_1se = cv$lambda.1se,
       beta = taxon_beta(cv$glmnet.fit$beta, idx, ncol(xx)), a0 = cv$glmnet.fit$a0,
       beta_tight = taxon_beta(ft$beta, idx, ncol(xx)), a0_tight = ft$a0)
}

# Binomial needs 10 samples per fold for the AUC: 100 samples, 15 taxa.
set.seed(7)
nB <- 100; pB <- 15
xB <- matrix(rgamma(nB * pB, 2), nB, pB) * rep(exp(rnorm(pB, sd = 1.5)), each = nB)
dimnames(xB) <- list(sprintf('S%03d', 1:nB), sprintf('T%02d', 1:pB))
LB <- log(xB)
ybB <- as.integer(1.2 * (LB[, 'T02'] - LB[, 'T05']) + 0.6 * (LB[, 'T09'] - LB[, 'T12']) +
                  rnorm(nB) > 0)
zB <- data.frame(age = round(rnorm(nB, 50, 10)), bmi = rnorm(nB, 25, 4),
                 row.names = rownames(xB))
c4m_bin <- list(x = xB, yb = ybB, z = zB)

c4m <- list(
  gaussian       = c4m_case(x, y, 'gaussian', seed = 11),
  gaussian_covar = c4m_case(x, y, 'gaussian', covar = z, seed = 12),
  binomial       = c4m_case(xB, factor(ybB), 'binomial', seed = 13),
  binomial_covar = c4m_case(xB, factor(ybB), 'binomial', covar = zB, seed = 14) )


saveRDS(
  list(small = small, big = big, c4m_bin = c4m_bin, zs = zs, zs_path = zs_path,
       glmnet_nozs = glmnet_nozs, c4m = c4m),
  'tests/tinytests/expected/coda_lasso.rds', version = 2 )
cat('size:', file.size('tests/tinytests/expected/coda_lasso.rds'), '\n')
