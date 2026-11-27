# Black-box probes of zCompositions' cmultRepl() (GPL; output only: its
# source was never read), which pinned down the formulas behind
# replace_zeros_czm(), replace_zeros_gbm(), replace_zeros_sq() and
# replace_zeros_bl(). See plans/zero_replacement.md, "Confirmed formulas".
#
# Run from the package root, with zCompositions 1.6.2 installed:
#
#     Rscript data-raw/expected/zero_replacement_probe.r
#
# It prints the largest relative difference between cmultRepl() and each
# candidate formula; the confirmed candidates differ by ~1e-16.

library(zCompositions)
stopifnot(as.character(packageVersion('zCompositions')) == '1.6.2')

cr <- function (X, ...) {
  as.matrix(suppressMessages(cmultRepl(
    X, output = 'prop', z.delete = FALSE, z.warning = 1,
    suppress.print = TRUE, ... )))
}

# The Bayesian-multiplicative posterior, for a given prior t (N x D) and
# strength s (length N): zeros get s t / (n + s); non-zero parts are scaled
# to make room.
bm <- function (X, t, s) {
  n <- rowSums(X)
  x <- X / n
  r <- x
  for (i in seq_len(nrow(X))) {
    z <- X[i, ] == 0
    imp <- s[i] * t[i, z] / (n[i] + s[i])
    r[i, z]  <- imp
    r[i, !z] <- x[i, !z] * (1 - sum(imp))
  }
  r
}

# The confirmed candidate for every method, with adjust.
candidate <- function (X, method, frac = 0.65, threshold = 0.5, adjust = TRUE) {
  N <- nrow(X); D <- ncol(X); n <- rowSums(X); x <- X / n
  # Pooled leave-one-out: the other samples' counts over their depths.
  tl <- (matrix(colSums(X), N, D, byrow = TRUE) - X) / (sum(n) - n)
  s <- switch(method, BL = rep(D, N), SQ = sqrt(n), GBM = 1 / exp(rowMeans(log(tl))))
  imp <- if (method == 'CZM') matrix(frac * threshold / n, N, D) else s * tl / (n + s)
  if (adjust) {
    xm <- x; xm[X == 0] <- Inf
    M <- matrix(apply(xm, 2, min), N, D, byrow = TRUE)
    big <- X == 0 & imp > M
    imp[big] <- frac * M[big]
  }
  r <- x
  for (i in seq_len(N)) {
    z <- X[i, ] == 0
    r[i, z]  <- imp[i, z]
    r[i, !z] <- x[i, !z] * (1 - sum(imp[i, z]))
  }
  r
}

rel <- function (a, b) max(abs(a - b) / abs(b))


# 1. The posterior form, through method = "user". ====

X <- rbind(c(10, 0, 5, 3, 2), c(4, 6, 0, 0, 10), c(1, 20, 7, 8, 4), c(30, 2, 3, 0, 15))
N <- nrow(X); D <- ncol(X); n <- rowSums(X)
set.seed(1)
tt <- matrix(runif(N * D), N); tt <- tt / rowSums(tt)
ss <- runif(N, 1, 5)
cat('user t, s vs posterior form:', rel(cr(X, method = 'user', t = tt, s = ss, adjust = FALSE), bm(X, tt, ss)), '\n')


# 2. t and s, separating the candidates (unequal depths 20, 20, 40, 50). ====

x   <- X / n
uni <- matrix(1 / D, N, D)
loo_mean   <- t(sapply(1:N, function (i) colMeans(x[-i, , drop = FALSE])))
loo_pooled <- (matrix(colSums(X), N, D, byrow = TRUE) - X) / (sum(n) - n)
gm <- function (t) exp(rowMeans(log(t)))

for (cs in list(
    list('BL',  't = 1/D,        s = D',       uni,        rep(D, N)),
    list('BL',  't = loo mean,   s = D',       loo_mean,   rep(D, N)),
    list('BL',  't = loo pooled, s = D',       loo_pooled, rep(D, N)),
    list('SQ',  't = 1/D,        s = sqrt(n)', uni,        sqrt(n)),
    list('SQ',  't = loo pooled, s = sqrt(n)', loo_pooled, sqrt(n)),
    list('GBM', 't = loo mean,   s = 1/g',     loo_mean,   1 / gm(loo_mean)),
    list('GBM', 't = loo pooled, s = 1/g',     loo_pooled, 1 / gm(loo_pooled)) ))
  cat(sprintf('%-4s %-28s %.3g\n', cs[[1]], cs[[2]],
              rel(cr(X, method = cs[[1]], adjust = FALSE), bm(X, cs[[3]], cs[[4]]))))


# 3. The full candidate on random tables. ====

set.seed(2)
worst <- 0; ran <- 0; adjusted <- c(BL = 0, SQ = 0, GBM = 0, CZM = 0)
for (trial in 1:200) {
  N <- sample(3:12, 1); D <- sample(3:15, 1)
  X <- matrix(rpois(N * D, sample(c(1, 3, 20), 1)) * sample(c(1, 1, 1, 5, 50), N * D, TRUE), N, D)
  X[sample(N * D, N * D / 3)] <- 0
  if (any(rowSums(X > 0) < 1) || any(colSums(X > 0) < 2) || !any(X == 0)) next
  ran <- ran + 1
  for (m in names(adjusted)) {
    if (max(abs(cr(X, method = m) - cr(X, method = m, adjust = FALSE))) > 0)
      adjusted[m] <- adjusted[m] + 1
    for (adj in c(TRUE, FALSE))
      for (fr in c(0.65, 0.3))
        for (th in c(0.5, 0.7))
          worst <- max(worst, rel(
            cr(X, method = m, adjust = adj, frac = fr, threshold = th),
            candidate(X, m, frac = fr, threshold = th, adjust = adj) ))
  }
}
cat('random tables:', ran, ' largest relative difference:', worst, '\n')
cat('tables where adjust changed the result:\n'); print(adjusted)


# 4. The adjust test is strict: 1/16 against a minimum of 1/16. ====

Xt <- rbind(c(4, 0, 4), c(8, 1, 7), c(2, 3, 3))
cat('tie, imputed (0.0625 = not adjusted):', cr(Xt, method = 'CZM', frac = 0.5, threshold = 1)[1, 2], '\n')


# 5. Edge cases. ====

X <- rbind(c(10, 0, 5, 3, 2), c(4, 6, 0, 0, 10), c(1, 20, 7, 8, 4), c(30, 2, 3, 0, 15))
edge <- function (label, X)
  for (m in c('BL', 'SQ', 'GBM', 'CZM')) {
    a <- tryCatch(cr(X, method = m),
                  error   = function (e) paste('error:', conditionMessage(e)),
                  warning = function (w) paste('warning:', conditionMessage(w)) )
    cat(sprintf('%-22s %-4s %s\n', label, m,
      if (is.character(a)) gsub('\\s+', ' ', a)
      else sprintf('candidate differs by %.3g', rel(a, candidate(X, m))) ))
  }
X1 <- X; X1[, 2] <- c(0, 0, 5, 0); edge('one positive value', X1)
X2 <- X; X2[, 2] <- 0;             edge('part zero everywhere', X2)
edge('one sample', X[1, , drop = FALSE])
edge('two samples', X[1:2, ])
edge('not whole numbers', X / 7)
edge('no zeros', X + 1)
