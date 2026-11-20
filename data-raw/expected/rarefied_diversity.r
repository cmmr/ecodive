# Writes tests/tinytests/expected/rarefied_diversity.rds: vegan::avgdist()
# (vegan 2.7-8 or later, GPL-2, output-only) on a small table, for
# test-rarefied_diversity.r. avgdist() averages Bray-Curtis dissimilarities
# over `iterations` rarefactions drawn with R's RNG (rrarefy()), so
# ecodive's beta_div_rarefied() can agree with it only within Monte Carlo
# error: the test checks every cell within 4 combined standard errors. The
# same draws are averaged twice, with meanfun = mean and meanfun = sd, to
# give vegan's own standard deviations.
#
# Run from the package root:
#
#     Rscript data-raw/expected/rarefied_diversity.r

library(vegan)
stopifnot(packageVersion('vegan') >= '2.7.8')

set.seed(2026)
x <- matrix(
  rpois(8 * 20, rep(rexp(20, 1 / 6), each = 8)), 8, 20,
  dimnames = list(sprintf('S%d', 1:8), sprintf('F%02d', 1:20)) )
x[, 1] <- x[, 1] + 30L
depth  <- as.integer(min(rowSums(x)))
iterations <- 5000L

set.seed(1)
avg_mean <- avgdist(x, sample = depth, iterations = iterations, meanfun = mean)
set.seed(1)
avg_sd   <- avgdist(x, sample = depth, iterations = iterations, meanfun = sd)

fx <- list(
  counts     = x,
  depth      = depth,
  iterations = iterations,
  mean       = avg_mean,
  sd         = avg_sd,
  vegan      = as.character(packageVersion('vegan')) )

saveRDS(fx, 'tests/tinytests/expected/rarefied_diversity.rds', version = 2)
