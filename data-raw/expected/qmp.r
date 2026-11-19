# Writes tests/tinytests/expected/qmp.rds: outputs of the authors' QMP script
# for test-qmp.r, so that the tests need no reference package.
#
# Run from the package root, with ecodive and phyloseq installed:
#
#     Rscript data-raw/expected/qmp.r
#
# data-raw/expected/qmp/QMP.R is rarefy_even_sampling_depth() from
# https://github.com/raeslab/QMP (master; BSD-3, (c) 2017 Raes Lab, see
# qmp/LICENSE). It rarefies with phyloseq::rarefy_even_depth() (phyloseq is
# AGPL: output only). The script is run unchanged, except that
# rarefy_even_depth() is wrapped to record the sample size it is asked for
# and the one it gives, and, for the mean over many draws, to vary its
# hard-coded rngseed. Rarefaction is random, so the tests compare:
#
#   * exactly: the sampling depths, each sample's number of reads after
#     rarefying (the targets), and the scaling;
#   * statistically: the mean over many seeds, which for both should be
#     x / n * load.

suppressPackageStartupMessages({ library(phyloseq); library(ecodive) })

# Only the function from the script (not its example call).
exprs <- parse(file.path('data-raw', 'expected', 'qmp', 'QMP.R'))
for (e in exprs)
  if (is.call(e) && as.character(e[[1]]) %in% c('<-', '=') &&
      identical(e[[2]], as.name('rarefy_even_sampling_depth')))
    eval(e, globalenv())

# Record what rarefy_even_depth() is asked for and what it gives. The
# script's function looks it up in the global environment first.
calls <- NULL
seed_override <- NULL
substitute_with <- NULL
rarefy_even_depth <- function (physeq, sample.size, rngseed, ...) {
  if (!is.null(seed_override)) rngseed <- seed_override
  r <- phyloseq::rarefy_even_depth(physeq, sample.size = sample.size, rngseed = rngseed, ...)
  if (!is.null(substitute_with)) r[] <- substitute_with[nrow(calls) %||% 0 + 1, ]
  calls <<- rbind(calls, data.frame(asked = sample.size, got = sum(r)))
  r
}
`%||%` <- function (a, b) if (is.null(a)) b else a

run <- function (x, load, ...) {
  calls <<- NULL
  q <- rarefy_even_sampling_depth(x, matrix(load, ncol = 1, dimnames = list(rownames(x), '#')))
  list(qmp = unclass(q), asked = calls$asked, got = calls$got)
}

cases <- list()

# 1. The script's own example.
a <- matrix(
  c(4,4,2,1,8,5,2,0,3,5,3,1,10,8,3,0,0,6,4,3), nrow = 5, ncol = 4, byrow = TRUE,
  dimnames = list(c('Sample A', 'Sample B', 'Sample C', 'Sample D', 'Sample E'),
                  c('taxa1', 'taxa2', 'taxa3', 'taxa4')) )
b <- c(10, 20, 34, 21, 12) * 100000
cases$example <- c(list(x = a, load = b), run(a, b))

# 2. A random table: non-integer targets, a few deep cells.
set.seed(2017)
x <- matrix(rpois(20 * 30, rexp(30, 1 / 6)), 20, 30,
  dimnames = list(sprintf('S%02d', 1:20), sprintf('T%02d', 1:30)))
x[sample(length(x), 150)] <- 0
x[1, 1] <- 400
L <- signif(exp(rnorm(20, 23, 0.8)), 4)
cases$random <- c(list(x = x, load = L), run(x, L))

# 3. A table where (n / L) * L for the smallest-depth sample comes out below
#    n in floating point, so a plain floor (the script's) loses a read.
for (k in 1:5000) { n1 <- 10 + k; L1 <- 0.1 * k + 0.07; if (n1 / L1 * L1 < n1) break }
x3 <- matrix(c(n1 - 5, 5, 0, 2 * n1 - 3, 3, 0, 4 * n1, 0, 7), 3, byrow = TRUE,
  dimnames = list(c('P', 'Q', 'R'), c('a', 'b', 'c')))
L3 <- c(L1, 2 * L1, 3 * L1)
cases$hair <- c(list(x = x3, load = L3), run(x3, L3))

# Mean over many draws (the script with rngseed = 1, 2, ...), on a small
# table, for a statistical comparison.
xs <- x[1:6, 1:8]
Ls <- L[1:6]
K  <- 400
acc <- 0
for (s in seq_len(K)) {
  seed_override <- s
  acc <- acc + run(xs, Ls)$qmp
}
seed_override <- NULL
mean_draws <- list(x = xs, load = Ls, K = K, mean = acc / K)

# With ecodive's rarefied counts substituted for phyloseq's, the script's
# output equals qmp_rarefy()'s (only the random draws differ).
q  <- qmp_rarefy(x, L, seed = 1)
tg <- attr(q, 'target')
r  <- as.matrix(rarefy(x, depth = tg, seed = 1, drop = FALSE, warn = FALSE))
calls <- NULL
substitute_with <- r
qs <- rarefy_even_sampling_depth(x, matrix(L, ncol = 1, dimnames = list(rownames(x), '#')))
substitute_with <- NULL
qq <- unclass(q); attributes(qq) <- attributes(qq)[c('dim', 'dimnames')]
cat('Script with ecodive\'s draws vs qmp_rarefy(): max relative difference',
    max(abs(qs - qq) / pmax(abs(qq), 1e-300)), '\n')
for (nm in names(cases))
  cat(nm, ': targets asked', paste(format(cases[[nm]]$asked), collapse = ' '),
      '\n   got', paste(cases[[nm]]$got, collapse = ' '), '\n')

saveRDS(
  list(cases = cases, mean_draws = mean_draws),
  file.path('tests', 'tinytests', 'expected', 'qmp.rds'),
  version = 2 )
