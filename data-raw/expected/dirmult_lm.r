# Writes tests/tinytests/expected/dirmult_lm.rds: inputs and reference
# results for test-dirmult_lm.r, so that the tests need no reference package.
#
# Run from the package root, with ecodive, ALDEx3 1.4.0 (CRAN, MIT; with
# lme4 and lmerTest for its mixed models) and ALDEx2 1.34.0 (Bioconductor
# 3.18, GPL >= 3: run for its outputs only) installed:
#
#     Rscript data-raw/expected/dirmult_lm.r
#
# * ALDEx3's aldex() on the same draws as ecodive: ALDEx3 draws its
#   compositions with rgamma(), samples then features, `nsample` at a time
#   (aldex.sampler(), rLogDirichlet()), so after set.seed() the same calls
#   regenerate them, and they are stored as the Gamma draws that
#   dirmult_lm_details() takes. Every aldex() call starts from the same
#   seed, so all of them draw the same compositions (checked). Scales:
#   ALDEx3's own clr.sm() draws (returned as logScale), given to ecodive as
#   `scale_draws`; and ecodive's own clr, tss and load scales (from
#   dirmult_lm_details()), given to ALDEx3 as a scale matrix. Two designs
#   (two groups; three groups and a numeric covariate) x three tests x
#   two p.adjust methods; and a random intercept (method = "lme4").
# * By distribution: ALDEx3 with its own 2000 draws (clr.sm, gamma 0.5) on
#   a deep table, and ALDEx2's aldex.glm() with its own 1000 clr draws, for
#   comparison with ecodive's own draws.

library(ecodive)
library(ALDEx3)   # ALDEx2 is used through ::, as both export aldex()

stopifnot(
  packageVersion('ALDEx3') == '1.4.0',
  packageVersion('ALDEx2') == '1.34.0' )

out <- file.path('tests', 'tinytests', 'expected', 'dirmult_lm.rds')


#--------------------------------------------------------
# A table (samples x features) and its metadata.
#--------------------------------------------------------
set.seed(20261124)
n   <- 15L
n_f <- 20L
meta <- data.frame(
  grp  = rep(c('a', 'b', 'c'), length.out = n),
  grp2 = rep(c('ctl', 'trt'), length.out = n),
  age  = round(runif(n, 20, 70)),
  id   = rep(sprintf('p%d', 1:5), each = 3),
  load = round(rlnorm(n, log(1e9), 0.6)),
  row.names = sprintf('S%02d', 1:n) )
mu  <- rlnorm(n_f, 3.5, 1.2)
fc  <- c(4, 0.25, 2, rep(1, n_f - 3))
lam <- outer(meta$grp2 == 'trt', log(fc)) + outer(meta$age - 45, rep(0.01, n_f))
counts <- matrix(
  rpois(n * n_f, exp(log(rep(mu, each = n)) + lam + rnorm(n * n_f, 0, 0.5))),
  n, n_f, dimnames = list(rownames(meta), sprintf('F%02d', 1:n_f)) )
counts[sample(length(counts), 25)] <- 0L
storage.mode(counts) <- 'integer'

n_draws <- 16L


# The Gamma draws behind aldex()'s compositions, after set.seed(seed).
regen <- function (x, S, seed) {
  set.seed(seed)
  G <- array(NA_real_, c(nrow(x), ncol(x), S))
  for (i in seq_len(nrow(x))) for (d in seq_len(ncol(x)))
    G[i, d, ] <- rgamma(S, x[i, d] + 0.5, 1)
  G
}
gamma_draws <- regen(counts, n_draws, 1)

aldex_on <- function (form, scale, test, adj, ...) {
  set.seed(1)
  r <- ALDEx3::aldex(
    t(counts), form, meta, nsample = n_draws, scale = scale, test = test,
    p.adjust.method = adj, n.cores = 1, ... )
  # The same compositions as gamma_draws.
  lc  <- aperm(r$logComp, c(2, 1, 3))
  lG  <- log(gamma_draws)
  lse <- apply(lG, c(1, 3), function (v) { m <- max(v); m + log(sum(exp(v - m))) })
  stopifnot(max(abs(lc - sweep(lG, c(1, 3), lse) / log(2))) < 1e-12)
  r
}

# The results, and with `draws` each draw's; the log scales when ALDEx3
# drew them (given = TRUE).
keep <- function (r, draws = FALSE, given = FALSE) {
  res <- list(
    estimate  = rowMeans(r$estimate, dims = 2),
    std.error = rowMeans(r$std.error, dims = 2),
    p.val     = r$p.val,
    p.val.adj = r$p.val.adj )
  if (given) res$logScale <- matrix(r$logScale, nrow(counts))
  if (draws) {
    res$draw_estimate  <- r$estimate
    res$draw_std.error <- r$std.error
  }
  res
}


#--------------------------------------------------------
# The cases on the same draws.
#--------------------------------------------------------
forms <- list(two = ~ grp2, three = ~ grp + age)
tests <- c(classical = 't', hc0 = 't.HC0', hc3 = 't.HC3')
adjs  <- c('BH', 'holm')

cases <- list()
add <- function (name, ...) cases[[name]] <<- list(...)

for (fn in names(forms)) for (se in names(tests)) for (adj in adjs) {

  form <- forms[[fn]]

  # ALDEx3's clr.sm() draws.
  r <- aldex_on(form, clr.sm, tests[[se]], adj, gamma = 0.5)
  add(
    paste(fn, se, adj, 'aldex_clr', sep = '/'), formula = form, se_type = se,
    p_adjust = adj, scale = 'given',
    ref = keep(r, draws = adj == 'BH' && se == 'hc3', given = TRUE) )

  # ecodive's own scale draws, given to ALDEx3.
  if (adj != 'BH') next
  for (sm in c('clr', 'tss', 'load')) {
    args <- list(
      counts, form, meta, n_draws = n_draws, se_type = se, p_adjust = adj,
      gamma_draws = gamma_draws, cpus = 1 )
    if (sm == 'load') { args$load <- 'load'; args$load_sd <- 0.4 }
    else              { args$scale_model <- sm; args$scale_sd <- 0.5 }
    e <- do.call(ecodive:::dirmult_lm_details, args)
    r <- aldex_on(form, e$draws$scale, tests[[se]], adj)
    add(
      paste(fn, se, adj, sm, sep = '/'), formula = form, se_type = se,
      p_adjust = adj, scale = sm, ref = keep(r) )
  }
}


#--------------------------------------------------------
# A random intercept: ALDEx3's method = "lme4", on
# ALDEx3's clr.sm() draws.
#--------------------------------------------------------
for (fn in c('two', 'three')) {
  form <- if (fn == 'two') ~ grp2 + (1 | id) else ~ grp + age + (1 | id)
  r <- aldex_on(form, clr.sm, 't', 'BH', gamma = 0.5, method = 'lme4')
  res <- keep(r, draws = TRUE, given = TRUE)
  # aldex() doesn't return the degrees of freedom: lmerTest's, for the
  # first four draws, fitted as ALDEx3's sr.mem() fits them.
  logW <- sweep(r$logComp, c(2, 3), r$logScale, `+`)
  res$draw_df <- array(NA_real_, c(nrow(r$X), n_f, 4L))
  for (d in 1:4) for (f in seq_len(n_f)) {
    dat   <- meta
    dat$y <- logW[f, , d]
    fit   <- suppressMessages(lmerTest::lmer(update(form, y ~ .), data = dat))
    res$draw_df[, f, d] <- coef(summary(fit))[, 'df']
  }
  add(
    paste(fn, 'lme4', sep = '/'), formula = form, se_type = 'classical',
    p_adjust = 'BH', scale = 'given', ref = res )
}


#--------------------------------------------------------
# By distribution: ALDEx3 and ALDEx2 with their own
# draws, on a deep table.
#--------------------------------------------------------
deep <- matrix(
  rpois(n * n_f, 2e4 * exp(log(rep(mu / mean(mu), each = n)) + lam + rnorm(n * n_f, 0, 0.5))),
  n, n_f, dimnames = dimnames(counts) )
storage.mode(deep) <- 'integer'

set.seed(2)
r3 <- ALDEx3::aldex(
  t(deep), ~ grp + age, meta, nsample = 2000, scale = clr.sm, gamma = 0.5,
  n.cores = 1, return.pars = c('X', 'estimate', 'std.error', 'p.val', 'p.val.adj') )
aldex3_own <- list(
  estimate    = rowMeans(r3$estimate, dims = 2),
  estimate_sd = apply(r3$estimate, c(1, 2), stats::sd),
  std.error   = rowMeans(r3$std.error, dims = 2),
  p.val       = r3$p.val,
  p.val.adj   = r3$p.val.adj )

set.seed(3)
mm  <- model.matrix(~ grp + age, meta)
clr <- ALDEx2::aldex.clr(t(deep), mm, mc.samples = 1000, denom = 'all', verbose = FALSE)
glm <- ALDEx2::aldex.glm(clr, mm)
aldex2_glm <- as.matrix(glm)


saveRDS(
  list(
    counts = counts, meta = meta, n_draws = n_draws, gamma_draws = gamma_draws,
    cases = cases, deep = deep, aldex3_own = aldex3_own, aldex2_glm = aldex2_glm,
    versions = c(
      ALDEx3   = as.character(packageVersion('ALDEx3')),
      ALDEx2   = as.character(packageVersion('ALDEx2')),
      lme4     = as.character(packageVersion('lme4')),
      lmerTest = as.character(packageVersion('lmerTest')),
      R        = as.character(getRversion()) )),
  out, version = 2 )
