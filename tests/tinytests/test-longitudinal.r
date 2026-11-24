# first_differences(), first_distances(), volatility() and nmit(), against
# q2-longitudinal 2026.10.0's own functions (fixtures from
# data-raw/expected/longitudinal.py) and direct R implementations.

E    <- readRDS(file.path('expected', 'longitudinal.rds'))
meta <- E$meta
rownames(meta) <- meta$sample
cnt  <- E$counts


# Rows of result `res` against q2's (sample, value) pairs `ref`: the same
# samples, and the largest difference of column `col`.
by_sample <- function (res, ref, col) {
  m <- match(ref$sample, res$sample)
  list(same = nrow(res) == nrow(ref) && !anyNA(m), diff = max(abs(res[[col]][m] - ref$value)))
}



#--------------------------------------------------
# first_differences(): q2's pairs and values, exactly.
#--------------------------------------------------
fd <- first_differences('shannon', 'subject', 'state', meta, replicates = 'drop')
r  <- by_sample(fd, E$fd$drop, 'difference')
expect_true(r$same)
expect_identical(r$diff, 0)
expect_inherits(fd, c('ecodive_first_differences', 'data.frame'))
expect_identical(names(fd), c('sample', 'subject', 'state_from', 'state', 'difference'))
expect_null(attr(fd, 'baseline'))

for (b in c('baseline', 'baseline7')) {
  bv <- if (b == 'baseline') 0 else 7
  x  <- first_differences('shannon', 'subject', 'state', meta, replicates = 'drop', baseline = bv)
  r  <- by_sample(x, E$fd[[b]], 'difference')
  expect_true(r$same, info = b)
  expect_identical(r$diff, 0, info = b)
  expect_true(all(x$state_from == bv), info = b)
  expect_identical(attr(x, 'baseline'), bv, info = b)
}

# Every feature of the count table at once, as relative abundances.
ft <- first_differences(cnt, 'subject', 'state', meta, replicates = 'drop', norm = 'percent')
expect_identical(names(ft), c('sample', 'subject', 'state_from', 'state', colnames(cnt)))
for (f in names(E$fd_table)) {
  r <- by_sample(ft, E$fd_table[[f]], f)
  expect_true(r$same, info = f)
  expect_true(r$diff < 1e-15, info = f)
}

# Ordered by pair of states, then subject; consecutive states only.
expect_false(is.unsorted(fd$state))
for (s in unique(fd$state))
  expect_false(is.unsorted(fd$subject[fd$state == s]), info = s)
st <- sort(unique(meta$state))
expect_identical(fd$state_from, st[match(fd$state, st) - 1L])

# The values: the later sample's minus the earlier one's.
expect_equal(fd$difference, meta[fd$sample, 'shannon'] -
  meta$shannon[match(paste(fd$subject, fd$state_from), paste(meta$subject, meta$state))])

# q2's skipped state: S03 has no day 7, so no (3, 7) or (7, 14) pair; S09
# (no day 0) none for (0, 3); S14 (one sample) none at all; the NA value
# (S08 at day 3) drops both of its pairs.
expect_identical(fd$state[fd$subject == 'S03'], c(3, 28))
expect_identical(fd$state[fd$subject == 'S09'], c(7, 14, 28))
expect_false('S14' %in% fd$subject)
expect_identical(fd$state[fd$subject == 'S08'], c(14, 28))
expect_identical(fd$state[fd$subject == 'S05'], c(3, 7))   # replicates dropped

# A vector with names, and in another order: the same pairs.
v  <- setNames(meta$shannon, meta$sample)
o  <- rev(seq_len(nrow(meta)))
x  <- first_differences(v[o], meta$subject[o], meta$state[o], replicates = 'drop')
expect_identical(x, fd)
x  <- first_differences(v[o], 'subject', 'state', meta, replicates = 'drop')
expect_identical(x, fd)

# Without names: the samples' positions.
x <- first_differences(unname(v), meta$subject, meta$state, replicates = 'drop')
expect_identical(x$sample, match(fd$sample, meta$sample))
expect_identical(x$difference, fd$difference)

# Count tables: margin, sparse types, other norms.
expect_identical(first_differences(t(cnt), 'subject', 'state', meta, replicates = 'drop', margin = 2L, norm = 'percent'), ft)
if (requireNamespace('Matrix', quietly = TRUE))
  expect_identical(first_differences(as(cnt * 1, 'CsparseMatrix'), 'subject', 'state', meta, replicates = 'drop', norm = 'percent'), ft)
if (requireNamespace('slam', quietly = TRUE))
  expect_identical(first_differences(slam::as.simple_triplet_matrix(cnt), 'subject', 'state', meta, replicates = 'drop', norm = 'percent'), ft)
x <- first_differences(cnt, 'subject', 'state', meta, replicates = 'drop')
expect_equal(x$F01, cnt[x$sample, 'F01'] -
  cnt[meta$sample[match(paste(x$subject, x$state_from), paste(meta$subject, meta$state))], 'F01'],
  check.attributes = FALSE)
x <- first_differences(cnt, 'subject', 'state', meta, replicates = 'drop', norm = 'clr', pseudocount = 1)
expect_true(all(is.finite(as.matrix(x[, colnames(cnt)]))))



#--------------------------------------------------
# Replicates.
#--------------------------------------------------
expect_error(
  first_differences('shannon', 'subject', 'state', meta),
  'Detected replicate samples for subject S05 at state 14 (2 samples)', fixed = TRUE )
expect_true(grepl('Detected replicate samples for individual (subject) S05 at state (state) 14.0', E$errors[['replicate']], fixed = TRUE))
expect_error(first_differences('shannon', 'subject', 'state', meta, replicates = 'x'), "`replicates` must be one of 'error', 'drop', 'random'")

reps <- meta$sample[meta$subject == 'S05' & meta$state == 14]
x1 <- first_differences('shannon', 'subject', 'state', meta, replicates = 'random')
expect_identical(x1, first_differences('shannon', 'subject', 'state', meta, replicates = 'r'))
expect_identical(nrow(x1), nrow(fd) + 2L)
i <- x1$subject == 'S05' & x1$state == 14
expect_true(sum(i) == 1 && x1$sample[i] %in% reps)
# The (14, 28) pair's earlier sample is a fresh pick, as q2's: different
# seeds pick both replicates for each side.
picks <- vapply(0:15, function (s) {
  y <- first_differences('shannon', 'subject', 'state', meta, replicates = 'random', seed = s)
  y$sample[y$subject == 'S05' & y$state == 14]
}, character(1))
expect_true(all(reps %in% picks))
d1 <- vapply(0:15, function (s) {
  y <- first_differences('shannon', 'subject', 'state', meta, replicates = 'random', seed = s)
  y$difference[y$subject == 'S05' & y$state == 28]
}, numeric(1))
expect_identical(length(unique(d1)), 2L)
# Otherwise the same as 'drop'.
expect_identical(as.list(x1[x1$subject != 'S05', ]), as.list(fd[fd$subject != 'S05', ]))



#--------------------------------------------------
# first_distances(): q2's distances, exactly.
#--------------------------------------------------
dm <- E$dm
fdd <- first_distances(dm, 'subject', 'state', meta, replicates = 'drop')
r   <- by_sample(fdd, E$fdist$drop, 'distance')
expect_true(r$same)
expect_identical(r$diff, 0)
expect_inherits(fdd, c('ecodive_first_distances', 'data.frame'))
expect_identical(names(fdd), c('sample', 'subject', 'state_from', 'state', 'distance'))
x <- first_distances(dm, 'subject', 'state', meta, replicates = 'drop', baseline = 0)
r <- by_sample(x, E$fdist$baseline, 'distance')
expect_true(r$same)
expect_identical(r$diff, 0)

# The same as indexing the full matrix.
m    <- as.matrix(dm)
from <- meta$sample[match(paste(fdd$subject, fdd$state_from), paste(meta$subject, meta$state))]
expect_identical(fdd$distance, unname(m[cbind(fdd$sample, from)]))
expect_identical(first_distances(m, 'subject', 'state', meta, replicates = 'drop'), fdd)

# Same pairs as first_differences().
# Same pairs as first_differences(), and the two that the NA value of
# shannon (S08 at day 3) drops there.
expect_true(all(fd$sample %in% fdd$sample))
extra <- fdd[!fdd$sample %in% fd$sample, ]
expect_identical(extra$subject, c('S08', 'S08'))
expect_identical(extra$state, c(3, 7))

# The title names the input's.
bc <- bray(cnt)
x  <- first_distances(bc, 'subject', 'state', meta, replicates = 'drop')
expect_identical(attr(x, 'display')$title, 'First Distances of Bray-Curtis Dissimilarity')
expect_identical(attr(fdd, 'display')$title, 'First Distances')
x <- first_distances(bray(cnt, norm = 'percent'), 'subject', 'state', meta, replicates = 'drop')
expect_true(max(abs(x$distance - fdd$distance)) < 1e-15)



#--------------------------------------------------
# Input checks (q2's errors where it has them).
#--------------------------------------------------
expect_true(grepl('state_column must contain at least two unique values', E$errors[['validate_one_state']]))
one <- meta[meta$state == 3, ]
expect_error(first_differences('shannon', 'subject', 'state', one), 'at least two distinct values')
expect_error(first_differences('shannon', 'subject', 'state', meta, replicates = 'drop', baseline = 5),
             '`baseline` must be one of the states: 5 is not in 0, 3, 7, 14, 28.', fixed = TRUE)
expect_true(grepl('baseline must be a valid state: 5.0 is not in', E$errors[['baseline']], fixed = TRUE))
expect_error(first_differences('shannon', 'subject', 'subject', meta), 'different columns')
expect_error(first_differences('shannon', 'subject', 'group', meta), '`states` must be numeric')
expect_error(first_differences('shannon', 'subject', 'nope', meta), "'nope' is not one")
expect_error(first_differences('nope', 'subject', 'state', meta), "'nope' is not one")
expect_error(first_differences('group', 'subject', 'state', meta), 'must be numeric')
expect_error(first_differences('shannon', state = 'state', data = meta), '`subjects` is required')
expect_error(first_differences('shannon', 'subject', data = meta), '`states` is required')
bad <- meta; bad$subject[3] <- NA
expect_error(first_differences('shannon', 'subject', 'state', bad), "Missing values (NA) in `data` column 'subject'", fixed = TRUE)
expect_error(first_differences(v, meta$subject, c(NA, meta$state[-1])), '`states` must not contain NAs')
expect_error(first_differences(v, meta$subject, meta$state[-1]), '`states` has 61 values')
expect_error(first_differences(letters, 1, 1), 'numeric vector')
expect_error(first_differences(c(1, 2), c('a', 'b'), c(1, 2)), 'No pairs of samples')
expect_error(first_differences(c(1, NA), c('a', 'a'), c(1, 2)), 'No pairs of samples')
expect_error(first_distances(dm, 'subject', 'state', meta[1:10, ]), 'have no row in `data`')
expect_error(first_distances(dist(1), 1, 1), 'at least 2 samples')
expect_error(volatility(cnt, 'subject', 'state', data = meta), 'numeric vector with one value per sample')

# Two samples of one subject.
x <- first_differences(c(1, 3), c('a', 'a'), c(5, 2))
expect_identical(x, structure(
  data.frame(sample = 1L, subject = 'a', state_from = 2, state = 5, difference = -2, stringsAsFactors = FALSE),
  class = c('ecodive_first_differences', 'data.frame') ))



#--------------------------------------------------
# nmit(): q2's distances.
#--------------------------------------------------
for (mth in c('kendall', 'pearson', 'spearman')) for (mn in c('frobenius', 'nuclear')) {
  key <- paste(mth, if (mn == 'frobenius') 'fro' else 'nuc', sep = '_')
  ref <- E$nmit[[key]]
  lab <- unname(E$nmit_ids[rownames(ref)])
  dimnames(ref) <- list(lab, lab)
  d <- suppressWarnings(nmit(cnt, 'subject', meta, method = mth, matrix_norm = mn))
  expect_identical(attr(d, 'Labels'), sort(lab), info = key)
  got <- as.matrix(d)[lab, lab]
  expect_true(max(abs(got - ref) / pmax(ref, 1)) < 1e-12, info = key)
}

d <- suppressWarnings(nmit(cnt, 'subject', meta))
expect_inherits(d, c('ecodive_dist', 'dist'))
expect_identical(attr(d, 'method'), 'nmit')
expect_identical(attr(d, 'display')$title, 'NMIT Distance')
expect_identical(attr(d, 'Size'), 14L)
expect_true(any(grepl('14 subjects', capture.output(print(d)))))

expect_warning(nmit(cnt, 'subject', meta), "2 subjects have fewer than 3 samples, so their correlations are only +/-1 or 0: 'S13' (2), 'S14' (1).", fixed = TRUE)

# Against R's cor(), norm() and svd(), on the subjects with 3 or more samples.
keep <- meta$sample[!meta$subject %in% c('S13', 'S14')]
rel  <- cnt[keep, ] / rowSums(cnt[keep, ])
sub  <- meta[keep, 'subject']
for (mth in c('kendall', 'pearson', 'spearman')) {
  C <- lapply(split(seq_along(sub), sub), function (i) {
    r <- suppressWarnings(stats::cor(rel[i, ], method = mth))
    # pandas: Kendall's diagonal is 1; Pearson's and Spearman's is
    # undefined (0) for a constant feature, where cor() gives 1.
    if (mth == 'kendall') diag(r) <- 1
    else diag(r)[apply(rel[i, ], 2, function (v) all(v == v[[1]]))] <- 0
    r[is.na(r)] <- 0
    r
  })
  S  <- names(C)
  fro <- nuc <- matrix(0, length(S), length(S), dimnames = list(S, S))
  for (a in S) for (b in S) {
    fro[a, b] <- norm(C[[a]] - C[[b]], 'F')
    nuc[a, b] <- sum(svd(C[[a]] - C[[b]])$d)
  }
  d1 <- nmit(cnt[keep, ], sub, method = mth)
  d2 <- nmit(cnt[keep, ], sub, method = mth, matrix_norm = 'nuclear')
  expect_equal(as.matrix(d1), fro, tolerance = 1e-12, info = mth, check.attributes = FALSE)
  expect_equal(as.matrix(d2), nuc, tolerance = 1e-12, info = mth, check.attributes = FALSE)
  expect_identical(rownames(as.matrix(d1)), S)
}

# The samples' order within subjects (and overall) doesn't matter.
o  <- sample.int(length(keep))
for (mth in c('kendall', 'pearson', 'spearman')) {
  d1 <- nmit(cnt[keep, ], sub, method = mth)
  d2 <- nmit(cnt[keep[o], ], sub[o], method = mth)
  expect_true(max(abs(d1 - d2)) < 1e-13, info = mth)
  if (mth == 'kendall') expect_identical(unclass(d1), unclass(d2), info = mth)
}

# Threads, input types and margins: identical.
base <- nmit(cnt[keep, ], sub, cpus = 1L)
for (cp in c(2L, 4L)) {
  expect_identical(nmit(cnt[keep, ], sub, cpus = cp), base, info = cp)
  expect_identical(nmit(cnt[keep, ], sub, matrix_norm = 'nuclear', method = 'pearson', cpus = cp),
                   nmit(cnt[keep, ], sub, matrix_norm = 'nuclear', method = 'pearson', cpus = 1L), info = cp)
}
expect_identical(nmit(t(cnt[keep, ]), sub, margin = 2L), base)
if (requireNamespace('Matrix', quietly = TRUE))
  expect_identical(nmit(as(cnt[keep, ] * 1, 'CsparseMatrix'), sub, cpus = 2L), base)
if (requireNamespace('slam', quietly = TRUE))
  expect_identical(nmit(slam::as.simple_triplet_matrix(cnt[keep, ]), sub, cpus = 2L), base)
# Relative abundances given as such, or norm = 'none' on them.
expect_equal(nmit(rel, sub, norm = 'none'), base, tolerance = 1e-14)

# A factor's levels order the subjects.
f <- factor(sub, levels = rev(sort(unique(sub))))
expect_identical(attr(nmit(cnt[keep, ], f), 'Labels'), levels(f))

# Larger: many subjects and features, threads.
set.seed(3)
big  <- matrix(rpois(240 * 60, 8), 240, 60)
bsub <- rep(sprintf('m%02d', 1:40), each = 6)
b1 <- nmit(big, bsub, method = 'spearman', cpus = 1L)
for (cp in c(2L, 3L, 4L)) expect_identical(nmit(big, bsub, method = 'spearman', cpus = cp), b1, info = cp)

# clust, then permanova() on the subjects.
dc <- nmit(cnt[keep, ], sub, clust = 'average')
expect_true(!is.null(attr(dc, 'order')))
subj <- unique(meta[keep, c('subject', 'group')])
rownames(subj) <- subj$subject
p <- permanova(base, 'group', subj, permutations = 99)
expect_inherits(p, 'ecodive_permanova')

expect_error(nmit(cnt, 'subject', meta, method = 'x'), "`method` must be one of 'kendall', 'pearson', 'spearman'")
expect_error(nmit(cnt, 'subject', meta, matrix_norm = 'x'), "`matrix_norm` must be one of")
expect_error(nmit(cnt, rep('a', nrow(cnt))), 'at least 2 subjects')
expect_error(nmit(cnt), '`subjects` is required')
bad <- cnt * 1; bad[1, 1] <- NA
expect_error(nmit(bad, 'subject', meta), 'must not contain NA')



#--------------------------------------------------
# volatility()
#--------------------------------------------------
vo <- volatility('shannon', 'subject', 'state', 'group', meta)
expect_inherits(vo, 'ecodive_volatility')
expect_identical(names(vo), c('means', 'global', 'trajectories'))
expect_identical(attr(vo, 'conf_level'), 0.95)
expect_identical(attr(vo, 'n_resamples'), 1000L)
expect_identical(attr(vo, 'display')$title, 'Volatility')

# Global limits: mean +/- k sd over the values.
y <- meta$shannon[!is.na(meta$shannon)]
g <- vo$global
expect_equal(unname(g[c('mean', 'sd')]), c(mean(y), sd(y)), tolerance = 1e-14)
expect_equal(unname(g[c('lower_3sd', 'lower_2sd', 'upper_2sd', 'upper_3sd')]),
             mean(y) + c(-3, -2, 2, 3) * sd(y), tolerance = 1e-14)

# Means and counts per group and state.
m <- vo$means
expect_identical(nrow(m), 10L)
expect_identical(levels(m$group), c('A', 'B'))
ok <- !is.na(meta$shannon)
key <- paste(meta$group, meta$state)[ok]
tm  <- tapply(meta$shannon[ok], key, mean)
expect_equal(m$mean, as.vector(tm[paste(m$group, m$state)]), tolerance = 1e-14)
expect_identical(m$n, as.integer(table(key)[paste(m$group, m$state)]))
expect_true(all(m$lower <= m$mean & m$mean <= m$upper))

# The interval: Vega's bootstrapCI() in R, on the same uniforms.
# Resample j takes values floor(u * n) + 1 of uniforms (j - 1) n + 1 .. j n.
boot_ci <- function (v, u, m, alpha) {
  n  <- length(v)
  mu <- colSums(matrix(v[floor(u * n) + 1], n, m)) / n
  stats::quantile(sort(mu), c(alpha / 2, 1 - alpha / 2), type = 7, names = FALSE)
}
for (k in seq_len(nrow(m))) {
  i <- which(ok & meta$group == m$group[k] & meta$state == m$state[k])
  u <- volatility_details(length(i), 1000L, 0L, k - 1L)
  expect_equal(c(m$lower[k], m$upper[k]), boot_ci(meta$shannon[i], u, 1000L, 0.05), tolerance = 1e-14, info = k)
}
v2 <- volatility('shannon', 'subject', 'state', 'group', meta, conf_level = 0.8, n_resamples = 50L, seed = 7L)
u  <- volatility_details(m$n[[3]], 50L, 7L, 2L)
i  <- which(ok & meta$group == 'A' & meta$state == 7)
expect_equal(c(v2$means$lower[3], v2$means$upper[3]), boot_ci(meta$shannon[i], u, 50L, 0.2), tolerance = 1e-14)
expect_identical(volatility('shannon', 'subject', 'state', 'group', meta, seed = 7L)$means[, 1:3], m[, 1:3])
expect_false(identical(volatility('shannon', 'subject', 'state', 'group', meta, seed = 7L)$means$lower, m$lower))
expect_identical(volatility('shannon', 'subject', 'state', 'group', meta), vo)

# Vega draws with Math.random(): its intervals match ecodive's in
# distribution. Coverage and width of 95% intervals for 400 samples of 8
# values from N(0, 1), by ecodive and by Vega's algorithm on R's generator.
# (One more sample at state 2, since volatility() needs two states.)
set.seed(11)
sim <- data.frame(y = rnorm(400 * 8), g = rep(sprintf('g%03d', 1:400), each = 8), t = 1)
sim <- rbind(sim, data.frame(y = 0, g = 'g001', t = 2))
vs  <- volatility('y', NULL, 't', 'g', sim)$means
vs  <- vs[vs$state == 1, ]
rr  <- t(vapply(split(sim$y[sim$t == 1], sim$g[sim$t == 1]), function (v)
  boot_ci(v, stats::runif(1000 * 8), 1000L, 0.05), numeric(2)))
expect_true(abs(mean(vs$lower < 0 & vs$upper > 0) - mean(rr[, 1] < 0 & rr[, 2] > 0)) < 0.05)
expect_true(abs(mean(vs$upper - vs$lower) / mean(rr[, 2] - rr[, 1]) - 1) < 0.02)

# Trajectories.
tr <- vo$trajectories
expect_identical(names(tr), c('sample', 'subject', 'group', 'state', 'value'))
expect_identical(nrow(tr), nrow(meta))
expect_identical(tr$value, meta[tr$sample, 'shannon'])
expect_false(is.unsorted(tr$subject))

# Without groups or subjects; vectors.
x <- volatility(setNames(meta$shannon, meta$sample), NULL, meta$state)
expect_null(x$trajectories)
expect_identical(levels(x$means$group), 'all')
expect_identical(x$global, vo$global)
x <- volatility('shannon', states = 'state', data = meta)
expect_null(x$trajectories)

# One value in a cell: a zero-width interval; no values: NA.
one <- volatility(c(1, 2, 3, NA), c('a', 'a', 'b', 'b'), c(1, 2, 1, 2))
expect_identical(one$means$n, c(2L, 1L))
expect_identical(c(one$means$lower[2], one$means$upper[2]), c(2, 2))
none <- volatility(c(1, NA), c('a', 'a'), c(1, 2))
expect_identical(none$means$n, c(1L, 0L))
expect_true(is.na(none$means$mean[2]) && is.na(none$means$lower[2]))
expect_identical(unname(none$global[['sd']]), 0)

expect_error(volatility('shannon', 'subject', 'state', data = meta, conf_level = 1), '`conf_level`')
expect_error(volatility('shannon', 'subject', 'state', data = meta, n_resamples = 0), '`n_resamples`')
expect_error(volatility('shannon', 'subject', 'subject', data = meta), 'different columns')
expect_error(volatility('shannon', 'subject', 'group', data = meta), '`states` must be numeric')
expect_error(volatility('shannon', 'subject', 'state', 'nope', meta), "'nope' is not one")



#--------------------------------------------------
# Printing and plots.
#--------------------------------------------------
out <- capture.output(print(vo))
expect_true(any(grepl('^Volatility', out)) && any(grepl('95% bootstrap', out)))
expect_true(any(grepl('First Differences', capture.output(print(fd)))))
expect_true(any(grepl('First Distances of Bray-Curtis', capture.output(print(first_distances(bc, 'subject', 'state', meta, replicates = 'drop'))))))

titles <- function (x, ...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control('enable')
  plot(x, ...)
  txt <- unlist(lapply(grDevices::recordPlot()[[1]], function (op) {
    a <- unlist(op[[2]], recursive = TRUE)
    a[vapply(a, is.character, logical(1))]
  }))
  gsub('\n', ' ', unique(unlist(txt)))
}

had_seed <- exists('.Random.seed', envir = globalenv())
seed0    <- if (had_seed) get('.Random.seed', envir = globalenv())
grDevices::pdf(NULL)
par0 <- graphics::par(no.readonly = TRUE)
drawn <- function (x, ...) {
  v <- withVisible(plot(x, ...))
  identical(v$value, x) && !v$visible
}
expect_true(drawn(vo))
expect_true(drawn(vo, spaghetti = FALSE, limits = FALSE, legend = NULL, col = c('red', 'blue')))
expect_true(drawn(volatility('shannon', NULL, 'state', data = meta)))
expect_true(drawn(fd))
expect_true(drawn(ft, feature = 'F03'))
expect_true(drawn(fdd))
expect_true(drawn(first_differences('shannon', 'subject', 'state', meta, replicates = 'drop', baseline = 0)))
# Only the plot's coordinates (usr, xaxp, yaxp) change, as with any plot.
par1 <- graphics::par(no.readonly = TRUE)
expect_identical(par1[setdiff(names(par1), c('usr', 'xaxp', 'yaxp'))], par0[setdiff(names(par0), c('usr', 'xaxp', 'yaxp'))])
expect_error(plot(ft, feature = 'nope'), '`feature` must be the name')
expect_error(plot(fd[0, ]), 'no rows')
expect_error(plot(vo, legend = 'middle'), '`legend`')
grDevices::dev.off()
expect_identical(exists('.Random.seed', envir = globalenv()), had_seed)
if (had_seed) expect_identical(get('.Random.seed', envir = globalenv()), seed0)

expect_true('Volatility' %in% titles(vo))
expect_true(all(c('state', 'shannon', 'group') %in% titles(vo)))
expect_true('First Differences' %in% titles(fd))
expect_true('Differences from state 0' %in% titles(first_differences('shannon', 'subject', 'state', meta, replicates = 'drop', baseline = 0)))
expect_true('First Distances of Bray-Curtis Dissimilarity' %in% titles(first_distances(bc, 'subject', 'state', meta, replicates = 'drop')))
expect_true('x' %in% titles(vo, main = 'x'))
