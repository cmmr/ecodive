# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Pairwise tests (`pairwise = TRUE`) for permanova(), anosim() and mrpp():
# each pair of groups is tested by running the same test on just those two
# groups' samples, as q2-diversity's beta_group_significance(pairwise=True)
# does (BSD-3, (c) 2016-2026 QIIME 2 development team), and the pairs'
# permutation p-values are adjusted for multiple testing with
# stats::p.adjust(). permdisp() has t tests of its own, as vegan's
# permutest(pairwise = TRUE), and only shares pairwise_method() and the
# printout.



#------------------------------------------------------------------------------
# Checks `pairwise` and `p_adjust`. Returns the p.adjust() method (matched as
# p.adjust() matches it, so 'bonf' is 'bonferroni'), or NULL when `pairwise`
# is FALSE. `p_given` is !missing(p_adjust) in the caller: p_adjust means
# nothing without pairwise = TRUE, so giving it then is an error.
#------------------------------------------------------------------------------
pairwise_method <- function (pairwise, p_adjust, p_given) {

  if (!isTRUE(pairwise) && !isFALSE(pairwise))
    stop('`pairwise` must be TRUE or FALSE.', call. = FALSE)

  if (!pairwise) {
    if (p_given)
      stop('`p_adjust` is only used with `pairwise = TRUE`.', call. = FALSE)
    return (NULL)
  }

  opts <- stats::p.adjust.methods
  m    <- if (is.character(p_adjust) && length(p_adjust) == 1 && !is.na(p_adjust))
            pmatch(p_adjust, opts)
          else NA_integer_

  if (is.na(m))
    stop(
      '`p_adjust` must be one of ', paste0("'", opts, "'", collapse = ', '),
      ' (see stats::p.adjust()).', call. = FALSE )

  opts[[m]]
}



#------------------------------------------------------------------------------
# The pairs' table.
#
#   dv        dist_vector() of the whole distance matrix.
#   grp       Each sample's group: a factor, aligned with dv's samples. Its
#             levels, in order, give the pairs, as utils::combn().
#   strata    Each sample's stratum (any values), or NULL.
#   perms     The `permutations` argument: a count, or a matrix with one
#             permutation of 1:n per row.
#   seed      The test's seed. Pair k's permutations come from
#             pair_seed(seed, k).
#   method      The p.adjust() method, from pairwise_method().
#   stat_names  The names of the statistics: the table's first columns.
#   test        function (dm, groups, strata, permutations, seed) running
#               the test on one pair's samples, and returning a list of
#               `stats` (the named statistics), `z`, `n_perm`, and
#               `p_value`, as pair_result().
#
# Returns a data frame with one row per pair, named "a-b": the statistics,
# z, n (the pair's samples), n_perm, p_value, and p_adjusted. A pair the test
# fails on gets NAs, with a warning; the test's own warnings are passed on,
# naming the pair.
#------------------------------------------------------------------------------
pairwise_table <- function (dv, grp, strata, perms, seed, method, stat_names, test) {

  lev <- levels(grp)
  if (length(lev) < 2) return (NULL)

  combos <- utils::combn(length(lev), 2)
  labs   <- paste(lev[combos[1, ]], lev[combos[2, ]], sep = '-')
  codes  <- as.integer(grp)

  rows <- lapply(seq_len(ncol(combos)), function (k) {

    idx <- which(codes %in% combos[, k])
    sub <- dist_subset(dv, idx)
    g   <- droplevels(grp[idx])
    s   <- if (is.null(strata)) NULL else strata[idx]
    p   <- if (is.matrix(perms)) induced_perms(perms, idx) else perms

    failed <- NULL
    res <- withCallingHandlers(
      tryCatch(
        test(sub, g, s, p, pair_seed(seed, k)),
        error = function (e) { failed <<- conditionMessage(e); NULL } ),
      warning = function (w) {
        warning('Pair ', labs[[k]], ': ', conditionMessage(w), call. = FALSE)
        invokeRestart('muffleWarning')
      })

    if (is.null(res)) {
      warning('Pair ', labs[[k]], ' could not be tested: ', failed, call. = FALSE)
      res <- list(
        stats   = stats::setNames(rep(NA_real_, length(stat_names)), stat_names),
        z       = NA_real_,
        n_perm  = if (is.matrix(perms)) nrow(perms) else as.integer(perms),
        p_value = NA_real_ )
    }

    c(as.list(res$stats),
      list(z = res$z, n = length(idx), n_perm = res$n_perm, p_value = res$p_value))
  })

  out <- do.call(rbind.data.frame, rows)
  out$n      <- as.integer(out$n)
  out$n_perm <- as.integer(out$n_perm)
  out$p_adjusted <- stats::p.adjust(out$p_value, method)
  rownames(out)  <- labs
  out
}



#------------------------------------------------------------------------------
# A dist object of the samples `idx` (increasing) of dist_vector() `dv`,
# taken straight from the packed lower triangle.
#------------------------------------------------------------------------------
dist_subset <- function (dv, idx) {

  m  <- length(idx)
  n  <- dv$n
  ij <- utils::combn(idx, 2)  # columns (i, j), i < j, in dist order
  i  <- ij[1, ]
  j  <- ij[2, ]
  # Element (i, j) of a dist of size n, i < j: column i starts after the
  # (i - 1) columns before it.
  at <- (i - 1) * (n - i / 2) + (j - i)

  labels <- if (is.null(dv$labels)) NULL else dv$labels[idx]
  structure(
    dv$dm[at],
    Size   = m,
    Labels = labels,
    Diag   = FALSE,
    Upper  = FALSE,
    class  = 'dist' )
}



#------------------------------------------------------------------------------
# The permutations of the samples `idx` (increasing) induced by a matrix of
# permutations of 1:n (one per row): each sample of `idx` is sent where the
# permutation's cycle through it first meets another sample of `idx`. This is
# a permutation of the pair's samples, uniformly distributed when the full
# one is, and it stays within strata when the full one does. Returns a
# matrix with one permutation of 1:length(idx) per row.
#------------------------------------------------------------------------------
induced_perms <- function (perms, idx) {

  n   <- ncol(perms)
  hit <- logical(n)
  hit[idx] <- TRUE

  if (nrow(perms) == 0) return (matrix(integer(0), 0, length(idx)))

  # All of a pair's samples, in every permutation, step together.
  rows <- rep(seq_len(nrow(perms)), each = length(idx))
  to   <- perms[cbind(rows, rep(idx, nrow(perms)))]
  out  <- !hit[to]
  while (any(out)) {
    to[out]  <- perms[cbind(rows[out], to[out])]
    out[out] <- !hit[to[out]]
  }

  pos <- integer(n)
  pos[idx] <- seq_along(idx)
  matrix(pos[to], nrow = nrow(perms), byrow = TRUE)
}



# The seed for pair k: see C_pair_seed() in src/permute.c.
pair_seed <- function (seed, k)
  .Call(C_pair_seed, as.integer(seed), as.integer(k))



#------------------------------------------------------------------------------
# Lines for test_details(): the pairs' table, naming the p.adjust() method.
#------------------------------------------------------------------------------
pairwise_details <- function (x, digits, heading = 'Pairwise comparisons') {

  if (is.null(x$pairwise)) return (character(0))

  c('',
    paste0(heading, if (x$p_adjust == 'none') ' (p-values not adjusted):'
                    else paste0(' (p-values adjusted by ', p_adjust_name(x$p_adjust), '):') ),
    utils::capture.output(print(x$pairwise, digits = digits)) )
}


p_adjust_name <- function (method) {
  switch(method,
    holm       = "Holm's method",
    hochberg   = "Hochberg's method",
    hommel     = "Hommel's method",
    bonferroni = 'Bonferroni',
    BH         = 'Benjamini-Hochberg',
    fdr        = 'Benjamini-Hochberg',
    BY         = 'Benjamini-Yekutieli',
    method )
}


# A pair's test result as pairwise_table() wants it: the statistic, then the
# fields named in `extra` (such as permanova()'s r2, as "R2").
pair_result <- function (res, extra = NULL, extra_names = extra) {
  stats <- c(res$statistic, stats::setNames(unlist(res[extra]), extra_names))
  list(stats = stats, z = res$z, n_perm = res$n_perm, p_value = res$p_value)
}
