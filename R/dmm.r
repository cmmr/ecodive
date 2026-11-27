# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Written from Holmes, Harris & Quince (2012) (clean room). Bioconductor's
# DirichletMultinomial (LGPL) and mothur (GPL) were only run, to compare
# outputs; their source was not read. See src/dmm.c.


#' Dirichlet-multinomial mixtures
#'
#' Clusters samples into community types with the Dirichlet-multinomial
#' mixture (DMM) model of Holmes, Harris & Quince (2012): each type is a
#' Dirichlet distribution of compositions, and each sample's counts are a
#' multinomial draw from a composition drawn from its type. The model is
#' fitted by expectation maximization (EM) for each number of types, and the
#' number of types is chosen by the Laplace approximation to the model
#' evidence.
#'
#' @param counts  A table of raw counts (whole numbers), such as a matrix,
#'        a `dgCMatrix` or a phyloseq object. Samples are rows unless
#'        `margin = 2`.
#'
#' @param n_types  The numbers of types to fit (each from 1 to the number of
#'        samples); a single value fixes it. Default `1:7`, cut to `1:n` for
#'        fewer than 7 samples.
#'
#' @param criterion  Which criterion chooses the number of types, the
#'        smallest value winning: `"laplace"` (the default), the negative
#'        log of the Laplace approximation to the model evidence, as in
#'        Holmes et al. (2012); `"aic"`; or `"bic"`. All three are reported
#'        in the `selection` attribute. See Details.
#'
#' @param n_starts  The number of random starts for each number of types;
#'        the fit with the largest posterior is kept. Default `1`.
#'
#' @param max_iters  The most EM iterations for each start. Default `1000`.
#'
#' @param tol  EM stops when the negative log posterior changes by less than
#'        `tol` times its value. Default `1e-10`.
#'
#' @param eta,nu  The shape and rate of the Gamma prior on every Dirichlet
#'        parameter. Defaults `0.1` and `0.1`, as in DirichletMultinomial.
#'
#' @param seed  An integer seed for the random starts. The same `seed`
#'        gives the same result. R's random number generator is not used.
#'        Default `0`.
#'
#' @param margin  If your samples are in the table's rows, set to `1L`
#'        (the default). If they are in columns, set to `2L`.
#'
#' @param cpus  How many parallel processing threads to use. Results are
#'        the same for any number. Default: all available cores.
#'
#' @details
#' **The model** (Holmes et al. 2012). Type \eqn{k} has weight \eqn{\pi_k}
#' and Dirichlet parameters \eqn{\alpha_{k1}, \dots, \alpha_{kS}}, with sum
#' \eqn{A_k}, its precision: large values mean the type's samples have
#' similar compositions, small values that they vary widely. The
#' probability of sample \eqn{i}'s counts \eqn{x_{i1}, \dots, x_{iS}}, with
#' total \eqn{J_i}, under type \eqn{k} is
#' \deqn{p(x_i | \alpha_k) = \frac{J_i!}{\prod_j x_{ij}!}
#'   \frac{\Gamma(A_k)}{\Gamma(A_k + J_i)}
#'   \prod_j \frac{\Gamma(\alpha_{kj} + x_{ij})}{\Gamma(\alpha_{kj})}.}
#' Every \eqn{\alpha_{kj}} has a Gamma(`eta`, `nu`) prior (shape and rate),
#' and is fitted on the log scale, \eqn{\lambda_{kj} = \log \alpha_{kj}}.
#'
#' **Fitting.** Each number of types starts from k-means clusters of the
#' samples' proportions (the best of 10 runs of k-means++ seeding and Lloyd's
#' iterations, drawn from the `seed`; one cluster for one type), with each
#' type's \eqn{\alpha} from its cluster's moments. EM then alternates: each
#' sample's posterior probability of each type, \eqn{z_{ik} \propto \pi_k
#' p(x_i | \alpha_k)}; then for each type the \eqn{\lambda_k} that maximise
#' \eqn{\sum_i z_{ik} \log p(x_i | \alpha_k)} plus the log prior, by
#' Newton's method, and \eqn{\pi_k = \frac1n \sum_i z_{ik}}. It stops when
#' the negative log posterior changes by less than `tol` times its value.
#' With `n_starts` above 1, k-means starts again from the next seeds, and
#' the fit with the smallest negative log posterior is kept. Each sample is
#' assigned to its most probable type. Types are numbered by decreasing
#' weight.
#'
#' **Choosing the number of types.** The `selection` attribute reports, for
#' each number of types \eqn{K}, with \eqn{M = KS + K - 1} parameters:
#' * `nle`: the negative log posterior at the fit,
#'   \eqn{-\log p(X | \hat\theta) - \log p(\hat\theta)}, including the
#'   multinomial coefficients and the prior's constants.
#' * `laplace`: the negative log of the Laplace approximation to the
#'   evidence, \eqn{\mathrm{nle} + \frac12 \log |H| - \frac M2 \log 2\pi},
#'   where \eqn{H} is the Hessian of the negative log posterior: one block
#'   per type, in \eqn{\lambda_k}, and the block of the multinomial in the
#'   weights.
#' * `aic`: \eqn{\mathrm{nle} + M}, and `bic`: \eqn{\mathrm{nle} + \frac M2
#'   \log n}. Both are on the scale of `nle`, half the usual AIC and BIC, as
#'   DirichletMultinomial reports them.
#'
#' The \eqn{K} with the smallest value of `criterion` is chosen. Holmes et al.
#' choose by the Laplace approximation, and argue that AIC and BIC can mislead
#' for these models.
#'
#' **Cost.** Each EM iteration takes time proportional to the number of
#' non-zero counts times the number of types, so sparse tables stay cheap.
#' The Hessian of each type is a diagonal plus a rank-one matrix, so Newton
#' steps and \eqn{\log|H|} need no \eqn{S \times S} matrix.
#'
#' **Counts.** The model is for raw counts: the multinomial step accounts for
#' each sample's sequencing depth. Proportions, rarefied averages or other
#' scaled values give a meaningless fit, and a warning. Aggregate to genus
#' (or another rank) first to keep the number of features moderate, as
#' Holmes et al. did.
#'
#' @section Differences from DirichletMultinomial:
#' * **The optimizer.** DirichletMultinomial (from Holmes et al.'s code)
#'   uses GSL's BFGS for each type's \eqn{\lambda}; `dmm()` uses Newton's
#'   method, and stops EM by `tol`. From the same starting weights both
#'   reach the same fit: the same posteriors, weights and negative log
#'   posterior (`dmm()`'s slightly lower). Small Dirichlet parameters (well
#'   under 0.01), where the posterior is nearly flat, can differ by up to
#'   about 1%, as DirichletMultinomial stops further from the optimum.
#' * **The start.** DirichletMultinomial starts from soft k-means with its
#'   own random numbers, so with more than one type the two may reach
#'   different local optima. Try several `n_starts`.
#' * **The Laplace approximation** matches DirichletMultinomial's for one
#'   type. For more types, its reported Hessian is not positive definite,
#'   and its intervals for the Dirichlet parameters are `NaN`; `dmm()` uses
#'   the Hessian of the negative log posterior described above, so
#'   `laplace` differs for more than one type, and may choose another
#'   number of types. `nle`, `aic` and `bic` match DirichletMultinomial's
#'   `goodnessOfFit()` (`NLE`, `AIC`, `BIC`).
#' * `alpha` and `profiles` are types \eqn{\times} features, where
#'   DirichletMultinomial's `fitted()` is features \eqn{\times} types.
#'   Types are ordered by decreasing weight.
#'
#' @return A data frame of class
#'   `c("ecodive_dmm", "ecodive_types", "data.frame")`, one row per sample
#'   (see [ecodive_types]), with columns `sample`, `type` (a factor, levels
#'   `"1"`, `"2"`, ... in decreasing weight) and `probability` (the posterior
#'   probability of the assigned type). Its attributes:
#'   * `types`: a data frame, one row per type, with `type`, `n_samples`,
#'     `top_feature` (the feature of largest Dirichlet mean), `weight`
#'     (\eqn{\pi_k}) and `precision` (\eqn{A_k}).
#'   * `profiles`: the types \eqn{\times} features matrix of Dirichlet means,
#'     \eqn{\alpha_{kj} / A_k}.
#'   * `selection`: a data frame, one row per number of types, with
#'     `n_types`, `laplace`, `aic`, `bic`, `nle`, `n_iters` (EM iterations),
#'     `converged` and `selected`.
#'   * `posterior`: the samples \eqn{\times} types matrix of posterior
#'     probabilities \eqn{z_{ik}}.
#'   * `alpha`: the types \eqn{\times} features matrix of Dirichlet
#'     parameters.
#'   * `criterion`: the criterion that chose the number of types.
#'   * `method`: `"dmm"`.
#'
#'   `plot()` draws the types' profiles or the selection (see
#'   [plot.ecodive_types()]); its title is "Dirichlet-Multinomial Mixtures".
#'
#' @section Same names in other packages:
#' `dmm::dmm()` (CRAN) fits dyadic mixed models, which are unrelated.
#' DirichletMultinomial calls this model `dmn()`.
#'
#' @references
#' Holmes, I., Harris, K., & Quince, C. (2012). Dirichlet multinomial
#' mixtures: generative models for microbial metagenomics. *PLoS ONE*, 7,
#' e30126. \doi{10.1371/journal.pone.0030126}
#'
#' Ding, T., & Schloss, P. D. (2014). Dynamics and associations of microbial
#' community types across the human body. *Nature*, 509, 357-360.
#' \doi{10.1038/nature13178}
#'
#' Arthur, D., & Vassilvitskii, S. (2007). k-means++: the advantages of
#' careful seeding. *Proceedings of the 18th ACM-SIAM Symposium on Discrete
#' Algorithms*, 1027-1035.
#'
#' @family community types
#' @seealso [ecodive_types], [plot.ecodive_types()], [enterotype()]
#' @export
#' @examples
#'     # Counts of 60 samples from three community types.
#'     set.seed(1)
#'     probs <- list(
#'       c(40, 5, 5, 5, 5, 1, 1), c(5, 40, 5, 1, 1, 5, 5), c(5, 5, 5, 20, 20, 1, 1) )
#'     counts <- t(sapply(rep(1:3, each = 20), function (k)
#'       rmultinom(1, 500, rgamma(7, probs[[k]])) ))
#'     colnames(counts) <- paste0('Taxon', 1:7)
#'
#'     dm <- dmm(counts, n_types = 1:5)
#'     dm
#'     table(dm$type, rep(1:3, each = 20))
#'
#'     plot(dm)
#'     plot(dm, 'selection')
#'
#'     # The types as groups in a test.
#'     permanova(bray(counts), dm$type)
dmm <- function (
    counts, n_types = 1:7, criterion = c('laplace', 'aic', 'bic'),
    n_starts = 1L, max_iters = 1000L, tol = 1e-10, eta = 0.1, nu = 0.1,
    seed = 0L, margin = 1L, cpus = n_cpus() ) {

  dmm_core(
    counts, n_types, criterion, n_starts, max_iters, tol, eta, nu, seed,
    margin, cpus, default_k = missing(n_types) )
}



# dmm(), with an internal hook: `init`, an n x K matrix of starting weights
# (one number of types), replacing the k-means starts (for comparisons with
# other software from shared starting values).
dmm_core <- function (
    counts, n_types = 1:7, criterion = c('laplace', 'aic', 'bic'),
    n_starts = 1L, max_iters = 1000L, tol = 1e-10, eta = 0.1, nu = 0.1,
    seed = 0L, margin = 1L, cpus = n_cpus(), default_k = FALSE, init = NULL ) {

  criterion <- match.arg(criterion)
  validate_counts()
  validate_margin()
  validate_cpus()

  if (is.logical(counts) && is.matrix(counts))
    storage.mode(counts) <- 'double'

  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  if (any(vals %% 1 != 0))
    warning(
      '`counts` has values that are not whole numbers. dmm() models raw ',
      'counts: proportions or scaled values give a meaningless fit.', call. = FALSE )
  remove('vals')

  smp     <- count_samples(counts, margin)
  n       <- smp$n
  samples <- if (is.null(smp$labels)) as.character(seq_len(n)) else smp$labels
  features <- count_samples(counts, 3L - margin)$labels
  S       <- count_samples(counts, 3L - margin)$n
  if (is.null(features)) features <- as.character(seq_len(S))

  if (n < 2)
    stop('dmm() needs at least 2 samples.', call. = FALSE)
  if (S < 2)
    stop('dmm() needs at least 2 features.', call. = FALSE)

  tot <- sample_totals(counts, margin)
  if (any(tot <= 0))
    stop(
      'Every sample needs a positive total; these have none: ',
      name_list(samples[tot <= 0]), '.', call. = FALSE )

  n_types   <- dmm_n_types(n_types, n, default_k)
  n_starts  <- dmm_count(n_starts,  'n_starts')
  max_iters <- dmm_count(max_iters, 'max_iters')
  validate_seed()
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol < 0)
    stop('`tol` must be a single non-negative number.', call. = FALSE)
  for (nm in c('eta', 'nu')) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v <= 0)
      stop('`', nm, '` must be a single positive number.', call. = FALSE)
  }

  if (!is.null(init)) {
    stopifnot(length(n_types) == 1, is.matrix(init), nrow(init) == n, ncol(init) == n_types)
    storage.mode(init) <- 'double'
  }


  #--------------------------------------------------
  # EM for each number of types (src/dmm.c).
  #--------------------------------------------------
  fits <- .Call(
    C_dmm, counts, margin, n_types, n_starts, max_iters, as.double(tol),
    seed, as.double(eta), as.double(nu), init, cpus )

  bad <- vapply(fits, function (f) f[[1]] == 0L, TRUE)
  if (any(bad))
    stop(
      'Could not start ', paste(n_types[bad], collapse = ', '),
      ' types: there are fewer distinct samples. Try smaller `n_types`.', call. = FALSE )

  M    <- n_types * S + n_types - 1
  nle  <- vapply(fits, `[[`, 0, 2)
  ldet <- vapply(fits, `[[`, 0, 3)

  selection <- data.frame(
    n_types   = n_types,
    laplace   = nle + ldet / 2 - M / 2 * log(2 * pi),
    aic       = nle + M,
    bic       = nle + M / 2 * log(n),
    nle       = nle,
    n_iters   = vapply(fits, `[[`, 0L, 4),
    converged = vapply(fits, `[[`, NA, 5) )

  crit <- selection[[criterion]]
  best <- which.min(replace(crit, is.nan(crit), Inf))
  if (length(best) == 0) best <- 1L # nocov
  selection$selected <- seq_along(n_types) == best

  if (!all(selection$converged))
    warning(
      'EM did not converge within ', max_iters, ' iterations for ',
      paste(n_types[!selection$converged], collapse = ', '), ' types.', call. = FALSE )
  if (anyNA(selection$laplace))
    warning(
      'The Laplace approximation is undefined (the Hessian is not positive ',
      'definite) for ', paste(n_types[is.na(selection$laplace)], collapse = ', '),
      ' types.', call. = FALSE )


  #--------------------------------------------------
  # The result.
  #--------------------------------------------------
  fit   <- fits[[best]]
  K     <- n_types[[best]]
  labs  <- as.character(seq_len(K))
  alpha <- fit[[7]]
  post  <- fit[[9]]
  dimnames(alpha) <- list(labs, features)
  dimnames(post)  <- list(samples, labs)

  A        <- rowSums(alpha)
  profiles <- alpha / A
  assign   <- max.col(post, ties.method = 'first')
  type     <- factor(assign, levels = seq_len(K), labels = labs)

  types <- data.frame(
    stringsAsFactors = FALSE,
    type        = factor(labs, levels = labs),
    n_samples   = tabulate(assign, K),
    top_feature = types_top_feature(profiles, K),
    weight      = fit[[8]],
    precision   = unname(A) )

  df <- data.frame(
    stringsAsFactors = FALSE,
    sample      = samples,
    type        = type,
    probability = post[cbind(seq_len(n), assign)] )

  new_ecodive_types(
    df, 'dmm', types, profiles, selection,
    posterior = post,
    alpha     = alpha,
    criterion = criterion )
}


# Checks `n_types` against `n` samples. The default 1:7 is cut to 1:n.
dmm_n_types <- function (n_types, n, default = FALSE) {
  if (default) return (seq_len(min(7L, n)))
  if (!is.numeric(n_types) || length(n_types) == 0 || anyNA(n_types) ||
      any(n_types %% 1 != 0))
    stop('`n_types` must be whole numbers.', call. = FALSE)
  n_types <- sort(unique(as.integer(n_types)))
  if (any(n_types < 1) || any(n_types > n))
    stop('`n_types` must be between 1 and ', n, ' (the number of samples).', call. = FALSE)
  n_types
}


dmm_count <- function (x, what) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 1 || x %% 1 != 0 ||
      x > .Machine$integer.max)
    stop('`', what, '` must be a single positive whole number.', call. = FALSE)
  as.integer(x)
}

