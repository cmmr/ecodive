# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# References as given by
# https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282


#' Alpha Diversity Wrapper Function
#' 
#' @inherit documentation
#' 
#' @param method   The name of an alpha diversity metric. One of `c('ace',
#'   'berger', 'brillouin', 'chao1', 'faith', 'fisher', 'inv_simpson',
#'   'margalef', 'mcintosh', 'menhinick', 'observed', 'shannon', 'simpson',
#'   'squares')`. Case-insensitive and partial name matching is supported.
#'   Programmatic access via `list_methods('alpha')`.
#'   
#' @param ...  Additional options to pass through to the called function. I.e.
#'   `cpus` or `tree`.
#' 
#' @return A numeric vector.
#' 
#' @details
#' 
#' ## Integer Count Requirements
#' 
#' A frequent and critical error in alpha diversity analysis is providing the
#' wrong type of data to a metric's formula. Some indices are mathematically
#' defined based on counts of individuals and require raw, integer abundance
#' data. Others are based on proportional abundances and can accept either
#' integer counts (which are then converted to proportions) or pre-normalized
#' proportional data. Using proportional data with a metric that requires
#' integer counts will return an error message.
#' 
#' | Requires Integer Counts Only | Can Use Proportional Data     |
#' | :--------------------------  | :---------------------------- |
#' | Chao1                        | Observed Features             |
#' | ACE                          | Shannon Index                 |
#' | Squares Richness Estimator   | Gini-Simpson Index            |
#' | Margalef's Index             | Inverse Simpson Index         |
#' | Menhinick's Index            | Berger-Parker Index           |
#' | Fisher's Alpha               | McIntosh Index                |
#' | Brillouin Index              | Faith's PD (presence/absence) |
#' 
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Shannon diversity values
#'     alpha_div(ex_counts, 'Shannon')
#'     
#'     # Chao1 diversity values
#'     alpha_div(ex_counts, 'c')
#'     
#'     # Faith PD values
#'     alpha_div(ex_counts, 'faith', tree = ex_tree)
#'     
#'     
alpha_div <- function (x, method, ...) {
  match_method(method, div = 'alpha')$func(x = x, ...)
}




