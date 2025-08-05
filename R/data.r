# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Example counts matrix
#' 
#' Genera found on four human body sites.
#' 
#' @format A matrix of 4 samples (columns) x 6 genera (rows).
#' @source Derived from The Human Microbiome Project dataset.
#' <https://commonfund.nih.gov/hmp>
"ex_counts"



#' Example phylogenetic tree
#' 
#' Companion tree for `ex_counts`.
#' 
#' @details
#' `ex_tree` encodes this tree structure:
#' ```
#'       +----------44---------- Haemophilus
#'   +-2-|
#'   |   +----------------68---------------- Bacteroides  
#'   |                      
#'   |             +---18---- Streptococcus
#'   |      +--12--|       
#'   |      |      +--11-- Staphylococcus
#'   +--11--|              
#'          |      +-----24----- Corynebacterium
#'          +--12--|
#'                 +--13-- Propionibacterium
#' ```
#' 
#' @format A `phylo` object.
"ex_tree"
