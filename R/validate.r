# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


NORM_PERCENT <- 1L
NORM_CLR     <- 2L
NORM_CHORD   <- 3L
NORM_BINARY  <- 4L


validate_args <- function () {
  env  <- parent.frame()
  args <- ls(env)
  
  # move counts, norm, pseudocount, and margin to head of the line
  args <- unique(c(intersect(c('counts', 'norm', 'pseudocount', 'margin'), args), sort(args)))
  
  for (arg in args)
    do.call(paste0('validate_', arg), list(env))
}


validate_alpha <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      if (!inherits(alpha, 'numeric'))
        alpha <- as.numeric(alpha)
      
      stopifnot(length(alpha) == 1)
      stopifnot(!is.na(alpha))
      stopifnot(alpha >= 0 && alpha <= 1)
    }),
    
    error = function (e) 
      stop(e$message, '\n`alpha` must be a single number between 0 and 1.')
  )
}


validate_counts <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      # Pull a tree from complex counts object.
      if (exists('tree', inherits = FALSE) && is.null(tree)) {
        
        if (inherits(counts, 'phyloseq')) {
          tree <- counts@phy_tree
        }
        else if (inherits(counts, 'rbiom')) {
          tree <- counts$tree
        }
        else if (inherits(counts, 'TreeSummarizedExperiment')) {
          tree <- counts@rowTree[[1]]
        }
        else {
          tree <- attr(counts, 'tree', exact = TRUE)
        }
      }
      
      # Derive matrix from simple vector or complex object.
      if (!inherits(counts, c('matrix', 'dgCMatrix', 'dgTMatrix', 'dgeMatrix', 'simple_triplet_matrix'))) {
        
        if (inherits(counts, 'rbiom')) {
          counts <- counts$counts # dgCMatrix
          margin <- 2L
        }
        
        else if (inherits(counts, 'phyloseq')) {
          margin <- ifelse(counts@otu_table@taxa_are_rows, 2L, 1L)
          counts <- counts@otu_table@.Data
        }
        
        else if (inherits(counts, 'TreeSummarizedExperiment')) {
          counts <- counts@assays@data[[1]]
          margin <- 2L
        }
        
        else if (inherits(counts, 'SummarizedExperiment')) {
          counts <- counts@assays@data[[1]]
          margin <- 2L
        }
        
        else if (is.vector(counts)) {
          counts <- matrix(data = counts, ncol = 1)
          margin <- 2L
        }
        
        else {
          counts <- as.matrix(counts)
        }
        
      }
      
      stopifnot(nrow(counts) > 0)
      stopifnot(ncol(counts) > 0)
    }),
    
    error = function (e) 
      stop(e$message, '\n`counts` must be a valid numeric matrix.')
  )
}


validate_cpus <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      stopifnot(is.numeric(cpus))
      stopifnot(length(cpus) == 1)
      stopifnot(!is.na(cpus))
      stopifnot(cpus > 0)
      stopifnot(cpus %% 1 == 0)
      
      if (!is.integer(cpus))
        cpus <- as.integer(cpus)
    }),
    
    error = function (e) 
      stop(e$message, '\n`cpus` must be an integer greater than 0.')
  )
}


validate_cutoff <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      stopifnot(is.numeric(cutoff))
      stopifnot(length(cutoff) == 1)
      stopifnot(!is.na(cutoff))
      stopifnot(cutoff > 0)
      stopifnot(cutoff %% 1 == 0)
      
      cutoff <- as.integer(cutoff)
      
    }),
    
    error = function (e) 
      stop(e$message, '\n`cutoff` must be a positive integer.')
  )
}


validate_depth <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      if (!is.null(depth)) {
        
        stopifnot(is.numeric(depth))
        stopifnot(length(depth) == 1)
        stopifnot(!is.na(depth))
        stopifnot(depth > 0)
        stopifnot(depth %% 1 == 0)
        
        depth <- as.integer(depth)
      }
      
    }),
    
    error = function (e) 
      stop(e$message, '\n`depth` must be a positive integer or NULL.')
  )
}


validate_digits <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      stopifnot(is.numeric(digits))
      stopifnot(length(digits) == 1)
      stopifnot(!is.na(digits))
      stopifnot(digits >= 0)
      stopifnot(digits <= 10)
      stopifnot(digits %% 1 == 0)
      
      digits <- as.integer(digits)
      
    }),
    
    error = function (e) 
      stop(e$message, '\n`digits` must be an integer between 0 and 10.')
  )
}


validate_drop <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      stopifnot(identical(drop, TRUE) || identical(drop, FALSE))
    }),
    error = function (e) 
      stop(e$message, '\n`drop` must be either TRUE or FALSE.')
  )
}


validate_newick <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      stopifnot(is.character(newick))
      stopifnot(length(newick) == 1)
      stopifnot(!is.na(newick))
      newick <- trimws(newick)
      stopifnot(nchar(newick) > 0)
    }),
    
    error = function (e) 
      stop(e$message, '\n`newick` must be a character string.')
  )
}


validate_pairs <- function (env = parent.frame()) {
  
  with(env, {
    if (ncol(counts) < 2)
      stop('`counts` must have at least two samples.')
  })
  
  tryCatch(
    with(env, {
      
      if (!is.null(pairs)) {
        
        n_samples   <- ncol(counts)
        n_distances <- n_samples * (n_samples - 1) / 2
        
        if (is.function(pairs))
          pairs <- local({
            m <- utils::combn(n_samples, 2)
            mapply(pairs, m[1,], m[2,])
          })
        
        if (is.logical(pairs)) {
          if (length(pairs) != n_distances)
            stop('logical vector `pairs` must have length ', n_distances)
          pairs <- which(pairs)
        }
        else if (is.numeric(pairs)) {
          if (is.double(pairs)) {
            if (any(pairs %% 1 > 0)) stop('non-integer values')
            pairs <- as.integer(pairs)
          }
          if (!all(pairs >= 1 & pairs <= n_distances))
            stop('expected `pairs` values between 1 and ', n_distances)
        }
        else {
          stop('cannot be ', typeof(pairs))
        }
        
        remove('n_samples', 'n_distances')
      }
    }),
    
    error = function (e) 
      stop(e$message, '\n`pairs` must be a numeric or logical vector.')
  )
}


validate_power <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      if (!inherits(power, 'numeric'))
        power <- as.numeric(power)
      
      stopifnot(length(power) == 1)
      stopifnot(!is.na(power))
    }),
    
    error = function (e) 
      stop(e$message, '\n`power` must be a single number.')
  )
}


validate_pseudocount <- function (env = parent.frame()) {
  with(env, {
    
    if (!identical(norm, NORM_CLR)) {
      pseudocount <- NULL
    }
    
    else {
      
      mtx_pkg <- get_matrix_package(counts)
      
      val_range <- switch(
        mtx_pkg,
        'base'   = range(counts),
        'slam'   = range(counts$v),
        'Matrix' = range(counts@x) )
      
      if (length(val_range) == 2 && !all(is.finite(val_range)))
        stop('`counts` contains non-finite values; cannot perform CLR normalization.')
      
      min_val <- val_range[1]
      max_val <- val_range[2]
      
      if (min_val < 0)
        stop('`counts` contains negative values; cannot perform CLR normalization.')
      
      if (max_val == 0)
        stop('`counts` contains only zeros; cannot perform CLR normalization.')
      
      if (min_val == max_val)
        stop('`counts` contains a single unique value; cannot perform CLR normalization.')
      
      
      if (!is.null(pseudocount)) {
        
        tryCatch(
          expr = {
            stopifnot(is.numeric(pseudocount))
            pseudocount <- as.double(pseudocount)
            stopifnot(length(pseudocount) == 1)
            stopifnot(!is.na(pseudocount))
            stopifnot(pseudocount > 0)
          },
          error = function (e) {
            stop('`pseudocount` must be a single positive number.')
        })
        
      }
      else {
        
        has_zeros <- min_val == 0
        
        # Check for zeros in sparse matrices
        if (!has_zeros && mtx_pkg != 'base')
          has_zeros <- switch(
            mtx_pkg,
            'slam'   = length(counts$v) < counts$nrow * counts$ncol,
            'Matrix' = length(counts@x) < prod(dim(counts)) )
        
        
        if (!has_zeros) {
          
          pseudocount <- NULL
          
        } else {
          
          # Calculate a "smart" default (half the smallest non-zero value)
          # This is generally safer than '1' for proportional data.
          pseudocount <- switch(
            mtx_pkg,
            'base'   = min(counts[counts > 0]),
            'slam'   = min(counts$v[counts$v > 0]),
            'Matrix' = min(counts@x[counts@x > 0]) )
          
          pseudocount <- pseudocount / 2
          
          warning(
            call. = FALSE, 
            paste0(
              "Zeros detected in data. Using pseudocount = ", format(pseudocount, digits = 3), "\n",
              "Calculated using formula: min(counts[counts > 0]) / 2.\n",
              "To suppress this warning, provide an explicit `pseudocount` argument." ))
        }
        
        remove('has_zeros')
      }
      
      remove('mtx_pkg', 'val_range', 'min_val')
    }
  })
}


validate_margin <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      margin <- as.integer(margin)
      stopifnot(identical(margin, 1L) || identical(margin, 2L))
    }),
    
    error = function (e) 
      stop(e$message, '\n`margin` must be 1 or 2.')
  )
}


validate_norm <- function (env = parent.frame()) {
  with(env, {
    
    if (!is.null(norm))
      norm <- switch(
        EXPR = match.arg(
          arg     = tolower(norm), 
          choices = c('none', 'percent', 'chord', 'binary', 'clr') ),
        'none'    = NULL,
        'percent' = NORM_PERCENT,
        'chord'   = NORM_CHORD,
        'binary'  = NORM_BINARY,
        'clr'     = NORM_CLR )
    
  })
}


validate_seed <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      stopifnot(is.numeric(seed))
      stopifnot(length(seed) == 1)
      stopifnot(!is.na(seed))
      stopifnot(seed >= 2**31 * -1)
      stopifnot(seed <= 2**31 - 1)
      stopifnot(seed %% 1 == 0)
      
      if (!is.integer(seed))
        seed <- as.integer(seed)
    }),
    
    error = function (e) 
      stop(e$message, '\n`seed` must be an integer between -2147483648 and 2147483647.')
  )
}


validate_times <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      if (!is.null(times)) {
        
        stopifnot(is.numeric(times))
        stopifnot(length(times) == 1)
        stopifnot(!is.na(times))
        stopifnot(times >= 0)
        stopifnot(times %% 1 == 0)
        
        if (!is.integer(times))
          times <- as.integer(times)
      }
    }),
    
    error = function (e) 
      stop(e$message, '\n`times` must be an integer greater than 0.')
  )
}


validate_tree <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      stopifnot(inherits(tree, 'phylo'))
      
      stopifnot(hasName(tree, 'edge'))
      stopifnot(is.matrix(tree$edge))
      stopifnot(ncol(tree$edge) == 2)
      if (typeof(tree$edge) != 'integer')
        tree$edge <- matrix(
          data     = as.integer(tree$edge), 
          nrow     = nrow(tree$edge), 
          ncol     = ncol(tree$edge),
          dimnames = dimnames(tree$edge) )
      
      stopifnot(hasName(tree, 'edge.length'))
      if (typeof(tree$edge.length) != 'double')
        tree$edge.length <- as.numeric(tree$edge.length)
      
      stopifnot(hasName(tree, 'tip.label'))
      
      if (margin == 1L) {
        
        stopifnot(!is.null(colnames(counts)))
        stopifnot(all(colnames(counts) %in% tree$tip.label))
        
        missing <- setdiff(tree$tip.label, colnames(counts))
        if (length(missing))
          counts <- cbind(
            counts, 
            matrix(
              data     = 0, 
              nrow     = nrow(counts),
              ncol     = length(missing), 
              dimnames = list(rownames(counts), missing) ))
        remove('missing')
        
        counts <- counts[,as.character(tree$tip.label),drop=FALSE]
      }
      else {
        
        stopifnot(!is.null(rownames(counts)))
        stopifnot(all(rownames(counts) %in% tree$tip.label))
        
        missing <- setdiff(tree$tip.label, rownames(counts))
        if (length(missing))
          counts <- rbind(
            counts, 
            matrix(
              data     = 0, 
              nrow     = length(missing), 
              ncol     = ncol(counts),
              dimnames = list(missing, rownames(counts)) ))
        remove('missing')
        
        counts <- counts[as.character(tree$tip.label),,drop=FALSE]
      }
    }),
    
    error = function (e) 
      stop(e$message, '\n`tree` is not a valid phylo object.')
  )
}


validate_underscores <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      
      if (!is.logical(underscores))
        underscores <- as.logical(underscores)
      
      stopifnot(length(underscores) == 1)
      stopifnot(!is.na(underscores))
    }),
    
    error = function (e) 
      stop(e$message, '\n`underscores` must be TRUE or FALSE.')
  )
}


validate_warn <- function (env = parent.frame()) {
  tryCatch(
    with(env, {
      stopifnot(identical(warn, TRUE) || identical(warn, FALSE))
    }),
    error = function (e) 
      stop(e$message, '\n`warn` must be either TRUE or FALSE.')
  )
}




assert_integer_counts <- function (env = parent.frame()) {
  with(env, {
    
    all_ints <- switch(
      get_matrix_package(counts),
      'base'   = all(counts   %% 1 == 0),
      'slam'   = all(counts$v %% 1 == 0),
      'Matrix' = all(counts@x %% 1 == 0) )
    
    if (!isTRUE(all_ints))
      stop('`counts` must be whole numbers (integers).')
    
    remove('all_ints')
  })
}

get_matrix_package <- function (counts) {
  
  if (is.matrix(counts)) {
    return ('base')
  } else if (inherits(counts, 'simple_triplet_matrix')) {
    return ('slam')
  } else {
    return ('Matrix')
  }
}

