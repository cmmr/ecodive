
validate_args <- function () {
  env  <- parent.frame()
  args <- ls(env)
  for (arg in args[order(args != 'counts')])
    do.call(paste0('validate_', arg), list(env))
}


validate_counts <- function (env) {
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
      
      # Pull counts matrix from complex counts object.
      if (!is.matrix(counts)) {
        
        if (inherits(counts, 'phyloseq')) {
          counts <- counts@otu_table
        }
        else if (inherits(counts, 'rbiom')) {
          counts <- counts$counts
        }
        else if (inherits(counts, 'TreeSummarizedExperiment')) {
          counts <- counts@assays@data[[1]]
        }
        else if (inherits(counts, 'SummarizedExperiment')) {
          counts <- counts@assays@data[[1]]
        }
        
        counts <- as.matrix(counts)
      }
      
      stopifnot(is.numeric(counts))
      stopifnot(length(dim(counts)) == 2)
      stopifnot(all(counts >= 0))
      stopifnot(nrow(counts) > 0)
      stopifnot(ncol(counts) > 0)
      
      if (typeof(counts) != 'double')
        counts <- matrix(
          data     = as.numeric(counts), 
          nrow     = nrow(counts), 
          ncol     = ncol(counts),
          dimnames = dimnames(counts) )
    }),
    
    error = function (e) 
      stop(e$message, '\n`counts` must be a valid numeric matrix.')
  )
}


validate_pairs <- function (env) {
  
  with(env, {
    if (ncol(counts) < 2)
      stop('`counts` must have at least two samples.')
  })
  
  tryCatch(
    with(env, {
      n_samples   <- ncol(counts)
      n_distances <- n_samples * (n_samples - 1) / 2
      
      if (is.null(pairs)) {
        pairs <- as.integer(seq_len(n_distances) - 1)
      }
      else {
        
        if (is.function(pairs))
          pairs <- local({
            m <- combn(n_samples, 2)
            mapply(pairs, m[1,], m[2,])
          })
        
        if (is.logical(pairs)) {
          stopifnot(exprObject = bquote(length(pairs) == .(n_distances)))
          pairs <- which(pairs)
        }
        else if (is.numeric(pairs)) {
          if (any(pairs %% 1 > 0)) stop('non-integer values')
          if (!all(pairs >= 1 & pairs <= n_distances))
            stop('expected `pairs` values between 1 and ', n_distances)
          pairs <- sort(unique(as.integer(pairs)))
        }
        else {
          stop('cannot be ', typeof(pairs))
        }
        
        pairs <- pairs - 1L
      }
      remove('n_samples', 'n_distances')
    }),
    
    error = function (e) 
      stop(e$message, '\n`pairs` must be a numeric or logical vector.')
  )
}


validate_weighted <- function (env) {
  tryCatch(
    with(env, {
      
      if (!is.logical(weighted))
        weighted <- as.logical(weighted)
      
      stopifnot(length(weighted) == 1)
      stopifnot(!is.na(weighted))
    }),
    
    error = function (e) 
      stop(e$message, '\n`weighted` must be TRUE or FALSE.')
  )
}



validate_newick <- function (env) {
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


validate_underscores <- function (env) {
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


validate_cpus <- function (env) {
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


validate_seed <- function (env) {
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


validate_times <- function (env) {
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


validate_depth <- function (env) {
  tryCatch(
    with(env, {
      
      stopifnot(is.numeric(depth))
      stopifnot(length(depth) == 1)
      stopifnot(!is.na(depth))
      stopifnot(depth > 0)
      stopifnot(depth %% 1 == 0 || depth < 1)
      
      if (depth %% 1 == 0) {
        stopifnot(depth <= max(colSums(counts)))
        depth <- as.integer(depth)
      }
      
    }),
    
    error = function (e) 
      stop(e$message, '\n`depth` must be a positive integer or be between 0 and 1.')
  )
}


validate_n_samples <- function (env) {
  tryCatch(
    with(env, {
      
      if (!is.null(n_samples)) {
      
        stopifnot(is.numeric(n_samples))
        stopifnot(length(n_samples) == 1)
        stopifnot(!is.na(n_samples))
        
        if (n_samples %% 1 == 0) {
          stopifnot(n_samples <= ncol(counts))
          stopifnot(n_samples > -ncol(counts))
          n_samples <- as.integer(n_samples)
        }
        else {
          stopifnot(n_samples > 0 && n_samples < 1)
        }
      }
      
    }),
    
    error = function (e) 
      stop(e$message, '\n`n_samples` must be an integer or be between 0 and 1.')
  )
}


validate_alpha <- function (env) {
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


validate_tree <- function (env) {
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
            dimnames = list(missing, colnames(counts)) ))
      remove('missing')
      
      counts <- counts[as.character(tree$tip.label),,drop=FALSE]
    }),
    
    error = function (e) 
      stop(e$message, '\n`tree` is not a valid phylo object.')
  )
}
