
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
      stop('`counts` must be a valid numeric matrix.\n', e$message)
  )
}


validate_pairs <- function (env) {
  
  with(env, {
    if (ncol(counts) < 2)
      stop('`counts` must have at least two samples.')
  })
  
  tryCatch(
    with(env, {
      n <- ncol(counts)
      n <- n * (n - 1) / 2
      
      if (is.null(pairs)) {
        pairs <- as.integer(seq_len(n) - 1)
      }
      else {
        
        if (is.function(pairs))
          pairs <- local({
            m <- combn(n, 2)
            mapply(pairs, m[1,], m[2,])
          })
        
        if (is.logical(pairs)) {
          stopifnot(exprObject = bquote(length(pairs) == .(n)))
          pairs <- which(pairs)
        }
        else if (is.numeric(pairs)) {
          if (any(pairs %% 1 > 0)) stop('non-integer values')
          stopifnot(all(pairs > 0))
          stopifnot(exprObject = bquote(all(pairs <= .(n))))
          pairs <- sort(unique(as.integer(pairs)))
        }
        else {
          stop('cannot be ', typeof(pairs))
        }
        
        pairs <- pairs - 1L
      }
      remove('n')
    }),
    
    error = function (e) 
      stop('`pairs` must be a numeric or logical vector.\n', e$message)
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
      stop('`weighted` must be TRUE or FALSE.\n', e$message)
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
      stop('`newick` must be a character string.\n', e$message)
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
      stop('`underscores` must be TRUE or FALSE.\n', e$message)
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
      stop('`cpus` must be an integer greater than 0.\n', e$message)
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
      stop('`seed` must be an integer between -2147483648 and 2147483647.\n', e$message)
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
      stop('`times` must be an integer greater than 0.\n', e$message)
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
      stop('`depth` must be a positive integer or be between 0 and 1.\n', e$message)
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
      stop('`n_samples` must be an integer or be between 0 and 1.\n', e$message)
  )
}


validate_alpha <- function (env) {
  tryCatch(
    with(env, {
      
      if (class(alpha) != 'numeric')
        alpha <- as.numeric(alpha)
      
      stopifnot(length(alpha) == 1)
      stopifnot(!is.na(alpha))
      stopifnot(alpha >= 0 && alpha <= 1)
    }),
    
    error = function (e) 
      stop('`alpha` must be a single number between 0 and 1.\n', e$message)
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
      stopifnot(all(tree$tip.label %in% rownames(counts)))
      counts <- counts[as.character(tree$tip.label),]
    }),
    
    error = function (e) 
      stop('`tree` is not a valid phylo object.\n', e$message)
  )
}
