# Main R file for the Heller-Heller-Gorfine test package

# The general test of independence (with or without handling of ties)
hhg.test = function(Dx, Dy, ties = T, w.sum = 0, w.max = 2, nr.perm = 10000, 
  total.nr.tests = 1, is.sequential = T, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
  nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.double(Dx) || !is.double(Dy) || !is.matrix(Dx) || !is.matrix(Dy) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != nrow(Dy) || nrow(Dy) != ncol(Dy)) {
    stop('Dx and Dy are expected to be square matrices of doubles, and must have the same number of rows/cols')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  
  if (ties) {
    test_type = .GENERAL_TEST
  } else {
    test_type = .NO_TIES_TEST
  }

  dummy.y = matrix(0, nrow(Dy), 1)
  extra_params = as.double(0)
  is_sequential = as.integer(is.sequential)

  wald = .configure.wald.sequential(total.nr.tests, is.sequential, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG', test_type, Dx, Dy, dummy.y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = T)
  return (ret)
}

# The univariate distribution-free test
xdp.test = function(x, y, variant = 'DDP', K = 3, correct.mi.bias = F, 
  w.sum = 0, w.max = 2) 
{
  if (variant == 'DDP') {
    if (K != as.integer(K) || K < 2) {
      stop('K must be an integer greater than 1')
    } else if (K == 2) {
      x.variant = 'spr.obs'
    } else if (K == 3) {
      x.variant = 'ppr.33.obs'
    } else if (K == 4) {
      x.variant = 'tpr.obs'
    } else {
      x.variant = 'ddp.obs'
    }
  } else if (variant == 'ADP') {
    if (K != as.integer(K) || K < 2) {
      stop('K must be an integer greater than 1')
    } else if (K == 2) {
      x.variant = 'spr.all'
    } else if (K == 3) {
      # One could use 'ppr.33.all' that has the same complexity, but in practice it is much slower
      x.variant = 'ddp.all'
    } else if (K == 4) {
      # tpr.all is too time consuming
      x.variant = 'ddp.all'
    } else {
      x.variant = 'ddp.all'
    }
  }

  ret = .hhg.test.udfree(x = x, y = y, variant = x.variant, K = K, correct.mi.bias = correct.mi.bias)
  
  if (((variant == 'DDP') && (K > 4)) || ((variant == 'ADP') && (K > 2))) {
    ret$max.chisq = NA
    ret$max.lr    = NA

    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
    
  return (ret)
}

# (I'm keeping the old interface because some big simulations rely on it)
.hhg.test.udfree = function(x, y, variant = 'ppr.33.obs', w.sum = 0, w.max = 2,
  nr.perm = 0, K = 3, correct.mi.bias = F, total.nr.tests = 1, 
  is.sequential = T, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
  nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.ordered(y)) {
    stop('y is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (K < 2 || K > length(x)) {
    stop('K should be strictly between 2 and length(x)')
  }

  if (variant == 'spr.obs') {
    test_type = .UDF_SPR_OBS
  } else if (variant == 'spr.all') {
    test_type = .UDF_SPR_ALL
  } else if (variant == 'ppr.22.obs') {
    test_type = .UDF_PPR_22_OBS
  } else if (variant == 'ppr.22.all') {
    test_type = .UDF_PPR_22_ALL
  } else if (variant == 'ppr.33.obs') {
    test_type = .UDF_PPR_33_OBS
  } else if (variant == 'ppr.33.all') {
    test_type = .UDF_PPR_33_ALL
  } else if (variant == 'tpr.obs') {
    test_type = .UDF_TPR_OBS
  } else if (variant == 'tpr.all') {
    test_type = .UDF_TPR_ALL
  } else if (variant == 'sppr.obs') {
    test_type = .UDF_SPPR_OBS
  } else if (variant == 'sppr.all') {
    test_type = .UDF_SPPR_ALL
  } else if (variant == 'ddp.obs') {
    test_type = .UDF_DDP_OBS
  } else if (variant == 'ddp.all') {
    test_type = .UDF_DDP_ALL
  } else {
  	stop('Unexpected variant specified.')
  }
  
  # Dx is used to store x
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  y = as.matrix(as.double(rank(y, ties.method = 'random')), nrow = length(y), ncol = 1)

  # For DDP, Dy is used to store x.sorted.by.y (first column), and y.sorted.by.x (second column)
  if (variant == 'ddp.obs' || variant == 'ddp.all') {
    Dy = cbind(Dx[order(y)], y[order(Dx)])
  } else {
    Dx = Dx - 1
    Dy = 0
    y = y - 1
  }
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)

  extra_params = as.double(c(K, correct.mi.bias))
  is_sequential = as.integer(is.sequential)

  wald = .configure.wald.sequential(total.nr.tests, is.sequential, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F)
  return (ret)
}
