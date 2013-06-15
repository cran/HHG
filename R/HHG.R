# Main R file for the Heller-Heller-Gorfine test package

.NO_TIES_TEST     = as.integer(2)
.GENERAL_TEST     = as.integer(3)
.UDF_SPR_OBS      = as.integer(4)
.UDF_SPR_ALL      = as.integer(5)
.UDF_PPR_22_OBS   = as.integer(6)
.UDF_PPR_22_ALL   = as.integer(7)
.UDF_PPR_33_OBS   = as.integer(8)
.UDF_PPR_33_ALL   = as.integer(9)
.UDF_TPR_OBS      = as.integer(10)
.UDF_TPR_ALL      = as.integer(11)
.UDF_SPPR_OBS     = as.integer(12)
.UDF_SPPR_ALL     = as.integer(13)
.UDF_DDP_OBS      = as.integer(14)
.UDF_DDP_ALL      = as.integer(15)

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
  if (!is.sequential) {
    total.nr.tests = 1 # just so that we set Wald parameters to *something* so we can pass them on to C
  }
  
  if (!is.null(total.nr.tests)) {
    if (total.nr.tests < 1) {
      stop('total.nr.tests should be at least 1')
    }
    
    alpha.hyp = 0.05 / max(1, log(total.nr.tests))
    alpha0 = 0.05
    beta0 = min(0.01, 0.05 / total.nr.tests)
    eps = 0.01
  } else if (any(is.null(c(alpha.hyp, alpha0, beta0, eps)))) {
    stop('Either total.nr.tests or all of {alpha.hyp, alpha0, beta0, eps} must be specified (i.e. be non-null).')
  }

  alpha_hyp = as.double(alpha.hyp)
  alpha0 = as.double(alpha0)
  beta0 = as.double(beta0)
  eps = as.double(eps)

  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG', test_type, Dx, Dy, dummy.y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  
  if (tables.wanted) {
    n = nrow(Dx)
    tbls = res[12 + (1:(4 * n^2))]
    tbls[tbls == -42] = NA
    tbls.df = data.frame(matrix(tbls, nrow = n^2, ncol = 4))
    names(tbls.df) = c('A11', 'A12', 'A21', 'A22')
  } else {
    tbls.df = NULL
  }
  
  if (perm.stats.wanted) {
    n = nrow(Dx)
    ps = data.frame(matrix(res[12 + 4 * n^2 * tables.wanted + (1:(6 * nr.perm))], nrow = nr.perm, ncol = 6))
    names(ps) = c('sc', 'sl', 'mc', 'ml', 'ht', 'edist')
  } else {
    ps = NULL
  }
  
  return (list(sum.chisq = res[7], sum.lr = res[8], 
               max.chisq = res[9], max.lr = res[10],
               perm.pval.hhg.sc = res[1], perm.pval.hhg.sl = res[2], 
               perm.pval.hhg.mc = res[3], perm.pval.hhg.ml = res[4],
               extras.edist = NA, extras.ht = res[11],
               extras.perm.pval.edist = NA, extras.perm.pval.ht = res[5],
               extras.hhg.tbls = tbls.df, extras.perm.stats = ps))
}

# The univariate distribution-free test
xdp.test = function(x, y, variant = 'ppr.33.obs', K = 3, correct.mi.bias = F) {
  res = .hhg.test.udfree(x = x, y = y, variant = variant, K = K, correct.mi.bias = correct.mi.bias)
  return (list(sum.chisq = res$sum.chisq, sum.lr = res$sum.lr, max.chisq = res$max.chisq, max.lr = res$max.lr))
}

# (I'm keeping the old interface because some big simulations rely on it)
.hhg.test.udfree = function(x, y, variant = 'ppr.33.obs', nr.perm = 0, K = 3, correct.mi.bias = F,
                           total.nr.tests = 1, is.sequential = T, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
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

  # Dy is used to store x.sorted.by.y (first column), and y.sorted.by.x (second column)
  if (variant == 'ddp.obs' || variant == 'ddp.all') {
    Dy = cbind(Dx[order(y)], y[order(x)])
  } else {
    Dx = Dx - 1
    Dy = 0
    y = y - 1
  }
  
  w.sum = 0
  w.max = 2

  extra_params = as.double(c(K, correct.mi.bias))
  
  is_sequential = as.integer(is.sequential)
  if (!is.sequential) {
    total.nr.tests = 1 # just so that we set Wald parameters to *something* so we can pass them on to C
  }
  
  if (!is.null(total.nr.tests)) {
    if (total.nr.tests < 1) {
      stop('total.nr.tests should be at least 1')
    }
    
    alpha_hyp = 0.05 / max(1, log(total.nr.tests))
    alpha0 = 0.05
    beta0 = min(0.01, 0.05 / total.nr.tests)
    eps = 0.01
  } else if (any(is.null(c(alpha.hyp, alpha0, beta0, eps)))) {
    stop('Either total.nr.tests or all of {alpha.hyp, alpha0, beta0, eps} must be specified (i.e. be non-null).')
  }

  alpha_hyp = as.double(alpha.hyp)
  alpha0 = as.double(alpha0)
  beta0 = as.double(beta0)
  eps = as.double(eps)

  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG', test_type, Dx, Dy, y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)

  if (tables.wanted) {
    n = nrow(Dx)
    tbls = res[12 + (1:(4 * n^2))]
    tbls[tbls == -42] = NA
    tbls.df = data.frame(matrix(tbls, nrow = n^2, ncol = 4))
    names(tbls.df) = c('A11', 'A12', 'A21', 'A22')
  } else {
    tbls.df = NULL
  }
  
  if (perm.stats.wanted) {
    n = nrow(Dx)
    ps = data.frame(matrix(res[12 + 4 * n^2 * tables.wanted + (1:(6 * nr.perm))], nrow = nr.perm, ncol = 6))
    names(ps) = c('sc', 'sl', 'mc', 'ml', 'ht', 'edist')
  } else {
    ps = NULL
  }
  
  return (list(sum.chisq = res[7], sum.lr = res[8], 
               max.chisq = res[9], max.lr = res[10],
               perm.pval.hhg.sc = res[1], perm.pval.hhg.sl = res[2], 
               perm.pval.hhg.mc = res[3], perm.pval.hhg.ml = res[4],
               extras.edist = res[12], extras.ht = res[11],
               extras.perm.pval.edist = res[6], extras.perm.pval.ht = res[5],
               extras.hhg.tbls = tbls.df, extras.perm.stats = ps))
}
