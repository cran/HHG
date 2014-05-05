# This file contains various constants and utility functions

.TWO_SAMPLE_TEST   = as.integer(0)
.K_SAMPLE_TEST     = as.integer(1)
.NO_TIES_TEST      = as.integer(2)
.GENERAL_TEST      = as.integer(3)
.UDF_SPR_OBS       = as.integer(4)
.UDF_SPR_ALL       = as.integer(5)
.UDF_PPR_22_OBS    = as.integer(6)
.UDF_PPR_22_ALL    = as.integer(7)
.UDF_PPR_33_OBS    = as.integer(8)
.UDF_PPR_33_ALL    = as.integer(9)
.UDF_TPR_OBS       = as.integer(10)
.UDF_TPR_ALL       = as.integer(11)
.UDF_SPPR_OBS      = as.integer(12)
.UDF_SPPR_ALL      = as.integer(13)
.UDF_DDP_OBS       = as.integer(14)
.UDF_DDP_ALL       = as.integer(15)


.configure.wald.sequential = function(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps) {
  if (!is.sequential) {
    total.nr.tests = 1 # just so that we set Wald parameters to *something* so we can pass them on to C (but their values will be ignored)
  }
  
  if (!is.null(total.nr.tests)) {
    if (total.nr.tests < 1) {
      stop('seq.total.nr.tests should be at least 1')
    }
    
    if (any(c(!is.null(alpha.hyp), !is.null(alpha0), !is.null(beta0), !is.null(eps)))) {
      warning('Wald sequential parameters (seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps) provided by user but are overwritten by defaults. See documentation on when and how to specify these arguments.')
    }
    
    alpha.hyp = 0.05 / max(1, log(total.nr.tests))
    alpha0 = 0.05
    beta0 = min(0.01, 0.05 / total.nr.tests)
    eps = 0.01
  } else if (any(c(is.null(alpha.hyp), is.null(alpha0), is.null(beta0), is.null(eps)))) {
    stop('Either seq.total.nr.tests or all of {seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps} must be specified (i.e. be non-null).')
  }
  
  return (list(alpha.hyp, alpha0, beta0, eps))
}

# Now organize the results produced by our C implementation into a convenient R structure
.organize.results = function(res, n, nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F) {
  has.grid = (grid.len > 0)
  variants = c('sc', 'sl', 'mc', 'ml')
  grid.variants = paste(variants, '.grid', sep = '')
  nr.variants = length(variants)
  base.nr.stats = 6
  total.nr.stats = ifelse(has.grid, base.nr.stats + nr.variants * grid.len, base.nr.stats)
  res.pvals.offset = 0
  res.obs.stats.offest = total.nr.stats
  res.tables.offset = total.nr.stats * 2
  res.tables.size = 4 * n^2 * tables.wanted
  res.perm.stats.offset = res.tables.offset + res.tables.size
  res.perm.stats.size = total.nr.stats * nr.perm * perm.stats.wanted
  
  ret = list()
  
  ret$sum.chisq = res[1 + res.obs.stats.offest]
  ret$sum.lr    = res[2 + res.obs.stats.offest]
  ret$max.chisq = res[3 + res.obs.stats.offest]
  ret$max.lr    = res[4 + res.obs.stats.offest]
  
  if (nr.perm > 0) {
    ret$perm.pval.hhg.sc = res[1]
    ret$perm.pval.hhg.sl = res[2]
    ret$perm.pval.hhg.mc = res[3]
    ret$perm.pval.hhg.ml = res[4]
  }
    
  if (extra.stats.wanted) {
    # NOTE: these are actually only provided for the 2-sample test, for now
    ret$extras.edist = res[6 + res.obs.stats.offest]
    ret$extras.ht = res[5 + res.obs.stats.offest]
    
    if (nr.perm > 0) {
      ret$extras.perm.pval.edist = res[6]
      ret$extras.perm.pval.ht = res[5]
    }
  }
  
  if (has.grid) {
    ret$obs.stats.grid = res[total.nr.stats + base.nr.stats + 1:(nr.variants * grid.len)]
    names(ret$obs.stats.grid) = paste(grid.variants, rep(1:grid.len, each = nr.variants), sep = '')
  }
  
  if (tables.wanted) {
    # FIXME: does this make sense in this context? we'll see
    tbls = res[res.tables.offset + (1:res.tables.size)]
    tbls[tbls == -42] = NA # ahem...
    tbls.df = data.frame(matrix(tbls, nrow = n^2, ncol = 4))
    names(tbls.df) = c('A11', 'A12', 'A21', 'A22')
    
    ret$extras.hhg.tbls = tbls.df
  }
  
  if ((nr.perm > 0) && perm.stats.wanted) {
    ps = data.frame(matrix(res[res.perm.stats.offset + (1:(nr.variants * nr.perm))], nrow = nr.perm))
    names(ps) = variants
    
    if (has.grid) {
      ofst = res.perm.stats.offset + base.nr.stats * nr.perm
      ps.grid = data.frame(matrix(res[ofst + (1:(nr.variants * grid.len * nr.perm))], nrow = nr.perm))
      names(ps.grid) = paste(grid.variants, rep(1:grid.len, each = nr.variants), sep = '')
      ps = cbind(ps, ps.grid)
    }

    ret$extras.perm.stats = ps
  }
  
  return (ret)
}
