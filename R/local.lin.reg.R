llr = function(X, max.second.derivative, bandwidth = NULL, Y = NULL, num.samples = rep(1, length(X)), threshold = 0, sigma.sq = NULL, change.derivative = TRUE, alpha = 0.95, max.window = max(abs(X - threshold)), num.bucket = 200, kernel = c("rectangular", "triangular"), minimization.target = c("mse", "ci.length"), use.homoskedatic.variance = FALSE) {
  
  kernel = match.arg(kernel)
  minimization.target = match.arg(minimization.target)
  
  # We compute our estimator based on a histogram summary of the data,
  # shifted such that the threshold is at 0. The breaks vector defines
  # the boundaries between buckets.
  xx = seq(-max.window, max.window, length.out = num.bucket)
  bin.width = xx[2] - xx[1]
  breaks = c(xx - bin.width/2, max(xx) + bin.width/2)
  
  # Construct a (weighted) histogram representing the X-values.
  inrange = which(abs(X - threshold) / max.window <= 1)
  bucket = cut(X[inrange] - threshold, breaks = breaks)
  bucket.map = Matrix::sparse.model.matrix(~bucket + 0, transpose = TRUE)
  X.counts = as.numeric(bucket.map %*% num.samples[inrange])
  
  # Actual mean of the data in the bucket; use this to make balance sharp
  X.mids = bucket.map %*% (X[inrange] * num.samples[inrange]) /
    X.counts
  
  # Naive initialization for sigma.sq if needed
  if (is.null(sigma.sq)) {
    if (is.null(Y)) {
      warning("Setting noise level to 1 as default...")
      sigma.sq = 1
    } else {
      regr.df = data.frame(X=X, W=X>=threshold, Y=Y)
      Y.fit = lm(Y ~ X * W, data = regr.df[inrange,], weights=num.samples[inrange])
      sigma.sq = sum((Y[inrange] - predict(Y.fit))^2 * num.samples[inrange]^2) /
        (sum(num.samples[inrange]) - 4)
    }
  }
  
  # The matrix M is used to intergrate a function over xx,
  # starting at 0. The matrix M2 integrates twice.
  M = outer(xx, xx, FUN=Vectorize(function(x1, x2) {
    if(x1 < 0 & x2 <= 0 & x1 < x2) {return(-bin.width)}
    if(x1 > 0 & x2 >= 0 & x1 > x2) {return(bin.width)}
    return(0)
  }))
  M2 = M %*% M
  
  if (!is.null(bandwidth)) {
    bw.vec = bandwidth
  } else {
    bw.vec = c(1:40) * max.window / 40
  }
  
  soln.vec = lapply(bw.vec, function(bw) {
    
    # only consider counts that occurs inside the bandwidth
    realized.idx = which((X.counts > 0) & (abs(xx) < bw))
    num.realized = length(realized.idx)
    
    signed.num.realized = min(sum(xx[realized.idx] > 0), sum(xx[realized.idx] < 0))
    
    if (signed.num.realized < 2) {
      return(list(max.mse=NA, max.bias=NA, homosk.plusminus=NA, gamma.xx=NA, realized.idx=NA))
    }
    
    # This optimizer learns bucket-wise gammas. Let k denote
    # the bucket index, n[k] the number of observations in
    # bucket k, and x[k] is the center of bucket k.
    #
    # We solve the following. Note that gamma[k] must be 0
    # if n[k] is 0, so we only optimize gamma over realized indices.
    #
    # argmin sum_k gamma[k]^2 * n[k]
    # subject to:
    #  sum_k n[k] gamma[k] = 0
    #  sum_k n[k] gamma[k] (2 W[k] - 1) = 2
    #  sum_k n[k] gamma[k] x[k] = 0
    
    if (kernel == "rectangular") {
      Dmat = diag(X.counts[realized.idx])
    } else if (kernel == "triangular") {
      Dmat = diag(X.counts[realized.idx] / (1 - abs(xx[realized.idx]) / bw))
    }
    dvec = rep(0, num.realized)
    
    Amat = cbind(X.counts[realized.idx],
                 X.counts[realized.idx] * sign(xx[realized.idx]),
                 X.counts[realized.idx] * X.mids[realized.idx])
    
    if(!change.derivative) {
      bvec = c(0, 2, 0)
      meq = 3
    } else {
      Amat = cbind(Amat,
                   X.counts[realized.idx] * pmax(X.mids[realized.idx], 0))
      bvec = c(0, 2, 0, 0)
      meq = 4
    }
    
    soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
    
    gamma.xx = rep(0, num.bucket)
    gamma.xx[realized.idx] = soln[1:num.realized]
    max.bias = max.second.derivative * sum(abs(t(M2) %*% (X.counts * gamma.xx)))
    sigma.hat.homosk = sqrt(sigma.sq * sum(X.counts[realized.idx] * gamma.xx[realized.idx]^2))
    max.mse = max.bias^2 + sigma.hat.homosk^2
    homosk.plusminus = get.plusminus(max.bias, sigma.hat.homosk, alpha)
    
    return(list(max.mse=max.mse,
               max.bias=max.bias,
               homosk.plusminus=homosk.plusminus,
               gamma.xx=gamma.xx,
               realized.idx=realized.idx))
  })
  
  # pick out the best soln
  if(minimization.target == "mse") {
    max.mse = unlist(sapply(soln.vec, function(vv) vv$max.mse))
    opt.idx = which.min(max.mse)
  } else {
    plusmin = unlist(sapply(soln.vec, function(vv) vv$homosk.plusminus))
    opt.idx = which.min(plusmin)
  }
  
  gamma.xx = soln.vec[[opt.idx]]$gamma.xx
  realized.idx = soln.vec[[opt.idx]]$realized.idx
  
  # Now map this x-wise function into a weight for each observation
  gamma = rep(0, length(X))
  gamma[inrange] = num.samples[inrange] * as.numeric(Matrix::t(bucket.map) %*% gamma.xx)
  
  # Compute the worst-case imbalance...
  max.bias = max.second.derivative * sum(abs(t(M2) %*% (X.counts * gamma.xx)))
  
  # If outcomes are provided, also compute confidence intervals for tau.
  if (!is.null(Y)) {
    
    # The point estimate
    tau.hat = sum(gamma * Y)
    
    if (use.homoskedatic.variance) {
      se.hat.tau = sqrt(sum(gamma^2 * sigma.sq / num.samples))
    } else {
      # A heteroskedaticity-robust variance estimate
      regr.df = data.frame(X=X, W=X>=threshold, Y=Y)
      Y.fit = lm(Y ~ X * W, data = regr.df[inrange,], weights=num.samples[inrange])
      Y.resid.sq = rep(0, length(Y))
      Y.resid.sq[inrange] = (Y[inrange] - predict(Y.fit))^2 *
        sum(num.samples[inrange]) / (sum(num.samples[inrange]) - 4)
      se.hat.tau = sqrt(sum(Y.resid.sq * gamma^2))
    }
    
    # Confidence intervals that account for both bias and variance
    tau.plusminus = get.plusminus(max.bias, se.hat.tau, alpha)
  } else {
    tau.hat = NULL
    se.hat.tau = sqrt(sigma.sq * sum(gamma^2))
    tau.plusminus = get.plusminus(max.bias, se.hat.tau, alpha)
  }
  
  ret = list(tau.hat=tau.hat,
             tau.plusminus=tau.plusminus,
             alpha=alpha,
             max.bias = max.bias,
             sampling.se=se.hat.tau,
             bandwidth=bw.vec[opt.idx],
             gamma=gamma,
             gamma.fun = data.frame(xx=xx[realized.idx] + threshold,
                                    gamma=gamma.xx[realized.idx]))
  class(ret) = "optrdd"
  return(ret)
}