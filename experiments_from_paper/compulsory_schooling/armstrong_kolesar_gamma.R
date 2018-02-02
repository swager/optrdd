ak.tau = function(X, max.second.derivative, Y = NULL, weights = rep(1, length(X)), threshold = 0, sigma.sq = NULL, change.derivative = TRUE, alpha = 0.95, max.window = max(abs(X - threshold)), num.bucket = 200) {
  
  # Naive initialization for sigma.sq if needed
  if (is.null(sigma.sq)) {
    if (is.null(Y)) {
      warning("Setting noise level to 1 as default...")
      sigma.sq = 1
    } else {
      Y.bar = sum(Y * weights) / sum(weights)
      sigma.sq = sum((Y - Y.bar)^2 * weights^2) / sum(weights)
    }
  }
  
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
  X.counts = as.numeric(bucket.map %*% weights[inrange])
  
  
  # only consider counts that occurs inside the bandwidth
  realized.idx = which(X.counts > 0)
  num.realized = length(realized.idx)
  
  # This optimizer learns bucket-wise gammas. Let k denote
  # the bucket index, n[k] the number of observations in
  # bucket k, and x[k] is the center of bucket k.
  #
  # We solve the following. Note that gamma[k] must be 0
  # if n[k] is 0, so we only optimize gamma over realized indices.
  #
  # argmin sum_k (gamma_+[k] - gamma_-[k])^2 * n[k] + B^2 I^2
  # subject to:
  #  sum_k n[k] (gamma_+[k] - gamma_-[k]) = 0
  #  sum_k n[k] (gamma_+[k] - gamma_-[k]) (2 W[k] - 1) = 2
  #  sum_k n[k] (gamma_+[k] - gamma_-[k]) x[k] = 0
  #  sum_k n[k] (gamma_+[k] - gamma_-[k]) (x[k])_+ = 0
  #  sum_k n[k] (gamma_+[k] + gamma_-[k]) x[k]^2 / 2 = I
  #  gamma_+, gamma_- >= 0
  
  Dmat =diag(c(sigma.sq * X.counts[realized.idx],
               sigma.sq * X.counts[realized.idx],
               max.second.derivative^2))
  dvec = rep(0, 2 * num.realized + 1)
  
  if(!change.derivative) {
    stop("Only implemented with derivate change at threshold.")
  }
  
  Amat = cbind(c(X.counts[realized.idx], -X.counts[realized.idx], 0),
               c(X.counts[realized.idx] * sign(xx[realized.idx]),
                 -X.counts[realized.idx] * sign(xx[realized.idx]), 0),
               c(X.counts[realized.idx] * xx[realized.idx],
                 -X.counts[realized.idx] * xx[realized.idx], 0),
               c(X.counts[realized.idx] *  pmax(xx[realized.idx], 0),
                 -X.counts[realized.idx] *  pmax(xx[realized.idx], 0), 0),
               c(X.counts[realized.idx] * xx[realized.idx]^2/2,
                 X.counts[realized.idx] * xx[realized.idx]^2/2, -1),
               diag(rep(1, 2 * num.realized + 1)))
  
  bvec = c(0, 2, 0, 0, 0, rep(rep(0, 2 * num.realized + 1)))
  meq = 5
  
  soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
  
  gamma.xx = rep(0, num.bucket)
  gamma.xx[realized.idx] = soln[1:num.realized] - soln[num.realized + (1:num.realized)]
  
  # Now map this x-wise function into a weight for each observation
  gamma = rep(0, length(X))
  gamma[inrange] = weights[inrange] * as.numeric(Matrix::t(bucket.map) %*% gamma.xx)
  
  max.bias = max.second.derivative * soln[1 + 2 * num.realized]
  
  # If outcomes are provided, also compute confidence intervals for tau.
  if (!is.null(Y)) {
    
    # The point estimate
    tau.hat = sum(gamma * Y)
    
    # A heteroskedaticity-robust variance estimate
    regr.df = data.frame(X=X, W=X>=threshold, Y=Y)
    Y.fit = lm(Y ~ X * W, data = regr.df[inrange,], weights=weights[inrange])
    Y.resid.sq = rep(0, length(Y))
    Y.resid.sq[inrange] = (Y[inrange] - predict(Y.fit))^2 * sum(weights[inrange]) / (sum(weights[inrange]) - 4)
    se.hat.tau = sqrt(sum(Y.resid.sq * gamma^2))
    
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
             gamma=gamma,
             gamma.fun = data.frame(xx=xx[realized.idx] + threshold,
                                    gamma=gamma.xx[realized.idx]))
  class(ret) = "optrdd"
  return(ret)
}