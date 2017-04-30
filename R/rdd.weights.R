optrdd = function(X, max.second.derivative, Y = NULL, num.samples = rep(1, length(X)), threshold = 0, sigma.sq = NULL, change.derivative = TRUE, alpha = 0.95, lambda.mult = 1, max.window = max(abs(X - threshold)), num.bucket = 200, use.homoskedatic.variance = FALSE) {

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
  
  realized.idx = which(X.counts > 0)
  num.realized = length(realized.idx)
  
  # The matrix M is used to intergrate a function over xx,
  # starting at 0. The matrix M2 integrates twice.
  M = outer(xx, xx, FUN=Vectorize(function(x1, x2) {
    if(x1 < 0 & x2 <= 0 & x1 < x2) {return(-bin.width)}
    if(x1 > 0 & x2 >= 0 & x1 > x2) {return(bin.width)}
    return(0)
  }))
  M2 = M %*% M
  
  # Given homoskedatisc errors with variance sigma.sq, and a bound
  # max.derivative on the second derivate, this choice of lambda minimizes
  # the worst-case MSE of the estimator for tau.
  # 
  # The factor lambda.mult can be used to tune the parameter choice (e.g.,
  # to make the CIs as short as possible).
  lambda = lambda.mult * max.second.derivative^2 / sigma.sq
  
  # This optimizer learns bucket-wise gammas. Let k denote
  # the bucket index, n[k] the number of observations in
  # bucket k, and x[k] is the center of bucket k.
  # W also have positive dummy variables such that
  # nu_+ + nu_- = abs(t(M2) %*% (n * gamma)); note that
  # the maximum error due to curvature is
  # ||t(M2) %*% (n * gamma)||_1.
  #
  # We solve the following. Note that gamma[k] must be 0
  # if n[k] is 0, so we only optimize gamma over realized indices.
  #
  # argmin sum_k gamma[k]^2 * n[k] + lambda z^2
  # subject to:
  #  sum_k n[k] gamma[k] = 0
  #  sum_k n[k] gamma[k] (2 W[k] - 1) = 2
  #  sum_k n[k] gamma[k] x[k] = 0
  #  sum_k nu_+[k] + nu_-[k] = z
  #  - t(M2) %*% (n * gamma) + nu_+ >= 0 (elem. wise)
  #  t(M2) %*% (n * gamma)[k] + nu_- >= 0 (elem. wise)
  #  nu_+, nu_- >= 0
  
  penalty.mat = diag(X.counts)[realized.idx,] %*% M2
  Dmat = diag(c(X.counts[realized.idx], rep(0.00000000001 * lambda, 2*num.bucket), lambda))
  dvec = rep(0, num.realized + 2* num.bucket + 1)
  Amat = cbind(c(X.counts[realized.idx], rep(0, 2*num.bucket + 1)),
               c(X.counts[realized.idx] * sign(xx[realized.idx]), rep(0, 2*num.bucket + 1)),
               c(X.counts[realized.idx] * xx[realized.idx], rep(0, 2*num.bucket + 1)),
               c(rep(0, num.realized), rep(1, 2*num.bucket), -1),
               rbind(-penalty.mat, diag(1, num.bucket), diag(0, num.bucket), 0),
               rbind(penalty.mat, diag(0, num.bucket), diag(1, num.bucket), 0),
               rbind(matrix(0, num.realized, num.bucket), diag(1, num.bucket), diag(0, num.bucket), 0),
               rbind(matrix(0, num.realized, num.bucket), diag(0, num.bucket), diag(1, num.bucket), 0))
  bvec = c(0, 2, 0, 0, rep(0, 4 * num.bucket))
  meq = 4
  
  # If we want to identify the CATE at the thresold, we also need to
  # add the constraint: sum_{k : x[k] >= 0} n[k] gamma[k] x[k] = 0
  if(change.derivative) {
    Amat = cbind(c(X.counts[realized.idx] * pmax(xx[realized.idx], 0), rep(0, 2*num.bucket + 1)), Amat)
    bvec = c(0, bvec)
    meq = meq + 1
  }
  
  # Solve the QP...
  soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
  
  # Extract the weighting function gamma(x), as a function of x
  gamma.xx = rep(0, num.bucket)
  gamma.xx[realized.idx] = soln[1:num.realized]
  
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
              gamma=gamma,
              gamma.fun = data.frame(xx=xx[realized.idx] + threshold,
                                     gamma=gamma.xx[realized.idx]))
  class(ret) = "optrdd"
  return(ret)
}

