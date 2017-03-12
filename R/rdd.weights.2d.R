optrdd.2d = function(X, max.second.derivative, Y = NULL, weights = NULL, threshold = c(0, 0), sigma.sq = NULL, change.derivative = TRUE, alpha = 0.95, lambda.mult = 1, max.window = c(max(abs(X[,1] - threshold[1])), max(abs(X[,2] - threshold[2]))), num.bucket = c(20, 20)) {
  
  if (ncol(X) != 2) { stop("The running variable must be bivariate") }
  
  if (!is.null(weights)) { stop("Weighted version of optrdd.2d not implemented.") }
  
  # Naive initialization for sigma.sq if needed
  if (is.null(sigma.sq)) {
    if (is.null(Y)) {
      warning("Setting noise level to 1 as default...")
      sigma.sq = 1
    } else {
      W = as.numeric((X[,1] > threshold[1]) & (X[,2] > threshold[2]))
      Y.hat = predict(lm(Y ~ X + W))
      sigma.sq = mean((Y - Y.hat)^2) * length(W) / (length(W) - 4)
    }
  }
  
  xx1 = seq(threshold[1] - max.window[1], threshold[1] + max.window[1], length.out = num.bucket[1])
  xx2 = seq(threshold[2] - max.window[2], threshold[2] + max.window[2], length.out = num.bucket[2])
  
  bin.width = c(xx1[2] - xx1[1], xx2[2] - xx2[1])
  xx12 = expand.grid(xx1, xx2)
  
  center.points = which(xx12[,1] %in% xx1[order(abs(xx1 - threshold[1]))[1:2]] &
                          xx12[,2] %in% xx2[order(abs(xx2 - threshold[2]))[1:2]])
  center.mat = matrix(0, length(center.points), nrow(xx12))
  for(iter in 1:length(center.points)) center.mat[iter, center.points[iter]] = 1
  
  # matrix that probes all local second derivatives, along axis and diagonals
  nabla = matrix(0, 4 * (num.bucket[1] - 2) * (num.bucket[2] - 2), nrow(xx12))
  curr.idx = 0
  for (i1 in 2:(num.bucket[1] - 1)) {
    for (i2 in 2:(num.bucket[2] - 1)) {
      nabla[curr.idx + 1,  (i1 - 1):(i1 + 1) + (i2 - 1) * num.bucket[1]] = c(1, -2, 1) / bin.width[1]^2
      nabla[curr.idx + 2, i1 + ((i2 - 2):i2) * num.bucket[1]] = c(1, -2, 1) / bin.width[2]^2
      nabla[curr.idx + 3,  (i1 - 1):(i1 + 1) + ((i2 - 2):i2) * num.bucket[1]] = c(1/2, -1, 1/2) / prod(bin.width)
      nabla[curr.idx + 4,  (i1 - 1):(i1 + 1) + (i2:(i2 - 2)) * num.bucket[1]] = c(1/2, -1, 1/2) / prod(bin.width)
      curr.idx = curr.idx + 4
    }
  }
  
  inrange = which(abs(X[,1] - threshold[1]) / max.window[1] < 1 & abs(X[,2] - threshold[2]) / max.window[2] < 1)
  X.inrange = X[inrange,]
  
  breaks1 = c(xx1 - bin.width[1]/2, max(xx1) + bin.width[1]/2)
  breaks2 = c(xx2 - bin.width[2]/2, max(xx2) + bin.width[2]/2)
  idx.1 = as.numeric(cut(X.inrange[,1], breaks = breaks1))
  idx.2 = as.numeric(cut(X.inrange[,2], breaks = breaks2))
  
  idx.to.bucket = sapply(1:nrow(X.inrange), function(iter) idx.1[iter] + (idx.2[iter] - 1) * num.bucket[1])
  X.counts = sapply(1:nrow(xx12), function(iter) sum(idx.to.bucket==iter))
  
  X.mids = t(sapply(1:nrow(xx12),function(iter) {
    ii = which(idx.to.bucket == iter)
    if(length(ii) == 0) return(as.numeric(xx12[iter,]))
    return(as.numeric(c(mean(X.inrange[ii,1]), mean(X.inrange[ii,2]))))
  }))
  
  realized.idx = which(X.counts > 0)
  num.realized = length(realized.idx)
  
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
  # nu_+ + nu_- = abs(t(Msq) %*% [row/col sums] %*% (n * gamma)); note that
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
  #  sum_k n[k] gamma[k] x[k,1] = 0
  #  sum_k n[k] gamma[k] x[k,2] = 0
  #  z >= 0
  #  
  #  and then, for a collection of fj...
  #  -sum_k n[k] gamma[k] fj[k] + fj >= 0
  
  treat = as.numeric(((xx12[realized.idx,1] - threshold[1]) > 0) & ((xx12[realized.idx,2] - threshold[2]) > 0))
  
  Dmat = diag(c(X.counts[realized.idx], lambda))
  dvec = rep(0, num.realized + 1)
  
  Amat = cbind(c(X.counts[realized.idx], 0),
               c(X.counts[realized.idx] * (2 * treat - 1), 0),
               c(X.counts[realized.idx] * X.mids[realized.idx,1], 0),
               c(X.counts[realized.idx] * X.mids[realized.idx,2], 0),
               c(X.counts[realized.idx] * xx12[realized.idx,1], 0),
               c(X.counts[realized.idx] * xx12[realized.idx,2], 0),
               c(rep(0, length(realized.idx)), 1))
  
  bvec = c(0, 2, 0, 0, 0, 0, 0)
  meq = 6
  
  if (change.derivative) {
    Amat = cbind(c(X.counts[realized.idx] * X.mids[realized.idx,1] * (2 * treat - 1), 0),
                 c(X.counts[realized.idx] * X.mids[realized.idx,2] * (2 * treat - 1), 0),
                 c(X.counts[realized.idx] * xx12[realized.idx,1] * (2 * treat - 1), 0),
                 c(X.counts[realized.idx] * xx12[realized.idx,2] * (2 * treat - 1), 0),
                 Amat)
    bvec = c(0, 0, 0, 0, bvec)
    meq = meq + 4
  }
  
  curr.value = 0
  
  for (iter in 1:100) {
    
    # Solve the problem while balancing the current dictionary of f-functions
    soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)
    
    # Find the function in the regularity class that yields the worst bias
    gamma.hat = rep(0, nrow(xx12))
    gamma.hat[realized.idx] = soln$solution[1:num.realized]
    objective.in = gamma.hat * X.counts
    
    worst.perturbation = lpSolve::lp(direction="max",
                            objective.in = objective.in,
                            const.mat = rbind(center.mat, nabla, -nabla),
                            const.dir = c(rep("=", nrow(center.mat)),
                                          rep("<=", 2*nrow(nabla))),
                            const.rhs =  c(rep(0, nrow(center.mat)),
                                           rep(1, 2*nrow(nabla))))
    ff = (worst.perturbation$solution - mean(worst.perturbation$solution))
    
    # Break once the worst bias over the {f}-class essentially matches
    # the actual worst-case bias over the regularity class
    claimed_worst_case_bias = soln$solution[num.realized + 1]
    actual_worst_case_bias = sum(X.counts[realized.idx] * ff[realized.idx] * gamma.hat[realized.idx])
    print(paste(claimed_worst_case_bias, " --- " , actual_worst_case_bias))
    if (claimed_worst_case_bias / actual_worst_case_bias > 0.99 & iter >= 3) break;
    
    image(xx1, xx2, matrix(ff, length(xx1), length(xx2)))
    
    # Add this function to the collection, so that the gamma_i balance it too
    Amat = cbind(Amat, c(-X.counts[realized.idx] * ff[realized.idx], 1))
    bvec = c(bvec, 0)
  }
  
  # Extract the weighting function gamma(x), as a function of x
  gamma.xx = rep(0, nrow(xx12))
  gamma.xx[realized.idx] = soln$solution[1:num.realized]
  
  gamma = rep(0, length(X))
  gamma[inrange] = gamma.hat[idx.to.bucket]
  
  imbalance = max.second.derivative^2 * sum(ff * gamma.hat * X.counts)^2
  
  ret = list(gamma=gamma, gamma.fun = data.frame(xx1=xx12[realized.idx,1], xx2=xx12[realized.idx,2], gamma=gamma.hat[realized.idx]), mean.var = sum(gamma^2), imbalance=imbalance)
  return(ret)
}