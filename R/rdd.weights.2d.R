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
      Y.hat = predict(lm(Y ~ X * W))
      sigma.sq = mean((Y - Y.hat)^2) * length(W) / (length(W) - 6)
    }
  }
  
  xx1 = seq(- max.window[1], max.window[1], length.out = num.bucket[1])
  xx2 = seq(- max.window[2], max.window[2], length.out = num.bucket[2])
  
  bin.width = c(xx1[2] - xx1[1], xx2[2] - xx2[1])
  xx12 = expand.grid(xx1, xx2)
  
  if (!change.derivative) {
    center.points = which(xx12[,1] %in% xx1[order(abs(xx1))[1:2]] &
                            xx12[,2] %in% xx2[order(abs(xx2))[1:2]])
  } else {
    center.points = which(xx12[,1] %in% xx1[order(abs(xx1))[1:4]] &
                            xx12[,2] %in% xx2[order(abs(xx2))[1:4]])
  }
  
  center.mat = matrix(0, length(center.points), nrow(xx12))
  for(iter in 1:length(center.points)) center.mat[iter, center.points[iter]] = 1
  
  # matrix that probes all local second derivatives, along axis and diagonals
  crit.1 = order(abs(xx1))[1:2]
  crit.2 = order(abs(xx2))[1:2]
  nabla = matrix(0, 4 * (num.bucket[1] - 2) * (num.bucket[2] - 2), nrow(xx12))
  curr.idx = 0
  for (i1 in 2:(num.bucket[1] - 1)) {
    for (i2 in 2:(num.bucket[2] - 1)) {
      # If "change.derivative", let f be different on different side of the boundary
      edge.1 = change.derivative & (i1 %in% crit.1 & i2 >= max(crit.2))
      edge.2 = change.derivative & (i2 %in% crit.2 & i1 >= max(crit.1))
      if (!edge.1) {
        nabla[curr.idx + 1,  (i1 - 1):(i1 + 1) + (i2 - 1) * num.bucket[1]] = c(1, -2, 1) / bin.width[1]^2
        curr.idx = curr.idx + 1
      }
      if (!edge.2) {
        nabla[curr.idx + 1, i1 + ((i2 - 2):i2) * num.bucket[1]] = c(1, -2, 1) / bin.width[2]^2
        curr.idx = curr.idx + 1
      }
      if (!(edge.1 | edge.2)) {
        nabla[curr.idx + 1,  (i1 - 1):(i1 + 1) + ((i2 - 2):i2) * num.bucket[1]] = c(1/2, -1, 1/2) / prod(bin.width)
        nabla[curr.idx + 2,  (i1 - 1):(i1 + 1) + (i2:(i2 - 2)) * num.bucket[1]] = c(1/2, -1, 1/2) / prod(bin.width)
        curr.idx = curr.idx + 2
      }
    }
  }
  
  # If change.derivate, throw away extra rows
  nabla = nabla[1:curr.idx,]
  
  # note the centering at the threshold
  inrange = which(abs(X[,1] - threshold[1]) / max.window[1] < 1 & abs(X[,2] - threshold[2]) / max.window[2] < 1)
  X.inrange = cbind(X[inrange, 1]  - threshold[1], X[inrange, 2]  - threshold[2])
  
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
  
  treat.all = as.numeric(((xx12[,1]) < 0) | ((xx12[,2]) < 0))
  treat = treat.all[realized.idx]
  
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
  worst = rep(0, nrow(xx12))
  
  lp.fail = FALSE
  
  for (iter in 1:100) {
    
    # Solve the problem while balancing the current dictionary of f-functions
    soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)
    
    # Find the function in the regularity class that yields the worst bias
    gamma.hat = rep(0, nrow(xx12))
    gamma.hat[realized.idx] = soln$solution[1:num.realized]
    objective.in = gamma.hat * X.counts
    
    nv = length(objective.in)
    
    if(!lp.fail) {
      SCL = max(rowSums(xx12^2))/2
      # Need to separate f into positive part and negative part, because
      # lpSolve only gives positive solutions
      worst.perturbation = lpSolve::lp(direction="max",
                                       c(objective.in, -objective.in),
                                       const.mat = rbind(
                                         cbind(center.mat, matrix(0, nrow(center.mat), nv)),
                                         cbind(matrix(0, nrow(center.mat), nv), center.mat),
                                         cbind(nabla, -nabla), cbind(-nabla, nabla), rep(1, 2 * nv)),
                                       const.dir = c(rep("=", 2 * nrow(center.mat)),
                                                     rep("<=", 2 * nrow(nabla) + 1)),
                                       const.rhs =  c(rep(0, 2 * nrow(center.mat)),
                                                      rep(1/SCL, 2 * nrow(nabla)), 2 * nv))
      worst = (worst.perturbation$solution[1:nv] - worst.perturbation$solution[nv + (1:nv)]) * SCL
      
      if (worst.perturbation$status != 0) {
        lp.fail = TRUE
        print("LP failed, used QP from now on...")
      }
    }
    
    if (lp.fail) {
      worst.pertubation = quadprog::solve.QP(diag(lambda / 100, nv),
                                             objective.in,
                                             t(rbind(center.mat, nabla, -nabla)),
                                             c(rep(0, nrow(center.mat)), rep(-1, 2 * nrow(nabla))),
                                             meq = 1)
      worst = worst.pertubation$solution
    }
    
    # Break once the worst bias over the {f}-class essentially matches
    # the actual worst-case bias over the regularity class
    claimed_worst_case_bias = soln$solution[num.realized + 1]
    actual_worst_case_bias = sum(X.counts[realized.idx] * worst[realized.idx] * gamma.hat[realized.idx])
    print(paste(claimed_worst_case_bias, " --- " , actual_worst_case_bias))
    if ((iter >= 3 & claimed_worst_case_bias / actual_worst_case_bias > 0.99) |
        (iter >= 10 & claimed_worst_case_bias / actual_worst_case_bias > 0.95)) break;
    
    image(xx1, xx2, matrix(worst, length(xx1), length(xx2)))
    abline(h=0)
    abline(v=0)
    
    # Add this function to the collection, so that the gamma_i balance it too
    Amat = cbind(Amat, c(-X.counts[realized.idx] * worst[realized.idx], 1))
    bvec = c(bvec, 0)
  }
  
  # Extract the weighting function gamma(x), as a function of x
  gamma.xx = rep(0, nrow(xx12))
  gamma.xx[realized.idx] = soln$solution[1:num.realized]
  
  gamma = rep(0, length(X))
  gamma[inrange] = gamma.hat[idx.to.bucket]
  
  max.bias = max.second.derivative * sum(worst * gamma.hat * X.counts)
  
  # If outcomes are provided, also compute confidence intervals for tau.
  if (!is.null(Y)) {
    
    # The point estimate
    tau.hat = sum(gamma * Y)
    
    # A heteroskedaticity-robust variance estimate
    regr.df = data.frame(X1=X.inrange[,1], X2=X.inrange[,2], W=treat.all[idx.to.bucket], Y=Y[inrange])
    Y.fit = lm(Y ~ W * (X1 + X2), data = regr.df)
    Y.resid.sq = rep(0, length(Y))
    Y.resid.sq[inrange] = (Y[inrange] - predict(Y.fit))^2 * length(inrange) / (length(inrange) - 6)
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
             gamma.fun = data.frame(xx1=xx12[realized.idx, 1] + threshold[1],
                                    xx2=xx12[realized.idx, 2] + threshold[2],
                                    gamma=gamma.xx[realized.idx]))
  class(ret) = "optrdd.2d"
  return(ret)
}