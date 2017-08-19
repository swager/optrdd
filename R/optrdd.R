#' Optimized regression discontinuity design
#'
#' Optimized estimation and inference of treamtment effects identified
#' via regression discontinuities
#'
#' @param X The running variables.
#' @param Y The outcomes. If null, only optimal weights are computed.
#' @param W Treatment assignments, typically of the form 1(X >= c).
#' @param max.second.derivative A bound on the second derivative of mu_w(x) = E[Y(w) | X = x].
#' @param center Point "c" at which CATE is to be estimated. If center = NULL, we estimate
#'               a weighted CATE, with weights chosen to minimize MSE (see Section 4.1 of paper).
#' @param sigma.sq The irreducible noise level. If null, estimated from the data.
#' @param alpha Coverage probability of confidence intervals.
#' @param lambda.mult Optional multplier that can be used to over- or under-penalize variance.
#' @param bin.width Bin width for discrete approximation.
#' @param use.homoskedatic.variance Whether confidence intervals should be built assuming homoskedasticity.
#'
#' @return A trained optrdd object.
#' @export
optrdd.new = function(X,
                  Y = NULL,
                  W,
                  max.second.derivative,
                  center = NULL,
                  sigma.sq = NULL,
                  alpha = 0.95,
                  lambda.mult = 1,
                  bin.width = NULL,
                  use.homoskedatic.variance = FALSE) {
  
    n = length(W)
    if (class(W) == "logical") W = as.numeric(W)
    if (!is.null(Y) & (length(Y) != n)) { stop("Y and W must have same length.") }
    if (is.null(dim(X))) { X = matrix(X, ncol = 1) }
    if (nrow(X) != n) { stop("The number of rows of X and the length of W must match") }
    
    nvar = ncol(X)
    if (nvar >= 3) { stop("Not yet implemented for 3 or more running variables.") }
    
    cate.at.pt = !is.null(center)
    if (is.null(center)) { center = colMeans(X) }
    
    # Naive initialization for sigma.sq if needed
    if (is.null(sigma.sq)) {
        if (is.null(Y)) {
            warning("Setting noise level to 1 as default...")
            sigma.sq = 1
        } else {
            Y.hat = predict(lm(Y ~ X * W))
            sigma.sq = mean((Y - Y.hat)^2) * length(W) / (length(W) - 2 - 2 * nvar)
        }
    }
    
    # Create discrete grid on which to optimize, and assign each training example
    # to a cell.
    
    if (nvar == 1) {
        
        if (is.null(bin.width)) {
            bin.width = (max(X[,1]) - min(X[,1])) / 400
        }
        breaks = seq(min(X[,1]) - bin.width/2, max(X[,1]) + bin.width, by = bin.width)
        xx.grid = breaks[-1] - bin.width/2
        idx.to.bucket = as.numeric(cut(X[,1], breaks = breaks))
        
    } else if (nvar == 2) {
        
        if (is.null(bin.width)) {
            bin.width = sqrt((max(X[,1]) - min(X[,1])) * (max(X[,2]) - min(X[,2])) / 1600)
        }
        breaks1 = seq(min(X[,1]) - bin.width/2, max(X[,1]) + bin.width, by = bin.width)
        breaks2 = seq(min(X[,2]) - bin.width/2, max(X[,2]) + bin.width, by = bin.width)
        xx1 = breaks1[-1] - bin.width/2
        xx2 = breaks2[-1] - bin.width/2
        xx.grid = expand.grid(xx1, xx2)
        
        idx.1 = as.numeric(cut(X[,1], breaks = breaks1))
        idx.2 = as.numeric(cut(X[,2], breaks = breaks2))
        idx.to.bucket = sapply(1:nrow(X), function(iter) idx.1[iter] + (idx.2[iter] - 1) * length(xx1))
        
    } else {
        stop("Not yet implemented for 3 or more running variables.")
    }
    
    # Define matrix of constraints
    # D2 is a raw curvature matrix, not accounting for bin.width
    if (nvar == 1) {
         D2 = Matrix::bandSparse(n=length(xx.grid)-2, m=length(xx.grid), k = c(0, 1, 2),
                    diag = list(rep(1, length(xx.grid)), rep(-2, length(xx.grid)))[c(1, 2, 1)])
         min.idx = max(which(xx.grid < center))
         center.idx = c(min.idx, min.idx + 1)
    } else {
        stop("Not yet implemented for 3 or more running variables.")
    }
    
    
    # Construct a (weighted) histogram representing the X and Y values.
    fact = factor(idx.to.bucket, levels = as.character(1:length(xx.grid)))
    bucket.map = Matrix::sparse.model.matrix(~fact + 0, transpose = TRUE)
    X.counts = as.numeric(bucket.map %*% rep(1, n))
    W.counts = as.numeric(bucket.map %*% W)
    
    realized.idx.0 = which(X.counts > W.counts)
    realized.idx.1 = which(W.counts > 0)
    num.realized.0 = length(realized.idx.0)
    num.realized.1 = length(realized.idx.1)
    
    # solution to opt problem is (G(0), G(1), lambda, f0, f1)
    num.lambda = 4 + cate.at.pt
    num.params = num.realized.0 + num.realized.1 + num.lambda + (1 + cate.at.pt) * ncol(D2)
    Dmat = Matrix::Diagonal(num.params,
                            c((X.counts[realized.idx.0] - W.counts[realized.idx.0]) / 2 / sigma.sq,
                              (W.counts[realized.idx.1]) / 2 / sigma.sq,
                              max.second.derivative^2 / 2,
                              rep(0, num.lambda - 1 + (1 + cate.at.pt) * ncol(D2)))
                            + 0.0000000001)
    dvec = c(rep(0, num.realized.0 + num.realized.1 + 1), -1, 1,
             rep(0, num.lambda - 3 + (1 + cate.at.pt) * ncol(D2)))
    Amat = Matrix::t(rbind(
             cbind(Matrix::Diagonal(num.realized.0, -1),
                   Matrix::Matrix(0, num.realized.0, num.realized.1),
                   0, 0, 1, xx.grid[realized.idx.0] - center,
                   if(cate.at.pt) { -xx.grid[realized.idx.0] + center } else { numeric() },
                   Matrix::Diagonal(length(xx.grid), 1)[realized.idx.0,],
                   if(cate.at.pt) { Matrix::Matrix(0, num.realized.0, length(xx.grid)) } else { numeric() }),
             cbind(Matrix::Matrix(0, num.realized.1, num.realized.0),
                   Matrix::Diagonal(num.realized.1, -1),
                   0, 1, 0, xx.grid[realized.idx.1] - center,
                   if(cate.at.pt) { xx.grid[realized.idx.1] - center } else { numeric() },
                   if(cate.at.pt) {
                       cbind(matrix(0, num.realized.1, length(xx.grid)),
                             Matrix::Diagonal(length(xx.grid), 1)[realized.idx.1,])
                   } else {
                       Matrix::Diagonal(length(xx.grid), 1)[realized.idx.1,]
                   }),
             c(rep(0, num.realized.0 + num.realized.1), 1, rep(0, num.lambda - 1 + (1 + cate.at.pt) * ncol(D2))),
             cbind(matrix(0, 2 * nrow(D2), num.realized.0 + num.realized.1),
                   bin.width^2 * max.second.derivative,
                   matrix(0, 2 * nrow(D2), num.lambda - 1),
                   rbind(D2, -D2),
                   if (cate.at.pt) { matrix(0, 2 * nrow(D2), length(xx.grid)) } else { numeric() }),
             if (cate.at.pt) {
                 cbind(matrix(0, 2 * nrow(D2), num.realized.0 + num.realized.1),
                       bin.width^2 * max.second.derivative,
                       matrix(0, 2 * nrow(D2), num.lambda - 1 + length(xx.grid)),
                       rbind(D2, -D2))
             }))
    
    meq = num.realized.0 + num.realized.1
    bvec = rep(0, ncol(Amat))
    
    soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=meq, factorized=FALSE)
    
    gamma.0 = rep(0, length(xx.grid))
    gamma.0[realized.idx.0] = - soln$solution[1:num.realized.0] / sigma.sq / 2
    gamma.1 = rep(0, length(xx.grid))
    gamma.1[realized.idx.1] = - soln$solution[num.realized.0 + 1:num.realized.1] / sigma.sq / 2
  
  # Now map this x-wise function into a weight for each observation
  gamma = rep(0, length(X))
  gamma[W==0] = as.numeric(Matrix::t(bucket.map[,W==0]) %*% gamma.0)
  gamma[W==1] = as.numeric(Matrix::t(bucket.map[,W==1]) %*% gamma.1)
  
  # Compute the worst-case imbalance...
  max.bias = soln$solution[num.realized.0 + num.realized.1 + 1] / (2 * max.second.derivative^2)
  
  # If outcomes are provided, also compute confidence intervals for tau.
  if (!is.null(Y)) {
    
    # The point estimate
    tau.hat = sum(gamma * Y)
    
    if (use.homoskedatic.variance) {
      se.hat.tau = sqrt(sum(gamma^2 * sigma.sq / num.samples))
    } else {
      # A heteroskedaticity-robust variance estimate
      regr.df = data.frame(X=X, W=W, Y=Y)
      Y.fit = lm(Y ~ X * W, data = regr.df)
      Y.resid.sq = (Y - predict(Y.fit))^2 * length(W) / (length(W) - 4)
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
             gamma.fun.0 = data.frame(xx=xx.grid[realized.idx.0],
                                     gamma=gamma.0[realized.idx.0]),
             gamma.fun.1 = data.frame(xx=xx.grid[realized.idx.1],
                                      gamma=gamma.1[realized.idx.1]))
  class(ret) = "optrdd"
  return(ret)
}

