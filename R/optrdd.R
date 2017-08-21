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
                      use.homoskedatic.variance = FALSE,
                      use.spline = TRUE,
                      spline.df = NULL) {
    
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
        num.bucket = length(xx.grid)
        
        xx.centered = matrix(xx.grid - center, ncol = 1)
        zeroed.idx = max(which(xx.centered < 0)) + c(0, 1)
        
        idx.to.bucket = as.numeric(cut(X[,1], breaks = breaks))
        
    } else if (nvar == 2) {
        
        if (is.null(bin.width)) {
            bin.width = sqrt((max(X[,1]) - min(X[,1])) * (max(X[,2]) - min(X[,2])) / 400)
        }
        breaks1 = seq(min(X[,1]) - bin.width/2, max(X[,1]) + bin.width, by = bin.width)
        breaks2 = seq(min(X[,2]) - bin.width/2, max(X[,2]) + bin.width, by = bin.width)
        xx1 = breaks1[-1] - bin.width/2
        xx2 = breaks2[-1] - bin.width/2
        xx.grid = expand.grid(xx1, xx2)
        num.bucket = nrow(xx.grid)
        
        xx.centered = t(t(xx.grid) - center)
        z1 = max(which(xx1 < 0)) + c(0, 1)
        z2 = max(which(xx2 < 0)) + c(0, 1)
        zeroed.idx = as.matrix(expand.grid(z1 - 1, z2 - 1)) %*% c(1, length(xx1))
        
        idx.1 = as.numeric(cut(X[,1], breaks = breaks1))
        idx.2 = as.numeric(cut(X[,2], breaks = breaks2))
        idx.to.bucket = sapply(1:nrow(X), function(iter) idx.1[iter] + (idx.2[iter] - 1) * length(xx1))
        
    } else {
        stop("Not yet implemented for 3 or more running variables.")
    }
    
    # Define matrix of constraints
    # D2 is a raw curvature matrix, not accounting for bin.width
    if (nvar == 1) {
        D2 = Matrix::bandSparse(n=num.bucket-2, m=num.bucket, k = c(0, 1, 2),
                                diag = list(rep(1, num.bucket), rep(-2, num.bucket))[c(1, 2, 1)])
    } else if (nvar == 2) {
        all.idx = expand.grid(1:length(xx1), 1:length(xx2))
        # remove corners
        all.idx = all.idx[!(all.idx[,1] %in% c(1, length(xx1)) & all.idx[,2] %in% c(1, length(xx2))),]
        D2.entries = do.call(rbind, sapply(1:nrow(all.idx), function(i12) {
            i1 = all.idx[i12,1]
            i2 = all.idx[i12,2]
            edge.1 = i1 %in% c(1, length(xx1))
            edge.2 = i2 %in% c(1, length(xx2))
            rbind(
                if (!edge.1) {
                    cbind(j=(i1 - 1):(i1 + 1) + (i2 - 1) * length(xx1),
                          x=c(1, -2, 1))
                } else { numeric() },
                if (!edge.2) {
                    cbind(j=i1 + ((i2 - 2):i2) * length(xx1),
                          x=c(1, -2, 1))
                } else { numeric() },
                if (!(edge.1 | edge.2)) {
                    cbind(j = c((i1 - 1):(i1 + 1) + ((i2 - 2):i2) * length(xx1),
                                (i1 - 1):(i1 + 1) + (i2:(i2 - 2)) * length(xx1)),
                          x = c(c(1/2, -1, 1/2), c(1/2, -1, 1/2)))
                } else { numeric() })
        }))
        D2.i = c(t(matrix(rep(1:(nrow(D2.entries)/3), 3), ncol = 3)))
        D2 = Matrix::sparseMatrix(i=D2.i, j=D2.entries[,1], x=D2.entries[,2])
    } else {
        stop("Not yet implemented for 3 or more running variables.")
    }
    
    # Construct a (weighted) histogram representing the X and Y values.
    fact = factor(idx.to.bucket, levels = as.character(1:num.bucket))
    bucket.map = Matrix::sparse.model.matrix(~fact + 0, transpose = TRUE)
    X.counts = as.numeric(bucket.map %*% rep(1, n))
    W.counts = as.numeric(bucket.map %*% W)
    
    realized.idx.0 = which(X.counts > W.counts)
    realized.idx.1 = which(W.counts > 0)
    num.realized.0 = length(realized.idx.0)
    num.realized.1 = length(realized.idx.1)
    
    # These matrices map non-parametric dual parameters f_w(x) to buckets
    # where gamma needs to be evaluated
    selector.0 = Matrix::Diagonal(num.bucket, 1)[realized.idx.0,]
    selector.1 = Matrix::Diagonal(num.bucket, 1)[realized.idx.1,]
    
    # We force centering.matrix %*% f_w(x) = 0, to ensure identification of f_w(x)
    centering.matrix = Matrix::sparseMatrix(dims = c(length(zeroed.idx), num.bucket),
                                            i = 1:length(zeroed.idx),
                                            j = zeroed.idx,
                                            x = rep(1, length(zeroed.idx)))
    
    if (use.spline) {
      if (is.null(spline.df)) {
        spline.df = 40 / nvar^2
      }
      if ((nvar == 1 && spline.df > nrow(xx.centered) / 2) |
          (nvar == 2 && spline.df > min(length(xx1), length(xx2)) * 0.7)) {
        use.spline = FALSE
      }
    }
    
    print(use.spline)
    
    if (use.spline) {
      
      if (nvar == 1) {
        basis.mat.raw = splines::ns(xx.grid, df=spline.df, intercept = TRUE)
        class(basis.mat.raw) = "matrix"
        basis.mat = Matrix::Matrix(basis.mat.raw)
      } else if (nvar == 2) {
        basis.mat.1 = splines::ns(xx.grid[,1], df=spline.df, intercept = TRUE)
        basis.mat.2 = splines::ns(xx.grid[,2], df=spline.df, intercept = TRUE)
        basis.mat = Matrix::sparse.model.matrix(~ basis.mat.1:basis.mat.2 + 0)
      } else {
        stop("Not yet implemented for 3 or more running variables.")
      }
      
      D2 = D2 %*% basis.mat
      selector.0 = selector.0 %*% basis.mat
      selector.1 = selector.1 %*% basis.mat
      centering.matrix = centering.matrix %*% basis.mat
      num.df = spline.df
    } else {
      basis.mat = Matrix::Diagonal(num.bucket, 1)
      num.df = num.bucket
    }
    
    # solution to opt problem is (G(0), G(1), lambda, f0, f1)
    num.lambda = 3 + nvar *(1 + cate.at.pt)
    num.params = num.realized.0 + num.realized.1 + num.lambda + (1 + cate.at.pt) * num.df
    Dmat.diagonal = c((X.counts[realized.idx.0] - W.counts[realized.idx.0]) / 2 / sigma.sq,
                      (W.counts[realized.idx.1]) / 2 / sigma.sq,
                      lambda.mult / max.second.derivative^2 / 2,
                      rep(0, num.lambda - 1 + (1 + cate.at.pt) * num.df))
    Dmat = Matrix::Diagonal(num.params, Dmat.diagonal + 0.000000000001)
    dvec = c(rep(0, num.realized.0 + num.realized.1 + 1), -1, 1,
             rep(0, num.lambda - 3 + (1 + cate.at.pt) * num.df))
    Amat = Matrix::t(rbind(
        cbind(Matrix::Diagonal(num.realized.0, -1),
              Matrix::Matrix(0, num.realized.0, num.realized.1),
              0, 0, 1, xx.centered[realized.idx.0,],
              if(cate.at.pt) { -xx.centered[realized.idx.0,] } else { numeric() },
              selector.0,
              if(cate.at.pt) { Matrix::Matrix(0, num.realized.0, num.df) } else { numeric() }),
        cbind(Matrix::Matrix(0, num.realized.1, num.realized.0),
              Matrix::Diagonal(num.realized.1, -1),
              0, 1, 0, xx.centered[realized.idx.1,],
              if(cate.at.pt) { xx.centered[realized.idx.1,] } else { numeric() },
              if(cate.at.pt) {
                  cbind(Matrix::Matrix(0, num.realized.1, num.df), selector.1)
              } else {
                  selector.1
              }),
        cbind(Matrix::Matrix(0, length(zeroed.idx), num.realized.0 + num.realized.1 + num.lambda),
              centering.matrix,
              Matrix::Matrix(0, length(zeroed.idx), as.numeric(cate.at.pt) * num.df)),
        if(cate.at.pt) {
          cbind(Matrix::Matrix(0, length(zeroed.idx), num.realized.0 + num.realized.1 + num.lambda + num.df),
                centering.matrix)
        } else { numeric() },
        c(rep(0, num.realized.0 + num.realized.1), 1, rep(0, num.lambda - 1 + (1 + cate.at.pt) * num.df)),
        cbind(Matrix::Matrix(0, 2 * nrow(D2), num.realized.0 + num.realized.1),
              bin.width^2,
              matrix(0, 2 * nrow(D2), num.lambda - 1),
              rbind(D2, -D2),
              if (cate.at.pt) { Matrix::Matrix(0, 2 * nrow(D2), num.df) } else { numeric() }),
        if (cate.at.pt) {
            cbind(Matrix::Matrix(0, 2 * nrow(D2), num.realized.0 + num.realized.1),
                  bin.width^2,
                  matrix(0, 2 * nrow(D2), num.lambda - 1 + num.df),
                  rbind(D2, -D2))
        }))
    
    meq = num.realized.0 + num.realized.1 + length(zeroed.idx) * (1 + cate.at.pt)
    bvec = rep(0, ncol(Amat))
    
    soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=meq, factorized=FALSE)
    
    gamma.0 = rep(0, num.bucket)
    gamma.0[realized.idx.0] = - soln$solution[1:num.realized.0] / sigma.sq / 2
    gamma.1 = rep(0, num.bucket)
    gamma.1[realized.idx.1] = - soln$solution[num.realized.0 + 1:num.realized.1] / sigma.sq / 2
    
    # Now map this x-wise function into a weight for each observation
    gamma = rep(0, length(W))
    gamma[W==0] = as.numeric(Matrix::t(bucket.map[,which(W==0)]) %*% gamma.0)
    gamma[W==1] = as.numeric(Matrix::t(bucket.map[,which(W==1)]) %*% gamma.1)
    
    # Compute the worst-case imbalance...
    t.hat = soln$solution[num.realized.0 + num.realized.1 + 1] / (2 * max.second.derivative^2)
    max.bias = max.second.derivative * t.hat
    
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

