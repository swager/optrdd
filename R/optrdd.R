#' Optimized regression discontinuity design
#'
#' Optimized estimation and inference of treamtment effects identified
#' via regression discontinuities
#'
#' @param X The running variables.
#' @param Y The outcomes. If null, only optimal weights are computed.
#' @param W Treatment assignments, typically of the form 1(X >= c).
#' @param max.second.derivative A bound on the second derivative of mu_w(x) = E[Y(w) | X = x].
#' @param estimation.point Point "c" at which CATE is to be estimated. If estimation.point = NULL,
#'                         we estimate a weighted CATE, with weights chosen to minimize MSE,
#'                         as in Section 4.1 of Imbens and Wager (2017).
#' @param sigma.sq The irreducible noise level. If null, estimated from the data.
#' @param alpha Coverage probability of confidence intervals.
#' @param lambda.mult Optional multplier that can be used to over- or under-penalize variance.
#' @param bin.width Bin width for discrete approximation.
#' @param num.bucket Number of bins for discrete approximation. Can only be used if bin.width = NULL.
#' @param use.homoskedatic.variance Whether confidence intervals should be built assuming homoskedasticity.
#' @param use.spline Whether non-parametric components should be modeled as quadratic splines
#'                   in order to reduce the number of optimization parameters, and potentially
#'                   improving computational performance.
#' @param spline.df Number of degrees of freedom (per running variable) used for spline computation.
#' @param try.elnet.for.sigma.sq Whether an elastic net on a spline basis should be used for estimating sigma^2.
#' @param optimizer Which optimizer to use? Mosek is a commercial solver, but free
#'                  academic licenses are available. Needs to be installed separately.
#'                  ECOS is an open-source interior-point solver for conic problems,
#'                  made available via the CVXR wrapper.
#'                  Quadprog is the default R solver; it may be slow on large problems, but
#'                  is very accurate on small problems.
#'                  SCS is an open-source "operator splitting" solver that implements a first order
#'                  method for solving very large cone programs to modest accuracy. The speed of SCS may
#'                  be helpful for prototyping; however, the results may be noticeably less accurate.
#'                  SCS is also accessed via the CVXR wrapper.
#'                  The option "auto" uses a heuristic to choose.
#' @param verbose whether the optimizer should print progress information
#'
#' @return A trained optrdd object.
#' 
#' @references Domahidi, A., Chu, E., & Boyd, S. (2013, July).
#' ECOS: An SOCP solver for embedded systems.
#' In Control Conference (ECC), 2013 European (pp. 3071-3076). IEEE.
#' 
#' @references Imbens, G., & Wager, S. (2017).
#' Optimized Regression Discontinuity Designs.
#' arXiv preprint arXiv:1705.01677.
#' 
#' @references O’Donoghue, B., Chu, E., Parikh, N., & Boyd, S. (2016).
#' Conic optimization via operator splitting and homogeneous self-dual embedding.
#' Journal of Optimization Theory and Applications, 169(3), 1042-1068.
#' 
#' @examples
#' # Simple regression discontinuity with discrete X
#' n = 4000; threshold = 0
#' X = sample(seq(-4, 4, by = 8/41.5), n, replace = TRUE)
#' W = as.numeric(X >= threshold)
#' Y = 0.4 * W + 1 / (1 + exp(2 * X)) + 0.2 * rnorm(n)
#' # using 0.4 for max.second.derivative would have been enough
#' out.1 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5, estimation.point = threshold)
#' print(out.1); plot(out.1, xlim = c(-1.5, 1.5))
#' 
#' # Now, treatment is instead allocated in a neighborhood of 0
#' thresh.low = -1; thresh.high = 1
#' W = as.numeric(thresh.low <= X & X <= thresh.high)
#' Y = 0.2 * (1 + X) * W + 1 / (1 + exp(2 * X)) + rnorm(n)
#' # This estimates CATE at specifically chosen points
#' out.2 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5, estimation.point = thresh.low)
#' print(out.2); plot(out.2, xlim = c(-2.5, 2.5))
#' out.3 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5, estimation.point = thresh.high)
#' print(out.3); plot(out.3, xlim = c(-2.5, 2.5))
#' # This estimates a weighted CATE, with lower variance
#' out.4 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5)
#' print(out.4); plot(out.4, xlim = c(-2.5, 2.5))
#'
#' \dontrun{
#' # RDD with multivariate running variable.
#' X = matrix(runif(n*2, -1, 1), n, 2)
#' W = as.numeric(X[,1] < 0 | X[,2] < 0)
#' Y = X[,1]^2/3 + W * (1 + X[,2]) + rnorm(n)
#' out.5 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 1)
#' print(out.5); plot(out.5)
#' out.6 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 1, estimation.point = c(0, 0.5))
#' print(out.6); plot(out.6)}
#' 
#' @export
optrdd = function(X,
                  Y = NULL,
                  W,
                  max.second.derivative,
                  estimation.point = NULL,
                  sigma.sq = NULL,
                  alpha = 0.95,
                  lambda.mult = 1,
                  bin.width = NULL,
                  num.bucket = NULL,
                  use.homoskedatic.variance = FALSE,
                  use.spline = TRUE,
                  spline.df = NULL,
                  try.elnet.for.sigma.sq = FALSE,
                  optimizer = c("auto", "mosek", "ECOS", "quadprog", "SCS"),
                  verbose = TRUE) {

    n = length(W)
    if (class(W) == "logical") W = as.numeric(W)
    if (!all(W %in% c(0, 1))) stop("The treatment assignment W must be binary.")
    if (!is.null(Y) & (length(Y) != n)) { stop("Y and W must have same length.") }
    if (is.null(dim(X))) { X = matrix(X, ncol = 1) }
    if (nrow(X) != n) { stop("The number of rows of X and the length of W must match") }
    if (length(max.second.derivative) != 1) { stop("max.second.derivative must be of length 1.") }
    if (!is.null(bin.width) & !is.null(num.bucket)) { stop("Only one of bin.width or num.bucket may be used.") }
    
    nvar = ncol(X)
    if (nvar >= 3) { stop("Not yet implemented for 3 or more running variables.") }
    
    cate.at.pt = !is.null(estimation.point)
    if (is.null(estimation.point)) { estimation.point = colMeans(X) }
    univariate.monotone = (nvar == 1) &&
        ((max(X[W==0, 1]) <= min(X[W==1, 1])) || (max(X[W==1, 1]) <= min(X[W==0, 1])))
    
    mosek_available = requireNamespace("Rmosek", quietly = TRUE)

    # Naive initialization for sigma.sq if needed
    if (is.null(sigma.sq)) {
        if (is.null(Y)) {
            warning("Setting noise level to 1 as default...")
            sigma.sq = 1
        } else {
            Y.hat = stats::predict(stats::lm(Y ~ X * W))
            sigma.sq = mean((Y - Y.hat)^2) * length(W) / (length(W) - 2 - 2 * nvar)
            if (try.elnet.for.sigma.sq){
                if(ncol(X) > 1) {
                    stop("Elastic net for sigma squared not implemented with more than 1 running variable.")
                }
                linear.params = 1 + 2 * ncol(X)
                elnet.df = 7
                ridge.mat = cbind(W, X, W * X, matrix(0, length(W), 2 * elnet.df))
                ridge.mat[W==0, linear.params + 1:elnet.df] = splines::ns(X[W==0,], df = elnet.df)
                ridge.mat[W==1, linear.params + elnet.df + 1:elnet.df] = splines::ns(X[W==1,], df = elnet.df)
                elnet = glmnet::cv.glmnet(ridge.mat, Y,
                                          penalty.factor = c(rep(0, linear.params),
                                                             rep(1, 2 * elnet.df)),
                                          keep = TRUE, alpha = 0.5)
                elnet.hat = elnet$fit.preval[,!is.na(colSums(elnet$fit.preval)),drop=FALSE][, elnet$lambda == elnet$lambda.1se]
                sigma.sq.elnet = mean((elnet.hat - Y)^2)
                sigma.sq = min(sigma.sq, sigma.sq.elnet)
            }
        }
    }

    optimizer = match.arg(optimizer)
    if (optimizer == "auto") {
        if (nvar == 1 &&
            use.spline &&
            (univariate.monotone || !cate.at.pt) &&
            length(unique(c(X))) <= 100 &&
            max.second.derivative / sigma.sq <= 4) {
            optimizer = "quadprog"
        } else {
            if (mosek_available) {
                optimizer = "mosek"
            } else {
                optimizer = "ECOS"
            }
        }
    }

    if (optimizer == "mosek") {
        if (!mosek_available) {
            optimizer = "ECOS"
            warning("The mosek optimizer is not installed; using ECOS instead.")
        }
        # if (!requireNamespace("Rmosek", quietly = TRUE)) {
        #     optimizer = "quadprog"
        #     if (nvar >= 2) {
        #         op = options("warn")
        #         on.exit(options(op))
        #         options(warn=1)
        #         warning(paste("The mosek optimizer is not installed; using quadprog instead.",
        #                       "This may be very slow with more than one running variable."))
        #     } else {
        #         warning(paste("The mosek optimizer is not installed; using quadprog instead."))
        #     }
        # }
    }

    if (optimizer == "SCS") {
        warning(paste("SCS is a fast an free optimizer, but doesn't solve the problem exactly.",
                      "Resulting confidence intervals may be needlessly long.",
                      "It is recommended to also try MOSEK."))
    }

    # Create discrete grid on which to optimize, and assign each training example
    # to a cell.
    if (nvar == 1) {
        
        if (is.null(bin.width)) {
            if (is.null(num.bucket)) { num.bucket = 2000 }
            bin.width = (max(X[,1]) - min(X[,1])) / num.bucket
        }
        breaks = seq(min(X[,1]) - bin.width/2, max(X[,1]) + bin.width, by = bin.width)
        xx.grid = breaks[-1] - bin.width/2
        num.bucket = length(xx.grid)
        
        xx.grid = matrix(xx.grid, ncol = 1)
        xx.centered = matrix(xx.grid - estimation.point, ncol = 1)
        zeroed.idx = max(which(xx.centered < 0)) + c(0, 1)
        
        idx.to.bucket = as.numeric(cut(X[,1], breaks = breaks))
        
    } else if (nvar == 2) {
        
        if (is.null(bin.width)) {
            if (is.null(num.bucket)) { 
                if (optimizer == "quadprog") {
                    num.bucket = 900
                    warning(paste("Using coarse discrete approximation of size 30x30",
                                  "to make quadprog run faster (i.e., num.bucket = 900.)"))
                } else {
                    num.bucket = 10000
                }
            }
            bin.width = sqrt((max(X[,1]) - min(X[,1])) * (max(X[,2]) - min(X[,2])) / num.bucket)
        }
        breaks1 = seq(min(X[,1]) - bin.width/2, max(X[,1]) + bin.width, by = bin.width)
        breaks2 = seq(min(X[,2]) - bin.width/2, max(X[,2]) + bin.width, by = bin.width)
        xx1 = breaks1[-1] - bin.width/2
        xx2 = breaks2[-1] - bin.width/2
        xx.grid = expand.grid(xx1, xx2)
        xx.centered = t(t(xx.grid) - estimation.point)
        num.bucket = nrow(xx.grid)

        z1 = max(which(xx1 < estimation.point[1])) + c(0, 1)
        z2 = max(which(xx2 < estimation.point[2])) + c(0, 1)
        zeroed.idx = as.matrix(expand.grid(z1 - 1, z2 - 1)) %*% c(1, length(xx1))
        
        idx.1 = as.numeric(cut(X[,1], breaks = breaks1))
        idx.2 = as.numeric(cut(X[,2], breaks = breaks2))
        idx.to.bucket = sapply(1:nrow(X), function(iter) idx.1[iter] + (idx.2[iter] - 1) * length(xx1))
        
    } else {
        stop("Not yet implemented for 3 or more running variables.")
    }
    
    # Define matrix of constraints.
    # D2 is a raw curvature matrix, not accounting for bin.width.
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
    
    # Computational performance can be improved by representing non-parametric
    # components in a quadratic spline basis.
    if (use.spline) {
        if (is.null(spline.df)) {
            if (optimizer == "mosek" || optimizer == "SCS" || optimizer == "ECOS") {
                spline.df = c(100, 45 - 15 * cate.at.pt)[nvar]
            } else {
                spline.df = c(40, 10)[nvar]
            }
        }
        if ((nvar == 1 && spline.df > nrow(xx.centered) / 2) |
            (nvar == 2 && spline.df > min(length(xx1), length(xx2)) * 0.7)) {
            use.spline = FALSE
        }
    }
    
    if (use.spline) {
        if (nvar == 1) {
            basis.mat.raw = splines::bs(xx.grid, degree = 2, df=spline.df, intercept = TRUE)
            class(basis.mat.raw) = "matrix"
            basis.mat = Matrix::Matrix(basis.mat.raw)
        } else if (nvar == 2) {
            basis.mat.1 = splines::bs(xx.grid[,1], degree = 2, df=spline.df, intercept = TRUE)
            basis.mat.2 = splines::bs(xx.grid[,2], degree = 2, df=spline.df, intercept = TRUE)
            basis.mat = Matrix::sparse.model.matrix(~ basis.mat.1:basis.mat.2 + 0)
        } else {
            stop("Not yet implemented for 3 or more running variables.")
        }
        
        D2 = D2 %*% basis.mat
        selector.0 = selector.0 %*% basis.mat
        selector.1 = selector.1 %*% basis.mat
        centering.matrix = centering.matrix %*% basis.mat
        num.df = ncol(basis.mat)
    } else {
        num.df = num.bucket
    }
    
    # If we are in one dimension and have a threshold RDD, it is
    # enough to estimate a single nuisance function f0.
    change.slope = cate.at.pt
    two.fun = cate.at.pt & !univariate.monotone
    
    # We now prepare inputs to a numerical optimizer. We seek to
    # solve a discretized version of equation (18) from Imbens and Wager (2017).
    # The parameters to the problem are ordered as (G(0), G(1), lambda, f0, f1).
    num.lambda = 3 + nvar * (1 + change.slope)
    num.params = num.realized.0 + num.realized.1 + num.lambda + (1 + two.fun) * num.df
    
    # The quadratic component of the objective is 1/2 sum_j Dmat.diagonal_j * params_j^2
    Dmat.diagonal = c((X.counts[realized.idx.0] - W.counts[realized.idx.0]) / 2 / sigma.sq,
                      (W.counts[realized.idx.1]) / 2 / sigma.sq,
                      lambda.mult / max.second.derivative^2 / 2,
                      rep(0, num.lambda - 1 + (1 + two.fun) * num.df))
    # The linear component of the objective is sum(dvec * params)
    dvec = c(rep(0, num.realized.0 + num.realized.1 + 1), 1, -1,
             rep(0, num.lambda - 3 + (1 + two.fun) * num.df))
    # We now impose constrains on Amat %*% params.
    # The first meq constrains are equality constraints (Amat %*% params = 0);
    # the remaining ones are inequalities (Amat %*% params >= 0).
    Amat = rbind(
        # Defines G(0) in terms of the other problem parameters (equality constraint)
        cbind(Matrix::Diagonal(num.realized.0, -1),
              Matrix::Matrix(0, num.realized.0, num.realized.1),
              0, 0, 1, xx.centered[realized.idx.0,],
              if(change.slope) { -xx.centered[realized.idx.0,] } else { numeric() },
              selector.0,
              if(two.fun) { Matrix::Matrix(0, num.realized.0, num.df) } else { numeric() }),
        # Defines G(1) in terms of the other problem parameters (equality constraint)
        cbind(Matrix::Matrix(0, num.realized.1, num.realized.0),
              Matrix::Diagonal(num.realized.1, -1),
              0, 1, 0, xx.centered[realized.idx.1,],
              if(change.slope) { xx.centered[realized.idx.1,] } else { numeric() },
              if(two.fun) {
                  cbind(Matrix::Matrix(0, num.realized.1, num.df), selector.1)
              } else {
                  selector.1
              }),
        # Ensure that f_w(c), f'_w(c) = 0
        cbind(Matrix::Matrix(0, length(zeroed.idx), num.realized.0 + num.realized.1 + num.lambda),
              centering.matrix,
              Matrix::Matrix(0, length(zeroed.idx), as.numeric(two.fun) * num.df)),
        if(two.fun) {
            cbind(Matrix::Matrix(0, length(zeroed.idx), num.realized.0 + num.realized.1 + num.lambda + num.df),
                  centering.matrix)
        } else { numeric() },
        # Bound the second derivative of f_0 by lambda_1
        cbind(Matrix::Matrix(0, 2 * nrow(D2), num.realized.0 + num.realized.1),
              bin.width^2,
              matrix(0, 2 * nrow(D2), num.lambda - 1),
              rbind(D2, -D2),
              if (two.fun) { Matrix::Matrix(0, 2 * nrow(D2), num.df) } else { numeric() }),
        # Bound the second derivative of f_1 by lambda_1
        if (two.fun) {
            cbind(Matrix::Matrix(0, 2 * nrow(D2), num.realized.0 + num.realized.1),
                  bin.width^2,
                  matrix(0, 2 * nrow(D2), num.lambda - 1 + num.df),
                  rbind(D2, -D2))
        })
    
    meq = num.realized.0 + num.realized.1 + length(zeroed.idx) * (1 + two.fun)
    bvec = rep(0, nrow(Amat))
    
    gamma.0 = rep(0, num.bucket)
    gamma.1 = rep(0, num.bucket)
    
    if (optimizer == "quadprog") {
        
        if (verbose) {
            print(paste0("Running quadrprog with problem of size: ",
                         dim(Amat)[1], " x ", dim(Amat)[2], "..."))
        }
        # For quadprog, we need Dmat to be positive definite, which is why we add a small number to the diagonal.
        # The conic implementation via mosek does not have this issue.
        soln = quadprog::solve.QP(Matrix::Diagonal(num.params, Dmat.diagonal + 0.000000001),
                                  -dvec,
                                  Matrix::t(Amat),
                                  bvec,
                                  meq=meq)
        
        gamma.0[realized.idx.0] = - soln$solution[1:num.realized.0] / sigma.sq / 2
        gamma.1[realized.idx.1] = - soln$solution[num.realized.0 + 1:num.realized.1] / sigma.sq / 2
        t.hat = soln$solution[num.realized.0 + num.realized.1 + 1] / (2 * max.second.derivative^2)
        
    } else if (optimizer == "SCS" || optimizer == "ECOS") {
      
        if (verbose && optimizer == "SCS") {
            print(paste0("Running CVXR/SCS with problem of size: ",
                         dim(Amat)[1], " x ", dim(Amat)[2], "..."))
        }
        
        if (verbose && optimizer == "ECOS") {
            print(paste0("Running CVXR/ECOS with problem of size: ",
                         dim(Amat)[1], " x ", dim(Amat)[2], "..."))
        }
        
        xx = CVXR::Variable(ncol(Amat))
        objective = sum(Dmat.diagonal/2 * xx^2 + dvec * xx)
        contraints = list(
            Amat[1:meq,] %*% xx == bvec[1:meq],
            Amat[(meq+1):nrow(Amat),] %*% xx >= bvec[(meq+1):nrow(Amat)]
        )
        cvx.problem = CVXR::Problem(CVXR::Minimize(objective), contraints)
        cvx.output = solve(cvx.problem, solver = optimizer, verbose = verbose)
        
        if (cvx.output$status != "optimal") {
            warning(paste0("CVXR returned with status: ",
                           cvx.output$status,
                           ". For better results, try another optimizer (MOSEK is recommended)."))
        }
        
        result = cvx.output$getValue(xx)
        gamma.0[realized.idx.0] = - result[1:num.realized.0] / sigma.sq / 2
        gamma.1[realized.idx.1] = - result[num.realized.0 + 1:num.realized.1] / sigma.sq / 2
        t.hat = result[num.realized.0 + num.realized.1 + 1] / (2 * max.second.derivative^2)
      
    } else if (optimizer == "mosek") {
        
        # We need to rescale our optimization parameters, such that Dmat has only
        # ones and zeros on the diagonal; i.e.,
        A.natural = Amat %*% Matrix::Diagonal(ncol(Amat), x=1/sqrt(Dmat.diagonal + as.numeric(Dmat.diagonal == 0)))
        
        mosek.problem <- list()
        
        # The A matrix relates parameters to constraints, via
        # blc <= A * { params } <= buc, and blx <= params <= bux
        # The conic fomulation adds two additional parameters to the problem, namely
        # a parameter "S" and "ONE", that are characterized by a second-order cone constraint
        # S * ONE >= 1/2 {params}' Dmat {params},
        # and the equality constraint ONE = 1
        mosek.problem$A <- cbind(A.natural, Matrix::Matrix(0, nrow(A.natural), 2))
        mosek.problem$bc <- rbind(blc = rep(0, nrow(A.natural)), buc = c(rep(0, meq), rep(Inf, nrow(A.natural) - meq)))
        mosek.problem$bx <- rbind(blx = c(rep(-Inf, ncol(A.natural)), 0, 1), bux = c(rep(Inf, ncol(A.natural)), Inf, 1))
        
        # This is the cone constraint
        mosek.problem$cones <- cbind(list("RQUAD", c(ncol(Amat) + 1, ncol(Amat) + 2, which(Dmat.diagonal != 0))))
        
        # We seek to minimize c * {params}
        mosek.problem$sense <- "min"
        mosek.problem$c <- c(dvec, 1, 0)
        
        #mosek.problem$dparam= list(BASIS_TOL_S=1.0e-9, BASIS_TOL_X=1.0e-9)
        
        if (verbose) {
            mosek.out = Rmosek::mosek(mosek.problem)
        } else {
            mosek.out = Rmosek::mosek(mosek.problem, opts=list(verbose=0))
        }
        
        if (mosek.out$response$code != 0) {
            warning(paste("MOSEK returned with status",
                          mosek.out$response$msg,
                          "For better results, try another optimizer."))
        }
        
        # We now also need to re-adjust for "natural" scaling
        gamma.0[realized.idx.0] = - mosek.out$sol$itr$xx[1:num.realized.0] / sigma.sq / 2 /
            sqrt(Dmat.diagonal[1:num.realized.0])
        gamma.1[realized.idx.1] = - mosek.out$sol$itr$xx[num.realized.0 + 1:num.realized.1] / sigma.sq / 2 /
            sqrt(Dmat.diagonal[num.realized.0 + 1:num.realized.1])
        t.hat = mosek.out$sol$itr$xx[num.realized.0 + num.realized.1 + 1] / (2 * max.second.derivative^2) /
            sqrt(Dmat.diagonal[num.realized.0 + num.realized.1 + 1])
        
    } else {
        stop("Optimizer choice not valid.")
    }
    
    # Now map the x-wise functions into a weight for each observation
    gamma = rep(0, length(W))
    gamma[W==0] = as.numeric(Matrix::t(bucket.map[,which(W==0)]) %*% gamma.0)
    gamma[W==1] = as.numeric(Matrix::t(bucket.map[,which(W==1)]) %*% gamma.1)
    
    # Patch up numerical inaccuracies
    gamma[W==0] = -gamma[W==0] / sum(gamma[W==0])
    gamma[W==1] = gamma[W==1] / sum(gamma[W==1])
    
    # Compute the worst-case imbalance...
    max.bias = max.second.derivative * t.hat
    
    # If outcomes are provided, also compute confidence intervals for tau.
    if (!is.null(Y)) {
        
        # The point estimate
        tau.hat = sum(gamma * Y)
        
        if (use.homoskedatic.variance) {
            se.hat.tau = sqrt(sum(gamma^2 * sigma.sq))
        } else {
            # A heteroskedaticity-robust variance estimate.
            # Weight the regression towards areas used in estimation, so we are not
            # too vulnerable to curvature effects far from the boundary.
            regr.df = data.frame(X=X, W=W, Y=Y, gamma.sq = gamma^2)
            regr.df = regr.df[regr.df$gamma.sq != 0,]
            Y.fit = stats::lm(Y ~ X * W, data = regr.df, weights = regr.df$gamma.sq)
            self.influence = stats::lm.influence(Y.fit)$hat
            Y.resid.sq = (regr.df$Y - stats::predict(Y.fit))^2
            se.hat.tau = sqrt(sum(Y.resid.sq * regr.df$gamma.sq / (1 - self.influence)))

            if (!try.elnet.for.sigma.sq & se.hat.tau^2 < 0.8 * sum(regr.df$gamma.sq) * sigma.sq) {
                warning(paste("Implicit estimate of sigma^2 may be too pessimistic,",
                              "resulting in valid but needlessly long confidence intervals.",
                              "Try the option `try.elnet.for.sigma.sq = TRUE` for potentially improved performance."))
            }
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
               gamma.fun.0 = data.frame(xx=xx.grid[realized.idx.0,],
                                        gamma=gamma.0[realized.idx.0]),
               gamma.fun.1 = data.frame(xx=xx.grid[realized.idx.1,],
                                        gamma=gamma.1[realized.idx.1]))
    class(ret) = "optrdd"
    return(ret)
}
