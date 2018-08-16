set.seed(1)

max.second.derivative = 1
K = 20
vv = 1.27 * 1:20
supp = 2 * (vv - trunc(vv)) - 1
prob = rexp(K)
prob = prob/sum(prob)

n = 2000

bucket = as.numeric(1:K %*% rmultinom(n, 1, prob))

X = supp[bucket]
threshold = 0
W = as.numeric(X >= threshold)
Y = 10 + X + rnorm(n) + W

# Test methods initially, and confirm gamma moments
rdd.free = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, verbose = FALSE, optimizer = "mosek")
rdd.cate = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose = FALSE, optimizer = "mosek")

test_that("optrdd gammas satisfy constraints", {
    tol = rdd.cate$gamma.fun.0[2, 1] - rdd.cate$gamma.fun.0[1, 1]
    expect_equal(sum(rdd.cate$gamma), 0)
    expect_equal(sum(rdd.cate$gamma * W), 1)
    expect_equal(sum(rdd.cate$gamma * X), 0, tolerance = tol)
    expect_equal(sum(rdd.cate$gamma * X * W), 0, tolerance = tol)
    expect_equal(sum(rdd.free$gamma), 0)
    expect_equal(sum(rdd.free$gamma * W), 1)
    expect_equal(sum(rdd.free$gamma * X), 0, tolerance = tol)
})

test_that("cate constraint hurts", {
    expect_true(rdd.cate$tau.plusminus > rdd.free$tau.plusminus)
})

# Check implementation against legacy implementation
test_that("results match legacy implementation", {
    skip_on_cran()
    source("../../baselines/old.optrdd.R")
    rdd.old = optrdd.primal(X=X, Y=Y, threshold = 0, max.second.derivative = max.second.derivative)
    expect_equal(rdd.cate$tau.hat, rdd.old$tau.hat, tolerance = rdd.cate$tau.plusminus)
    expect_equal(rdd.cate$tau.plusminus, rdd.old$tau.plusminus, tolerance = 0.01)
})

# Test optimization strategies
rdd.free.raw = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, bin.width = 0.05, use.spline = FALSE, verbose = FALSE)
rdd.cate.raw = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, bin.width = 0.05, use.spline = FALSE, verbose = FALSE)
rdd.free.qp = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, optimizer = "quadprog", verbose = FALSE)
rdd.cate.qp = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, optimizer = "quadprog", verbose = FALSE)
rdd.free.mk = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, optimizer = "mosek", verbose = FALSE)
rdd.cate.mk = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, optimizer = "mosek", verbose = FALSE)
rdd.free.ecos = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, optimizer = "ECOS", verbose = FALSE)
rdd.cate.ecos = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, optimizer = "ECOS", verbose = FALSE)
expect_warning(rdd.free.scs <- optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, optimizer = "SCS", verbose = FALSE))
expect_warning(rdd.cate.scs <- optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, optimizer = "SCS", verbose = FALSE))

test_that("optimization strategies are equivalent", {
    expect_equal(rdd.free$tau.hat, rdd.free.raw$tau.hat, tolerance = rdd.free$tau.plusminus)
    expect_equal(rdd.free$tau.plusminus, rdd.free.raw$tau.plusminus, tolerance = 0.01)
    expect_equal(rdd.free$tau.hat, rdd.free.qp$tau.hat, tolerance = 0.01)
    expect_equal(rdd.free$tau.plusminus, rdd.free.qp$tau.plusminus, tolerance = 0.01)
    expect_equal(rdd.free$tau.hat, rdd.free.mk$tau.hat, tolerance = 0.001)
    expect_equal(rdd.free$tau.plusminus, rdd.free.mk$tau.plusminus, tolerance = 0.001)
    expect_equal(rdd.free$tau.hat, rdd.free.ecos$tau.hat, tolerance = 0.001)
    expect_equal(rdd.free$tau.plusminus, rdd.free.ecos$tau.plusminus, tolerance = 0.001)
    
    # SCS is not actually quite the same...
    expect_equal(rdd.free$tau.hat, rdd.free.scs$tau.hat, tolerance = 0.05)
    expect_equal(rdd.free$tau.plusminus, rdd.free.scs$tau.plusminus, tolerance = 0.05)
    
    expect_equal(rdd.cate$tau.hat, rdd.cate.raw$tau.hat, tolerance = rdd.cate$tau.plusminus)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.raw$tau.plusminus, tolerance = 0.05)
    expect_equal(rdd.cate$tau.hat, rdd.cate.qp$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.qp$tau.plusminus, tolerance = 0.05)
    expect_equal(rdd.cate$tau.hat, rdd.cate.mk$tau.hat, tolerance = 0.001)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.mk$tau.plusminus, tolerance = 0.001)
    expect_equal(rdd.cate$tau.hat, rdd.cate.ecos$tau.hat, tolerance = 0.005)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.ecos$tau.plusminus, tolerance = 0.005)
    
    # SCS is not actually quite the same...
    expect_equal(rdd.cate$tau.hat, rdd.cate.scs$tau.hat, tolerance = 0.05)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.scs$tau.plusminus, tolerance = 0.12)
})

# Test sigma square estimation for optrdd
rdd.fixed = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose=FALSE,
                   sigma.sq=1, use.homoskedatic.variance=FALSE)
rdd.homosk = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose=FALSE,
                    sigma.sq=1, use.homoskedatic.variance=TRUE)

test_that("oprdd gets variance almost right", {
    expect_equal(rdd.cate$tau.hat, rdd.fixed$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.fixed$tau.plusminus, tolerance = 0.05)
    expect_equal(rdd.cate$tau.hat, rdd.homosk$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.homosk$tau.plusminus, tolerance = 0.05)
})


# Test bias-adjusted confidence interval function

test_that("test plusminus function", {
    max.bias = 1
    se = 2
    alpha = 0.95
    pm = get.plusminus(max.bias, se, alpha)
    err = pnorm(-(pm + max.bias)/se) + pnorm(-(pm - max.bias)/se)
    expect_equal(alpha + err, 1, tolerance = 10^(-5))
    
    pm2 = get.plusminus(0, 1, 0.9)
    expect_equal(pm2, qnorm(0.95), tolerance = 10^(-5))
})

# Test 2d optrdd
X.2d = cbind(X, runif(n, -1, 1))
W = X.2d[, 1] < 0 | X.2d[, 2] < 0

MOSEK = requireNamespace("Rmosek", quietly = TRUE)
rdd.2d.free.mk = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                        verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "mosek")
rdd.2d.cate.mk = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                        verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "mosek")
# For quadprog, make the problem easier...
rdd.2d.free.qp = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                        verbose = FALSE, spline.df = 6, bin.width = 0.2, optimizer = "quadprog")
rdd.2d.cate.qp = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                        verbose = FALSE, spline.df = 6, bin.width = 0.2, optimizer = "quadprog")


rdd.2d.free.ecos = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "ECOS")
rdd.2d.cate.ecos = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "ECOS")
expect_warning(rdd.2d.free.scs <- optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                     verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "SCS"))
expect_warning(rdd.2d.cate.scs <- optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                     verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "SCS"))


test_that("2d-optrdd gammas satisfy constraints with mosek", {
    tol = 0.01
    expect_equal(sum(rdd.2d.free.mk$gamma), 0)
    expect_equal(sum(rdd.2d.free.mk$gamma * W), 1)
    expect_equal(sum(rdd.2d.free.mk$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.free.mk$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.mk$gamma), 0)
    expect_equal(sum(rdd.2d.cate.mk$gamma * W), 1)
    expect_equal(sum(rdd.2d.cate.mk$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.mk$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.mk$gamma * W * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.mk$gamma * W * X.2d[, 2]), 0, tolerance = tol)
})

test_that("2d-optrdd gammas satisfy constraints with quadprog", {
    tol = 0.05 # note the somewhate loose tolerance
    expect_equal(sum(rdd.2d.free.qp$gamma), 0)
    expect_equal(sum(rdd.2d.free.qp$gamma * W), 1)
    expect_equal(sum(rdd.2d.free.qp$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.free.qp$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.qp$gamma), 0)
    expect_equal(sum(rdd.2d.cate.qp$gamma * W), 1)
    expect_equal(sum(rdd.2d.cate.qp$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.qp$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.qp$gamma * W * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.qp$gamma * W * X.2d[, 2]), 0, tolerance = tol)
})

test_that("2d-optrdd gammas satisfy constraints with ECOS", {
    tol = 0.01
    expect_equal(sum(rdd.2d.free.ecos$gamma), 0)
    expect_equal(sum(rdd.2d.free.ecos$gamma * W), 1)
    expect_equal(sum(rdd.2d.free.ecos$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.free.ecos$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.ecos$gamma), 0)
    expect_equal(sum(rdd.2d.cate.ecos$gamma * W), 1)
    expect_equal(sum(rdd.2d.cate.ecos$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.ecos$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.ecos$gamma * W * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.ecos$gamma * W * X.2d[, 2]), 0, tolerance = tol)
})

test_that("2d-optrdd gammas satisfy constraints with SCS", {
    tol = 0.01
    expect_equal(sum(rdd.2d.free.scs$gamma), 0)
    expect_equal(sum(rdd.2d.free.scs$gamma * W), 1)
    expect_equal(sum(rdd.2d.free.scs$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.free.scs$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.scs$gamma), 0)
    expect_equal(sum(rdd.2d.cate.scs$gamma * W), 1)
    expect_equal(sum(rdd.2d.cate.scs$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.scs$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.scs$gamma * W * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate.scs$gamma * W * X.2d[, 2]), 0, tolerance = tol)
})

test_that("MOSEK/ECOS give same answer on 2d problem", {
    expect_equal(rdd.2d.free.mk$tau.hat, rdd.2d.free.ecos$tau.hat, tolerance = 0.001)
    expect_equal(rdd.2d.free.mk$tau.plusminus, rdd.2d.free.ecos$tau.plusminus, tolerance = 0.001)
    expect_equal(rdd.2d.cate.mk$tau.hat, rdd.2d.cate.ecos$tau.hat, tolerance = 0.001)
    expect_equal(rdd.2d.cate.mk$tau.plusminus, rdd.2d.cate.ecos$tau.plusminus, tolerance = 0.001)
})

test_that("MOSEK/SCS give same answer on 2d problem", {
    # note the looser tolerance
    expect_equal(rdd.2d.free.mk$tau.hat, rdd.2d.free.scs$tau.hat, tolerance = 0.05)
    expect_equal(rdd.2d.free.mk$tau.plusminus, rdd.2d.free.scs$tau.plusminus, tolerance = 0.01)
    expect_equal(rdd.2d.cate.mk$tau.hat, rdd.2d.cate.scs$tau.hat, tolerance = 0.05)
    expect_equal(rdd.2d.cate.mk$tau.plusminus, rdd.2d.cate.scs$tau.plusminus, tolerance = 0.01)
})

rdd.2d.free.raw.mk = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                         use.spline = FALSE, bin.width = 0.05, verbose = FALSE, optimizer = "mosek")
rdd.2d.free.raw.ecos = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                             use.spline = FALSE, bin.width = 0.05, verbose = FALSE, optimizer = "ECOS")
expect_warning(rdd.2d.free.raw.scs <- optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                                             use.spline = FALSE, bin.width = 0.05, verbose = FALSE, optimizer = "SCS"))

test_that("Spline approximation doesn't affect MOSEK", {
    expect_equal(rdd.2d.free.mk$tau.hat, rdd.2d.free.raw.mk$tau.hat, tolerance = rdd.2d.free.mk$tau.plusminus)
    expect_equal(rdd.2d.free.mk$tau.plusminus, rdd.2d.free.raw.mk$tau.plusminus, tolerance = 0.01)
})

test_that("Spline approximation doesn't affect ECOS", {
    expect_equal(rdd.2d.free.ecos$tau.hat, rdd.2d.free.raw.ecos$tau.hat, tolerance = rdd.2d.free.ecos$tau.plusminus)
    expect_equal(rdd.2d.free.ecos$tau.plusminus, rdd.2d.free.raw.ecos$tau.plusminus, tolerance = 0.01)
})


test_that("Spline approximation doesn't affect SCS", {
    expect_equal(rdd.2d.free.scs$tau.hat, rdd.2d.free.raw.scs$tau.hat, tolerance = rdd.2d.free.scs$tau.plusminus)
    expect_equal(rdd.2d.free.scs$tau.plusminus, rdd.2d.free.raw.scs$tau.plusminus, tolerance = 0.07)
})


test_that("baseline local linear regression implementation works", {
    
    skip_on_cran()
    source("../../baselines/local.lin.reg.R")
    rectangle = llr(X, Y = Y, max.second.derivative, kernel = "rectangular", 
                    minimization.target = "mse", max.window = 1)
    triangle = llr(X, Y = Y, max.second.derivative, kernel = "triangular", 
                   minimization.target = "mse", max.window = 1)
    
    half.bucket = min(rectangle$gamma.fun[-1, 1] -
                          rectangle$gamma.fun[-nrow(rectangle$gamma.fun), 1])
    expect_equal(sum(rectangle$gamma), 0)
    expect_equal(sum(rectangle$gamma * (X > 0)), 1)
    expect_equal(sum(rectangle$gamma * X), 0, tolerance = half.bucket)
    expect_equal(sum(rectangle$gamma * X * (X > 0)), 0, tolerance = half.bucket)
    
    half.bucket = min(triangle$gamma.fun[-1, 1] -
                          triangle$gamma.fun[-nrow(triangle$gamma.fun), 1])
    expect_equal(sum(triangle$gamma), 0)
    expect_equal(sum(triangle$gamma * (X > 0)), 1)
    expect_equal(sum(triangle$gamma * X), 0, tolerance = half.bucket)
    expect_equal(sum(triangle$gamma * X * (X > 0)), 0, tolerance = half.bucket)
    
    expect_true(rdd.cate$max.bias^2 + rdd.cate$sampling.se^2 <
                    triangle$max.bias^2 + triangle$sampling.se^2)
    expect_true(triangle$max.bias^2 + triangle$sampling.se^2 < 
                    rectangle$max.bias^2 + rectangle$sampling.se^2)
    
    # Test sigma square estimation for llr
    triangle.fixed = llr(X, Y = Y, max.second.derivative, sigma.sq = 1, max.window = 1, 
                         use.homoskedatic.variance = TRUE, kernel = "triangular",
                         minimization.target = "mse")
    
    expect_equal(triangle$tau.hat, triangle.fixed$tau.hat, tolerance = 0.05)
    expect_equal(triangle$tau.plusminus, triangle.fixed$tau.plusminus, 
                 tolerance = 0.05)
})

