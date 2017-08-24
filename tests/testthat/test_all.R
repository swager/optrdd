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
rdd.free = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, verbose = FALSE)
rdd.cate = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose = FALSE)

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

test_that("optimization strategies are equivalent", {
    expect_equal(rdd.free$tau.hat, rdd.free.raw$tau.hat, tolerance = rdd.free$tau.plusminus)
    expect_equal(rdd.free$tau.plusminus, rdd.free.raw$tau.plusminus, tolerance = 0.01)
    expect_equal(rdd.free$tau.hat, rdd.free.qp$tau.hat, tolerance = 0.01)
    expect_equal(rdd.free$tau.plusminus, rdd.free.qp$tau.plusminus, tolerance = 0.01)
    expect_equal(rdd.free$tau.hat, rdd.free.mk$tau.hat, tolerance = 0.01)
    expect_equal(rdd.free$tau.plusminus, rdd.free.mk$tau.plusminus, tolerance = 0.01)
    
    expect_equal(rdd.cate$tau.hat, rdd.cate.raw$tau.hat, tolerance = rdd.cate$tau.plusminus)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.raw$tau.plusminus, tolerance = 0.05)
    expect_equal(rdd.cate$tau.hat, rdd.cate.qp$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.qp$tau.plusminus, tolerance = 0.05)
    expect_equal(rdd.cate$tau.hat, rdd.cate.mk$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.mk$tau.plusminus, tolerance = 0.05)
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
if (MOSEK) {
    rdd.2d.free = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 20, bin.width = 0.05)
    rdd.2d.cate = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 20, bin.width = 0.05)
} else { # if MOSEK isn't installed, make problem easier
    rdd.2d.free = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 6, bin.width = 0.2)
    rdd.2d.cate = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 6, bin.width = 0.2)
}


test_that("2d-optrdd gammas satisfy constraints", {
    tol = 0.05
    expect_equal(sum(rdd.2d.free$gamma), 0)
    expect_equal(sum(rdd.2d.free$gamma * W), 1)
    expect_equal(sum(rdd.2d.free$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.free$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate$gamma), 0)
    expect_equal(sum(rdd.2d.cate$gamma * W), 1)
    expect_equal(sum(rdd.2d.cate$gamma * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate$gamma * X.2d[, 2]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate$gamma * W * X.2d[, 1]), 0, tolerance = tol)
    expect_equal(sum(rdd.2d.cate$gamma * W * X.2d[, 2]), 0, tolerance = tol)
})

test_that("optimization strategies are equivalent", {
    skip_if_not(MOSEK) # too slow without MOSEK
    rdd.2d.free.raw = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                             use.spline = FALSE, bin.width = 0.05, verbose = FALSE)
    expect_equal(rdd.2d.free$tau.hat, rdd.2d.free.raw$tau.hat, tolerance = rdd.2d.free$tau.plusminus)
    expect_equal(rdd.2d.free$tau.plusminus, rdd.2d.free.raw$tau.plusminus, tolerance = 0.01)
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

