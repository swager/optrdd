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
rdd.free.qp = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, verbose = FALSE, optimizer = "quadprog")
rdd.cate.qp = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose = FALSE, optimizer = "quadprog")

rdd.free.ecos = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, verbose = FALSE, optimizer = "ECOS")
rdd.cate.ecos = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose = FALSE, optimizer = "ECOS")

test_that("optrdd gammas satisfy constraints with quadprog", {
    tol = rdd.cate.qp$gamma.fun.0[2, 1] - rdd.cate.qp$gamma.fun.0[1, 1]
    expect_equal(sum(rdd.cate.qp$gamma), 0)
    expect_equal(sum(rdd.cate.qp$gamma * W), 1)
    expect_equal(sum(rdd.cate.qp$gamma * X), 0, tolerance = tol)
    expect_equal(sum(rdd.cate.qp$gamma * X * W), 0, tolerance = tol)
    expect_equal(sum(rdd.free.qp$gamma), 0)
    expect_equal(sum(rdd.free.qp$gamma * W), 1)
    expect_equal(sum(rdd.free.qp$gamma * X), 0, tolerance = tol)
})

test_that("optrdd gammas satisfy constraints with ECOS", {
    tol = rdd.cate.ecos$gamma.fun.0[2, 1] - rdd.cate.ecos$gamma.fun.0[1, 1]
    expect_equal(sum(rdd.cate.ecos$gamma), 0)
    expect_equal(sum(rdd.cate.ecos$gamma * W), 1)
    expect_equal(sum(rdd.cate.ecos$gamma * X), 0, tolerance = tol)
    expect_equal(sum(rdd.cate.ecos$gamma * X * W), 0, tolerance = tol)
    expect_equal(sum(rdd.free.ecos$gamma), 0)
    expect_equal(sum(rdd.free.ecos$gamma * W), 1)
    expect_equal(sum(rdd.free.ecos$gamma * X), 0, tolerance = tol)
})

test_that("ECOS and quadprog optimizers match", {
    expect_equal(rdd.cate.qp$tau.hat, rdd.cate.ecos$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate.qp$tau.plusminus, rdd.cate.ecos$tau.plusminus, tolerance = 0.01)
})

test_that("cate constraint hurts", {
    expect_true(rdd.cate.qp$tau.plusminus > rdd.free.qp$tau.plusminus)
    expect_true(rdd.cate.ecos$tau.plusminus > rdd.free.ecos$tau.plusminus)
})


# Test optimization strategies

rdd.free = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, optimizer = "quadprog", verbose = FALSE)
rdd.cate = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, optimizer = "quadprog", verbose = FALSE)
rdd.free.ecos = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, optimizer = "ECOS", verbose = FALSE)
rdd.cate.ecos = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, optimizer = "ECOS", verbose = FALSE)
rdd.free.scs = optrdd(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, optimizer = "SCS", verbose = FALSE)
rdd.cate.scs = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, optimizer = "SCS", verbose = FALSE)

test_that("optimization strategies are equivalent", {

    expect_equal(rdd.free$tau.hat, rdd.free.ecos$tau.hat, tolerance = 0.01)
    expect_equal(rdd.free$tau.plusminus, rdd.free.ecos$tau.plusminus, tolerance = 0.01)
    # SCS is too inaccurate, so this test gets flaky
    #expect_equal(rdd.free$tau.hat, rdd.free.scs$tau.hat, tolerance = 0.2)
    #expect_equal(rdd.free$tau.plusminus, rdd.free.scs$tau.plusminus, tolerance = 0.2)
    
    expect_equal(rdd.cate$tau.hat, rdd.cate.ecos$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.cate.ecos$tau.plusminus, tolerance = 0.01)
    #expect_equal(rdd.cate$tau.hat, rdd.cate.scs$tau.hat, tolerance = 0.2)
    #expect_equal(rdd.cate$tau.plusminus, rdd.cate.scs$tau.plusminus, tolerance = 0.2)
})

# Test sigma square estimation for optrdd
rdd.fixed = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose=FALSE,
                   sigma.sq=1, use.homoskedatic.variance=FALSE, optimizer = "ECOS")
rdd.homosk = optrdd(X=X, Y=Y, W=W, estimation.point = threshold, max.second.derivative = max.second.derivative, verbose=FALSE,
                    sigma.sq=1, use.homoskedatic.variance=TRUE, optimizer = "ECOS")

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

# For quadprog, make the problem easier...
rdd.2d.free.qp = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                        verbose = FALSE, spline.df = 6, bin.width = 0.2, optimizer = "quadprog")
rdd.2d.cate.qp = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                        verbose = FALSE, spline.df = 6, bin.width = 0.2, optimizer = "quadprog")


rdd.2d.free.ecos = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "ECOS")
rdd.2d.cate.ecos = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                         verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "ECOS")
rdd.2d.free.scs = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                     verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "SCS")
rdd.2d.cate.scs = optrdd(X=X.2d, Y=Y, W=W, estimation.point = c(0, 0), max.second.derivative = max.second.derivative,
                     verbose = FALSE, spline.df = 20, bin.width = 0.05, optimizer = "SCS")

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
    tol = 0.03
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
    tol = 0.03
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

test_that("ECOS/SCS give same answer on 2d problem", {
    # note the looser tolerance
    expect_equal(rdd.2d.free.ecos$tau.hat, rdd.2d.free.scs$tau.hat, tolerance = 0.1)
    expect_equal(rdd.2d.free.ecos$tau.plusminus, rdd.2d.free.scs$tau.plusminus, tolerance = 0.05)
    expect_equal(rdd.2d.cate.ecos$tau.hat, rdd.2d.cate.scs$tau.hat, tolerance = 0.1)
    expect_equal(rdd.2d.cate.ecos$tau.plusminus, rdd.2d.cate.scs$tau.plusminus, tolerance = 0.05)
})

rdd.2d.free.raw.ecos = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
                             use.spline = FALSE, bin.width = 0.05, verbose = FALSE, optimizer = "ECOS")
test_that("Spline approximation doesn't affect ECOS", {
    expect_equal(rdd.2d.free.ecos$tau.hat, rdd.2d.free.raw.ecos$tau.hat, tolerance = rdd.2d.free.ecos$tau.plusminus)
    expect_equal(rdd.2d.free.ecos$tau.plusminus, rdd.2d.free.raw.ecos$tau.plusminus, tolerance = 0.01)
})

# SCS is again a little flaky
# rdd.2d.free.raw.scs = optrdd(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative,
#                              use.spline = FALSE, bin.width = 0.05, verbose = FALSE, optimizer = "SCS")
# test_that("Spline approximation doesn't affect SCS", {
#     expect_equal(rdd.2d.free.scs$tau.hat, rdd.2d.free.raw.scs$tau.hat, tolerance = rdd.2d.free.scs$tau.plusminus)
#     expect_equal(rdd.2d.free.scs$tau.plusminus, rdd.2d.free.raw.scs$tau.plusminus, tolerance = 0.1)
# })

