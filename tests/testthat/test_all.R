set.seed(1)

max.second.derivative = 1
K = 20
supp = runif(K, -1, 1)
prob = rexp(K)
prob = prob/sum(prob)

n = 2000

bucket = as.numeric(1:K %*% rmultinom(n, 1, prob))

X = supp[bucket]
threshold = 0
W = as.numeric(X >= threshold)
Y = 10 + 20 * X + rnorm(n) + W

# Test methods initially, and confirm gamma moments

rdd.old = optrdd(X, Y=Y, max.second.derivative, max.window = 1)
rdd.free = optrdd.new(X=X, Y=Y, W=W, max.second.derivative = max.second.derivative, verbose = FALSE)
rdd.cate = optrdd.new(X=X, Y=Y, W=W, center = threshold, max.second.derivative = max.second.derivative, verbose = FALSE)

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


# Test sigma square estimation for optrdd
rdd.fixed = optrdd.new(X=X, Y=Y, W=W, center = threshold, max.second.derivative = max.second.derivative, verbose=FALSE,
                   sigma.sq=1, use.homoskedatic.variance=FALSE)
rdd.homosk = optrdd.new(X=X, Y=Y, W=W, center = threshold, max.second.derivative = max.second.derivative, verbose=FALSE,
                       sigma.sq=1, use.homoskedatic.variance=TRUE)

test_that("oprdd gets variance almost right", {
    expect_equal(rdd.cate$tau.hat, rdd.fixed$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.fixed$tau.plusminus, tolerance = 0.01)
    expect_equal(rdd.cate$tau.hat, rdd.homosk$tau.hat, tolerance = 0.01)
    expect_equal(rdd.cate$tau.plusminus, rdd.homosk$tau.plusminus, tolerance = 0.05)
})


# Test bias-adjusted confidence interval function

max.bias = 1
se = 2
alpha = 0.95
pm = get.plusminus(max.bias, se, alpha)
err = pnorm(-(pm + max.bias)/se) + pnorm(-(pm - max.bias)/se)

test_that("test plusminus function", {
    expect_equal(alpha + err, 1, tolerance = 10^(-5))
})

# Test 2d optrdd
X.2d = cbind(X, runif(n, -1, 1))
W = X.2d[, 1] < 0 | X.2d[, 2] < 0

rdd.2d.old = optrdd.2d(X=X.2d, Y=Y, max.second.derivative)
rdd.2d.free = optrdd.new(X=X.2d, Y=Y, W=W, max.second.derivative = max.second.derivative, verbose = FALSE)
rdd.2d.cate = optrdd.new(X=X.2d, Y=Y, W=W, center = c(0, 0), max.second.derivative = max.second.derivative, verbose = FALSE)

test_that("2d-optrdd gammas satisfy constraints", {
    tol = rdd.2d.cate$gamma.fun.0[2,1] - rdd.2d.cate$gamma.fun.0[1,1]
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

rdd.2d.center = optrdd.2d(X.2d, Y = Y, max.second.derivative, estimate.cate.at.point = TRUE, 
    center.treated.sample = TRUE)
test_that("centered 2d-optrdd gammas satisfy constraints", {
    expect_equal(sum(rdd.2d.center$gamma), 0)
    expect_equal(sum(rdd.2d.center$gamma * W), 1)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 1]), 0)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 2]), 0)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 1] * W), 0)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 2] * W), 0)
})




# Test baseline methods
rectangle = llr(X, Y = Y, max.second.derivative, kernel = "rectangular", 
                minimization.target = "mse", max.window = 1)

test_that("rectangular kernel gammas satisfy constraints", {
    half.bucket = min(rectangle$gamma.fun[-1, 1] -
                          rectangle$gamma.fun[-nrow(rectangle$gamma.fun), 1])/2
    expect_equal(sum(rectangle$gamma), 0)
    expect_equal(sum(rectangle$gamma * (X > 0)), 1)
    expect_equal(sum(rectangle$gamma * X), 0, tolerance = half.bucket)
    expect_equal(sum(rectangle$gamma * X * (X > 0)), 0, tolerance = half.bucket)
})


triangle = llr(X, Y = Y, max.second.derivative, kernel = "triangular", 
               minimization.target = "mse", max.window = 1)

test_that("triangular kernel gammas satisfy constraints", {
    half.bucket = min(triangle$gamma.fun[-1, 1] -
                          triangle$gamma.fun[-nrow(triangle$gamma.fun), 1])/2
    expect_equal(sum(triangle$gamma), 0)
    expect_equal(sum(triangle$gamma * (X > 0)), 1)
    expect_equal(sum(triangle$gamma * X), 0, tolerance = half.bucket)
    expect_equal(sum(triangle$gamma * X * (X > 0)), 0, tolerance = half.bucket)
})


test_that("relative performance of methods is as expected", {
    expect_true(rdd$tau.plusminus < triangle$tau.plusminus)
    expect_true(triangle$tau.plusminus < rectangle$tau.plusminus)
})


# Test aggregation for llr

triangle.agg = llr(X.agg, Y = Y.agg, max.second.derivative, kernel = "triangular", 
                   minimization.target = "mse", max.window = 1, num.samples = n.agg)

test_that("aggregation for llr is roughly the same", {
    expect_equal(triangle$tau.hat, triangle.agg$tau.hat, tolerance = 0.1)
    expect_equal(triangle$tau.plusminus, triangle.agg$tau.plusminus, tolerance = 0.1)
})

# Test sigma square estimation for llr

# If we do not estimate variance, then aggregation shouldn't do
# anything
triangle.fixed = llr(X, Y = Y, max.second.derivative, sigma.sq = 1, max.window = 1, 
                     use.homoskedatic.variance = TRUE, kernel = "triangular", minimization.target = "mse")
triangle.agg.fixed = llr(X.agg, Y = Y.agg, max.second.derivative, sigma.sq = 1, 
                         max.window = 1, num.samples = n.agg, use.homoskedatic.variance = TRUE, 
                         kernel = "triangular", minimization.target = "mse")

test_that("llr gets variance almost right", {
    expect_equal(triangle$tau.hat, triangle.fixed$tau.hat, tolerance = 0.01)
    expect_equal(triangle$tau.plusminus, triangle.fixed$tau.plusminus, 
                 tolerance = 0.01)
})

test_that("aggregation for llr is exact when variance is known", {
    expect_equal(triangle.fixed$tau.hat, triangle.agg.fixed$tau.hat)
    expect_equal(triangle.fixed$tau.plusminus, triangle.agg.fixed$tau.plusminus)
})


