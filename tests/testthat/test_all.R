set.seed(1)

max.second.derivative = 5
K = 1000
supp = rnorm(K)
prob = rexp(K)
prob = prob/sum(prob)

n = 2000

bucket = as.numeric(1:K %*% rmultinom(n, 1, prob))
X = supp[bucket]
Y = 10 + 20 * X + rnorm(n)

# Test methods initially, and confirm gamma moments

rdd = optrdd(X, Y = Y, max.second.derivative, max.window = 1)

test_that("optrdd gammas satisfy constraints", {
    half.bucket = min(rdd$gamma.fun[-1, 1] -
                        rdd$gamma.fun[-nrow(rdd$gamma.fun), 1])/2
    expect_equal(sum(rdd$gamma), 0)
    expect_equal(sum(rdd$gamma * (X > 0)), 1)
    expect_equal(sum(rdd$gamma * X), 0, tolerance = half.bucket)
    expect_equal(sum(rdd$gamma * X * (X > 0)), 0, tolerance = half.bucket)
})

rectangle = llr(X, Y = Y, max.second.derivative, kernel = "rectangular", 
    minimization.target = "mse", max.window = 1)

test_that("rectangular kernel gammas satisfy constraints", {
    expect_equal(sum(rectangle$gamma), 0)
    expect_equal(sum(rectangle$gamma * (X > 0)), 1)
    expect_equal(sum(rectangle$gamma * X), 0)
    expect_equal(sum(rectangle$gamma * X * (X > 0)), 0)
})


triangle = llr(X, Y = Y, max.second.derivative, kernel = "triangular", 
    minimization.target = "mse", max.window = 1)

test_that("triangular kernel gammas satisfy constraints", {
    expect_equal(sum(triangle$gamma), 0)
    expect_equal(sum(triangle$gamma * (X > 0)), 1)
    expect_equal(sum(triangle$gamma * X), 0)
    expect_equal(sum(triangle$gamma * X * (X > 0)), 0)
})


test_that("relative performance of methods is as expected", {
    expect_true(rdd$tau.plusminus < triangle$tau.plusminus)
    expect_true(triangle$tau.plusminus < rectangle$tau.plusminus)
})

# Test aggregation for optrdd

X.agg = unique(X)
Y.agg = sapply(X.agg, function(x) mean(Y[X == x]))
n.agg = sapply(X.agg, function(x) sum(X == x))

rdd.agg = optrdd(X.agg, Y = Y.agg, max.second.derivative, max.window = 1, 
    num.samples = n.agg)

test_that("aggregation for optrdd is roughly the same", {
    expect_equal(rdd$tau.hat, rdd.agg$tau.hat, tolerance = 0.1)
    expect_equal(rdd$tau.plusminus, rdd.agg$tau.plusminus, tolerance = 0.1)
})

# Test sigma square estimation for optrdd

# If we do not estimate variance, then aggregation shouldn't do
# anything
rdd.fixed = optrdd(X, Y = Y, max.second.derivative, sigma.sq = 1, max.window = 1, 
    use.homoskedatic.variance = TRUE)
rdd.agg.fixed = optrdd(X.agg, Y = Y.agg, max.second.derivative, sigma.sq = 1, 
    max.window = 1, num.samples = n.agg, use.homoskedatic.variance = TRUE)

test_that("oprdd gets variance almost right", {
    expect_equal(rdd$tau.hat, rdd.fixed$tau.hat, tolerance = 0.01)
    expect_equal(rdd$tau.plusminus, rdd.fixed$tau.plusminus, tolerance = 0.01)
})

test_that("aggregation for optrdd is exact when variance is known", {
    expect_equal(rdd.fixed$tau.hat, rdd.agg.fixed$tau.hat)
    expect_equal(rdd.fixed$tau.plusminus, rdd.agg.fixed$tau.plusminus)
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

max.window = c(1, 1)
X.2d = cbind(X, rnorm(n))
W = X.2d[, 1] < 0 | X.2d[, 2] < 0

rdd.2d = optrdd.2d(X.2d, Y = Y, max.second.derivative, max.window = max.window)
test_that("2d-optrdd gammas satisfy constraints", {
    expect_equal(sum(rdd.2d$gamma), 0)
    expect_equal(sum(rdd.2d$gamma * W), 1)
    expect_equal(sum(rdd.2d$gamma * X.2d[, 1]), 0)
    expect_equal(sum(rdd.2d$gamma * X.2d[, 2]), 0)
})

rdd.2d.center = optrdd.2d(X.2d, Y = Y, max.second.derivative, estimate.cate.at.point = TRUE, 
    center.treated.sample = TRUE, max.window = max.window)
test_that("centered 2d-optrdd gammas satisfy constraints", {
    expect_equal(sum(rdd.2d.center$gamma), 0)
    expect_equal(sum(rdd.2d.center$gamma * W), 1)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 1]), 0)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 2]), 0)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 1] * W), 0)
    expect_equal(sum(rdd.2d.center$gamma * X.2d[, 2] * W), 0)
})


