set.seed(1)

ludwig.miller = Vectorize(function(x) {
    if (x < 0) {
        mu = 3.71 + 2.3 * x + 3.28 * x^2 + 1.45 * x^3 + 0.23 * x^4 + 0.03 * x^5
    } else {
        mu = 0.26  + 18.49 * x - 54.81 * x^2 + 74.30 * x^3 - 45.02 * x^4 + 9.83 * x^5
    }
    mu
})
B.ludwig.miller = 2 * 54.81
tau.ludwig.miller = - 3.71 + 0.26

fff = ludwig.miller
Bmax = B.ludwig.miller
tau.true = tau.ludwig.miller

n = 5000
sigma = 0.1295

X = 2 * rbeta(n, 2, 4) - 1
Y = fff(X) + sigma * rnorm(n)
W = as.numeric(X > 0)

test_that("Elastic net is recommended when appropriate.", {
    # Since there is a lot of curvature, we need to use an elastic net to
    # get good CIs (and avoid a warning)
    expect_warning(run1 <- optrdd(X, Y, W, Bmax, optimizer = "mosek", verbose = FALSE))
    expect_silent(run2 <- optrdd(X, Y, W, Bmax, optimizer = "mosek", verbose = FALSE, try.elnet.for.sigma.sq = TRUE))
    expect_true(run2$tau.plusminus < 0.97 * run1$tau.plusminus)
    expect_true(abs(run2$tau.hat - tau.true) < run2$tau.plusminus)
})

test_that("Optimizer warnings are emitted.", {
    # ECOS is unable to solve accurately here. There should be a warning.
    expect_warning(run3 <- optrdd(X, Y, W, Bmax, optimizer = "ECOS", verbose = FALSE, try.elnet.for.sigma.sq = TRUE))
    # SCS sometimes doesn't find the true optimum. There should be a warning about this.
    expect_warning(run4 <- optrdd(X, Y, W, Bmax, optimizer = "SCS", verbose = FALSE, try.elnet.for.sigma.sq = TRUE))
})