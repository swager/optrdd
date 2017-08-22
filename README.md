# optrdd
Optimized inference in regression discontinuity designs, as proposed by Imbens and Wager (2017).

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("swager/optrdd")
```

This package currently works with two optimizers: `mosek` an `quadprog`.
Mosek is a commercial interior point solver that needs to be installed separately,
while quadprog is a standard `R` optimization library.
Both optimizers appear to work well with a univariate running variable;
however, with multi-dimensional running variables, we strongly recommend
[installing mosek](INSTALLING_MOSEK.md).

Example usage:

```R
library(optrdd)

# Simple regression discontinuity with discrete X
n = 4000; threshold = 0
X = sample(seq(-4, 4, by = 8/41.5), n, replace = TRUE)
W = as.numeric(X >= threshold)
Y = 0.4 * W + 1 / (1 + exp(2 * X)) + 0.2 * rnorm(n)
# using 0.4 for max.second.derivative would have been enough
out.1 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5, estimation.point = threshold)
print(out.1); plot(out.1, xlim = c(-1.5, 1.5))

# Now, treatment is instead allocated in a neighborhood of 0
thresh.low = -1; thresh.high = 1
W = as.numeric(thresh.low <= X & X <= thresh.high)
Y = 0.2 * (1 + X) * W + 1 / (1 + exp(2 * X)) + rnorm(n)
# This estimates CATE at specifically chosen points
out.2 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5, estimation.point = thresh.low)
print(out.2); plot(out.2, xlim = c(-2.5, 2.5))
out.3 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5, estimation.point = thresh.high)
print(out.3); plot(out.3, xlim = c(-2.5, 2.5))
# This estimates a weighted CATE, with lower variance
out.4 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 0.5)
print(out.4); plot(out.4, xlim = c(-2.5, 2.5))

# RDD with multivariate running variable. Warning: slow without mosek.
X = matrix(runif(n*2, -1, 1), n, 2)
W = as.numeric(X[,1] < 0 | X[,2] < 0)
Y = X[,1]^2/3 + W * (1 + X[,2])
out.5 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 1)
print(out.5); plot(out.5)
out.6 = optrdd(X=X, Y=Y, W=W, max.second.derivative = 1, estimation.point = c(0, 0.5))
print(out.6); plot(out.6)
```

#### References
Guido Imbens and Stefan Wager.
<b>Optimized Regression Discontinuity Designs.</b> 2017.
[<a href="https://arxiv.org/pdf/1705.01677.pdf">arxiv</a>]