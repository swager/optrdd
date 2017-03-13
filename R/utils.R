get.plusminus = function(max.bias, sampling.se, alpha = 0.95) {
  rel.bias = max.bias / sampling.se
  zz = uniroot(function(z) pnorm(rel.bias - z) + pnorm(-rel.bias - z) - 0.05, c(0, rel.bias - qnorm((1 - alpha) / 2)))$root
  zz * sampling.se
}

plot.optrdd = function(obj) {
  plot(obj$gamma.fun)
  abline(h=0, lty = 3)
}

print.optrdd = function(obj) {
  if (!is.null(obj$tau.hat)) {
    print(paste0(100 * obj$alpha, "% CI for tau: ", signif(obj$tau.hat, 2), " +/- ", signif(obj$tau.plusminus, 2)))
  } else {
    print(paste0(100 * obj$alpha, "% CI for tau: [point estimate] +/- ", signif(obj$tau.plusminus, 2)))
  }
}

print.optrdd.2d = function(obj) {
  print.optrdd(obj)
}
  