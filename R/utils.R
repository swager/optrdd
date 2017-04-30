#' Bias-adjusted Gaussian confidence intervals.
#'
#' @param max.bias Worst-case bias of estimate.
#' @param sampling.se Sampling error of estimate.
#' @param alpha Coverage probability of confidence interval.
#'
#' @return Half-width of confidence interval.
#' @export
get.plusminus = function(max.bias, sampling.se, alpha = 0.95) {
    rel.bias = max.bias/sampling.se
    zz = uniroot(function(z) pnorm(rel.bias - z) +
                   pnorm(-rel.bias - z) - 0.05,
                 c(0, rel.bias - qnorm((1 - alpha)/3)))$root
    zz * sampling.se
}

#' @export
plot.optrdd = function(obj) {
    plot(obj$gamma.fun)
    abline(h = 0, lty = 3)
}

#' @export
print.optrdd = function(obj) {
    if (!is.null(obj$tau.hat)) {
        print(paste0(100 * obj$alpha, "% CI for tau: ",
                     signif(obj$tau.hat, 2), " +/- ", signif(obj$tau.plusminus, 2)))
    } else {
        print(paste0(100 * obj$alpha, "% CI for tau: [point estimate] +/- ", 
            signif(obj$tau.plusminus, 2)))
    }
}

#' @export
summary.optrdd = function(obj) {
    unlist(obj)[1:7]
}

#' @export
print.optrdd.2d = function(obj) {
    print.optrdd(obj)
}

#' @export
plot.optrdd.2d = function(obj, xlab = "x1", ylab = "x2") {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
        stop("RColorBrewer needed for this function to work. Please install it.", 
            call. = FALSE)
    }
    gamma.xx = -obj$gamma.fun[, 3]
    cidx = 51 + round(50 * gamma.xx/max(abs(gamma.xx)))
    hc = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(101)
    plot(obj$gamma.fun[, 1:2], col = hc[cidx], pch = 16, xlab = xlab, 
        ylab = ylab)
    segments(2 * max(obj$gamma.fun[, 1]), 0, 0, 0)
    segments(0, 0, 0, 2 * max(obj$gamma.fun[, 2]))
    points(obj$tau.center[1], obj$tau.center[2], pch = 4, 
        cex = 1.5, lwd = 3)
}

#' @export
summary.optrdd.2d = function(obj) {
    summary.optrdd(obj)
}
