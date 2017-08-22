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
print.llr = function(obj) {
    print.optrdd(obj)
}

#' @export
plot.llr = function(obj) {
    plot(obj$gamma.fun)
    abline(h = 0, lty = 3)
}

#' @export
summary.optrdd = function(obj) {
    unlist(obj)[1:length(obj)]
}

plot.optrdd = function(obj) {
    nvar = dim(obj$gamma.fun.0)[2] - 1
    if (nvar == 1) {
        
        all.x = c(obj$gamma.fun.0[,1], obj$gamma.fun.1[,1])
        xrng = range(all.x)
        yrng = range(obj$gamma)
        plot(NA, NA, xlab = "x", ylab = "gamma", xlim = xrng, ylim = yrng)
        if (length (unique(all.x) > 40)) {
            points(obj$gamma.fun.0, col = 4, pch = 16, cex = 1.5)
            points(obj$gamma.fun.1, col = 2, pch = 16, cex = 1.5)
        } else {
            lines(obj$gamma.fun.0, col = 4, lwd = 2)
            lines(obj$gamma.fun.1, col = 2, lwd = 2)
        }
        abline(h=0, lty = 3)
        
    } else if (nvar == 2) {
        
        if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
            stop("RColorBrewer needed for this function to work. Please install it.", 
                 call. = FALSE)
        }
        gamma.all = c(obj$gamma.fun.0[, 3], obj$gamma.fun.1[, 3])
        cidx = 51 + round(50 * gamma.all/max(abs(gamma.all)))
        hc = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(101)
        
        x1rng = range(obj$gamma.fun.0[, 1], obj$gamma.fun.1[, 1])
        x2rng = range(obj$gamma.fun.0[, 2], obj$gamma.fun.1[, 2])
        
        plot(NA, NA, xlim = x1rng, ylim = x2rng, pch = 16, xlab = "x1", ylab = "x2")
        points(obj$gamma.fun.0[, 1:2], col = hc[cidx[1:nrow(obj$gamma.fun.0)]], pch = 16)
        points(obj$gamma.fun.1[, 1:2], col = hc[cidx[nrow(obj$gamma.fun.0) + 1:nrow(obj$gamma.fun.1)]], pch = 1, lwd = 2)
    
    } else {
        stop("Corrupted object.")
    }
}

#' @export
plot.optrdd.2d = function(obj, xlab = "x1", ylab = "x2") {

}

#' @export
summary.optrdd.2d = function(obj) {
    summary.optrdd(obj)
}
