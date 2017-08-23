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
    zz = stats::uniroot(function(z) stats::pnorm(rel.bias - z) +
                            stats::pnorm(-rel.bias - z) + alpha - 1,
                        c(0, rel.bias - stats::qnorm((1 - alpha)/3)))$root
    zz * sampling.se
}

#' @export
print.optrdd = function(x, ...) {
    if (!is.null(x$tau.hat)) {
        print(paste0(100 * x$alpha, "% CI for tau: ",
                     signif(x$tau.hat, 2), " +/- ", signif(x$tau.plusminus, 2)))
    } else {
        print(paste0(100 * x$alpha, "% CI for tau: [point estimate] +/- ", 
                     signif(x$tau.plusminus, 2)))
    }
}

#' @export
summary.optrdd = function(object, ...) {
    unlist(object)[1:5]
}

#' @export
plot.optrdd = function(x, ...) {
    
    nvar = dim(x$gamma.fun.0)[2] - 1
    args = list(...)
    
    if (nvar == 1) {
        
        all.x = c(x$gamma.fun.0[,1], x$gamma.fun.1[,1])
        if (!"xlim" %in% names(args)) {
            args$xlim = range(all.x)
        }
        if (!"ylim" %in% names(args)) {
            args$ylim = range(x$gamma)
        }
        if (!"xlab" %in% names(args)) {
            args$xlab = "x"
        }
        if (!"ylab" %in% names(args)) {
            args$ylab = "gamma"
        }
        args$x = NA; args$y = NA
        do.call(graphics::plot, args)
        if (length (unique(all.x) > 40)) {
            graphics::points(x$gamma.fun.0, col = 4, pch = 16, cex = 1.5)
            graphics::points(x$gamma.fun.1, col = 2, pch = 16, cex = 1.5)
        } else {
            graphics::lines(x$gamma.fun.0, col = 4, lwd = 2)
            graphics::lines(x$gamma.fun.1, col = 2, lwd = 2)
        }
        graphics::abline(h=0, lty = 3)
        
    } else if (nvar == 2) {
        
        if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
            stop("RColorBrewer needed for this function to work. Please install it.", 
                 call. = FALSE)
        }
        gamma.all = c(x$gamma.fun.0[, 3], x$gamma.fun.1[, 3])
        cidx = 51 + round(50 * gamma.all/max(abs(gamma.all)))
        hc = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(101)
        
        x1rng = range(x$gamma.fun.0[, 1], x$gamma.fun.1[, 1])
        x2rng = range(x$gamma.fun.0[, 2], x$gamma.fun.1[, 2])
        
        if (!"xlim" %in% names(args)) {
            args$xlim = x1rng
        }
        if (!"ylim" %in% names(args)) {
            args$ylim = x2rng
        }
        if (!"xlab" %in% names(args)) {
            args$xlab = "x1"
        }
        if (!"ylab" %in% names(args)) {
            args$ylab = "x2"
        }
        args$x = NA; args$y = NA
        do.call(graphics::plot, args)
        graphics::points(x$gamma.fun.0[, 1:2],
                         col = hc[cidx[1:nrow(x$gamma.fun.0)]],
                         pch = 16)
        graphics::points(x$gamma.fun.1[, 1:2],
                         col = hc[cidx[nrow(x$gamma.fun.0) + 1:nrow(x$gamma.fun.1)]],
                         pch = 1, lwd = 2)
    
    } else {
        stop("Corrupted object.")
    }
}
