rm(list = ls())

set.seed(2)

detach("package:optrdd", unload=TRUE)
library(optrdd)
source("~/git/optrdd/baselines/local.lin.reg.R")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

COLS = RColorBrewer::brewer.pal(9, "Set1")

max.second.derivative = 5
K = 40
supp = rnorm(K)
prob = rexp(K)
prob = prob/sum(prob)

nn.all = round(1000 * exp(((-3):12 / 3) * log(3)))
nn.plot = 1000 * 3^(0:3)

maxerr.list = lapply(nn.all, function(n) {
#maxerr.list = lapply(nn.plot, function(n) {
  
  bucket = as.numeric(1:K %*% rmultinom(n, 1, prob))
  X = supp[bucket]

  X.short=X[abs(X) <= 1]  
  W.short=as.numeric(X.short > 0)
  rdd = optrdd(X=X.short, W=W.short, max.second.derivative=max.second.derivative, estimation.point = 0, sigma.sq=1, num.bucket = 600)
  rectangle = llr(X.short, max.second.derivative, sigma.sq=1, kernel="rectangular", num.bucket = 600)
  triangle = llr(X.short, max.second.derivative, sigma.sq=1, kernel="triangular", num.bucket = 600)
  
  maxerr = data.frame(bias=c(rdd$max.bias, rectangle$max.bias, triangle$max.bias),
                      se=c(rdd$sampling.se, rectangle$sampling.se, triangle$sampling.se))
  rownames(maxerr) = c("rdd", "rectangle", "triangle")
  
  if (n %in% nn.plot) {
    
    X.g = rnorm(n)
    X.g = X.g[abs(X.g) <= 1]
    W.g = as.numeric(X.g > 0)
    rdd.g = optrdd(X=X.g, W=W.g, max.second.derivative=max.second.derivative, sigma.sq=1, estimation.point = 0, bin.width = 2/600)
    
    pdf(paste0("plots/first_", n, ".pdf"))
    plot(NA, NA, ylim = c(-5.5, 5.5), xlim = c(-0.7, 0.7), xlab = "", ylab = "", yaxt='n', cex.axis=1.5)
    points(rdd$gamma.fun.0[,1], rdd$gamma.fun.0[,2] * n^(4/5), pch = 16, cex=2, col=COLS[1])
    points(rdd$gamma.fun.1[,1], rdd$gamma.fun.1[,2] * n^(4/5), pch = 16, cex=2, col=COLS[1])
    lines(rdd.g$gamma.fun.0[,1], rdd.g$gamma.fun.0[,2] * n^(4/5), lwd = 3, col=COLS[2])
    lines(rdd.g$gamma.fun.1[,1], rdd.g$gamma.fun.1[,2] * n^(4/5), lwd=3, col=COLS[2])
    abline(h=0, lty = 2, lwd = 3)
    abline(v=0, lty = 2, lwd = 1)
    dev.off()
  }
  
  return(maxerr)
  
})

out = Reduce(cbind, lapply(maxerr.list, function(mm) rowSums(mm^2)))

pdf("plots/mse_plot.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(nn.all, out[2,]/out[1,], log="x", ylim=range(c(1, range(out[2,]/out[1,]), range(out[3,]/out[1,]))), xlab="n", ylab = "Relative Excess Error", pch = 24, col = COLS[5], bg = COLS[5], cex = 1.5)
points(nn.all, out[3,]/out[1,], col = COLS[4], pch = 25, cex = 1.5, bg=COLS[4])
abline(h=1, lty = 2, lwd = 2, col = 1)
abline(h=1.05, lty = 4, lwd = 1)
abline(h=1.1, lty = 4, lwd = 1)
abline(h=1.15, lty = 4, lwd = 1)
abline(h=1.2, lty = 4, lwd = 1)
legend("topleft", c("Rectangular Kernel", "Triangular Kernel"), pch = c(24, 25), col = COLS[5:4], pt.bg = COLS[5:4], cex = 1.5, bg="white")
par=pardef
dev.off()

pdf("plots/intro_pmf.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim = c(-0.7, 0.7), ylim = c(0, 0.1), xlab = "X", ylab = "Probability Mass")
for(iter in 1:K) {
  if(abs(supp[iter]) < 0.75) {
    segments(supp[iter], 0, supp[iter], prob[iter], lwd = 4, col = COLS[7])
  }
}
abline(h=0, lwd = 4)
abline(v=0, lty = 2)
par=pardef
dev.off()



