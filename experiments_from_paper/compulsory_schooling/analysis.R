set.seed(1234)

rm(list = ls())

#library(foreign)
detach("package:optrdd", unload=TRUE)
library(optrdd)
source("~/git/optrdd/baselines/local.lin.reg.R")

# print more
print.optrdd = function (obj) 
{
  if (!is.null(obj$tau.hat)) {
    print(paste0(100 * obj$alpha, "% CI for tau: $", round(obj$tau.hat, 4),
                 " \\pm ", round(obj$tau.plusminus, 4), "$"))
  }
  else {
    print(paste0(100 * obj$alpha, "% CI for tau: [point estimate] +/- ", 
                 round(obj$tau.plusminus, 3)))
  }
}
print.llr = function(obj) {
    print.optrdd(obj)
}


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# filtering as in kolesar & rothe
# cghs = foreign::read.dta("combined general household survey.dta")
# cghs$yearat14 <- cghs$yobirth+14
# d <- within(cghs, {
#   k1 <- yearat14 >= 35 & yearat14 <= 65 & age <= 64 &
#     agelfted >= 10 & !is.na(agelfted);
#   k2 <- k1 & !is.na(earn) & nireland==0;
#   k3 <- k2 & abs(yearat14-47) <= 6
#   k4 <- k2 & abs(yearat14-47) <= 3
#   learn <- log(earn)
# })
# d$x <- d$yearat14
# 
# ## save subset of GB data
# data <- d[d$k2, c("learn", "x", "agelfted", "sex", "datyear")]
# write.csv(data, "uk_analysis_sample.csv")
data <- read.csv("uk_analysis_sample.csv")

threshold <- 46.99
max.window <- 12.1
datsub <- data[which(data$x - threshold <= max.window),]
W = as.numeric(datsub$x > threshold)

# repro K & R
rect.003 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.003,
    threshold = threshold, num.bucket = 400,
    minimization.target = "ci.length")
rect.006 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.006,
               threshold = threshold, num.bucket = 400,
               minimization.target = "ci.length")
rect.012 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.012,
               threshold = threshold, num.bucket = 400,
               minimization.target = "ci.length")
rect.03 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.03,
    threshold = threshold, num.bucket = 400,
    minimization.target = "ci.length")

# now with triangular kernel
tri.003 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.003, kernel = "triangular",
           threshold = threshold, num.bucket = 400,
           minimization.target = "ci.length")
tri.006 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.006, kernel = "triangular",
              threshold = threshold, num.bucket = 400,
              minimization.target = "ci.length")
tri.012 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.012, kernel = "triangular",
              threshold = threshold, num.bucket = 400,
              minimization.target = "ci.length")
tri.03 = llr(datsub$x, Y=datsub$learn, max.second.derivative = 0.03, kernel = "triangular",
              threshold = threshold, num.bucket = 400,
             minimization.target = "ci.length")

opt.003 = optrdd(datsub$x, Y=datsub$learn, W=W, max.second.derivative = 0.003,
                 estimation.point = threshold, num.bucket = 400)
opt.006 = optrdd(datsub$x, Y=datsub$learn, max.second.derivative = 0.006,
                 estimation.point = threshold, W=W, num.bucket = 400)
opt.012 = optrdd(datsub$x, Y=datsub$learn, max.second.derivative = 0.012,
                 estimation.point = threshold, W=W, num.bucket = 400)
opt.03 = optrdd(datsub$x, Y=datsub$learn, max.second.derivative = 0.03,
                estimation.point = threshold, W=W, num.bucket = 400)

rect.003
rect.03
tri.003
tri.03
opt.003
opt.03

cat(paste(substr(print(rect.003), 17, 35),
  substr(print(tri.003), 17, 35),
  substr(print(opt.003), 17, 35), sep = " & "))

cat(paste(substr(print(rect.006), 17, 35),
          substr(print(tri.006), 17, 35),
          substr(print(opt.006), 17, 35), sep = " & "))

cat(paste(substr(print(rect.012), 17, 35),
          substr(print(tri.012), 17, 35),
          substr(print(opt.012), 17, 35), sep = " & "))

cat(paste(substr(print(rect.03), 17, 35),
          substr(print(tri.03), 17, 35),
          substr(print(opt.03), 17, 35), sep = " & "))

#
# This is used to pick the curvature bound B. Note that this regression
# should not be interpreted as an estimate of the RDD parameter.
#

W = datsub$x >= threshold
X = (datsub$x - threshold)
X2 = (X - threshold)^2
Y = datsub$learn
RDF = data.frame(W=W, X=X, X2=X2, Y=Y)

summary(lm(Y ~ W * (X + X2), data = RDF))


COLS = RColorBrewer::brewer.pal(9, "Set1")
pdf("oreopoulos_weights.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim=(47 + 6 * c(-1, 1)), pch = 4, cex=2, xlab = "Year turned 14", ylab = "gamma", ylim = range(c(opt.012$gamma.fun.0[,2], opt.012$gamma.fun.1[,2], tri.012$gamma.fun[,2], rect.012$gamma.fun[,2])))
points(opt.012$gamma.fun.0[,1], opt.012$gamma.fun.0[,2], pch = 4, lwd=3, col = COLS[2], cex=2)
points(opt.012$gamma.fun.1[,1], opt.012$gamma.fun.1[,2], pch = 4, lwd=3, col = COLS[2], cex=2)
lines(rect.012$gamma.fun[rect.012$gamma.fun[,1] < threshold,1], rect.012$gamma.fun[rect.012$gamma.fun[,1] < threshold,2], lwd = 3, col = COLS[5])
lines(rect.012$gamma.fun[rect.012$gamma.fun[,1] >= threshold,1], rect.012$gamma.fun[rect.012$gamma.fun[,1] >= threshold,2], lwd = 3, col = COLS[5])
points(tri.012$gamma.fun, pch = 3, cex=2, col = COLS[4], lwd = 3)
abline(h = 0, lty = 3, lwd = 2)
legend("bottomright", c("Optimized", "Tri Kenrel", "Rect Kernel"), lwd = 3, col = COLS[c(2, 4, 5)], pch = c(4, 3, NA), lty = c(NA, NA, 1), cex = 2)
par=pardef
dev.off()


source("armstrong_kolesar_gamma.R")
ak.003 = ak.tau(datsub$x, Y=datsub$learn, max.second.derivative = 0.003,
                threshold = threshold, num.bucket = 400)
ak.006 = ak.tau(datsub$x, Y=datsub$learn, max.second.derivative = 0.006,
                threshold = threshold, num.bucket = 400)
ak.012 = ak.tau(datsub$x, Y=datsub$learn, max.second.derivative = 0.012,
                threshold = threshold, num.bucket = 400)
ak.03 = ak.tau(datsub$x, Y=datsub$learn, max.second.derivative = 0.03,
                threshold = threshold, num.bucket = 400)


ak = sapply(list(ak.003, ak.006, ak.012, ak.03), print)




