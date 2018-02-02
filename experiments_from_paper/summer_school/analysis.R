rm(list = ls())

set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#detach("package:optrdd", unload=TRUE)
library(optrdd)

library(RColorBrewer)

data.all = read.csv("ssextract_karthik.csv")[,-1]
hist(data.all$mdcut01)
hist(data.all$rdcut01)

THRESH = 40

data = data.all[pmax(abs(data.all$mdcut01), abs(data.all$rdcut01)) <= THRESH,]
data$itt = data$mdcut01 > 0 & data$rdcut01 > 0

table(data$mdcut01 > 0, data$rdcut01 > 0)
round(mean(data$ssatyn01), 3)
round(mean(data$ssatyn01[data$mdcut01 < 0 & data$rdcut01 < 0]), 3)
round(mean(data$ssatyn01[data$mdcut01 > 0 & data$rdcut01 < 0]), 3)
round(mean(data$ssatyn01[data$mdcut01 < 0 & data$rdcut01 > 0]), 3)
round(mean(data$ssatyn01[data$mdcut01 > 0 & data$rdcut01 > 0]), 3)

data.sm = data.all[pmax(abs(data.all$mdcut01), abs(data.all$rdcut01)) < 20,]

data.sm$pass = data.sm$mdcut01 > 0 & data.sm$rdcut01 > 0
summary(lm(zmscr02 ~ mdcut01 + rdcut01 + pass, data = data.sm))

X = as.numeric(data$mdcut01)
Y = as.numeric(data$zmscr02)

uu = unique(X)
yy = sapply(uu, function(uuu) mean(Y[X == uuu]))
plot(uu, yy)

X1 = as.numeric(data$mdcut01) / THRESH
X2 = as.numeric(data$rdcut01) / THRESH
X = cbind(X1, X2)
Y.math = as.numeric(data$zmscr02)
Y.reading = as.numeric(data$zrscr02)
W = as.numeric(X[,1] <= 0 | X[,2] <= 0)

threshold = c(0, 0)
max.window = c(1, 1)
num.bucket = c(40, 40)

# Guess at max second derivative
DF = data.frame(Y=Y.reading, X1=X1, X1.2=(X1 - mean(X1))^2, X2=X2, X2.2=(X2 - mean(X2))^2, X12=(X1 - mean(X1))*(X2 - mean(X2)), W=as.numeric(X1 < 0 | X2 < 0))
lmb = coef(lm(Y ~ W * ., data = DF))
M0.curv = matrix(c(2*lmb[4], lmb[7], lmb[7], 2*lmb[6]), 2, 2)
M1.curv = M0.curv + matrix(c(2*lmb[9], lmb[12], lmb[12], 2*lmb[11]), 2, 2)
svd(M0.curv)$d
svd(M1.curv)$d

# Biggest curvature effects:
# 
# For math, among treated (i.e., summer school) sample, curvature of -0.2 in the
# (8, 5) direction (i.e., summer school maybe doesn't help good students,
# esp. students good at math?)
# 
# For reading, among controls (no summer school) sample, curvature of +0.46
# in the (1, 2) direction (i.e., good readers improve on their own?).

subjects = c("math", "reading")
max.derivs = c(0.5, 1)
cate.at.pts = c(TRUE, FALSE)
#centers = c(TRUE, FALSE)

curr.idx = 1
summaries = list()

for (subject in subjects) {
  for (max.second.derivative in max.derivs) {
    #for (center in centers) {
      for (cate.at.pt in cate.at.pts) {
        
        center = cate.at.pt
        if (!center & cate.at.pt) next;
        
        if (subject == "math") {
          Y = Y.math
        } else {
          Y = Y.reading
        }
        
        if (cate.at.pt) {
            estimation.point = threshold
        } else {
            estimation.point = NULL
        }
        gamma = optrdd(X = X, Y = Y, W = W,
                        max.second.derivative = max.second.derivative,
                        estimation.point = estimation.point)
        print(gamma)
        
        
        
        pdf(paste0("output/gamma_", subject, "_B_", 10 * max.second.derivative,
                   "_cate_", cate.at.pt, "_center_", center, ".pdf"))
        #plot(gamma, xlab="math score", ylab="reading score")
        x=gamma
        gamma.all = c(x$gamma.fun.0[, 3], x$gamma.fun.1[, 3])
        cidx = 51 + round(50 * gamma.all/max(abs(gamma.all)))
        hc = (grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu")))(101)[101:1]
        x1rng = range(x$gamma.fun.0[, 1], x$gamma.fun.1[, 1])
        x2rng = range(x$gamma.fun.0[, 2], x$gamma.fun.1[, 2])
        
        pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        plot(NA, NA, xlim = x1rng, ylim = x2rng,
             xlab="math score", ylab="reading score")
        points(x$gamma.fun.0[, 1], x$gamma.fun.0[, 2],
               col = hc[cidx[1:length(x$gamma.fun.0[, 3])]], pch = 10, lwd = 1.5)
        points(x$gamma.fun.1[, 1], x$gamma.fun.1[, 2],
               col = hc[cidx[length(x$gamma.fun.0[, 3]) + 1:length(x$gamma.fun.1[, 3])]], pch = 16, lwd = 1.5)
        segments(0, 0, 0, 2, lwd = 2)
        segments(0, 0, 2, 0, lwd = 2)
        
        if (cate.at.pt) {
            points(estimation.point[1], estimation.point[2], lwd = 4, cex = 1.5, pch = 4)
        } else {
            middle = colSums(X[W == 1,] * gamma$gamma[W==1])
            points(middle[1], middle[2], lwd = 4, cex = 1.25, pch = 5)
        }
        
        par = pardef
        dev.off()
        
        
        
        save(gamma, file=paste0("output/object_", subject, "_B_", max.second.derivative,
                   "_cate_", cate.at.pt, "_center_", center, ".RData"))
        
        summaries[[curr.idx]] = c(subject=subject,
                                  max.second.derivative=max.second.derivative,
                                  cate.at.pt=cate.at.pt,
                                  center=cate.at.pt,
                                  summary(gamma))
        curr.idx = curr.idx + 1
      }
    #}
  }
}

result_summaries = data.frame(Reduce(rbind, summaries))
write.csv(result_summaries, file="output/result_summaries.csv")
