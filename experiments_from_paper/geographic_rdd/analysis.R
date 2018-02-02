rm(list=ls())

set.seed(1234)

library(optrdd)
library(foreign)
library(RColorBrewer)
library(xtable)
library(glmnet)
library(splines)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pointsALL <- read.dbf("KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf")
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X

data.voting <- read.dta("KeeleTitiunik2014-PA-replication-files/Data/Voters/Voters_Final.dta",convert.dates = FALSE, convert.factors = FALSE, convert.underscore = FALSE)
data.voting.nona = data.voting[!is.na(data.voting$treat),]

plot(data.voting.nona$longitude,
     data.voting.nona$latitude,
     col=data.voting.nona$treat + 1,  pch = ".", cex = 3)
lines(pointsALL$longitude, pointsALL$latitude, lwd = 4)

WV = data.voting.nona$treat
XV = (pi * 6371 /180) * cbind(
    A=data.voting.nona$longitude - mean(data.voting.nona$longitude),
    B=data.voting.nona$latitude - mean(data.voting.nona$latitude))

hc = ((grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu")))(101))[101:1]
purple = ((grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Purples")))(101))

compute_curvature = function(xgrid.pred, nA, nB, centers, binw) {
    curvs0 = sapply(centers, function(idx) {
        c((xgrid.pred[idx + 1] + xgrid.pred[idx - 1] - 2 * xgrid.pred[idx]) / binw^2,
          (xgrid.pred[idx + nA] + xgrid.pred[idx - nA] - 2 * xgrid.pred[idx]) / binw^2,
          (xgrid.pred[idx + nA + 1] + xgrid.pred[idx - nA - 1] - 2 * xgrid.pred[idx]) / binw^2 / 2,
          (xgrid.pred[idx + nA - 1] + xgrid.pred[idx - nA + 1] - 2 * xgrid.pred[idx]) / binw^2 / 2)
    })
    curvs = apply(curvs0, 2, function(iii) max(abs(iii)))
    quantile(curvs, 0.95, na.rm=TRUE)
}

get_curvature = function(xx, outcome, ww, make_plot = TRUE, binw = 0.1) {
    
    xx = data.frame(xx)
    yy = data.voting.nona[,outcome]
    
    names(xx) = c("A", "B")
    
    gridA = seq(min(xx$A) - 1.5 * binw, max(xx$A) + 1.5 * binw, by = binw)
    gridB = seq(min(xx$B) - 1.5 * binw, max(xx$B) + 1.5 * binw, by = binw)
    xgrid = data.frame(expand.grid(A=gridA, B = gridB))
    xspl.all = model.matrix(~ 0 + ns(A, df = 7) * ns(B, df = 7),
                            data = rbind(xx, xgrid))
    xspl = xspl.all[1:nrow(xx),]
    
    fit0 = cv.glmnet(xspl[ww==0,], yy[ww==0], alpha = 0)
    fit1 = cv.glmnet(xspl[ww==1,], yy[ww==1], alpha = 0)
    
    if (make_plot) {
        pred = rep(NA, nrow(xx))
        pred[ww==0] = predict(fit0, xspl[ww==0,], s="lambda.1se")
        pred[ww==1] = predict(fit1, xspl[ww==1,], s="lambda.1se")
        print(paste(outcome, range(pred)))
        cidx = 1 + round(100 * (pred - min(pred))/(max(abs(pred)) - min(pred)))
        pdf(paste0("output/raw_pred_", outcome, ".pdf"))
        pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        plot(NA, NA, xlim = range(data.voting.nona$longitude), ylim = range(data.voting.nona$latitude),
             xlab = "longitude", ylab = "latitude")
        points(data.voting.nona$longitude, data.voting.nona$latitude, col = purple[cidx], pch = 16)
        lines(pointsALL$longitude, pointsALL$latitude, lwd = 2, lty = 1)
        par = pardef
        dev.off()
    }

    xgrid.spl = xspl.all[nrow(xx) + 1:nrow(xgrid),]
    xgrid.pred.0 = predict(fit0, xgrid.spl, s="lambda.1se")
    xgrid.pred.1 = predict(fit1, xgrid.spl, s="lambda.1se")
    
    nA = length(gridA)
    nB = length(gridB)
    
    bucketA = as.numeric(cut(xx[,1], gridA))
    bucketB = as.numeric(cut(xx[,2], gridB))
    bucket = bucketA + nA * bucketB
    
    c.hat.0 = compute_curvature(xgrid.pred.0, nA, nB, centers = bucket[ww == 0], binw = binw)
    c.hat.1 = compute_curvature(xgrid.pred.1, nA, nB, centers = bucket[ww == 1], binw = binw)
    max(c.hat.0, c.hat.1)
}

#
# Run analysis
#

#outcomes = c("e2008p", "age", "black", "hisp", "dem", "female")
outcomes = c("e2008p", "age", "black", "dem", "female")

BOUNDARY_PT = 37
estimation.point = (pi * 6371 /180) *
    c(pointsALL$longitude[BOUNDARY_PT] - mean(data.voting.nona$longitude),
      pointsALL$latitude[BOUNDARY_PT] - mean(data.voting.nona$latitude))

res.all =  lapply(outcomes, function(outcome) {
    YV = data.voting.nona[,outcome]
    max.curv = get_curvature(XV, outcome, WV)
    print(max.curv)
    out = optrdd(X=XV, Y=YV, W=WV, max.second.derivative = max.curv)
    out.pt = optrdd(X=XV, Y=YV, W=WV, max.second.derivative = max.curv, estimation.point = estimation.point)
    return(list(out, out.pt))
})

#
# Make results table
#

res.avg = lapply(res.all, function(x) x[[1]])
res.pt = lapply(res.all, function(x) x[[2]])

parsed = cbind(c("WATE", rep("", 4)), outcomes,
               t(sapply(res.avg, function(out) {
                   c(sprintf("%.3f", round(out$tau.hat, 3)),
                     paste0("(", sprintf("%.3f", round(out$tau.hat - out$tau.plusminus, 3)),
                            ", ", sprintf("%.3f", round(out$tau.hat + out$tau.plusminus, 3))
                            , ")"),
                     sprintf("%.3f", round(out$max.bias, 3)), 
                     sprintf("%.3f", round(out$sampling.se, 3)))
               })))

parsed.pt = cbind(c("CATE", rep("", 4)), outcomes,
                      t(sapply(res.pt, function(out) {
                          c(sprintf("%.3f", round(out$tau.hat, 3)),
                            paste0("(", sprintf("%.3f", round(out$tau.hat - out$tau.plusminus, 3)),
                                   ", ", sprintf("%.3f", round(out$tau.hat + out$tau.plusminus, 3))
                                   , ")"),
                            sprintf("%.3f", round(out$max.bias, 3)), 
                            sprintf("%.3f", round(out$sampling.se, 3)))
                      })))

tab.all = rbind(parsed, parsed.pt)
colnames(tab.all) = c("", "outcome", "point estimate", "conf interval", "max bias", "sampling std err")
xtab = xtable(tab.all, align=c("r", "|", "r", "r", "|", "c", "c", "c", "c", "|"))
print(xtab, include.rownames = FALSE,
      hline.after = c(-1, 0, 5, 10),
      file="output/geo_out.tex")


#
# Make gamma plots
#

main.res = res.all[[1]]

out = main.res[[1]]
gamma.all = c(out$gamma.fun.0[, 3], out$gamma.fun.1[, 3])
cidx = 51 + round(50 * gamma.all/max(abs(gamma.all)))

lon.all = c(out$gamma.fun.0[,1], out$gamma.fun.1[,1]) /
    (pi * 6371 /180) + mean(data.voting.nona$longitude)
lat.all = c(out$gamma.fun.0[,2], out$gamma.fun.1[,2]) /
    (pi * 6371 /180) + mean(data.voting.nona$latitude)

control = 1:nrow(out$gamma.fun.0)
treat = nrow(out$gamma.fun.0) + 1:nrow(out$gamma.fun.1)

pdf("output/geo_WATE.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim = range(lon.all), ylim = range(lat.all),
     xlab = "longitude", ylab = "latitude")
points(lon.all[control], lat.all[control], col = hc[cidx[control]], pch = 10, lwd = 1.5)
points(lon.all[treat], lat.all[treat], col = hc[cidx[treat]], pch = 16, lwd = 1.5)
lines(pointsALL$longitude, pointsALL$latitude, lwd = 2, lty = 1)
par = pardef
dev.off()


out.pt = main.res[[2]]

gamma.all.pt = c(out.pt$gamma.fun.0[, 3], out.pt$gamma.fun.1[, 3])
cidx.pt = 51 + round(50 * gamma.all.pt/max(abs(gamma.all.pt)))

lon.all.pt = c(out.pt$gamma.fun.0[,1], out.pt$gamma.fun.1[,1]) /
    (pi * 6371 /180) + mean(data.voting.nona$longitude)
lat.all.pt = c(out.pt$gamma.fun.0[,2], out.pt$gamma.fun.1[,2]) /
    (pi * 6371 /180) + mean(data.voting.nona$latitude)

control.pt = 1:nrow(out.pt$gamma.fun.0)
treat.pt = nrow(out.pt$gamma.fun.0) + 1:nrow(out.pt$gamma.fun.1)

pdf("output/geo_CATE.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim = range(lon.all.pt), ylim = range(lat.all.pt),
     xlab = "longitude", ylab = "latitude")
points(lon.all.pt[control.pt], lat.all.pt[control.pt],
       col = hc[cidx.pt[control.pt]], pch = 10, lwd = 1.5)
points(lon.all.pt[treat.pt], lat.all.pt[treat.pt],
       col = hc[cidx.pt[treat.pt]], pch = 16, lwd = 1.5)
lines(pointsALL$longitude, pointsALL$latitude, lwd = 2, lty = 1)
points(pointsALL$longitude[BOUNDARY_PT], pointsALL$latitude[BOUNDARY_PT],
       lwd = 4, cex = 1.5, pch = 4)
par = pardef
dev.off()
